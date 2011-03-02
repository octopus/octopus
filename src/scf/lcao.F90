!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module lcao_m
  use batch_m
  use blacs_proc_grid_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use h_sys_output_m
  use lalg_adv_m
  use lalg_basic_m
  use lapack_m
  use loct_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2
  use parser_m
  use profiling_m
  use simul_box_m
  use scalapack_m
  use solids_m
  use species_m
  use species_pot_m
  use states_m
  use states_block_m
  use states_calc_m
  use states_dim_m
  use states_io_m
  use submesh_m
  use system_m
  use v_ks_m
  use varinfo_m

  implicit none

  private
  public ::            &
    lcao_t,            &
    lcao_init,         &
    lcao_wf,           &
    lcao_run,          &
    lcao_end,          &
    lcao_is_available, &
    lcao_num_orbitals

  integer, public, parameter ::     &
    LCAO_START_NONE    = 0, &
    LCAO_START_STATES  = 2, &
    LCAO_START_FULL    = 3

  type lcao_t
    private
    logical           :: initialized !< are k, s and v1 matrices filled?
    integer           :: norbs !< number of orbitals
    integer           :: maxorbs
    integer, pointer  :: atom(:)
    integer, pointer  :: level(:)
    integer, pointer  :: ddim(:)
    logical           :: alternative
    logical           :: derivative
    
    ! For the alternative LCAO
    logical             :: keep_orb     !< Whether we keep orbitals in memory.
    FLOAT,   pointer    :: radius(:)    !< The localization radius of each atom orbitals
    FLOAT               :: lapdist      !< This is the extra distance that the Laplacian adds to the localization radius.
    integer             :: mult         !< The number of basis per atomic function (with derivatives is 2, 1 otherwise).
    integer             :: maxorb       !< The maximum value of the orbitals over all atoms.
    integer             :: nbasis       !< The total number of basis functions.
    ! The following functions map between a basis index and atom/orbital index
    integer, pointer    :: basis_atom(:) !< The atom that corresponds to a certain basis index
    integer, pointer    :: basis_orb(:)  !< The orbital that corresponds to a certain basis index
    integer, pointer    :: atom_orb_basis(:, :) !< The basis index that coorrespond to a certain
    integer, pointer    :: norb_atom(:)  !< The number of orbitals per atom including mult.
    logical             :: parallel      !< Whether the LCAO is done in parallel
    integer             :: lsize(1:2)
    integer             :: nproc(1:2)
    integer             :: myroc(1:2)
    integer             :: desc(1:BLACS_DLEN)
    logical, pointer    :: calc_atom(:)
    FLOAT               :: diag_tol
  end type lcao_t
  
  type(profile_t), save :: prof_orbitals

contains

  ! ---------------------------------------------------------
  subroutine lcao_init(this, gr, geo, st)
    type(lcao_t),         intent(out)   :: this
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(in)    :: st

    integer :: ia, n, ii, jj, maxj, idim

    PUSH_SUB(lcao_init)

    ! nullify everything so we can check for associated pointers when deallocating
    nullify(this%atom)
    nullify(this%level)
    nullify(this%ddim)

    this%initialized = .true.

    !%Variable LCAOAlternative
    !%Type logical
    !%Default false
    !%Section SCF::LCAO
    !%Description
    !% If this variable is set, the LCAO procedure will use an
    !% alternative (and experimental) implementation. It is faster for
    !% large systems and parallel in states.
    !%End
    call parse_logical(datasets_check('LCAOAlternative'), .false., this%alternative)

    if(.not. this%alternative) then

      ! count the number of orbitals available
      maxj = 0
      this%maxorbs = 0
      do ia = 1, geo%natoms
        maxj = max(maxj, species_niwfs(geo%atom(ia)%spec) )
        this%maxorbs = this%maxorbs + species_niwfs(geo%atom(ia)%spec)
      end do

      this%maxorbs = this%maxorbs*st%d%dim

      if(this%maxorbs < st%nst) then
        this%initialized = .false.
        write(message(1),'(a)') 'Cannot do LCAO initial calculation because there are not enough atomic orbitals.'
        write(message(2),'(a,i6,a,i6,a)') 'Required: ', st%nst, '. Available: ', this%maxorbs, '.'
        call write_warning(2)
        POP_SUB(lcao_init)
        return
      end if

      ! generate tables to know which indices each atomic orbital has

      SAFE_ALLOCATE( this%atom(1:this%maxorbs))
      SAFE_ALLOCATE(this%level(1:this%maxorbs))
      SAFE_ALLOCATE( this%ddim(1:this%maxorbs))

      ! Each atom provides niwfs pseudo-orbitals (this number is given in
      ! geo%atom(ia)%spec%niwfs for atom number ia). This number is
      ! actually multiplied by two in case of spin-unrestricted or spinors
      ! calculations.
      !
      ! The pseudo-orbitals are placed in order in the following way (Natoms
      ! is the total number of atoms).
      !
      ! n = 1 => first orbital of atom 1,
      ! n = 2 => first orbital of atom 2.
      ! n = 3 => first orbital of atom 3.
      ! ....
      ! n = Natoms => first orbital of atom Natoms
      ! n = Natoms + 1 = > second orbital of atom 1
      ! ....
      !
      ! If at some point in this loop an atom pseudo cannot provide the corresponding
      ! orbital (because the niws orbitals have been exhausted), it moves on to the following
      ! atom.
      !
      ! In the spinors case, it changes a bit:
      !
      ! n = 1 => first spin-up orbital of atom 1, assigned to the spin-up component of the spinor.
      ! n = 2 => first spin-down orbital of atom 1, assigned to the spin-down component of the spinor.
      ! n = 3 => first spin-up orbital of atom 2, assigned to the spin-up component of the spinor.

      ii = 1
      do jj = 1, maxj
        do ia = 1, geo%natoms
          do idim = 1,st%d%dim
            if(jj > species_niwfs(geo%atom(ia)%spec) ) cycle

            this%atom(ii) = ia
            this%level(ii) = jj
            this%ddim(ii) = idim

            ii = ii + 1
          end do
        end do
      end do

      ASSERT(ii - 1 == this%maxorbs)

      !%Variable LCAODimension
      !%Type integer
      !%Default 0
      !%Section SCF::LCAO
      !%Description
      !% Before starting the SCF cycle, an initial LCAO calculation can be performed
      !% in order to obtain reasonable initial guesses for spin-orbitals and densities.
      !% For this purpose, the code calculates a number of atomic orbitals -- this
      !% number depends on the given species. The default dimension for the LCAO basis
      !% set will be the sum of all these numbers, unless this dimension is larger than
      !% twice the number of required orbitals for the full calculation. 
      !%
      !% This dimension however can be changed by making use of this
      !% variable. Note that <tt>LCAODimension</tt> cannot be smaller than the
      !% number of orbitals needed in the full calculation -- if
      !% <tt>LCAODimension</tt> is smaller, it will be silently increased to meet
      !% this requirement. In the same way, if <tt>LCAODimension</tt> is larger
      !% than the available number of atomic orbitals, it will be
      !% reduced. If you want to use the largest possible number, set
      !% <tt>LCAODimension</tt> to a negative number.
      !%End
      call parse_integer(datasets_check('LCAODimension'), 0, n)

      if(n > 0 .and. n <= st%nst) then
        this%norbs = st%nst
      else if(n > st%nst .and. n <= this%maxorbs) then
        this%norbs = n
      else if(n == 0) then
        this%norbs = min(this%maxorbs, 2*st%nst)
      else
        this%norbs = this%maxorbs
      end if

      ASSERT(this%norbs >= st%nst)
      ASSERT(this%norbs <= this%maxorbs)
      
      nullify(this%radius)
      nullify(this%basis_atom)
      nullify(this%basis_orb)
      nullify(this%atom_orb_basis)
      nullify(this%norb_atom)
      nullify(this%calc_atom)
    else
      call lcao2_init()
    end if

    POP_SUB(lcao_init)

  contains
    subroutine lcao2_init()
      integer :: iatom, iorb, norbs
      FLOAT   :: maxradius
      integer :: ibasis
#ifdef HAVE_SCALAPACK
      integer :: jatom, jorb, jbasis, ilbasis, jlbasis, proc(1:2), info, nbl
#endif
      PUSH_SUB(lcao_init.lcao2_init)

      message(1) = "Info: Using LCAO alternative implementation."
      call write_info(1)

      call messages_experimental('LCAO alternative implementation')

      !%Variable LCAOKeepOrbitals
      !%Type logical
      !%Default yes
      !%Section SCF::LCAO
      !%Description
      !% If set to yes (the default) Octopus keeps atomic arbitals in
      !% memory during the LCAO procedure. If set to no, the orbitals
      !% are generated each time that they are needed, increasing
      !% computational time but saving memory.
      !%
      !% When set to yes, Octopus prints the amount of memory per node
      !% that is required to store the orbitals.
      !%
      !%End
      call parse_logical(datasets_check('LCAOKeepOrbitals'), .true., this%keep_orb)

      !%Variable LCAOExtraOrbitals
      !%Type logical
      !%Default false
      !%Section SCF::LCAO
      !%Description
      !% (experimental) If this variable is set to yes, the LCAO
      !% procedure will add an extra set of numerical orbitals (by
      !% using the derivative of the radial part of the original
      !% orbitals).
      !%End
      call parse_logical(datasets_check('LCAOExtraOrbitals'), .false., this%derivative)

      if(this%derivative) call messages_experimental('LCAO extra orbitals')

      !%Variable LCAODiagTol
      !%Type float
      !%Default 1e-10
      !%Section SCF::LCAO
      !%Description
      !% The tolerance for the diagonalization of the LCAO Hamiltonian. The default is 1e-10.
      !%End
      call parse_float(datasets_check('LCAODiagTol'), CNST(1e-10), this%diag_tol)

      if(this%derivative) then
        this%mult = 2
      else
        this%mult = 1
      end if

      SAFE_ALLOCATE(this%norb_atom(1:geo%natoms))

      this%maxorb = 0
      this%nbasis = 0
      do iatom = 1, geo%natoms
        this%norb_atom(iatom) = this%mult*species_niwfs(geo%atom(iatom)%spec)
        this%maxorb = max(this%maxorb, species_niwfs(geo%atom(iatom)%spec))
        this%nbasis = this%nbasis + species_niwfs(geo%atom(iatom)%spec)
      end do

      this%maxorb = this%maxorb*this%mult
      this%nbasis = this%nbasis*this%mult

      SAFE_ALLOCATE(this%basis_atom(1:this%nbasis))
      SAFE_ALLOCATE(this%basis_orb(1:this%nbasis))
      SAFE_ALLOCATE(this%atom_orb_basis(1:geo%natoms, 1:this%maxorb))

      ! Initialize the mapping between indexes

      ibasis = 0
      do iatom = 1, geo%natoms
        norbs = species_niwfs(geo%atom(iatom)%spec)

        do iorb = 1, this%mult*norbs
          ibasis = ibasis + 1
          this%atom_orb_basis(iatom, iorb) = ibasis
          this%basis_atom(ibasis) = iatom
          this%basis_orb(ibasis) = iorb
        end do
      end do

      ! this is determined by the stencil we are using and the spacing
      this%lapdist = maxval(abs(gr%mesh%idx%enlarge)*gr%mesh%spacing)

      ! calculate the radius of each orbital
      SAFE_ALLOCATE(this%radius(1:geo%natoms))

      do iatom = 1, geo%natoms
        norbs = species_niwfs(geo%atom(iatom)%spec)

        maxradius = M_ZERO
        do iorb = 1, norbs
          maxradius = max(maxradius, species_get_iwf_radius(geo%atom(iatom)%spec, iorb, is = 1))
        end do

        if(this%derivative) maxradius = maxradius + this%lapdist

        maxradius = min(maxradius, M_TWO*maxval(gr%mesh%sb%lsize(1:gr%mesh%sb%dim)))

        this%radius(iatom) = maxradius
      end do

      SAFE_ALLOCATE(this%calc_atom(1:geo%natoms))
      this%calc_atom = .true.

      ! initialize parallel data
#ifndef HAVE_SCALAPACK
      this%parallel = .false.
#else
      this%parallel = (st%parallel_in_states .or. gr%mesh%parallel_in_domains) &
        .and. .not. blacs_proc_grid_null(st%dom_st_proc_grid)

      if(this%parallel) then      
        nbl = min(16, this%nbasis)

        ! The size of the distributed matrix in each node
        this%lsize(1) = max(1, numroc(this%nbasis, nbl, st%dom_st_proc_grid%myrow, 0, st%dom_st_proc_grid%nprow))
        this%lsize(2) = max(1, numroc(this%nbasis, nbl, st%dom_st_proc_grid%mycol, 0, st%dom_st_proc_grid%npcol))

        this%nproc(1) = st%dom_st_proc_grid%nprow
        this%nproc(2) = st%dom_st_proc_grid%npcol
        this%myroc(1) = st%dom_st_proc_grid%myrow
        this%myroc(2) = st%dom_st_proc_grid%mycol

        call descinit(this%desc(1), this%nbasis, this%nbasis, nbl, nbl, 0, 0, &
          st%dom_st_proc_grid%context, this%lsize(1), info)

        ASSERT(info == 0)

        this%calc_atom = .false.
        do iatom = 1, geo%natoms
          ibasis = this%atom_orb_basis(iatom, 1)

          do jatom = 1, geo%natoms
            jbasis = this%atom_orb_basis(jatom, 1)

            do iorb = 1, this%norb_atom(iatom)
              do jorb = 1, this%norb_atom(jatom)
                call lcao_local_index(this,  ibasis - 1 + iorb,  jbasis - 1 + jorb, &
                  ilbasis, jlbasis, proc(1), proc(2))

                this%calc_atom(this%basis_atom(jbasis)) = &
                  this%calc_atom(this%basis_atom(jbasis)) .or. proc(2) == this%myroc(2)

              end do
            end do

          end do
        end do

      end if
#endif

      POP_SUB(lcao_init.lcao2_init)
    end subroutine lcao2_init

  end subroutine lcao_init


  ! ---------------------------------------------------------
  subroutine lcao_run(sys, hm, st_start)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    integer, optional,   intent(in)    :: st_start ! use for unoccupied-states run

    integer :: lcao_start_default, lcao_start
    type(lcao_t) :: lcao
    integer :: s1, s2, k1, k2, is, ik, ip, idim
    logical :: lcao_done

    PUSH_SUB(lcao_run)

    if (.not. present(st_start)) then
      call guess_density(sys%gr%fine%mesh, sys%gr%sb, sys%geo, sys%st%qtot, sys%st%d%nspin, &
        sys%st%d%spin_channels, sys%st%rho)

      ! set up Hamiltonian (we do not call system_h_setup here because we do not want to
      ! overwrite the guess density)
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)

      ! get the effective potential (we don`t need the eigenvalues yet)
      call v_ks_calc(sys%ks, hm, sys%st, calc_eigenval=.false., calc_berry=.false.)
      ! eigenvalues have nevertheless to be initialized to something
      sys%st%eigenval = M_ZERO
    else
      call v_ks_calc(sys%ks, hm, sys%st, calc_eigenval=.true.)

      if(st_start .gt. sys%st%nst) then ! nothing to be done in LCAO
        POP_SUB(lcao_run)
        return
      endif
    endif

    ! The initial LCAO calculation is done by default if we have pseudopotentials.
    ! Otherwise, it is not the default value and has to be enforced in the input file.
    lcao_start_default = LCAO_START_FULL
    if(sys%geo%only_user_def) lcao_start_default = LCAO_START_NONE
    
    !%Variable LCAOStart
    !%Type integer
    !%Section SCF
    !%Description
    !% Before starting a SCF calculation, <tt>Octopus</tt> can perform
    !% a LCAO calculation. These can provide <tt>Octopus</tt> with a good set
    !% of initial wavefunctions and with a new guess for the density.
    !% (Up to the current version, only a minimal basis set is used.)
    !% The default is <tt>lcao_full</tt> unless all species are user-defined, in which case
    !% the default is <tt>lcao_none</tt>.
    !%Option lcao_none 0
    !% Do not perform a LCAO calculation before the SCF cycle. Instead use random wavefunctions.
    !%Option lcao_states 2
    !% Do a LCAO calculation before the SCF cycle and use the resulting wavefunctions as 
    !% initial wavefunctions without changing the guess density.
    !% This will speed up the convergence of the eigensolver during the first SCF iterations.
    !%Option lcao_full 3
    !% Do a LCAO calculation before the SCF cycle and use the LCAO wavefunctions to build a new
    !% guess density and a new KS potential.
    !% Using the LCAO density as a new guess density may improve the convergence, but can
    !% also slow it down or yield wrong results (especially for spin-polarized calculations).
    !%End
    call parse_integer(datasets_check('LCAOStart'), lcao_start_default, lcao_start)
    if(.not.varinfo_valid_option('LCAOStart', lcao_start)) call input_error('LCAOStart')
    call messages_print_var_option(stdout, 'LCAOStart', lcao_start)

    lcao_done = .false.
    if (lcao_start /= LCAO_START_NONE) then
      call lcao_init(lcao, sys%gr, sys%geo, sys%st)

      ! after initialized, can check that LCAO is possible
      if(lcao_is_available(lcao)) then
        lcao_done = .true.
        
        if(present(st_start)) then
          call lcao_wf(lcao, sys%st, sys%gr, sys%geo, hm, start=st_start)
        else
          call lcao_wf(lcao, sys%st, sys%gr, sys%geo, hm)
        endif

        if (.not. present(st_start)) then
          !Just populate again the states, so that the eigenvalues are properly written
          call states_fermi(sys%st, sys%gr%mesh)
          call states_write_eigenvalues(stdout, sys%st%nst, sys%st, sys%gr%sb)

          ! Update the density and the Hamiltonian
          if (lcao_start == LCAO_START_FULL) call system_h_setup(sys, hm)
        endif
      endif

      call lcao_end(lcao)

    endif

    if(.not. lcao_done .and. .not. present(st_start)) then

      ! FIXME: the following initialization is wrong when not all
      ! wavefunctions are calculated by the Lippmann-Schwinger
      ! equation.
      ! Use free states as initial wavefunctions.
      if(sys%gr%ob_grid%open_boundaries) then
        ASSERT(sys%st%ob_nst .eq. sys%st%nst)
        ASSERT(sys%st%ob_d%nik .eq. sys%st%d%nik)
        s1 = sys%st%st_start
        s2 = sys%st%st_end
        k1 = sys%st%d%kpt%start
        k2 = sys%st%d%kpt%end
        ! the following copying does NOT ALWAYS work, especially for large numbers of k2
        !sys%st%zpsi(1:sys%gr%mesh%np, :, s1:s2, k1:k2) = sys%st%zphi(1:sys%gr%mesh%np, :, s1:s2, k1:k2)
        ! so do it the stupid and slow way
        forall (ik = k1:k2, is = s1:s2, idim = 1:sys%st%d%dim, ip = 1:sys%gr%mesh%np)
          sys%st%zpsi(ip, idim, is, ik) = sys%st%zphi(ip, idim, is, ik)
        end forall
      else
        ! Randomly generate the initial wavefunctions.
        call states_generate_random(sys%st, sys%gr%mesh)
        message(1) = "Orthogonalizing random wavefunctions."
        call write_info(1)
        call states_orthogonalize(sys%st, sys%gr%mesh)
        call v_ks_calc(sys%ks, hm, sys%st, calc_eigenval=.true.) ! get potentials
        call states_fermi(sys%st, sys%gr%mesh)                           ! occupations
      end if

    end if

    POP_SUB(lcao_run)
  end subroutine lcao_run

  ! ---------------------------------------------------------
  subroutine lcao_end(this)
    type(lcao_t), intent(inout) :: this

    PUSH_SUB(lcao_end)

    SAFE_DEALLOCATE_P(this%calc_atom)
    SAFE_DEALLOCATE_P(this%norb_atom)
    SAFE_DEALLOCATE_P(this%basis_atom)
    SAFE_DEALLOCATE_P(this%basis_orb)
    SAFE_DEALLOCATE_P(this%atom_orb_basis)
    SAFE_DEALLOCATE_P(this%radius)
    SAFE_DEALLOCATE_P(this%atom)
    SAFE_DEALLOCATE_P(this%level)
    SAFE_DEALLOCATE_P(this%ddim)

    this%initialized = .false.
    POP_SUB(lcao_end)
  end subroutine lcao_end


  ! ---------------------------------------------------------
  subroutine lcao_wf(this, st, gr, geo, hm, start)
    type(lcao_t),        intent(inout) :: this
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(in)    :: geo
    type(hamiltonian_t), intent(in)    :: hm
    integer, optional,   intent(in)    :: start

    integer :: start_
    type(profile_t), save :: prof

    ASSERT(this%initialized)

    call profiling_in(prof, "LCAO")
    PUSH_SUB(lcao_wf)

    start_ = 1
    if(present(start)) start_ = start

    if(this%alternative) then
      if (states_are_real(st)) then
        call dlcao_wf2(this, st, gr, geo, hm, start_)
      else
        call zlcao_wf2(this, st, gr, geo, hm, start_)
      end if
    else
      if (states_are_real(st)) then
        call dlcao_wf(this, st, gr, geo, hm, start_)
      else
        call zlcao_wf(this, st, gr, geo, hm, start_)
      end if
    end if
    POP_SUB(lcao_wf)
    call profiling_out(prof)
  end subroutine lcao_wf


  ! ---------------------------------------------------------
  logical function lcao_is_available(this) result(available)
    type(lcao_t), intent(in) :: this

    PUSH_SUB(lcao_is_available)
    available = this%initialized

    POP_SUB(lcao_is_available)
  end function lcao_is_available


  ! ---------------------------------------------------------
  integer function lcao_num_orbitals(this) result(norbs)
    type(lcao_t), intent(in) :: this

    PUSH_SUB(lcao_num_orbitals)
    norbs = this%norbs

    POP_SUB(lcao_num_orbitals)
  end function lcao_num_orbitals

  ! ---------------------------------------------------------

  subroutine lcao_local_index(this, ig, jg, il, jl, prow, pcol)
    type(lcao_t), intent(in)  :: this
    integer,      intent(in)  :: ig
    integer,      intent(in)  :: jg
    integer,      intent(out) :: il
    integer,      intent(out) :: jl
    integer,      intent(out) :: prow
    integer,      intent(out) :: pcol
    
    ! no PUSH_SUB, called too often
#ifdef HAVE_SCALAPACK
    call infog2l(ig, jg, this%desc(1), this%nproc(1), this%nproc(2), this%myroc(1), this%myroc(2), &
      il, jl, prow, pcol)
#else
    il = ig
    jl = jg
    prow = 0
    pcol = 0
#endif

  end subroutine lcao_local_index

  ! --------------------------------------------------------- 
  
  !> This function generates the set of an atomic orbitals for an atom
  !! and stores it in the batch orbitalb. It can be called when the
  !! orbitals are already store. In that case it does not do anything.
  subroutine lcao_get_orbital(orbitalb, sphere, st, geo, ispin, iatom, norbs)
    type(batch_t),     intent(inout) :: orbitalb
    type(submesh_t),   intent(in)    :: sphere
    type(states_t),    intent(in)    :: st
    type(geometry_t),  intent(in)    :: geo
    integer,           intent(in)    :: ispin
    integer,           intent(in)    :: iatom
    integer,           intent(in)    :: norbs

    integer :: iorb

    PUSH_SUB(lcao_get_orbital)

    if(.not. batch_is_ok(orbitalb)) then

      call profiling_in(prof_orbitals, "LCAO_ORBITALS")

      ! allocate memory
      call dbatch_new(orbitalb, 1, norbs, sphere%ns)
      
      ! generate the orbitals
      do iorb = 1, norbs
        if(iorb > species_niwfs(geo%atom(iatom)%spec)) then
          call species_get_orbital_submesh(geo%atom(iatom)%spec, sphere, iorb - species_niwfs(geo%atom(iatom)%spec), &
            st%d%dim, ispin, geo%atom(iatom)%x, orbitalb%states(iorb)%dpsi(:, 1), derivative = .true.)
        else
          call species_get_orbital_submesh(geo%atom(iatom)%spec, sphere, iorb, &
            st%d%dim, ispin, geo%atom(iatom)%x, orbitalb%states(iorb)%dpsi(:, 1))
        end if
      end do
 
      call profiling_out(prof_orbitals)
    end if

    POP_SUB(lcao_get_orbital)

  end subroutine lcao_get_orbital

  ! ---------------------------------------------------------

  !> This function deallocates a set of an atomic orbitals for an
  !! atom. It can be called when the batch is empty, in that case it
  !! does not do anything.
  subroutine lcao_end_orbital(orbitalb)
    type(batch_t),   intent(inout) :: orbitalb

    PUSH_SUB(lcao_end_orbital)

    if(batch_is_ok(orbitalb)) then
      call dbatch_delete(orbitalb)
    end if

    POP_SUB(lcao_end_orbital)

  end subroutine lcao_end_orbital

  ! ---------------------------------------------------------

#include "undef.F90"
#include "real.F90"
#include "lcao_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lcao_inc.F90"


end module lcao_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
