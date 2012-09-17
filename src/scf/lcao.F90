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
  use atom_m
  use batch_m
  use blacs_proc_grid_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use output_m
  use lalg_adv_m
  use lalg_basic_m
  use lapack_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2
  use parser_m
  use profiling_m
  use ps_m
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
    integer           :: mode
    logical           :: initialized !< are k, s and v1 matrices filled?
    integer           :: norbs   !< number of orbitals used
    integer           :: maxorbs !< largest number of orbitals that could be used
    integer, pointer  :: atom(:)
    integer, pointer  :: level(:)
    integer, pointer  :: ddim(:)
    logical           :: alternative
    logical           :: derivative
    integer, pointer  :: cst(:, :)
    integer, pointer  :: ck(:, :)
    real(4), pointer  :: dbuff(:, :, :, :)
    real(8), pointer  :: zbuff(:, :, :, :)

    !> For the alternative LCAO
    logical             :: keep_orb     !< Whether we keep orbitals in memory.
    FLOAT,   pointer    :: radius(:)    !< The localization radius of each atom orbitals
    FLOAT               :: lapdist      !< This is the extra distance that the Laplacian adds to the localization radius.
    integer             :: mult         !< The number of basis per atomic function (with derivatives is 2, 1 otherwise).
    integer             :: maxorb       !< The maximum value of the orbitals over all atoms.
    integer             :: nbasis       !< The total number of basis functions.
    !> The following functions map between a basis index and atom/orbital index
    integer, pointer    :: basis_atom(:) !< The atom that corresponds to a certain basis index
    integer, pointer    :: basis_orb(:)  !< The orbital that corresponds to a certain basis index
    integer, pointer    :: atom_orb_basis(:, :) !< The basis index that corresponds to a certain atom and orbital
    integer, pointer    :: norb_atom(:)  !< The number of orbitals per atom including mult.
    logical             :: parallel      !< Whether the LCAO is done in parallel
    integer             :: lsize(1:2)
    integer             :: nproc(1:2)
    integer             :: myroc(1:2)
    integer             :: desc(1:BLACS_DLEN)
    logical, pointer    :: calc_atom(:)
    FLOAT               :: diag_tol
    type(submesh_t), pointer :: sphere(:)
    type(batch_t),   pointer :: orbitals(:)
  end type lcao_t
  
  type(profile_t), save :: prof_orbitals

  integer, parameter :: INITRHO_PARAMAGNETIC  = 1, &
                        INITRHO_FERROMAGNETIC = 2, &
                        INITRHO_RANDOM        = 3, &
                        INITRHO_USERDEF       = 77

contains

  ! ---------------------------------------------------------
  subroutine lcao_init(this, gr, geo, st)
    type(lcao_t),         intent(out)   :: this
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(in)    :: st

    integer :: ia, n, ii, jj, maxj, idim
    integer :: mode_default

    PUSH_SUB(lcao_init)

    ! nullify everything so we can check for associated pointers when deallocating
    nullify(this%atom)
    nullify(this%level)
    nullify(this%ddim)
    nullify(this%cst)
    nullify(this%ck)
    nullify(this%dbuff)
    nullify(this%zbuff)

    nullify(this%radius)
    nullify(this%basis_atom)
    nullify(this%basis_orb)
    nullify(this%atom_orb_basis)
    nullify(this%norb_atom)
    nullify(this%calc_atom)
    nullify(this%sphere)
    nullify(this%orbitals)

    this%initialized = .true.

    ! The initial LCAO calculation is done by default if we have pseudopotentials.
    ! Otherwise, it is not the default value and has to be enforced in the input file.
    mode_default = LCAO_START_FULL
    if(geo%only_user_def) mode_default = LCAO_START_NONE
    
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
    call parse_integer(datasets_check('LCAOStart'), mode_default, this%mode)
    if(.not.varinfo_valid_option('LCAOStart', this%mode)) call input_error('LCAOStart')

    call messages_print_var_option(stdout, 'LCAOStart', this%mode)

    if(this%mode == LCAO_START_NONE) then
      POP_SUB(lcao_init)
      return
    end if
    
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
        call messages_warning(2)
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
        ! n was made too small
        this%norbs = st%nst
      else if(n > st%nst .and. n <= this%maxorbs) then
        ! n is a reasonable value
        this%norbs = n
      else if(n == 0) then
        ! using the default
        this%norbs = min(this%maxorbs, 2*st%nst)
      else
        ! n was negative, or greater than maxorbs
        this%norbs = this%maxorbs
      end if

      ASSERT(this%norbs >= st%nst)
      ASSERT(this%norbs <= this%maxorbs)
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

      call messages_write('Info: Using LCAO alternative implementation.')
      call messages_info()

      call messages_experimental('LCAO alternative implementation')

      !%Variable LCAOKeepOrbitals
      !%Type logical
      !%Default yes
      !%Section SCF::LCAO
      !%Description
      !% If set to yes (the default) Octopus keeps atomic orbitals in
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

      ! Initialize the mapping between indices

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
    integer, optional,   intent(in)    :: st_start !< use for unoccupied-states run

    type(lcao_t) :: lcao
    integer :: s1, s2, k1, k2, is, ik, ip, idim
    logical :: lcao_done
    type(profile_t), save :: prof

    PUSH_SUB(lcao_run)

    if(sys%ks%theory_level == CLASSICAL) then
      POP_SUB(lcao_run) 
      return
    end if

    if (present(st_start)) then
      ! If we are doing unocc calculation, do not mess with the correct eigenvalues
      ! of the occupied states.
      call v_ks_calc(sys%ks, hm, sys%st, sys%geo, calc_eigenval=.not. present(st_start))

      ASSERT(st_start <= sys%st%nst)
      if(st_start .gt. sys%st%nst) then ! nothing to be done in LCAO
        POP_SUB(lcao_run)
        return
      endif
    endif

    call profiling_in(prof, 'LCAO_RUN')

    call lcao_init(lcao, sys%gr, sys%geo, sys%st)

    call lcao_init_orbitals(lcao, sys%st, sys%gr, sys%geo, hm, start = st_start)

    if (.not. present(st_start)) then
      call lcao_guess_density(lcao, sys%st, sys%gr, sys%gr%sb, sys%geo, sys%st%qtot, sys%st%d%nspin, &
        sys%st%d%spin_channels, sys%st%rho)
      
      ! set up Hamiltonian (we do not call system_h_setup here because we do not want to
      ! overwrite the guess density)
      message(1) = 'Info: Setting up Hamiltonian.'
      call messages_info(1)
    
      ! get the effective potential (we don`t need the eigenvalues yet)
      call v_ks_calc(sys%ks, hm, sys%st, sys%geo, calc_eigenval=.false., calc_berry=.false.)
      ! eigenvalues have nevertheless to be initialized to something
      sys%st%eigenval = M_ZERO
      
    end if


    lcao_done = .false.

    ! after initialized, can check that LCAO is possible
    if(lcao_is_available(lcao)) then
      lcao_done = .true.
      
      call lcao_wf(lcao, sys%st, sys%gr, sys%geo, hm, start = st_start)
      
      if (.not. present(st_start)) then
        call states_fermi(sys%st, sys%gr%mesh)
        call states_write_eigenvalues(stdout, sys%st%nst, sys%st, sys%gr%sb)
        
        ! Update the density and the Hamiltonian
        if (lcao%mode == LCAO_START_FULL) call system_h_setup(sys, hm, calc_eigenval = .false.)
      endif
    else
      if(.not. present(st_start)) call init_states(sys%st, sys%gr%mesh, sys%geo)
    end if
    
    if(.not. lcao_done) then

      if(.not. sys%gr%ob_grid%open_boundaries) then

        ! Randomly generate the initial wavefunctions.
        call states_generate_random(sys%st, sys%gr%mesh, ist_start_ = st_start)

        call messages_write('Orthogonalizing random wavefunctions.')
        call messages_info()

        call states_orthogonalize(sys%st, sys%gr%mesh)
        ! If we are doing unocc calculation, do not mess with the correct eigenvalues and occupations
        ! of the occupied states.
        call v_ks_calc(sys%ks, hm, sys%st, sys%geo, calc_eigenval=.not. present(st_start)) ! get potentials
        if(.not. present(st_start)) call states_fermi(sys%st, sys%gr%mesh) ! occupations

      else
  
        ! FIXME: the following initialization is wrong when not all
        ! wavefunctions are calculated by the Lippmann-Schwinger
        ! equation.
        ! Use free states as initial wavefunctions.

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
      end if

    end if

    call lcao_end(lcao)

    call profiling_out(prof)
    POP_SUB(lcao_run)

  contains

    subroutine init_states(st, mesh, geo)
      type(states_t),   intent(inout) :: st
      type(mesh_t),     intent(inout) :: mesh
      type(geometry_t), intent(in)    :: geo

      integer :: iatom, iorb, maxorbs, ist, idim, iqn, ispin
      FLOAT, allocatable :: aorbital(:)

      PUSH_SUB(lcao_run.init_states)

      SAFE_ALLOCATE(aorbital(1:mesh%np))

      maxorbs = 0
      do iatom = 1, geo%natoms
        maxorbs = max(maxorbs, species_niwfs(geo%atom(iatom)%spec))
      end do

      ist = 0
      do iorb = 1, maxorbs
        do iatom = 1, geo%natoms
          if(iorb > species_niwfs(geo%atom(iatom)%spec)) cycle

          INCR(ist, 1)
          if(ist < st%st_start .or. ist > st%st_end) cycle

          do ispin = 1, st%d%spin_channels ! we have to iterate over spinor dimensions or spin
            idim = min(st%d%dim, ispin)

            call species_get_orbital(geo%atom(iatom)%spec, mesh, iorb, 1, geo%atom(iatom)%x, aorbital)

            do iqn = st%d%kpt%start, st%d%kpt%end
              if(st%d%ispin == SPIN_POLARIZED .and. ispin /= states_dim_get_spin_index(st%d, iqn)) cycle

              call states_set_state(st, mesh, idim, ist, iqn, aorbital)

            end do
          end do

        end do
      end do

      SAFE_DEALLOCATE_A(aorbital)
      POP_SUB(lcao_run.init_states)
    end subroutine init_states

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
    SAFE_DEALLOCATE_P(this%sphere)
    SAFE_DEALLOCATE_P(this%orbitals)

    SAFE_DEALLOCATE_P(this%atom)
    SAFE_DEALLOCATE_P(this%level)
    SAFE_DEALLOCATE_P(this%ddim)
    SAFE_DEALLOCATE_P(this%cst)
    SAFE_DEALLOCATE_P(this%ck)
    SAFE_DEALLOCATE_P(this%dbuff)
    SAFE_DEALLOCATE_P(this%zbuff)

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
        call dlcao_alt_wf(this, st, gr, geo, hm, start_)
      else
        call zlcao_alt_wf(this, st, gr, geo, hm, start_)
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

    available = this%initialized .and. this%mode /= LCAO_START_NONE

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
  !! orbitals are already stored. In that case it does not do anything.
  subroutine lcao_alt_get_orbital(orbitalb, sphere, geo, ispin, iatom, norbs)
    type(batch_t),     intent(inout) :: orbitalb
    type(submesh_t),   intent(in)    :: sphere
    type(geometry_t),  intent(in)    :: geo
    integer,           intent(in)    :: ispin
    integer,           intent(in)    :: iatom
    integer,           intent(in)    :: norbs

    integer :: iorb

    PUSH_SUB(lcao_alt_get_orbital)

    if(.not. batch_is_ok(orbitalb)) then

      call profiling_in(prof_orbitals, "LCAO_ORBITALS")

      ! allocate memory
      call dbatch_new(orbitalb, 1, norbs, sphere%np)
      
      ! generate the orbitals
      do iorb = 1, norbs
        if(iorb > species_niwfs(geo%atom(iatom)%spec)) then
          call species_get_orbital_submesh(geo%atom(iatom)%spec, sphere, iorb - species_niwfs(geo%atom(iatom)%spec), &
            ispin, geo%atom(iatom)%x, orbitalb%states(iorb)%dpsi(:, 1), derivative = .true.)
        else
          call species_get_orbital_submesh(geo%atom(iatom)%spec, sphere, iorb, &
            ispin, geo%atom(iatom)%x, orbitalb%states(iorb)%dpsi(:, 1))
        end if
      end do
 
      call profiling_out(prof_orbitals)
    end if

    POP_SUB(lcao_alt_get_orbital)

  end subroutine lcao_alt_get_orbital

  ! ---------------------------------------------------------

  !> This function deallocates a set of an atomic orbitals for an
  !! atom. It can be called when the batch is empty, in that case it
  !! does not do anything.
  subroutine lcao_alt_end_orbital(orbitalb)
    type(batch_t),   intent(inout) :: orbitalb

    PUSH_SUB(lcao_alt_end_orbital)

    if(batch_is_ok(orbitalb)) then
      call batch_delete(orbitalb)
    end if

    POP_SUB(lcao_alt_end_orbital)

  end subroutine lcao_alt_end_orbital

  ! ---------------------------------------------------------

  subroutine lcao_atom_density(this, st, gr, sb, geo, iatom, spin_channels, rho)
    type(lcao_t),      intent(inout) :: this
    type(states_t),    intent(in)    :: st
    type(grid_t),      intent(in)    :: gr
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(in)    :: geo
    integer,           intent(in)    :: iatom
    integer,           intent(in)    :: spin_channels
    FLOAT,             intent(inout) :: rho(:, :) !< (gr[%fine]%mesh%np, spin_channels)
    
    FLOAT, allocatable :: dorbital(:, :)
    CMPLX, allocatable :: zorbital(:, :)
    FLOAT, allocatable :: factors(:)
    FLOAT :: factor, aa
    integer :: iorb, ip, ii, ll, mm
    type(ps_t), pointer :: ps
    logical :: use_stored_orbitals

    PUSH_SUB(lcao_atom_density)

    rho = M_ZERO

    use_stored_orbitals = species_is_ps(geo%atom(iatom)%spec) &
      .and. states_are_real(st) .and. spin_channels == 1 .and. lcao_is_available(this) &
      .and. st%d%dim == 1 .and. .not. gr%have_fine_mesh

    ps => species_ps(geo%atom(iatom)%spec)

    ! we can use the orbitals we already have calculated
    if(use_stored_orbitals) then
      ASSERT(.not. gr%have_fine_mesh)

      if(.not. this%alternative) then
        
        if(states_are_real(st)) then
          SAFE_ALLOCATE(dorbital(1:gr%mesh%np, 1:st%d%dim))
        else
          SAFE_ALLOCATE(zorbital(1:gr%mesh%np, 1:st%d%dim))
        end if
        
        do iorb = 1, this%norbs
          if(iatom /= this%atom(iorb)) cycle
          
          call species_iwf_ilm(geo%atom(iatom)%spec, this%level(iorb), 1, ii, ll, mm)
          factor = ps%conf%occ(ii, 1)/(CNST(2.0)*ll + CNST(1.0))
         
          if(states_are_real(st)) then
            call dget_ao(this, st, gr%mesh, geo, iorb, 1, dorbital, use_psi = .true.)
            !$omp parallel do
            do ip = 1, gr%mesh%np
              rho(ip, 1) = rho(ip, 1) + factor*dorbital(ip, 1)**2
            end do
          else
            call zget_ao(this, st, gr%mesh, geo, iorb, 1, zorbital, use_psi = .true.)
            !$omp parallel do
            do ip = 1, gr%mesh%np
              rho(ip, 1) = rho(ip, 1) + factor*abs(zorbital(ip, 1))**2
            end do
          end if
          
        end do

        SAFE_DEALLOCATE_A(dorbital)
        SAFE_DEALLOCATE_A(zorbital)

      else

        call lcao_alt_get_orbital(this%orbitals(iatom), this%sphere(iatom), geo, 1, iatom, this%norb_atom(iatom))

        SAFE_ALLOCATE(factors(1:this%norb_atom(iatom)))

        do iorb = 1, this%norb_atom(iatom)
          call species_iwf_ilm(geo%atom(iatom)%spec, iorb, 1, ii, ll, mm)
          factors(iorb) = ps%conf%occ(ii, 1)/(CNST(2.0)*ll + CNST(1.0))
        end do

        !$omp parallel do private(ip, aa)
        do ip = 1, this%sphere(iatom)%np
          aa = CNST(0.0)
          do iorb = 1, this%norb_atom(iatom)
            aa = aa + factors(iorb)*this%orbitals(iatom)%states_linear(iorb)%dpsi(ip)**2
          end do
          rho(this%sphere(iatom)%map(ip), 1) = aa
        end do

        SAFE_DEALLOCATE_A(factors)

      end if

    else
      call species_atom_density(gr%fine%mesh, sb, geo%atom(iatom), spin_channels, rho)
    end if

    POP_SUB(lcao_atom_density)
  end subroutine lcao_atom_density

  ! ---------------------------------------------------------
  !> builds a density which is the sum of the atomic densities
  subroutine lcao_guess_density(this, st, gr, sb, geo, qtot, nspin, spin_channels, rho)
    type(lcao_t),      intent(inout) :: this
    type(states_t),    intent(in)    :: st
    type(grid_t),      intent(in)    :: gr
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(in)    :: geo
    FLOAT,             intent(in)    :: qtot  !< the total charge of the system
    integer,           intent(in)    :: nspin, spin_channels
    FLOAT,             intent(out)   :: rho(:, :)

    integer :: ia, is, idir, gmd_opt
    integer, save :: iseed = 321
    type(block_t) :: blk
    FLOAT :: rr, rnd, phi, theta, mag(1:3), lmag, n1, n2
    FLOAT, allocatable :: atom_rho(:,:)
    logical :: parallelized_in_atoms


    PUSH_SUB(lcao_guess_density)

    parallelized_in_atoms = .false.

    if (spin_channels == 1) then
      gmd_opt = INITRHO_PARAMAGNETIC
    else
      !%Variable GuessMagnetDensity
      !%Type integer
      !%Default ferromagnetic
      !%Section SCF
      !%Description
      !% The guess density for the SCF cycle is just the sum of all the atomic densities.
      !% When performing spin-polarized or non-collinear-spin calculations this option sets 
      !% the guess magnetization density.
      !%
      !% For anti-ferromagnetic configurations, the <tt>user_defined</tt> option should be used.
      !%
      !% Note that if the <tt>paramagnetic</tt> option is used, the final ground state will also be
      !% paramagnetic, but the same is not true for the other options.
      !%Option paramagnetic 1
      !% Magnetization density is zero.
      !%Option ferromagnetic 2
      !% Magnetization density is the sum of the atomic magnetization densities.
      !%Option random 3
      !% Each atomic magnetization density is randomly rotated.
      !%Option user_defined 77
      !% The atomic magnetization densities are rotated so that the magnetization 
      !% vector has the same direction as a vector provided by the user. In this case,
      !% the <tt>AtomsMagnetDirection</tt> block has to be set.
      !%End
      call parse_integer(datasets_check('GuessMagnetDensity'), INITRHO_FERROMAGNETIC, gmd_opt)
      if(.not.varinfo_valid_option('GuessMagnetDensity', gmd_opt)) call input_error('GuessMagnetDensity')
      call messages_print_var_option(stdout, 'GuessMagnetDensity', gmd_opt)
    end if

    rho = M_ZERO
    select case (gmd_opt)
    case (INITRHO_PARAMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:gr%fine%mesh%np, 1:1))

      parallelized_in_atoms = .true.

      do ia = geo%atoms_dist%start, geo%atoms_dist%end
        call lcao_atom_density(this, st, gr, sb, geo, ia, 1, atom_rho)
        rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + atom_rho(1:gr%fine%mesh%np, 1)
      end do

      if (spin_channels == 2) then
        rho(1:gr%fine%mesh%np, 1) = M_HALF*rho(1:gr%fine%mesh%np, 1)
        rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 1)
      end if

    case (INITRHO_FERROMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:gr%fine%mesh%np, 1:2))

      parallelized_in_atoms = .true.

      atom_rho = M_ZERO
      rho = M_ZERO
      do ia = geo%atoms_dist%start, geo%atoms_dist%end
        call lcao_atom_density(this, st, gr, sb, geo, ia, 2, atom_rho(1:gr%fine%mesh%np, 1:2))
        rho(1:gr%fine%mesh%np, 1:2) = rho(1:gr%fine%mesh%np, 1:2) + atom_rho(1:gr%fine%mesh%np, 1:2)
      end do

    case (INITRHO_RANDOM) ! Randomly oriented spins
      SAFE_ALLOCATE(atom_rho(1:gr%fine%mesh%np, 1:2))
      do ia = 1, geo%natoms
        call lcao_atom_density(this, st, gr, sb, geo, ia, 2, atom_rho)

        if (nspin == 2) then
          call quickrnd(iseed, rnd)
          rnd = rnd - M_HALF
          if (rnd > M_ZERO) then
            rho(1:gr%fine%mesh%np, 1:2) = rho(1:gr%fine%mesh%np, 1:2) + atom_rho(1:gr%fine%mesh%np, 1:2)
          else
            rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + atom_rho(1:gr%fine%mesh%np, 2)
            rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 2) + atom_rho(1:gr%fine%mesh%np, 1)
          end if
        elseif (nspin == 4) then
          call quickrnd(iseed, phi)
          call quickrnd(iseed, theta)
          phi = phi*M_TWO*M_PI
          theta = theta*M_PI
          rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + cos(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 1) &
            + sin(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 2)
          rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 2) + sin(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 1) &
            + cos(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 2)
          rho(1:gr%fine%mesh%np, 3) = rho(1:gr%fine%mesh%np, 3) + cos(theta/M_TWO)*sin(theta/M_TWO)*cos(phi)* &
            (atom_rho(1:gr%fine%mesh%np, 1) - atom_rho(1:gr%fine%mesh%np, 2))
          rho(1:gr%fine%mesh%np, 4) = rho(1:gr%fine%mesh%np, 4) - cos(theta/M_TWO)*sin(theta/M_TWO)*sin(phi)* &
            (atom_rho(1:gr%fine%mesh%np, 1) - atom_rho(1:gr%fine%mesh%np, 2))
        end if
      end do

    case (INITRHO_USERDEF) ! User-defined
      
      !%Variable AtomsMagnetDirection
      !%Type block
      !%Section Hamiltonian
      !%Description
      !% This option is only used when <tt>GuessMagnetDensity</tt> is
      !% set to <tt>user_defined</tt>. It provides a direction for the
      !% magnetization vector of each atom when building the guess
      !% density. In order to do that, the user should specify the
      !% coordinates of a vector that has the desired direction and
      !% norm.  Note that it is necessary to maintain the ordering in
      !% which the species were defined in the coordinates
      !% specifications.
      !%
      !% For spin-polarized calculations, the vectors should have only
      !% one component; for non-collinear-spin calculations, they
      !% should have three components.
      !%End
      if(parse_block(datasets_check('AtomsMagnetDirection'), blk) < 0) then
        message(1) = "AtomsMagnetDirection block is not defined."
        call messages_fatal(1)
      end if

      if (parse_block_n(blk) /= geo%natoms) then
        message(1) = "AtomsMagnetDirection block has the wrong number of rows."
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(atom_rho(1:gr%fine%mesh%np, 1:2))
      do ia = 1, geo%natoms
        !Read from AtomsMagnetDirection block 
        if (nspin == 2) then
          call parse_block_float(blk, ia-1, 0, mag(1))
          lmag = abs(mag(1))
        elseif (nspin == 4) then
          do idir = 1, 3
            call parse_block_float(blk, ia-1, idir-1, mag(idir))
            if (abs(mag(idir)) < CNST(1.0e-20)) mag(idir) = M_ZERO
          end do
          lmag = sqrt(dot_product(mag(1:3), mag(1:3)))
        end if

        !Get atomic density
        call lcao_atom_density(this, st, gr, sb, geo, ia, 2, atom_rho)

        !Scale magnetization density
        n1 = dmf_integrate(gr%fine%mesh, atom_rho(:, 1))
        n2 = dmf_integrate(gr%fine%mesh, atom_rho(:, 2))
        if (lmag > n1 + n2) then
          mag = mag*(n1 + n2)/lmag
          lmag = n1 + n2
        elseif (lmag == M_ZERO) then
          if (n1 - n2 == M_ZERO) then
            rho(1:gr%fine%mesh%np, 1:2) = rho(1:gr%fine%mesh%np, 1:2) + atom_rho(1:gr%fine%mesh%np, 1:2)
          else
            atom_rho(:, 1) = (atom_rho(:, 1) + atom_rho(:, 2))/M_TWO
            rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + atom_rho(1:gr%fine%mesh%np, 1)
            rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 2) + atom_rho(1:gr%fine%mesh%np, 1)
          end if
          cycle
        end if
        if (n1 - n2 /= lmag .and. n2 /= M_ZERO) then
          if (n1 - n2 < lmag) then
            atom_rho(:, 1) = atom_rho(:, 1) + (lmag - n1 + n2)/M_TWO/n2*atom_rho(:, 2)
            atom_rho(:, 2) = (n1 + n2 - lmag)/M_TWO/n2*atom_rho(:, 2)
          elseif (n1 - n2 > lmag) then
            atom_rho(:, 2) = atom_rho(:, 2) + (n1 - n2 - lmag)/M_TWO/n1*atom_rho(:, 1)
            atom_rho(:, 1) = (lmag + n1 + n2)/M_TWO/n1*atom_rho(:, 1)
          end if
        end if

        !Rotate magnetization density
        if (nspin == 2) then
          if (mag(1) > M_ZERO) then
            rho(1:gr%fine%mesh%np, 1:2) = rho(1:gr%fine%mesh%np, 1:2) + atom_rho(1:gr%fine%mesh%np, 1:2)
          else
            rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + atom_rho(1:gr%fine%mesh%np, 2)
            rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 2) + atom_rho(1:gr%fine%mesh%np, 1)
          end if

        elseif (nspin == 4) then
          theta = acos(mag(3)/lmag)
          if (mag(1) == M_ZERO) then
            if (mag(2) == M_ZERO) then
              phi = M_ZERO
            elseif (mag(2) < M_ZERO) then
              phi = M_PI*M_TWOTHIRD
            elseif (mag(2) > M_ZERO) then
              phi = M_PI*M_HALF
            end if
          else
            if (mag(2) < M_ZERO) then
              phi = M_TWO*M_PI - acos(mag(1)/sin(theta)/lmag)
            elseif (mag(2) >= M_ZERO) then
              phi = acos(mag(1)/sin(theta)/lmag)
            end if
          end if

          rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + cos(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 1) &
            + sin(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 2)
          rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 2) + sin(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 1) &
               + cos(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 2)
          rho(1:gr%fine%mesh%np, 3) = rho(1:gr%fine%mesh%np, 3) + cos(theta/M_TWO)*sin(theta/M_TWO)*cos(phi)* &
            (atom_rho(1:gr%fine%mesh%np, 1) - atom_rho(1:gr%fine%mesh%np, 2))
          rho(1:gr%fine%mesh%np, 4) = rho(1:gr%fine%mesh%np, 4) - cos(theta/M_TWO)*sin(theta/M_TWO)*sin(phi)* &
            (atom_rho(1:gr%fine%mesh%np, 1) - atom_rho(1:gr%fine%mesh%np, 2))
        end if
      end do

      call parse_block_end(blk)

    end select


#ifdef HAVE_MPI
    if(geo%atoms_dist%parallel .and. parallelized_in_atoms) then
      do is = 1, spin_channels
        atom_rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, is)
        call MPI_Allreduce(atom_rho(1, 1), rho(1, is), gr%fine%mesh%np, &
          MPI_FLOAT, MPI_SUM, geo%atoms_dist%mpi_grp%comm, mpi_err)
      end do
    end if
#endif

    ! we now renormalize the density (necessary if we have a charged system)
    rr = M_ZERO
    do is = 1, spin_channels
      rr = rr + dmf_integrate(gr%fine%mesh, rho(:, is))
    end do

    write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', rr
    call messages_info(1)

    rr = qtot / rr
    rho = rr * rho
    rr = M_ZERO
    do is = 1, spin_channels
      rr = rr + dmf_integrate(gr%fine%mesh, rho(:, is))
    end do

    write(message(1),'(a,f13.6)')'Info: Renormalized total charge = ', rr
    call messages_info(1)

    SAFE_DEALLOCATE_A(atom_rho)
    POP_SUB(lcao_guess_density)
  end subroutine lcao_guess_density

  ! ---------------------------------------------------------

  subroutine lcao_init_orbitals(this, st, gr, geo, hm, start)
    type(lcao_t),        intent(inout) :: this
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(in)    :: geo
    type(hamiltonian_t), intent(in)    :: hm
    integer, optional,   intent(in)    :: start

    if(.not. lcao_is_available(this)) return
  
    PUSH_SUB(lcao_init_orbitals)
        
    if(.not. this%alternative) then
      if(states_are_real(st)) then
        call dinit_orbitals(this, st, gr, geo, hm, start)
      else
        call zinit_orbitals(this, st, gr, geo, hm, start)
      end if
    else
      if(states_are_real(st)) then
        call dlcao_alt_init_orbitals(this, st, gr, geo, hm, start)
      else
        call zlcao_alt_init_orbitals(this, st, gr, geo, hm, start)
      end if

    end if

    POP_SUB(lcao_init_orbitals)
  end subroutine lcao_init_orbitals
  
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
