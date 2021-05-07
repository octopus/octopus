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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module lcao_oct_m
  use atom_oct_m
  use atomic_orbital_oct_m
  use batch_oct_m
  use blacs_oct_m
  use blacs_proc_grid_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use lapack_oct_m
  use loct_oct_m
  use magnetic_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use ps_oct_m
  use quickrnd_oct_m
  use simul_box_oct_m
  use scalapack_oct_m
  use smear_oct_m
  use space_oct_m
  use species_oct_m
  use species_pot_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use states_elec_io_oct_m
  use submesh_oct_m
  use symmetrizer_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use varinfo_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::             &
    lcao_t,             &
    lcao_init,          &
    lcao_init_orbitals, &
    lcao_wf,            &
    lcao_run,           &
    lcao_end,           &
    lcao_is_available,  &
    lcao_num_orbitals

  type lcao_t
    private
    integer                 :: mode
    logical                 :: complex_ylms !< whether to use real or complex Ylms
    logical                 :: initialized !< are k, s and v1 matrices filled?
    integer                 :: norbs   !< number of orbitals used
    integer                 :: maxorbs !< largest number of orbitals that could be used
    integer,    allocatable :: atom(:)
    integer,    allocatable :: level(:)
    integer,    allocatable :: ddim(:)
    logical                 :: alternative
    logical                 :: derivative
    integer,    allocatable :: cst(:, :)
    integer,    allocatable :: ck(:, :)
    real(4),    allocatable :: dbuff(:, :, :, :) !< single-precision buffer
    complex(4), allocatable :: zbuff(:, :, :, :) !< single-precision buffer
    logical                 :: initialized_orbitals
    FLOAT                   :: orbital_scale_factor
 
    !> For optimization purposes
    logical,    allocatable :: is_empty(:) !< True if the submesh contains zero point.
                                           !< This occurs in domain parallelization if the atom is not
                                           !< in the local domain

    !> For the alternative LCAO
    logical                 :: keep_orb     !< Whether we keep orbitals in memory.
    FLOAT,      allocatable :: radius(:)    !< The localization radius of each atom orbitals
    FLOAT                   :: lapdist      !< This is the extra distance that the Laplacian adds to the localization radius.
    integer                 :: mult         !< The number of basis orbitals per atomic function (with derivatives is 2, 1 otherwise).
    integer                 :: maxorb       !< The maximum value of the orbitals over all atoms.
    !> The following functions map between a basis index and atom/orbital index
    integer, allocatable    :: basis_atom(:) !< The atom that corresponds to a certain basis index
    integer, allocatable    :: basis_orb(:)  !< The orbital that corresponds to a certain basis index
    integer, allocatable    :: atom_orb_basis(:, :) !< The basis index that corresponds to a certain atom and orbital
    integer, allocatable    :: norb_atom(:)  !< The number of orbitals per atom including mult.
    logical                 :: parallel      !< Whether the LCAO is done in parallel
    integer                 :: lsize(1:2)
    integer                 :: nproc(1:2)
    integer                 :: myroc(1:2)
    integer                 :: desc(1:BLACS_DLEN)
    logical, allocatable    :: calc_atom(:)
    FLOAT                   :: diag_tol
    type(submesh_t), allocatable :: sphere(:)
    type(batch_t),   allocatable :: orbitals(:)
    logical,         allocatable :: is_orbital_initialized(:) !< array to store which orbitals are already initialized
  end type lcao_t
  
  type(profile_t), save :: prof_orbitals

  integer, parameter :: INITRHO_PARAMAGNETIC  = 1, &
                        INITRHO_FERROMAGNETIC = 2, &
                        INITRHO_RANDOM        = 3, &
                        INITRHO_USERDEF       = 77

contains

  ! ---------------------------------------------------------
  subroutine lcao_init(this, namespace, space, gr, ions, st)
    type(lcao_t),         intent(out) :: this
    type(namespace_t),    intent(in)  :: namespace
    type(space_t),        intent(in)  :: space
    type(grid_t),         intent(in)  :: gr
    type(ions_t),         intent(in)  :: ions
    type(states_elec_t),  intent(in)  :: st

    integer :: ia, n, iorb, jj, maxj, idim
    integer :: ii, ll, mm, norbs
    integer :: mode_default
    FLOAT   :: max_orb_radius, maxradius
    integer :: iunit_o

    PUSH_SUB(lcao_init)

    this%initialized = .true.

    ! initialization, in case we leave this routine before LCAOAlternative is parsed
    this%alternative = .false.

    ! The initial LCAO calculation is done by default if we have species representing atoms.
    ! Otherwise, it is not the default value and has to be enforced in the input file.
    mode_default = OPTION__LCAOSTART__LCAO_STATES
    if(ions%only_user_def) mode_default = OPTION__LCAOSTART__LCAO_NONE
    
    !%Variable LCAOStart
    !%Type integer
    !%Section SCF::LCAO
    !%Description
    !% Before starting a SCF calculation, <tt>Octopus</tt> can perform
    !% a linear combination of atomic orbitals (LCAO) calculation.
    !% These can provide <tt>Octopus</tt> with a good set
    !% of initial wavefunctions and with a new guess for the density.
    !% (Up to the current version, only a minimal basis set is used.)
    !% The default is <tt>lcao_states</tt> if at least one species representing an atom is present.
    !% The default is <tt>lcao_none</tt> if all species are <tt>species_user_defined</tt>,
    !% <tt>species_charge_density</tt>, <tt>species_from_file</tt>, or <tt>species_jellium_slab</tt>.
    !%
    !% The initial guess densities for LCAO are taken from the atomic orbitals for pseudopotential species;
    !% from the natural charge density for <tt>species_charge_density</tt>, <tt>species_point</tt>,
    !% <tt>species_jellium</tt>, and <tt>species_jellium_slab</tt>;
    !% or uniform for <tt>species_full_delta</tt>, <tt>species_full_gaussian</tt>,
    !% <tt>species_user_defined</tt>, or <tt>species_from_file</tt>.
    !% Pseudopotential species use the pseudo-wavefunctions as orbitals, full-potential atomic species
    !% (<tt>species_full_delta</tt> and <tt>species_full_gaussian</tt>) use hydrogenic wavefunctions, and
    !% others use harmonic-oscillator wavefunctions.
    !%
    !% Note: Some pseudopotential files (CPI, FHI for example) do not
    !% contain full information about the orbitals. In this case,
    !% Octopus generates the starting density from the normalized
    !% square root of the local potential. If no orbitals are
    !% available at all from the pseudopotential files, Octopus will
    !% not be able to perform an LCAO and the initial states will be
    !% randomized.
    !%
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
    call parse_variable(namespace, 'LCAOStart', mode_default, this%mode)
    if(.not.varinfo_valid_option('LCAOStart', this%mode)) call messages_input_error(namespace, 'LCAOStart')

    call messages_print_var_option(stdout, 'LCAOStart', this%mode)

    if(this%mode == OPTION__LCAOSTART__LCAO_NONE) then
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
    !% large systems and parallel in states. It is not working for spinors, however.
    !%End
    call parse_variable(namespace, 'LCAOAlternative', .false., this%alternative)
    ! DAS: For spinors, you will always get magnetization in (1, 0, 0) direction, and the
    ! eigenvalues will be incorrect. This is due to confusion between spins and spinors in the code.
    if(st%d%ispin == SPINORS .and. this%alternative) then
      message(1) = "LCAOAlternative is not working for spinors."
      call messages_fatal(1)
    end if
    if (space%is_periodic() .and. this%alternative) then
      call messages_experimental("LCAOAlternative in periodic systems")
      ! specifically, if you get the message about submesh radius > box size, results will probably be totally wrong.
    end if

    !%Variable LCAOComplexYlms
    !%Type logical
    !%Default false
    !%Section SCF::LCAO
    !%Description
    !% If set to true, and using complex states, complex spherical harmonics will be used, <i>i.e.</i>
    !% with <math>e^{\pm i m \phi}</math>.
    !% If false, real spherical harmonics with <math>\sin(m \phi)</math> or <math>\cos(m \phi)</math> are used.
    !% This variable will make it more likely to get states that are eigenvectors of the <math>L_z</math>
    !% operator, with a definite angular momentum.
    !%End

    if(states_are_complex(st)) then
      call parse_variable(namespace, 'LCAOComplexYlms', .false., this%complex_ylms)
    else
      this%complex_ylms = .false.
    end if

    if(debug%info .and. mpi_grp_is_root(mpi_world)) then
      call io_mkdir('debug/lcao', namespace)
      iunit_o = io_open('debug/lcao/orbitals', namespace, action='write')
      write(iunit_o,'(7a6)') 'iorb', 'atom', 'level', 'i', 'l', 'm', 'spin'
    end if

    if(.not. this%alternative) then

      !%Variable LCAOScaleFactor
      !%Type float
      !%Default 1.0
      !%Section SCF::LCAO
      !%Description
      !% The coordinates of the atomic orbitals used by the LCAO
      !% procedure will be rescaled by the value of this variable. 1.0 means no rescaling.
      !%End
      call parse_variable(namespace, 'LCAOScaleFactor', CNST(1.0), this%orbital_scale_factor)
      call messages_print_var_value(stdout, 'LCAOScaleFactor', this%orbital_scale_factor)

      !%Variable LCAOMaximumOrbitalRadius
      !%Type float
      !%Default 20.0 a.u.
      !%Section SCF::LCAO
      !%Description
      !% The LCAO procedure will ignore orbitals that have an
      !% extent greater that this value.
      !%End
      call parse_variable(namespace, 'LCAOMaximumOrbitalRadius', CNST(20.0), max_orb_radius, unit = units_inp%length)
      call messages_print_var_value(stdout, 'LCAOMaximumOrbitalRadius', max_orb_radius, units_out%length)

      ! count the number of orbitals available
      maxj = 0
      this%maxorbs = 0
      do ia = 1, ions%natoms
        maxj = max(maxj, species_niwfs(ions%atom(ia)%species) )
        this%maxorbs = this%maxorbs + species_niwfs(ions%atom(ia)%species)
      end do

      this%maxorbs = this%maxorbs*st%d%dim

      if(this%maxorbs == 0) then
        call messages_write('The are no atomic orbitals available, cannot do LCAO.')
        call messages_warning()
        this%mode = OPTION__LCAOSTART__LCAO_NONE
        POP_SUB(lcao_init)
        return
      end if
      
      ! generate tables to know which indices each atomic orbital has

      SAFE_ALLOCATE( this%atom(1:this%maxorbs))
      SAFE_ALLOCATE(this%level(1:this%maxorbs))
      SAFE_ALLOCATE( this%ddim(1:this%maxorbs))

      SAFE_ALLOCATE(this%is_empty(1:this%maxorbs))
      this%is_empty = .false.

      ! this is determined by the stencil we are using and the spacing
      this%lapdist = maxval(abs(gr%mesh%idx%enlarge)*gr%mesh%spacing)

      ! calculate the radius of each orbital
      SAFE_ALLOCATE(this%radius(1:ions%natoms))

      do ia = 1, ions%natoms
        norbs = species_niwfs(ions%atom(ia)%species)

        maxradius = M_ZERO
        do iorb = 1, norbs
          call species_iwf_ilm(ions%atom(ia)%species, iorb, 1, ii, ll, mm)
          maxradius = max(maxradius, species_get_iwf_radius(ions%atom(ia)%species, ii, is = 1))
        end do

        maxradius = min(maxradius, M_TWO*maxval(gr%sb%lsize(1:gr%sb%dim)))

        this%radius(ia) = maxradius
      end do


      ! Each atom provides niwfs pseudo-orbitals (this number is given in
      ! ions%atom(ia)%species%niwfs for atom number ia). This number is
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

      iorb = 1
      do jj = 1, maxj
        do ia = 1, ions%natoms
          do idim = 1,st%d%dim
            if(jj > species_niwfs(ions%atom(ia)%species) ) cycle
            call species_iwf_ilm(ions%atom(ia)%species, jj, idim, ii, ll, mm)            
            if(this%orbital_scale_factor*species_get_iwf_radius(ions%atom(ia)%species, ii, is = 1) >= max_orb_radius) cycle

            this%atom(iorb) = ia
            this%level(iorb) = jj
            this%ddim(iorb) = idim

            if(debug%info .and. mpi_grp_is_root(mpi_world)) then
              write(iunit_o,'(7i6)') iorb, this%atom(iorb), this%level(iorb), ii, ll, mm, this%ddim(iorb)
            end if

            iorb = iorb + 1
          end do
        end do
      end do

      if(debug%info .and. mpi_grp_is_root(mpi_world)) then
        call io_close(iunit_o)
      end if

      ! some orbitals might have been removed because of their radii
      if(this%maxorbs /= iorb - 1) then
        call messages_write('Info: ')
        call messages_write(this%maxorbs - iorb + 1)
        call messages_write(' of ')
        call messages_write(this%maxorbs)
        call messages_write(' orbitals cannot be used for the LCAO calculation,')
        call messages_new_line()
        call messages_write('      their radii exceeds LCAOMaximumOrbitalRadius (')
        call messages_write(max_orb_radius, units = units_out%length, fmt = '(f6.1)')
        call messages_write(').')
        call messages_warning()

        this%maxorbs = iorb - 1
      end if

      if(this%maxorbs < st%nst) then
        call messages_write('Cannot do LCAO for all states because there are not enough atomic orbitals.')
        call messages_new_line()

        call messages_write('Required: ')
        call messages_write(st%nst)
        call messages_write('. Available: ')
        call messages_write(this%maxorbs)
        call messages_write('. ')
        call messages_write(st%nst - this%maxorbs)
        call messages_write(' orbitals will be randomized.')
        call messages_warning()
      end if

      !%Variable LCAODimension
      !%Type integer
      !%Section SCF::LCAO
      !%Description
      !% (Only applies if <tt>LCAOAlternative = no</tt>.)
      !% Before starting the SCF cycle, an initial LCAO calculation can be performed
      !% in order to obtain reasonable initial guesses for spin-orbitals and densities.
      !% For this purpose, the code calculates a number of atomic orbitals.
      !% The number available for a species described by a pseudopotential is all the
      !% orbitals up the maximum angular momentum in the pseudopotential, minus any orbitals that
      !% are found to be unbound. For non-pseudopotential species, the number is equal to
      !% twice the valence charge.
      !% The default dimension for the LCAO basis
      !% set will be the sum of all these numbers, or twice the number of required orbitals
      !% for the full calculation, whichever is less.
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
      call parse_variable(namespace, 'LCAODimension', 0, n)

      if(n > 0 .and. n <= st%nst .and. st%nst <= this%maxorbs) then
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

      ASSERT(this%norbs <= this%maxorbs)

      SAFE_ALLOCATE(this%cst(1:this%norbs, 1:st%d%spin_channels))
      SAFE_ALLOCATE(this%ck(1:this%norbs, 1:st%d%spin_channels))
      this%initialized_orbitals = .false.
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
      !% Only applies if <tt>LCAOAlternative = true</tt>.
      !% If set to yes (the default) Octopus keeps atomic orbitals in
      !% memory during the LCAO procedure. If set to no, the orbitals
      !% are generated each time that they are needed, increasing
      !% computational time but saving memory.
      !%
      !% When set to yes, Octopus prints the amount of memory per node
      !% that is required to store the orbitals.
      !%
      !%End
      call parse_variable(namespace, 'LCAOKeepOrbitals', .true., this%keep_orb)

      !%Variable LCAOExtraOrbitals
      !%Type logical
      !%Default false
      !%Section SCF::LCAO
      !%Description
      !% Only applies if <tt>LCAOAlternative = true</tt>, and all species are pseudopotentials.
      !% (experimental) If this variable is set to yes, the LCAO
      !% procedure will add an extra set of numerical orbitals (by
      !% using the derivative of the radial part of the original
      !% orbitals). Note that this corresponds roughly to adding orbitals
      !% with higher principal quantum numbers, but the same angular momentum.
      !% This option may cause problems for unoccupied states since you may miss
      !% some lower-lying states which correspond to higher angular momenta instead
      !% of higher principal quantum number.
      !%End
      call parse_variable(namespace, 'LCAOExtraOrbitals', .false., this%derivative)

      ! DAS: if you calculate the Na atom this way, spin-polarized, with just one unoccupied state,
      ! you will obtain states (up and down) which are actually the 10th states if you start with
      ! random wavefunctions! We really need to implement taking the derivative of the angular part
      ! instead to be sure of getting decent results!

      if(this%derivative) then
        call messages_experimental('LCAO extra orbitals')

        if(st%nst * st%smear%el_per_state > st%qtot) then
          message(1) = "Lower-lying empty states may be missed with LCAOExtraOrbitals."
          call messages_warning(1)
        end if
      end if

      !%Variable LCAODiagTol
      !%Type float
      !%Default 1e-10
      !%Section SCF::LCAO
      !%Description
      !% Only applies if <tt>LCAOAlternative = true</tt>.
      !% The tolerance for the diagonalization of the LCAO Hamiltonian.
      !%End
      call parse_variable(namespace, 'LCAODiagTol', CNST(1e-10), this%diag_tol)

      if(this%derivative) then
        this%mult = 2
      else
        this%mult = 1
      end if

      SAFE_ALLOCATE(this%sphere(1:ions%natoms))
      SAFE_ALLOCATE(this%orbitals(1:ions%natoms))
      SAFE_ALLOCATE(this%is_orbital_initialized(1:ions%natoms))
      this%is_orbital_initialized = .false.

      SAFE_ALLOCATE(this%norb_atom(1:ions%natoms))

      this%maxorb = 0
      this%norbs = 0
      do iatom = 1, ions%natoms
        this%norb_atom(iatom) = this%mult*species_niwfs(ions%atom(iatom)%species)
        this%maxorb = max(this%maxorb, species_niwfs(ions%atom(iatom)%species))
        this%norbs = this%norbs + species_niwfs(ions%atom(iatom)%species)
      end do

      this%maxorb = this%maxorb*this%mult
      this%norbs = this%norbs*this%mult

      SAFE_ALLOCATE(this%basis_atom(1:this%norbs))
      SAFE_ALLOCATE(this%basis_orb(1:this%norbs))
      SAFE_ALLOCATE(this%atom_orb_basis(1:ions%natoms, 1:this%maxorb))

      ! Initialize the mapping between indices

      ibasis = 0
      do iatom = 1, ions%natoms
        norbs = species_niwfs(ions%atom(iatom)%species)

        do iorb = 1, this%mult*norbs
          ibasis = ibasis + 1
          this%atom_orb_basis(iatom, iorb) = ibasis
          this%basis_atom(ibasis) = iatom
          this%basis_orb(ibasis) = iorb

          ! no stored spin index in alternative mode
          if(debug%info .and. mpi_grp_is_root(mpi_world)) then
            call species_iwf_ilm(ions%atom(iatom)%species, iorb, 1, ii, ll, mm)
            write(iunit_o,'(7i6)') ibasis, iatom, iorb, ii, ll, mm, 1
          end if
        end do
      end do

      if(debug%info .and. mpi_grp_is_root(mpi_world)) then
        call io_close(iunit_o)
      end if

      ! this is determined by the stencil we are using and the spacing
      this%lapdist = maxval(abs(gr%mesh%idx%enlarge)*gr%mesh%spacing)

      ! calculate the radius of each orbital
      SAFE_ALLOCATE(this%radius(1:ions%natoms))

      do iatom = 1, ions%natoms
        norbs = species_niwfs(ions%atom(iatom)%species)

        maxradius = M_ZERO
        do iorb = 1, norbs
          call species_iwf_ilm(ions%atom(iatom)%species, iorb, 1, ii, ll, mm)
          maxradius = max(maxradius, species_get_iwf_radius(ions%atom(iatom)%species, ii, is = 1))
        end do

        if(this%derivative) maxradius = maxradius + this%lapdist

        maxradius = min(maxradius, M_TWO*maxval(gr%sb%lsize(1:gr%sb%dim)))

        this%radius(iatom) = maxradius
      end do

      SAFE_ALLOCATE(this%calc_atom(1:ions%natoms))
      this%calc_atom = .true.

      ! initialize parallel data
#ifndef HAVE_SCALAPACK
      this%parallel = .false.
#else
      this%parallel = (st%parallel_in_states .or. gr%mesh%parallel_in_domains) &
        .and. .not. blacs_proc_grid_null(st%dom_st_proc_grid)

      if(this%parallel) then      
        nbl = min(16, this%norbs)

        ! The size of the distributed matrix in each node
        this%lsize(1) = max(1, numroc(this%norbs, nbl, st%dom_st_proc_grid%myrow, 0, st%dom_st_proc_grid%nprow))
        this%lsize(2) = max(1, numroc(this%norbs, nbl, st%dom_st_proc_grid%mycol, 0, st%dom_st_proc_grid%npcol))

        this%nproc(1) = st%dom_st_proc_grid%nprow
        this%nproc(2) = st%dom_st_proc_grid%npcol
        this%myroc(1) = st%dom_st_proc_grid%myrow
        this%myroc(2) = st%dom_st_proc_grid%mycol

        call descinit(this%desc(1), this%norbs, this%norbs, nbl, nbl, 0, 0, &
          st%dom_st_proc_grid%context, this%lsize(1), info)

        if(info /= 0) then
          write(message(1), '(a,i6)') 'descinit for BLACS failed with error code ', info
          call messages_fatal(1)
        end if

        this%calc_atom = .false.
        do iatom = 1, ions%natoms
          ibasis = this%atom_orb_basis(iatom, 1)

          do jatom = 1, ions%natoms
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
  subroutine lcao_run(namespace, space, gr, ions, st, ks, hm, st_start, lmm_r)
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(in)    :: ions
    type(states_elec_t),      intent(inout) :: st
    type(v_ks_t),             intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,        optional, intent(in)    :: st_start !< use for unoccupied-states run
    FLOAT,          optional, intent(in)    :: lmm_r !< used only if not present(st_start)

    type(lcao_t) :: lcao
    integer :: st_start_random
    logical :: lcao_done
    type(profile_t), save :: prof

    PUSH_SUB(lcao_run)

    if (present(st_start)) then
      ! If we are doing unocc calculation, do not mess with the correct eigenvalues
      ! of the occupied states.
      call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval=.not. present(st_start), calc_current=.false.)

      ASSERT(st_start >= 1)
      if (st_start > st%nst) then ! nothing to be done in LCAO
        POP_SUB(lcao_run)
        return
      end if
    end if

    call profiling_in(prof, 'LCAO_RUN')

    call lcao_init(lcao, namespace, space, gr, ions, st)

    call lcao_init_orbitals(lcao, st, gr, ions, start = st_start)

    if (.not. present(st_start)) then
      call lcao_guess_density(lcao, namespace, st, gr, hm, gr%sb, ions, st%qtot, st%d%ispin, st%d%spin_channels, st%rho)

      if (st%d%ispin > UNPOLARIZED) then
        ASSERT(present(lmm_r))
        call write_magnetic_moments(stdout, gr%mesh, st, ions, gr%der%boundaries, lmm_r)
      end if

      ! set up Hamiltonian (we do not call v_ks_h_setup here because we do not want to
      ! overwrite the guess density)
      message(1) = 'Info: Setting up Hamiltonian.'
      call messages_info(1)

      ! get the effective potential (we don`t need the eigenvalues yet)
      call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval=.false., calc_current=.false., calc_energy=.false.)
      ! eigenvalues have nevertheless to be initialized to something
      if(st%smear%method == SMEAR_SEMICONDUCTOR .and. lcao_is_available(lcao)) then
        st%eigenval = M_HUGE
      else
        !For smearing functions with finite temperature, we cannot set them to M_HUGE
        st%eigenval = M_ZERO
      end if

    end if

    lcao_done = .false.

    ! after initialized, can check that LCAO is possible
    if(lcao_is_available(lcao)) then
      lcao_done = .true.

      if(present(st_start)) then
        write(message(1),'(a,i8,a)') 'Performing LCAO for states ', st_start, ' and above'
        call messages_info(1)
      end if

      call lcao_wf(lcao, st, gr, ions, hm, namespace, start = st_start)

      if (.not. present(st_start)) then
        call states_elec_fermi(st, namespace, gr%mesh)
        call states_elec_write_eigenvalues(stdout, min(st%nst, lcao%norbs), st, space, hm%kpoints)

        ! Update the density and the Hamiltonian
        if (lcao%mode == OPTION__LCAOSTART__LCAO_FULL) then
          call v_ks_h_setup(namespace, space, gr, ions, st, ks, hm, calc_eigenval = .false., calc_current=.false.)
          if (st%d%ispin > UNPOLARIZED) then
            ASSERT(present(lmm_r))
            call write_magnetic_moments(stdout, gr%mesh, st, ions, gr%der%boundaries, lmm_r)
          end if
        end if
      end if
    end if

    if (.not. lcao_done .or. lcao%norbs < st%nst) then

      if(lcao_done) then
        st_start_random = lcao%norbs + 1
      else
        st_start_random = 1
      end if
      if(present(st_start)) st_start_random = max(st_start, st_start_random)

      if(st_start_random > 1) then
        write(message(1),'(a,i8,a)') 'Generating random wavefunctions for states ', st_start_random, ' and above'
        call messages_info(1)
      end if

      ! Randomly generate the initial wavefunctions.
      call states_elec_generate_random(st, gr%mesh, hm%kpoints, ist_start_ = st_start_random, normalized = .false.)

      call messages_write('Orthogonalizing wavefunctions.')
      call messages_info()
      call states_elec_orthogonalize(st, namespace, gr%mesh)

      if(.not. lcao_done) then
        ! If we are doing unocc calculation, do not mess with the correct eigenvalues and occupations
        ! of the occupied states.
        call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval=.not. present(st_start), calc_current=.false.) ! get potentials
        if(.not. present(st_start)) then
          call states_elec_fermi(st, namespace, gr%mesh) ! occupations
        end if

      end if

    else if (present(st_start)) then

      if(st_start > 1) then
        call messages_write('Orthogonalizing wavefunctions.')
        call messages_info()
        call states_elec_orthogonalize(st, namespace, gr%mesh)
      end if

    end if

    call lcao_end(lcao)


    call profiling_out(prof)
    POP_SUB(lcao_run)
  end subroutine lcao_run

  ! ---------------------------------------------------------
  subroutine lcao_end(this)
    type(lcao_t), intent(inout) :: this

    PUSH_SUB(lcao_end)

    SAFE_DEALLOCATE_A(this%calc_atom)
    SAFE_DEALLOCATE_A(this%norb_atom)
    SAFE_DEALLOCATE_A(this%basis_atom)
    SAFE_DEALLOCATE_A(this%basis_orb)
    SAFE_DEALLOCATE_A(this%atom_orb_basis)
    SAFE_DEALLOCATE_A(this%radius)
    SAFE_DEALLOCATE_A(this%sphere)
    SAFE_DEALLOCATE_A(this%orbitals)

    SAFE_DEALLOCATE_A(this%atom)
    SAFE_DEALLOCATE_A(this%level)
    SAFE_DEALLOCATE_A(this%ddim)
    SAFE_DEALLOCATE_A(this%cst)
    SAFE_DEALLOCATE_A(this%ck)
    SAFE_DEALLOCATE_A(this%dbuff)
    SAFE_DEALLOCATE_A(this%zbuff)

    this%initialized = .false.
    POP_SUB(lcao_end)
  end subroutine lcao_end


  ! ---------------------------------------------------------
  subroutine lcao_wf(this, st, gr, ions, hm, namespace, start)
    type(lcao_t),             intent(inout) :: this
    type(states_elec_t),      intent(inout) :: st
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(in)    :: ions
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(namespace_t),        intent(in)    :: namespace
    integer, optional,        intent(in)    :: start

    integer :: start_
    type(profile_t), save :: prof

    ASSERT(this%initialized)

    call profiling_in(prof, "LCAO")
    PUSH_SUB(lcao_wf)

    start_ = 1
    if(present(start)) start_ = start

    if(this%alternative) then
      if (states_are_real(st)) then
        call dlcao_alt_wf(this, st, gr, ions, hm, namespace, start_)
      else
        call zlcao_alt_wf(this, st, gr, ions, hm, namespace, start_)
      end if
    else
      if (states_are_real(st)) then
        call dlcao_wf(this, st, gr, ions, hm, namespace, start_)
      else
        call zlcao_wf(this, st, gr, ions, hm, namespace, start_)
      end if
    end if
    POP_SUB(lcao_wf)
    call profiling_out(prof)
  end subroutine lcao_wf


  ! ---------------------------------------------------------
  logical function lcao_is_available(this) result(available)
    type(lcao_t), intent(in) :: this

    PUSH_SUB(lcao_is_available)

    available = this%initialized .and. this%mode /= OPTION__LCAOSTART__LCAO_NONE

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

  !> This function deallocates a set of an atomic orbitals for an
  !! atom. It can be called when the batch is empty, in that case it
  !! does not do anything.
  subroutine lcao_alt_end_orbital(this, iatom)
    type(lcao_t),   intent(inout) :: this
    integer,        intent(in)    :: iatom

    PUSH_SUB(lcao_alt_end_orbital)

    if(this%is_orbital_initialized(iatom)) then
      call this%orbitals(iatom)%end()
      this%is_orbital_initialized(iatom) = .false.
    end if

    POP_SUB(lcao_alt_end_orbital)

  end subroutine lcao_alt_end_orbital

  ! ---------------------------------------------------------

  subroutine lcao_atom_density(this, namespace, st, mesh, sb, ions, iatom, spin_channels, rho)
    type(lcao_t),             intent(inout) :: this
    type(namespace_t),        intent(in)    :: namespace
    type(states_elec_t),      intent(in)    :: st
    type(mesh_t),             intent(in)    :: mesh
    type(simul_box_t),        intent(in)    :: sb
    type(ions_t),     target, intent(in)    :: ions
    integer,                  intent(in)    :: iatom
    integer,                  intent(in)    :: spin_channels
    FLOAT,                    intent(inout) :: rho(:, :) !< (gr%mesh%np, spin_channels)
    
    FLOAT, allocatable :: dorbital(:, :)
    CMPLX, allocatable :: zorbital(:, :)
    FLOAT, allocatable :: factors(:)
    FLOAT :: factor, aa
    integer :: iorb, ip, ii, ll, mm
    type(ps_t), pointer :: ps
    logical :: use_stored_orbitals

    PUSH_SUB(lcao_atom_density)

    rho = M_ZERO

    use_stored_orbitals = species_is_ps(ions%atom(iatom)%species) &
      .and. states_are_real(st) .and. spin_channels == 1 .and. lcao_is_available(this) &
      .and. st%d%dim == 1 .and. .not. ions%space%is_periodic()

    ps => species_ps(ions%atom(iatom)%species)

    ! we can use the orbitals we already have calculated
    if(use_stored_orbitals) then
      !There is no periodic copies here, so this will not work for periodic systems
      ASSERT(.not. ions%space%is_periodic())

      if(.not. this%alternative) then
        
        if(states_are_real(st)) then
          SAFE_ALLOCATE(dorbital(1:mesh%np, 1:st%d%dim))
        else
          SAFE_ALLOCATE(zorbital(1:mesh%np, 1:st%d%dim))
        end if
        
        do iorb = 1, this%norbs
          if(iatom /= this%atom(iorb)) cycle
          
          call species_iwf_ilm(ions%atom(iatom)%species, this%level(iorb), 1, ii, ll, mm)
          factor = ps%conf%occ(ii, 1)/(CNST(2.0)*ll + CNST(1.0))
         
          if(states_are_real(st)) then
            call dget_ao(this, st, mesh, ions, iorb, 1, dorbital, use_psi = .true.)
            !$omp parallel do
            do ip = 1, mesh%np
              rho(ip, 1) = rho(ip, 1) + factor*dorbital(ip, 1)**2
            end do
          else
            call zget_ao(this, st, mesh, ions, iorb, 1, zorbital, use_psi = .true.)
            !$omp parallel do
            do ip = 1, mesh%np
              rho(ip, 1) = rho(ip, 1) + factor*abs(zorbital(ip, 1))**2
            end do
          end if
          
        end do

        SAFE_DEALLOCATE_A(dorbital)
        SAFE_DEALLOCATE_A(zorbital)

      else

        ! for simplicity, always use real ones here.
        call dlcao_alt_get_orbital(this, this%sphere(iatom), ions, 1, iatom, this%norb_atom(iatom))

        ! the extra orbitals with the derivative are not relevant here, hence we divide by this%mult
        SAFE_ALLOCATE(factors(1:this%norb_atom(iatom)/this%mult))

        do iorb = 1, this%norb_atom(iatom)/this%mult
          call species_iwf_ilm(ions%atom(iatom)%species, iorb, 1, ii, ll, mm)
          factors(iorb) = ps%conf%occ(ii, 1)/(CNST(2.0)*ll + CNST(1.0))
        end do

        !$omp parallel do private(ip, aa, iorb) 
        do ip = 1, this%sphere(iatom)%np
          aa = CNST(0.0)
          do iorb = 1, this%norb_atom(iatom)/this%mult
            aa = aa + factors(iorb)*this%orbitals(iatom)%dff_linear(ip, iorb)**2
          end do
          !Due to the mapping function, more than one task could write to the same point in the array
          !$omp critical
          rho(this%sphere(iatom)%map(ip), 1) = aa
          !$omp end critical
        end do

        SAFE_DEALLOCATE_A(factors)

      end if

    else
      call species_atom_density(mesh, ions%space, namespace, sb, ions%atom(iatom), spin_channels, rho)
    end if

    POP_SUB(lcao_atom_density)
  end subroutine lcao_atom_density

  ! ---------------------------------------------------------
  !> builds a density which is the sum of the atomic densities
  subroutine lcao_guess_density(this, namespace, st, gr, hm, sb, ions, qtot, ispin, spin_channels, rho)
    type(lcao_t),             intent(inout) :: this
    type(namespace_t),        intent(in)    :: namespace
    type(states_elec_t),      intent(in)    :: st
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(simul_box_t),        intent(in)    :: sb
    type(ions_t),             intent(in)    :: ions
    FLOAT,                    intent(in)    :: qtot  !< the total charge of the system
    integer,                  intent(in)    :: ispin, spin_channels
    FLOAT,                    intent(out)   :: rho(:, :)

    integer :: ia, is, idir, gmd_opt, ip
    integer, save :: iseed = 321
    type(block_t) :: blk
    FLOAT :: rr, rnd, phi, theta, mag(1:3), lmag, n1, n2
    FLOAT, allocatable :: atom_rho(:,:)
    logical :: parallelized_in_atoms
    type(symmetrizer_t) :: symmetrizer


    PUSH_SUB(lcao_guess_density)

    parallelized_in_atoms = .false.

    if (spin_channels == 1) then
      gmd_opt = INITRHO_PARAMAGNETIC
    else
      !%Variable GuessMagnetDensity
      !%Type integer
      !%Default ferromagnetic
      !%Section SCF::LCAO
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
      call parse_variable(namespace, 'GuessMagnetDensity', INITRHO_FERROMAGNETIC, gmd_opt)
      if(.not.varinfo_valid_option('GuessMagnetDensity', gmd_opt)) call messages_input_error(namespace, 'GuessMagnetDensity')
      call messages_print_var_option(stdout, 'GuessMagnetDensity', gmd_opt)
    end if

    if(parse_is_defined(namespace, 'GuessMagnetDensity') .and. (hm%theory_level == HARTREE_FOCK &
        .or. hm%theory_level == GENERALIZED_KOHN_SHAM_DFT)) then
      message(1) = "GuessMagnetDensity cannot be used for Hartree-Fock and generalized Kohn-Sham calculation."
      message(2) = "Please perform a LDA or GGA calculation first and restart from this calculation."
      call messages_fatal(2)
    end if

    rho = M_ZERO
    select case (gmd_opt)
    case (INITRHO_PARAMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:gr%mesh%np, 1:spin_channels))

      parallelized_in_atoms = .true.

      do ia = ions%atoms_dist%start, ions%atoms_dist%end
        call lcao_atom_density(this, namespace, st, gr%mesh, sb, ions, ia, spin_channels, atom_rho)
        rho(1:gr%mesh%np, 1:spin_channels) = rho(1:gr%mesh%np, 1:spin_channels) + &
                                                  atom_rho(1:gr%mesh%np, 1:spin_channels)
      end do

      if (spin_channels == 2) then
        rho(1:gr%mesh%np, 1) = M_HALF*(sum(rho(1:gr%mesh%np, 1:2), dim=2))
        rho(1:gr%mesh%np, 2) = rho(1:gr%mesh%np, 1)
      end if

    case (INITRHO_FERROMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:gr%mesh%np, 1:2))

      parallelized_in_atoms = .true.

      rho = M_ZERO
      do ia = ions%atoms_dist%start, ions%atoms_dist%end
        call lcao_atom_density(this, namespace, st, gr%mesh, sb, ions, ia, 2, atom_rho(1:gr%mesh%np, 1:2))
        rho(1:gr%mesh%np, 1:2) = rho(1:gr%mesh%np, 1:2) + atom_rho(1:gr%mesh%np, 1:2)
      end do

    case (INITRHO_RANDOM) ! Randomly oriented spins
      SAFE_ALLOCATE(atom_rho(1:gr%mesh%np, 1:2))
      do ia = 1, ions%natoms
        call lcao_atom_density(this, namespace, st, gr%mesh, sb, ions, ia, 2, atom_rho)

        if (ispin == SPIN_POLARIZED) then
          call quickrnd(iseed, rnd)
          rnd = rnd - M_HALF
          if (rnd > M_ZERO) then
            rho(1:gr%mesh%np, 1:2) = rho(1:gr%mesh%np, 1:2) + atom_rho(1:gr%mesh%np, 1:2)
          else
            rho(1:gr%mesh%np, 1) = rho(1:gr%mesh%np, 1) + atom_rho(1:gr%mesh%np, 2)
            rho(1:gr%mesh%np, 2) = rho(1:gr%mesh%np, 2) + atom_rho(1:gr%mesh%np, 1)
          end if
        elseif (ispin == SPINORS) then
          call quickrnd(iseed, phi)
          call quickrnd(iseed, theta)
          phi = phi*M_TWO*M_PI
          theta = theta*M_PI*M_HALF
          
          call accumulate_rotated_density(gr%mesh, rho, atom_rho, theta, phi)
        end if
      end do

    case (INITRHO_USERDEF) ! User-defined
      
      !%Variable AtomsMagnetDirection
      !%Type block
      !%Section SCF::LCAO
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
      !% should have three components. If the norm of the vector is greater
      !% than the number of valence electrons in the atom, it will be rescaled
      !% to this number, which is the maximum possible magnetization.
      !%End
      if(parse_block(namespace, 'AtomsMagnetDirection', blk) < 0) then
        message(1) = "AtomsMagnetDirection block is not defined."
        call messages_fatal(1)
      end if

      if (parse_block_n(blk) /= ions%natoms) then
        message(1) = "AtomsMagnetDirection block has the wrong number of rows."
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(atom_rho(1:gr%mesh%np, 1:2))
      do ia = 1, ions%natoms
        !Read from AtomsMagnetDirection block 
        if (ispin == SPIN_POLARIZED) then
          call parse_block_float(blk, ia-1, 0, mag(1))
          mag(2:3) = M_ZERO !Else, this is unitialized and lead to a FPE in the case (lmag > n1+n2) 
          lmag = abs(mag(1))
        elseif (ispin == SPINORS) then
          do idir = 1, 3
            call parse_block_float(blk, ia-1, idir-1, mag(idir))
            if (abs(mag(idir)) < CNST(1.0e-20)) mag(idir) = M_ZERO
          end do
          lmag = sqrt(dot_product(mag(1:3), mag(1:3)))
        end if

        !Get atomic density
        call lcao_atom_density(this, namespace, st, gr%mesh, sb, ions, ia, 2, atom_rho)

        !Scale magnetization density
        n1 = dmf_integrate(gr%mesh, atom_rho(:, 1))
        n2 = dmf_integrate(gr%mesh, atom_rho(:, 2))
 
        if (lmag > n1 + n2) then
          mag = mag*(n1 + n2)/lmag
          lmag = n1 + n2
        elseif (abs(lmag) <= M_EPSILON) then
          if (abs(n1 - n2) <= M_EPSILON) then
           call lalg_axpy(gr%mesh%np, 2, M_ONE, atom_rho, rho)
          else
            !$omp parallel do simd
            do ip = 1, gr%mesh%np
              atom_rho(ip, 1) = M_HALF*(atom_rho(ip, 1) + atom_rho(ip, 2))
              rho(ip, 1) = rho(ip, 1) + atom_rho(ip, 1)
              rho(ip, 2) = rho(ip, 2) + atom_rho(ip, 1)
            end do
          end if
          cycle
        end if

        if (n1 - n2 /= lmag .and. n2 /= M_ZERO) then
          if (n1 - n2 < lmag) then
            call lalg_axpy(gr%mesh%np, (lmag - n1 + n2)/M_TWO/n2, atom_rho(:, 2), atom_rho(:, 1))
            call lalg_scal(gr%mesh%np, (n1 + n2 - lmag)/M_TWO/n2, atom_rho(:, 2))
          elseif (n1 - n2 > lmag) then
            call lalg_axpy(gr%mesh%np, (n1 - n2 - lmag)/M_TWO/n1, atom_rho(:, 1), atom_rho(:, 2))
            call lalg_scal(gr%mesh%np, (n1 + n2 + lmag)/M_TWO/n1, atom_rho(:, 1))
          end if
        end if

        !Rotate magnetization density
        if (ispin == SPIN_POLARIZED) then

          if (mag(1) > M_ZERO) then
            call lalg_axpy(gr%mesh%np, 2, M_ONE, atom_rho, rho)
          else
            call lalg_axpy(gr%mesh%np, M_ONE, atom_rho(:,2), rho(:,1))
            call lalg_axpy(gr%mesh%np, M_ONE, atom_rho(:,1), rho(:,2))
          end if

        elseif (ispin == SPINORS) then
          theta = acos(mag(3)/lmag)
          if (abs(mag(1)) <= M_EPSILON) then
            if (abs(mag(2)) <= M_EPSILON) then
              phi = M_ZERO
            elseif (mag(2) < M_ZERO) then
              phi = M_PI*CNST(3.0/2.0)
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
          theta = M_HALF*theta
          call accumulate_rotated_density(gr%mesh, rho, atom_rho, theta, phi)
        end if
      end do

      call parse_block_end(blk)

    end select


#ifdef HAVE_MPI
    if(ions%atoms_dist%parallel .and. parallelized_in_atoms) then
      ! NOTE: if random or user_defined are made parallelized in atoms, below should be st%d%nspin instead of spin_channels
      do is = 1, spin_channels
        call lalg_copy(gr%mesh%np, rho(:,is), atom_rho(:,1))
        call MPI_Allreduce(atom_rho(1, 1), rho(1, is), gr%mesh%np, &
          MPI_FLOAT, MPI_SUM, ions%atoms_dist%mpi_grp%comm, mpi_err)
      end do
    end if
#endif

    ! we now renormalize the density (necessary if we have a charged system)
    rr = M_ZERO
    do is = 1, spin_channels
      rr = rr + dmf_integrate(gr%mesh, rho(:, is), reduce = .false.)
    end do
    if(gr%mesh%parallel_in_domains) then
      call gr%mesh%allreduce(rr)
    end if

    write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', rr
    call messages_info(1)

    if (abs(rr) > M_EPSILON) then ! We only renormalize if the density is not zero
      rr = qtot / rr
      rho = rr * rho
    end if
    rr = M_ZERO
    do is = 1, spin_channels
      rr = rr + dmf_integrate(gr%mesh, rho(:, is), reduce = .false.)
    end do
    if(gr%mesh%parallel_in_domains) then
      call gr%mesh%allreduce(rr)
    end if

    write(message(1),'(a,f13.6)')'Info: Renormalized total charge = ', rr
    call messages_info(1)

    if(st%symmetrize_density) then
      call symmetrizer_init(symmetrizer, gr%mesh, gr%symm)

      do is = 1, st%d%nspin
        call dsymmetrizer_apply(symmetrizer, gr%mesh, field = rho(:, is), &
                                 symmfield = atom_rho(:, 1))
        call lalg_copy(gr%mesh%np, atom_rho(:, 1), rho(:, is))
      end do

      call symmetrizer_end(symmetrizer)
    end if

    SAFE_DEALLOCATE_A(atom_rho)
    POP_SUB(lcao_guess_density)
  end subroutine lcao_guess_density

  ! ---------------------------------------------------------
  subroutine accumulate_rotated_density(mesh, rho, atom_rho, theta, phi)
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(inout) :: rho(:,:)
    FLOAT,               intent(in)    :: atom_rho(:,:)
    FLOAT,               intent(in)    :: theta, phi

    integer :: ip

    PUSH_SUB(accumulate_rotated_density)

    !$omp parallel do simd
    do ip = 1, mesh%np
      rho(ip, 1) = rho(ip, 1) + cos(theta)**2*atom_rho(ip, 1) + sin(theta)**2*atom_rho(ip, 2)
      rho(ip, 2) = rho(ip, 2) + sin(theta)**2*atom_rho(ip, 1) + cos(theta)**2*atom_rho(ip, 2)
      rho(ip, 3) = rho(ip, 3) + cos(theta)*sin(theta)*cos(phi)*(atom_rho(ip, 1)-atom_rho(ip, 2))
      rho(ip, 4) = rho(ip, 4) - cos(theta)*sin(theta)*sin(phi)*(atom_rho(ip, 1)-atom_rho(ip, 2))
    end do

    POP_SUB(accumulate_rotated_density)
  end subroutine accumulate_rotated_density

  ! ---------------------------------------------------------

  subroutine lcao_init_orbitals(this, st, gr, ions, start)
    type(lcao_t),        intent(inout) :: this
    type(states_elec_t), intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    type(ions_t),        intent(in)    :: ions
    integer, optional,   intent(in)    :: start

    if(.not. lcao_is_available(this)) return
  
    PUSH_SUB(lcao_init_orbitals)
        
    if(.not. this%alternative) then
      if(states_are_real(st)) then
        call dinit_orbitals(this, st, gr, ions, start)
      else
        call zinit_orbitals(this, st, gr, ions, start)
      end if
    else
      if(states_are_real(st)) then
        call dlcao_alt_init_orbitals(this, st, gr, ions, start)
      else
        call zlcao_alt_init_orbitals(this, st, gr, ions, start)
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


end module lcao_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
