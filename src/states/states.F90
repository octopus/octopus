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

module states_oct_m
  use accel_oct_m
  use blacs_proc_grid_oct_m
  use boundaries_oct_m
  use calc_mode_par_oct_m
  use comm_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use derivatives_oct_m
  use distributed_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use kpoints_oct_m
  use loct_oct_m
  use loct_pointer_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use modelmb_particles_oct_m
  use mpi_oct_m
  use multicomm_oct_m
#ifdef HAVE_OPENMP
  use omp_lib
#endif
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use states_group_oct_m
  use states_dim_oct_m
  use symmetrizer_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private

  public ::                           &
    states_t,                         &
    states_priv_t,                    &
    states_init,                      &
    states_look,                      &
    states_densities_init,            &
    states_exec_init,                 &
    states_allocate_wfns,             &
    states_allocate_current,          &
    states_deallocate_wfns,           &
    states_null,                      &
    states_end,                       &
    states_copy,                      &
    states_generate_random,           &
    states_fermi,                     &
    states_eigenvalues_sum,           &
    states_calc_quantities,           &
    state_is_local,                   &
    state_kpt_is_local,               &
    states_distribute_nodes,          &
    states_wfns_memory,               &
    states_are_complex,               &
    states_are_real,                  &
    states_set_complex,               &
    states_get_state,                 &
    states_set_state,                 &
    states_get_points,                &
    states_pack,                      &
    states_unpack,                    &
    states_are_packed,                &
    states_write_info,                &
    states_set_zero,                  &
    states_block_min,                 &
    states_block_max,                 &
    states_block_size,                &
    states_count_pairs,               &
    occupied_states,                  &
    states_type,                      &
    states_set_phase

  type states_priv_t
    private
    type(type_t) :: wfs_type              !< real (TYPE_FLOAT) or complex (TYPE_CMPLX) wavefunctions
  end type states_priv_t

  type states_t
    ! Components are public by default
    type(states_dim_t)       :: d
    type(states_priv_t)      :: priv                  !< the private components
    integer                  :: nst                   !< Number of states in each irreducible subspace
    integer                  :: nst_conv              !< Number of states to be converged for unocc calc.

    logical                  :: only_userdef_istates  !< only use user-defined states as initial states in propagation
     
    type(states_group_t)     :: group

    !> used for the user-defined wavefunctions (they are stored as formula strings)
    !! (st%d%dim, st%nst, st%d%nik)
    character(len=1024), allocatable :: user_def_states(:,:,:)

    !> the densities and currents (after all we are doing DFT :)
    FLOAT, pointer :: rho(:,:)         !< rho(gr%mesh%np_part, st%d%nspin)
    FLOAT, pointer :: current(:, :, :) !<   current(gr%mesh%np_part, gr%sb%dim, st%d%nspin)

    !> k-point resolved current
    FLOAT, pointer :: current_kpt(:,:,:) !< current(gr%mesh%np_part, gr%sb%dim, kpt_start:kpt_end)


    FLOAT, pointer :: rho_core(:)      !< core charge for nl core corrections

    !> It may be required to "freeze" the deepest orbitals during the evolution; the density
    !! of these orbitals is kept in frozen_rho. It is different from rho_core.
    FLOAT, pointer :: frozen_rho(:, :)

    logical        :: calc_eigenval
    logical        :: uniform_occ   !< .true. if occupations are equal for all states: no empty states, and no smearing
    
    FLOAT, pointer :: eigenval(:,:) !< obviously the eigenvalues
    logical        :: fixed_occ     !< should the occupation numbers be fixed?
    logical        :: restart_fixed_occ !< should the occupation numbers be fixed by restart?
    logical        :: restart_reorder_occs !< used for restart with altered occupation numbers
    FLOAT, pointer :: occ(:,:)      !< the occupation numbers
    logical, private :: fixed_spins   !< In spinors mode, the spin direction is set
                                    !< for the initial (random) orbitals.
    FLOAT, pointer :: spin(:, :, :)

    FLOAT          :: qtot          !< (-) The total charge in the system (used in Fermi)
    FLOAT          :: val_charge    !< valence charge

    logical        :: fromScratch
    type(smear_t)  :: smear         ! smearing of the electronic occupations

    type(modelmb_particle_t) :: modelmbparticles
    integer,pointer:: mmb_nspindown(:,:) !< number of down spins in the selected Young diagram for each type and state
    integer,pointer:: mmb_iyoung(:,:)    !< index of the selected Young diagram for each type and state
    FLOAT, pointer :: mmb_proj(:)        !< projection of the state onto the chosen Young diagram

    !> This is stuff needed for the parallelization in states.
    logical                     :: parallel_in_states !< Am I parallel in states?
    type(mpi_grp_t)             :: mpi_grp            !< The MPI group related to the parallelization in states.
    type(mpi_grp_t)             :: dom_st_mpi_grp     !< The MPI group related to the domains-states "plane".
    type(mpi_grp_t)             :: st_kpt_mpi_grp     !< The MPI group related to the states-kpoints "plane".
    type(mpi_grp_t)             :: dom_st_kpt_mpi_grp !< The MPI group related to the domains-states-kpoints "cube".
#ifdef HAVE_SCALAPACK
    type(blacs_proc_grid_t)     :: dom_st_proc_grid   !< The BLACS process grid for the domains-states plane
#endif
    type(distributed_t)         :: dist
    logical                     :: scalapack_compatible !< Whether the states parallelization uses ScaLAPACK layout
    integer                     :: lnst               !< Number of states on local node.
    integer                     :: st_start, st_end   !< Range of states processed by local node.
    integer, pointer            :: node(:)            !< To which node belongs each state.
    type(multicomm_all_pairs_t), private :: ap        !< All-pairs schedule.

    logical                     :: symmetrize_density
    logical, private            :: packed

    integer                     :: randomization      !< Method used to generate random states
  end type states_t

  !> Method used to generate random states
  integer, public, parameter :: &
    PAR_INDEPENDENT = 1,              &
    PAR_DEPENDENT   = 2


  interface states_get_state
    module procedure dstates_get_state1, zstates_get_state1, dstates_get_state2, zstates_get_state2
    module procedure dstates_get_state3, zstates_get_state3, dstates_get_state4, zstates_get_state4
  end interface states_get_state

  interface states_set_state
    module procedure dstates_set_state1, zstates_set_state1, dstates_set_state2, zstates_set_state2
    module procedure dstates_set_state3, zstates_set_state3, dstates_set_state4, zstates_set_state4
  end interface states_set_state

  interface states_get_points
    module procedure dstates_get_points1, zstates_get_points1, dstates_get_points2, zstates_get_points2 
  end interface states_get_points

contains

  ! ---------------------------------------------------------
  subroutine states_null(st)
    type(states_t), intent(inout) :: st

    PUSH_SUB(states_null)

    call states_dim_null(st%d)
    call states_group_null(st%group)
    call distributed_nullify(st%dist)
    
    st%d%orth_method = 0
    call modelmb_particles_nullify(st%modelmbparticles)
    nullify(st%mmb_nspindown)
    nullify(st%mmb_iyoung)
    nullify(st%mmb_proj)

    st%priv%wfs_type = TYPE_FLOAT ! By default, calculations use real wavefunctions

    nullify(st%rho, st%current)
    nullify(st%current_kpt)
    nullify(st%rho_core, st%frozen_rho)
    nullify(st%eigenval, st%occ, st%spin)

    st%parallel_in_states = .false.
#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_nullify(st%dom_st_proc_grid)
#endif
    nullify(st%node)
    nullify(st%ap%schedule)

    st%packed = .false.

    POP_SUB(states_null)
  end subroutine states_null


  ! ---------------------------------------------------------
  subroutine states_init(st, parser, gr, geo)
    type(states_t), target, intent(inout) :: st
    type(parser_t),         intent(in)    :: parser
    type(grid_t),           intent(in)    :: gr
    type(geometry_t),       intent(in)    :: geo

    FLOAT :: excess_charge
    integer :: nempty, ntot, default, nthreads
    integer :: nempty_conv
    logical :: force

    PUSH_SUB(states_init)

    st%fromScratch = .true. ! this will be reset if restart_read is called
    call states_null(st)


    !%Variable SpinComponents
    !%Type integer
    !%Default unpolarized
    !%Section States
    !%Description
    !% The calculations may be done in three different ways: spin-restricted (TD)DFT (<i>i.e.</i>, doubly
    !% occupied "closed shells"), spin-unrestricted or "spin-polarized" (TD)DFT (<i>i.e.</i> we have two
    !% electronic systems, one with spin up and one with spin down), or making use of two-component
    !% spinors.
    !%Option unpolarized 1
    !% Spin-restricted calculations.
    !%Option polarized 2
    !%Option spin_polarized 2
    !% (Synonym <tt>polarized</tt>.) Spin-unrestricted, also known as spin-DFT, SDFT. This mode will double the number of
    !% wavefunctions necessary for a spin-unpolarized calculation.
    !%Option non_collinear 3
    !%Option spinors 3
    !% (Synonym: <tt>non_collinear</tt>.) The spin-orbitals are two-component spinors. This effectively allows the spin-density to
    !% be oriented non-collinearly: <i>i.e.</i> the magnetization vector is allowed to take different
    !% directions at different points. This vector is always in 3D regardless of <tt>Dimensions</tt>.
    !%End
    call parse_variable(parser, 'SpinComponents', UNPOLARIZED, st%d%ispin)
    if(.not.varinfo_valid_option('SpinComponents', st%d%ispin)) call messages_input_error('SpinComponents')
    call messages_print_var_option(stdout, 'SpinComponents', st%d%ispin)
    ! Use of spinors requires complex wavefunctions.
    if (st%d%ispin == SPINORS) st%priv%wfs_type = TYPE_CMPLX

    if(st%d%ispin /= UNPOLARIZED .and. gr%sb%kpoints%use_time_reversal) then
      message(1) = "Time reversal symmetry is only implemented for unpolarized spins."
      message(2) = "Use KPointsUseTimeReversal = no."
      call messages_fatal(2)
    end if
      

    !%Variable ExcessCharge
    !%Type float
    !%Default 0.0
    !%Section States
    !%Description
    !% The net charge of the system. A negative value means that we are adding
    !% electrons, while a positive value means we are taking electrons
    !% from the system.
    !%End
    call parse_variable(parser, 'ExcessCharge', M_ZERO, excess_charge)

    !%Variable CalcEigenvalues
    !%Type logical
    !%Default yes
    !%Section SCF
    !%Description
    !% (Experimental) When this variable is set to <tt>no</tt>,
    !% Octopus will not calculate the eigenvalues or eigenvectors of
    !% the Hamiltonian. Instead, Octopus will obtain the occupied
    !% subspace. The advantage that calculation can be made faster by
    !% avoiding subspace diagonalization and other calculations.
    !%
    !% This mode cannot be used with unoccupied states.    
    !%End
    call parse_variable(parser, 'CalcEigenvalues', .true., st%calc_eigenval)
    if(.not. st%calc_eigenval) call messages_experimental('CalcEigenvalues = .false.')
    
    !%Variable TotalStates
    !%Type integer
    !%Default 0
    !%Section States
    !%Description
    !% This variable sets the total number of states that Octopus will
    !% use. This is normally not necessary since by default Octopus
    !% sets the number of states to the minimum necessary to hold the
    !% electrons present in the system. (This default behavior is
    !% obtained by setting <tt>TotalStates</tt> to 0).
    !%
    !% If you want to add some unoccupied states, probably it is more convenient to use the variable
    !% <tt>ExtraStates</tt>.
    !%End
    call parse_variable(parser, 'TotalStates', 0, ntot)
    if (ntot < 0) then
      write(message(1), '(a,i5,a)') "Input: '", ntot, "' is not a valid value for TotalStates."
      call messages_fatal(1)
    end if

    !%Variable ExtraStates
    !%Type integer
    !%Default 0
    !%Section States
    !%Description
    !% The number of states is in principle calculated considering the minimum
    !% numbers of states necessary to hold the electrons present in the system.
    !% The number of electrons is
    !% in turn calculated considering the nature of the species supplied in the
    !% <tt>Species</tt> block, and the value of the <tt>ExcessCharge</tt> variable.
    !% However, one may command <tt>Octopus</tt> to use more states, which is necessary if one wants to
    !% use fractional occupational numbers, either fixed from the beginning through
    !% the <tt>Occupations</tt> block or by prescribing
    !% an electronic temperature with <tt>Smearing</tt>, or in order to calculate
    !% excited states (including with <tt>CalculationMode = unocc</tt>).
    !%End
    call parse_variable(parser, 'ExtraStates', 0, nempty)
    if (nempty < 0) then
      write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid value for ExtraStates."
      message(2) = '(0 <= ExtraStates)'
      call messages_fatal(2)
    end if

    if(ntot > 0 .and. nempty > 0) then
      message(1) = 'You cannot set TotalStates and ExtraStates at the same time.'
      call messages_fatal(1)
    end if

    !%Variable ExtraStatesToConverge
    !%Type integer
    !%Default 0
    !%Section States
    !%Description
    !% Only for unocc calculations.
    !% Specifies the number of extra states that will be considered for reaching the convergence.
    !% Together with <tt>ExtraStates</tt>, one can have some more states which will not be
    !% considered for the convergence criteria, thus making the convergence of the
    !% unocc calculation faster.
    !% By default, all extra states need to be converged.
    !%End
    call parse_variable(parser, 'ExtraStatesToConverge', nempty, nempty_conv)
    if (nempty < 0) then
      write(message(1), '(a,i5,a)') "Input: '", nempty_conv, "' is not a valid value for ExtraStatesToConverge."
      message(2) = '(0 <= ExtraStatesToConverge)'
      call messages_fatal(2)
    end if

    if(nempty_conv > nempty) then
      message(1) = 'You cannot set ExtraStatesToConverge to an higer value than ExtraStates.'
      call messages_fatal(1)
    end if

    ! For non-periodic systems this should just return the Gamma point
    call states_choose_kpoints(st%d, gr%sb)

    call geometry_val_charge(geo, st%val_charge)

    st%qtot = -(st%val_charge + excess_charge)

    if(st%qtot < -M_EPSILON) then
      write(message(1),'(a,f12.6,a)') 'Total charge = ', st%qtot, ' < 0'
      message(2) = 'Check Species and ExcessCharge.'
      call messages_fatal(2, only_root_writes = .true.)
    endif

    select case(st%d%ispin)
    case(UNPOLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%d%nspin = 1
      st%d%spin_channels = 1
    case(SPIN_POLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%d%nspin = 2
      st%d%spin_channels = 2
    case(SPINORS)
      st%d%dim = 2
      st%nst = int(st%qtot)
      if(st%nst < st%qtot) st%nst = st%nst + 1
      st%d%nspin = 4
      st%d%spin_channels = 2
    end select
    
    if(ntot > 0) then
      if(ntot < st%nst) then
        message(1) = 'TotalStates is smaller than the number of states required by the system.'
        call messages_fatal(1)
      end if

      st%nst = ntot
    end if

    st%nst_conv = st%nst + nempty_conv
    st%nst = st%nst + nempty
    if(st%nst == 0) then
      message(1) = "Cannot run with number of states = zero."
      call messages_fatal(1)
    end if

    !%Variable StatesBlockSize
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% Some routines work over blocks of eigenfunctions, which
    !% generally improves performance at the expense of increased
    !% memory consumption. This variable selects the size of the
    !% blocks to be used. If OpenCl is enabled, the default is 32;
    !% otherwise it is max(4, 2*nthreads).
    !%End

    nthreads = 1
#ifdef HAVE_OPENMP
    !$omp parallel
    !$omp master
    nthreads = omp_get_num_threads()
    !$omp end master
    !$omp end parallel
#endif    

    if(accel_is_enabled()) then
      default = 32
    else
      default = max(4, 2*nthreads)
    end if

    if(default > pad_pow2(st%nst)) default = pad_pow2(st%nst)

    ASSERT(default > 0)

    call parse_variable(parser, 'StatesBlockSize', default, st%d%block_size)
    if(st%d%block_size < 1) then
      call messages_write("The variable 'StatesBlockSize' must be greater than 0.")
      call messages_fatal()
    end if

    st%d%block_size = min(st%d%block_size, st%nst)
    conf%target_states_block_size = st%d%block_size

    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)

    ! Periodic systems require complex wavefunctions
    ! but not if it is Gamma-point only
    if(simul_box_is_periodic(gr%sb)) then
      if(.not. (kpoints_number(gr%sb%kpoints) == 1 .and. kpoints_point_is_gamma(gr%sb%kpoints, 1))) then
        st%priv%wfs_type = TYPE_CMPLX
      end if
    end if

    !%Variable OnlyUserDefinedInitialStates
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% If true, then only user-defined states from the block <tt>UserDefinedStates</tt>
    !% will be used as initial states for a time-propagation. No attempt is made
    !% to load ground-state orbitals from a previous ground-state run.
    !%End
    call parse_variable(parser, 'OnlyUserDefinedInitialStates', .false., st%only_userdef_istates)

    ! we now allocate some arrays
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    st%occ      = M_ZERO
    ! allocate space for formula strings that define user-defined states
    if(parse_is_defined(parser, 'UserDefinedStates') .or. parse_is_defined(parser, 'OCTInitialUserdefined') &
         .or. parse_is_defined(parser, 'OCTTargetUserdefined')) then
      SAFE_ALLOCATE(st%user_def_states(1:st%d%dim, 1:st%nst, 1:st%d%nik))
      ! initially we mark all 'formulas' as undefined
      st%user_def_states(1:st%d%dim, 1:st%nst, 1:st%d%nik) = 'undefined'
    end if

    if(st%d%ispin == SPINORS) then
      SAFE_ALLOCATE(st%spin(1:3, 1:st%nst, 1:st%d%nik))
    else
      nullify(st%spin)
    end if

    !%Variable StatesRandomization
    !%Type integer
    !%Default par_independent
    !%Section States
    !%Description
    !% The randomization of states can be done in two ways: 
    !% i) a parallelisation independent way (default), where the random states are identical, 
    !% irrespectively of the number of tasks and 
    !% ii) a parallelisation dependent way, which can prevent linear dependency
    !%  to occur for large systems.
    !%Option par_independent 1
    !% Parallelisation-independent randomization of states.
    !%Option par_dependent 2
    !% The randomization depends on the number of taks used in the calculation.
    !%End
    call parse_variable(parser, 'StatesRandomization', PAR_INDEPENDENT, st%randomization)


    call states_read_initial_occs(st, parser, excess_charge, gr%sb%kpoints)
    call states_read_initial_spins(st, parser)

    st%st_start = 1
    st%st_end = st%nst
    st%lnst = st%nst
    SAFE_ALLOCATE(st%node(1:st%nst))
    st%node(1:st%nst) = 0

    call mpi_grp_init(st%mpi_grp, -1)
    st%parallel_in_states = .false.

    call distributed_nullify(st%d%kpt, st%d%nik)

    call modelmb_particles_init(st%modelmbparticles, parser, gr)
    if (st%modelmbparticles%nparticle > 0) then
      ! FIXME: check why this is not initialized properly in the test, or why it is written out when not initialized
      SAFE_ALLOCATE(st%mmb_nspindown(1:st%modelmbparticles%ntype_of_particle, 1:st%nst))
      st%mmb_nspindown(:,:) = -1
      SAFE_ALLOCATE(st%mmb_iyoung(1:st%modelmbparticles%ntype_of_particle, 1:st%nst))
      st%mmb_iyoung(:,:) = -1
      SAFE_ALLOCATE(st%mmb_proj(1:st%nst))
      st%mmb_proj(:) = M_ZERO
    end if

    !%Variable SymmetrizeDensity
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% When enabled the density is symmetrized. Currently, this can
    !% only be done for periodic systems. (Experimental.)
    !%End
    call parse_variable(parser, 'SymmetrizeDensity', gr%sb%kpoints%use_symmetries, st%symmetrize_density)
    call messages_print_var_value(stdout, 'SymmetrizeDensity', st%symmetrize_density)

#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_nullify(st%dom_st_proc_grid)
#endif

    !%Variable ForceComplex
    !%Type logical
    !%Default no
    !%Section Execution::Debug
    !%Description
    !% Normally <tt>Octopus</tt> determines automatically the type necessary
    !% for the wavefunctions. When set to yes this variable will
    !% force the use of complex wavefunctions.
    !%
    !% Warning: This variable is designed for testing and
    !% benchmarking and normal users need not use it.
    !%End
    call parse_variable(parser, 'ForceComplex', .false., force)

    if(force) call states_set_complex(st)

    st%packed = .false.

    POP_SUB(states_init)
  end subroutine states_init

  !> Reads the 'states' file in the restart directory, and finds out
  !! the nik, dim, and nst contained in it.
  ! ---------------------------------------------------------
  subroutine states_look(restart, nik, dim, nst, ierr)
    type(restart_t), intent(in)  :: restart
    integer,         intent(out) :: nik
    integer,         intent(out) :: dim
    integer,         intent(out) :: nst
    integer,         intent(out) :: ierr

    character(len=256) :: lines(3)
    character(len=20)   :: char
    integer :: iunit

    PUSH_SUB(states_look)

    ierr = 0

    iunit = restart_open(restart, 'states')
    call restart_read(restart, iunit, lines, 3, ierr)
    if (ierr == 0) then
      read(lines(1), *) char, nst
      read(lines(2), *) char, dim
      read(lines(3), *) char, nik
    end if
    call restart_close(restart, iunit)

    POP_SUB(states_look)
  end subroutine states_look

  ! ---------------------------------------------------------
  !> Reads from the input file the initial occupations, if the
  !! block "Occupations" is present. Otherwise, it makes an initial
  !! guess for the occupations, maybe using the "Smearing"
  !! variable.
  !!
  !! The resulting occupations are placed on the st\%occ variable. The
  !! boolean st\%fixed_occ is also set to .true., if the occupations are
  !! set by the user through the "Occupations" block; false otherwise.
  subroutine states_read_initial_occs(st, parser, excess_charge, kpoints)
    type(states_t),  intent(inout) :: st
    type(parser_t),  intent(in)    :: parser
    FLOAT,           intent(in)    :: excess_charge
    type(kpoints_t), intent(in)    :: kpoints

    integer :: ik, ist, ispin, nspin, ncols, nrows, el_per_state, icol, start_pos, spin_n
    type(block_t) :: blk
    FLOAT :: rr, charge
    logical :: integral_occs, unoccupied_states
    FLOAT, allocatable :: read_occs(:, :)
    FLOAT :: charge_in_block

    PUSH_SUB(states_read_initial_occs)

    !%Variable RestartFixedOccupations
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% Setting this variable will make the restart proceed as
    !% if the occupations from the previous calculation had been set via the <tt>Occupations</tt> block,
    !% <i>i.e.</i> fixed. Otherwise, occupations will be determined by smearing.
    !%End
    call parse_variable(parser, 'RestartFixedOccupations', .false., st%restart_fixed_occ)
    ! we will turn on st%fixed_occ if restart_read is ever called

    !%Variable Occupations
    !%Type block
    !%Section States
    !%Description
    !% The occupation numbers of the orbitals can be fixed through the use of this
    !% variable. For example:
    !%
    !% <tt>%Occupations
    !% <br>&nbsp;&nbsp;2 | 2 | 2 | 2 | 2
    !% <br>%</tt>
    !%
    !% would fix the occupations of the five states to 2. There can be
    !% at most as many columns as states in the calculation. If there are fewer columns
    !% than states, then the code will assume that the user is indicating the occupations
    !% of the uppermost states where all lower states have full occupation (i.e. 2 for spin-unpolarized
    !% calculations, 1 otherwise) and all higher states have zero occupation. The first column
    !% will be taken to refer to the lowest state such that the occupations would be consistent
    !% with the correct total charge. For example, if there are 8 electrons and 10 states (from
    !% <tt>ExtraStates = 6</tt>), then an abbreviated specification
    !%
    !% <tt>%Occupations
    !% <br>&nbsp;&nbsp;1 | 0 | 1
    !% <br>%</tt>
    !%
    !% would be equivalent to a full specification
    !%
    !% <tt>%Occupations
    !% <br>&nbsp;&nbsp;2 | 2 | 2 | 1 | 0 | 1 | 0 | 0 | 0 | 0
    !% <br>%</tt>
    !%
    !% This is an example of use for constrained density-functional theory,
    !% crudely emulating a HOMO->LUMO+1 optical excitation.
    !% The number of rows should be equal
    !% to the number of k-points times the number of spins. For example, for a finite system
    !% with <tt>SpinComponents == spin_polarized</tt>,
    !% this block should contain two lines, one for each spin channel.
    !% All rows must have the same number of columns.
    !%
    !% The <tt>Occupations</tt> block is useful for the ground state of highly symmetric
    !% small systems (like an open-shell atom), to fix the occupation numbers
    !% of degenerate states in order to help <tt>octopus</tt> to converge. This is to
    !% be used in conjuction with <tt>ExtraStates</tt>. For example, to calculate the
    !% carbon atom, one would do:
    !%
    !% <tt>ExtraStates = 2
    !% <br>%Occupations
    !% <br>&nbsp;&nbsp;2 | 2/3 | 2/3 | 2/3
    !% <br>%</tt>
    !%
    !% If you want the calculation to be spin-polarized (which makes more sense), you could do:
    !%
    !% <tt>ExtraStates = 2
    !% <br>%Occupations
    !% <br>&nbsp;&nbsp; 2/3 | 2/3 | 2/3
    !% <br>&nbsp;&nbsp; 0   |   0 |   0
    !% <br>%</tt>
    !%
    !% Note that in this case the first state is absent, the code will calculate four states
    !% (two because there are four electrons, plus two because <tt>ExtraStates</tt> = 2), and since
    !% it finds only three columns, it will occupy the first state with one electron for each
    !% of the spin options.
    !%
    !% If the sum of occupations is not equal to the total charge set by <tt>ExcessCharge</tt>,
    !% an error message is printed.
    !% If <tt>FromScratch = no</tt> and <tt>RestartFixedOccupations = yes</tt>,
    !% this block will be ignored.
    !%End

    integral_occs = .true.

    occ_fix: if(parse_block(parser, 'Occupations', blk) == 0) then
      ! read in occupations
      st%fixed_occ = .true.

      ncols = parse_block_cols(blk, 0)
      if(ncols > st%nst) then
        message(1) = "Too many columns in block Occupations."
        call messages_warning(1)
        call messages_input_error("Occupations")
      end if

      nrows = parse_block_n(blk)
      if(nrows /= st%d%nik) then
        message(1) = "Wrong number of rows in block Occupations."
        call messages_warning(1)
        call messages_input_error("Occupations")
      end if

      do ik = 1, st%d%nik - 1
        if(parse_block_cols(blk, ik) /= ncols) then
          message(1) = "All rows in block Occupations must have the same number of columns."
          call messages_warning(1)
          call messages_input_error("Occupations")
        end if
      end do

      ! Now we fill all the "missing" states with the maximum occupation.
      if(st%d%ispin == UNPOLARIZED) then
        el_per_state = 2
      else
        el_per_state = 1
      end if

      SAFE_ALLOCATE(read_occs(1:ncols, 1:st%d%nik))

      charge_in_block = M_ZERO
      do ik = 1, st%d%nik
        do icol = 1, ncols
          call parse_block_float(blk, ik - 1, icol - 1, read_occs(icol, ik))
          charge_in_block = charge_in_block + read_occs(icol, ik) * st%d%kweights(ik)
        end do
      end do
      
      spin_n = 2
      select case(st%d%ispin)
        case(UNPOLARIZED) 
          spin_n = 2 
        case(SPIN_POLARIZED)
          spin_n = 2
        case(SPINORS)
          spin_n = 1         
      end select 
      
      start_pos = int((st%qtot - charge_in_block)/spin_n)

      if(start_pos + ncols > st%nst) then
        message(1) = "To balance charge, the first column in block Occupations is taken to refer to state"
        write(message(2),'(a,i6,a)') "number ", start_pos, " but there are too many columns for the number of states."
        write(message(3),'(a,i6,a)') "Solution: set ExtraStates = ", start_pos + ncols - st%nst
        call messages_fatal(3)
      end if

      do ik = 1, st%d%nik
        do ist = 1, start_pos
          st%occ(ist, ik) = el_per_state
        end do
      end do

      do ik = 1, st%d%nik
        do ist = start_pos + 1, start_pos + ncols
          st%occ(ist, ik) = read_occs(ist - start_pos, ik)
          integral_occs = integral_occs .and. &
            abs((st%occ(ist, ik) - el_per_state) * st%occ(ist, ik))  <=  M_EPSILON
        end do
      end do

      do ik = 1, st%d%nik
        do ist = start_pos + ncols + 1, st%nst
          st%occ(ist, ik) = M_ZERO
        end do
      end do

      call parse_block_end(blk)

      SAFE_DEALLOCATE_A(read_occs)

    else
      st%fixed_occ = .false.
      integral_occs = .false.

      ! first guess for occupation...paramagnetic configuration
      rr = M_ONE
      if(st%d%ispin == UNPOLARIZED) rr = M_TWO

      st%occ  = M_ZERO
      st%qtot = -(st%val_charge + excess_charge)

      nspin = 1
      if(st%d%nspin == 2) nspin = 2

      do ik = 1, st%d%nik, nspin
        charge = M_ZERO
        do ist = 1, st%nst
          do ispin = ik, ik + nspin - 1
            st%occ(ist, ispin) = min(rr, -(st%val_charge + excess_charge) - charge)
            charge = charge + st%occ(ist, ispin)
          end do
        end do
      end do

    end if occ_fix

    !%Variable RestartReorderOccs
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% Consider doing a ground-state calculation, and then restarting with new occupations set
    !% with the <tt>Occupations</tt> block, in an attempt to populate the orbitals of the original
    !% calculation. However, the eigenvalues may reorder as the density changes, in which case the
    !% occupations will now be referring to different orbitals. Setting this variable to yes will
    !% try to solve this issue when the restart data is being read, by reordering the occupations
    !% according to the order of the expectation values of the restart wavefunctions.
    !%End
    if(st%fixed_occ) then
      call parse_variable(parser, 'RestartReorderOccs', .false., st%restart_reorder_occs)
    else
      st%restart_reorder_occs = .false.
    end if

    call smear_init(st%smear, parser, st%d%ispin, st%fixed_occ, integral_occs, kpoints)

    unoccupied_states = (st%d%ispin /= SPINORS .and. st%nst*2 > st%qtot) .or. (st%d%ispin == SPINORS .and. st%nst > st%qtot)
    
    if(.not. smear_is_semiconducting(st%smear) .and. .not. st%smear%method == SMEAR_FIXED_OCC) then
      if(.not. unoccupied_states) then
        call messages_write('Smearing needs unoccupied states (via ExtraStates or TotalStates) to be useful.')
        call messages_warning()
      end if
    end if

    ! sanity check
    charge = M_ZERO
    do ist = 1, st%nst
      charge = charge + sum(st%occ(ist, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
    end do
    if(abs(charge - st%qtot) > CNST(1e-6)) then
      message(1) = "Initial occupations do not integrate to total charge."
      write(message(2), '(6x,f12.6,a,f12.6)') charge, ' != ', st%qtot
      call messages_fatal(2, only_root_writes = .true.)
    end if

    st%uniform_occ = smear_is_semiconducting(st%smear) .and. .not. unoccupied_states

    if(.not. st%calc_eigenval .and. .not. st%uniform_occ) then
      call messages_write('Calculation of the eigenvalues is required with unoccupied states', new_line = .true.)
      call messages_write('or smearing.')
      call messages_fatal()
    end if
    
    POP_SUB(states_read_initial_occs)
  end subroutine states_read_initial_occs


  ! ---------------------------------------------------------
  !> Reads, if present, the "InitialSpins" block. This is only
  !! done in spinors mode; otherwise the routine does nothing. The
  !! resulting spins are placed onto the st\%spin pointer. The boolean
  !! st\%fixed_spins is set to true if (and only if) the InitialSpins
  !! block is present.
  subroutine states_read_initial_spins(st, parser)
    type(states_t), intent(inout) :: st
    type(parser_t), intent(in)    :: parser

    integer :: i, j
    type(block_t) :: blk

    PUSH_SUB(states_read_initial_spins)

    st%fixed_spins = .false.
    if(st%d%ispin /= SPINORS) then
      POP_SUB(states_read_initial_spins)
      return
    end if

    !%Variable InitialSpins
    !%Type block
    !%Section States
    !%Description
    !% The spin character of the initial random guesses for the spinors can
    !% be fixed by making use of this block. Note that this will not "fix" the
    !% the spins during the calculation (this cannot be done in spinors mode, in
    !% being able to change the spins is why the spinors mode exists in the first
    !% place).
    !%
    !% This block is meaningless and ignored if the run is not in spinors mode
    !% (<tt>SpinComponents = spinors</tt>).
    !%
    !% The structure of the block is very simple: each column contains the desired
    !% <math>\left< S_x \right>, \left< S_y \right>, \left< S_z \right> </math> for each spinor.
    !% If the calculation is for a periodic system
    !% and there is more than one <i>k</i>-point, the spins of all the <i>k</i>-points are
    !% the same.
    !%
    !% For example, if we have two spinors, and we want one in the <math>S_x</math> "down" state,
    !% and another one in the <math>S_x</math> "up" state:
    !%
    !% <tt>%InitialSpins
    !% <br>&nbsp;&nbsp;&nbsp; 0.5 | 0.0 | 0.0
    !% <br>&nbsp;&nbsp; -0.5 | 0.0 | 0.0
    !% <br>%</tt>
    !%
    !% WARNING: if the calculation is for a system described by pseudopotentials (as
    !% opposed to user-defined potentials or model systems), this option is
    !% meaningless since the random spinors are overwritten by the atomic orbitals.
    !%
    !% This constraint must be fulfilled:
    !% <br><math> \left< S_x \right>^2 + \left< S_y \right>^2 + \left< S_z \right>^2 = \frac{1}{4} </math>
    !%End
    spin_fix: if(parse_block(parser, 'InitialSpins', blk)==0) then
      do i = 1, st%nst
        do j = 1, 3
          call parse_block_float(blk, i-1, j-1, st%spin(j, i, 1))
        end do
        if( abs(sum(st%spin(1:3, i, 1)**2) - M_FOURTH) > CNST(1.0e-6)) call messages_input_error('InitialSpins')
      end do
      call parse_block_end(blk)
      st%fixed_spins = .true.
      do i = 2, st%d%nik
        st%spin(:, :, i) = st%spin(:, :, 1)
      end do
    end if spin_fix

    POP_SUB(states_read_initial_spins)
  end subroutine states_read_initial_spins


  ! ---------------------------------------------------------
  !> Allocates the KS wavefunctions defined within a states_t structure.
  subroutine states_allocate_wfns(st, mesh, wfs_type, alloc_Left)
    type(states_t),         intent(inout)   :: st
    type(mesh_t),           intent(in)      :: mesh
    type(type_t), optional, intent(in)      :: wfs_type
    logical,      optional, intent(in)      :: alloc_Left !< allocate an additional set of wfs to store left eigenstates

    PUSH_SUB(states_allocate_wfns)

    if (present(wfs_type)) then
      ASSERT(wfs_type == TYPE_FLOAT .or. wfs_type == TYPE_CMPLX)
      st%priv%wfs_type = wfs_type
    end if

    call states_init_block(st, mesh)
    call states_set_zero(st)

    POP_SUB(states_allocate_wfns)
  end subroutine states_allocate_wfns

  !---------------------------------------------------------------------
  !> Initializes the data components in st that describe how the states
  !! are distributed in blocks:
  !!
  !! st\%nblocks: this is the number of blocks in which the states are divided. Note that
  !!   this number is the total number of blocks, regardless of how many are actually stored
  !!   in each node.
  !! block_start: in each node, the index of the first block.
  !! block_end: in each node, the index of the last block.
  !!   If the states are not parallelized, then block_start is 1 and block_end is st\%nblocks.
  !! st\%iblock(1:st\%nst, 1:st\%d\%nik): it points, for each state, to the block that contains it.
  !! st\%block_is_local(): st\%block_is_local(ib) is .true. if block ib is stored in the running node.
  !! st\%block_range(1:st\%nblocks, 1:2): Block ib contains states fromn st\%block_range(ib, 1) to st\%block_range(ib, 2)
  !! st\%block_size(1:st\%nblocks): Block ib contains a number st\%block_size(ib) of states.
  !! st\%block_initialized: it should be .false. on entry, and .true. after exiting this routine.
  !!
  !! The set of batches st\%psib(1:st\%nblocks) contains the blocks themselves.
  subroutine states_init_block(st, mesh, verbose)
    type(states_t),           intent(inout) :: st
    type(mesh_t),             intent(in)    :: mesh
    logical, optional,        intent(in)    :: verbose

    integer :: ib, iqn, ist
    logical :: same_node, verbose_
    integer, allocatable :: bstart(:), bend(:)

    PUSH_SUB(states_init_block)

    SAFE_ALLOCATE(bstart(1:st%nst))
    SAFE_ALLOCATE(bend(1:st%nst))
    SAFE_ALLOCATE(st%group%iblock(1:st%nst, 1:st%d%nik))

    st%group%iblock = 0

    verbose_ = optional_default(verbose, .true.)

    ! count and assign blocks
    ib = 0
    st%group%nblocks = 0
    bstart(1) = 1
    do ist = 1, st%nst
      INCR(ib, 1)

      st%group%iblock(ist, st%d%kpt%start:st%d%kpt%end) = st%group%nblocks + 1

      same_node = .true.
      if(st%parallel_in_states .and. ist /= st%nst) then
        ! We have to avoid that states that are in different nodes end
        ! up in the same block
        same_node = (st%node(ist + 1) == st%node(ist))
      end if

      if(ib == st%d%block_size .or. ist == st%nst .or. .not. same_node) then
        ib = 0
        INCR(st%group%nblocks, 1)
        bend(st%group%nblocks) = ist
        if(ist /= st%nst) bstart(st%group%nblocks + 1) = ist + 1
      end if
    end do

    SAFE_ALLOCATE(st%group%psib(1:st%group%nblocks, st%d%kpt%start:st%d%kpt%end))

    SAFE_ALLOCATE(st%group%block_is_local(1:st%group%nblocks, st%d%kpt%start:st%d%kpt%end))
    st%group%block_is_local = .false.
    st%group%block_start  = -1
    st%group%block_end    = -2  ! this will make that loops block_start:block_end do not run if not initialized

    do ib = 1, st%group%nblocks
      if(bstart(ib) >= st%st_start .and. bend(ib) <= st%st_end) then
        if(st%group%block_start == -1) st%group%block_start = ib
        st%group%block_end = ib
        do iqn = st%d%kpt%start, st%d%kpt%end
          st%group%block_is_local(ib, iqn) = .true.

          if (states_are_real(st)) then
            call batch_init(st%group%psib(ib, iqn), st%d%dim, bend(ib) - bstart(ib) + 1)
            call dbatch_allocate(st%group%psib(ib, iqn), bstart(ib), bend(ib), mesh%np_part, mirror = st%d%mirror_states)
          else
            call batch_init(st%group%psib(ib, iqn), st%d%dim, bend(ib) - bstart(ib) + 1)
            call zbatch_allocate(st%group%psib(ib, iqn), bstart(ib), bend(ib), mesh%np_part, mirror = st%d%mirror_states)
          end if
          
        end do
      end if
    end do

    SAFE_ALLOCATE(st%group%block_range(1:st%group%nblocks, 1:2))
    SAFE_ALLOCATE(st%group%block_size(1:st%group%nblocks))
    
    st%group%block_range(1:st%group%nblocks, 1) = bstart(1:st%group%nblocks)
    st%group%block_range(1:st%group%nblocks, 2) = bend(1:st%group%nblocks)
    st%group%block_size(1:st%group%nblocks) = bend(1:st%group%nblocks) - bstart(1:st%group%nblocks) + 1

    st%group%block_initialized = .true.

    SAFE_ALLOCATE(st%group%block_node(1:st%group%nblocks))

    ASSERT(associated(st%node))
    ASSERT(all(st%node >= 0) .and. all(st%node < st%mpi_grp%size))
    
    do ib = 1, st%group%nblocks
      st%group%block_node(ib) = st%node(st%group%block_range(ib, 1))
      ASSERT(st%group%block_node(ib) == st%node(st%group%block_range(ib, 2)))
    end do
    
    if(verbose_) then
      call messages_write('Info: Blocks of states')
      call messages_info()
      do ib = 1, st%group%nblocks
        call messages_write('      Block ')
        call messages_write(ib, fmt = 'i8')
        call messages_write(' contains ')
        call messages_write(st%group%block_size(ib), fmt = 'i8')
        call messages_write(' states')
        if(st%group%block_size(ib) > 0) then
          call messages_write(':')
          call messages_write(st%group%block_range(ib, 1), fmt = 'i8')
          call messages_write(' - ')
          call messages_write(st%group%block_range(ib, 2), fmt = 'i8')
        end if
        call messages_info()
      end do
    end if
    
!!$!!!!DEBUG
!!$    ! some debug output that I will keep here for the moment
!!$    if(mpi_grp_is_root(mpi_world)) then
!!$      print*, "NST       ", st%nst
!!$      print*, "BLOCKSIZE ", st%d%block_size
!!$      print*, "NBLOCKS   ", st%group%nblocks
!!$
!!$      print*, "==============="
!!$      do ist = 1, st%nst
!!$        print*, st%node(ist), ist, st%group%iblock(ist, 1)
!!$      end do
!!$      print*, "==============="
!!$
!!$      do ib = 1, st%group%nblocks
!!$        print*, ib, bstart(ib), bend(ib)
!!$      end do
!!$
!!$    end if
!!$!!!!ENDOFDEBUG

    SAFE_DEALLOCATE_A(bstart)
    SAFE_DEALLOCATE_A(bend)
    POP_SUB(states_init_block)
  end subroutine states_init_block


  ! ---------------------------------------------------------
  !> Deallocates the KS wavefunctions defined within a states_t structure.
  subroutine states_deallocate_wfns(st)
    type(states_t), intent(inout) :: st

    integer :: ib, iq

    PUSH_SUB(states_deallocate_wfns)

    if (st%group%block_initialized) then
       do ib = 1, st%group%nblocks
          do iq = st%d%kpt%start, st%d%kpt%end
            if(st%group%block_is_local(ib, iq)) then
              call batch_end(st%group%psib(ib, iq))
            end if
          end do
       end do

       SAFE_DEALLOCATE_P(st%group%psib)

       SAFE_DEALLOCATE_P(st%group%iblock)
       SAFE_DEALLOCATE_P(st%group%block_range)
       SAFE_DEALLOCATE_P(st%group%block_size)
       SAFE_DEALLOCATE_P(st%group%block_is_local)
       SAFE_DEALLOCATE_A(st%group%block_node)
       st%group%block_initialized = .false.
    end if

    POP_SUB(states_deallocate_wfns)
  end subroutine states_deallocate_wfns


  ! ---------------------------------------------------------
  subroutine states_densities_init(st, gr, geo)
    type(states_t), target, intent(inout) :: st
    type(grid_t),           intent(in)    :: gr
    type(geometry_t),       intent(in)    :: geo

    FLOAT :: fsize

    PUSH_SUB(states_densities_init)

    SAFE_ALLOCATE(st%rho(1:gr%fine%mesh%np_part, 1:st%d%nspin))
    st%rho = M_ZERO

    if(geo%nlcc) then
      SAFE_ALLOCATE(st%rho_core(1:gr%fine%mesh%np))
      st%rho_core(:) = M_ZERO
    end if

    fsize = gr%mesh%np_part*CNST(8.0)*st%d%block_size

    call messages_write('Info: states-block size = ')
    call messages_write(fsize, fmt = '(f10.1)', align_left = .true., units = unit_megabytes, print_units = .true.)
    call messages_info()

    POP_SUB(states_densities_init)
  end subroutine states_densities_init

  !---------------------------------------------------------------------
  subroutine states_allocate_current(st, gr)
    type(states_t), target, intent(inout) :: st
    type(grid_t),           intent(in)    :: gr

    PUSH_SUB(states_allocate_current)
    
    if(.not. associated(st%current)) then
      SAFE_ALLOCATE(st%current(1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:st%d%nspin))
      st%current = M_ZERO
    end if

    if(.not. associated(st%current_kpt)) then
      SAFE_ALLOCATE(st%current_kpt(1:gr%mesh%np_part,1:gr%mesh%sb%dim,st%d%kpt%start:st%d%kpt%end))
      st%current_kpt = M_ZERO
    end if

    POP_SUB(states_allocate_current)
  end subroutine states_allocate_current

  !---------------------------------------------------------------------
  !> This subroutine: (i) Fills in the block size (st\%d\%block_size);
  !! (ii) Finds out whether or not to pack the states (st\%d\%pack_states);
  !! (iii) Finds out the orthogonalization method (st\%d\%orth_method).
  subroutine states_exec_init(st, parser, mc)
    type(states_t),    intent(inout) :: st
    type(parser_t),    intent(in)    :: parser
    type(multicomm_t), intent(in)    :: mc

    integer :: default
    logical :: defaultl

    PUSH_SUB(states_exec_init)

    !%Variable StatesPack
    !%Type logical
    !%Section Execution::Optimization
    !%Description
    !% When set to yes, states are stored in packed mode, which improves
    !% performance considerably. Not all parts of the code will profit from
    !% this, but should nevertheless work regardless of how the states are
    !% stored.
    !%
    !% If OpenCL is used and this variable is set to yes, Octopus
    !% will store the wave-functions in device (GPU) memory. If
    !% there is not enough memory to store all the wave-functions,
    !% execution will stop with an error.
    !%
    !% See also the related <tt>HamiltonianApplyPacked</tt> variable.
    !%
    !% The default is yes except when using OpenCL.
    !%End

    defaultl = .true.
    if(accel_is_enabled()) then
      defaultl = .false.
    end if
    call parse_variable(parser, 'StatesPack', defaultl, st%d%pack_states)

    call messages_print_var_value(stdout, 'StatesPack', st%d%pack_states)

    !%Variable StatesMirror
    !%Type logical
    !%Section Execution::Optimization
    !%Description
    !% When this is enabled, Octopus keeps a copy of the states in
    !% main memory. This speeds up calculations when working with
    !% GPUs, as the memory does not to be copied back, but consumes
    !% more main memory.
    !%
    !% The default is false, except when acceleration is enabled and
    !% StatesPack is disabled.
    !%End

    defaultl = .false.
    if(accel_is_enabled() .and. .not. st%d%pack_states) then
      defaultl = .true.
    end if
    call parse_variable(parser, 'StatesMirror', defaultl, st%d%mirror_states)

    call messages_print_var_value(stdout, 'StatesMirror', st%d%mirror_states)

    !%Variable StatesOrthogonalization
    !%Type integer
    !%Section SCF::Eigensolver
    !%Description
    !% The full orthogonalization method used by some
    !% eigensolvers. The default is <tt>cholesky_serial</tt>, except with state
    !% parallelization, the default is <tt>cholesky_parallel</tt>.
    !%Option cholesky_serial 1
    !% Cholesky decomposition implemented using
    !% BLAS/LAPACK. Can be used with domain parallelization but not
    !% state parallelization.
    !%Option cholesky_parallel 2
    !% Cholesky decomposition implemented using
    !% ScaLAPACK. Compatible with states parallelization.
    !%Option cgs 3
    !% Classical Gram-Schmidt (CGS) orthogonalization.
    !% Can be used with domain parallelization but not state parallelization.
    !% The algorithm is defined in Giraud et al., Computers and Mathematics with Applications 50, 1069 (2005).
    !%Option mgs 4
    !% Modified Gram-Schmidt (MGS) orthogonalization.
    !% Can be used with domain parallelization but not state parallelization.
    !% The algorithm is defined in Giraud et al., Computers and Mathematics with Applications 50, 1069 (2005).
    !%Option drcgs 5
    !% Classical Gram-Schmidt orthogonalization with double-step reorthogonalization.
    !% Can be used with domain parallelization but not state parallelization.
    !% The algorithm is taken from Giraud et al., Computers and Mathematics with Applications 50, 1069 (2005). 
    !% According to this reference, this is much more precise than CGS or MGS algorithms. The MGS version seems not to improve much the stability and would require more communications over the domains.
    !%End

    default = OPTION__STATESORTHOGONALIZATION__CHOLESKY_SERIAL
#ifdef HAVE_SCALAPACK
    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then
      default = OPTION__STATESORTHOGONALIZATION__CHOLESKY_PARALLEL
    end if
#endif
    
    call parse_variable(parser, 'StatesOrthogonalization', default, st%d%orth_method)

    if(.not.varinfo_valid_option('StatesOrthogonalization', st%d%orth_method)) call messages_input_error('StatesOrthogonalization')
    call messages_print_var_option(stdout, 'StatesOrthogonalization', st%d%orth_method)

    !%Variable StatesCLDeviceMemory
    !%Type float
    !%Section Execution::Optimization
    !%Default -512
    !%Description
    !% This variable selects the amount of OpenCL device memory that
    !% will be used by Octopus to store the states. 
    !%
    !% A positive number smaller than 1 indicates a fraction of the total
    !% device memory. A number larger than one indicates an absolute
    !% amount of memory in megabytes. A negative number indicates an
    !% amount of memory in megabytes that would be subtracted from
    !% the total device memory.
    !%End
    call parse_variable(parser, 'StatesCLDeviceMemory', CNST(-512.0), st%d%cl_states_mem)

    POP_SUB(states_exec_init)
  end subroutine states_exec_init


  ! ---------------------------------------------------------
  subroutine states_copy(stout, stin, exclude_wfns, exclude_eigenval)
    type(states_t), target, intent(inout) :: stout
    type(states_t),         intent(in)    :: stin
    logical, optional,      intent(in)    :: exclude_wfns !< do not copy wavefunctions, densities, node
    logical, optional,      intent(in)    :: exclude_eigenval !< do not copy eigenvalues, occ, spin

    logical :: exclude_wfns_

    PUSH_SUB(states_copy)

    exclude_wfns_ = optional_default(exclude_wfns, .false.)

    call states_null(stout)

    call states_dim_copy(stout%d, stin%d)

    call modelmb_particles_copy(stout%modelmbparticles, stin%modelmbparticles)
    if (stin%modelmbparticles%nparticle > 0) then
      call loct_pointer_copy(stout%mmb_nspindown, stin%mmb_nspindown)
      call loct_pointer_copy(stout%mmb_iyoung, stin%mmb_iyoung)
      call loct_pointer_copy(stout%mmb_proj, stin%mmb_proj)
    end if

    stout%priv%wfs_type = stin%priv%wfs_type
    stout%nst           = stin%nst

    stout%only_userdef_istates = stin%only_userdef_istates

    if(.not. exclude_wfns_) call loct_pointer_copy(stout%rho, stin%rho)

    stout%calc_eigenval = stin%calc_eigenval
    stout%uniform_occ = stin%uniform_occ
    
    if(.not. optional_default(exclude_eigenval, .false.)) then
      call loct_pointer_copy(stout%eigenval, stin%eigenval)
      call loct_pointer_copy(stout%occ, stin%occ)
      call loct_pointer_copy(stout%spin, stin%spin)
    end if

    ! the call to init_block is done at the end of this subroutine
    ! it allocates iblock, psib, block_is_local
    stout%group%nblocks = stin%group%nblocks

    call loct_allocatable_copy(stout%user_def_states, stin%user_def_states)

    call loct_pointer_copy(stout%current, stin%current)
    call loct_pointer_copy(stout%current_kpt, stin%current_kpt)
 
    call loct_pointer_copy(stout%rho_core, stin%rho_core)
    call loct_pointer_copy(stout%frozen_rho, stin%frozen_rho)

    stout%fixed_occ = stin%fixed_occ
    stout%restart_fixed_occ = stin%restart_fixed_occ

    stout%fixed_spins = stin%fixed_spins

    stout%qtot       = stin%qtot
    stout%val_charge = stin%val_charge

    call smear_copy(stout%smear, stin%smear)

    stout%parallel_in_states = stin%parallel_in_states
    call mpi_grp_copy(stout%mpi_grp, stin%mpi_grp)
    stout%dom_st_kpt_mpi_grp = stin%dom_st_kpt_mpi_grp
    stout%st_kpt_mpi_grp     = stin%st_kpt_mpi_grp
    call loct_pointer_copy(stout%node, stin%node)

#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_copy(stin%dom_st_proc_grid, stout%dom_st_proc_grid)
#endif

    call distributed_copy(stin%dist, stout%dist)
    
    stout%scalapack_compatible = stin%scalapack_compatible

    stout%lnst       = stin%lnst
    stout%st_start   = stin%st_start
    stout%st_end     = stin%st_end

    if(stin%parallel_in_states) call multicomm_all_pairs_copy(stout%ap, stin%ap)

    stout%symmetrize_density = stin%symmetrize_density

    if(.not. exclude_wfns_) call states_group_copy(stin%d,stin%group, stout%group)

    stout%packed = stin%packed

    stout%randomization = stin%randomization

    POP_SUB(states_copy)
  end subroutine states_copy


  ! ---------------------------------------------------------
  subroutine states_end(st)
    type(states_t), intent(inout) :: st

    PUSH_SUB(states_end)

    call states_dim_end(st%d)

    if (st%modelmbparticles%nparticle > 0) then
      SAFE_DEALLOCATE_P(st%mmb_nspindown)
      SAFE_DEALLOCATE_P(st%mmb_iyoung)
      SAFE_DEALLOCATE_P(st%mmb_proj)
    end if
    call modelmb_particles_end(st%modelmbparticles)

    ! this deallocates dpsi, zpsi, psib, iblock, iblock
    call states_deallocate_wfns(st)

    SAFE_DEALLOCATE_A(st%user_def_states)

    SAFE_DEALLOCATE_P(st%rho)
    SAFE_DEALLOCATE_P(st%eigenval)

    SAFE_DEALLOCATE_P(st%current)
    SAFE_DEALLOCATE_P(st%current_kpt)
    SAFE_DEALLOCATE_P(st%rho_core)
    SAFE_DEALLOCATE_P(st%frozen_rho)
    SAFE_DEALLOCATE_P(st%occ)
    SAFE_DEALLOCATE_P(st%spin)

#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_end(st%dom_st_proc_grid)
#endif

    call distributed_end(st%dist)

    SAFE_DEALLOCATE_P(st%node)

    if(st%parallel_in_states) then
      SAFE_DEALLOCATE_P(st%ap%schedule)
    end if

    POP_SUB(states_end)
  end subroutine states_end

  ! ---------------------------------------------------------
  !> generate a hydrogen s-wavefunction around a random point
  subroutine states_generate_random(st, mesh, sb, ist_start_, ist_end_, ikpt_start_, ikpt_end_, normalized)
    type(states_t),    intent(inout) :: st
    type(mesh_t),      intent(in)    :: mesh
    type(simul_box_t), intent(in)    :: sb
    integer, optional, intent(in)    :: ist_start_
    integer, optional, intent(in)    :: ist_end_
    integer, optional, intent(in)    :: ikpt_start_
    integer, optional, intent(in)    :: ikpt_end_
    logical, optional, intent(in)    :: normalized !< whether generate states should have norm 1, true by default
    
    integer :: ist, ik, id, ist_start, ist_end, jst, ikpt_start, ikpt_end
    CMPLX   :: alpha, beta
    FLOAT, allocatable :: dpsi(:,  :)
    CMPLX, allocatable :: zpsi(:,  :), zpsi2(:)
    integer :: ikpoint, ip

    PUSH_SUB(states_generate_random)
 
    ist_start = optional_default(ist_start_, 1)
    ist_end = optional_default(ist_end_, st%nst)
    ikpt_start = optional_default(ikpt_start_, 1)
    ikpt_end = optional_default(ikpt_end_, st%d%nik)

    SAFE_ALLOCATE(dpsi(1:mesh%np, 1:st%d%dim))
    if (states_are_complex(st)) then
      SAFE_ALLOCATE(zpsi(1:mesh%np, 1:st%d%dim))
    end if

    select case(st%d%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)

      do ik = ikpt_start, ikpt_end
        ikpoint = states_dim_get_kpoint_index(st%d, ik)
        do ist = ist_start, ist_end
          if (states_are_real(st).or.kpoints_point_is_gamma(sb%kpoints, ikpoint)) then
            if(st%randomization == PAR_INDEPENDENT) then
              call dmf_random(mesh, dpsi(:, 1), mesh%vp%xlocal-1, normalized = normalized)
            else
              call dmf_random(mesh, dpsi(:, 1), normalized = normalized)
            end if
            if(.not. state_kpt_is_local(st, ist, ik)) cycle
            if(states_are_complex(st)) then !Gamma point
              forall(ip=1:mesh%np) 
                zpsi(ip,1) = cmplx(dpsi(ip,1), M_ZERO)
              end forall
              call states_set_state(st, mesh, ist,  ik, zpsi)
            else
              call states_set_state(st, mesh, ist,  ik, dpsi)
            end if
          else
            if(st%randomization == PAR_INDEPENDENT) then
              call zmf_random(mesh, zpsi(:, 1), mesh%vp%xlocal-1, normalized = normalized)
            else
              call zmf_random(mesh, zpsi(:, 1), normalized = normalized)
            end if
            if(.not. state_kpt_is_local(st, ist, ik)) cycle
            call states_set_state(st, mesh, ist,  ik, zpsi)
          end if
        end do
      end do

    case(SPINORS)

      ASSERT(states_are_complex(st))

      if(st%fixed_spins) then

        do ik = ikpt_start, ikpt_end
          ikpoint = states_dim_get_kpoint_index(st%d, ik)
          do ist = ist_start, ist_end
            if(kpoints_point_is_gamma(sb%kpoints, ikpoint)) then
              if(st%randomization == PAR_INDEPENDENT) then
                call dmf_random(mesh, dpsi(:, 1), mesh%vp%xlocal-1, normalized = normalized)
              else
                call dmf_random(mesh, dpsi(:, 1), normalized = normalized)
                if(.not. state_kpt_is_local(st, ist, ik)) cycle
              end if
              forall(ip=1:mesh%np)
                zpsi(ip,1) = cmplx(dpsi(ip,1), M_ZERO)
              end forall
              call states_set_state(st, mesh, ist,  ik, zpsi)
            else
              if(st%randomization == PAR_INDEPENDENT) then
                call zmf_random(mesh, zpsi(:, 1), mesh%vp%xlocal-1, normalized = normalized)
              else
                call zmf_random(mesh, zpsi(:, 1), normalized = normalized)
                if(.not. state_kpt_is_local(st, ist, ik)) cycle
              end if
            end if
            if(.not. state_kpt_is_local(st, ist, ik)) cycle
            ! In this case, the spinors are made of a spatial part times a vector [alpha beta]^T in
            ! spin space (i.e., same spatial part for each spin component). So (alpha, beta)
            ! determines the spin values. The values of (alpha, beta) can be be obtained
            ! with simple formulae from <Sx>, <Sy>, <Sz>.
            !
            ! Note that here we orthonormalize the orbital part. This ensures that the spinors
            ! are untouched later in the general orthonormalization, and therefore the spin values
            ! of each spinor remain the same.
            SAFE_ALLOCATE(zpsi2(1:mesh%np))
            do jst = ist_start, ist - 1
              call states_get_state(st, mesh, 1, jst, ik, zpsi2)
              zpsi(1:mesh%np, 1) = zpsi(1:mesh%np, 1) - zmf_dotp(mesh, zpsi(:, 1), zpsi2)*zpsi2(1:mesh%np)
            end do
            SAFE_DEALLOCATE_A(zpsi2)
            
            call zmf_normalize(mesh, 1, zpsi)
            zpsi(1:mesh%np, 2) = zpsi(1:mesh%np, 1)

            alpha = TOCMPLX(sqrt(M_HALF + st%spin(3, ist, ik)), M_ZERO)
            beta  = TOCMPLX(sqrt(M_ONE - abs(alpha)**2), M_ZERO)
            if(abs(alpha) > M_ZERO) then
              beta = TOCMPLX(st%spin(1, ist, ik) / abs(alpha), st%spin(2, ist, ik) / abs(alpha))
            end if
            zpsi(1:mesh%np, 1) = alpha*zpsi(1:mesh%np, 1)
            zpsi(1:mesh%np, 2) = beta*zpsi(1:mesh%np, 2)
            call states_set_state(st, mesh, ist,  ik, zpsi)
          end do
        end do
      else
        do ik = ikpt_start, ikpt_end
          do ist = ist_start, ist_end
            do id = 1, st%d%dim
              if(st%randomization == PAR_INDEPENDENT) then
                call zmf_random(mesh, zpsi(:, id), mesh%vp%xlocal-1, normalized = normalized)
              else
                call zmf_random(mesh, zpsi(:, id), normalized = normalized)
              end if
            end do
            if(.not. state_kpt_is_local(st, ist, ik)) cycle
            call states_set_state(st, mesh, ist,  ik, zpsi)
          end do
        end do
      end if

    end select

    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(zpsi)

    POP_SUB(states_generate_random)
  end subroutine states_generate_random

  ! ---------------------------------------------------------
  subroutine states_fermi(st, mesh)
    type(states_t), intent(inout) :: st
    type(mesh_t),   intent(in)    :: mesh

    !> Local variables.
    integer            :: ist, ik
    FLOAT              :: charge
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(states_fermi)

    call smear_find_fermi_energy(st%smear, st%eigenval, st%occ, st%qtot, &
      st%d%nik, st%nst, st%d%kweights)

    call smear_fill_occupations(st%smear, st%eigenval, st%occ, &
      st%d%nik, st%nst)
        
    ! check if everything is OK
    charge = M_ZERO
    do ist = 1, st%nst
      charge = charge + sum(st%occ(ist, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
    end do
    if(abs(charge-st%qtot) > CNST(1e-6)) then
      message(1) = 'Occupations do not integrate to total charge.'
      write(message(2), '(6x,f12.8,a,f12.8)') charge, ' != ', st%qtot
      call messages_warning(2)
      if(charge < M_EPSILON) then
        message(1) = "There don't seem to be any electrons at all!"
        call messages_fatal(1)
      end if
    end if

    if(st%d%ispin == SPINORS) then
      ASSERT(states_are_complex(st))
      
      st%spin(:,:,:) = M_ZERO
      
      SAFE_ALLOCATE(zpsi(1:mesh%np, st%d%dim))
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call states_get_state(st, mesh, ist, ik, zpsi)
          st%spin(1:3, ist, ik) = state_spin(mesh, zpsi)
        end do
      end do
      SAFE_DEALLOCATE_A(zpsi)

#if defined(HAVE_MPI)        
        if(st%parallel_in_states .or. st%d%kpt%parallel) then
          call comm_allreduce(st%st_kpt_mpi_grp%comm, st%spin)
        end if
#endif      
            
    end if

    POP_SUB(states_fermi)
  end subroutine states_fermi


  ! ---------------------------------------------------------
  !> function to calculate the eigenvalues sum using occupations as weights
  function states_eigenvalues_sum(st, alt_eig) result(tot)
    type(states_t),  intent(in) :: st
    FLOAT, optional, intent(in) :: alt_eig(st%st_start:, st%d%kpt%start:) !< (:st%st_end, :st%d%kpt%end)
    FLOAT                       :: tot

    integer :: ik

    PUSH_SUB(states_eigenvalues_sum)

    tot = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      if(present(alt_eig)) then
        tot = tot + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik) * &
          alt_eig(st%st_start:st%st_end, ik))
      else
        tot = tot + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik) * &
          st%eigenval(st%st_start:st%st_end, ik))
      end if
    end do

    if(st%parallel_in_states .or. st%d%kpt%parallel) call comm_allreduce(st%st_kpt_mpi_grp%comm, tot)

    POP_SUB(states_eigenvalues_sum)
  end function states_eigenvalues_sum


  ! ---------------------------------------------------------
  subroutine states_distribute_nodes(st, mc)
    type(states_t),    intent(inout) :: st
    type(multicomm_t), intent(in)    :: mc

    PUSH_SUB(states_distribute_nodes)

    ! Defaults.
    st%node(:)            = 0
    st%st_start           = 1
    st%st_end             = st%nst
    st%lnst               = st%nst
    st%parallel_in_states = .false.
    call mpi_grp_init(st%mpi_grp, mc%group_comm(P_STRATEGY_STATES))
    call mpi_grp_init(st%dom_st_kpt_mpi_grp, mc%dom_st_kpt_comm)
    call mpi_grp_init(st%dom_st_mpi_grp, mc%dom_st_comm)
    call mpi_grp_init(st%st_kpt_mpi_grp, mc%st_kpt_comm)

#ifdef HAVE_SCALAPACK
    !%Variable ScaLAPACKCompatible
    !%Type logical
    !%Section Execution::Parallelization
    !%Description
    !% Whether to use a layout for states parallelization which is compatible with ScaLAPACK.
    !% The default is yes for <tt>CalculationMode = gs, unocc, go</tt> without k-point parallelization,
    !% and no otherwise. (Setting to other than default is experimental.)
    !% The value must be yes if any ScaLAPACK routines are called in the course of the run;
    !% it must be set by hand for <tt>td</tt> with <tt>TDDynamics = bo</tt>.
    !% This variable has no effect unless you are using states parallelization and have linked ScaLAPACK.
    !% Note: currently, use of ScaLAPACK is not compatible with task parallelization (<i>i.e.</i> slaves).
    !%End
    call parse_variable(parser, 'ScaLAPACKCompatible', &
      calc_mode_par_scalapack_compat() .and. .not. st%d%kpt%parallel, st%scalapack_compatible)
    if((calc_mode_par_scalapack_compat() .and. .not. st%d%kpt%parallel) .neqv. st%scalapack_compatible) &
      call messages_experimental('Setting ScaLAPACKCompatible to other than default')

    if(st%scalapack_compatible) then
      if(multicomm_have_slaves(mc)) &
        call messages_not_implemented("ScaLAPACK usage with task parallelization (slaves)")
      call blacs_proc_grid_init(st%dom_st_proc_grid, st%dom_st_mpi_grp)
    else
      call blacs_proc_grid_nullify(st%dom_st_proc_grid)
    end if
#else
    st%scalapack_compatible = .false.
#endif

    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then

#ifdef HAVE_MPI
      call multicomm_create_all_pairs(st%mpi_grp, st%ap)
#endif

      if(st%nst < st%mpi_grp%size) then
        message(1) = "Have more processors than necessary"
        write(message(2),'(i4,a,i4,a)') st%mpi_grp%size, " processors and ", st%nst, " states."
        call messages_fatal(2)
      end if

      call distributed_init(st%dist, st%nst, st%mpi_grp%comm, "states", scalapack_compat = st%scalapack_compatible)

      st%st_start = st%dist%start
      st%st_end   = st%dist%end
      st%lnst     = st%dist%nlocal
      st%node(1:st%nst) = st%dist%node(1:st%nst)
      st%parallel_in_states = st%dist%parallel

    end if

    POP_SUB(states_distribute_nodes)
  end subroutine states_distribute_nodes


  ! ---------------------------------------------------------
  subroutine states_set_complex(st)
    type(states_t),    intent(inout) :: st

    PUSH_SUB(states_set_complex)
    st%priv%wfs_type = TYPE_CMPLX

    POP_SUB(states_set_complex)
  end subroutine states_set_complex

  ! ---------------------------------------------------------
  pure logical function states_are_complex(st) result (wac)
    type(states_t),    intent(in) :: st

    wac = (st%priv%wfs_type == TYPE_CMPLX)

  end function states_are_complex


  ! ---------------------------------------------------------
  pure logical function states_are_real(st) result (war)
    type(states_t),    intent(in) :: st

    war = (st%priv%wfs_type == TYPE_FLOAT)

  end function states_are_real

  ! ---------------------------------------------------------


  pure type(type_t) function states_type(st)
    type(states_t),    intent(in) :: st
    
    states_type = st%priv%wfs_type
    
  end function states_type
  
  
  ! ---------------------------------------------------------
  !
  !> This function can calculate several quantities that depend on
  !! derivatives of the orbitals from the states and the density.
  !! The quantities to be calculated depend on the arguments passed.
  subroutine states_calc_quantities(der, st, nlcc, &
    kinetic_energy_density, paramagnetic_current, density_gradient, density_laplacian, gi_kinetic_energy_density)
    type(derivatives_t),     intent(in)    :: der
    type(states_t),          intent(in)    :: st
    logical,                 intent(in)    :: nlcc
    FLOAT, optional, target, intent(out)   :: kinetic_energy_density(:,:)       !< The kinetic energy density.
    FLOAT, optional, target, intent(out)   :: paramagnetic_current(:,:,:)       !< The paramagnetic current.
    FLOAT, optional,         intent(out)   :: density_gradient(:,:,:)           !< The gradient of the density.
    FLOAT, optional,         intent(out)   :: density_laplacian(:,:)            !< The Laplacian of the density.
    FLOAT, optional,         intent(out)   :: gi_kinetic_energy_density(:,:)    !< The gauge-invariant kinetic energy density.

    FLOAT, pointer :: jp(:, :, :)
    FLOAT, pointer :: tau(:, :)
    CMPLX, allocatable :: wf_psi(:,:), gwf_psi(:,:,:), wf_psi_conj(:,:), lwf_psi(:,:)
    FLOAT, allocatable :: abs_wf_psi(:), abs_gwf_psi(:)
    CMPLX, allocatable :: psi_gpsi(:)
    CMPLX   :: c_tmp
    integer :: is, ik, ist, i_dim, st_dim, ii
    FLOAT   :: ww, kpoint(1:MAX_DIM)
    logical :: something_to_do
    FLOAT, allocatable :: symm(:, :)
    type(symmetrizer_t) :: symmetrizer
    type(profile_t), save :: prof

    call profiling_in(prof, "STATES_CALC_QUANTITIES")

    PUSH_SUB(states_calc_quantities)

    something_to_do = present(kinetic_energy_density) .or. present(gi_kinetic_energy_density) .or. &
      present(paramagnetic_current) .or. present(density_gradient) .or. present(density_laplacian)
    ASSERT(something_to_do)

    SAFE_ALLOCATE( wf_psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE( wf_psi_conj(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gwf_psi(1:der%mesh%np, 1:der%mesh%sb%dim, 1:st%d%dim))
    SAFE_ALLOCATE(abs_wf_psi(1:der%mesh%np))
    SAFE_ALLOCATE(abs_gwf_psi(1:der%mesh%np))
    SAFE_ALLOCATE(psi_gpsi(1:der%mesh%np))
    if(present(density_laplacian)) then
      SAFE_ALLOCATE(lwf_psi(1:der%mesh%np, 1:st%d%dim))
    end if

    nullify(tau)
    if(present(kinetic_energy_density)) tau => kinetic_energy_density

    nullify(jp)
    if(present(paramagnetic_current)) jp => paramagnetic_current

    ! for the gauge-invariant kinetic energy density we need the
    ! current and the kinetic energy density
    if(present(gi_kinetic_energy_density)) then
      if(.not. present(paramagnetic_current) .and. states_are_complex(st)) then
        SAFE_ALLOCATE(jp(1:der%mesh%np, 1:der%mesh%sb%dim, 1:st%d%nspin))
      end if
      if(.not. present(kinetic_energy_density)) then
        SAFE_ALLOCATE(tau(1:der%mesh%np, 1:st%d%nspin))
      end if
    end if

    if(associated(tau)) tau = M_ZERO
    if(associated(jp)) jp = M_ZERO
    if(present(density_gradient)) density_gradient(:,:,:) = M_ZERO
    if(present(density_laplacian)) density_laplacian(:,:) = M_ZERO
    if(present(gi_kinetic_energy_density)) gi_kinetic_energy_density = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end

      kpoint(1:der%mesh%sb%dim) = kpoints_get_point(der%mesh%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))
      is = states_dim_get_spin_index(st%d, ik)

      do ist = st%st_start, st%st_end

        ww = st%d%kweights(ik)*st%occ(ist, ik)
        if(abs(ww) <= M_EPSILON) cycle

        ! all calculations will be done with complex wavefunctions
        call states_get_state(st, der%mesh, ist, ik, wf_psi)

        do st_dim = 1, st%d%dim
          call boundaries_set(der%boundaries, wf_psi(:, st_dim))
        end do

        ! calculate gradient of the wavefunction
        do st_dim = 1, st%d%dim
          call zderivatives_grad(der, wf_psi(:,st_dim), gwf_psi(:,:,st_dim), set_bc = .false.)
        end do

        ! calculate the Laplacian of the wavefunction
        if (present(density_laplacian)) then
          do st_dim = 1, st%d%dim
            call zderivatives_lapl(der, wf_psi(:,st_dim), lwf_psi(:,st_dim), set_bc = .false.)
          end do
        end if

        !We precompute some quantites, to avoid to compute it many times
        wf_psi_conj(1:der%mesh%np, 1:st%d%dim) = conjg(wf_psi(1:der%mesh%np,1:st%d%dim))
        abs_wf_psi(1:der%mesh%np) = real(wf_psi_conj(1:der%mesh%np, 1)*wf_psi(1:der%mesh%np, 1))

        if(present(density_laplacian)) then
          density_laplacian(1:der%mesh%np, is) = density_laplacian(1:der%mesh%np, is) + &
               ww*M_TWO*real(wf_psi_conj(1:der%mesh%np, 1)*lwf_psi(1:der%mesh%np, 1))
          if(st%d%ispin == SPINORS) then
            density_laplacian(1:der%mesh%np, 2) = density_laplacian(1:der%mesh%np, 2) + &
                 ww*M_TWO*real(wf_psi_conj(1:der%mesh%np, 2)*lwf_psi(1:der%mesh%np, 2))
            density_laplacian(1:der%mesh%np, 3) = density_laplacian(1:der%mesh%np, 3) + &
                 ww*real (lwf_psi(1:der%mesh%np, 1)*wf_psi_conj(1:der%mesh%np, 2) + &
                 wf_psi(1:der%mesh%np, 1)*conjg(lwf_psi(1:der%mesh%np, 2)))
            density_laplacian(1:der%mesh%np, 4) = density_laplacian(1:der%mesh%np, 4) + &
                 ww*aimag(lwf_psi(1:der%mesh%np, 1)*wf_psi_conj(1:der%mesh%np, 2) + &
                 wf_psi(1:der%mesh%np, 1)*conjg(lwf_psi(1:der%mesh%np, 2)))
          end if
        end if
        
        do i_dim = 1, der%mesh%sb%dim

          !We precompute some quantites, to avoid to compute it many times
          psi_gpsi(1:der%mesh%np) = wf_psi_conj(1:der%mesh%np, 1)*gwf_psi(1:der%mesh%np,i_dim,1)
          abs_gwf_psi(1:der%mesh%np) = real(conjg(gwf_psi(1:der%mesh%np, i_dim, 1))*gwf_psi(1:der%mesh%np, i_dim, 1))

          if(present(density_gradient)) &
               density_gradient(1:der%mesh%np, i_dim, is) = density_gradient(1:der%mesh%np, i_dim, is) &
                      + ww*M_TWO*real(psi_gpsi(1:der%mesh%np))
          if(present(density_laplacian)) &
               density_laplacian(1:der%mesh%np, is) = density_laplacian(1:der%mesh%np, is)             &
                      + ww*M_TWO*abs_gwf_psi(1:der%mesh%np)

          if(associated(jp)) then
            if (.not.(states_are_real(st))) then
              jp(1:der%mesh%np, i_dim, is) = jp(1:der%mesh%np, i_dim, is) + &
                    ww*aimag(psi_gpsi(1:der%mesh%np)) &
                  - ww*abs_wf_psi(1:der%mesh%np)*kpoint(i_dim)
            else
              jp(1:der%mesh%np, i_dim, is) = M_ZERO
            end if
          end if

          if (associated(tau)) then
            tau (1:der%mesh%np, is)   = tau (1:der%mesh%np, is)        + &
                 ww*(abs_gwf_psi(1:der%mesh%np) + abs(kpoint(i_dim))**2*abs_wf_psi(1:der%mesh%np)  &
                     - M_TWO*aimag(psi_gpsi(1:der%mesh%np))*kpoint(i_dim))
          end if

          if(st%d%ispin == SPINORS) then
            if(present(density_gradient)) then
              density_gradient(1:der%mesh%np, i_dim, 2) = density_gradient(1:der%mesh%np, i_dim, 2) + &
                   ww*M_TWO*real(wf_psi_conj(1:der%mesh%np, 2)*gwf_psi(1:der%mesh%np, i_dim, 2))
              density_gradient(1:der%mesh%np, i_dim, 3) = density_gradient(1:der%mesh%np, i_dim, 3) + ww* &
                   real (gwf_psi(1:der%mesh%np, i_dim, 1)*wf_psi_conj(1:der%mesh%np, 2) + &
                   wf_psi(1:der%mesh%np, 1)*conjg(gwf_psi(1:der%mesh%np, i_dim, 2)))
              density_gradient(1:der%mesh%np, i_dim, 4) = density_gradient(1:der%mesh%np, i_dim, 4) + ww* &
                   aimag(gwf_psi(1:der%mesh%np, i_dim, 1)*wf_psi_conj(1:der%mesh%np, 2) + &
                   wf_psi(1:der%mesh%np, 1)*conjg(gwf_psi(1:der%mesh%np, i_dim, 2)))
            end if

            if(present(density_laplacian)) then
              density_laplacian(1:der%mesh%np, 2) = density_laplacian(1:der%mesh%np, 2)         + &
                   ww*M_TWO*real(conjg(gwf_psi(1:der%mesh%np, i_dim, 2))*gwf_psi(1:der%mesh%np, i_dim, 2))
              density_laplacian(1:der%mesh%np, 3) = density_laplacian(1:der%mesh%np, 3)         + &
                   ww*M_TWO*real (gwf_psi(1:der%mesh%np, i_dim, 1)*conjg(gwf_psi(1:der%mesh%np, i_dim, 2)))
              density_laplacian(1:der%mesh%np, 4) = density_laplacian(1:der%mesh%np, 4)         + &
                   ww*M_TWO*aimag(gwf_psi(1:der%mesh%np, i_dim, 1)*conjg(gwf_psi(1:der%mesh%np, i_dim, 2)))
            end if

            ! the expression for the paramagnetic current with spinors is
            !     j = ( jp(1)             jp(3) + i jp(4) )
            !         (-jp(3) + i jp(4)   jp(2)           )
            if(associated(jp)) then
              jp(1:der%mesh%np, i_dim, 2) = jp(1:der%mesh%np, i_dim, 2) + &
                   ww*aimag(wf_psi_conj(1:der%mesh%np, 2)*gwf_psi(1:der%mesh%np, i_dim, 2))
              do ii = 1, der%mesh%np
                c_tmp = wf_psi_conj(ii, 1)*gwf_psi(ii, i_dim, 2) - wf_psi(ii, 2)*conjg(gwf_psi(ii, i_dim, 1))
                jp(ii, i_dim, 3) = jp(ii, i_dim, 3) + ww* real(c_tmp)
                jp(ii, i_dim, 4) = jp(ii, i_dim, 4) + ww*aimag(c_tmp)
              end do
            end if

            ! the expression for the paramagnetic current with spinors is
            !     t = ( tau(1)              tau(3) + i tau(4) )
            !         ( tau(3) - i tau(4)   tau(2)            )
            if(associated(tau)) then
              tau (1:der%mesh%np, 2) = tau (1:der%mesh%np, 2) + ww*abs(gwf_psi(1:der%mesh%np, i_dim, 2))**2
              do ii = 1, der%mesh%np
                c_tmp = conjg(gwf_psi(ii, i_dim, 1))*gwf_psi(ii, i_dim, 2)
                tau(ii, 3) = tau(ii, 3) + ww* real(c_tmp)
                tau(ii, 4) = tau(ii, 4) + ww*aimag(c_tmp)
              end do
            end if

            ASSERT(.not. present(gi_kinetic_energy_density))

          end if !SPINORS

        end do

      end do
    end do

    SAFE_DEALLOCATE_A(wf_psi_conj)
    SAFE_DEALLOCATE_A(abs_wf_psi)
    SAFE_DEALLOCATE_A(abs_gwf_psi)
    SAFE_DEALLOCATE_A(psi_gpsi)

    if(.not. present(gi_kinetic_energy_density)) then
      if(.not. present(paramagnetic_current)) then
        SAFE_DEALLOCATE_P(jp)
      end if
      if(.not. present(kinetic_energy_density)) then
        SAFE_DEALLOCATE_P(tau)
      end if
    end if

    if(st%parallel_in_states .or. st%d%kpt%parallel) call reduce_all(st%st_kpt_mpi_grp)

    ! We have to symmetrize everything as they are calculated from the
    ! wavefunctions.
    ! This must be done before compute the gauge-invariant kinetic energy density 
    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symm(1:der%mesh%np, 1:der%mesh%sb%dim))
      call symmetrizer_init(symmetrizer, der%mesh)
      do is = 1, st%d%nspin
        if(associated(tau)) then
          call dsymmetrizer_apply(symmetrizer, der%mesh%np, field = tau(:, is), symmfield = symm(:,1), &
            suppress_warning = .true.)
          tau(1:der%mesh%np, is) = symm(1:der%mesh%np,1)
        end if

        if(present(density_laplacian)) then
          call dsymmetrizer_apply(symmetrizer, der%mesh%np, field = density_laplacian(:, is), symmfield = symm(:,1), &
            suppress_warning = .true.)
          density_laplacian(1:der%mesh%np, is) = symm(1:der%mesh%np,1)
        end if

        if(associated(jp)) then 
          call dsymmetrizer_apply(symmetrizer, der%mesh%np, field_vector = jp(:, :, is), symmfield_vector = symm, &
            suppress_warning = .true.)
          jp(1:der%mesh%np, 1:der%mesh%sb%dim, is) = symm(1:der%mesh%np, 1:der%mesh%sb%dim)
        end if
 
        if(present(density_gradient)) then
          call dsymmetrizer_apply(symmetrizer, der%mesh%np, field_vector = density_gradient(:, :, is), &
            symmfield_vector = symm, suppress_warning = .true.)
          density_gradient(1:der%mesh%np, 1:der%mesh%sb%dim, is) = symm(1:der%mesh%np, 1:der%mesh%sb%dim)
        end if   
      end do
      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symm) 
    end if


    if(associated(st%rho_core) .and. nlcc .and. (present(density_laplacian) .or. present(density_gradient))) then
       forall(ii=1:der%mesh%np)
         wf_psi(ii, 1) = st%rho_core(ii)/st%d%spin_channels
       end forall

       call boundaries_set(der%boundaries, wf_psi(:, 1))

       if(present(density_gradient)) then
         ! calculate gradient of the NLCC
         call zderivatives_grad(der, wf_psi(:,1), gwf_psi(:,:,1), set_bc = .false.)
         do is = 1, st%d%spin_channels
           density_gradient(1:der%mesh%np, 1:der%mesh%sb%dim, is) = density_gradient(1:der%mesh%np, 1:der%mesh%sb%dim, is) + &
                                                                    gwf_psi(1:der%mesh%np, 1:der%mesh%sb%dim,1)
         end do
       end if

       ! calculate the Laplacian of the wavefunction
       if (present(density_laplacian)) then
         call zderivatives_lapl(der, wf_psi(:,1), lwf_psi(:,1), set_bc = .false.)

         do is = 1, st%d%spin_channels
           density_laplacian(1:der%mesh%np, is) = density_laplacian(1:der%mesh%np, is) + lwf_psi(1:der%mesh%np, 1)
         end do
      end if
    end if

    SAFE_DEALLOCATE_A(wf_psi)
    SAFE_DEALLOCATE_A(gwf_psi)
    SAFE_DEALLOCATE_A(lwf_psi)


    !We compute the gauge-invariant kinetic energy density
    if(present(gi_kinetic_energy_density) .and. st%d%ispin /= SPINORS) then
      do is = 1, st%d%nspin
        ASSERT(associated(tau))
        gi_kinetic_energy_density(1:der%mesh%np, is) = tau(1:der%mesh%np, is)
        if(states_are_complex(st)) then
          ASSERT(associated(jp))
          do ii = 1, der%mesh%np
            if(st%rho(ii, is) < CNST(1.0e-7)) cycle
            gi_kinetic_energy_density(ii, is) = &
              gi_kinetic_energy_density(ii, is) - sum(jp(ii,1:der%mesh%sb%dim, is)**2)/st%rho(ii, is)
          end do
        end if
      end do
    end if

    if(.not. present(kinetic_energy_density)) then
      SAFE_DEALLOCATE_P(tau)
    end if
    if(.not. present(paramagnetic_current)) then
      SAFE_DEALLOCATE_P(jp)
    end if


    POP_SUB(states_calc_quantities)

    call profiling_out(prof)

  contains

    subroutine reduce_all(grp)
      type(mpi_grp_t), intent(in)  :: grp

      PUSH_SUB(states_calc_quantities.reduce_all)

      if(associated(tau)) call comm_allreduce(grp%comm, tau, dim = (/der%mesh%np, st%d%nspin/))

      if (present(density_laplacian)) call comm_allreduce(grp%comm, density_laplacian, dim = (/der%mesh%np, st%d%nspin/))

      do is = 1, st%d%nspin
        if(associated(jp)) call comm_allreduce(grp%comm, jp(:, :, is), dim = (/der%mesh%np, der%mesh%sb%dim/))

        if(present(density_gradient)) &
          call comm_allreduce(grp%comm, density_gradient(:, :, is), dim = (/der%mesh%np, der%mesh%sb%dim/))
      end do

      POP_SUB(states_calc_quantities.reduce_all)
    end subroutine reduce_all

  end subroutine states_calc_quantities


  ! ---------------------------------------------------------
  function state_spin(mesh, f1) result(spin)
    type(mesh_t), intent(in) :: mesh
    CMPLX,        intent(in) :: f1(:, :)
    FLOAT                    :: spin(1:3)

    CMPLX :: z

    PUSH_SUB(state_spin)

    z = zmf_dotp(mesh, f1(:, 1) , f1(:, 2))

    spin(1) = M_TWO*dble(z)
    spin(2) = M_TWO*aimag(z)
    spin(3) = zmf_nrm2(mesh, f1(:, 1))**2 - zmf_nrm2(mesh, f1(:, 2))**2
    spin = M_HALF*spin ! spin is half the sigma matrix.

    POP_SUB(state_spin)
  end function state_spin

  ! ---------------------------------------------------------
  logical function state_is_local(st, ist)
    type(states_t), intent(in) :: st
    integer,        intent(in) :: ist

    PUSH_SUB(state_is_local)

    state_is_local = ist >= st%st_start.and.ist <= st%st_end

    POP_SUB(state_is_local)
  end function state_is_local

  ! ---------------------------------------------------------
  logical function state_kpt_is_local(st, ist, ik)
    type(states_t), intent(in) :: st
    integer,        intent(in) :: ist
    integer,        intent(in) :: ik

    PUSH_SUB(state_kpt_is_local)

    state_kpt_is_local = ist >= st%st_start .and. ist <= st%st_end .and. &
      ik >= st%d%kpt%start .and. ik <= st%d%kpt%end

    POP_SUB(state_kpt_is_local)
  end function state_kpt_is_local


  ! ---------------------------------------------------------

  real(8) function states_wfns_memory(st, mesh) result(memory)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: mesh

    PUSH_SUB(states_wfns_memory)
    memory = 0.0_8

    ! orbitals
    memory = memory + REAL_PRECISION*dble(mesh%np_part_global)*st%d%dim*dble(st%nst)*st%d%kpt%nglobal

    POP_SUB(states_wfns_memory)
  end function states_wfns_memory

  ! ---------------------------------------------------------

  subroutine states_pack(st, copy)
    type(states_t),    intent(inout) :: st
    logical, optional, intent(in)    :: copy

    integer :: iqn, ib
    integer(8) :: max_mem, mem

    PUSH_SUB(states_pack)

    ! nothing to do, already packed
    if (st%packed) then
      POP_SUB(states_pack)
      return
    end if

    st%packed = .true.

    if(accel_is_enabled()) then
      max_mem = accel_global_memory_size()
      
      if(st%d%cl_states_mem > CNST(1.0)) then
        max_mem = int(st%d%cl_states_mem, 8)*(1024_8)**2
      else if(st%d%cl_states_mem < CNST(0.0)) then
        max_mem = max_mem + int(st%d%cl_states_mem, 8)*(1024_8)**2
      else
        max_mem = int(st%d%cl_states_mem*real(max_mem, REAL_PRECISION), 8)
      end if
    else
      max_mem = HUGE(max_mem)
    end if

    mem = 0
    qnloop: do iqn = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end

        mem = mem + batch_pack_size(st%group%psib(ib, iqn))

        if(mem > max_mem) then
          call messages_write('Not enough CL device memory to store all states simultaneously.', new_line = .true.)
          call messages_write('Only ')
          call messages_write(ib - st%group%block_start)
          call messages_write(' of ')
          call messages_write(st%group%block_end - st%group%block_start + 1)
          call messages_write(' blocks will be stored in device memory.', new_line = .true.)
          call messages_warning()
          exit qnloop
        end if
        
        call batch_pack(st%group%psib(ib, iqn), copy)
      end do
    end do qnloop

    POP_SUB(states_pack)
  end subroutine states_pack

  ! ------------------------------------------------------------

  subroutine states_unpack(st, copy)
    type(states_t),    intent(inout) :: st
    logical, optional, intent(in)    :: copy

    integer :: iqn, ib

    PUSH_SUB(states_unpack)

    if(st%packed) then
      st%packed = .false.

      do iqn = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          if(batch_is_packed(st%group%psib(ib, iqn))) call batch_unpack(st%group%psib(ib, iqn), copy)
        end do
      end do
    end if

    POP_SUB(states_unpack)
  end subroutine states_unpack

  ! -----------------------------------------------------------

  subroutine states_write_info(st)
    type(states_t),    intent(in) :: st

    PUSH_SUB(states_write_info)

    call messages_print_stress(stdout, "States")

    write(message(1), '(a,f12.3)') 'Total electronic charge  = ', st%qtot
    write(message(2), '(a,i8)')    'Number of states         = ', st%nst
    write(message(3), '(a,i8)')    'States block-size        = ', st%d%block_size
    call messages_info(3)

    call messages_print_stress(stdout)

    POP_SUB(states_write_info)
  end subroutine states_write_info
 
  ! -----------------------------------------------------------

  logical pure function states_are_packed(st) result(packed)
    type(states_t),    intent(in) :: st

    packed = st%packed
  end function states_are_packed

  ! ------------------------------------------------------------

  subroutine states_set_zero(st)
    type(states_t),    intent(inout) :: st

    integer :: iqn, ib

    PUSH_SUB(states_set_zero)

    do iqn = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call batch_set_zero(st%group%psib(ib, iqn))
      end do
    end do
    
    POP_SUB(states_set_zero)
  end subroutine states_set_zero

  ! ------------------------------------------------------------

  integer pure function states_block_min(st, ib) result(range)
    type(states_t),    intent(in) :: st
    integer,           intent(in) :: ib
    
    range = st%group%block_range(ib, 1)
  end function states_block_min

  ! ------------------------------------------------------------

  integer pure function states_block_max(st, ib) result(range)
    type(states_t),    intent(in) :: st
    integer,           intent(in) :: ib
    
    range = st%group%block_range(ib, 2)
  end function states_block_max

  ! ------------------------------------------------------------

  integer pure function states_block_size(st, ib) result(size)
    type(states_t),    intent(in) :: st
    integer,           intent(in) :: ib
    
    size = st%group%block_size(ib)
  end function states_block_size

  ! ---------------------------------------------------------
  !> number of occupied-unoccipied pairs for Casida
  subroutine states_count_pairs(st, parser, n_pairs, n_occ, n_unocc, is_included, is_frac_occ)
    type(states_t),    intent(in)  :: st
    type(parser_t),    intent(in)  :: parser
    integer,           intent(out) :: n_pairs
    integer,           intent(out) :: n_occ(:)   !< nik
    integer,           intent(out) :: n_unocc(:) !< nik
    logical, pointer,  intent(out) :: is_included(:,:,:) !< (max(n_occ), max(n_unocc), st%d%nik)
    logical,           intent(out) :: is_frac_occ !< are there fractional occupations?

    integer :: ik, ist, ast, n_filled, n_partially_filled, n_half_filled
    character(len=80) :: nst_string, default, wfn_list
    FLOAT :: energy_window

    PUSH_SUB(states_count_pairs)

    is_frac_occ = .false.
    do ik = 1, st%d%nik
      call occupied_states(st, ik, n_filled, n_partially_filled, n_half_filled)
      if(n_partially_filled > 0 .or. n_half_filled > 0) is_frac_occ = .true.
      n_occ(ik) = n_filled + n_partially_filled + n_half_filled
      n_unocc(ik) = st%nst - n_filled
      ! when we implement occupations, partially occupied levels need to be counted as both occ and unocc.
    end do

    !%Variable CasidaKSEnergyWindow
    !%Type float
    !%Section Linear Response::Casida
    !%Description
    !% An alternative to <tt>CasidaKohnShamStates</tt> for specifying which occupied-unoccupied
    !% transitions will be used: all those whose eigenvalue differences are less than this
    !% number will be included. If a value less than 0 is supplied, this criterion will not be used.
    !%End

    call parse_variable(parser, 'CasidaKSEnergyWindow', -M_ONE, energy_window, units_inp%energy)

    !%Variable CasidaKohnShamStates
    !%Type string
    !%Section Linear Response::Casida
    !%Default all states
    !%Description
    !% The calculation of the excitation spectrum of a system in the Casida frequency-domain
    !% formulation of linear-response time-dependent density functional theory (TDDFT)
    !% implies the use of a basis set of occupied/unoccupied Kohn-Sham orbitals. This
    !% basis set should, in principle, include all pairs formed by all occupied states,
    !% and an infinite number of unoccupied states. In practice, one has to truncate this
    !% basis set, selecting a number of occupied and unoccupied states that will form the
    !% pairs. These states are specified with this variable. If there are, say, 15 occupied
    !% states, and one sets this variable to the value "10-18", this means that occupied
    !% states from 10 to 15, and unoccupied states from 16 to 18 will be considered.
    !%
    !% This variable is a string in list form, <i>i.e.</i> expressions such as "1,2-5,8-15" are
    !% valid. You should include a non-zero number of unoccupied states and a non-zero number
    !% of occupied states.
    !%End

    n_pairs = 0
    SAFE_ALLOCATE(is_included(maxval(n_occ), minval(n_occ) + 1:st%nst , st%d%nik))
    is_included(:,:,:) = .false.

    if(energy_window < M_ZERO) then
      write(nst_string,'(i6)') st%nst
      write(default,'(a,a)') "1-", trim(adjustl(nst_string))
      call parse_variable(parser, 'CasidaKohnShamStates', default, wfn_list)

      write(message(1),'(a,a)') "Info: States that form the basis: ", trim(wfn_list)
      call messages_info(1)

      ! count pairs
      n_pairs = 0
      do ik = 1, st%d%nik
        do ast = n_occ(ik) + 1, st%nst
          if(loct_isinstringlist(ast, wfn_list)) then
            do ist = 1, n_occ(ik)
              if(loct_isinstringlist(ist, wfn_list)) then
                n_pairs = n_pairs + 1
                is_included(ist, ast, ik) = .true.
              end if
            end do
          end if
        end do
      end do

    else ! using CasidaKSEnergyWindow

      write(message(1),'(a,f12.6,a)') "Info: including transitions with energy < ", &
        units_from_atomic(units_out%energy, energy_window), trim(units_abbrev(units_out%energy))
      call messages_info(1)

      ! count pairs
      n_pairs = 0
      do ik = 1, st%d%nik
        do ast = n_occ(ik) + 1, st%nst
          do ist = 1, n_occ(ik)
            if(st%eigenval(ast, ik) - st%eigenval(ist, ik) < energy_window) then
              n_pairs = n_pairs + 1
              is_included(ist, ast, ik) = .true.
            end if
          end do
        end do
      end do

    end if

    POP_SUB(states_count_pairs)
  end subroutine states_count_pairs

  ! ---------------------------------------------------------
  !> Returns information about which single-particle orbitals are
  !! occupied or not in a _many-particle_ state st:
  !!   n_filled are the number of orbitals that are totally filled
  !!            (the occupation number is two, if ispin = UNPOLARIZED,
  !!            or it is one in the other cases).
  !!   n_half_filled is only meaningful if ispin = UNPOLARIZED. It 
  !!            is the number of orbitals where there is only one 
  !!            electron in the orbital.
  !!   n_partially_filled is the number of orbitals that are neither filled,
  !!            half-filled, nor empty.
  !! The integer arrays filled, partially_filled and half_filled point
  !!   to the indices where the filled, partially filled and half_filled
  !!   orbitals are, respectively.
  subroutine occupied_states(st, ik, n_filled, n_partially_filled, n_half_filled, &
                             filled, partially_filled, half_filled)
    type(states_t),    intent(in)  :: st
    integer,           intent(in)  :: ik
    integer,           intent(out) :: n_filled, n_partially_filled, n_half_filled
    integer, optional, intent(out) :: filled(:), partially_filled(:), half_filled(:)

    integer :: ist
    FLOAT, parameter :: M_THRESHOLD = CNST(1.0e-6)

    PUSH_SUB(occupied_states)

    if(present(filled))           filled(:) = 0
    if(present(partially_filled)) partially_filled(:) = 0
    if(present(half_filled))      half_filled(:) = 0
    n_filled = 0
    n_partially_filled = 0
    n_half_filled = 0

    select case(st%d%ispin)
    case(UNPOLARIZED)
      do ist = 1, st%nst
        if(abs(st%occ(ist, ik) - M_TWO) < M_THRESHOLD) then
          n_filled = n_filled + 1
          if(present(filled)) filled(n_filled) = ist
        elseif(abs(st%occ(ist, ik) - M_ONE) < M_THRESHOLD) then
          n_half_filled = n_half_filled + 1
          if(present(half_filled)) half_filled(n_half_filled) = ist
        elseif(st%occ(ist, ik) > M_THRESHOLD ) then
          n_partially_filled = n_partially_filled + 1
          if(present(partially_filled)) partially_filled(n_partially_filled) = ist
        elseif(abs(st%occ(ist, ik)) > M_THRESHOLD ) then
          write(message(1),*) 'Internal error in occupied_states: Illegal occupation value ', st%occ(ist, ik)
          call messages_fatal(1)
         end if
      end do
    case(SPIN_POLARIZED, SPINORS)
      do ist = 1, st%nst
        if(abs(st%occ(ist, ik)-M_ONE) < M_THRESHOLD) then
          n_filled = n_filled + 1
          if(present(filled)) filled(n_filled) = ist
        elseif(st%occ(ist, ik) > M_THRESHOLD ) then
          n_partially_filled = n_partially_filled + 1
          if(present(partially_filled)) partially_filled(n_partially_filled) = ist
        elseif(abs(st%occ(ist, ik)) > M_THRESHOLD ) then
          write(message(1),*) 'Internal error in occupied_states: Illegal occupation value ', st%occ(ist, ik)
          call messages_fatal(1)
         end if
      end do
    end select

    POP_SUB(occupied_states)
  end subroutine occupied_states


  ! ------------------------------------------------------------
subroutine states_set_phase(st_d, psi, phase, np, conjugate)
  type(states_dim_t),intent(in)    :: st_d
  CMPLX,          intent(inout)    :: psi(:, :)
  CMPLX,             intent(in)    :: phase(:)
  integer,           intent(in)    :: np
  logical,           intent(in)    :: conjugate

  integer :: idim, ip

  PUSH_SUB(states_set_phase)

  if(conjugate) then
    ! Apply the phase that contains both the k-point and vector-potential terms.
    do idim = 1, st_d%dim
      !$omp parallel do
      do ip = 1, np
        psi(ip, idim) = conjg(phase(ip))*psi(ip, idim)
      end do
      !$omp end parallel do
    end do
  else
    ! Apply the conjugate of the phase that contains both the k-point and vector-potential terms.
    do idim = 1, st_d%dim
      !$omp parallel do
      do ip = 1, np
        psi(ip, idim) = phase(ip)*psi(ip, idim)
      end do
      !$omp end parallel do
    end do
  end if

  POP_SUB(states_set_phase)

end subroutine  states_set_phase

  
#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
