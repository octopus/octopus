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

module states_m
  use crystal_m
  use blas_m
  use datasets_m
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use hardware_m
  use io_function_m
  use io_m
  use lalg_basic_m
  use loct_m
  use loct_parser_m
  use math_m
  use messages_m
  use mesh_function_m
  use mesh_m
  use mpi_m
  use mpi_lib_m
  use multicomm_m
  use profiling_m
  use simul_box_m
  use string_m
  use units_m
  use varinfo_m

  implicit none

  private
  public ::                           &
    states_t,                         &
    states_dim_t,                     &
    states_init,                      &
    states_look,                      &
    states_read_user_def_orbitals,    &
    states_densities_init,            &
    states_allocate_wfns,             &
    states_deallocate_wfns,           &
    states_null,                      &
    states_end,                       &
    states_copy,                      &
    states_dim_copy,                  &
    states_dim_end,                   &
    states_generate_random,           &
    states_orthogonalize,             &
    states_fermi,                     &
    states_eigenvalues_sum,           &
    states_write_eigenvalues,         &
    states_write_dos,                 &
    states_write_bands,               &
    states_write_fermi_energy,        &
    states_degeneracy_matrix,         &
    states_spin_channel,              &
    states_dens_accumulate,           &
    states_dens_reduce,               &
    states_calc_dens,                 &
    states_calc_tau_jp_gn,            &
    state_is_local,                   &
    kpoints_write_info,               &
    kpoint_is_gamma,                  &
    wfs_are_complex,                  &
    wfs_are_real,                     &
    states_dump,                      &
    rotate_states,                    &
    states_freeze_orbitals,           &
    assignment(=)


  public ::                         &
    dstates_gram_schmidt,           &
    zstates_gram_schmidt,           &
    dstates_gram_schmidt_full,      &
    zstates_gram_schmidt_full,      &
    dstates_dotp,                   &
    zstates_dotp,                   &
    dstates_nrm2,                   &
    zstates_nrm2,                   &
    dstates_normalize_orbital,      &
    zstates_normalize_orbital,      &
    dstates_residue,                &
    zstates_residue,                &
    dstates_calc_momentum,          &
    zstates_calc_momentum,          &
    dstates_angular_momentum,       &
    zstates_angular_momentum,       &
    dstates_matrix,                 &
    zstates_matrix,                 &
    dstates_linear_combination,     &
    zstates_linear_combination,     &
    states_distribute_nodes

  type states_dim_t
    integer :: dim                  ! Dimension of the state (one or two for spinors)
    integer :: nik                  ! Number of irreducible subspaces
    integer :: nik_axis(MAX_DIM)    ! Number of kpoints per axis
    integer :: ispin                ! spin mode (unpolarized, spin polarized, spinors)
    integer :: nspin                ! dimension of rho (1, 2 or 4)
    integer :: spin_channels        ! 1 or 2, wether spin is or not considered.
    logical :: cdft                 ! Are we using Current-DFT or not?
    FLOAT, pointer :: kpoints(:,:)  ! obviously the kpoints
    FLOAT, pointer :: kweights(:)   ! weights for the kpoint integrations
  end type states_dim_t

  type states_t
    type(states_dim_t) :: d
    integer :: nst                  ! Number of states in each irreducible subspace

    integer :: wfs_type             ! real (M_REAL) or complex (M_CMPLX) wavefunctions
    ! pointers to the wavefunctions 
    logical :: only_userdef_istates ! only use user defined states initial states in propagation
    FLOAT, pointer :: dpsi(:,:,:,:) ! dpsi(sys%NP_PART, st%d%dim, st%nst, st%d%nik)
    CMPLX, pointer :: zpsi(:,:,:,:) ! zpsi(sys%NP_PART, st%d%dim, st%nst, st%d%nik)

    ! used for the user defined wavefunctions (they are stored as formula strings)
    character(len=1024), pointer :: user_def_states(:,:,:) ! (st%d%dim, st%nst, st%d%nik)

    ! the densities and currents (after all we are doing DFT :)
    FLOAT, pointer :: rho(:,:)      ! rho(gr%m%np_part, st%d%nspin)
    FLOAT, pointer :: j(:,:,:)      !   j(gr%m%np_part, gr%sb%dim, st%d%nspin)

    logical        :: nlcc          ! do we have non-linear core corrections
    FLOAT, pointer :: rho_core(:)   ! core charge for nl core corrections

    FLOAT, pointer :: eigenval(:,:) ! obviously the eigenvalues
    logical        :: fixed_occ     ! should the occupation numbers be fixed?
    FLOAT, pointer :: occ(:,:)      ! the occupation numbers
    logical        :: fixed_spins   ! In spinors mode, the spin direction is set
                                    ! for the initial (random) orbitals.
    FLOAT, pointer :: spin(:, :, :)
    FLOAT, pointer :: momentum(:, :, :)

    FLOAT :: qtot                   ! (-) The total charge in the system (used in Fermi)
    FLOAT :: val_charge             ! valence charge

    FLOAT :: el_temp                ! electronic temperature for the Fermi function
    FLOAT :: ef                     ! the fermi energy

    ! This is stuff needed for the parallelization in states.
    logical                     :: parallel_in_states ! Am I parallel in states?
    type(mpi_grp_t)             :: mpi_grp            ! The MPI group related to the parallelization in states.
    type(mpi_grp_t)             :: dom_st             ! The MPI group related to the domain-states "plane".
    integer                     :: lnst               ! Number of states on local node.
    integer                     :: st_start, st_end   ! Range of states processed by local node.
    integer, pointer            :: node(:)            ! To which node belongs each state.
    integer, pointer            :: st_range(:, :)     ! Node r manages states st_range(1, r) to
                                                      ! st_range(2, r) for r = 0, ..., mpi_grp%size-1,
                                                      ! i. e. st_start = st_range(1, r) and
                                                      ! st_end = st_range(2, r) on node r.
    integer, pointer            :: st_num(:)          ! Number of states on node r, i. e.
                                                      ! st_num(r) = st_num(2, r)-st_num(1, r).
    type(multicomm_all_pairs_t) :: ap                 ! All-pairs schedule.
  end type states_t



  ! Parameters...
  integer, public, parameter ::     &
    UNPOLARIZED    = 1,             &
    SPIN_POLARIZED = 2,             &
    SPINORS        = 3

  interface assignment (=)
    module procedure states_copy
  end interface

contains

  ! ---------------------------------------------------------
  subroutine states_null(st)
    type(states_t), intent(inout) :: st
    call push_sub('states.states_null')

    nullify(st%dpsi, st%zpsi, st%rho, st%j, st%rho_core, st%eigenval)
    nullify(st%occ, st%spin, st%momentum, st%node, st%user_def_states)
    nullify(st%d%kpoints, st%d%kweights)
    nullify(st%st_range, st%st_num)

    ! By default, calculations use real wave-functions
    st%wfs_type = M_REAL

    call pop_sub()
  end subroutine states_null


  ! ---------------------------------------------------------
  subroutine states_init(st, gr, geo)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo

    FLOAT :: excess_charge
    integer :: nempty

    call push_sub('states.states_init')

    call states_null(st)

    !%Variable SpinComponents
    !%Type integer
    !%Default unpolarized
    !%Section States
    !%Description
    !% The calculations may be done in three different ways: spin-restricted (TD)DFT (i.e., doubly
    !% occupied "closed shells"), spin-unsrestricted or "spin-polarized" (TD)DFT (i.e. we have two
    !% electronic systes, one with spin up and one with spin down), or making use of two-component
    !% spinors.
    !%Option unpolarized 1
    !% Spin-restricted calculations.
    !%Option polarized 2
    !%Option spin_polarized 2
    !% Spin unrestricted, also know as spin-DFT, SDFT. This mode will double the number of wave
    !% functions necessary for a spin-unpolarised calculation.
    !%Option non_collinear 3
    !%Option spinors 3
    !% The spin-orbitals are two-component spinors. This effectively allows the spin-density to
    !% arrange non-collinearly - i.e. the magnetization vector is allowed to take different
    !% directions in different points.
    !%End
    call loct_parse_int(check_inp('SpinComponents'), UNPOLARIZED, st%d%ispin)
    if(.not.varinfo_valid_option('SpinComponents', st%d%ispin)) call input_error('SpinComponents')
    call messages_print_var_option(stdout, 'SpinComponents', st%d%ispin)
    ! Use of Spinors requires complex wave-functions
    if (st%d%ispin == SPINORS) st%wfs_type = M_CMPLX


    !%Variable ExcessCharge
    !%Type float
    !%Default 0.0
    !%Section States
    !%Description
    !% The net charge of the system. A negative value means that we are adding 
    !% electrons, while a positive value means we are taking electrons
    !% from the system.
    !%End
    call loct_parse_float(check_inp('ExcessCharge'), M_ZERO, excess_charge)


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
    !% However, one may command <tt>octopus</tt> to put more states, which is necessary if one wants to
    !% use fractional occupational numbers, either fixed from the origin through
    !% the <tt>Occupations</tt> block or by prescribing
    !% an electronic temperature with <tt>ElectronicTemperature</tt>.
    !%
    !% Note that this number is unrelated to <tt>CalculationMode == unocc</tt>.
    !%End
    call loct_parse_int(check_inp('ExtraStates'), 0, nempty)
    if (nempty < 0) then
      write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid ExtraStates"
      message(2) = '(0 <= ExtraStates)'
      call write_fatal(2)
    end if

    call geometry_val_charge(geo, st%val_charge)
    st%qtot = -(st%val_charge + excess_charge)

    select case(st%d%ispin)
    case(UNPOLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty
      st%d%nspin = 1
      st%d%spin_channels = 1
    case(SPIN_POLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty
      st%d%nik = st%d%nik*2
      st%d%nspin = 2
      st%d%spin_channels = 2
    case(SPINORS)
      st%d%dim = 2
      st%nst = int(st%qtot)
      if(st%nst < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty
      st%d%nspin = 4
      st%d%spin_channels = 2
    end select

    ! current
    call loct_parse_logical(check_inp('CurrentDFT'), .false., st%d%cdft)
    if (st%d%cdft) then
      ! Use of CDFT requires complex wave-functions
      st%wfs_type = M_CMPLX

      if(st%d%ispin == SPINORS) then
        message(1) = "Sorry, Current DFT not working yet for spinors"
        call write_fatal(1)
      end if
      message(1) = "Info: Using Current DFT"
      call write_info(1)
    end if

    ! For non-periodic systems this should just return the Gamma point
    call states_choose_kpoints(st%d, gr%sb, geo)

    ! Periodic systems require complex wave-functions
    if(simul_box_is_periodic(gr%sb)) st%wfs_type = M_CMPLX

    ! Transport calculations require complex wave-functions.
    if(calc_mode.eq.M_TD_TRANSPORT) then
      st%wfs_type = M_CMPLX
    end if

    !%Variable OnlyUserDefinedInitialStates
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% If true, then only user defined states from the block UserDefinedStates
    !% will be used as initial states for a time propagation. No attempt is made
    !% to load ground state orbitals from a previous ground state run.
    !%End
    call loct_parse_logical(check_inp('OnlyUserDefinedInitialStates'), .false., st%only_userdef_istates)

    ! we now allocate some arrays
    ALLOCATE(st%occ     (st%nst, st%d%nik),      st%nst*st%d%nik)
    ALLOCATE(st%eigenval(st%nst, st%d%nik),      st%nst*st%d%nik)
    ALLOCATE(st%momentum(3, st%nst, st%d%nik), 3*st%nst*st%d%nik)
    st%occ = M_ZERO
    st%eigenval = M_ZERO
    st%momentum = M_ZERO
    ! allocate space for formula strings that define user defined states
    ALLOCATE(st%user_def_states(st%d%dim, st%nst, st%d%nik), st%d%dim*st%nst*st%d%nik)
    if(st%d%ispin == SPINORS) then
      ALLOCATE(st%spin(3, st%nst, st%d%nik), st%nst*st%d%nik*3)
    else
      nullify(st%spin)
    end if

    ! initially we mark all 'formulas' as undefined
    st%user_def_states(1:st%d%dim, 1:st%nst, 1:st%d%nik) = 'undefined'

    call states_read_initial_occs(st, excess_charge)
    call states_read_initial_spins(st)

    st%st_start = 1
    st%st_end = st%nst
    st%lnst = st%nst
    ALLOCATE(st%node(st%nst), st%nst)
    st%node(1:st%nst) = 0

    call mpi_grp_init(st%mpi_grp, -1)
    st%parallel_in_states = .false.

    nullify(st%dpsi, st%zpsi)

    call pop_sub()
  end subroutine states_init


  ! ---------------------------------------------------------
  ! Reads from the input file the initial occupations, if the
  ! block "Occupations" is present. Otherwise, it makes an initial
  ! guess for the occupations, maybe using the "ElectronicTemperature"
  ! variable.
  ! The resulting occupations are placed on the st%occ variable. The
  ! boolean st%fixed_occ is also set to .true., if the occupations are
  ! set by the user through the "Occupations" block; false otherwise.
  subroutine states_read_initial_occs(st, excess_charge)
    type(states_t), intent(inout) :: st
    FLOAT, intent(in) :: excess_charge

    integer :: i, j, ncols
    type(block_t) :: blk
    FLOAT :: r

    call push_sub('states.states_read_initial_occs')
    !%Variable Occupations
    !%Type block
    !%Section States
    !%Description
    !% The occupation numbers of the orbitals can be fixed through the use of this
    !% variable. For example:
    !%
    !% <tt>%Occupations
    !% <br>&nbsp;&nbsp;2.0 | 2.0 | 2.0 | 2.0 | 2.0
    !% <br>%</tt>
    !%
    !% would fix the occupations of the five states to <i>2.0</i>. There must be,
    !% at most, as many columns as states in the calculation. If there are less columns
    !% than states, then the code will assume that the user is indicating the occupations
    !% of the uppermost states, assigning maximum occupation (i.e. 2 for spin-unpolarized
    !% calculations, 1 otherwise) to the lower states. If <tt>SpinComponents == polarized</tt>
    !% this block should contain two lines, one for each spin channel.
    !% This variable is very useful when dealing with highly symmetric small systems
    !% (like an open shell atom), for it allows us to fix the occupation numbers
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
    !% <br>&nbsp;&nbsp; 1/3 | 1/3 | 1/3
    !% <br>&nbsp;&nbsp; 0   |   0 |   0
    !% <br>%</tt>
    !%
    !% Note that in this case the first state is absent; the code will calculate four states
    !% (two because there are four electrons, plus two because ExtraStates = 2), and since
    !% it finds only three columns, it will occupy the first state with one electron for each
    !% of the spin options.
    !%End



    occ_fix: if(loct_parse_block(check_inp('Occupations'), blk)==0) then
      ! read in occupations
      st%fixed_occ = .true.

      ! Reads the number of columns in the first row. This assumes that all rows
      ! have the same column number; otherwise the code will stop with an error.
      ncols = loct_parse_block_cols(blk, 0)
      if(ncols > st%nst) then
        call input_error("Occupations")
      end if
      ! Now we fill al the "missing" states with the maximum occupation.
      do i = 1, st%d%nik
        do j = 1, st%nst - ncols
          if(st%d%ispin == UNPOLARIZED) then
            st%occ(j, i) = M_TWO
          else
            st%occ(j, i) = M_ONE
          end if
        end do
      end do
      do i = 1, st%d%nik
        do j = st%nst - ncols + 1, st%nst 
          call loct_parse_block_float(blk, i-1, j-1-(st%nst-ncols), st%occ(j, i))
        end do
      end do
      call loct_parse_block_end(blk)

    else
      st%fixed_occ = .false.

      ! first guess for occupation...paramagnetic configuration
      if(st%d%ispin == UNPOLARIZED) then
        r = M_TWO
      else
        r = M_ONE
      end if
      st%occ  = M_ZERO
      st%qtot = M_ZERO

      do j = 1, st%nst
        do i = 1, st%d%nik
          st%occ(j, i) = min(r, -(st%val_charge + excess_charge) - st%qtot)
          st%qtot = st%qtot + st%occ(j, i)

        end do
      end do

      !%Variable ElectronicTemperature
      !%Type float
      !%Default 0.0
      !%Section States
      !%Description
      !% If <tt>Occupations</tt> is not set, <tt>ElectronicTemperature</tt> is the
      !% temperature in the Fermi-Dirac function used to distribute the electrons
      !% among the existing states.
      !%End
      call loct_parse_float(check_inp('ElectronicTemperature'), M_ZERO, st%el_temp)
      st%el_temp = st%el_temp * units_inp%energy%factor
    end if occ_fix

    call pop_sub()
  end subroutine states_read_initial_occs


  ! ---------------------------------------------------------
  ! Reads, if present, the "InitialSpins" block. This is only
  ! done in spinors mode; otherwise the routine does nothing. The
  ! resulting spins are placed onto the st%spin pointer. The boolean
  ! st%fixed_spins is set to true if (and only if) the InitialSpins
  ! block is present.
  subroutine states_read_initial_spins(st)
    type(states_t), intent(inout) :: st
    integer :: i, j
    type(block_t) :: blk

    call push_sub('states.states_read_initial_spins')

    st%fixed_spins = .false.
    if(st%d%ispin .ne. SPINORS) then
      call pop_sub(); return
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
    !% (i.e. SpinComponents = spinors). 
    !%
    !% The structure of the block is very simple: each column contains the desired
    !% <Sx>, <Sy>, <Sz> for each spinor. If the calculation is for a periodic system
    !% and there are more than one k point, the spins of all the k-point space are
    !% the same.
    !%
    !% For example, if we have two spinors, and we want one in the Sx "down" state,
    !% and another one in the Sx "up" state:
    !%
    !% <tt>%InitialSpins
    !% <br>&nbsp;&nbsp;  0.5 | 0.0 | 0.0
    !% <br>&nbsp;&nbsp; -0.5 | 0.0 | 0.0
    !% <br>%</tt>
    !%
    !% WARNING: if the calculation is for a system described by pseudopotentials (as
    !% opposed to using user defined potentials or model systems), this option is
    !% meaningless since the random spinors are overwritten by the atomic orbitals.
    !%
    !% There are a couple of physical constrains that have to be fulfilled:
    !%
    !% (A) | <S_i> | <= 1/2
    !%
    !% (B) <S_x>^2 + <S_y>^2 + <S_z>^2 = 1/4
    !%
    !%End
    spin_fix: if(loct_parse_block(check_inp('InitialSpins'), blk)==0) then
      do i = 1, st%nst
        do j = 1, 3
          call loct_parse_block_float(blk, i-1, j-1, st%spin(j, i, 1))
        end do
        ! This checks (B).
        if( abs(sum(st%spin(1:3, i, 1)**2) - M_FOURTH) > CNST(1.0e-6)) call input_error('InitialSpins')
      end do
      call loct_parse_block_end(blk)
      ! This checks (A). In fact (A) follows from (B), so maybe this is not necessary...
      if(any(abs(st%spin(:, :, :)) > M_HALF)) then
        call input_error('InitialSpins')
      end if
      st%fixed_spins = .true.
      do i = 2, st%d%nik
        st%spin(:, :, i) = st%spin(:, :, 1)     
      end do
    end if spin_fix

    call pop_sub()
  end subroutine states_read_initial_spins


  ! ---------------------------------------------------------
  ! Allocates the KS wavefunctions defined within an states_t
  ! structure.
  subroutine states_allocate_wfns(st, m, wfs_type)
    type(states_t), intent(inout) :: st
    type(mesh_t),    intent(in)    :: m
    integer, optional, intent(in) :: wfs_type

    integer :: n, ik, ist, idim
    logical :: force

    call push_sub('states.states_allocate_wfns')

    if (present(wfs_type)) then
      ASSERT(wfs_type == M_REAL .or. wfs_type == M_CMPLX)
      st%wfs_type = wfs_type
    end if

    !%Variable ForceComplex
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% Normally Octopus determines automatically the type necessary
    !% for the wave functions. When set to yes this variable will
    !% force the use of complex wavefunctions. 
    !%
    !% Warning: This variable is designed for testing and
    !% benchmarching and normal users need not use it.
    !%
    !%End
    call loct_parse_logical(check_inp('ForceComplex'), .false., force)
    
    if(force) st%wfs_type = M_CMPLX

    n = m%np_part * st%d%dim * st%lnst * st%d%nik
    if (st%wfs_type == M_REAL) then
      ALLOCATE(st%dpsi(m%np_part, st%d%dim, st%st_start:st%st_end, st%d%nik), n)

      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            !$omp parallel workshare
            st%dpsi(1:m%np, idim, ist, ik) = M_ZERO
            st%dpsi(m%np+1:m%np_part, idim, ist, ik) = M_ZERO
            !$omp end parallel workshare
          end do
        end do
      end do

    else
      ALLOCATE(st%zpsi(m%np_part, st%d%dim, st%st_start:st%st_end, st%d%nik), n)

      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            !$omp parallel workshare
            st%zpsi(1:m%np, idim, ist, ik) = M_Z0
            st%zpsi(m%np+1:m%np_part, idim, ist, ik) = M_Z0
            !$omp end parallel workshare
          end do
        end do
      end do

    end if

    call pop_sub()
  end subroutine states_allocate_wfns


  ! ---------------------------------------------------------
  ! Deallocates the KS wavefunctions defined within an states_t
  ! structure.
  subroutine states_deallocate_wfns(st)
    type(states_t), intent(inout) :: st

    call push_sub('states.states_deallocate_wfns')

    if (st%wfs_type == M_REAL) then
      deallocate(st%dpsi); nullify(st%dpsi)
    else
      deallocate(st%zpsi); nullify(st%zpsi)
    end if

    call pop_sub()
  end subroutine states_deallocate_wfns

  ! ---------------------------------------------------------
  ! This routine transforms the orbitals of state "st", according
  ! to the transformation matrix "u".
  !
  ! Each row of u contains the coefficients of the new orbitals
  ! in terms of the old ones.
  ! ---------------------------------------------------------
  subroutine rotate_states(mesh, st, stin, u)
    type(mesh_t),      intent(in)    :: mesh
    type(states_t),    intent(inout) :: st
    type(states_t),    intent(in)    :: stin
    CMPLX,             intent(in)    :: u(:, :)

    integer :: ik

    call push_sub('states.rotate_states')

    if(st%wfs_type == M_REAL) then
      do ik = 1, st%d%nik
        call lalg_gemm(mesh%np_part*st%d%dim, st%nst, stin%nst, M_ONE, stin%dpsi(:, :, 1:stin%nst, ik), &
                       transpose(real(u(:, :), REAL_PRECISION)), M_ZERO, st%dpsi(:, :, :, ik))
      end do
    else
      do ik = 1, st%d%nik
        call lalg_gemm(mesh%np_part*st%d%dim, st%nst, stin%nst, M_z1, stin%zpsi(:, :, 1:stin%nst, ik), &
                       transpose(u(:, :)), M_z0, st%zpsi(:, :, :, ik))
      end do
    end if

    call pop_sub()
  end subroutine rotate_states


  ! ---------------------------------------------------------
  ! the routine reads formulas for user defined wavefunctions 
  ! from the input file and fills the respective orbitals
  subroutine states_read_user_def_orbitals(mesh, st)
    type(mesh_t),      intent(in) :: mesh
    type(states_t), intent(inout) :: st    

    type(block_t) :: blk
    integer   :: ip, id, is, ik, nstates, state_from, ierr, ncols
    integer   :: ib, idim, inst, inik, normalize
    FLOAT     :: x(MAX_DIM), r, psi_re, psi_im
    character(len=150) :: filename

    integer, parameter ::      &
      state_from_formula  = 1, &
      state_from_file     = 0, &
      normalize_yes       = 1, &
      normalize_no        = 0

    call push_sub('td.read_user_def_states')

    !%Variable UserDefinedStates
    !%Type block
    !%Section States
    !%Description
    !% Instead of using the ground state as initial state for
    !% time propagations it might be interesting in some cases 
    !% to specify alternative states. Similar to user defined
    !% potentials this block allows to specify formulas for
    !% the orbitals at t=0.
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedStates
    !% <br>&nbsp;&nbsp; 1 | 1 | 1 | formula | "exp(-r^2)*exp(-i*0.2*x)" | normalize_yes
    !% <br>%</tt>
    !%
    !% The first column specifies the component of the spinor, 
    !% the second column the number of the state and the third 
    !% contains kpoint and spin quantum numbers. Column four
    !% indicates that column five should be interpreted as formula
    !% for the correspondig orbital.
    !% 
    !% Alternatively, if column four states file the state will
    !% be read from the file given in column five.
    !%
    !% <tt>%UserDefinedStates
    !% <br>&nbsp;&nbsp; 1 | 1 | 1 | file | "/path/to/file" | normalize_no
    !% <br>%</tt>
    !% 
    !% Octopus reads first the ground state orbitals from
    !% the restart/gs directory. Only the states that are
    !% specified in the above block will be overwritten with
    !% the given analytical expression for the orbital.
    !%
    !% The sixth (optional) column indicates if Octopus should renormalize the orbital
    !% or not. The default, i. e. no sixth column given, is to renormalize.
    !%
    !%Option file 0
    !% Read initsial orbital from file
    !%Option formula 1
    !% Calculate initial orbital by given analytic expression
    !%Option normalize_yes 1
    !% Normalize orbitals (default)
    !%Option normalize_no 0
    !% Do not normalize orbitals
    !%End
    if(loct_parse_block(check_inp('UserDefinedStates'), blk) == 0) then

      call messages_print_stress(stdout, trim('Substitution of orbitals'))

      ! find out how many lines (i.e. states) the block has
      nstates = loct_parse_block_n(blk)

      ! read all lines
      do ib = 1, nstates
        ! Check that number of columns is five or six.
        ncols = loct_parse_block_cols(blk, ib-1)
        if(ncols.lt.5.or.ncols.gt.6) then
          message(1) = 'Each line in the UserDefinedStates block must have'
          message(2) = 'five or six columns.'
          call write_fatal(2)
        end if
        
        call loct_parse_block_int(blk, ib-1, 0, idim)
        call loct_parse_block_int(blk, ib-1, 1, inst)
        call loct_parse_block_int(blk, ib-1, 2, inik)

        ! Calculate from expression or read from file?
        call loct_parse_block_int(blk, ib-1, 3, state_from)

        ! loop over all states
        do id = 1, st%d%dim
          do is = 1, st%nst
            do ik = 1, st%d%nik

              ! does the block entry match and is this node responsible?
              if(.not.(id.eq.idim .and. is.eq.inst .and. ik.eq.inik    &
                .and. st%st_start.le.is .and. st%st_end.ge.is) ) cycle

              select case(state_from)

              case(state_from_formula)
                ! parse formula string
                call loct_parse_block_string(                            &
                  blk, ib-1, 4, st%user_def_states(id, is, ik))

                write(message(1), '(a,3i5)') 'Substituting state of orbital with k, ist, dim = ', ik, is, id
                write(message(2), '(2a)') '  with the expression:'
                write(message(3), '(2a)') '  ',trim(st%user_def_states(id, is, ik))
                call write_info(3)

                ! convert to C string
                call conv_to_C_string(st%user_def_states(id, is, ik))

                ! fill states with user defined formulas
                do ip = 1, mesh%np
                  x = mesh%x(ip, :)
                  r = sqrt(sum(x(:)**2))

                  ! parse user defined expressions
                  call loct_parse_expression(psi_re, psi_im,             &
                    x(1), x(2), x(3), r, M_ZERO, st%user_def_states(id, is, ik))
                  ! fill state
                  st%zpsi(ip, id, is, ik) = psi_re + M_zI*psi_im
                end do

              case(state_from_file)
                ! The input format can be coded in column four now. As it is
                ! not used now, we just say "file".
                ! Read the filename.
                call loct_parse_block_string(blk, ib-1, 4, filename)

                write(message(1), '(a,3i5)') 'Substituting state of orbital with k, ist, dim = ', ik, is, id
                write(message(2), '(2a)') '  with data from file:'
                write(message(3), '(2a)') '  ',trim(filename)
                call write_info(3)

                ! finally read the state
                call zinput_function(filename, mesh, st%zpsi(:, id, is, ik), ierr, .true.)

              case default
                message(1) = 'Wrong entry in UserDefinedStates, column 4.'
                message(2) = 'You may state "formula" or "file" here.'
                call write_fatal(2)
              end select

              ! normalize orbital
              if(loct_parse_block_cols(blk, ib-1).eq.6) then
                call loct_parse_block_int(blk, ib-1, 5, normalize)
              else
                normalize = 1
              end if
              select case(normalize)
              case(normalize_no)
              case(normalize_yes)
                call zstates_normalize_orbital(mesh, st%d%dim, st%zpsi(:,:, is, ik))
              case default
                message(1) = 'The sixth column in UserDefinedStates may either be'
                message(2) = '"normalize_yes" or "normalize_no"'
                call write_fatal(2)
              end select
            end do
          end do
        end do

      end do

      call loct_parse_block_end(blk)
      call messages_print_stress(stdout)

    else
      message(1) = '"UserDefinesStates" has to be specified as block.'
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine states_read_user_def_orbitals


  ! ---------------------------------------------------------
  subroutine states_densities_init(st, gr, geo)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo
    call push_sub('states.states_densities_init')

    ! allocate arrays for charge and current densities
    ALLOCATE(st%rho(NP_PART, st%d%nspin), NP_PART*st%d%nspin)
    ALLOCATE(st%j(NP_PART, NDIM, st%d%nspin), NP_PART*NDIM*st%d%nspin)
    st%rho  = M_ZERO
    st%j    = M_ZERO
    st%nlcc = geo%nlcc
    if(st%nlcc) then
      ALLOCATE(st%rho_core(gr%m%np), gr%m%np)
      st%rho_core(:) = M_ZERO
    end if

    call pop_sub()
  end subroutine states_densities_init


  ! ---------------------------------------------------------
  subroutine states_dim_copy(dout, din)
    type(states_dim_t), intent(out) :: dout
    type(states_dim_t), intent(in)  :: din
    integer :: i

    call push_sub('states.states_dim_copy')

    dout%dim            = din%dim
    dout%nik            = din%nik
    dout%nik_axis       = din%nik_axis
    dout%ispin          = din%ispin
    dout%nspin          = din%nspin
    dout%spin_channels  = din%spin_channels
    dout%cdft           = din%cdft
    if(associated(din%kpoints)) then
      i = size(din%kpoints, 1)*size(din%kpoints, 2)
      ALLOCATE(dout%kpoints(size(din%kpoints, 1), size(din%kpoints, 2)), i)
      dout%kpoints = din%kpoints
    end if
    if(associated(din%kweights)) then
      i = size(din%kweights, 1)
      ALLOCATE(dout%kweights(size(din%kweights, 1)), i)
      dout%kweights = din%kweights
    end if

    call pop_sub()
  end subroutine states_dim_copy
  ! ---------------------------------------------------------



  ! ---------------------------------------------------------
  subroutine states_dim_end(d)
    type(states_dim_t), intent(inout) :: d
    if(associated(d%kpoints)) then
      deallocate(d%kpoints); nullify(d%kpoints)
    end if

    if(associated(d%kweights)) then
      deallocate(d%kweights); nullify(d%kweights)
    end if
  end subroutine states_dim_end



  ! ---------------------------------------------------------
  subroutine states_copy(stout, stin)
    type(states_t), intent(inout) :: stout
    type(states_t), intent(in)    :: stin

    integer :: i, j, k, l

    call push_sub('states.states_copy')

    call states_null(stout)

    stout%wfs_type = stin%wfs_type
    call states_dim_copy(stout%d, stin%d)
    stout%nst        = stin%nst
    stout%qtot       = stin%qtot
    stout%val_charge = stin%val_charge
    stout%el_temp    = stin%el_temp
    stout%ef         = stin%ef
    stout%parallel_in_states = stin%parallel_in_states
    stout%lnst       = stin%lnst
    stout%st_start   = stin%st_start
    stout%st_end     = stin%st_end
    if(associated(stin%dpsi)) then
      i = size(stin%dpsi, 1)*stin%d%dim*(stin%st_end-stin%st_start+1)*stin%d%nik
      ALLOCATE(stout%dpsi(size(stin%dpsi, 1), stin%d%dim, stin%st_start:stin%st_end, stin%d%nik), i)
      do k = 1, stin%d%nik
        do j = stin%st_start, stin%st_end
          stout%dpsi(:, :, j, k) = stin%dpsi(:, :, j, k)
        end do
      end do
    end if
    if(associated(stin%zpsi)) then
      i = size(stin%zpsi, 1)*stin%d%dim*(stin%st_end-stin%st_start+1)*stin%d%nik
      ALLOCATE(stout%zpsi(size(stin%zpsi, 1), stin%d%dim, stin%st_start:stin%st_end, stin%d%nik), i)
      do k = 1, stin%d%nik
        do j = stin%st_start, stin%st_end
          stout%zpsi(:, :, j, k) = stin%zpsi(:, :, j, k)
        end do
      end do
    end if
    if(associated(stin%user_def_states)) then
      j = size(stin%user_def_states, 1)
      k = size(stin%user_def_states, 2)
      l = size(stin%user_def_states, 3)
      i = j*k*l
      ALLOCATE(stout%user_def_states(j, k, l), i)
      stout%user_def_states = stin%user_def_states
    end if
    if(associated(stin%rho)) then
      i = size(stin%rho, 1)*size(stin%rho, 2)
      ALLOCATE(stout%rho(size(stin%rho, 1), size(stin%rho, 2)), i)
      stout%rho = stin%rho
    end if
    if(associated(stin%j)) then
      i = size(stin%j, 1)*size(stin%j, 2)*size(stin%j, 3)
      ALLOCATE(stout%j(size(stin%j, 1), size(stin%j, 2), size(stin%j, 3)), i)
      do j = 1, size(stin%j, 3)
        do k = 1, size(stin%j, 2)
          stout%j(:, k, j) = stin%j(:, k, j)
        end do
      end do
    end if
    stout%nlcc = stin%nlcc
    if(associated(stin%rho_core)) then
      i = size(stin%rho_core, 1)
      ALLOCATE(stout%rho_core(size(stin%rho_core, 1)), i)
      stout%rho_core = stin%rho_core
    end if
    if(associated(stin%eigenval)) then
      i = (stin%st_end-stin%st_start)*stin%d%nik
      ALLOCATE(stout%eigenval(stin%st_start:stin%st_end, stin%d%nik), i)
      stout%eigenval = stin%eigenval
    end if
    stout%fixed_occ = stin%fixed_occ
    if(associated(stin%occ)) then
      i = size(stin%occ, 1)*size(stin%occ, 2)
      ALLOCATE(stout%occ(size(stin%occ, 1), size(stin%occ, 2)), i)
      stout%occ = stin%occ
    end if
    stout%fixed_spins = stin%fixed_spins
    if(associated(stin%spin)) then
      i = size(stin%spin, 1)*size(stin%spin, 2)*size(stin%spin, 3)
      ALLOCATE(stout%spin(size(stin%spin, 1), size(stin%spin, 2), size(stin%spin, 3)), i)
      stout%spin = stin%spin
    end if
    if(associated(stin%momentum)) then
      i = size(stin%momentum, 1)*size(stin%momentum, 2)*size(stin%momentum, 3)
      ALLOCATE(stout%momentum(size(stin%momentum, 1), size(stin%momentum, 2), size(stin%momentum, 3)), i)
      stout%momentum = stin%momentum
    end if
    if(associated(stin%node)) then
      i = size(stin%node)
      ALLOCATE(stout%node(size(stin%node)), i)
      stout%node = stin%node
    end if
    call mpi_grp_copy(stout%mpi_grp, stin%mpi_grp)
    if(associated(stin%st_range)) then
      i = size(stin%st_range, 1)*size(stin%st_range, 2)
      ALLOCATE(stout%st_range(2,0:stin%mpi_grp%size), i)
      stout%st_range = stin%st_range
    end if
    if(associated(stin%st_num)) then
      i = size(stin%st_num, 1)
      ALLOCATE(stout%st_num(0:stin%mpi_grp%size), i)
      stout%st_num = stin%st_num
    end if

    call pop_sub()
  end subroutine states_copy


  ! ---------------------------------------------------------
  subroutine states_end(st)
    type(states_t), intent(inout) :: st

    call push_sub('states.states_end')

    if(associated(st%rho)) then
      deallocate(st%rho); nullify(st%rho)
    end if

    if(associated(st%occ)) then
      deallocate(st%occ); nullify(st%occ)
    end if

    if(associated(st%eigenval)) then
      deallocate(st%eigenval); nullify(st%eigenval)
    end if

    if(associated(st%momentum)) then
      deallocate(st%momentum); nullify(st%momentum)
    end if

    if(associated(st%node)) then
      deallocate(st%node); nullify(st%node)
    end if

    if(associated(st%j)) then
      deallocate(st%j)
      nullify(st%j)
    end if

    if(associated(st%rho_core)) then
      deallocate(st%rho_core)
      nullify(st%rho_core)
    end if

    if(st%d%ispin==SPINORS .and. associated(st%spin)) then
      deallocate(st%spin); nullify(st%spin)
    end if

    if(associated(st%dpsi)) then
      deallocate(st%dpsi); nullify(st%dpsi)
    end if

    if(associated(st%zpsi)) then
      deallocate(st%zpsi); nullify(st%zpsi)
    end if

    call states_dim_end(st%d)

    if(associated(st%user_def_states)) then
      deallocate(st%user_def_states); nullify(st%user_def_states)
    end if

    if(st%parallel_in_states) then
      deallocate(st%st_range)
      nullify(st%st_range)
      deallocate(st%st_num)
      nullify(st%st_num)
      deallocate(st%ap%schedule)
      nullify(st%ap%schedule)
    end if

    call pop_sub()
  end subroutine states_end


  ! ---------------------------------------------------------
  ! Calculates the new density out the wavefunctions and
  ! occupations...
  subroutine states_dens_accumulate(st, np, rho, ist)
    type(states_t), intent(in)    :: st
    integer,        intent(in)    :: np
    FLOAT,          intent(inout) :: rho(:,:)
    integer,        intent(in)    :: ist

    integer :: ip, ik, sp
    CMPLX   :: c
    type(profile_t), save :: prof

    call push_sub('states.states_dens_accumulate')
    call profiling_in(prof, "CALC_DENSITY")

    sp = 1
    if(st%d%ispin == SPIN_POLARIZED) sp = 2

    do ik = 1, st%d%nik, sp

      if (st%wfs_type == M_REAL) then
        do ip = 1, np
          rho(ip, 1) = rho(ip, 1) + st%d%kweights(ik)*st%occ(ist, ik)*st%dpsi(ip, 1, ist, ik)**2
        end do
      else
        do ip = 1, np
          rho(ip, 1) = rho(ip, 1) + st%d%kweights(ik)*st%occ(ist, ik)*&
               (real(st%zpsi(ip, 1, ist, ik), REAL_PRECISION)**2 + aimag(st%zpsi(ip, 1, ist, ik))**2)
        end do
      end if

      select case(st%d%ispin)
      case(SPIN_POLARIZED)
        if (st%wfs_type == M_REAL) then
          do ip = 1, np
            rho(ip, 2) = rho(ip, 2) + st%d%kweights(ik + 1)*st%occ(ist, ik + 1)*(st%dpsi(ip, 1, ist, ik + 1))**2
          end do
        else
          do ip = 1, np
            rho(ip, 2) = rho(ip, 2) + st%d%kweights(ik + 1)*st%occ(ist, ik + 1)*&
                 (real(st%zpsi(ip, 1, ist, ik + 1), REAL_PRECISION)**2 + aimag(st%zpsi(ip, 1, ist, ik + 1))**2)
          end do
        end if
      case(SPINORS) ! in this case wave-functions are always complex
        do ip = 1, np
          rho(ip, 2) = rho(ip, 2) + st%d%kweights(ik)*st%occ(ist, ik)*&
               (real(st%zpsi(ip, 2, ist, ik), REAL_PRECISION)**2 + aimag(st%zpsi(ip, 2, ist, ik))**2)
          
          c = st%d%kweights(ik)*st%occ(ist, ik)*st%zpsi(ip, 1, ist, ik)*conjg(st%zpsi(ip, 2, ist, ik))
          rho(ip, 3) = rho(ip, 3) + real(c, REAL_PRECISION)
          rho(ip, 4) = rho(ip, 4) + aimag(c)
        end do
      end select
      
    end do

    call profiling_out(prof)

    call pop_sub()
  end subroutine states_dens_accumulate

  subroutine states_dens_reduce(st, np, rho)
    type(states_t), intent(in)    :: st
    integer,        intent(in)    :: np
    FLOAT,          intent(inout) :: rho(:,:)

#ifdef HAVE_MPI
    integer :: ispin
    FLOAT,  allocatable :: reduce_rho(:)
    type(profile_t), save :: reduce_prof
#endif

    call push_sub('states.states_dens_reduce')

#ifdef HAVE_MPI
    ! reduce density
    if(st%parallel_in_states) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
      ALLOCATE(reduce_rho(1:np), np)
      do ispin = 1, st%d%nspin
        call MPI_Allreduce(rho(1, ispin), reduce_rho(1), np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
        call lalg_copy(np, reduce_rho, rho(:, ispin))
      end do
      deallocate(reduce_rho)
      call profiling_out(reduce_prof)
    end if
#endif

    call pop_sub()
  end subroutine states_dens_reduce

  subroutine states_calc_dens(st, np, rho)
    type(states_t), intent(in)  :: st
    integer,        intent(in)  :: np
    FLOAT,          intent(out) :: rho(:,:)

    integer :: ispin, ist

    call push_sub('states.states_calc_dens')

    do ispin = 1, st%d%nspin
      !$omp parallel workshare
      rho(1:np, ispin) = M_ZERO
      !$omp end parallel workshare
    end do

    do ist = st%st_start, st%st_end
      call states_dens_accumulate(st, np, rho, ist)
    end do

    call states_dens_reduce(st, np, rho)

    call pop_sub()
  end subroutine states_calc_dens

  ! ---------------------------------------------------------
  ! generate a hydrogen s-wavefunction around a random point
  subroutine states_generate_random(st, m, ist_start_, ist_end_)
    type(states_t),    intent(inout) :: st
    type(mesh_t),      intent(in)    :: m
    integer, optional, intent(in)    :: ist_start_, ist_end_

    integer :: ist, ik, id, ist_start, ist_end, j, seed
    CMPLX   :: alpha, beta

    call push_sub('states.states_generate_random')

    ist_start = 1
    if(present(ist_start_)) ist_start = ist_start_
    ist_end = st%nst
    if(present(ist_end_)) ist_end = ist_end_

    if(st%parallel_in_states) then
      seed = st%mpi_grp%rank
    else
      seed = 0
    end if

    select case(st%d%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)

      do ik = 1, st%d%nik
        do ist = ist_start, ist_end
          if (st%wfs_type == M_REAL) then
            call dmf_random(m, st%dpsi(:, 1, ist, ik), seed)
          else
            call zmf_random(m, st%zpsi(:, 1, ist, ik), seed)
          end if
          st%eigenval(ist, ik) = M_ZERO
        end do
      end do

    case(SPINORS)

      ASSERT(st%wfs_type == M_CMPLX)

      if(st%fixed_spins) then
        do ik = 1, st%d%nik
          do ist = ist_start, ist_end
            call zmf_random(m, st%zpsi(:, 1, ist, ik))
            ! In this case, the spinors are made of a spatial part times a vector [alpha beta]^T in 
            ! spin space (i.e., same spatial part for each spin component). So (alpha, beta)
            ! determines the spin values. The values of (alpha, beta) can be be obtained
            ! with simple formulae from <Sx>, <Sy>, <Sz>.
            !
            ! Note that here we orthonormalize the orbital part. This ensures that the spinors
            ! are untouched later in the general orthonormalization, and therefore the spin values
            ! of each spinor remain the same.
            do j = ist_start, ist - 1
              st%zpsi(:, 1, ist, ik) = st%zpsi(:, 1, ist, ik) - &
                                       zmf_dotp(m, st%zpsi(:, 1, ist, ik), st%zpsi(:, 1, j, ik)) * &
                                       st%zpsi(:, 1, j, ik)
            end do
            st%zpsi(:, 1, ist, ik) = st%zpsi(:, 1, ist, ik) / zmf_nrm2(m, st%zpsi(:, 1, ist, ik))
            st%zpsi(:, 2, ist, ik) = st%zpsi(:, 1, ist, ik)
            alpha = cmplx(sqrt(M_HALF + st%spin(3, ist, ik)), M_ZERO, REAL_PRECISION)
            beta  = cmplx(sqrt(M_ONE - abs(alpha)**2), M_ZERO, REAL_PRECISION)
            if(abs(alpha) > M_ZERO) then
              beta = cmplx(st%spin(1, ist, ik) / abs(alpha), st%spin(2, ist, ik) / abs(alpha), REAL_PRECISION)
            end if
            st%zpsi(:, 1, ist, ik) = alpha * st%zpsi(:, 1, ist, ik)
            st%zpsi(:, 2, ist, ik) = beta * st%zpsi(:, 2, ist, ik)
            st%eigenval(ist, ik) = M_ZERO
          end do
        end do
      else
        do ik = 1, st%d%nik
          do ist = ist_start, ist_end
            do id = 1, st%d%dim
              call zmf_random(m, st%zpsi(:, id, ist, ik))
            end do
            st%eigenval(ist, ik) = M_ZERO
          end do
        end do
      end if

    end select

    call pop_sub()
  end subroutine states_generate_random


  ! ---------------------------------------------------------
  subroutine states_orthogonalize(st, m, ik_, start_)
    type(states_t),    intent(inout) :: st
    type(mesh_t),      intent(in)    :: m
    integer, optional, intent(in)    :: ik_, start_

    integer :: ik, ik_start, ik_end
    integer :: start

    start = 1
    if(present(start_)) start = start_
    if(present(ik_)) then
      ik_start = ik_
      ik_end   = ik_
    else
      ik_start = 1
      ik_end   = st%d%nik
    end if
    
    do ik = ik_start, ik_end
      if (st%wfs_type == M_REAL) then
        call dstates_gram_schmidt_full(st, st%nst, m, st%d%dim, st%dpsi(:,:,:,ik), start)
      else
        call zstates_gram_schmidt_full(st, st%nst, m, st%d%dim, st%zpsi(:,:,:,ik), start)
      end if
    end do

  end subroutine states_orthogonalize


  ! ---------------------------------------------------------
  subroutine states_fermi(st, m)
    type(states_t), intent(inout) :: st
    type(mesh_t),   intent(in)    :: m

    ! Local variables.
    integer            :: ie, ik, iter
    integer, parameter :: nitmax = 200
    FLOAT              :: drange, t, emin, emax, sumq
    FLOAT, parameter   :: tol = CNST(1.0e-10)
    logical            :: conv
#if defined(HAVE_MPI)
    integer            :: j
    integer            :: tmp
    FLOAT, allocatable :: lspin(:, :) ! To exchange spin.
#endif

    call push_sub('states.fermi')

    if(st%fixed_occ) then ! nothing to do
      ! Calculate magnetizations...
      if(st%d%ispin == SPINORS) then
        do ik = 1, st%d%nik
          do ie = st%st_start, st%st_end
            if (st%wfs_type == M_REAL) then
              write(message(1),'(a)') 'Internal error in states_fermi'
              call write_fatal(1)
            else
              st%spin(1:3, ie, ik) = state_spin(m, st%zpsi(:, :, ie, ik))
            end if
          end do
#if defined(HAVE_MPI)
          if(st%parallel_in_states) then
            ALLOCATE(lspin(3, st%lnst), 3*st%lnst)
            lspin = st%spin(1:3, st%st_start:st%st_end, ik)
            do j = 1, 3
              call lmpi_gen_alltoallv(st%lnst, lspin(j, :), tmp, st%spin(j, :, ik), st%mpi_grp)
            end do
            deallocate(lspin)
          end if
#endif
        end do
      end if
      call pop_sub()
      return
    end if

    ! Initializations
    emin = minval(st%eigenval)
    emax = maxval(st%eigenval)

    if(st%d%ispin == SPINORS) then
      sumq = real(st%nst, REAL_PRECISION)
    else
      sumq = M_TWO*st%nst
    end if

    t = max(st%el_temp, CNST(1.0e-6))
    st%ef = emax

    conv = .true.
    if (abs(sumq - st%qtot) > tol) conv = .false.
    if (conv) then ! all orbitals are full; nothing to be done
      st%occ = M_TWO/st%d%spin_channels!st%d%nspin
      ! Calculate magnetizations...
      if(st%d%ispin == SPINORS) then
        do ik = 1, st%d%nik
          do ie = st%st_start, st%st_end
            st%spin(1:3, ie, ik) = state_spin(m, st%zpsi(:, :, ie, ik))
          end do
#if defined(HAVE_MPI)
          if(st%parallel_in_states) then
            ALLOCATE(lspin(3, st%lnst), 3*st%lnst)
            lspin = st%spin(1:3, st%st_start:st%st_end, ik)
            do j = 1, 3
              call lmpi_gen_alltoallv(st%lnst, lspin(j, :), tmp, st%spin(j, :, ik), st%mpi_grp)
            end do
            deallocate(lspin)
          end if
#endif
        end do
      end if
      call pop_sub()
      return
    end if

    if (sumq < st%qtot) then ! not enough states
      message(1) = 'Fermi: Not enough states'
      write(message(2),'(6x,a,f12.6,a,f12.6)')'(total charge = ', st%qtot, &
        ' max charge = ', sumq
      call write_fatal(2)
    end if

    drange = t*sqrt(-log(tol*CNST(.01)))

    emin = emin - drange
    emax = emax + drange

    do iter = 1, nitmax
      st%ef = M_HALF*(emin + emax)
      sumq  = M_ZERO

      do ik = 1, st%d%nik
        do ie = 1, st%nst
          sumq = sumq + st%d%kweights(ik)/st%d%spin_channels * & !st%d%nspin * &
            stepf((st%eigenval(ie, ik) - st%ef)/t)
        end do
      end do

      conv = .true.
      if(abs(sumq - st%qtot) > tol) conv = .false.
      if(conv) exit

      if(sumq <= st%qtot ) emin = st%ef
      if(sumq >= st%qtot ) emax = st%ef
    end do

    if(iter == nitmax) then
      message(1) = 'Fermi: did not converge'
      call write_fatal(1)
    end if

    do ik = 1, st%d%nik
      do ie = 1, st%nst
        st%occ(ie, ik) = stepf((st%eigenval(ie, ik) - st%ef)/t)/st%d%spin_channels!st%d%nspin
      end do
    end do

    ! Calculate magnetizations...
    if(st%d%ispin == SPINORS) then
      do ik = 1, st%d%nik
        do ie = st%st_start, st%st_end
          st%spin(1:3, ie, ik) = state_spin(m, st%zpsi(:, :, ie, ik))
        end do
#if defined(HAVE_MPI)
        if(st%parallel_in_states) then
          ALLOCATE(lspin(3, st%lnst), 3*st%lnst)
          lspin = st%spin(1:3, st%st_start:st%st_end, ik)
          do j = 1, 3
            call lmpi_gen_alltoallv(st%lnst, lspin(j, :), tmp, st%spin(j, :, ik), st%mpi_grp)
          end do
          deallocate(lspin)
        end if
#endif
      end do
    end if

    call pop_sub()
  end subroutine states_fermi


  ! ---------------------------------------------------------
  ! function to calculate the eigenvalues sum using occupations as weights
  function states_eigenvalues_sum(st, x) result(e)
    type(states_t), intent(in)  :: st
    FLOAT                       :: e
    FLOAT, optional, intent(in) :: x(st%st_start:st%st_end, 1:st%d%nik)

    integer :: ik
#ifdef HAVE_MPI
    FLOAT   :: s
#endif

    call push_sub('states.states_eigenvalues_sum')

    e = M_ZERO
    do ik = 1, st%d%nik
      if(present(x)) then
        e = e + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik)* &
          x(st%st_start:st%st_end, ik))
      else
        e = e + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik)* &
          st%eigenval(st%st_start:st%st_end, ik))
      end if
    end do

#ifdef HAVE_MPI
    if(st%parallel_in_states) then
      call MPI_Allreduce(e, s, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      e = s
    end if
#endif

    call pop_sub()
  end function states_eigenvalues_sum


  ! ---------------------------------------------------------
  subroutine states_write_eigenvalues(iunit, nst, st, sb, error)
    integer,           intent(in) :: iunit, nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb
    FLOAT,             intent(in), optional :: error(nst, st%d%nik)

    integer ik, j, ns, is
    FLOAT :: o
    character(len=80) tmp_str(MAX_DIM), cspin

    call push_sub('states.states_write_eigenvalues')

    ns = 1
    if(st%d%nspin == 2) ns = 2

    message(1) = 'Eigenvalues [' // trim(units_out%energy%abbrev) // ']'
    call write_info(1, iunit)
    if (st%d%nik > ns) then
      message(1) = 'Kpoints [' // trim(units_out%length%abbrev) // '^-1]'
      call write_info(1, iunit)
    end if

    if(.not.mpi_grp_is_root(mpi_world)) return

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then
        write(message(1), '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
          st%d%kpoints(1, ik)*units_out%length%factor, ',',            &
          st%d%kpoints(2, ik)*units_out%length%factor, ',',            &
          st%d%kpoints(3, ik)*units_out%length%factor, ')'
        call write_info(1, iunit)
      end if

      if(present(error)) then
        if(st%d%ispin .eq. SPINORS) then
          write(message(1), '(a4,1x,a5,1x,a12,1x,a12,2x,a4,4x,a4,4x,a4,5x,a5)')   &
            '#st',' Spin',' Eigenvalue', 'Occupation ', '<Sx>', '<Sy>', '<Sz>', 'Error'
        else
          write(message(1), '(a4,1x,a5,1x,a12,4x,a12,1x,a10)')   &
            '#st',' Spin',' Eigenvalue', 'Occupation ', 'Error'
        end if
      else
        if(st%d%ispin .eq. SPINORS) then
          write(message(1), '(a4,1x,a5,1x,a12,1x,a12,2x,a4,4x,a4,4x,a4)')   &
            '#st',' Spin',' Eigenvalue', 'Occupation ', '<Sx>', '<Sy>', '<Sz>'
        else
          write(message(1), '(a4,1x,a5,1x,a12,4x,a12,1x)')       &
            '#st',' Spin',' Eigenvalue', 'Occupation '
        end if
      end if
      call write_info(1, iunit)

      do j = 1, nst
        do is = 0, ns-1
          if(j > st%nst) then
            o = M_ZERO
          else
            o = st%occ(j, ik+is)
          end if

          if(is.eq.0) cspin = 'up'
          if(is.eq.1) cspin = 'dn'
          if(st%d%ispin.eq.UNPOLARIZED.or.st%d%ispin.eq.SPINORS) cspin = '--'

          write(tmp_str(1), '(i4,3x,a2)') j, trim(cspin)
          if(simul_box_is_periodic(sb)) then
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,3x,4f5.2)') &
                (st%eigenval(j, ik)-st%ef)/units_out%energy%factor, o, st%spin(1:3, j, ik)
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                (st%eigenval(j, ik+is))/units_out%energy%factor, o
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik), ')'
            end if
          else
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,5x,f5.2,3x,3f8.4)') &
                st%eigenval(j, ik)/units_out%energy%factor, o, st%spin(1:3, j, ik)
              if(present(error)) write(tmp_str(3), '(a3,es7.1,a1)')'  (', error(j, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                st%eigenval(j, ik+is)/units_out%energy%factor, o
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik), ')'
            end if
          end if
          if(present(error)) then
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))//trim(tmp_str(3))
          else
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))
          end if
          call write_info(1, iunit)
        end do
      end do
    end do

    call pop_sub()
  end subroutine states_write_eigenvalues


  ! ---------------------------------------------------------
  subroutine states_write_bands(dir, nst, st, sb)
    character(len=*),  intent(in) :: dir    
    integer,           intent(in) :: nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb

    integer :: i, ik, j, ns, is
    integer, allocatable :: iunit(:)
    FLOAT   :: factor(MAX_DIM)
    logical :: grace_mode, gnuplot_mode
    character(len=80) :: filename    

    call push_sub('states.states_write_bands')

    if(.not.mpi_grp_is_root(mpi_world)) return

    !%Variable OutputBandsGnuplotMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in gnuplot friendly format
    !%End
    call loct_parse_logical(check_inp('OutputBandsGnuplotMode'), .true., gnuplot_mode)

    !%Variable OutputBandsGraceMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in grace friendly format
    !%End
    call loct_parse_logical(check_inp('OutputBandsGraceMode'), .false., grace_mode)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    ALLOCATE(iunit(0:ns-1), ns)

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do i = 1, MAX_DIM
      factor(i) = M_ONE
      if (sb%klattice(i,i) /= M_ZERO) factor(i) = sb%klattice(i,i)
    end do

    if (gnuplot_mode) then
      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i1.1,a)') 'bands-gp-', is+1,'.dat'
        else
          write(filename, '(a)') 'bands-gp.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a, i6)') '# kx ky kz (unscaled), kx ky kz (scaled), bands:', nst
      end do

      ! output bands in gnuplot format
      do j = 1, nst
        do ik = 1, st%d%nik, ns
          do is = 0, ns-1
            write(iunit(is), '(1x,6f14.8,3x,f14.8)')            &
              st%d%kpoints(1:MAX_DIM, ik+is),                   & ! unscaled
              st%d%kpoints(1:MAX_DIM, ik+is)/factor(1:MAX_DIM), & ! scaled
              st%eigenval(j, ik+is)/units_out%energy%factor
          end do
        end do
        do is = 0, ns-1
          write(iunit(is), '(a)') ''
        end do
      end do
      do is = 0, ns-1
        call io_close(iunit(is))
      end do
    end if

    if (grace_mode) then
      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i1.1,a)') 'bands-grace-', is+1,'.dat'
        else
          write(filename, '(a)') 'bands-grace.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a, i6)') '# kx ky kz (unscaled), kx ky kz (scaled), bands:', nst
      end do

      ! output bands in xmgrace format, i.e.:
      ! k_x, k_y, k_z, e_1, e_2, ..., e_n
      do ik = 1, st%d%nik, ns
        do is = 0, ns-1
          write(iunit(is), '(1x,6f14.8,3x,16384f14.8)')         &
            st%d%kpoints(1:MAX_DIM, ik+is),                     & ! unscaled
            st%d%kpoints(1:MAX_DIM, ik+is)/factor(1:MAX_DIM),   & ! scaled
            (st%eigenval(j, ik+is)/units_out%energy%factor, j = 1, nst)
        end do
      end do
      do is = 0, ns-1
        call io_close(iunit(is))
      end do        
    end if

    deallocate(iunit)

    call pop_sub()
  end subroutine states_write_bands


  ! ---------------------------------------------------------
  subroutine states_write_dos(dir, st)
    character(len=*), intent(in) :: dir
    type(states_t),   intent(in) :: st

    integer :: ie, ik, ist, epoints, is, ns
    integer, allocatable :: iunit(:)
    FLOAT   :: emin, emax, de, gamma, energy
    FLOAT   :: evalmax, evalmin, tdos, eextend
    FLOAT, allocatable :: dos(:,:,:)
    character(len=64)  :: filename

    call push_sub('states.states_write_dos')

    evalmin = minval(st%eigenval/units_out%energy%factor)
    evalmax = maxval(st%eigenval/units_out%energy%factor)
    ! we extend the energy mesh by this amount
    eextend  = (evalmax - evalmin) / M_FOUR

    !%Variable DOSEnergyMin
    !%Type float
    !%Default 
    !%Section Output
    !%Description
    !% Lower bound for the energy mesh of the DOS
    !%End
    call loct_parse_float(check_inp('DOSEnergyMin'), evalmin - eextend, emin)

    !%Variable DOSEnergyMax
    !%Type float
    !%Default 
    !%Section Output
    !%Description
    !% Upper bound for the energy mesh of the DOS
    !%End
    call loct_parse_float(check_inp('DOSEnergyMax'), evalmax + eextend, emax)

    !%Variable DOSEnergyPoints
    !%Type integer
    !%Default 500
    !%Section Output
    !%Description
    !% Determines how many energy points octopus should use for 
    !% the DOS energy grid
    !%End
    call loct_parse_int(check_inp('DOSEnergyPoints'), 500, epoints)

    !%Variable DOSGamma
    !%Type float
    !%Default 
    !%Section Output
    !%Description
    !% Determines the width of the Lorentzian which is used to sum 
    !% up the DOS sum
    !%End
    call loct_parse_float(check_inp('DOSGamma'), &
      CNST(0.008)/units_out%energy%factor, gamma)

    ! spacing for energy mesh
    de = (emax - emin) / (epoints - 1)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    ! space for state dependent DOS
    ALLOCATE(dos(epoints, st%nst, 0:ns-1), epoints*st%nst*ns)
    ALLOCATE(iunit(0:ns-1), ns)    

    ! compute band/spin resolved density of states
    do ist = 1, st%nst

      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i4.4,a,i1.1,a)') 'dos-', ist, '-', is+1,'.dat'
        else
          write(filename, '(a,i4.4,a)') 'dos-', ist, '.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a)') '# energy, band resolved DOS'
      end do

      do ie = 1, epoints
        energy = emin + (ie - 1) * de
        dos(ie, ist, :) = M_ZERO
        ! sum up Lorentzians
        do ik = 1, st%d%nik, ns
          do is = 0, ns-1
            dos(ie, ist, is) = dos(ie, ist, is) + st%d%kweights(ik+is) * M_ONE/M_Pi * &
              gamma / ( (energy - st%eigenval(ist, ik+is)/units_out%energy%factor)**2 + gamma**2 )
          end do
        end do
        do is = 0, ns-1
          write(message(1), '(2f12.6)') energy, dos(ie, ist, is)
          call write_info(1, iunit(is))
        end do
      end do

      do is = 0, ns-1
        call io_close(iunit(is))
      end do
    end do

    ! for spin polarized calculations also output spin resolved tdos
    if(st%d%nspin .gt. 1) then    
      do is = 0, ns-1
        write(filename, '(a,i1.1,a)') 'total-dos-', is+1,'.dat'
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a)') '# energy, total DOS (spin resolved)'

        do ie = 1, epoints
          energy = emin + (ie - 1) * de
          tdos = M_ZERO
          do ist = 1, st%nst
            tdos = tdos + dos(ie, ist, is)
          end do
          write(message(1), '(2f12.6)') energy, tdos
          call write_info(1, iunit(is))
        end do

        call io_close(iunit(is))
      end do
    end if


    iunit(0) = io_open(trim(dir)//'/'//'total-dos.dat', action='write')    

    ! compute total density of states
    do ie = 1, epoints
      energy = emin + (ie - 1) * de
      tdos = M_ZERO
      do ist = 1, st%nst
        do is = 0, ns-1
          tdos = tdos + dos(ie, ist, is)
        end do
      end do
      write(message(1), '(2f12.6)') energy, tdos
      call write_info(1, iunit(0))
    end do

    call io_close(iunit(0))

    deallocate(iunit, dos)

    call pop_sub()
  end subroutine states_write_dos


  ! ---------------------------------------------------------
  subroutine states_write_fermi_energy(dir, st, m, sb)
    character(len=*),  intent(in) :: dir
    type(states_t), intent(inout) :: st
    type(mesh_t),      intent(in) :: m
    type(simul_box_t), intent(in) :: sb

    integer :: iunit, i
    FLOAT :: scale, maxdos
    FLOAT :: factor(MAX_DIM)

    call push_sub('states.states_write_fermi_energy')

    call states_fermi(st, m)

    iunit = io_open(trim(dir)//'/'//'bands-efermi.dat', action='write')    

    scale = units_out%energy%factor

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do i = 1, MAX_DIM
      factor(i) = M_ONE
      if (sb%klattice(i,i) /= M_ZERO) factor(i) = sb%klattice(i,i)
    end do

    ! write fermi energy in a format that can be used together 
    ! with bands.dat
    write(message(1), '(a)') '# Fermi energy in a format compatible with bands-gp.dat'

    write(message(2), '(7f12.6)')          &
      minval(st%d%kpoints(1,:)),           &
      minval(st%d%kpoints(2,:)),           &
      minval(st%d%kpoints(3,:)),           &
      minval(st%d%kpoints(1,:)/factor(1)), &
      minval(st%d%kpoints(2,:)/factor(2)), &
      minval(st%d%kpoints(3,:)/factor(3)), &
      st%ef/scale

    ! Gamma point
    write(message(3), '(7f12.6)')          &
      (M_ZERO, i = 1, 6),                  &
      st%ef/scale

    write(message(4), '(7f12.6)')          &
      maxval(st%d%kpoints(1,:)),           &
      maxval(st%d%kpoints(2,:)),           &
      maxval(st%d%kpoints(3,:)),           &
      maxval(st%d%kpoints(1,:)/factor(1)), &
      maxval(st%d%kpoints(2,:)/factor(2)), &
      maxval(st%d%kpoints(3,:)/factor(3)), &
      st%ef/scale

    call write_info(4, iunit)
    call io_close(iunit)

    ! now we write the same information so that it can be used 
    ! together with total-dos.dat
    iunit = io_open(trim(dir)//'/'//'total-dos-efermi.dat', action='write')    

    write(message(1), '(a)') '# Fermi energy in a format compatible with total-dos.dat'    

    ! this is the maximum that tdos can reach
    maxdos = sum(st%d%kweights) * st%nst

    write(message(2), '(4f12.6)') st%ef/scale, M_ZERO
    write(message(3), '(4f12.6)') st%ef/scale, maxdos

    call write_info(3, iunit)
    call io_close(iunit)

    call pop_sub()
  end subroutine states_write_fermi_energy


  ! -------------------------------------------------------
  subroutine states_degeneracy_matrix(st)
    type(states_t), intent(in) :: st

    integer :: is, js, inst, inik, dsize, iunit
    integer, allocatable :: eindex(:,:), sindex(:)
    integer, allocatable :: degeneracy_matrix(:, :)
    FLOAT,   allocatable :: eigenval_sorted(:)
    FLOAT :: degen_thres, evis, evjs

    call push_sub('states.states_degeneracy_matrix')

    ALLOCATE(eigenval_sorted(st%nst*st%d%nik),   st%nst*st%d%nik)
    ALLOCATE(         sindex(st%nst*st%d%nik),   st%nst*st%d%nik)
    ALLOCATE(      eindex(2, st%nst*st%d%nik), 2*st%nst*st%d%nik)
    dsize = st%nst*st%d%nik * st%nst*st%d%nik
    ALLOCATE(degeneracy_matrix(st%nst*st%d%nik, st%nst*st%d%nik), dsize)

    ! convert double index "inst, inik" to single index "is"
    ! and keep mapping array
    is = 1
    do inst = 1, st%nst
      do inik = 1, st%d%nik
        eigenval_sorted(is) = st%eigenval(inst, inik)        
        eindex(1, is) = inst
        eindex(2, is) = inik
        is = is + 1
      end do
    end do

    ! sort eigenvalues
    call sort(eigenval_sorted, sindex)

    !%Variable DegeneracyMatrixThreshold
    !%Type float
    !%Default 1e-5
    !%Section States
    !%Description
    !% A state j with energy E_j will be considered degenerate with a state
    !% with energy E_i, if  E_i - threshold < E_j < E_i + threshold.
    !%End
    call loct_parse_float(check_inp('DegeneracyMatrixThreshold'), CNST(1e-5), degen_thres)    

    ! setup degeneracy matrix. the matrix summarizes the degeneracy relations 
    ! among the states
    degeneracy_matrix = 0

    do is = 1, st%nst*st%d%nik
      do js = 1, st%nst*st%d%nik

        ! a state is always degenerate to itself
        if ( is.eq.js ) cycle

        evis = st%eigenval(eindex(1, sindex(is)), eindex(2, sindex(is)))
        evjs = st%eigenval(eindex(1, sindex(js)), eindex(2, sindex(js)))

        ! is evjs in the "evis plus minus threshold" bracket?
        if( (evjs.gt.evis - degen_thres).and.(evjs.lt.evis + degen_thres) ) then
          ! mark forward scattering states with +1 and backward scattering
          ! states with -1
          degeneracy_matrix(is, js) = &
            sign(M_ONE, st%momentum(1, eindex(1, sindex(js)), eindex(2, sindex(js))))
        end if

      end do
    end do

    if(mpi_grp_is_root(mpi_world)) then

      ! write matrix to "restart/gs" directory
      iunit = io_open(trim(tmpdir)//'gs/degeneracy_matrix', action='write', is_tmp = .true.)

      write(iunit, '(a)') '# index  kx ky kz  eigenvalue  degeneracy matrix'

      do is = 1, st%nst*st%d%nik
        write(iunit, '(i6,4e24.16,32767i3)') is, st%d%kpoints(:, eindex(2, sindex(is))), &
          eigenval_sorted(is), (degeneracy_matrix(is, js), js = 1, st%nst*st%d%nik)
      end do

      call io_close(iunit)

      ! write index vectors to "restart/gs" directory
      iunit = io_open(trim(tmpdir)//'gs/index_vectors', action='write', is_tmp = .true.)    

      write(iunit, '(a)') '# index  sindex  eindex1 eindex2'

      do is = 1, st%nst*st%d%nik
        write(iunit,'(4i6)') is, sindex(is), eindex(1, sindex(is)), eindex(2, sindex(is))
      end do

      call io_close(iunit)
    end if

    deallocate(eigenval_sorted, sindex, eindex)
    deallocate(degeneracy_matrix)

    call pop_sub()
  end subroutine states_degeneracy_matrix


  ! -------------------------------------------------------
  integer function states_spin_channel(ispin, ik, dim)
    integer, intent(in) :: ispin, ik, dim

    ASSERT(ispin >= UNPOLARIZED .or. ispin <= SPINORS)
    ASSERT(ik > 0)
    ASSERT(dim==1 .or. dim==2)
    ASSERT(.not.(ispin.ne.3 .and. dim==2))

    select case(ispin)
    case(1); states_spin_channel = 1
    case(2); states_spin_channel = mod(ik+1, 2)+1
    case(3); states_spin_channel = dim
    case default; states_spin_channel = -1
    end select

  end function states_spin_channel


  ! ---------------------------------------------------------
  subroutine states_distribute_nodes(st, mc)
    type(states_t),    intent(inout) :: st
    type(multicomm_t), intent(in)    :: mc

#ifdef HAVE_MPI
    integer :: sn, sn1, r, j, k
#endif

    call push_sub('states.states_distribute_nodes')

    ! Defaults.
    st%node(:)            = 0
    st%st_start           = 1
    st%st_end             = st%nst
    st%lnst               = st%nst
    st%parallel_in_states = .false.
    call mpi_grp_init(st%mpi_grp, -1)

#if defined(HAVE_MPI)
    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then
      st%parallel_in_states = .true.
      call mpi_grp_init(st%mpi_grp, mc%group_comm(P_STRATEGY_STATES))
      call mpi_grp_init(st%dom_st, mc%dom_st_comm)
      call multicomm_create_all_pairs(st%mpi_grp, st%ap)

     if(st%nst < st%mpi_grp%size) then
       message(1) = "Have more processors than necessary"
       write(message(2),'(i4,a,i4,a)') st%mpi_grp%size, " processors and ", st%nst, " states."
       call write_fatal(2)
     end if

     ALLOCATE(st%st_range(2, 0:st%mpi_grp%size-1), 2*st%mpi_grp%size)
     ALLOCATE(st%st_num(0:st%mpi_grp%size-1), st%mpi_grp%size)

     call multicomm_divide_range(st%nst, st%mpi_grp%size, st%st_range(1, :), st%st_range(2, :), st%st_num)

     do k = 0, st%mpi_grp%size - 1
       write(message(1),'(a,i4,a,i7,a,i7)') &
            'Info: Nodes in states-group ', k, ' will manage states', st%st_range(1, k), " - ", st%st_range(2, k)
       call write_info(1)
       if(st%mpi_grp%rank .eq. k) then
         st%st_start = st%st_range(1, k)
         st%st_end   = st%st_range(2, k)
         st%lnst     = st%st_num(k)
       endif
     end do

     sn  = st%nst/st%mpi_grp%size
     sn1 = sn + 1
     r  = mod(st%nst, st%mpi_grp%size)
     do j = 1, r
       st%node((j-1)*sn1+1:j*sn1) = j - 1
     end do
     k = sn1*r
     call MPI_Barrier(st%mpi_grp%comm, mpi_err)
     do j = 1, st%mpi_grp%size - r
       st%node(k+(j-1)*sn+1:k+j*sn) = r + j - 1
     end do
   end if
#endif

    call pop_sub()
  end subroutine states_distribute_nodes


  ! ---------------------------------------------------------
  logical function wfs_are_complex(st) result (wac)
    type(states_t),    intent(in) :: st
    wac = (st%wfs_type == M_CMPLX)
  end function wfs_are_complex


  ! ---------------------------------------------------------
  logical function wfs_are_real(st) result (war)
    type(states_t),    intent(in) :: st
    war = (st%wfs_type == M_REAL)
  end function wfs_are_real


  ! ---------------------------------------------------------
  subroutine states_dump(st, iunit)
    type(states_t), intent(in) :: st
     integer,       intent(in) :: iunit

     call push_sub('states.states_dump')
     
     write(iunit, '(a20,1i10)')  'nst=                ', st%nst
     write(iunit, '(a20,1i10)')  'dim=                ', st%d%dim
     write(iunit, '(a20,1i10)')  'nik=                ', st%d%nik

     call pop_sub()
  end subroutine states_dump


  ! ---------------------------------------------------------
  subroutine states_calc_tau_jp_gn(gr, st, tau, jp, grho)
    type(grid_t),    intent(inout) :: gr
    type(states_t),  intent(inout) :: st
    FLOAT, optional, intent(out)   ::  tau(:,:)    ! (NP, st%d%nspin)
    FLOAT, optional, intent(out)   ::   jp(:,:,:)  ! (NP, NDIM, st%d%nspin)
    FLOAT, optional, intent(out)   :: grho(:,:,:)  ! (NP, NDIM, st%d%nspin)

    CMPLX, allocatable :: wf_psi(:,:), gwf_psi(:,:,:)
    CMPLX   :: c_tmp
    integer :: sp, is, ik, ik_tmp, ist, i_dim, st_dim, ii
    FLOAT   :: ww

#if defined(HAVE_MPI)
    FLOAT, allocatable :: tmp_reduce(:)
    integer :: mpi_err
#endif

    ALLOCATE(wf_psi(NP_PART, st%d%dim),  NP_PART*st%d%dim)
    ALLOCATE(gwf_psi(NP, NDIM, st%d%dim), NP*NDIM*st%d%dim)   

    sp = 1
    if(st%d%ispin == SPIN_POLARIZED) sp = 2

    ASSERT(present( tau).or.present(  jp).or.present(grho))

    if(present( tau))  tau(:,:)   = M_ZERO
    if(present(  jp))   jp(:,:,:) = M_ZERO
    if(present(grho)) grho(:,:,:) = M_ZERO

    do is = 1, sp
      do ik_tmp = 1, st%d%nik, sp
        ik = ik_tmp + is - 1

        do ist = st%st_start, st%st_end

          ! all calculations will be done with complex wave-functions
          if (st%wfs_type == M_REAL) then
            wf_psi(:,:) = cmplx(st%dpsi(:,:, ist, ik), KIND=REAL_PRECISION)
          else
            wf_psi(:,:) = st%zpsi(:,:, ist, ik)
          end if

          ! calculate gradient of the wave-function
          do st_dim = 1, st%d%dim
            call zf_gradient(gr%sb, gr%f_der, wf_psi(:,st_dim), gwf_psi(:,:,st_dim))
          end do

          ww = st%d%kweights(ik)*st%occ(ist, ik)
          do i_dim = 1, NDIM
            if(present(grho)) &
              grho(1:NP, i_dim, is) = grho(1:NP, i_dim, is) + &
                ww*M_TWO*real(conjg(wf_psi(1:NP, 1))*gwf_psi(1:NP, i_dim, 1))
            if(present(  jp)) &
              jp  (1:NP, i_dim, is) = jp  (1:NP, i_dim, is) + &
                ww*aimag(conjg(wf_psi(1:NP, 1))*gwf_psi(1:NP, i_dim, 1))
            if(present( tau)) then
              tau (1:NP, is)        = tau (1:NP, is)        + &
                ww*abs(gwf_psi(1:NP, i_dim, 1))**2
            end if

            if(st%d%ispin == SPINORS) then
              if(present(grho)) then
                grho(1:NP, i_dim, 2) = grho(1:NP, i_dim, 2) + &
                  ww*M_TWO*real(conjg(wf_psi(1:NP, 2))*gwf_psi(1:NP, i_dim, 2))
                grho(1:NP, i_dim, 3) = grho(1:NP, i_dim, 3) + ww* &
                  real (gwf_psi(1:NP, i_dim, 1)*conjg(wf_psi(1:NP, 2)) + &
                    wf_psi(1:NP, 1)*conjg(gwf_psi(1:NP, i_dim, 2)))
                grho(1:NP, i_dim, 4) = grho(1:NP, i_dim, 4) + ww* &
                  aimag(gwf_psi(1:NP, i_dim, 1)*conjg(wf_psi(1:NP, 2)) + &
                    wf_psi(1:NP, 1)*conjg(gwf_psi(1:NP, i_dim, 2)))
              end if
            
              ! the expression for the paramagnetic current with spinors is
              !     j = ( jp(1)             jp(3) + i jp(4) ) 
              !         (-jp(3) + i jp(4)   jp(2)           )
              if(present(  jp)) then
                jp  (1:NP, i_dim, 2) = jp  (1:NP, i_dim, 2) + &
                  ww*aimag(conjg(wf_psi(1:NP, 2))*gwf_psi(1:NP, i_dim, 2))
                do ii = 1, NP
                  c_tmp = conjg(wf_psi(ii, 1))*gwf_psi(ii, i_dim, 2) - wf_psi(ii, 2)*conjg(gwf_psi(ii, i_dim, 1))
                  jp(ii, i_dim, 3) = jp(ii, i_dim, 3) + ww* real(c_tmp)
                  jp(ii, i_dim, 4) = jp(ii, i_dim, 4) + ww*aimag(c_tmp)
                end do
              end if

              ! the expression for the paramagnetic current with spinors is
              !     t = ( tau(1)              tau(3) + i tau(4) ) 
              !         ( tau(3) - i tau(4)   tau(2)            )
              if(present( tau)) then
                tau (1:NP, 2)        = tau (1:NP, 2)        + ww*abs(gwf_psi(1:NP, i_dim, 2))**2
                do ii = 1, NP
                  c_tmp = conjg(gwf_psi(ii, i_dim, 1))*gwf_psi(ii, i_dim, 2)
                  tau(ii, 3) = tau(ii, 3) + ww* real(c_tmp)
                  tau(ii, 4) = tau(ii, 4) + ww*aimag(c_tmp)
                end do
              end if
              
            end if
          end do

        end do
      end do
    end do

    deallocate(wf_psi, gwf_psi)

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      ALLOCATE(tmp_reduce(1:NP), NP)

      do is = 1, st%d%nspin
        if(present(tau)) then
          call MPI_Allreduce(tau(1, is), tmp_reduce(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
          tau(1:NP, is) = tmp_reduce(1:NP)       
        end if

        do i_dim = 1, NDIM
          if(present(jp)) then
            call MPI_Allreduce(jp(1, i_dim, is), tmp_reduce(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
            jp(1:NP, i_dim, is) = tmp_reduce(1:NP)
          end if

          if(present(grho)) then
            call MPI_Allreduce(grho(1, i_dim, is), tmp_reduce(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
            grho(1:NP, i_dim, is) = tmp_reduce(1:NP)
          end if
        end do

      end do
      deallocate(tmp_reduce)
    end if
#endif            
  end subroutine states_calc_tau_jp_gn


  ! ---------------------------------------------------------
  function state_spin(m, f1) result(s)
    FLOAT, dimension(3) :: s
    type(mesh_t), intent(in) :: m
    CMPLX,  intent(in) :: f1(:, :)

    CMPLX :: z

    call push_sub('states.zstate_spin')

    z = zmf_dotp(m, f1(:, 1) , f1(:, 2))

    s(1) = M_TWO * z
    s(2) = M_TWO * aimag(z)
    s(3) = zmf_dotp(m, f1(:, 1), f1(:, 1)) - zmf_dotp(m, f1(:, 2), f1(:, 2))
    s = s * M_HALF ! spin is half the sigma matrix.

    call pop_sub()
  end function state_spin


  ! ---------------------------------------------------------
  ! Reads the state stored in directory "dir", and finds out
  ! the kpoints, dim, and nst contained in it.
  ! ---------------------------------------------------------
  subroutine states_look (dir, m, kpoints, dim, nst, ierr)
    character(len=*), intent(in)    :: dir
    type(mesh_t),     intent(in)    :: m
    integer,          intent(out)   :: kpoints, dim, nst, ierr

    character(len=256) :: line
    character(len=12)  :: filename
    character(len=1)   :: char
    integer :: iunit, iunit2, err, i, ist, idim, ik
    FLOAT :: occ, eigenval

    call push_sub('states.states_look')

    ierr = 0
    iunit  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp = .true., grp = m%mpi_grp)
    if(iunit < 0) then
      ierr = -1
      return
    end if
    iunit2 = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp = .true., grp = m%mpi_grp)
    if(iunit2 < 0) then
      call io_close(iunit, grp = m%mpi_grp)
      ierr = -1
      return
    end if

    ! Skip two lines.
    call iopar_read(m%mpi_grp, iunit, line, err); call iopar_read(m%mpi_grp, iunit, line, err)
    call iopar_read(m%mpi_grp, iunit2, line, err); call iopar_read(m%mpi_grp, iunit2, line, err)

    kpoints = 1
    dim = 1
    nst = 1
    do
      call iopar_read(m%mpi_grp, iunit, line, i)
      read(line, '(a)') char
      if(i.ne.0.or.char=='%') exit
      read(line, *) ik, char, ist, char, idim, char, filename
      if(ik > kpoints) kpoints = ik
      if(idim == 2)    dim     = 2
      if(ist>nst)      nst     = ist
      call iopar_read(m%mpi_grp, iunit2, line, err)
      read(line, *) occ, char, eigenval
    end do

    call io_close(iunit, grp = m%mpi_grp)
    call io_close(iunit2, grp = m%mpi_grp)
    call pop_sub()
  end subroutine states_look


  ! ---------------------------------------------------------
  logical function state_is_local(st, ist)
    type(states_t), intent(in) :: st
    integer,        intent(in) :: ist

    call push_sub('states.state_is_local')

    state_is_local = ist.ge.st%st_start.and.ist.le.st%st_end

    call pop_sub()
  end function state_is_local


  ! ---------------------------------------------------------
  subroutine states_freeze_orbitals(st, gr, mc, n)
    type(states_t), intent(inout) :: st
    type(grid_t),   intent(in)    :: gr
    type(multicomm_t), intent(in) :: mc
    integer,        intent(in)    :: n

    integer :: ispin, ist, ik
    FLOAT, allocatable :: rho(:, :)
    type(states_t) :: staux

    call push_sub('states.states_freeze_orbitals')

    if(n >= st%nst) then
      write(message(1),'(a)') 'Attempting to freeze a number of orbitals which is larger or equal to'
      write(message(2),'(a)') 'the total number. The program has to stop.'
      call write_fatal(2)
    end if

    ! We will put the frozen density into st%rho_core. We will put the total density, summing up
    ! the possible spin-up an spin-down contributions. This could be refined later...
    if(.not.associated(st%rho_core)) then
      ALLOCATE(st%rho_core(gr%m%np), gr%m%np)
      st%rho_core(:) = M_ZERO
    end if

    ALLOCATE(rho(gr%m%np, st%d%nspin), gr%m%np*st%d%nspin)
    do ist = st%st_start, st%st_end
      if(ist > n) cycle
      call states_dens_accumulate(st, gr%m%np, rho, ist)
    end do
    call states_dens_reduce(st, gr%m%np, rho)

    do ispin = 1, st%d%nspin
      st%rho_core(:) = st%rho_core(:) + rho(:, ispin)
    end do

    call states_copy(staux, st)

    st%nst = st%nst - n

    call states_deallocate_wfns(st)
    call states_distribute_nodes(st, mc)
    call states_allocate_wfns(st, gr%m, M_CMPLX)

#if defined(HAVE_MPI) 
    if(staux%parallel_in_states) then

      do ik = 1, st%d%nik
        do ist = staux%st_start, staux%st_end
          if(.not.state_is_local(st, ist-n) then
            call MPI_Send(staux%zpsi(:, :, ist, ik), m%np_part*st%d%dim, R_MPITYPE, staux%node(ist), &
              0, st%mpi_grp%comm, mpi_err)
          end if
   
        end do
      end do

      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end

          if(state_is_local(staux, n+ist)) then
            st%zpsi(:, :, ist, ik) = staux%zpsi(:, :, n + ist, ik)
          else
            call MPI_Recv(st%zpsi(1, 1, ist, ik), m%np_part*st%d%dim, R_MPITYPE, st%node(n+ist), 0, st%mpi_grp%comm, mpi_err)
          end if

        end do
      end do

   else
     do ik = 1, st%d%nik
       do ist = st%st_start, st%st_end
         st%zpsi(:, :, ist, ik) = staux%zpsi(:, :, n + ist, ik)
       end do
     end do
   end if

#else

   do ik = 1, st%d%nik
     do ist = st%st_start, st%st_end
       st%zpsi(:, :, ist, ik) = staux%zpsi(:, :, n + ist, ik)
     end do
   end do

#endif

    call states_end(staux)
    deallocate(rho)
    call pop_sub()
  end subroutine states_freeze_orbitals


#include "states_kpoints.F90"

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
