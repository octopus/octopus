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
  use batch_m
  use blas_m
  use derivatives_m
  use calc_mode_m
  use crystal_m
  use distributed_m
  use datasets_m
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
  use modelmb_particles_m
  use mpi_m
  use mpi_lib_m
  use multicomm_m
  use ob_green_m
  use profiling_m
  use simul_box_m
  use smear_m
  use states_dim_m
  use units_m
  use varinfo_m
  use lalg_adv_m

  implicit none

  private

  public ::                           &
    states_t,                         &
    states_init,                      &
    states_look,                      &
    states_densities_init,            &
    states_allocate_wfns,             &
    states_deallocate_wfns,           &
    states_null,                      &
    states_end,                       &
    states_copy,                      &
    states_generate_random,           &
    states_fermi,                     &
    states_eigenvalues_sum,           &
    states_allocate_free_states,      &
    states_deallocate_free_states,    &
    states_write_eigenvalues,         &
    states_write_dos,                 &
    states_write_tpa,                 &
    states_write_bands,               &
    states_write_fermi_energy,        &
    states_spin_channel,              &
    states_dens_accumulate,           &
    states_dens_accumulate_batch,     &
    states_dens_reduce,               &
    states_calc_dens,                 &
    states_calc_tau_jp_gn,            &
    state_is_local,                   &
    states_dump,                      &
    states_distribute_nodes,          &
    states_wfns_memory,               &
    states_freeze_orbitals,           &
    states_total_density,             &
    states_init_green

  public ::                           &
    states_are_complex,               &
    states_are_real

  type states_t
    type(states_dim_t) :: d
    type(modelmb_particle_t) :: modelmbparticles
    integer :: nst                  ! Number of states in each irreducible subspace

    integer :: wfs_type             ! real (M_REAL) or complex (M_CMPLX) wavefunctions
    ! pointers to the wavefunctions 
    logical :: only_userdef_istates ! only use user-defined states as initial states in propagation
    FLOAT, pointer :: dpsi(:,:,:,:) ! dpsi(sys%gr%mesh%np_part, st%d%dim, st%nst, st%d%nik)
    CMPLX, pointer :: zpsi(:,:,:,:) ! zpsi(sys%gr%mesh%np_part, st%d%dim, st%nst, st%d%nik)

    logical            :: open_boundaries
    CMPLX, pointer     :: zphi(:, :, :, :)  ! Free states for open-boundary calculations.
    CMPLX, pointer     :: ob_intf_psi(:, :, :, :, :, :) ! (np, 2, st%d%dim, st%nst, st%d%nik, nleads)
    FLOAT, pointer     :: ob_rho(:, :, :)   ! Density of the lead unit cells.
    FLOAT, pointer     :: ob_eigenval(:, :) ! Eigenvalues of free states.
    type(states_dim_t) :: ob_d              ! Dims. of the unscattered systems.
    integer            :: ob_nst            ! nst of the unscattered systems.
    integer            :: ob_ncs            ! No. of continuum states of open system.
                                            ! ob_ncs = ob_nst*st%ob_d%nik / st%d%nik
    FLOAT, pointer     :: ob_occ(:, :)      ! occupations
    CMPLX, pointer     :: ob_green(:, :, :, :, :, :) ! (np, np, nspin, ncs, nik, nleads) Green`s function of the leads.

    ! used for the user-defined wavefunctions (they are stored as formula strings)
    character(len=1024), pointer :: user_def_states(:,:,:) ! (st%d%dim, st%nst, st%d%nik)

    ! the densities and currents (after all we are doing DFT :)
    FLOAT, pointer :: rho(:,:)      ! rho(gr%mesh%np_part, st%d%nspin)
    FLOAT, pointer :: current(:, :, :)      !   current(gr%mesh%np_part, gr%sb%dim, st%d%nspin)

    logical        :: nlcc          ! do we have non-linear core corrections
    FLOAT, pointer :: rho_core(:)   ! core charge for nl core corrections

    ! It may be required to "freeze" the deepest orbitals during the evolution; the density
    ! of these orbitals is kept in frozen_rho. It is different from rho_core.
    FLOAT, pointer :: frozen_rho(:, :)

    FLOAT, pointer :: eigenval(:,:) ! obviously the eigenvalues
    logical        :: fixed_occ     ! should the occupation numbers be fixed?
    FLOAT, pointer :: occ(:,:)      ! the occupation numbers
    logical        :: fixed_spins   ! In spinors mode, the spin direction is set
                                    ! for the initial (random) orbitals.
    FLOAT, pointer :: spin(:, :, :)

    FLOAT          :: qtot          ! (-) The total charge in the system (used in Fermi)
    FLOAT          :: val_charge    ! valence charge

    type(smear_t)  :: smear         ! smearing of the electronic occupations

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

    logical                     :: np_size            ! whether the states were allocated with size mesh%np instead of size np
  end type states_t

contains

  ! ---------------------------------------------------------
  subroutine states_null(st)
    type(states_t), intent(inout) :: st
    call push_sub('states.states_null')

    nullify(st%dpsi, st%zpsi, st%zphi, st%rho, st%current, st%rho_core, st%frozen_rho, st%eigenval)
    nullify(st%ob_intf_psi, st%ob_rho, st%ob_eigenval, st%ob_occ, st%ob_green)
    nullify(st%occ, st%spin, st%node, st%user_def_states)
    nullify(st%d%kpoints, st%d%kweights)
    nullify(st%st_range, st%st_num)

    ! By default, calculations use real wave-functions
    st%wfs_type = M_REAL

    call modelmb_particles_nullify(st%modelmbparticles)

    call pop_sub()
  end subroutine states_null


  ! ---------------------------------------------------------
  subroutine states_init(st, gr, geo)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo

    FLOAT :: excess_charge
    integer :: nempty, ierr, il
    integer, allocatable :: ob_k(:), ob_st(:), ob_d(:)

    call push_sub('states.states_init')

    call states_null(st)

    !%Variable SpinComponents
    !%Type integer
    !%Default unpolarized
    !%Section States
    !%Description
    !% The calculations may be done in three different ways: spin-restricted (TD)DFT (i.e., doubly
    !% occupied "closed shells"), spin-unrestricted or "spin-polarized" (TD)DFT (i.e. we have two
    !% electronic systes, one with spin up and one with spin down), or making use of two-component
    !% spinors.
    !%Option unpolarized 1
    !% Spin-restricted calculations.
    !%Option polarized 2
    !%Option spin_polarized 2
    !% Spin unrestricted, also known as spin-DFT, SDFT. This mode will double the number of wave
    !% functions necessary for a spin-unpolarized calculation.
    !%Option non_collinear 3
    !%Option spinors 3
    !% The spin-orbitals are two-component spinors. This effectively allows the spin-density to
    !% arrange non-collinearly - i.e. the magnetization vector is allowed to take different
    !% directions in different points.
    !%End
    call loct_parse_int(datasets_check('SpinComponents'), UNPOLARIZED, st%d%ispin)
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
    call loct_parse_float(datasets_check('ExcessCharge'), M_ZERO, excess_charge)


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
    !% use fractional occupational numbers, either fixed from the beginning through
    !% the <tt>Occupations</tt> block or by prescribing
    !% an electronic temperature with <tt>Smearing</tt>.
    !%
    !% Note that this number is unrelated to <tt>CalculationMode == unocc</tt>.
    !%End
    call loct_parse_int(datasets_check('ExtraStates'), 0, nempty)
    if (nempty < 0) then
      write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid ExtraStates"
      message(2) = '(0 <= ExtraStates)'
      call write_fatal(2)
    end if

    ! For non-periodic systems this should just return the Gamma point
    call states_choose_kpoints(st%d, gr%sb, geo)

    call geometry_val_charge(geo, st%val_charge)
    
    if(gr%sb%open_boundaries) then
      excess_charge = -st%val_charge
    end if

    st%qtot = -(st%val_charge + excess_charge)

    nullify(st%ob_intf_psi)
    st%open_boundaries = .false.
    ! When doing open-boundary calculations the number of free states is
    ! determined by the previous periodic calculation.
    st%open_boundaries = gr%sb%open_boundaries
    if(gr%sb%open_boundaries) then
      SAFE_ALLOCATE( ob_k(1:NLEADS))
      SAFE_ALLOCATE(ob_st(1:NLEADS))
      SAFE_ALLOCATE( ob_d(1:NLEADS))
      do il = 1, NLEADS
        call states_look(trim(gr%sb%lead_restart_dir(il))//'/gs', mpi_world, &
          ob_k(il), ob_d(il), ob_st(il), ierr, .true.)
        if(ierr.ne.0) then
          message(1) = 'Could not read the number of states of the periodic calculation'
          message(2) = 'from '//trim(gr%sb%lead_restart_dir(il))//'/gs.'
          call write_fatal(2)
        end if
      end do
      if(ob_k(LEFT).ne.ob_k(RIGHT).or. &
        ob_st(LEFT).ne.ob_st(LEFT).or. &
        ob_d(LEFT).ne.ob_d(RIGHT)) then
        message(1) = 'The number of states for the left and right leads are not equal.'
        call write_fatal(1)
      end if
      st%ob_d%dim = ob_d(LEFT)
      st%ob_nst   = ob_st(LEFT)
      st%ob_d%nik = ob_k(LEFT)
      st%d%nik = st%ob_d%nik
      SAFE_DEALLOCATE_A(ob_d)
      SAFE_DEALLOCATE_A(ob_st)
      SAFE_DEALLOCATE_A(ob_k)
      call distributed_nullify(st%ob_d%kpt, 0)
      if((st%d%ispin.eq.UNPOLARIZED.and.st%ob_d%dim.ne.1) .or.   &
        (st%d%ispin.eq.SPIN_POLARIZED.and.st%ob_d%dim.ne.1) .or. &
        (st%d%ispin.eq.SPINORS.and.st%ob_d%dim.ne.2)) then
        message(1) = 'The spin type of the leads calculation from '//gr%sb%lead_restart_dir(LEFT)
        message(2) = 'and SpinComponents of the current run do not match.'
        call write_fatal(2)
      end if
      ! If the system is spin-polarized one half of the free states
      ! goes to the spin-up k-index, the other half to the spin-down
      ! k-index, we therefore divide by st%d%nik.
      st%ob_ncs = st%ob_d%nik*st%ob_nst / st%d%nik
      st%ob_ncs = 1
      SAFE_DEALLOCATE_P(st%d%kpoints)
      SAFE_DEALLOCATE_P(st%d%kweights)
      SAFE_ALLOCATE( st%d%kpoints(1:MAX_DIM, 1:st%d%nik))
      SAFE_ALLOCATE(st%d%kweights(1:st%d%nik))
      st%d%kpoints  = M_ZERO
      st%d%kweights = M_ZERO
      st%d%kweights(1) = M_ONE
      SAFE_ALLOCATE(st%ob_d%kpoints(1:MAX_DIM, 1:st%ob_d%nik))
      SAFE_ALLOCATE(st%ob_d%kweights(1:st%ob_d%nik))
      SAFE_ALLOCATE(st%ob_eigenval(1:st%ob_nst, 1:st%ob_d%nik))
      SAFE_ALLOCATE(st%ob_occ(1:st%ob_nst, 1:st%ob_d%nik))
      st%ob_d%kpoints  = M_ZERO
      st%ob_d%kweights = M_ZERO
      st%ob_eigenval   = huge(st%ob_eigenval)
      st%ob_occ        = M_ZERO
      call read_ob_eigenval_and_occ()
    else
      st%ob_nst   = 0
      st%ob_ncs   = 0
      st%ob_d%nik = 0
      st%ob_d%dim = 0
    end if

    select case(st%d%ispin)
    case(UNPOLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty + st%ob_ncs
      st%d%nspin = 1
      st%d%spin_channels = 1
    case(SPIN_POLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty + st%ob_ncs
      ! st%d%nik = st%d%nik*2
      st%d%nspin = 2
      st%d%spin_channels = 2
    case(SPINORS)
      st%d%dim = 2
      st%nst = int(st%qtot)
      if(st%nst < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty + st%ob_ncs
      st%d%nspin = 4
      st%d%spin_channels = 2
    end select

    !%Variable StatesSaveMemory
    !%Type logical
    !%Default false
    !%Section Execution::Optimization
    !%Description
    !% (experimental) If set to yes, the wave functions will require
    !% less memory, at the expense of increased computational
    !% time. The default is no.
    !%End
    call loct_parse_logical(datasets_check('StatesSaveMemory'), .false., st%np_size)

    ! FIXME: For now, open-boundary calculations are only possible for
    ! continuum states, i.e. for those states treated by the Lippmann-
    ! Schwinger approach during SCF.
    if(gr%sb%open_boundaries) then
      if(st%nst.ne.st%ob_nst .or. st%d%nik.ne.st%ob_d%nik) then
        message(1) = 'Open-boundary calculations for possibly bound states'
        message(2) = 'are not possible yet. You have to match your number'
        message(3) = 'of states to the number of free states of your previous'
        message(4) = 'periodic run.'
        write(message(5), '(a,i5,a)') 'Your finite system contributes ', st%nst-st%ob_nst, ' states,'
        write(message(6), '(a,i5,a)') 'while your periodic calculation had ', st%ob_nst, ' states.'
        write(message(7), '(a,i5,a)') 'Your finite system contributes ', st%d%nik-st%ob_d%nik, ' k-points,'
        write(message(8), '(a,i5,a)') 'while your periodic calculation had ', st%ob_d%nik, ' k-points.'
        call write_fatal(8)
      end if
    end if

    ! current
    call loct_parse_logical(datasets_check('CurrentDFT'), .false., st%d%cdft)
    if (st%d%cdft) then
      ! Use of CDFT requires complex wave-functions
      st%wfs_type = M_CMPLX

      if(st%d%ispin == SPINORS) then
        message(1) = "Sorry, current DFT not working yet for spinors"
        call write_fatal(1)
      end if
      message(1) = "Info: Using current DFT"
      call write_info(1)
    end if

    ! Periodic systems require complex wave-functions
    ! but not if it is gamma-point only
    if(simul_box_is_periodic(gr%sb)) then
      if(.not. (st%d%nik == 1 .and. kpoint_is_gamma(st%d, 1))) then
        st%wfs_type = M_CMPLX
      endif
    endif

    ! Calculations with open boundaries require complex wavefunctions.
    if(gr%sb%open_boundaries) then
      st%wfs_type = M_CMPLX
    end if

    !%Variable OnlyUserDefinedInitialStates
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% If true, then only user-defined states from the block UserDefinedStates
    !% will be used as initial states for a time propagation. No attempt is made
    !% to load ground-state orbitals from a previous ground-state run.
    !%End
    call loct_parse_logical(datasets_check('OnlyUserDefinedInitialStates'), .false., st%only_userdef_istates)

    ! we now allocate some arrays
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)
    st%occ      = M_ZERO
    ! allocate space for formula strings that define user-defined states
    SAFE_ALLOCATE(st%user_def_states(1:st%d%dim, 1:st%nst, 1:st%d%nik))
    if(st%d%ispin == SPINORS) then
      SAFE_ALLOCATE(st%spin(1:3, 1:st%nst, 1:st%d%nik))
    else
      nullify(st%spin)
    end if

    ! initially we mark all 'formulas' as undefined
    st%user_def_states(1:st%d%dim, 1:st%nst, 1:st%d%nik) = 'undefined'

    call states_read_initial_occs(st, excess_charge)
    call states_read_initial_spins(st)

    nullify(st%zphi)

    st%st_start = 1
    st%st_end = st%nst
    st%lnst = st%nst
    SAFE_ALLOCATE(st%node(1:st%nst))
    st%node(1:st%nst) = 0

    call mpi_grp_init(st%mpi_grp, -1)
    st%parallel_in_states = .false.

    nullify(st%dpsi, st%zpsi)

    call distributed_nullify(st%d%kpt, st%d%nik)

    call modelmb_particles_init (st%modelmbparticles,gr)

    call pop_sub()

  contains

    subroutine read_ob_eigenval_and_occ()
      integer            :: occs, jk, ist, ik, err
      FLOAT              :: flt, eigenval, occ
      character          :: char
      character(len=256) :: restart_dir, line, chars

      call push_sub('states.read_ob_eigenval_and_occ')

      restart_dir = trim(gr%sb%lead_restart_dir(LEFT))//'/gs'

      occs = io_open(trim(restart_dir)//'/occs', action='read', is_tmp=.true., grp=mpi_world)
      if(occs.lt.0) then
        message(1) = 'Could not read '//trim(restart_dir)//'/occs.'
        call write_fatal(1)
      end if

      ! Skip two lines.
      call iopar_read(mpi_world, occs, line, err)
      call iopar_read(mpi_world, occs, line, err)

      jk = 1
      do
        ! Check for end of file.
        call iopar_read(mpi_world, occs, line, err)

        read(line, '(a)') char
        if(char.eq.'%') exit
        call iopar_backspace(mpi_world, occs)

        ! Extract eigenvalue.
        call iopar_read(mpi_world, occs, line, err)
        read(line, *) occ, char, eigenval, char, flt, char, flt, char, flt, char, &
          flt, chars, ik, char, ist
        if(occ.gt.CNST(1e-5)) then
          if(st%d%ispin.eq.SPIN_POLARIZED) then
            if(is_spin_up(ik)) then
!              st%ob_eigenval(jst, SPIN_UP) = eigenval
!              st%ob_occ(jst, SPIN_UP)      = occ
            else
!              st%ob_eigenval(jst, SPIN_DOWN) = eigenval
!              st%ob_occ(jst, SPIN_DOWN)      = occ
            end if
          else
            st%ob_eigenval(1, jk) = eigenval
            st%ob_occ(1, jk)      = occ
          end if
          jk = jk + 1
        end if
      end do

      call io_close(occs)

      call pop_sub()
    end subroutine read_ob_eigenval_and_occ
  end subroutine states_init


  ! ---------------------------------------------------------
  ! Allocate free states.
  subroutine states_allocate_free_states(st, gr)
    type(states_t), intent(inout) :: st
    type(grid_t),   intent(in)    :: gr

    call push_sub('states.states_allocate_free_states')

    ! FIXME: spin-polarized free states ignored.
    if(gr%sb%open_boundaries) then
      SAFE_ALLOCATE(st%zphi(1:gr%mesh%np_part, 1:st%ob_d%dim, 1:st%ob_nst, 1:st%ob_d%nik))
      SAFE_ALLOCATE(st%ob_rho(1:gr%mesh%lead_unit_cell(LEFT)%np, 1:st%d%nspin, 1:NLEADS))
      st%zphi = M_z0
    else
      nullify(st%zphi)
    end if

    call pop_sub()
  end subroutine states_allocate_free_states


  ! ---------------------------------------------------------
  ! Deallocate free states.
  subroutine states_deallocate_free_states(st, gr)
    type(states_t), intent(inout) :: st
    type(grid_t),   intent(in)    :: gr

    call push_sub('states.states_deallocate_free_states')

    if(gr%sb%open_boundaries) then
      SAFE_DEALLOCATE_P(st%zphi)
      SAFE_DEALLOCATE_P(st%ob_rho)
    end if

    call pop_sub()
  end subroutine states_deallocate_free_states


  ! ---------------------------------------------------------
  ! Reads from the input file the initial occupations, if the
  ! block "Occupations" is present. Otherwise, it makes an initial
  ! guess for the occupations, maybe using the "Smearing"
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
    !% would fix the occupations of the five states to <i>2.0</i>. There can be
    !% at most as many columns as states in the calculation. If there are fewer columns
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
    !% an error message is printed. If FromScratch = no, this block is ignored when restart information
    !% is read, and the previous occupations are used.
    !%End

    if(st%open_boundaries) then
      st%fixed_occ = .true.
      st%occ  = st%ob_occ
      st%qtot = M_ZERO
      do i = 1, st%nst
        st%qtot = st%qtot + sum(st%occ(i, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
      end do

    else
      occ_fix: if(loct_parse_block(datasets_check('Occupations'), blk)==0) then
        ! read in occupations
        st%fixed_occ = .true.

        ! Reads the number of columns in the first row. This assumes that all rows
        ! have the same column number; otherwise the code will stop with an error.
        ncols = loct_parse_block_cols(blk, 0)
        if(ncols > st%nst) then
          call input_error("Occupations")
        end if
        ! Now we fill all the "missing" states with the maximum occupation.
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
      end if occ_fix
    end if

    call smear_init(st%smear, st%d%ispin, st%fixed_occ)

    ! sanity check
    r = M_ZERO
    do i = 1, st%nst
      r = r + sum(st%occ(i, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
    end do
    if(abs(r - st%qtot) > CNST(1e-6)) then
      message(1) = "Occupations do not integrate to total charge!"
      write(message(2), '(6x,f12.6,a,f12.6)') r, ' != ', st%qtot
      call write_warning(2)
    end if

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
    spin_fix: if(loct_parse_block(datasets_check('InitialSpins'), blk)==0) then
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
  subroutine states_allocate_wfns(st, mesh, wfs_type)
    type(states_t),    intent(inout) :: st
    type(mesh_t),      intent(in)    :: mesh
    integer, optional, intent(in)    :: wfs_type

    integer :: np, ik, ist, idim, st1, st2, k1, k2, size
    logical :: force

    call push_sub('states.states_allocate_wfns')

    if (present(wfs_type)) then
      ASSERT(wfs_type == M_REAL .or. wfs_type == M_CMPLX)
      st%wfs_type = wfs_type
    end if

    !%Variable ForceComplex
    !%Type logical
    !%Default no
    !%Section Execution::Debug
    !%Description
    !% Normally Octopus determines automatically the type necessary
    !% for the wavefunctions. When set to yes this variable will
    !% force the use of complex wavefunctions. 
    !%
    !% Warning: This variable is designed for testing and
    !% benchmarking and normal users need not use it.
    !%
    !%End
    call loct_parse_logical(datasets_check('ForceComplex'), .false., force)

    if(force) st%wfs_type = M_CMPLX

    st1 = st%st_start; st2 = st%st_end
    k1 = st%d%kpt%start; k2 = st%d%kpt%end
    
    if(st%np_size) then
      size = mesh%np
    else
      size = mesh%np_part
    end if

    if (states_are_real(st)) then
      SAFE_ALLOCATE(st%dpsi(1:size, 1:st%d%dim, st1:st2, k1:k2))

      do ik = k1, k2
        do ist = st1, st2
          do idim = 1, st%d%dim
            st%dpsi(1:size, idim, ist, ik) = M_ZERO
          end do
        end do
      end do

    else
      SAFE_ALLOCATE(st%zpsi(1:size, 1:st%d%dim, st1:st2, k1:k2))

      do ik = k1, k2
        do ist = st1, st2
          do idim = 1, st%d%dim
            st%zpsi(1:size, idim, ist, ik) = M_Z0
           end do
        end do
      end do
    end if

    if(calc_mode_is(CM_TD).and.st%open_boundaries) then
      np = mesh%lead_unit_cell(LEFT)%np
      SAFE_ALLOCATE(st%ob_intf_psi(1:np, 1:2, 1:st%d%dim, st1:st2, k1:k2, 1:NLEADS))

      st%ob_intf_psi = M_z0
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
      SAFE_DEALLOCATE_P(st%dpsi)
    else
      SAFE_DEALLOCATE_P(st%zpsi)
    end if

    if(st%open_boundaries.and.calc_mode_is(CM_TD)) then
      SAFE_DEALLOCATE_P(st%ob_intf_psi)
    end if

    call pop_sub()
  end subroutine states_deallocate_wfns


  ! ---------------------------------------------------------
  subroutine states_densities_init(st, gr, geo, mc)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo
    type(multicomm_t), intent(in)    :: mc

    call push_sub('states.states_densities_init')

    ! allocate arrays for charge and current densities
    SAFE_ALLOCATE(st%rho(1:gr%mesh%np_part, 1:st%d%nspin))
    st%rho  = M_ZERO
    if(st%d%cdft) then
      SAFE_ALLOCATE(st%current(1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:st%d%nspin))
      st%current = M_ZERO
    end if
    st%nlcc = geo%nlcc
    if(st%nlcc) then
      SAFE_ALLOCATE(st%rho_core(1:gr%mesh%np))
      st%rho_core(:) = M_ZERO
    end if

    !%Variable StatesBlockSize
    !%Type integer
    !%Default 4
    !%Section Execution::Optimization
    !%Description
    !% Some routines work over blocks of eigenfunctions, this
    !% generally improves performance at the expense of increased
    !% memory consumption. This variable selects the size of the
    !% blocks to be used, the default size is max(4, 2*nthreads).
    !%End

    ! This has to be here as it requires mc%nthreads that is not available in states_init
    call loct_parse_int(datasets_check('StatesBlockSize'), max(4, 2*mc%nthreads), st%d%block_size)
    if(st%d%block_size < 1) then
      message(1) = "Error: The variable 'StatesBlockSize' must be greater than 0."
      call write_fatal(1)
    end if
    
    st%d%block_size = min(st%d%block_size, st%nst)

    !%Variable StatesWindowSize
    !%Type integer
    !%Default nst
    !%Section Execution::Optimization
    !%Description
    !% (experimental) The orthogonalization of the wave-functions can
    !% be done by windows, this value selects the size of the
    !% windows. The default size is total number of states, that
    !% disables window orthogonalization.
    !%End

    ! This has to be here as it requires mc%nthreads that is not available in states_init
    call loct_parse_int(datasets_check('StatesWindowSize'), st%nst, st%d%window_size)
    if(st%d%block_size < 1) then
      message(1) = "Error: The variable 'StatesBlockSize' must be greater than 0."
      call write_fatal(1)
    end if
    
    st%d%window_size = min(st%d%window_size, st%nst)

    call pop_sub()
  end subroutine states_densities_init


  ! ---------------------------------------------------------
  subroutine states_copy(stout, stin)
    type(states_t), intent(inout) :: stout
    type(states_t), intent(in)    :: stin

    call push_sub('states.states_copy')

    call states_null(stout)

    stout%wfs_type   = stin%wfs_type
    call states_dim_copy(stout%d, stin%d)
    stout%nst        = stin%nst
    stout%qtot       = stin%qtot
    stout%val_charge = stin%val_charge
    call smear_copy(stout%smear, stin%smear)
    stout%parallel_in_states = stin%parallel_in_states
    stout%lnst       = stin%lnst
    stout%st_start   = stin%st_start
    stout%st_end     = stin%st_end
    stout%np_size    = stin%np_size
    call loct_pointer_copy(stout%dpsi, stin%dpsi)
    call loct_pointer_copy(stout%zpsi, stin%zpsi)
    call loct_pointer_copy(stout%user_def_states, stin%user_def_states)

    call loct_pointer_copy(stout%rho, stin%rho)
    call loct_pointer_copy(stout%current, stin%current)
    stout%nlcc = stin%nlcc
    call loct_pointer_copy(stout%rho_core, stin%rho_core)
    call loct_pointer_copy(stout%frozen_rho, stin%frozen_rho)
    call loct_pointer_copy(stout%eigenval, stin%eigenval)
    stout%fixed_occ = stin%fixed_occ
    call loct_pointer_copy(stout%occ, stin%occ)
    stout%fixed_spins = stin%fixed_spins
    call loct_pointer_copy(stout%spin, stin%spin)
    call loct_pointer_copy(stout%node, stin%node)
    call mpi_grp_copy(stout%mpi_grp, stin%mpi_grp)
    call loct_pointer_copy(stout%st_range, stin%st_range)
    call loct_pointer_copy(stout%st_num, stin%st_num)

    call modelmb_particles_copy(stout%modelmbparticles, stin%modelmbparticles)

    if(stin%parallel_in_states) call multicomm_all_pairs_copy(stout%ap, stin%ap)

    stout%open_boundaries = stin%open_boundaries
    ! Some of the "open boundaries" variables are not copied.

    call pop_sub()
  end subroutine states_copy


  ! ---------------------------------------------------------
  subroutine states_end(st)
    type(states_t), intent(inout) :: st

    call push_sub('states.states_end')
    
    SAFE_DEALLOCATE_P(st%dpsi)
    SAFE_DEALLOCATE_P(st%zpsi)
    SAFE_DEALLOCATE_P(st%user_def_states)

    SAFE_DEALLOCATE_P(st%rho)
    SAFE_DEALLOCATE_P(st%current)
    SAFE_DEALLOCATE_P(st%rho_core)
    SAFE_DEALLOCATE_P(st%frozen_rho)
    SAFE_DEALLOCATE_P(st%eigenval)

    SAFE_DEALLOCATE_P(st%occ)
    SAFE_DEALLOCATE_P(st%spin)

    SAFE_DEALLOCATE_P(st%node)
    SAFE_DEALLOCATE_P(st%st_range)
    SAFE_DEALLOCATE_P(st%st_num)

    if(st%parallel_in_states) then
      SAFE_DEALLOCATE_P(st%ap%schedule)
    end if

    SAFE_DEALLOCATE_P(st%zphi)

    call states_dim_end(st%d)
    if(st%open_boundaries) then
      call states_dim_end(st%ob_d)
      SAFE_DEALLOCATE_P(st%ob_eigenval)
      SAFE_DEALLOCATE_P(st%ob_occ)
      SAFE_DEALLOCATE_P(st%ob_green)
    end if

    SAFE_DEALLOCATE_P(st%user_def_states)

    call modelmb_particles_end(st%modelmbparticles)

    call pop_sub()
  end subroutine states_end


  ! ---------------------------------------------------------
  ! Calculates the new density out the wavefunctions and
  ! occupations...
  subroutine states_dens_accumulate(st, np, rho, ist, ik)
    type(states_t), intent(in)    :: st
    integer,        intent(in)    :: np
    FLOAT,          intent(inout) :: rho(:,:)
    integer,        intent(in)    :: ist
    integer,        intent(in)    :: ik
    
    integer :: ip, ispin
    CMPLX   :: c
    type(profile_t), save :: prof

#ifndef USE_OMP
    call push_sub('states.states_dens_accumulate')
#endif
    call profiling_in(prof, "CALC_DENSITY")

    ispin = states_dim_get_spin_index(st%d, ik)
    
    if (st%wfs_type == M_REAL) then
      do ip = 1, np
        rho(ip, ispin) = rho(ip, ispin) + st%d%kweights(ik)*st%occ(ist, ik)*st%dpsi(ip, 1, ist, ik)**2
      end do
    else
      do ip = 1, np
        rho(ip, ispin) = rho(ip, ispin) + st%d%kweights(ik)*st%occ(ist, ik)*&
             (real(st%zpsi(ip, 1, ist, ik), REAL_PRECISION)**2 + aimag(st%zpsi(ip, 1, ist, ik))**2)
      end do
    end if
    
    if(st%d%ispin == SPINORS) then ! in this case wave-functions are always complex
      do ip = 1, np
        rho(ip, 2) = rho(ip, 2) + st%d%kweights(ik)*st%occ(ist, ik)*&
             (real(st%zpsi(ip, 2, ist, ik), REAL_PRECISION)**2 + aimag(st%zpsi(ip, 2, ist, ik))**2)
        
        c = st%d%kweights(ik)*st%occ(ist, ik)*st%zpsi(ip, 1, ist, ik)*conjg(st%zpsi(ip, 2, ist, ik))
        rho(ip, 3) = rho(ip, 3) + real(c, REAL_PRECISION)
        rho(ip, 4) = rho(ip, 4) + aimag(c)
      end do
    end if
    
    call profiling_out(prof)

#ifndef USE_OMP
    call pop_sub()
#endif
  end subroutine states_dens_accumulate


  subroutine states_dens_accumulate_batch(np, rho, st, psib, ik)
    integer,        intent(in)    :: np
    FLOAT,          intent(inout) :: rho(:,:)
    type(states_t), intent(in)    :: st
    type(batch_t),  intent(in)    :: psib
    integer,        intent(in)    :: ik
    
    integer :: ist, ist2, ip, ispin
    CMPLX   :: c
    type(profile_t), save :: prof

#ifndef USE_OMP
    call push_sub('states.states_dens_accumulate_batch')
#endif
    call profiling_in(prof, "CALC_DENSITY")

    ispin = states_dim_get_spin_index(st%d, ik)
    
    if (st%wfs_type == M_REAL) then
      do ist = 1, psib%nst
        ist2 = psib%states(ist)%ist
        forall(ip = 1:np)
          rho(ip, ispin) = rho(ip, ispin) + st%d%kweights(ik)*st%occ(ist2, ik)*&
            psib%states(ist)%dpsi(ip, 1)**2
        end forall
      end do
    else
      do ist = 1, psib%nst
        ist2 = psib%states(ist)%ist
        forall(ip = 1:np)
          rho(ip, ispin) = rho(ip, ispin) + st%d%kweights(ik)*st%occ(ist2, ik)* ( &
            real (psib%states(ist)%zpsi(ip, 1), REAL_PRECISION)**2 + &
            aimag(psib%states(ist)%zpsi(ip, 1))**2)
        end forall
      end do
    end if
    
    if(st%d%ispin == SPINORS) then ! in this case wave-functions are always complex
     do ist = 1, psib%nst
        ist2 = psib%states(ist)%ist
        do ip = 1, np
          rho(ip, 2) = rho(ip, 2) + st%d%kweights(ik)*st%occ(ist2, ik)* ( &
            real (psib%states(ist)%zpsi(ip, 2), REAL_PRECISION)**2 + &
            aimag(psib%states(ist)%zpsi(ip, 2)**2))
        
          c = st%d%kweights(ik)*st%occ(ist, ik)* &
            psib%states(ist)%zpsi(ip, 1)*conjg(psib%states(ist)%zpsi(ip, 2))
          rho(ip, 3) = rho(ip, 3) + real(c, REAL_PRECISION)
          rho(ip, 4) = rho(ip, 4) + aimag(c)
        end do
      end do
    end if
    
    call profiling_out(prof)

#ifndef USE_OMP
    call pop_sub()
#endif
  end subroutine states_dens_accumulate_batch


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
    ! reduce over states
    if(st%parallel_in_states) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
#ifndef HAVE_MPI2
      SAFE_ALLOCATE(reduce_rho(1:np))
#endif
      do ispin = 1, st%d%nspin
#ifndef HAVE_MPI2
        call blas_copy(np, rho(1, ispin), 1, reduce_rho(1), 1)
#endif
        call MPI_Allreduce(MPI_IN_PLACE_OR(reduce_rho(1)), rho(1, ispin), np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      end do
      SAFE_DEALLOCATE_A(reduce_rho)
      call profiling_out(reduce_prof)
    end if

    ! reduce over kpoints
    if(st%d%kpt%parallel) then
      call profiling_in(reduce_prof, "DENSITY_REDUCE")
      SAFE_ALLOCATE(reduce_rho(1:np))
      do ispin = 1, st%d%nspin
        call MPI_Allreduce(rho(1, ispin), reduce_rho(1), np, MPI_FLOAT, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
        call blas_copy(np, reduce_rho(1), 1, rho(1, ispin), 1)
      end do
      SAFE_DEALLOCATE_A(reduce_rho)
      call profiling_out(reduce_prof)
    end if
#endif

    call pop_sub()
  end subroutine states_dens_reduce


  ! ---------------------------------------------------------
  ! Computes the density from the orbitals in st. If rho is
  ! present, the density is placed there; if it is not present,
  ! the density is placed in st%rho.
  ! ---------------------------------------------------------
  subroutine states_calc_dens(st, np, rho)
    type(states_t), intent(in)  :: st
    integer,        intent(in)  :: np
    FLOAT, optional, target, intent(out) :: rho(:,:)

    integer :: ik, ist

    FLOAT, pointer :: dens(:, :)

    call push_sub('states.states_calc_dens')

    if(present(rho)) then
      dens => rho
    else
      dens => st%rho
    end if

    dens(1:np, 1:st%d%nspin) = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        call states_dens_accumulate(st, np, dens, ist, ik)
      end do
    end do

    call states_dens_reduce(st, np, dens)

    nullify(dens)
    call pop_sub()
  end subroutine states_calc_dens

  ! ---------------------------------------------------------
  ! generate a hydrogen s-wavefunction around a random point
  subroutine states_generate_random(st, mesh, ist_start_, ist_end_)
    type(states_t),    intent(inout) :: st
    type(mesh_t),      intent(in)    :: mesh
    integer, optional, intent(in)    :: ist_start_, ist_end_

    integer :: ist, ik, id, ist_start, ist_end, j, seed
    CMPLX   :: alpha, beta

    call push_sub('states.states_generate_random')

    ist_start = st%st_start
    if(present(ist_start_)) ist_start = max(ist_start, ist_start_)
    ist_end = st%st_end
    if(present(ist_end_)) ist_end = min(ist_end, ist_end_)

    if(st%parallel_in_states) then
      seed = st%mpi_grp%rank
    else
      seed = 0
    end if
    if(st%d%kpt%parallel) then
      seed = st%d%kpt%mpi_grp%rank
    else
      seed = 0
    end if

    select case(st%d%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = ist_start, ist_end
          if (states_are_real(st)) then
            call dmf_random(mesh, st%dpsi(:, 1, ist, ik), seed)
          else
            call zmf_random(mesh, st%zpsi(:, 1, ist, ik), seed)
          end if
          st%eigenval(ist, ik) = M_ZERO
        end do
      end do

    case(SPINORS)

      ASSERT(st%wfs_type == M_CMPLX)

      if(st%fixed_spins) then
        
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = ist_start, ist_end
            call zmf_random(mesh, st%zpsi(:, 1, ist, ik))
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
                                       zmf_dotp(mesh, st%zpsi(:, 1, ist, ik), st%zpsi(:, 1, j, ik)) * &
                                       st%zpsi(:, 1, j, ik)
            end do
            st%zpsi(:, 1, ist, ik) = st%zpsi(:, 1, ist, ik) / zmf_nrm2(mesh, st%zpsi(:, 1, ist, ik))
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
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = ist_start, ist_end
            do id = 1, st%d%dim
              call zmf_random(mesh, st%zpsi(:, id, ist, ik))
            end do
            st%eigenval(ist, ik) = M_ZERO
          end do
        end do
      end if

    end select

    call pop_sub()
  end subroutine states_generate_random

  ! ---------------------------------------------------------
  subroutine states_fermi(st, mesh)
    type(states_t), intent(inout) :: st
    type(mesh_t),   intent(in)    :: mesh

    ! Local variables.
    integer            :: ist, ik
    FLOAT              :: charge
#if defined(HAVE_MPI)
    integer            :: j
    integer            :: tmp
    FLOAT, allocatable :: lspin(:, :) ! To exchange spin.
#endif

    call push_sub('states.states_fermi')

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
      message(1) = 'Occupations do not integrate to total charge!'
      write(message(2), '(6x,f12.8,a,f12.8)') charge, '.ne.', st%qtot
      call write_warning(2)
    end if

    if(st%d%ispin == SPINORS) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          if (st%wfs_type == M_REAL) then
            write(message(1),'(a)') 'Internal error in states_fermi'
            call write_fatal(1)
          else
            st%spin(1:3, ist, ik) = state_spin(mesh, st%zpsi(:, :, ist, ik))
          end if
        end do
#if defined(HAVE_MPI)
        if(st%parallel_in_states) then
          SAFE_ALLOCATE(lspin(1:3, 1:st%lnst))
          lspin = st%spin(1:3, st%st_start:st%st_end, ik)
          do j = 1, 3
            call lmpi_gen_allgatherv(st%lnst, lspin(j, :), tmp, st%spin(j, :, ik), st%mpi_grp)
          end do
          SAFE_DEALLOCATE_A(lspin)
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
    do ik = st%d%kpt%start, st%d%kpt%end
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

    if(st%d%kpt%parallel) then
      call MPI_Allreduce(e, s, 1, MPI_FLOAT, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
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
    FLOAT :: occ
    character(len=80) tmp_str(MAX_DIM), cspin

    call push_sub('states.states_write_eigenvalues')

    ns = 1
    if(st%d%nspin == 2) ns = 2

    message(1) = 'Eigenvalues [' // trim(units_abbrev(units_out%energy)) // ']'
    call write_info(1, iunit)
    if (st%d%nik > ns) then
      message(1) = 'Kpoints [' // trim(units_abbrev(unit_one/units_out%length)) //']'
      call write_info(1, iunit)
    end if

    if(.not.mpi_grp_is_root(mpi_world)) then
      call pop_sub(); return
    end if

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then
        write(message(1), '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
          units_from_atomic(unit_one/units_out%length, st%d%kpoints(1, ik)), ',', &
          units_from_atomic(unit_one/units_out%length, st%d%kpoints(2, ik)), ',',            &
          units_from_atomic(unit_one/units_out%length, st%d%kpoints(3, ik)), ')'
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
            occ = M_ZERO
          else
            occ = st%occ(j, ik+is)
          end if

          if(is.eq.0) cspin = 'up'
          if(is.eq.1) cspin = 'dn'
          if(st%d%ispin.eq.UNPOLARIZED.or.st%d%ispin.eq.SPINORS) cspin = '--'

          write(tmp_str(1), '(i4,3x,a2)') j, trim(cspin)
          if(simul_box_is_periodic(sb)) then
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,3x,4f5.2)') &
                units_from_atomic(units_out%energy, st%eigenval(j, ik) - st%smear%e_fermi), occ, st%spin(1:3, j, ik)
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                units_from_atomic(units_out%energy, st%eigenval(j, ik+is)), occ
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik), ')'
            end if
          else
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,5x,f5.2,3x,3f8.4)') &
                units_from_atomic(units_out%energy, st%eigenval(j, ik)), occ, st%spin(1:3, j, ik)
              if(present(error)) write(tmp_str(3), '(a3,es7.1,a1)')'  (', error(j, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                units_from_atomic(units_out%energy, st%eigenval(j, ik+is)), occ
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

    if(.not.mpi_grp_is_root(mpi_world)) return

    call push_sub('states.states_write_bands')

    !%Variable OutputBandsGnuplotMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in Gnuplot-friendly format.
    !%End
    call loct_parse_logical(datasets_check('OutputBandsGnuplotMode'), .true., gnuplot_mode)

    !%Variable OutputBandsGraceMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in Grace-friendly format.
    !%End
    call loct_parse_logical(datasets_check('OutputBandsGraceMode'), .false., grace_mode)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    SAFE_ALLOCATE(iunit(0:ns-1))

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
              st%d%kpoints(1:sb%dim, ik+is),                   & ! unscaled
              st%d%kpoints(1:sb%dim, ik+is)/factor(1:sb%dim), & ! scaled
              st%eigenval(j, ik + is)/units_out%energy%factor
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

    SAFE_DEALLOCATE_A(iunit)

    call pop_sub()
  end subroutine states_write_bands

  ! ---------------------------------------------------------
  subroutine states_write_tpa(dir, gr, st)
    character(len=*), intent(in) :: dir
    type(grid_t),     intent(in) :: gr
    type(states_t),   intent(in) :: st

    type(block_t) :: blk
    integer       :: ncols, icoord, ist, ik, tpa_initialst, tpa_initialk
    integer       :: iunit

    FLOAT, allocatable  :: ff(:)
    CMPLX, allocatable  :: cff(:)
    FLOAT, allocatable  :: osc(:)
    FLOAT               :: transition_energy, osc_strength, dsf

    FLOAT, parameter    :: M_THRESHOLD = CNST(1.0e-6)
    logical             :: use_qvector = .false.
    FLOAT, allocatable  :: qvector(:)

    call push_sub('states.states_write_tpa')

    ! find the orbital with half occupation
    tpa_initialst = -1
    do ist = 1, st%nst
      do ik = 1, st%d%nik
        if (abs(st%occ(ist,ik)-0.5) .lt. M_THRESHOLD) then
          tpa_initialst = ist
          tpa_initialk  = ik
        end if
      end do
    end do

    ! make sure that half occupancy was found
    if(tpa_initialst.eq.-1) then
      if(mpi_grp_is_root(mpi_world)) then

        message(1) = 'No orbital with half occupancy found. TPA output is not written.'
        call write_warning(1)

        call pop_sub()
        return

      end if
    end if

    !%Variable MomentumTransfer
    !%Type block
    !%Section States
    !%Description
    !% Momentum-transfer vector q to be used when calculating matrix elements
    !% <f|exp(iq.r)|i>. This enables the calculation of the dynamic structure factor,
    !% which is closely related to generalized oscillator strengths.
    !% If the vector is not given but TPA output is requested (Output=TPA),
    !% only the oscillator strengths are written in the output file.
    !% For example, to use q = ( 0.1 , 0.2 , 0.3 ), set
    !% <tt>%MomentumTransfer
    !% <br>&nbsp;&nbsp; 0.1 | 0.2 | 0.3
    !% <br>%</tt>
    !%End
    if(loct_parse_block(datasets_check('MomentumTransfer'),blk)==0) then

      ! check if input makes sense
      ncols = loct_parse_block_cols(blk, 0)

      if(ncols .ne. gr%mesh%sb%dim ) then ! wrong size

        if(mpi_grp_is_root(mpi_world)) then
          message(1) = 'Inconsistent size of momentum-transfer vector. It will not be used in the TPA calculation.'
          call write_warning(1)
        end if

      else ! correct size

        use_qvector = .true.
        SAFE_ALLOCATE(qvector(1:gr%mesh%sb%dim))

        do icoord = 1,gr%mesh%sb%dim    !for x,y,z
          call loct_parse_block_float(blk, 0, icoord-1, qvector(icoord))
          qvector(icoord) = qvector(icoord) / units_inp%length%factor 
        end do

      end if

    end if

    ! calculate the matrix elements

    SAFE_ALLOCATE(ff(1:gr%mesh%np))
    if(use_qvector) then
      SAFE_ALLOCATE(cff(1:gr%mesh%np))
    end if
    SAFE_ALLOCATE(osc(1:gr%mesh%sb%dim))

    ! root writes output to file

    if(mpi_grp_is_root(mpi_world)) then

      iunit = io_open(trim(dir)//'/'//trim('tpa_xas'), action='write')    

      ! header
      if(use_qvector) then
        write (message(1),'(a1,a30,3(es14.5,1x),a1)') '#', ' momentum-transfer vector : (', &
                                                     & qvector(:)*units_out%length%factor,')'
        select case(gr%mesh%sb%dim)
          case(1); write(message(2), '(a1,4(a15,1x))') '#', 'E' , '<x>', '<f>', 'S(q,omega)'
          case(2); write(message(2), '(a1,5(a15,1x))') '#', 'E' , '<x>', '<y>', '<f>', 'S(q,omega)'
          case(3); write(message(2), '(a1,6(a15,1x))') '#', 'E' , '<x>', '<y>', '<z>', '<f>', 'S(q,omega)'
        end select
        call write_info(2,iunit)
      else
        select case(gr%mesh%sb%dim)
          case(1); write(message(1), '(a1,3(a15,1x))') '#', 'E' , '<x>', '<f>'
          case(2); write(message(1), '(a1,4(a15,1x))') '#', 'E' , '<x>', '<y>', '<f>'
          case(3); write(message(1), '(a1,5(a15,1x))') '#', 'E' , '<x>', '<y>', '<z>', '<f>'
        end select
        call write_info(1,iunit)
      end if

    end if

    ! loop through every state
    do ist = 1,st%nst

      ! final states are the unoccupied ones
      if (abs(st%occ(ist,tpa_initialk)) .lt. M_THRESHOLD) then

        osc_strength=M_ZERO;
        transition_energy=st%eigenval(ist,tpa_initialk)-st%eigenval(tpa_initialst,tpa_initialk)

        ! dipole matrix elements <f|x|i> etc. -> oscillator strengths
        do icoord=1,gr%mesh%sb%dim    ! for x,y,z

          ff(1:gr%mesh%np) = st%dpsi(1:gr%mesh%np,1,tpa_initialst,tpa_initialk) * &
                       &  gr%mesh%x(1:gr%mesh%np,icoord)                        * &
                       &  st%dpsi(1:gr%mesh%np,1,ist,tpa_initialk)
          osc(icoord)  = dmf_integrate(gr%mesh, ff)
          osc_strength = osc_strength + 2.0/real(gr%mesh%sb%dim)*transition_energy*abs(osc(icoord))**2.0

        end do

        ! matrix elements <f|exp(iq.r)|i> -> dynamic structure factor
        if (use_qvector) then

          cff(1:gr%mesh%np) = cmplx(st%dpsi(1:gr%mesh%np,1,tpa_initialst,tpa_initialk)) * &
                       &   cmplx(st%dpsi(1:gr%mesh%np,1,ist,tpa_initialk))
          do icoord=1,gr%mesh%sb%dim    ! for x,y,z
            cff(1:gr%mesh%np) = cff(1:gr%mesh%np) * exp(M_zI*gr%mesh%x(1:gr%mesh%np,icoord)*qvector(icoord))
          end do

          dsf = abs(zmf_integrate(gr%mesh, cff))**2.0
        end if

        ! write oscillator strengths (+ dynamic structure factor if qvector if given) into file
        if(mpi_grp_is_root(mpi_world)) then

          if(use_qvector) then
            write(message(1), '(1x,6(es15.8,1x))') units_from_atomic(units_out%energy, transition_energy), osc(:), osc_strength, &
                                                   units_from_atomic(unit_one/units_out%energy, dsf)
          else
            write(message(1), '(1x,6(es15.8,1x))') units_from_atomic(units_out%energy, transition_energy), osc(:), osc_strength
          endif

          call write_info(1,iunit)

        end if

      end if

    end do

    ! finally close the file
    if(mpi_grp_is_root(mpi_world)) then
      call io_close(iunit)
    end if

    SAFE_DEALLOCATE_A(ff);
    if(use_qvector) then
      SAFE_DEALLOCATE_A(cff);
    end if
    SAFE_DEALLOCATE_A(osc);
    if (use_qvector) then
      SAFE_DEALLOCATE_A(qvector)
    end if
  
    call pop_sub()
 
  end subroutine states_write_tpa

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
    call loct_parse_float(datasets_check('DOSEnergyMin'), evalmin - eextend, emin)

    !%Variable DOSEnergyMax
    !%Type float
    !%Default 
    !%Section Output
    !%Description
    !% Upper bound for the energy mesh of the DOS
    !%End
    call loct_parse_float(datasets_check('DOSEnergyMax'), evalmax + eextend, emax)

    !%Variable DOSEnergyPoints
    !%Type integer
    !%Default 500
    !%Section Output
    !%Description
    !% Determines how many energy points octopus should use for 
    !% the DOS energy grid
    !%End
    call loct_parse_int(datasets_check('DOSEnergyPoints'), 500, epoints)

    !%Variable DOSGamma
    !%Type float
    !%Default 
    !%Section Output
    !%Description
    !% Determines the width of the Lorentzian which is used to sum 
    !% up the DOS sum
    !%End
    call loct_parse_float(datasets_check('DOSGamma'), &
      CNST(0.008)/units_out%energy%factor, gamma)

    ! spacing for energy mesh
    de = (emax - emin) / (epoints - 1)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    ! space for state dependent DOS
    SAFE_ALLOCATE(dos(1:epoints, 1:st%nst, 0:ns-1))
    SAFE_ALLOCATE(iunit(0:ns-1))    

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

    ! for spin-polarized calculations also output spin-resolved tdos
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

    SAFE_DEALLOCATE_A(iunit)
    SAFE_DEALLOCATE_A(dos)

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

    ! write Fermi energy in a format that can be used together 
    ! with bands.dat
    write(message(1), '(a)') '# Fermi energy in a format compatible with bands-gp.dat'

    write(message(2), '(7f12.6)')          &
      minval(st%d%kpoints(1,:)),           &
      minval(st%d%kpoints(2,:)),           &
      minval(st%d%kpoints(3,:)),           &
      minval(st%d%kpoints(1,:)/factor(1)), &
      minval(st%d%kpoints(2,:)/factor(2)), &
      minval(st%d%kpoints(3,:)/factor(3)), &
      st%smear%e_fermi/scale

    ! Gamma point
    write(message(3), '(7f12.6)')          &
      (M_ZERO, i = 1, 6),                  &
      st%smear%e_fermi/scale

    write(message(4), '(7f12.6)')          &
      maxval(st%d%kpoints(1,:)),           &
      maxval(st%d%kpoints(2,:)),           &
      maxval(st%d%kpoints(3,:)),           &
      maxval(st%d%kpoints(1,:)/factor(1)), &
      maxval(st%d%kpoints(2,:)/factor(2)), &
      maxval(st%d%kpoints(3,:)/factor(3)), &
      st%smear%e_fermi/scale

    call write_info(4, iunit)
    call io_close(iunit)

    ! now we write the same information so that it can be used 
    ! together with total-dos.dat
    iunit = io_open(trim(dir)//'/'//'total-dos-efermi.dat', action='write')    

    write(message(1), '(a)') '# Fermi energy in a format compatible with total-dos.dat'    

    ! this is the maximum that tdos can reach
    maxdos = sum(st%d%kweights) * st%nst

    write(message(2), '(4f12.6)') st%smear%e_fermi/scale, M_ZERO
    write(message(3), '(4f12.6)') st%smear%e_fermi/scale, maxdos

    call write_info(3, iunit)
    call io_close(iunit)

    call pop_sub()
  end subroutine states_write_fermi_energy
  ! ---------------------------------------------------------


  ! -------------------------------------------------------
  integer pure function states_spin_channel(ispin, ik, dim)
    integer, intent(in) :: ispin, ik, dim

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

     SAFE_ALLOCATE(st%st_range(1:2, 0:st%mpi_grp%size-1))
     SAFE_ALLOCATE(st%st_num(0:st%mpi_grp%size-1))

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
  logical function states_are_complex(st) result (wac)
    type(states_t),    intent(in) :: st
    wac = (st%wfs_type == M_CMPLX)
  end function states_are_complex


  ! ---------------------------------------------------------
  logical function states_are_real(st) result (war)
    type(states_t),    intent(in) :: st
    war = (st%wfs_type == M_REAL)
  end function states_are_real


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
    FLOAT, optional, intent(out)   ::  tau(:,:)    ! (gr%mesh%np, st%d%nspin)
    FLOAT, optional, intent(out)   ::   jp(:,:,:)  ! (gr%mesh%np, gr%mesh%sb%dim, st%d%nspin)
    FLOAT, optional, intent(out)   :: grho(:,:,:)  ! (gr%mesh%np, gr%mesh%sb%dim, st%d%nspin)

    CMPLX, allocatable :: wf_psi(:,:), gwf_psi(:,:,:)
    CMPLX   :: c_tmp
    integer :: sp, is, ik, ik_tmp, ist, i_dim, st_dim, ii
    FLOAT   :: ww

#if defined(HAVE_MPI)
    FLOAT, allocatable :: tmp_reduce(:)
    integer :: mpi_err
#endif

    SAFE_ALLOCATE( wf_psi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gwf_psi(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%dim))   

    sp = 1
    if(st%d%ispin == SPIN_POLARIZED) sp = 2

    ASSERT(present( tau).or.present(  jp).or.present(grho))

    if(present( tau))  tau(:,:)   = M_ZERO
    if(present(  jp))   jp(:,:,:) = M_ZERO
    if(present(grho)) grho(:,:,:) = M_ZERO

    do is = 1, sp
      do ik_tmp = st%d%kpt%start, st%d%kpt%end, sp
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
            call zderivatives_grad(gr%der, wf_psi(:,st_dim), gwf_psi(:,:,st_dim))
          end do

          ww = st%d%kweights(ik)*st%occ(ist, ik)

          do i_dim = 1, gr%mesh%sb%dim
            if(present(grho)) &
              grho(1:gr%mesh%np, i_dim, is) = grho(1:gr%mesh%np, i_dim, is) + &
                ww*M_TWO*real(conjg(wf_psi(1:gr%mesh%np, 1))*gwf_psi(1:gr%mesh%np, i_dim, 1))
            if(present(  jp)) &
              jp  (1:gr%mesh%np, i_dim, is) = jp  (1:gr%mesh%np, i_dim, is) + &
                ww*aimag(conjg(wf_psi(1:gr%mesh%np, 1))*gwf_psi(1:gr%mesh%np, i_dim, 1))
            if(present( tau)) then
              tau (1:gr%mesh%np, is)        = tau (1:gr%mesh%np, is)        + &
                ww*abs(gwf_psi(1:gr%mesh%np, i_dim, 1))**2
            end if

            if(st%d%ispin == SPINORS) then
              if(present(grho)) then
                grho(1:gr%mesh%np, i_dim, 2) = grho(1:gr%mesh%np, i_dim, 2) + &
                  ww*M_TWO*real(conjg(wf_psi(1:gr%mesh%np, 2))*gwf_psi(1:gr%mesh%np, i_dim, 2))
                grho(1:gr%mesh%np, i_dim, 3) = grho(1:gr%mesh%np, i_dim, 3) + ww* &
                  real (gwf_psi(1:gr%mesh%np, i_dim, 1)*conjg(wf_psi(1:gr%mesh%np, 2)) + &
                    wf_psi(1:gr%mesh%np, 1)*conjg(gwf_psi(1:gr%mesh%np, i_dim, 2)))
                grho(1:gr%mesh%np, i_dim, 4) = grho(1:gr%mesh%np, i_dim, 4) + ww* &
                  aimag(gwf_psi(1:gr%mesh%np, i_dim, 1)*conjg(wf_psi(1:gr%mesh%np, 2)) + &
                    wf_psi(1:gr%mesh%np, 1)*conjg(gwf_psi(1:gr%mesh%np, i_dim, 2)))
              end if
            
              ! the expression for the paramagnetic current with spinors is
              !     j = ( jp(1)             jp(3) + i jp(4) ) 
              !         (-jp(3) + i jp(4)   jp(2)           )
              if(present(  jp)) then
                jp  (1:gr%mesh%np, i_dim, 2) = jp  (1:gr%mesh%np, i_dim, 2) + &
                  ww*aimag(conjg(wf_psi(1:gr%mesh%np, 2))*gwf_psi(1:gr%mesh%np, i_dim, 2))
                do ii = 1, gr%mesh%np
                  c_tmp = conjg(wf_psi(ii, 1))*gwf_psi(ii, i_dim, 2) - wf_psi(ii, 2)*conjg(gwf_psi(ii, i_dim, 1))
                  jp(ii, i_dim, 3) = jp(ii, i_dim, 3) + ww* real(c_tmp)
                  jp(ii, i_dim, 4) = jp(ii, i_dim, 4) + ww*aimag(c_tmp)
                end do
              end if

              ! the expression for the paramagnetic current with spinors is
              !     t = ( tau(1)              tau(3) + i tau(4) ) 
              !         ( tau(3) - i tau(4)   tau(2)            )
              if(present( tau)) then
                tau (1:gr%mesh%np, 2) = tau (1:gr%mesh%np, 2) + ww*abs(gwf_psi(1:gr%mesh%np, i_dim, 2))**2
                do ii = 1, gr%mesh%np
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

    SAFE_DEALLOCATE_A(wf_psi)
    SAFE_DEALLOCATE_A(gwf_psi)

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      call reduce_all(st%mpi_grp)
    end if
    if(st%d%kpt%parallel) then
      call reduce_all(st%d%kpt%mpi_grp)
    end if

  contains 

    subroutine reduce_all(grp)
      type(mpi_grp_t), intent(in)  :: grp

      call push_sub('states.reduce_all')

      SAFE_ALLOCATE(tmp_reduce(1:gr%mesh%np))

      do is = 1, st%d%nspin
        if(present(tau)) then
          call MPI_Allreduce(tau(1, is), tmp_reduce(1), gr%mesh%np, MPI_FLOAT, MPI_SUM, grp%comm, mpi_err)
          tau(1:gr%mesh%np, is) = tmp_reduce(1:gr%mesh%np)       
        end if

        do i_dim = 1, gr%mesh%sb%dim
          if(present(jp)) then
            call MPI_Allreduce(jp(1, i_dim, is), tmp_reduce(1), gr%mesh%np, MPI_FLOAT, MPI_SUM, &
                 grp%comm, mpi_err)
            jp(1:gr%mesh%np, i_dim, is) = tmp_reduce(1:gr%mesh%np)
          end if

          if(present(grho)) then
            call MPI_Allreduce(grho(1, i_dim, is), tmp_reduce(1), gr%mesh%np, MPI_FLOAT, MPI_SUM, &
                 grp%comm, mpi_err)
            grho(1:gr%mesh%np, i_dim, is) = tmp_reduce(1:gr%mesh%np)
          end if
        end do

      end do
      SAFE_DEALLOCATE_A(tmp_reduce)

      call pop_sub()
    end subroutine reduce_all

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
  subroutine states_look(dir, mpi_grp, kpoints, dim, nst, ierr, only_occupied)
    character(len=*), intent(in)    :: dir
    type(mpi_grp_t),  intent(in)    :: mpi_grp
    integer,          intent(out)   :: kpoints, dim, nst, ierr
    logical, intent(in), optional   :: only_occupied

    character(len=256) :: line
    character(len=12)  :: filename
    character(len=1)   :: char
    integer :: iunit, iunit2, err, i, ist, idim, ik, counter_kpoints
    FLOAT :: occ, eigenval
    logical :: only_occupied_

    call push_sub('states.states_look')

    only_occupied_ = .false.
    if(present(only_occupied)) only_occupied_ = only_occupied
    ierr = 0
    iunit  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(iunit < 0) then
      ierr = -1
      call pop_sub
      return
    end if
    iunit2 = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(iunit2 < 0) then
      call io_close(iunit, grp = mpi_grp)
      ierr = -1
      call pop_sub()
      return
    end if

    ! Skip two lines.
    call iopar_read(mpi_grp, iunit, line, err); call iopar_read(mpi_grp, iunit, line, err)
    call iopar_read(mpi_grp, iunit2, line, err); call iopar_read(mpi_grp, iunit2, line, err)

    kpoints = 1
    dim = 1
    nst = 1
    counter_kpoints = 0

    do
      call iopar_read(mpi_grp, iunit, line, i)
      read(line, '(a)') char
      if(i.ne.0.or.char=='%') exit
      read(line, *) ik, char, ist, char, idim, char, filename
      if(ik > kpoints) kpoints = ik
      if(idim == 2)    dim     = 2
      if(ist>nst)      nst     = ist
      call iopar_read(mpi_grp, iunit2, line, err)
      read(line, *) occ, char, eigenval
      if(occ.gt.CNST(1e-5)) counter_kpoints = counter_kpoints + 1
    end do
    if(only_occupied_) kpoints = counter_kpoints

    call io_close(iunit, grp = mpi_grp)
    call io_close(iunit2, grp = mpi_grp)
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

    integer :: ist, ik
    type(states_t) :: staux

    call push_sub('states.states_freeze_orbitals')

    if(n >= st%nst) then
      write(message(1),'(a)') 'Attempting to freeze a number of orbitals which is larger or equal to'
      write(message(2),'(a)') 'the total number. The program has to stop.'
      call write_fatal(2)
    end if

    if(.not.associated(st%frozen_rho)) then
      SAFE_ALLOCATE(st%frozen_rho(1:gr%mesh%np, 1:st%d%dim))
      st%frozen_rho = M_ZERO
    end if

    st%frozen_rho = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        if(ist > n) cycle
        call states_dens_accumulate(st, gr%mesh%np, st%frozen_rho, ist, ik)
      end do
    end do
    call states_dens_reduce(st, gr%mesh%np, st%frozen_rho)

    call states_copy(staux, st)

    st%nst = st%nst - n

    call states_deallocate_wfns(st)
    call states_distribute_nodes(st, mc)
    call states_allocate_wfns(st, gr%mesh, M_CMPLX)

#if defined(HAVE_MPI) 

    if(staux%parallel_in_states) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = staux%st_start, staux%st_end
          if(ist <= n) cycle
          if(.not.state_is_local(st, ist-n)) then
            call mpi_send(staux%zpsi(1, 1, ist, ik), gr%mesh%np_part*st%d%dim, MPI_CMPLX, staux%node(ist), &
              ist, st%mpi_grp%comm, mpi_err)

            call mpi_recv(st%zpsi(1, 1, ist-n, ik), gr%mesh%np_part*st%d%dim, MPI_CMPLX, st%node(ist-n), &
              ist, st%mpi_grp%comm, mpi_err)
          else
            st%zpsi(:, :, ist-n, ik) = staux%zpsi(:, :, ist, ik)
          end if
   
        end do
      end do
   else
     do ik = st%d%kpt%start, st%d%kpt%end
       do ist = staux%st_start, staux%st_end
         if(ist <= n) cycle
         st%zpsi(:, :, ist-n, ik) = staux%zpsi(:, :, ist, ik)
       end do
     end do
   end if

#else

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%zpsi(:, :, ist, ik) = staux%zpsi(:, :, n + ist, ik)
      end do
    end do

#endif

    SAFE_DEALLOCATE_P(st%occ)
    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_ALLOCATE(st%occ     (1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)
    st%occ      = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        st%occ(ist, ik) = staux%occ(n+ist, ik)
        st%eigenval(ist, ik) = staux%eigenval(n+ist, ik)
      end do
    end do

    call states_end(staux)
    call pop_sub()
  end subroutine states_freeze_orbitals


  ! ---------------------------------------------------------
  ! this routine calculates the total electronic density,
  ! which is the sum of the part coming from the orbitals, the
  ! non-linear core corrections and the frozen orbitals
  subroutine states_total_density(st, mesh, rho)
    type(states_t), intent(in)  :: st
    type(mesh_t),   intent(in)  :: mesh
    FLOAT,          intent(out) :: rho(:,:)

    integer :: is, ip

    call push_sub('states.states_total_density')

    forall(is = 1:st%d%nspin, ip = 1:mesh%np)
      rho(ip, is) = st%rho(ip, is)
    end forall

    if(associated(st%rho_core)) then
      forall(is = 1:st%d%spin_channels, ip = 1:mesh%np)
        rho(ip, is) = rho(ip, is) + st%rho_core(ip)/st%d%nspin
      end forall
    end if

    ! Add, if it exists, the frozen density from the inner orbitals.
    if(associated(st%frozen_rho)) then
      forall(is = 1:st%d%spin_channels, ip = 1:mesh%np)
        rho(ip, is) = rho(ip, is) + st%frozen_rho(ip, is)
      end forall
    end if

    call pop_sub()
  end subroutine states_total_density


  ! ---------------------------------------------------------
  real(8) function states_wfns_memory(st, mesh) result(memory)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: mesh
    
    memory = 0.0_8

    ! orbitals
    memory = memory + REAL_PRECISION*dble(mesh%np_part_global)*st%d%dim*dble(st%nst)*st%d%kpt%nglobal

  end function states_wfns_memory


  ! ---------------------------------------------------------
  ! initialize the surface Green`s functions of the leads
  subroutine states_init_green(st, gr, nspin, d_ispin, diag, offdiag)
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: nspin
    integer,             intent(in)    :: d_ispin
    CMPLX,               intent(in)    :: diag(:, :, :, :)      ! Diagonal block of the lead Hamiltonian.
    CMPLX,               intent(in)    :: offdiag(:, :, :)      ! Off-diagonal block of the lead Hamiltonian.

    character(len=1), allocatable  :: ln(:)
    character(len=2)      :: spin
    character(len=256)    :: fmt, fname_real, fname_imag
    FLOAT                 :: energy
    integer  :: np, ik, ist, il, ispin, s1, s2, k1, k2
    integer  :: green_real, green_imag, irow

    call push_sub('states.states_init_green')

    if(calc_mode_is(CM_GS)) then
      SAFE_ALLOCATE(ln(1:NLEADS))
      np = gr%intf(LEFT)%np
      ln(LEFT)  = 'L'; ln(RIGHT) = 'R'
      ! Calculate Green`s function of the leads.
      ! FIXME: For spinors, this calculation is almost certainly wrong.
      ASSERT(st%ob_nst == st%nst)
      ASSERT(st%ob_d%nik == st%d%nik)
      s1 = st%st_start; s2 = st%st_end
      k1 = st%d%kpt%start; k2 = st%d%kpt%end
      SAFE_ALLOCATE(st%ob_green(1:np, 1:np, 1:nspin, s1:s2, k1:k2, 1:NLEADS))
      call messages_print_stress(stdout, "Lead Green's functions")
      message(1) = ' st#     k#  Spin  Lead     Energy'
      call write_info(1)
#ifdef HAVE_MPI 
      ! wait for all processors to finish 
      if(st%d%kpt%parallel) then 
        call MPI_Barrier(st%d%kpt%mpi_grp%comm, mpi_err) 
      end if 
#endif
      do ik = k1, k2
        do ist = s1, s2
          energy = st%ob_eigenval(ist, ik)
          do il = 1, NLEADS
            do ispin = 1, nspin
              select case(d_ispin)
              case(UNPOLARIZED)
                spin = '--'
              case(SPIN_POLARIZED)
                if(is_spin_up(ik)) then
                  spin = 'up'
                else
                  spin = 'dn'
                end if
                ! This is nonsense, but at least all indices are present.
              case(SPINORS)
                if(ispin.eq.1) then
                  spin = 'up'
                else
                  spin = 'dn'
                end if
              end select
              write(message(1), '(i4,3x,i4,3x,a2,5x,a1,1x,f12.6)') ist, ik, spin, ln(il), energy
              call write_info(1)
              ! TODO magnetic gs
              call lead_green(energy, diag(:, :, ispin, il), offdiag(:, :, il), &
                  np, st%ob_green(:, :, ispin, ist, ik, il), .true.)

              ! Write the entire Green`s function to a file.
              if(in_debug_mode) then
                call io_mkdir('debug/open_boundaries')
                write(fname_real, '(3a,i4.4,a,i3.3,a,i1.1,a)') 'debug/open_boundaries/green-', &
                  trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', ispin, '.real'
                write(fname_imag, '(3a,i4.4,a,i3.3,a,i1.1,a)') 'debug/open_boundaries/green-', &
                  trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', ispin, '.imag'
                green_real = io_open(fname_real, action='write', grp=st%d%kpt%mpi_grp, is_tmp=.false.)
                green_imag = io_open(fname_imag, action='write', grp=st%d%kpt%mpi_grp, is_tmp=.false.)

                write(fmt, '(a,i6,a)') '(', np, 'e14.4)'
                do irow = 1, np
                  write(green_real, fmt) real(st%ob_green(irow, :, ispin, ist, ik, il))
                  write(green_imag, fmt) aimag(st%ob_green(irow, :, ispin, ist, ik, il))
                end do
                call io_close(green_real); call io_close(green_imag)
              end if
            end do
          end do
        end do
      end do
      call messages_print_stress(stdout)
      SAFE_DEALLOCATE_A(ln)
    end if

    call pop_sub()
  end subroutine states_init_green

end module states_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
