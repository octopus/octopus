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
!! $Id$

#include "global.h"

module td_m
  use batch_m
  use calc_mode_m
  use cmplxscl_m
  use datasets_m
  use density_m
  use energy_calc_m
  use epot_m
  use forces_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use ground_state_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use ion_dynamics_m
  use kick_m
  use lasers_m
  use lalg_basic_m
  use loct_m
  use loct_math_m
  use math_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use PES_m
  use profiling_m
  use projector_m
  use restart_m
  use scf_m
  use simul_box_m
  use species_m
  use spectrum_m
  use states_m
  use states_calc_m
  use states_dim_m
  use states_io_m
  use states_restart_m
  use system_m
  use propagator_m
  use td_write_m
  use types_m
  use unit_m
  use unit_system_m
  use v_ks_m
  use varinfo_m

  implicit none

  private
  public ::               &
    td_t,                 &
    td_run,               &
    td_run_init,          &
    td_init,              &
    td_end

  !> Parameters.
  integer, parameter :: &
       EHRENFEST = 1,   &
       BO        = 2

  type td_t
    type(propagator_t)   :: tr             !< contains the details of the time-evolution
    type(scf_t)          :: scf
    type(ion_dynamics_t) :: ions
    FLOAT                :: dt             !< time step
    integer              :: max_iter       !< maximum number of iterations to perform
    integer              :: iter           !< the actual iteration
    logical              :: recalculate_gs !< Recalculate ground-state along the evolution.

    type(pes_t)          :: pesv
    type(gauge_force_t)  :: gauge_force

    FLOAT                :: mu
    integer              :: dynamics
    integer              :: energy_update_iter
  end type td_t


contains

  subroutine td_run_init()

    PUSH_SUB(td_run_init)

    call calc_mode_set_parallelization(P_STRATEGY_STATES, default = .true.)

    POP_SUB(td_run_init)
  end subroutine td_run_init

  ! ---------------------------------------------------------
  subroutine td_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch

    type(td_t)                :: td
    type(td_write_t)          :: write_handler
    type(grid_t),     pointer :: gr   ! some shortcuts
    type(states_t),   pointer :: st
    type(geometry_t), pointer :: geo
    logical                   :: stopping, cmplxscl
    integer                   :: iter, ierr, scsteps
    real(8)                   :: etime
    type(profile_t),     save :: prof
    type(restart_t)           :: restart_load, restart_dump

    PUSH_SUB(td_run)

    cmplxscl = hm%cmplxscl%space

    ! some shortcuts
    gr  => sys%gr
    geo => sys%geo
    st  => sys%st

    if(simul_box_is_periodic(gr%mesh%sb)) call messages_experimental('Time propagation for periodic systems')

    call td_init(td, sys, hm)

    ! Allocate wavefunctions during time-propagation
    if(td%dynamics == EHRENFEST) then
      !complex wfs are required for Ehrenfest
      call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX, alloc_Left = cmplxscl)
      if(st%open_boundaries) then
        ASSERT(associated(gr%ob_grid%lead))
        call states_allocate_intf_wfns(st, gr%ob_grid%lead(:)%mesh)
      end if
    else
      call states_allocate_wfns(st, gr%mesh, alloc_Left = cmplxscl)
    end if

    ! Calculate initial value of the gauge vector field
    call gauge_field_init(hm%ep%gfield, gr%sb)

    call init_wfs()

    if(td%iter >= td%max_iter) then
      call end_()
      POP_SUB(td_run)
      return
    endif

    ! Calculate initial forces and kinetic energy
    if(ion_dynamics_ions_move(td%ions)) then
      if(td%iter > 0) then
        call td_read_coordinates()
        call hamiltonian_epot_generate(hm, gr, geo, st, time = td%iter*td%dt)
      end if

      call forces_calculate(gr, geo, hm, st, td%iter*td%dt, td%dt)

      geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)
    end if

    call td_write_init(write_handler, gr, st, hm, geo, &
         ion_dynamics_ions_move(td%ions), gauge_field_is_applied(hm%ep%gfield), hm%ep%kick, td%iter, td%max_iter, td%dt)

    if(td%iter == 0) call td_run_zero_iter()

    if (gauge_field_is_applied(hm%ep%gfield)) call gauge_field_get_force(gr, geo, hm%ep%proj, hm%phase, st, td%gauge_force)

    !call td_check_trotter(td, sys, h)
    td%iter = td%iter + 1

    call restart_init(restart_dump, RESTART_TD, RESTART_TYPE_DUMP, st%dom_st_kpt_mpi_grp, mesh=gr%mesh)
    if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) then
      ! We will also use the TD restart directory as temporary storage during the time propagation
      call restart_init(restart_load, RESTART_TD, RESTART_TYPE_DUMP, st%dom_st_kpt_mpi_grp, &
                        mesh=gr%mesh)
    end if

    call messages_print_stress(stdout, "Time-Dependent Simulation")
    call print_header()

    if(td%pesv%calc_rc .or. td%pesv%calc_mask .and. fromScratch) then
      call pes_init_write(td%pesv,gr%mesh,st)
    end if

    if(st%d%pack_states .and. hamiltonian_apply_packed(hm, gr%mesh)) call states_pack(st)

    etime = loct_clock()
    ! This is the time-propagation loop. It starts at t=0 and finishes
    ! at td%max_iter*dt. The index i runs from 1 to td%max_iter, and
    ! step "iter" means propagation from (iter-1)*dt to iter*dt.
    propagation: do iter = td%iter, td%max_iter

      stopping = clean_stop(sys%mc%master_comm)
      call profiling_in(prof, "TIME_STEP")

      if(iter > 1) then
        if( ((iter-1)*td%dt <= hm%ep%kick%time) .and. (iter*td%dt > hm%ep%kick%time) ) then
          if(.not. cmplxscl) then
            call kick_apply(gr, st, td%ions, geo, hm%ep%kick)
          else
            call kick_apply(gr, st, td%ions, geo, hm%ep%kick, hm%cmplxscl%theta)
          end if
          call td_write_kick(gr, hm, sys%outp, geo, iter)
        end if
      end if

      !Apply mask absorbing boundaries
      if(hm%ab == MASK_ABSORBING) call zvmask(gr, hm, st) 

      ! time iterate the system, one time step.
      select case(td%dynamics)
      case(EHRENFEST)
        call propagator_dt(sys%ks, hm, gr, st, td%tr, iter*td%dt, td%dt, td%mu, td%max_iter, iter, td%ions, geo, &
          gauge_force = td%gauge_force, scsteps = scsteps, &
          update_energy = (mod(iter, td%energy_update_iter) == 0) .or. (iter == td%max_iter) )
      case(BO)
        call propagator_dt_bo(td%scf, gr, sys%ks, st, hm, td%gauge_force, geo, sys%mc, sys%outp, iter, td%dt, td%ions, scsteps)
      end select

      !Photoelectron stuff 
      if(td%pesv%calc_rc .or. td%pesv%calc_mask ) &
        call pes_calc(td%pesv, gr%mesh, st, mod(iter, sys%outp%output_interval), td%dt, iter)

      call td_write_iter(write_handler, gr, st, hm, geo, hm%ep%kick, td%dt, iter)

      ! write down data
      call check_point()

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()

      call profiling_out(prof)
      if (stopping) exit

    end do propagation

    if(st%d%pack_states .and. hamiltonian_apply_packed(hm, gr%mesh)) call states_unpack(st)
    
    call restart_end(restart_dump)
    if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) call restart_end(restart_load)
    call td_write_end(write_handler)
    call end_()

#ifdef HAVE_MPI
    ! wait for all processors to finish
    if(st%parallel_in_states) then
      call MPI_Barrier(st%mpi_grp%comm, mpi_err)
    end if
#endif

    POP_SUB(td_run)

  contains

    subroutine print_header()

      if(.not.cmplxscl) then
        write(message(1), '(a7,1x,a14,a14,a10,a17)') 'Iter ', 'Time ', 'Energy ', 'SC Steps', 'Elapsed Time '
      else
        write(message(1), '(a7,1x,a14,a14,a14,a10,a17)') &
          'Iter ', 'Time ', 'Re(Energy) ','Im(Energy) ', 'SC Steps', 'Elapsed Time '
      end if

      call messages_info(1)
      call messages_print_stress(stdout)

    end subroutine print_header

    ! ---------------------------------------------------------
    subroutine check_point()
      PUSH_SUB(td_run.check_point)

      if(.not. cmplxscl) then
        write(message(1), '(i7,1x,2f14.6,i10,f14.3)') iter, &
          units_from_atomic(units_out%time, iter*td%dt), &
          units_from_atomic(units_out%energy, hm%energy%total + geo%kinetic_energy), &
          scsteps, loct_clock() - etime
      else
        write(message(1), '(i7,1x,3f14.6,i10,f14.3)') iter, &
          units_from_atomic(units_out%time, iter*td%dt), &
          units_from_atomic(units_out%energy, hm%energy%total + geo%kinetic_energy), &
          units_from_atomic(units_out%energy, hm%energy%Imtotal), &
          scsteps, loct_clock() - etime
      end if

      call messages_info(1)
      etime = loct_clock()

      if((sys%outp%output_interval > 0 .and. mod(iter, sys%outp%output_interval) == 0) .or. &
         iter == td%max_iter .or. stopping) then ! output
        call td_write_data(write_handler, gr, st, hm, sys%ks, sys%outp, geo, iter, td%dt)
      end if

      if (mod(iter, sys%outp%restart_write_interval) == 0 .or. iter == td%max_iter .or. stopping) then ! restart
        !if(iter == td%max_iter) sys%outp%iter = ii - 1
        call td_dump(restart_dump, gr, st, hm, td, iter, ierr)
        if (ierr /= 0) then
          message(1) = "Unable to write time-dependent restart information."
          call messages_warning(1)
        end if

        call pes_output(td%pesv, gr%mesh, st, iter, sys%outp, td%dt, gr, geo)

        if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) then
          call messages_print_stress(stdout, 'Recalculating the ground state.')
          fromScratch = .false.
          call ground_state_run(sys, hm, fromScratch)
          call states_load(restart_load, st, gr, ierr, iter=iter)
          if (ierr /= 0) then
            message(1) = "Unable to load TD states."
            call messages_fatal(1)
          end if
          call messages_print_stress(stdout, "Time-dependent simulation proceeds")
          call print_header()
        end if
      end if

      POP_SUB(td_run.check_point)
    end subroutine check_point

   ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(td_run.end_)

      ! free memory
      call states_deallocate_wfns(st)
      call ion_dynamics_end(td%ions)
      call td_end(td)
      if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) call restart_end(restart_load)

      POP_SUB(td_run.end_)
    end subroutine end_

    ! ---------------------------------------------------------
    subroutine init_wfs()

      integer :: ierr, ist, jst, freeze_orbitals
      FLOAT :: x
      type(block_t) :: blk
      type(states_t) :: stin
      CMPLX, allocatable :: rotation_matrix(:, :)
      logical :: freeze_hxc
      type(restart_t) :: restart

      PUSH_SUB(td_run.init_wfs)

      if (.not. fromscratch) then
        call restart_init(restart, RESTART_TD, RESTART_TYPE_LOAD, st%dom_st_kpt_mpi_grp, mesh=gr%mesh)
        call td_load(restart, gr, st, hm, td, ierr)
        if(ierr /= 0) then
          fromScratch = .true.
          td%iter = 0
          message(1) = "Unable to read time-dependent restart information: Starting from scratch"
          call messages_warning(1)
        end if
        call restart_end(restart)
      end if

      if (td%iter >= td%max_iter) then
        message(1) = "All requested iterations have already been done. Use FromScratch = yes if you want to redo them."
        call messages_info(1)
        POP_SUB(td_run.init_wfs)
        return
      end if

      if (fromScratch) then
        call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, st%dom_st_kpt_mpi_grp, mesh=gr%mesh, exact=.true.)

        if(.not. st%only_userdef_istates) then
          call states_load(restart, st, gr, ierr, label = ": gs", kpts_dont_matter = st%open_boundaries)
          if (ierr /= 0) then
            message(1) = 'Unable to read ground-state wavefunctions.'
            call messages_fatal(1)
          end if
          ! extract the interface wave function
          if(st%open_boundaries) call states_get_ob_intf(st, gr)
        end if

        ! check if we should deploy user-defined wavefunctions.
        ! according to the settings in the input file the routine
        ! overwrites orbitals that were read from restart/gs
        if(parse_isdef(datasets_check('UserDefinedStates')) /= 0) then
          call states_read_user_def_orbitals(gr%mesh, st)
        end if
        

        !%Variable TransformStates
        !%Type block
        !%Default no
        !%Section States
        !%Description
        !% Before starting the <tt>td</tt> calculation, the initial states (that are
        !% read from the <tt>restart/gs</tt> directory, which should have been
        !% generated in a previous ground-state calculation) can be "transformed"
        !% among themselves. The block <tt>TransformStates</tt> gives the transformation matrix
        !% to be used. The number of rows of the matrix should equal the number
        !% of the states present in the time-dependent calculation (the independent
        !% spin and <i>k</i>-point subspaces are all transformed equally); the number of
        !% columns should be equal to the number of states present in the
        !% <tt>restart/gs</tt> directory. This number may be different: for example,
        !% one could have run previously in <tt>unocc</tt> mode in order to obtain unoccupied
        !% Kohn-Sham states, and therefore <tt>restart/gs</tt> will contain more states.
        !% These states can be used in the transformation.
        !%
        !% Note that the code will not check the orthonormality of the new states!
        !%
        !% Each line provides the coefficients of the new states, in terms of
        !% the old ones. The coefficients are complex, but the imaginary part will be
        !% ignored for real wavefunctions.
        !%End
        if(parse_isdef(datasets_check('TransformStates')) /= 0) then
          if(parse_block(datasets_check('TransformStates'), blk) == 0) then
            call states_copy(stin, st)
            SAFE_DEALLOCATE_P(stin%zpsi)
            call states_look_and_load(restart, stin, gr)
            ! FIXME: rotation matrix should be R_TYPE
            SAFE_ALLOCATE(rotation_matrix(1:st%nst, 1:stin%nst))
            rotation_matrix = M_z0
            do ist = 1, st%nst
              do jst = 1, parse_block_cols(blk, ist-1)
                call parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
              end do
            end do
            if(states_are_real(st)) then
              call dstates_rotate(gr%mesh, st, stin, real(rotation_matrix, REAL_PRECISION))
            else
              call zstates_rotate(gr%mesh, st, stin, rotation_matrix)
            endif
            SAFE_DEALLOCATE_A(rotation_matrix)
            call states_end(stin)
          else
            message(1) = '"TransformStates" has to be specified as block.'
            call messages_info(1)
            call input_error('TransformStates')
          end if
        end if

        call restart_end(restart)
      end if



      !%Variable TDFreezeOrbitals
      !%Type integer
      !%Default 0
      !%Section Time-Dependent
      !%Description
      !% You have the possibility of "freezing" a number of orbitals during a time-propagation.
      !% The Hartree and exchange-correlation potential due to these orbitals (which
      !% will be the lowest-energy ones) will be added during the propagation, but the orbitals
      !% will not be propagated.
      !%
      !% <b>WARNING: NOT TESTED YET.</b>
      !%Option sae -1
      !% Single-active-electron approximation. This option is only valid for time-dependent
      !% calculations (<tt>CalculationMode = td</tt>). Also, the nuclei should not move.
      !% The idea is that all orbitals except the last one are frozen. The orbitals are to
      !% be read from a previous ground-state calculation. The active orbital is then treated
      !% as independent (whether it contains one electron or two) -- although it will
      !% feel the Hartree and exchange-correlation potentials from the ground-state electronic
      !% configuration.
      !%
      !% It is almost equivalent to setting <tt>TDFreezeOrbitals = N-1</tt>, where <tt>N</tt> is the number
      !% of orbitals, but not completely.
      !%End
      call parse_integer(datasets_check('TDFreezeOrbitals'), 0, freeze_orbitals)

      if(freeze_orbitals /= 0) call messages_experimental('TDFreezeOrbitals')

      if(.not. cmplxscl) then
        call density_calc(st, gr, st%rho)
      else
        call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
      end if

      if(freeze_orbitals > 0) then
        ! In this case, we first freeze the orbitals, then calculate the Hxc potential.
        call states_freeze_orbitals(st, gr, sys%mc, freeze_orbitals)
        write(message(1),'(a,i4,a,i4,a)') 'Info: The lowest', freeze_orbitals, &
          ' orbitals have been frozen.', st%nst, ' will be propagated.'
        call messages_info(1)
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
      elseif(freeze_orbitals < 0) then
        ! This means SAE approximation. We calculate the Hxc first, then freeze all
        ! orbitals minus one.
        write(message(1),'(a)') 'Info: The single-active-electron approximation will be used.'
        call messages_info(1)
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
        call states_freeze_orbitals(st, gr, sys%mc, n = st%nst-1)
        call v_ks_freeze_hxc(sys%ks)
        call density_calc(st, gr, st%rho)
      else
        ! Normal run.
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
      end if

      !%Variable TDFreezeHXC
      !%Type logical
      !%Default no
      !%Section Time-Dependent
      !%Description
      !% The electrons are evolved as independent particles feeling the Hartree and 
      !% exchange-correlation potentials from the ground-state electronic configuration.
      !%End
      call parse_logical(datasets_check('TDFreezeHXC'), .false., freeze_hxc)
      if(freeze_hxc) then 
        write(message(1),'(a)') 'Info: Freezing Hartree and exchange-correlation potentials.'
        call messages_info(1)
        call v_ks_freeze_hxc(sys%ks)
      end if

      x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
      if(st%parallel_in_states) then
        call MPI_Bcast(x, 1, MPI_FLOAT, 0, st%mpi_grp%comm, mpi_err)
      end if
#endif
      call hamiltonian_span(hm, minval(gr%mesh%spacing(1:gr%mesh%sb%dim)), x)
      ! initialize Fermi energy
      call states_fermi(st, gr%mesh)
      call energy_calc_total(hm, gr, st)

      POP_SUB(td_run.init_wfs)
    end subroutine init_wfs


    ! ---------------------------------------------------------
    subroutine td_run_zero_iter()
      PUSH_SUB(td_run.td_run_zero_iter)

      if (gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_init_vec_pot(hm%ep%gfield, gr%sb, st)
        call hamiltonian_update(hm, gr%mesh, time = td%dt*td%iter)
      end if

      call td_write_iter(write_handler, gr, st, hm, geo, hm%ep%kick, td%dt, 0)

      ! I apply the delta electric field *after* td_write_iter, otherwise the
      ! dipole matrix elements in write_proj are wrong
      if(hm%ep%kick%time  ==  M_ZERO) then
        if(.not. cmplxscl) then
          call kick_apply(gr, st, td%ions, geo, hm%ep%kick)
        else
          call kick_apply(gr, st, td%ions, geo, hm%ep%kick, hm%cmplxscl%theta)
        end if
        call td_write_kick(gr, hm, sys%outp, geo, 0)
      end if
      call propagator_run_zero_iter(hm, gr, td%tr)
      if (sys%outp%output_interval > 0) then
        call td_write_data(write_handler, gr, st, hm, sys%ks, sys%outp, geo, 0)
      end if

      POP_SUB(td_run.td_run_zero_iter)
    end subroutine td_run_zero_iter


    ! ---------------------------------------------------------
    !> reads the pos and vel from coordinates file
    subroutine td_read_coordinates() 
      integer :: iatom, iter, iunit
      PUSH_SUB(td_run.td_read_coordinates)

      call io_assign(iunit)
      open(unit = iunit, file = io_workpath('td.general/coordinates'), action='read', status='old')

      if(iunit < 0) then
        message(1) = "Could not open file '"//trim(io_workpath('td.general/coordinates'))//"'."
        message(2) = "Starting simulation from initial geometry."
        call messages_warning(2)
        POP_SUB(td_run.td_read_coordinates)
        return
      end if

      call io_skip_header(iunit)
      do iter = 0, td%iter - 1
        read(iunit, *) ! skip previous iterations... sorry, but no portable seek in Fortran
      end do
      read(iunit, '(28x)', advance='no') ! skip the time index.

      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%x(1:gr%mesh%sb%dim)
        geo%atom(iatom)%x(:) = units_to_atomic(units_inp%length, geo%atom(iatom)%x(:))
      end do
      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%v(1:gr%mesh%sb%dim)
        geo%atom(iatom)%v(:) = units_to_atomic(units_inp%velocity, geo%atom(iatom)%v(:))
      end do
      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%f(1:gr%mesh%sb%dim)
        geo%atom(iatom)%f(:) = units_to_atomic(units_inp%force, geo%atom(iatom)%f(:))
      end do

      call io_close(iunit)
      POP_SUB(td_run.td_read_coordinates)
    end subroutine td_read_coordinates

  end subroutine td_run


  ! ---------------------------------------------------------
  subroutine td_dump(restart, gr, st, hm, td, iter, ierr)
    type(restart_t),     intent(in)  :: restart
    type(grid_t),        intent(in)  :: gr
    type(states_t),      intent(in)  :: st
    type(hamiltonian_t), intent(in)  :: hm
    type(td_t),          intent(in)  :: td
    integer,             intent(in)  :: iter
    integer,             intent(out) :: ierr

    logical :: cmplxscl
    integer :: err, err2

    PUSH_SUB(td_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(td_dump)
      return
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Writing td restart."
      call messages_info(1)
    end if

    ! first write resume file
    call states_dump(restart, st, gr, err, iter=iter)
    if (err /= 0) ierr = ierr + 1

    cmplxscl = st%cmplxscl%space      
    call vksinterp_dump(td%tr%vksold, restart, gr, cmplxscl, st%d%nspin, err2)
    if (err2 /= 0) ierr = ierr + 2

    call pes_dump(restart, td%pesv, st, err)
    if (err /= 0) ierr = ierr + 4

    ! Gauge field restart
    if (gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_dump(restart, hm%ep%gfield, ierr)
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Writing td restart done."
      call messages_info(1)
    end if

    POP_SUB(td_dump)
  end subroutine td_dump

  ! ---------------------------------------------------------
  subroutine td_load(restart, gr, st, hm, td, ierr)
    type(restart_t),     intent(in)    :: restart
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: hm
    type(td_t),          intent(inout) :: td
    integer,             intent(out)   :: ierr

    logical :: cmplxscl
    integer :: err, err2
    PUSH_SUB(td_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(td_load)
      return
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Reading td restart."
      call messages_info(1)
    end if

    ! Read states
    call states_load(restart, st, gr, err, iter=td%iter, read_left = st%have_left_states, label = ": td")
    if (err /= 0) then
      ierr = ierr + 1
    else if (st%open_boundaries) then
      ! extract the interface wave function
      call states_get_ob_intf(st, gr)
    end if

    ! read potential from previous interactions
    cmplxscl = st%cmplxscl%space
    call vksinterp_load(td%tr%vksold, restart, gr, cmplxscl, st%d%nspin, err2)

    if (err2 /= 0) ierr = ierr + 2

    ! read PES restart
    if (td%pesv%calc_rc .or. td%pesv%calc_mask) then
      call pes_load(restart, td%pesv, st, err)
      if (err /= 0) ierr = ierr + 4
    end if

    ! Gauge field restart
    if (gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_load(restart, hm%ep%gfield, err)
      if (err /= 0) then
        ierr = ierr + 8
      else
        call hamiltonian_update(hm, gr%mesh, time = td%dt*td%iter)
      end if
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Reading td restart done."
      call messages_info(1)
    end if

    POP_SUB(td_load)
  end subroutine td_load


#include "td_init_inc.F90"

end module td_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
