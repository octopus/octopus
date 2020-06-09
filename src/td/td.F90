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

module td_oct_m
  use boundaries_oct_m
  use boundary_op_oct_m
  use calc_mode_par_oct_m
  use current_oct_m
  use classical_particle_oct_m
  use density_oct_m
  use energy_calc_oct_m
  use epot_oct_m
  use forces_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use ground_state_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use kick_oct_m
  use lasers_oct_m
  use lda_u_oct_m
  use lda_u_io_oct_m
  use linked_list_oct_m
  use loct_oct_m
  use maxwell_boundary_op_oct_m
  use messages_oct_m
  use modelmb_exchange_syms_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pes_oct_m
  use poisson_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_oct_m
  use propagator_base_oct_m
  use restart_oct_m
  use scdm_oct_m
  use scf_oct_m
  use scissor_oct_m
  use simul_box_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_restart_oct_m
  use system_oct_m
  use system_abst_oct_m
  use td_write_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use varinfo_oct_m
  use walltimer_oct_m
  use xc_oct_m

  implicit none

  private
  public ::               &
    td_t,                 &
    td_run,               &
    td_run_init,          &
    td_init,              &
    td_end,               &
    multisys_td_run,      &
    transform_states

  !> Parameters.
  integer, parameter :: &
    EHRENFEST = 1,   &
    BO        = 2

  type td_t
    private
    type(propagator_t),   public :: tr             !< contains the details of the time-evolution
    type(scf_t)                  :: scf
    type(ion_dynamics_t), public :: ions
    FLOAT,                public :: dt             !< time step
    integer,              public :: max_iter       !< maximum number of iterations to perform
    integer,              public :: iter           !< the actual iteration
    logical                      :: recalculate_gs !< Recalculate ground-state along the evolution.

    type(pes_t)                  :: pesv

    FLOAT,                public :: mu
    integer                      :: dynamics
    integer                      :: energy_update_iter
    FLOAT                        :: scissor

    logical                      :: freeze_occ
    logical                      :: freeze_u
  end type td_t


contains

  subroutine td_run_init()

    PUSH_SUB(td_run_init)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .true.)

    POP_SUB(td_run_init)
  end subroutine td_run_init

  ! ---------------------------------------------------------

  subroutine td_init(td, sys)
    type(td_t),       intent(inout) :: td
    type(system_t),   intent(inout) :: sys

    integer :: default
    FLOAT   :: spacing, default_dt, propagation_time

    PUSH_SUB(td_init)

    if (sys%hm%pcm%run_pcm) call messages_experimental("PCM for CalculationMode = td")

    if(sys%gr%sb%kpoints%use_symmetries) call messages_experimental("KPoints symmetries with CalculationMode = td")

    call ion_dynamics_init(td%ions, sys%namespace, sys%geo)

    if(ion_dynamics_ions_move(td%ions)) then
      if(sys%gr%sb%kpoints%use_symmetries) then
        message(1) = "KPoints symmetries cannot be used with moving ions."
        message(2) = "Please set KPointsSymmetries = no."
        call messages_fatal(2, namespace=sys%namespace)
      end if
      if(sys%st%symmetrize_density) then
        message(1) = "Symmetrization of the density cannot be used with moving ions."
        message(2) = "Please set SymmetrizeDensity = no."
        call messages_fatal(2, namespace=sys%namespace)
      end if
    end if

    td%iter = 0

    !%Variable TDIonicTimeScale
    !%Type float
    !%Default 1.0
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable defines the factor between the timescale of ionic
    !% and electronic movement. It allows reasonably fast
    !% Born-Oppenheimer molecular-dynamics simulations based on
    !% Ehrenfest dynamics. The value of this variable is equivalent to
    !% the role of <math>\mu</math> in Car-Parrinello. Increasing it
    !% linearly accelerates the time step of the ion
    !% dynamics, but also increases the deviation of the system from the
    !% Born-Oppenheimer surface. The default is 1, which means that both
    !% timescales are the same. Note that a value different than 1
    !% implies that the electrons will not follow physical behaviour.
    !%
    !% According to our tests, values around 10 are reasonable, but it
    !% will depend on your system, mainly on the width of the gap.
    !%
    !% Important: The electronic time step will be the value of
    !% <tt>TDTimeStep</tt> divided by this variable, so if you have determined an
    !% optimal electronic time step (that we can call <i>dte</i>), it is
    !% recommended that you define your time step as:
    !%
    !% <tt>TDTimeStep</tt> = <i>dte</i> * <tt>TDIonicTimeScale</tt>
    !%
    !% so you will always use the optimal electronic time step
    !% (<a href=http://arxiv.org/abs/0710.3321>more details</a>).
    !%End
    call parse_variable(sys%namespace, 'TDIonicTimeScale', CNST(1.0), td%mu)

    if (td%mu <= M_ZERO) then
       write(message(1),'(a)') 'Input: TDIonicTimeScale must be positive.'
       call messages_fatal(1)
    end if

    call messages_print_var_value(stdout, 'TDIonicTimeScale', td%mu)


    !%Variable TDTimeStep
    !%Type float
    !%Section Time-Dependent::Propagation
    !%Description
    !% The time-step for the time propagation. For most propagators you
    !% want to use the largest value that is possible without the
    !% evolution becoming unstable.
    !%
    !% The default value is the maximum value that we have found
    !% empirically that is stable for the spacing <math>h</math>:
    !% <math>dt = 0.0426 - 0.207 h + 0.808 h^2</math>
    !% (from parabolic fit to Fig. 4 of http://dx.doi.org/10.1021/ct800518j,
    !% probably valid for 3D systems only).
    !% However, you might need to adjust this value.
    !%End

    spacing = minval(sys%gr%mesh%spacing(1:sys%gr%sb%dim))
    default_dt = CNST(0.0426) - CNST(0.207)*spacing + CNST(0.808)*spacing**2
    default_dt = default_dt*td%mu

    call parse_variable(sys%namespace, 'TDTimeStep', default_dt, td%dt, unit = units_inp%time)

    if (td%dt <= M_ZERO) then
       write(message(1),'(a)') 'Input: TDTimeStep must be positive.'
       call messages_fatal(1, namespace=sys%namespace)
    end if

    call messages_print_var_value(stdout, 'TDTimeStep', td%dt, unit = units_out%time)


    if(parse_is_defined(sys%namespace, 'TDMaxSteps') .and. parse_is_defined(sys%namespace, 'TDPropagationTime')) then
       call messages_write('You cannot set TDMaxSteps and TDPropagationTime at the same time')
       call messages_fatal(namespace=sys%namespace)
    end if

    !%Variable TDPropagationTime
    !%Type float
    !%Section Time-Dependent::Propagation
    !%Description
    !% The length of the time propagation. You cannot set this variable
    !% at the same time as <tt>TDMaxSteps</tt>. By default this variable will
    !% not be used.
    !%
    !% The units for this variable are <math>\hbar</math>/Hartree (or <math>\hbar</math>/eV if you
    !% selected <tt>ev_angstrom</tt> as input units). The approximate conversions to
    !% femtoseconds are 1 fs = 41.34 <math>\hbar</math>/Hartree = 1.52 <math>\hbar</math>/eV.
    !%End
    call parse_variable(sys%namespace, 'TDPropagationTime', CNST(-1.0), propagation_time, unit = units_inp%time)

    call messages_obsolete_variable(sys%namespace, 'TDMaximumIter', 'TDMaxSteps')

    !%Variable TDMaxSteps
    !%Type integer
    !%Default 1500
    !%Section Time-Dependent::Propagation
    !%Description
    !% Number of time-propagation steps that will be performed. You
    !% cannot use this variable together with <tt>TDPropagationTime</tt>.
    !%End
    default = 1500
    if(propagation_time > CNST(0.0)) default = nint(propagation_time/td%dt)
    call parse_variable(sys%namespace, 'TDMaxSteps', default, td%max_iter)

    if(propagation_time <= CNST(0.0)) propagation_time = td%dt*td%max_iter

    call messages_print_var_value(stdout, 'TDPropagationTime', propagation_time, unit = units_out%time)
    call messages_print_var_value(stdout, 'TDMaxSteps', td%max_iter)

    if(td%max_iter < 1) then
       write(message(1), '(a,i6,a)') "Input: '", td%max_iter, "' is not a valid value for TDMaxSteps."
       message(2) = '(TDMaxSteps <= 1)'
       call messages_fatal(2, namespace=sys%namespace)
    end if

    td%iter = 0

    td%dt = td%dt/td%mu

    ! now the photoelectron stuff
    call pes_init(td%pesv, sys%namespace, sys%gr%mesh, sys%gr%sb, sys%st, sys%outp%restart_write_interval, sys%hm, td%max_iter, &
         td%dt)

    !%Variable TDDynamics
    !%Type integer
    !%Default ehrenfest
    !%Section Time-Dependent::Propagation
    !%Description
    !% Type of dynamics to follow during a time propagation.
    !% For BO, you must set <tt>MoveIons = yes</tt>.
    !%Option ehrenfest 1
    !% Ehrenfest dynamics.
    !%Option bo 2
    !% Born-Oppenheimer (Experimental).
    !%End

    call parse_variable(sys%namespace, 'TDDynamics', EHRENFEST, td%dynamics)
    if(.not.varinfo_valid_option('TDDynamics', td%dynamics)) call messages_input_error(sys%namespace, 'TDDynamics')
    call messages_print_var_option(stdout, 'TDDynamics', td%dynamics)
    if(td%dynamics .ne. EHRENFEST) then
       if(.not.ion_dynamics_ions_move(td%ions)) call messages_input_error(sys%namespace, 'TDDynamics')
    end if

    !%Variable RecalculateGSDuringEvolution
    !%Type logical
    !%Default no
    !%Section Time-Dependent::Propagation
    !%Description
    !% In order to calculate some information about the system along the
    !% evolution (e.g. projection onto the ground-state KS determinant,
    !% projection of the TDKS spin-orbitals onto the ground-state KS
    !% spin-orbitals), the ground-state KS orbitals are needed. If the
    !% ionic potential changes -- that is, the ions move -- one may want
    !% to recalculate the ground state. You may do this by setting this
    !% variable.
    !%
    !% The recalculation is not done every time step, but only every
    !% <tt>RestartWriteInterval</tt> time steps.
    !%End
    call parse_variable(sys%namespace, 'RecalculateGSDuringEvolution', .false., td%recalculate_gs)
    if( sys%hm%lda_u_level /= DFT_U_NONE .and. td%recalculate_gs) &
         call messages_not_implemented("DFT+U with RecalculateGSDuringEvolution=yes", namespace=sys%namespace)

    !%Variable TDScissor
    !%Type float
    !%Default 0.0
    !%Section Time-Dependent
    !%Description
    !% (experimental) If set, a scissor operator will be applied in the
    !% Hamiltonian, shifting the excitation energies by the amount
    !% specified. By default, it is not applied.
    !%End
    call parse_variable(sys%namespace, 'TDScissor', CNST(0.0), td%scissor)
    td%scissor = units_to_atomic(units_inp%energy, td%scissor)
    call messages_print_var_value(stdout, 'TDScissor', td%scissor)

    call propagator_init(sys%gr, sys%namespace, sys%st, td%tr, &
         ion_dynamics_ions_move(td%ions) .or. gauge_field_is_applied(sys%hm%ep%gfield), family_is_mgga_with_exc(sys%ks%xc))

    if(sys%hm%ep%no_lasers>0.and.mpi_grp_is_root(mpi_world)) then
       call messages_print_stress(stdout, "Time-dependent external fields", namespace=sys%namespace)
       call laser_write_info(sys%hm%ep%lasers, stdout, td%dt, td%max_iter)
       call messages_print_stress(stdout, namespace=sys%namespace)
    end if

    !%Variable TDEnergyUpdateIter
    !%Type integer
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable controls after how many iterations Octopus
    !% updates the total energy during a time-propagation run. For
    !% iterations where the energy is not updated, the last calculated
    !% value is reported. If you set this variable to 1, the energy
    !% will be calculated in each step. The default value is 10,
    !% unless the ions move, in which case the default is 1.
    !%End

    default = 10
    if(ion_dynamics_ions_move(td%ions)) default = 1
    call parse_variable(sys%namespace, 'TDEnergyUpdateIter', default, td%energy_update_iter)

    if(ion_dynamics_ions_move(td%ions) .and. td%energy_update_iter /= 1) then
       call messages_experimental('TDEnergyUpdateIter /= 1 when moving ions')
    end if

    if(sys%gr%der%boundaries%spiralBC .and. sys%hm%ep%reltype == SPIN_ORBIT) then
       message(1) = "Generalized Bloch theorem cannot be used with spin-orbit coupling."
       call messages_fatal(1, namespace=sys%namespace)
    end if

    if(sys%gr%der%boundaries%spiralBC .and. &
         any(abs(sys%hm%ep%kick%easy_axis(1:2)) > M_EPSILON)) then
       message(1) = "Generalized Bloch theorem cannot be used for an easy axis along the z direction."
       call messages_fatal(1, namespace=sys%namespace)
    end if

    POP_SUB(td_init)
  end subroutine td_init

  ! ---------------------------------------------------------
  
  subroutine td_end(td)
    type(td_t), intent(inout) :: td

    PUSH_SUB(td_end)

    call pes_end(td%pesv)
    call propagator_end(td%tr)  ! clean the evolution method
    call ion_dynamics_end(td%ions)

    if(td%dynamics == BO) call scf_end(td%scf)

    POP_SUB(td_end)
  end subroutine td_end

  ! ---------------------------------------------------------
  
  subroutine td_run(sys, fromScratch)
    type(system_t), target, intent(inout) :: sys
    logical,                intent(inout) :: fromScratch

    type(td_t)                   :: td
    type(td_write_t)             :: write_handler
    type(grid_t),        pointer :: gr   ! some shortcuts
    type(states_elec_t), pointer :: st
    type(geometry_t),    pointer :: geo
    logical                      :: stopping
#ifdef HAVE_MPI
    logical                      :: stopping_tmp
#endif
    integer                      :: iter, ierr, scsteps
    FLOAT                        :: etime
    type(profile_t),        save :: prof
    type(restart_t)              :: restart_load, restart_dump

    PUSH_SUB(td_run)


    call td_init(td, sys)

    ! some shortcuts
    gr  => sys%gr
    geo => sys%geo
    st  => sys%st

    ! Allocate wavefunctions during time-propagation
    if(td%dynamics == EHRENFEST) then
       !Note: this is not really clean to do this
       if(sys%hm%lda_u_level /= DFT_U_NONE .and. states_are_real(st)) then
          call lda_u_end(sys%hm%lda_u)
          !complex wfs are required for Ehrenfest
          call states_elec_allocate_wfns(st, gr%mesh, TYPE_CMPLX, packed=.true.)
          call lda_u_init(sys%hm%lda_u, sys%namespace, sys%hm%lda_u_level, gr, geo, st, sys%hm%psolver)
       else
          !complex wfs are required for Ehrenfest
          call states_elec_allocate_wfns(st, gr%mesh, TYPE_CMPLX, packed=.true.)
       end if
    else
       call states_elec_allocate_wfns(st, gr%mesh, packed=.true.)
       call scf_init(td%scf, sys%namespace, sys%gr, sys%geo, sys%st, sys%mc, sys%hm, sys%ks)
    end if

    if(sys%hm%scdm_EXX) then
       call scdm_init(sys%hm%scdm, sys%namespace, st, gr%der, sys%hm%psolver%cube, operate_on_scdm = .true.)
       ! make sure scdm is constructed as soon as it is needed
       scdm_is_local = .false.
    end if

    if (gauge_field_is_applied(sys%hm%ep%gfield)) then
       !if the gauge field is applied, we need to tell v_ks to calculate the current
       call v_ks_calculate_current(sys%ks, .true.)

       ! initialize the vector field and update the hamiltonian
       call gauge_field_init_vec_pot(sys%hm%ep%gfield, gr%sb, st)
       call hamiltonian_elec_update(sys%hm, gr%mesh, sys%namespace, time = td%dt*td%iter)
    end if

    call init_wfs()

    if(td%iter >= td%max_iter) then
       call end_()
       POP_SUB(td_run)
       return
    end if

    ! Calculate initial forces and kinetic energy
    if(ion_dynamics_ions_move(td%ions)) then
       if(td%iter > 0) then
          call td_read_coordinates()
          call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, gr, geo, st, time = td%iter*td%dt)
       end if

       call forces_calculate(gr, sys%namespace, geo, sys%hm, st, sys%ks, t = td%iter*td%dt, dt = td%dt)

       geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)
    else
       if(bitand(sys%outp%what, OPTION__OUTPUT__FORCES) /= 0) then
          call forces_calculate(gr, sys%namespace, geo, sys%hm, st, sys%ks, t = td%iter*td%dt, dt = td%dt)
       end if
    end if

    call td_write_init(write_handler, sys%namespace, sys%outp, gr, st, sys%hm, geo, sys%ks, &
         ion_dynamics_ions_move(td%ions), gauge_field_is_applied(sys%hm%ep%gfield), &
         sys%hm%ep%kick, td%iter, td%max_iter, td%dt, sys%mc)

    if(td%scissor > M_EPSILON) then
       call scissor_init(sys%hm%scissor, sys%namespace, st, gr, sys%hm%d, td%scissor, sys%mc)
    end if

    if(td%iter == 0) call td_run_zero_iter()

    if (gauge_field_is_applied(sys%hm%ep%gfield)) call gauge_field_get_force(sys%hm%ep%gfield, gr, st)

    !call td_check_trotter(td, sys, h)
    td%iter = td%iter + 1

    call restart_init(restart_dump, sys%namespace, RESTART_TD, RESTART_TYPE_DUMP, sys%mc, ierr, mesh=gr%mesh)
    if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) then
       ! We will also use the TD restart directory as temporary storage during the time propagation
       call restart_init(restart_load, sys%namespace, RESTART_TD, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh)
    end if

    call messages_print_stress(stdout, "Time-Dependent Simulation", namespace=sys%namespace)
    call print_header(sys%namespace)

    if(td%pesv%calc_spm .or. td%pesv%calc_mask .and. fromScratch) then
       call pes_init_write(td%pesv,gr%mesh,st, sys%namespace)
    end if

    if(st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call st%pack()

    etime = loct_clock()
    ! This is the time-propagation loop. It starts at t=0 and finishes
    ! at td%max_iter*dt. The index i runs from 1 to td%max_iter, and
    ! step "iter" means propagation from (iter-1)*dt to iter*dt.
    propagation: do iter = td%iter, td%max_iter

       stopping = clean_stop(sys%mc%master_comm) .or. walltimer_alarm()

#ifdef HAVE_MPI
       call MPI_Allreduce(stopping, stopping_tmp, 1, MPI_LOGICAL, MPI_LOR, sys%mc%master_comm, mpi_err)
       stopping = stopping_tmp
#endif      

       call profiling_in(prof, "TIME_STEP")

       if(iter > 1) then
          if( ((iter-1)*td%dt <= sys%hm%ep%kick%time) .and. (iter*td%dt > sys%hm%ep%kick%time) ) then
             if( .not.sys%hm%pcm%localf ) then
                call kick_apply(gr%mesh, st, td%ions, geo, sys%hm%ep%kick, sys%hm%psolver)
             else
                call kick_apply(gr%mesh, st, td%ions, geo, sys%hm%ep%kick, sys%hm%psolver, pcm = sys%hm%pcm)
             end if
             call td_write_kick(sys%outp, sys%namespace, gr%mesh, sys%hm%ep%kick, geo, iter)
             !We activate the sprial BC only after the kick,
             !to be sure that the first iteration corresponds to the ground state
             if(gr%der%boundaries%spiralBC) gr%der%boundaries%spiral = .true.
          end if
       end if

       ! in case use scdm localized states for exact exchange and request a new localization
       if(sys%hm%scdm_EXX) scdm_is_local = .false.

       ! time iterate the system, one time step.
       select case(td%dynamics)
       case(EHRENFEST)
          call propagator_dt(sys%ks, sys%namespace, sys%hm, gr, st, td%tr, iter*td%dt, td%dt, &
               td%energy_update_iter*td%mu, iter, td%ions, geo, sys%outp, scsteps = scsteps, &
               update_energy = (mod(iter, td%energy_update_iter) == 0) .or. (iter == td%max_iter) )
       case(BO)
          call propagator_dt_bo(td%scf, sys%namespace, gr, sys%ks, st, sys%hm, geo, sys%mc, sys%outp, iter, td%dt, &
               td%ions, scsteps)
       end select

       !Apply mask absorbing boundaries
       if(sys%hm%bc%abtype == MASK_ABSORBING) call zvmask(gr%mesh, sys%hm, st)

       !Photoelectron stuff
       if (td%pesv%calc_spm .or. td%pesv%calc_mask .or. td%pesv%calc_flux) then
          call pes_calc(td%pesv, sys%namespace, gr%mesh, st, td%dt, iter, gr, sys%hm, stopping)
       end if

       call td_write_iter(write_handler, sys%namespace, sys%outp, gr, st, sys%hm, geo, sys%hm%ep%kick, td%dt, iter)

       ! write down data
       call check_point()

       ! check if debug mode should be enabled or disabled on the fly
       call io_debug_on_the_fly(sys%namespace)

       call profiling_out(prof)
       if (stopping) exit

    end do propagation

    if(st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call st%unpack()

    call restart_end(restart_dump)
    if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) call restart_end(restart_load)
    call td_write_end(write_handler)
    call end_()

    POP_SUB(td_run)

  contains

    subroutine print_header(namespace)
    type(namespace_t),    intent(in)    :: namespace

      write(message(1), '(a7,1x,a14,a14,a10,a17)') 'Iter ', 'Time ', 'Energy ', 'SC Steps', 'Elapsed Time '

      call messages_info(1)
      call messages_print_stress(stdout, namespace%get())

    end subroutine print_header

    ! ---------------------------------------------------------
    subroutine check_point()
      PUSH_SUB(td_run.check_point)

      write(message(1), '(i7,1x,2f14.6,i10,f14.3)') iter, &
        units_from_atomic(units_out%time, iter*td%dt), &
        units_from_atomic(units_out%energy, sys%hm%energy%total + geo%kinetic_energy), &
        scsteps, loct_clock() - etime

      call messages_info(1)
      etime = loct_clock()

      if((sys%outp%output_interval > 0 .and. mod(iter, sys%outp%output_interval) == 0) .or. &
        iter == td%max_iter .or. stopping) then ! output
        ! TODO this now overwrites wf inside st. If this is not wanted need to add an optional overwrite=no flag
        if (st%modelmbparticles%nparticle > 0) then
          call modelmb_sym_all_states (gr, st)
        end if
        call td_write_output(sys%namespace, gr, st, sys%hm, sys%ks, sys%outp, geo, iter, td%dt)
      end if

      if (mod(iter, sys%outp%restart_write_interval) == 0 .or. iter == td%max_iter .or. stopping) then ! restart
        !if(iter == td%max_iter) sys%outp%iter = ii - 1
        call td_write_data(write_handler)
        call td_dump(restart_dump, sys%namespace, gr, st, sys%hm, td, iter, ierr)
        if (ierr /= 0) then
          message(1) = "Unable to write time-dependent restart information."
          call messages_warning(1, namespace=sys%namespace)
        end if

        call pes_output(td%pesv, sys%namespace, gr%mesh, st, iter, sys%outp, td%dt, gr, geo)

        if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) then
          call messages_print_stress(stdout, 'Recalculating the ground state.', namespace=sys%namespace)
          fromScratch = .false.
          call states_elec_deallocate_wfns(sys%st)
          call ground_state_run(sys, fromScratch)
          call states_elec_allocate_wfns(sys%st, gr%mesh, packed=.true.)
          call td_load(restart_load, sys%namespace, gr, st, sys%hm, td, ierr)
          if (ierr /= 0) then
            message(1) = "Unable to load TD states."
            call messages_fatal(1, namespace=sys%namespace)
          end if
          call density_calc(st, gr, st%rho)
          call v_ks_calc(sys%ks, sys%namespace, sys%hm, st, sys%geo, calc_eigenval=.true., time = iter*td%dt, calc_energy=.true.)
          call forces_calculate(gr, sys%namespace, geo, sys%hm, st, sys%ks, t = iter*td%dt, dt = td%dt)
          call messages_print_stress(stdout, "Time-dependent simulation proceeds", namespace=sys%namespace)
          call print_header(sys%namespace)
        end if
      end if

      POP_SUB(td_run.check_point)
    end subroutine check_point

    ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(td_run.end_)

      ! free memory
      call states_elec_deallocate_wfns(st)
      call td_end(td)
      if (ion_dynamics_ions_move(td%ions) .and. td%recalculate_gs) call restart_end(restart_load)

      POP_SUB(td_run.end_)
    end subroutine end_

    ! ---------------------------------------------------------
    subroutine init_wfs()

      integer :: ierr, freeze_orbitals
      FLOAT :: x
      logical :: freeze_hxc, freeze_occ, freeze_u
      type(restart_t) :: restart, restart_frozen

      PUSH_SUB(td_run.init_wfs)

      !%Variable TDFreezeOrbitals
      !%Type integer
      !%Default 0
      !%Section Time-Dependent
      !%Description
      !% (Experimental) You have the possibility of "freezing" a number of orbitals during a time-propagation.
      !% The Hartree and exchange-correlation potential due to these orbitals (which
      !% will be the lowest-energy ones) will be added during the propagation, but the orbitals
      !% will not be propagated.
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
      call parse_variable(sys%namespace, 'TDFreezeOrbitals', 0, freeze_orbitals)

      if(freeze_orbitals /= 0) then
        call messages_experimental('TDFreezeOrbitals')

        if(sys%hm%lda_u_level /= DFT_U_NONE) then
          call messages_not_implemented('TDFreezeOrbitals with DFT+U', namespace=sys%namespace)
        end if
      end if

      if (.not. fromscratch) then
        !We redistribute the states before the restarting
        if(freeze_orbitals > 0) then
          call states_elec_freeze_redistribute_states(st, sys%namespace, sys%gr, sys%mc, freeze_orbitals)
        end if

        call restart_init(restart, sys%namespace, RESTART_TD, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh)
        if(ierr == 0) &
          call td_load(restart, sys%namespace, gr, st, sys%hm, td, ierr)
        if(ierr /= 0) then
          fromScratch = .true.
          td%iter = 0
          message(1) = "Unable to read time-dependent restart information: Starting from scratch"
          call messages_warning(1, namespace=sys%namespace)
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
        call restart_init(restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh, exact=.true.)

        if(.not. st%only_userdef_istates) then
          if(ierr == 0) call states_elec_load(restart, sys%namespace, st, gr, ierr, label = ": gs")
          if (ierr /= 0) then
            message(1) = 'Unable to read ground-state wavefunctions.'
            call messages_fatal(1, namespace=sys%namespace)
          end if
        end if

        ! check if we should deploy user-defined wavefunctions.
        ! according to the settings in the input file the routine
        ! overwrites orbitals that were read from restart/gs
        if(parse_is_defined(sys%namespace, 'UserDefinedStates')) call states_elec_read_user_def_orbitals(gr%mesh, sys%namespace, st)

        call transform_states(st, sys%namespace, restart, gr)
        call restart_end(restart)
      end if

      !We activate the sprial BC only after the kick, 
      !to be sure that the first iteration corresponds to the ground state
      if(gr%der%boundaries%spiralBC) then
        if((td%iter-1)*td%dt > sys%hm%ep%kick%time .and. gr%der%boundaries%spiralBC) then
          gr%der%boundaries%spiral = .true.
        end if
        sys%hm%hm_base%spin => st%spin
        !We fill st%spin. In case of restart, we read it in td_load
        if(fromScratch) call states_elec_fermi(st, sys%namespace, gr%mesh) 
      end if

      ! Initialize the occupation matrices and U for LDA+U
      ! This must be called before parsing TDFreezeOccupations and TDFreezeU
      ! in order that the code does properly the initialization.
      call lda_u_update_occ_matrices(sys%hm%lda_u, sys%namespace, gr%mesh, st, sys%hm%hm_base, sys%hm%energy)

      if(freeze_orbitals > 0) then
        if(fromScratch) then
          ! In this case, we first freeze the orbitals, then calculate the Hxc potential.
          call states_elec_freeze_orbitals(st, sys%namespace, gr, sys%mc, freeze_orbitals, family_is_mgga(sys%ks%xc_family))
        else
          call restart_init(restart, sys%namespace, RESTART_TD, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh)
          if(ierr == 0) &
            call td_load_frozen(restart, gr, st, sys%hm, ierr)
          if(ierr /= 0) then
            td%iter = 0
            message(1) = "Unable to read frozen restart information."
            call messages_fatal(1, namespace=sys%namespace)
          end if
          call restart_end(restart)
        end if
        write(message(1),'(a,i4,a,i4,a)') 'Info: The lowest', freeze_orbitals, &
          ' orbitals have been frozen.', st%nst, ' will be propagated.'
        call messages_info(1)
        call states_elec_freeze_adjust_qtot(st)
        call density_calc(st, gr, st%rho)
        call v_ks_calc(sys%ks, sys%namespace, sys%hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
      elseif(freeze_orbitals < 0) then
        ! This means SAE approximation. We calculate the Hxc first, then freeze all
        ! orbitals minus one.
        write(message(1),'(a)') 'Info: The single-active-electron approximation will be used.'
        call messages_info(1)
        call v_ks_calc(sys%ks, sys%namespace, sys%hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
        if(fromScratch) then
          call states_elec_freeze_orbitals(st, sys%namespace, gr, sys%mc, st%nst-1, family_is_mgga(sys%ks%xc_family))
        else
           call messages_not_implemented("TDFreezeOrbials < 0 with FromScratch=no", namespace=sys%namespace)
        end if
        call v_ks_freeze_hxc(sys%ks)
        call density_calc(st, gr, st%rho)
      else
        ! Normal run.
        call density_calc(st, gr, st%rho)
        call v_ks_calc(sys%ks, sys%namespace, sys%hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
      end if

      !%Variable TDFreezeHXC
      !%Type logical
      !%Default no
      !%Section Time-Dependent
      !%Description
      !% The electrons are evolved as independent particles feeling the Hartree and 
      !% exchange-correlation potentials from the ground-state electronic configuration.
      !%End
      call parse_variable(sys%namespace, 'TDFreezeHXC', .false., freeze_hxc)
      if(freeze_hxc) then 
        write(message(1),'(a)') 'Info: Freezing Hartree and exchange-correlation potentials.'
        call messages_info(1)
        call v_ks_freeze_hxc(sys%ks)

        !In this case we should reload GS wavefunctions 
        if(.not.fromScratch) then
          call messages_not_implemented("TDFreezeHXC with FromScratch=no", namespace=sys%namespace)
        end if
      end if

      x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
      if(st%parallel_in_states) then
        call MPI_Bcast(x, 1, MPI_FLOAT, 0, st%mpi_grp%comm, mpi_err)
      end if
#endif
      call sys%hm%update_span(minval(gr%mesh%spacing(1:gr%mesh%sb%dim)), x)
      ! initialize Fermi energy
      ! Note: we already have calculated the spin above.
      call states_elec_fermi(st, sys%namespace, gr%mesh, compute_spin = .not. gr%der%boundaries%spiralBC)
      call energy_calc_total(sys%namespace, sys%hm, gr, st)

      !%Variable TDFreezeDFTUOccupations
      !%Type logical
      !%Default no
      !%Section Time-Dependent
      !%Description
      !% The occupation matrices than enters in the LDA+U potential
      !% are not evolved during the time evolution.
      !%End
      call parse_variable(sys%namespace, 'TDFreezeDFTUOccupations', .false., freeze_occ)
      if(freeze_occ) then
        write(message(1),'(a)') 'Info: Freezing DFT+U occupation matrices that enters in the DFT+U potential.'
        call messages_info(1)
        call lda_u_freeze_occ(sys%hm%lda_u)

        !In this case we should reload GS wavefunctions 
        if(sys%hm%lda_u_level /= DFT_U_NONE .and..not.fromScratch) then 
          call restart_init(restart_frozen, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh)
          call lda_u_load(restart_frozen, sys%hm%lda_u, st, sys%hm%energy%dft_u, ierr, occ_only = .true.)
          call restart_end(restart_frozen)
        end if
      end if

      !%Variable TDFreezeU
      !%Type logical
      !%Default no
      !%Section Time-Dependent
      !%Description
      !% The effective U of LDA+U is not evolved during the time evolution.
      !%End
      call parse_variable(sys%namespace, 'TDFreezeU', .false., freeze_u)
      if(freeze_u) then
        write(message(1),'(a)') 'Info: Freezing the effective U of DFT+U.'
        call messages_info(1)
        call lda_u_freeze_u(sys%hm%lda_u)

        !In this case we should reload GS wavefunctions
        if(sys%hm%lda_u_level == DFT_U_ACBN0 .and. .not.fromScratch) then
          call restart_init(restart_frozen, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh)
          call lda_u_load(restart_frozen, sys%hm%lda_u, st, sys%hm%energy%dft_u, ierr, u_only = .true.)
          call restart_end(restart_frozen)    
          write(message(1),'(a)') 'Loaded GS effective U of DFT+U'
          call messages_info(1)
          call lda_u_write_U(sys%hm%lda_u, stdout)
          call lda_u_write_V(sys%hm%lda_u, stdout)
        end if
      end if

      POP_SUB(td_run.init_wfs)
    end subroutine init_wfs


    ! ---------------------------------------------------------
    subroutine td_run_zero_iter()
      PUSH_SUB(td_run.td_run_zero_iter)

      call td_write_iter(write_handler, sys%namespace, sys%outp, gr, st, sys%hm, geo, sys%hm%ep%kick, td%dt, 0)

      ! I apply the delta electric field *after* td_write_iter, otherwise the
      ! dipole matrix elements in write_proj are wrong
      if(abs(sys%hm%ep%kick%time)  <=  M_EPSILON) then
        if( .not.sys%hm%pcm%localf ) then
          call kick_apply(gr%mesh, st, td%ions, geo, sys%hm%ep%kick, sys%hm%psolver)
        else
          call kick_apply(gr%mesh, st, td%ions, geo, sys%hm%ep%kick, sys%hm%psolver, pcm = sys%hm%pcm)
        end if
        call td_write_kick(sys%outp, sys%namespace, gr%mesh, sys%hm%ep%kick, geo, 0)

        !We activate the sprial BC only after the kick 
        if(gr%der%boundaries%spiralBC) then
          gr%der%boundaries%spiral = .true.
        end if
      end if
      call propagator_run_zero_iter(sys%hm, gr, td%tr)
      if (sys%outp%output_interval > 0) then
        call td_write_data(write_handler)
        call td_write_output(sys%namespace, gr, st, sys%hm, sys%ks, sys%outp, geo, 0)
      end if

      POP_SUB(td_run.td_run_zero_iter)
    end subroutine td_run_zero_iter


    ! ---------------------------------------------------------
    !> reads the pos and vel from coordinates file
    subroutine td_read_coordinates() 
      integer :: iatom, iter, iunit
      PUSH_SUB(td_run.td_read_coordinates)

      call io_assign(iunit)
      iunit = io_open(io_workpath('td.general/coordinates', sys%namespace), sys%namespace, action='read', status='old')

      if(iunit < 0) then
        message(1) = "Could not open file '"//trim(io_workpath('td.general/coordinates', sys%namespace))//"'."
        message(2) = "Starting simulation from initial geometry."
        call messages_warning(2, namespace=sys%namespace)
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
        geo%atom(iatom)%x(:) = units_to_atomic(units_out%length, geo%atom(iatom)%x(:))
      end do
      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%v(1:gr%mesh%sb%dim)
        geo%atom(iatom)%v(:) = units_to_atomic(units_out%velocity, geo%atom(iatom)%v(:))
      end do
      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%f(1:gr%mesh%sb%dim)
        geo%atom(iatom)%f(:) = units_to_atomic(units_out%force, geo%atom(iatom)%f(:))
      end do

      call io_close(iunit)

      POP_SUB(td_run.td_read_coordinates)
    end subroutine td_read_coordinates

  end subroutine td_run


  ! ---------------------------------------------------------
  subroutine multisys_td_run(systems, fromScratch)
    type(multisystem_t), intent(inout) :: systems
    logical,             intent(inout) :: fromScratch

    integer :: it, internal_loop
    integer, parameter :: MAX_PROPAGATOR_STEPS = 1000
    FLOAT :: smallest_algo_dt

    PUSH_SUB(multisys_td_run)

    call messages_write('Info: Running Multi-System time evolution')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    ! this should eventually be moved up to run.F90 when all systems
    ! are derived classes from system_abst
    call systems%init_interactions()

    ! Initialize all propagators and find the smallest time-step
    smallest_algo_dt = CNST(1e10)
    call systems%init_propagator(smallest_algo_dt)

    ! Initialize all the clocks
    call systems%init_clocks(smallest_algo_dt)

    ! Set initial conditions
    call systems%initial_conditions(.true.)

    call messages_print_stress(stdout, "Multi-system propagation", namespace=systems%namespace)

    call systems%propagation_start()

    ! The full TD loop
    it = 0
    do while (.not. systems%has_reached_final_propagation_time())

      it = it + 1

      internal_loop = 1
      do while (.not. systems%propagation_step_is_done() .and. internal_loop < MAX_PROPAGATOR_STEPS)
        call systems%dt_operation()
        internal_loop = internal_loop + 1
      end do
      call systems%propagation_step_finish(it)

      write (message(1), '(a)') repeat ('-', 71)
      call messages_info(1)
    end do

    call systems%propagation_finish()

    POP_SUB(multisys_td_run)
  end subroutine multisys_td_run


  ! ---------------------------------------------------------
  subroutine transform_states(st, namespace, restart, gr, prefix)
    type(states_elec_t),        intent(inout) :: st
    type(namespace_t),          intent(in)    :: namespace
    type(restart_t),            intent(inout) :: restart
    type(grid_t),               intent(in)    :: gr
    character(len=*), optional, intent(in)    :: prefix

    type(states_elec_t) :: stin
    type(block_t) :: blk
    CMPLX, allocatable :: rotation_matrix(:,:), psi(:, :)
    integer :: ist, jst, ncols, iqn
    character(len=256) :: block_name
    
    PUSH_SUB(transform_states)

    block_name = trim(optional_default(prefix, "")) // "TransformStates"
    
    !%Variable TransformStates
    !%Type block
    !%Default no
    !%Section States
    !%Description
    !% Before starting the <tt>td</tt> calculation, the initial states (that are
    !% read from the <tt>restart/gs</tt> directory, which should have been
    !% generated in a previous ground-state calculation) can be "transformed"
    !% among themselves. The block <tt>TransformStates</tt> gives the transformation matrix
    !% to be used. The number of rows and columns of the matrix should equal the number
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
    !% Note: This variable cannot be used when parallel in states.
    !%End
    if(parse_is_defined(namespace, trim(block_name))) then
      if(parse_block(namespace, trim(block_name), blk) == 0) then
        if(st%parallel_in_states) &
          call messages_not_implemented(trim(block_name) // " parallel in states", namespace=namespace)
        if(parse_block_n(blk) /= st%nst) then
          message(1) = "Number of rows in block " // trim(block_name) // " must equal number of states in this calculation."
          call messages_fatal(1, namespace=namespace)
        end if
        call states_elec_copy(stin, st, exclude_wfns = .true.)
        call states_elec_look_and_load(restart, namespace, stin, gr)

        ! FIXME: rotation matrix should be R_TYPE
        SAFE_ALLOCATE(rotation_matrix(1:stin%nst, 1:stin%nst))
        SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))

        rotation_matrix = M_z0
        do ist = 1, stin%nst
          rotation_matrix(ist, ist) = CNST(1.0)
        end do

        do ist = 1, st%nst
          ncols = parse_block_cols(blk, ist-1)
          if(ncols /= stin%nst) then
            write(message(1),'(a,i6,a,i6,3a,i6,a)') "Number of columns (", ncols, ") in row ", ist, " of block ", &
              trim(block_name), " must equal number of states (", stin%nst, ") read from gs restart."
            call messages_fatal(1, namespace=namespace)
          end if
          do jst = 1, stin%nst
            call parse_block_cmplx(blk, ist - 1, jst - 1, rotation_matrix(jst, ist))
          end do
        end do

        call parse_block_end(blk)

        do iqn = st%d%kpt%start, st%d%kpt%end
          if(states_are_real(st)) then
            call states_elec_rotate(stin, namespace, gr%mesh, TOFLOAT(rotation_matrix), iqn)
          else
            call states_elec_rotate(stin, namespace, gr%mesh, rotation_matrix, iqn)
          end if

          do ist = st%st_start, st%st_end 
            call states_elec_get_state(stin, gr%mesh, ist, iqn, psi)
            call states_elec_set_state(st, gr%mesh, ist, iqn, psi)
          end do

        end do

        SAFE_DEALLOCATE_A(rotation_matrix)
        SAFE_DEALLOCATE_A(psi)

        call states_elec_end(stin)

      else
        call messages_input_error(namespace, trim(block_name), '"' // trim(block_name) // '" has to be specified as block.')
      end if
      
    end if

    POP_SUB(transform_states)
  end subroutine transform_states

  ! ---------------------------------------------------------
  subroutine td_dump(restart, namespace, gr, st, hm, td, iter, ierr)
    type(restart_t),          intent(in)  :: restart
    type(namespace_t),        intent(in)  :: namespace
    type(grid_t),             intent(in)  :: gr
    type(states_elec_t),      intent(in)  :: st
    type(hamiltonian_elec_t), intent(in)  :: hm
    type(td_t),               intent(in)  :: td
    integer,                  intent(in)  :: iter
    integer,                  intent(out) :: ierr

    integer :: err, err2

    PUSH_SUB(td_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(td_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing td restart."
      call messages_info(1)
    end if

    ! first write resume file
    call states_elec_dump(restart, st, gr, err, iter=iter)
    if (err /= 0) ierr = ierr + 1

    call states_elec_dump_rho(restart, st, gr, ierr, iter=iter)
    if (err /= 0) ierr = ierr + 1 

    if(hm%lda_u_level /= DFT_U_NONE) then
      call lda_u_dump(restart, hm%lda_u, st, ierr)
      if (err /= 0) ierr = ierr + 1
    end if

    call potential_interpolation_dump(td%tr%vksold, restart, gr, st%d%nspin, err2)
    if (err2 /= 0) ierr = ierr + 2

    call pes_dump(td%pesv, namespace, restart, st, gr%mesh, err)
    if (err /= 0) ierr = ierr + 4

    ! Gauge field restart
    if (gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_dump(restart, hm%ep%gfield, ierr)
    end if

    if(gr%der%boundaries%spiralBC) then
      call states_elec_dump_spin(restart, st, err)
      if(err /= 0) ierr = ierr + 8
    end if

    if(associated(st%frozen_rho)) then
      call states_elec_dump_frozen(restart, st, gr, ierr)
    end if

    if (debug%info) then
      message(1) = "Debug: Writing td restart done."
      call messages_info(1)
    end if

    POP_SUB(td_dump)
  end subroutine td_dump

  ! ---------------------------------------------------------
  subroutine td_load(restart, namespace, gr, st, hm, td, ierr)
    type(restart_t),     intent(in)    :: restart
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(in)    :: gr
    type(states_elec_t), intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(td_t),          intent(inout) :: td
    integer,             intent(out)   :: ierr

    integer :: err, err2
    PUSH_SUB(td_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(td_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading td restart."
      call messages_info(1)
    end if

    ! Read states
    call states_elec_load(restart, namespace, st, gr, err, iter=td%iter, label = ": td")
    if (err /= 0) then
      ierr = ierr + 1
    end if

    ! read potential from previous interactions
    call potential_interpolation_load(td%tr%vksold, namespace, restart, gr, st%d%nspin, err2)

    if (err2 /= 0) ierr = ierr + 2

    ! read PES restart
    if (td%pesv%calc_spm .or. td%pesv%calc_mask .or. td%pesv%calc_flux) then
      call pes_load(td%pesv, namespace, restart, st, err)
      if (err /= 0) ierr = ierr + 4
    end if

    ! Gauge field restart
    if (gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_load(restart, hm%ep%gfield, err)
      if (err /= 0) then
        ierr = ierr + 8
      else
        call hamiltonian_elec_update(hm, gr%mesh, namespace, time = td%dt*td%iter)
      end if
    end if

    if(gr%der%boundaries%spiralBC) then
      call states_elec_load_spin(restart, st, err)
      !To ensure back compatibility, if the file is not present, we use the 
      !current states to get the spins
      if(err /= 0) call states_elec_fermi(st, namespace, gr%mesh)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading td restart done."
      call messages_info(1)
    end if

    POP_SUB(td_load)
  end subroutine td_load

  ! ---------------------------------------------------------
  subroutine td_load_frozen(restart, gr, st, hm, ierr)
    type(restart_t),     intent(in)    :: restart
    type(grid_t),        intent(in)    :: gr
    type(states_elec_t), intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,             intent(out)   :: ierr

    PUSH_SUB(td_load_frozen)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(td_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading td frozen restart."
      call messages_info(1)
    end if

    SAFE_ALLOCATE(st%frozen_rho(1:gr%mesh%np,1:st%d%nspin))
    if(family_is_mgga(hm%xc%family)) then
      SAFE_ALLOCATE(st%frozen_tau(1:gr%mesh%np,1:st%d%nspin))
      SAFE_ALLOCATE(st%frozen_gdens(1:gr%mesh%np,1:gr%sb%dim,1:st%d%nspin))
      SAFE_ALLOCATE(st%frozen_ldens(1:gr%mesh%np,1:st%d%nspin))
    end if

    call states_elec_load_frozen(restart, st, gr, ierr)

    if (debug%info) then
      message(1) = "Debug: Reading td frozen restart done."
      call messages_info(1)
    end if

    POP_SUB(td_load_frozen)
  end subroutine td_load_frozen


end module td_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
