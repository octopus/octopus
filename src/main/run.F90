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

module run_oct_m
  use accel_oct_m
  use casida_oct_m
  use em_resp_oct_m
  use fft_oct_m
  use geom_opt_oct_m
  use global_oct_m
  use ground_state_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_mxll_oct_m
  use invert_ks_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use mpi_debug_oct_m
  use memory_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use opt_control_oct_m
  use parser_oct_m
  use pcm_oct_m
  use phonons_fd_oct_m
  use phonons_lr_oct_m
  use poisson_oct_m
  use kdotp_oct_m
  use profiling_oct_m
  use pulpo_oct_m
  use restart_oct_m
  use static_pol_oct_m
  use system_oct_m
  use system_mxll_oct_m
  use td_oct_m
  use test_oct_m
  use unit_system_oct_m
  use unocc_oct_m
  use varinfo_oct_m
  use vdw_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                      &
    run

  integer :: calc_mode_id

  integer, parameter :: LR = 1, FD = 2

  integer, public, parameter ::   &
    CM_NONE               =   0,  &
    CM_GS                 =   1,  &
    CM_UNOCC              =   2,  &
    CM_TD                 =   3,  &
    CM_GEOM_OPT           =   5,  &
    CM_OPT_CONTROL        =   7,  &
    CM_LR_POL             =   8,  &
    CM_CASIDA             =   9,  &
    CM_VDW                =  11,  &
    CM_PHONONS_LR         =  12,  &
    CM_ONE_SHOT           =  14,  &
    CM_KDOTP              =  15,  &
    CM_DUMMY              =  17,  &
    CM_INVERTKDS          =  18,  &
    CM_TEST               =  19,  &
    CM_MAXWELL_FREE       =  20,  &
    CM_PULPO_A_FEIRA      =  99

contains

  ! ---------------------------------------------------------
  integer function get_resp_method(namespace)
    type(namespace_t),    intent(in)    :: namespace

    PUSH_SUB(get_resp_method)
    
    !%Variable ResponseMethod
    !%Type integer
    !%Default sternheimer
    !%Section Linear Response
    !%Description
    !% Some response properties can be calculated either via
    !% Sternheimer linear response or by using finite
    !% differences. You can use this variable to select how you want
    !% them to be calculated, it applies to <tt>em_resp</tt> and <tt>vib_modes</tt>
    !% calculation modes. By default, the Sternheimer linear-response
    !% technique is used.
    !%Option sternheimer 1
    !% The linear response is obtained by solving a self-consistent
    !% Sternheimer equation for the variation of the orbitals. This
    !% is the recommended method.
    !%Option finite_differences 2
    !% Properties are calculated as a finite-differences derivative of
    !% the energy obtained by several ground-state calculations. This
    !% method, slow and limited only to static response, is kept
    !% mainly because it is simple and useful for testing purposes.
    !%End
    
    call parse_variable(namespace, 'ResponseMethod', LR, get_resp_method)

    if(.not.varinfo_valid_option('ResponseMethod', get_resp_method)) then
      call messages_input_error('ResponseMethod')
    end if

    POP_SUB(get_resp_method)
  end function get_resp_method
  
  ! ---------------------------------------------------------
  subroutine run(namespace, cm)
    type(namespace_t), intent(in) :: namespace
    integer,           intent(in) :: cm

    type(linked_list_t) :: systems
    type(list_iterator_t) :: iter
    class(*), pointer :: sys
    type(profile_t), save :: calc_mode_prof
    logical :: fromScratch

    PUSH_SUB(run)

    calc_mode_id = cm

    call messages_print_stress(stdout, "Calculation Mode")
    call messages_print_var_option(stdout, "CalculationMode", calc_mode_id)
    call messages_print_stress(stdout)

    call calc_mode_init()

    if(calc_mode_id == CM_PULPO_A_FEIRA) then
      call pulpo_print()
      POP_SUB(run)
      return
    end if

    call restart_module_init(namespace)

    call accel_init(mpi_world, namespace)

    ! initialize FFTs
    call fft_all_init(namespace)

    call unit_system_init(namespace)

    if(calc_mode_id == CM_TEST) then
      call test_run(namespace)
      call fft_all_end()
#ifdef HAVE_MPI
      call mpi_debug_statistics()
#endif
      POP_SUB(run)
      return
    end if

    ! Initialize systems
    call multisystem_init(systems, namespace)

    ! Loop over systems
    call iter%start(systems)
    do while (iter%has_next())
      sys => iter%get_next()
      select type (sys)

      type is (system_mxll_t)
        select case(calc_mode_id)
        case (CM_TD)
           call td_run(sys, fromScratch)
        case default
           message(1) = "Maxwell systems currently support only TD propagation"
           call messages_info(1)
        end select

      type is (system_t)
        if (sys%hm%pcm%run_pcm) then
          select case (calc_mode_id)
          case (CM_GS)
            if (sys%hm%pcm%epsilon_infty /= sys%hm%pcm%epsilon_0 .and. sys%hm%pcm%tdlevel /= PCM_TD_EQ) then
              message(1) = 'Non-equilbrium PCM is not active in a time-independent run.'
              message(2) = 'You set epsilon_infty /= epsilon_0, but epsilon_infty is not relevant for CalculationMode = gs.'
              message(3) = 'By definition, the ground state is in equilibrium with the solvent.'
              message(4) = 'Therefore, the only relevant dielectric constant is the static one.'
              message(5) = 'Nevertheless, the dynamical PCM response matrix is evaluated for benchamarking purposes.'
              call messages_warning(5)
            end if
          case (CM_TD)
            call messages_experimental("PCM for CalculationMode = td")
          case default
            call messages_not_implemented("PCM for CalculationMode /= gs or td")
          end select

          if ( (sys%mc%par_strategy /= P_STRATEGY_SERIAL).and.(sys%mc%par_strategy /= P_STRATEGY_STATES) ) then
            call messages_experimental('Parallel in domain calculations with PCM')
          end if
        end if

        call messages_print_stress(stdout, 'Approximate memory requirements')
        call memory_run(sys)
        call messages_print_stress(stdout)

        if(calc_mode_id /= CM_DUMMY) then
          message(1) = "Info: Generating external potential"
          call messages_info(1)
          call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%gr, sys%geo, sys%st)
          message(1) = "      done."
          call messages_info(1)
        end if

        if(sys%ks%theory_level /= INDEPENDENT_PARTICLES) then
          call poisson_async_init(sys%hm%psolver, sys%mc)
          ! slave nodes do not call the calculation routine
          if(multicomm_is_slave(sys%mc))then
            !for the moment we only have one type of slave
            call poisson_slave_work(sys%hm%psolver)
          end if
        end if

        if(.not. multicomm_is_slave(sys%mc)) then
          call messages_write('Info: Octopus initialization completed.', new_line = .true.)
          call messages_write('Info: Starting calculation mode.')
          call messages_info()

          !%Variable FromScratch
          !%Type logical
          !%Default false
          !%Section Execution
          !%Description
          !% When this variable is set to true, <tt>Octopus</tt> will perform a
          !% calculation from the beginning, without looking for restart
          !% information.
          !%End

          call parse_variable(namespace, 'FromScratch', .false., fromScratch)

          call profiling_in(calc_mode_prof, "CALC_MODE")

          select case(calc_mode_id)
          case(CM_GS)
            call ground_state_run(sys, fromScratch)
          case(CM_UNOCC)
            call unocc_run(sys, fromScratch)
          case(CM_TD)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = td")
            call td_run(sys, fromScratch)
          case(CM_LR_POL)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = em_resp")
            select case(get_resp_method(sys%namespace))
            case(FD)
              call static_pol_run(sys, fromScratch)
            case(LR)
              call em_resp_run(sys, fromScratch)
            end select
          case(CM_VDW)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = vdw")
            call vdW_run(sys, fromScratch)
          case(CM_GEOM_OPT)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = go")
            call geom_opt_run(sys, fromScratch)
          case(CM_PHONONS_LR)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = vib_modes")
            select case(get_resp_method(sys%namespace))
            case(FD)
              call phonons_run(sys)
            case(LR)
              call phonons_lr_run(sys, fromscratch)
            end select
          case(CM_OPT_CONTROL)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = opt_control")
            call opt_control_run(sys)
          case(CM_CASIDA)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = casida")
            call casida_run(sys, fromScratch)
          case(CM_ONE_SHOT)
            message(1) = "CalculationMode = one_shot is obsolete. Please use gs with MaximumIter = 0."
            call messages_fatal(1)
          case(CM_KDOTP)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = kdotp")
            call kdotp_lr_run(sys, fromScratch)
          case(CM_DUMMY)
          case(CM_INVERTKDS)
            if(sys%gr%sb%kpoints%use_symmetries) &
              call messages_experimental("KPoints symmetries with CalculationMode = invert_ks")
            call invert_ks_run(sys)
          case(CM_PULPO_A_FEIRA)
            ASSERT(.false.) !this is handled before, if we get here, it is an error
          end select

          call profiling_out(calc_mode_prof)
        end if

        if(sys%ks%theory_level /= INDEPENDENT_PARTICLES) call poisson_async_end(sys%hm%psolver, sys%mc)

      class default
        message(1) = "Unknow system type."
        call messages_fatal(1)
      end select
    end do

    ! Finalize systems
    call multisystem_end(systems)

    call fft_all_end()

    call accel_end()


#ifdef HAVE_MPI
    call mpi_debug_statistics()
#endif

    POP_SUB(run)

  contains

    subroutine calc_mode_init()

      PUSH_SUB(calc_mode_init)

      select case(calc_mode_id)
      case(CM_GS, CM_GEOM_OPT, CM_UNOCC)
        call ground_state_run_init()
      case(CM_TD)
        call td_run_init()
      case(CM_CASIDA)
        call casida_run_init()
      end select

      POP_SUB(calc_mode_init)
    end subroutine calc_mode_init

  end subroutine run

end module run_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
