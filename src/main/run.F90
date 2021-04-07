!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2020 M. Oliveira
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
  use external_potential_oct_m
  use fft_oct_m
  use geom_opt_oct_m
  use global_oct_m
  use ground_state_oct_m
  use interactions_factory_oct_m
  use interaction_partner_oct_m
  use invert_ks_oct_m
  use messages_oct_m
  use mpi_debug_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multisystem_basic_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use opt_control_oct_m
  use parser_oct_m
  use phonons_fd_oct_m
  use phonons_lr_oct_m
  use poisson_oct_m
  use kdotp_oct_m
  use profiling_oct_m
  use pulpo_oct_m
  use restart_oct_m
  use static_pol_oct_m
  use system_factory_oct_m
  use system_oct_m
  use td_oct_m
  use test_oct_m
  use time_dependent_oct_m
  use unit_system_oct_m
  use unocc_oct_m
  use varinfo_oct_m
  use vdw_oct_m

  implicit none

  private
  public ::                      &
    run

  integer, parameter :: LR = 1, FD = 2

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
      call messages_input_error(namespace, 'ResponseMethod')
    end if

    POP_SUB(get_resp_method)
  end function get_resp_method
  
  ! ---------------------------------------------------------
  subroutine run(namespace, calc_mode_id)
    type(namespace_t), intent(in) :: namespace
    integer,           intent(in) :: calc_mode_id

    type(partner_list_t) :: partners
    class(system_t), pointer :: systems
    type(system_factory_t) :: system_factory
    type(interactions_factory_t) :: interactions_factory
    type(profile_t), save :: calc_mode_prof
    logical :: from_scratch
    integer :: iunit_out
    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner

    PUSH_SUB(run)

    call messages_print_stress(stdout, "Calculation Mode")
    call messages_print_var_option(stdout, "CalculationMode", calc_mode_id)
    call messages_print_stress(stdout)

    call calc_mode_init()

    if (calc_mode_id == OPTION__CALCULATIONMODE__RECIPE) then
      call pulpo_print()
      POP_SUB(run)
      return
    end if

    call restart_module_init(namespace)

    call unit_system_init(namespace)

    call accel_init(mpi_world, namespace)

    ! initialize FFTs
    call fft_all_init(namespace)

    if (calc_mode_id == OPTION__CALCULATIONMODE__TEST) then
      call test_run(namespace)
      call fft_all_end()
#ifdef HAVE_MPI
      call mpi_debug_statistics()
#endif
      POP_SUB(run)
      return
    end if

    ! Create systems
    if (parse_is_defined(namespace, "Systems")) then
      ! We are running in multi-system mode
      systems => multisystem_basic_t(namespace, system_factory)
    else
      ! Fall back to old behaviour
      systems => electrons_t(namespace, generate_epot = calc_mode_id /= OPTION__CALCULATIONMODE__DUMMY)
    end if

    ! initialize everything that needs parallelization
    call systems%init_parallelization(mpi_world)

    ! Create list of partners
    select type (systems)
    class is (multisystem_basic_t)
      ! Systems are also partners
      partners = systems%list

      ! Add external potentials to partners list
      call load_external_potentials(partners, namespace)

    type is (electrons_t)
      call partners%add(systems)
    end select

    ! Create and initialize interactions
    call interactions_factory%create_interactions(systems, partners)
    call systems%init_all_interactions()

    select type (systems)
    class is (multisystem_basic_t)
      ! Write the interaction graph as a DOT graph for debug
      if ( (debug%interaction_graph .or. debug%interaction_graph_full) .and. mpi_grp_is_root(mpi_world)) then
        iunit_out = io_open('debug/interaction_graph.dot', systems%namespace, action='write')
        write(iunit_out, '(a)') 'digraph {'
        call systems%write_interaction_graph(iunit_out, debug%interaction_graph_full)
        write(iunit_out, '(a)') '}'
        call io_close(iunit_out)
      end if
    end select

    if (.not. systems%process_is_slave()) then
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
      call parse_variable(namespace, 'FromScratch', .false., from_scratch)

      call profiling_in(calc_mode_prof, "CALC_MODE")

      select case (calc_mode_id)
      case (OPTION__CALCULATIONMODE__GS)
        call ground_state_run(systems, from_scratch)
      case (OPTION__CALCULATIONMODE__UNOCC)
        call unocc_run(systems, from_scratch)
      case (OPTION__CALCULATIONMODE__TD)
        call time_dependent_run(systems, from_scratch)
      case (OPTION__CALCULATIONMODE__GO)
        call geom_opt_run(systems, from_scratch)
      case (OPTION__CALCULATIONMODE__OPT_CONTROL)
        call opt_control_run(systems)
      case (OPTION__CALCULATIONMODE__EM_RESP)
        select case(get_resp_method(namespace))
        case(FD)
          call static_pol_run(systems, from_scratch)
        case(LR)
          call em_resp_run(systems, from_scratch)
        end select
      case (OPTION__CALCULATIONMODE__CASIDA)
        call casida_run(systems, from_scratch)
      case (OPTION__CALCULATIONMODE__VDW)
        call vdW_run(systems, from_scratch)
      case (OPTION__CALCULATIONMODE__VIB_MODES)
        select case(get_resp_method(namespace))
        case(FD)
          call phonons_run(systems)
        case(LR)
          call phonons_lr_run(systems, from_scratch)
        end select
      case (OPTION__CALCULATIONMODE__ONE_SHOT)
        message(1) = "CalculationMode = one_shot is obsolete. Please use gs with MaximumIter = 0."
        call messages_fatal(1)
      case (OPTION__CALCULATIONMODE__KDOTP)
        call kdotp_lr_run(systems, from_scratch)
      case (OPTION__CALCULATIONMODE__DUMMY)
      case (OPTION__CALCULATIONMODE__INVERT_KS)
        call invert_ks_run(systems)
      case (OPTION__CALCULATIONMODE__RECIPE)
        ASSERT(.false.) !this is handled before, if we get here, it is an error
      end select

      call profiling_out(calc_mode_prof)
    end if

    select type (systems)
    class is (multisystem_basic_t)
      !Deallocate the external potentials
      call iter%start(partners)
      do while (iter%has_next())
        select type(ptr => iter%get_next())
        class is(external_potential_t)
          partner => ptr
          SAFE_DEALLOCATE_P(partner)
        end select
      end do
    end select

    ! Finalize systems
    SAFE_DEALLOCATE_P(systems)

    call fft_all_end()

    call accel_end()

#ifdef HAVE_MPI
    call mpi_debug_statistics()
#endif

    POP_SUB(run)

  contains

    subroutine calc_mode_init()

      PUSH_SUB(calc_mode_init)

      select case (calc_mode_id)
      case (OPTION__CALCULATIONMODE__GS, OPTION__CALCULATIONMODE__GO, OPTION__CALCULATIONMODE__UNOCC)
        call ground_state_run_init()
      case (OPTION__CALCULATIONMODE__TD)
        call td_run_init()
      case (OPTION__CALCULATIONMODE__CASIDA)
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
