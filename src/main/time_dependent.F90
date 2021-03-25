!! Copyright (C) 2019-2020 M. Oliveira, H. Appel, N. Tancogne-Dejean
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

module time_dependent_oct_m
  use electrons_oct_m
  use global_oct_m
  use messages_oct_m
  use multisystem_basic_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use td_oct_m
  use unit_system_oct_m

  implicit none

  private
  public :: time_dependent_run

contains

  ! ---------------------------------------------------------
  subroutine time_dependent_run(system, from_scratch)
    class(*), intent(inout) :: system
    logical,  intent(inout) :: from_scratch

    PUSH_SUB(time_dependent_run)

    select type (system)
    class is (multisystem_basic_t)
      call time_dependent_run_multisystem(system, from_scratch)
    type is (electrons_t)
      call time_dependent_run_legacy(system, from_scratch)
    end select

    POP_SUB(time_dependent_run)
  end subroutine time_dependent_run

  ! ---------------------------------------------------------
  subroutine time_dependent_run_multisystem(systems, from_scratch)
    class(multisystem_basic_t), intent(inout) :: systems
    logical,              intent(in)    :: from_scratch

    integer :: internal_loop
    integer, parameter :: MAX_PROPAGATOR_STEPS = 1000
    FLOAT :: final_time

    PUSH_SUB(time_dependent_run_multisystem)

    call messages_write('Info: Running Multi-System time evolution')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    ! Get final propagation time from input
    ! This variable is also defined (and properly documented) in td/td.F90.
    ! This is temporary, until all the propagators are moved to the new framework.
    call parse_variable(systems%namespace, 'TDPropagationTime', CNST(-1.0), final_time, unit = units_inp%time)
    if (final_time <= M_ZERO) then
      call messages_input_error(systems%namespace, 'TDPropagationTime', 'must be greater than zero')
    end if
    call messages_print_var_value(stdout, 'TDPropagationTime', final_time)

    ! Initialize all propagators
    call systems%init_propagator()

    ! Set initial conditions
    call systems%initial_conditions(from_scratch)

    call messages_print_stress(stdout, "Multi-system propagation", namespace=systems%namespace)


    call systems%propagation_start()

    call multisystem_debug_start_log()

    ! The full TD loop
    do while (.not. systems%has_reached_final_propagation_time(final_time))
      internal_loop = 1
      do while (internal_loop < MAX_PROPAGATOR_STEPS)
        ! Do one algorithmic operation
        call systems%dt_operation()

        ! Exit loop if a propagation step is done
        if (systems%prop%step_is_done()) exit

        internal_loop = internal_loop + 1
      end do

      write (message(1), '(a)') repeat ('-', 71)
      call messages_info(1)

    end do

    call systems%propagation_finish()

    POP_SUB(time_dependent_run_multisystem)
  end subroutine time_dependent_run_multisystem

  ! ---------------------------------------------------------
  subroutine time_dependent_run_legacy(electrons, from_scratch)
    class(electrons_t), intent(inout) :: electrons
    logical,            intent(inout) :: from_scratch

    PUSH_SUB(time_dependent_run_legacy)

    call td_init(electrons%td, electrons%namespace, electrons%space, electrons%gr, electrons%geo, electrons%st, electrons%ks, &
      electrons%hm, electrons%outp)
    call td_init_run(electrons%td, electrons%namespace, electrons%mc, electrons%gr, electrons%geo, electrons%st, electrons%ks, &
      electrons%hm, electrons%outp, electrons%space, from_scratch)
    call td_run(electrons%td, electrons%namespace, electrons%mc, electrons%gr, electrons%geo, electrons%st, electrons%ks, &
      electrons%hm, electrons%outp, electrons%space, from_scratch)
    call td_end_run(electrons%td, electrons%st, electrons%hm)
    call td_end(electrons%td)

    POP_SUB(time_dependent_run_legacy)
  end subroutine time_dependent_run_legacy

end module time_dependent_oct_m
