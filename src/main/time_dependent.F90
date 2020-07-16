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
  use multisystem_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use td_oct_m

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
    class is (multisystem_t)
      call time_dependent_run_multisystem(system, from_scratch)
    type is (electrons_t)
      call time_dependent_run_legacy(system, from_scratch)
    end select

    POP_SUB(time_dependent_run)
  end subroutine time_dependent_run

  ! ---------------------------------------------------------
  subroutine time_dependent_run_multisystem(systems, from_scratch)
    class(multisystem_t), intent(inout) :: systems
    logical,              intent(in)    :: from_scratch

    integer :: it, internal_loop
    integer, parameter :: MAX_PROPAGATOR_STEPS = 1000
    FLOAT :: smallest_algo_dt

    PUSH_SUB(time_dependent_run_multisystem)

    call messages_write('Info: Running Multi-System time evolution')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    ! Initialize all propagators and find the smallest time-step
    smallest_algo_dt = CNST(1e10)
    call systems%init_propagator(smallest_algo_dt)

    ! Initialize all the clocks
    call systems%init_clocks(smallest_algo_dt)

    ! Set initial conditions
    call systems%initial_conditions(from_scratch)

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

    POP_SUB(time_dependent_run_multisystem)
  end subroutine time_dependent_run_multisystem

  ! ---------------------------------------------------------
  subroutine time_dependent_run_legacy(electrons, from_scratch)
    class(electrons_t), intent(inout) :: electrons
    logical,            intent(inout) :: from_scratch

    PUSH_SUB(time_dependent_run_legacy)

    call td_run(electrons%namespace, electrons%mc, electrons%gr, electrons%geo, electrons%st, electrons%ks, electrons%hm, &
      electrons%outp, from_scratch)

    POP_SUB(time_dependent_run_legacy)
  end subroutine time_dependent_run_legacy

end module time_dependent_oct_m
