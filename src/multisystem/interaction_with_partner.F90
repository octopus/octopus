!! Copyright (C) 2020 M. Oliveira, Heiko Appel
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

module interaction_with_partner_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use profiling_oct_m
  implicit none

  private
  public ::                      &
    interaction_with_partner_t,  &
    interaction_with_partner_end

  !> Some interactions involve two systems. In this case the interaction is a
  !! unidirectional relationship between those two systems. One of the systems
  !! owns the interaction and feels its effects. The other system is referred to
  !! as the interaction partner.
  type, extends(interaction_t), abstract :: interaction_with_partner_t
    private
    class(interaction_partner_t), public, pointer :: partner

    integer,              public :: n_partner_quantities  !< Number of quantities needed from the partner
    integer, allocatable, public :: partner_quantities(:) !< Identifiers of the quantities needed from the partner
  contains
    procedure :: update => interaction_with_partner_update
  end type interaction_with_partner_t

contains

  ! ---------------------------------------------------------
  logical function interaction_with_partner_update(this, requested_time) result(updated)
    class(interaction_with_partner_t), intent(inout) :: this
    class(clock_t),                    intent(in)    :: requested_time

    logical :: allowed_to_update
    type(event_handle_t) :: debug_handle

    PUSH_SUB(interaction_with_partner_update)

    ! We should only try to update the interaction if it is not yet at the requested time
    ASSERT(.not. (this%clock == requested_time))

    debug_handle = multisystem_debug_write_event_in(event = event_function_call_t("interaction_with_partner_update"), &
                                                   extra="target: "//trim(this%label)//"-"//trim(this%partner%namespace%get()), &
                                                   interaction_clock = this%clock,      &
                                                   partner_clock = this%partner%clock,  & 
                                                   requested_clock = requested_time)

    allowed_to_update = this%partner%update_exposed_quantities(requested_time, this)

    if (allowed_to_update) then
      ! We can now compute the interaction from the updated quantities
      call this%calculate()

      ! Update was successful, so set new interaction time
      updated = .true.
      call this%clock%set_time(requested_time)
      call multisystem_debug_write_marker(event = event_clock_update_t( "interaction", &
                                                    trim(this%label)//"-"//trim(this%partner%namespace%get()), & 
                                                    this%clock, "set") )

      if (debug%info) then
        write(message(1), '(a,a,a,a,a)') "Debug: -- Updated '", trim(this%label), "' interaction with '", &
          trim(this%partner%namespace%get()), "'"
        write(message(2), '(a,f16.6,a,f16.6)') "Debug: ---- Requested time is ", requested_time%time(), &
          " and partner time is ", this%partner%clock%time()
        call messages_info(2)
      end if
    else
      if (debug%info) then
        write(message(1), '(a,a,a,a,a)') "Debug: -- Cannot update yet the '", trim(this%label), "' interaction with '", &
          trim(this%partner%namespace%get()), "'"
        write(message(2), '(a,f16.6,a,f16.6)') "Debug: ---- Requested time is ", requested_time%time(), &
          " and partner time is ", this%partner%clock%time()
        call messages_info(2)
      end if
      updated = .false.
    end if

    call multisystem_debug_write_event_out(debug_handle, update = updated,   &
                                            interaction_clock = this%clock,  &
                                            partner_clock = this%partner%clock, &
                                            requested_clock = requested_time)

    POP_SUB(interaction_with_partner_update)
  end function interaction_with_partner_update

  ! ---------------------------------------------------------
  subroutine interaction_with_partner_end(this)
    class(interaction_with_partner_t), intent(inout) :: this

    PUSH_SUB(interaction_with_partner_end)

    SAFE_DEALLOCATE_A(this%partner_quantities)
    nullify(this%partner)
    call interaction_end(this)

    POP_SUB(interaction_with_partner_end)
  end subroutine interaction_with_partner_end

end module interaction_with_partner_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
