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

module interaction_gravity_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use messages_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use system_abst_oct_m

  implicit none

  private
  public ::                &
    interaction_gravity_t

  type, extends(interaction_abst_t) :: interaction_gravity_t
    integer :: dim

    FLOAT :: force(MAX_DIM)

    class(system_abst_t), pointer :: partner

    FLOAT, pointer :: system_mass
    FLOAT, pointer :: system_pos(:)

    FLOAT, pointer :: partner_mass
    FLOAT, pointer :: partner_pos(:)

  contains
    procedure :: update => interaction_gravity_update
    procedure :: end => interaction_gravity_end
    final :: interaction_gravity_finalize
  end type interaction_gravity_t

  interface interaction_gravity_t
    module procedure interaction_gravity_init
  end interface interaction_gravity_t

contains

  ! ---------------------------------------------------------

  function interaction_gravity_init(dim, partner) result(this)
    integer,                      intent(in) :: dim
    class(system_abst_t), target, intent(in) :: partner
    class(interaction_gravity_t), pointer    :: this

    PUSH_SUB(interaction_gravity_init)

    SAFE_ALLOCATE(this)

    this%dim = dim
    this%partner => partner

    !Gravity interaction needs only one quantity from each system, which is the position
    this%n_system_quantities = 2
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%system_quantities(this%n_system_quantities))
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = MASS
    this%partner_quantities(1) = POSITION
    this%partner_quantities(2) = MASS

    call partner%set_pointers_to_interaction(this)

    POP_SUB(interaction_gravity_init)
  end function interaction_gravity_init

  ! ---------------------------------------------------------
  logical function interaction_gravity_update(this, clock) result(updated)
    class(interaction_gravity_t), intent(inout) :: this
    class(clock_t),               intent(in)    :: clock

    logical :: allowed_to_update
    FLOAT, parameter :: GG = CNST(6.67430e-11)
    FLOAT :: dist3

    PUSH_SUB(interaction_gravity_update)

    !The interaction has already been updated to the desired time
    if (this%clock == clock) then
      if (debug%info) then
        write(message(1), '(a,a)') "Debug: -- Interaction already up-to-date with ", trim(this%partner%namespace%get())
        call messages_info(1)
      end if
      updated = .true.
      POP_SUB(interaction_gravity_update)
      return
    end if

    allowed_to_update = this%partner%update_exposed_quantities(clock, this%n_partner_quantities, this%partner_quantities)

    if (allowed_to_update) then
      !We can now compute the interaction from the updated pointers
      ASSERT(associated(this%partner_pos))
      ASSERT(associated(this%system_pos))
      ASSERT(associated(this%partner_mass))
      ASSERT(associated(this%system_mass))

      dist3 = sum((this%partner_pos(1:this%dim) - this%system_pos(1:this%dim))**2)**(M_THREE/M_TWO)

      this%force(1:this%dim) = (this%partner_pos(1:this%dim) - this%system_pos(1:this%dim)) &
        / dist3 * (GG * this%system_mass * this%partner_mass)

      ! Update was successful, so set new interaction time
      updated = .true.
      call this%clock%set_time(clock)

      if (debug%info) then
        write(message(1), '(a,a)') "Debug: -- Updated interaction with ", trim(this%partner%namespace%get())
        write(message(2), '(a,i3,a,i3)') "Debug: ---- clocks are ", clock%get_tick(), " and ", this%partner%clock%get_tick()
        call messages_info(2)
      end if
    else
      if (debug%info) then
        write(message(1), '(a,a)') "Debug: -- Cannot update yet the interaction with ", trim(this%partner%namespace%get())
        write(message(2), '(a,i3,a,i3)') "Debug: ---- clocks are ", clock%get_tick(), " and ", this%partner%clock%get_tick()
        call messages_info(2)
      end if
      updated = .false.
    end if

    POP_SUB(interaction_gravity_update)
  end function interaction_gravity_update

  ! ---------------------------------------------------------
  subroutine interaction_gravity_end(this)
    class(interaction_gravity_t), intent(inout) :: this

    PUSH_SUB(interaction_gravity_end)

    nullify(this%partner)
    this%force = M_ZERO
    nullify(this%system_mass)
    nullify(this%system_pos)
    nullify(this%partner_mass)
    nullify(this%partner_pos)

    call interaction_abst_end(this)

    POP_SUB(interaction_gravity_end)
  end subroutine interaction_gravity_end

  ! ---------------------------------------------------------
  subroutine interaction_gravity_finalize(this)
    type(interaction_gravity_t), intent(inout) :: this

    PUSH_SUB(interaction_gravity_finalize)

    call this%end()

    POP_SUB(interaction_gravity_finalize)
  end subroutine interaction_gravity_finalize

end module interaction_gravity_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
