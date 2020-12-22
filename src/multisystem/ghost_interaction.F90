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

module ghost_interaction_oct_m
  use global_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                &
    ghost_interaction_t

  type, extends(interaction_with_partner_t) :: ghost_interaction_t
  contains
    procedure :: calculate => ghost_interaction_calculate
    final :: ghost_interaction_finalize
  end type ghost_interaction_t

  interface ghost_interaction_t
    module procedure ghost_interaction_init
  end interface ghost_interaction_t

contains

  ! ---------------------------------------------------------
  function ghost_interaction_init(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(ghost_interaction_t),           pointer       :: this

    PUSH_SUB(ghost_interaction_init)

    SAFE_ALLOCATE(this)

    this%label = "ghost"

    this%partner => partner

    ! A ghost interaction does not require any quantity
    this%n_system_quantities = 0
    this%n_partner_quantities = 0

    POP_SUB(ghost_interaction_init)
  end function ghost_interaction_init

  ! ---------------------------------------------------------
  subroutine ghost_interaction_calculate(this)
    class(ghost_interaction_t), intent(inout) :: this

    PUSH_SUB(ghost_interaction_calculate)

    ! A ghost interaction does not do anything

    POP_SUB(ghost_interaction_calculate)
  end subroutine ghost_interaction_calculate

  ! ---------------------------------------------------------
  subroutine ghost_interaction_finalize(this)
    type(ghost_interaction_t), intent(inout) :: this

    PUSH_SUB(ghost_interaction_finalize)

    call interaction_with_partner_end(this)

    POP_SUB(ghost_interaction_finalize)
  end subroutine ghost_interaction_finalize

end module ghost_interaction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
