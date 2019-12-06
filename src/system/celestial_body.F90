!! Copyright (C) 2019 N. Tancogne-Dejean
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

module celestial_body_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use system_abst_oct_m

  implicit none

  private
  public ::               &
    celestial_body_t

  type, extends(system_abst_t) :: celestial_body_t
    private

    FLOAT :: mass
    FLOAT :: pos(1:MAX_DIM)
    FLOAT :: vel(1:MAX_DIM)
    FLOAT :: acc(1:MAX_DIM)

    class(propagator_abst_t), pointer :: prop

  contains
    procedure :: do_td_operation => celestial_body_do_td
    procedure :: pull_interaction => celestial_body_pull
    procedure :: get_needed_quantity => celestial_body_needed_quantity
    procedure(celestial_body_init) :: init
    procedure :: set_propagator => celestial_body_set_prop
  end type celestial_body_t

contains

  subroutine init(sys)
    class(celestial_body_t),  intent(in) :: sys

    PUSH_SUB(celestial_body_init)


    POP_SUB(celestial_body_init)

  end subroutine init

  ! ---------------------------------------------------------
  integer function celestial_body_needed_quantity(sys)
    class(celestial_body_t),  intent(in) :: sys

    PUSH_SUB(celestial_body_needed_quantity)

    ce_needed_quantity = FORCE

    POP_SUB(celestial_body_needed_quantity)

  end function celestial_body_needed_quantity

  ! ---------------------------------------------------------
  subroutine celestial_body_set_prop(sys, prop)
    class(celestial_body_t),          intent(inout) :: sys
    class(propagator_abst_t), target, intent(in)    :: prop

    PUSH_SUB(celestial_body_set_prop)

    sys%prop => prop

    POP_SUB(celestial_body_set_prop)

  end subroutine celestial_body_set_prop

end module celestial_body_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
