!! Copyright (C)  2020 N. Tancogne-Dejean
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

module energy_criteria_oct_m
  use convergence_criteria_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use unit_oct_m

  implicit none

  private
  public ::                  &
    energy_criteria_t

  type, extends(convergence_criteria_t) :: energy_criteria_t
    private

  contains
    final     :: energy_criteria_end 
  end type energy_criteria_t

  interface energy_criteria_t
    procedure energy_criteria_constructor
  end interface energy_criteria_t


contains

  ! ---------------------------------------------------------
  function energy_criteria_constructor(tol_abs, tol_rel, unit) result(crit)
    FLOAT,                       intent(in) :: tol_abs
    FLOAT,                       intent(in) :: tol_rel
    type(unit_t),        target, intent(in) :: unit
    class(energy_criteria_t),  pointer :: crit

    PUSH_SUB(energy_criteria_constructor)

    SAFE_ALLOCATE(crit)

    crit%tol_abs = tol_abs
    crit%tol_rel = tol_rel
    crit%unit => unit
    crit%label = 'en'

    POP_SUB(energy_criteria_constructor)
  end function energy_criteria_constructor

  ! ---------------------------------------------------------
  subroutine energy_criteria_end(this) 
    type(energy_criteria_t),   intent(inout) :: this

    PUSH_SUB(energy_criteria_end)

    call convergence_criteria_end(this)

    POP_SUB(energy_criteria_end)
  end subroutine energy_criteria_end

end module energy_criteria_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
