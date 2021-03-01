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

module eigenval_criteria_oct_m
  use convergence_criteria_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use unit_oct_m

  implicit none

  private
  public ::                  &
    eigenval_criteria_t

  type, extends(convergence_criteria_t) :: eigenval_criteria_t
    private

  contains
    final     :: eigenval_criteria_end
  end type eigenval_criteria_t

  interface eigenval_criteria_t
    procedure eigenval_criteria_constructor
  end interface eigenval_criteria_t


contains

  ! ---------------------------------------------------------
  function eigenval_criteria_constructor(tol_abs, tol_rel, unit) result(crit)
    FLOAT,                       intent(in) :: tol_abs
    FLOAT,                       intent(in) :: tol_rel
    type(unit_t),        target, intent(in) :: unit
    class(eigenval_criteria_t),  pointer :: crit

    PUSH_SUB(eigenval_criteria_constructor)

    SAFE_ALLOCATE(crit)

    crit%tol_abs = tol_abs
    crit%tol_rel = tol_rel
    crit%unit => unit
    crit%label = 'evsum'

    POP_SUB(eigenval_criteria_constructor)
  end function eigenval_criteria_constructor

  ! ---------------------------------------------------------
  subroutine eigenval_criteria_end(this)
    type(eigenval_criteria_t),   intent(inout) :: this

    PUSH_SUB(eigenval_criteria_end)

    call convergence_criteria_end(this)

    POP_SUB(eigenval_criteria_end)
  end subroutine eigenval_criteria_end

end module eigenval_criteria_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
