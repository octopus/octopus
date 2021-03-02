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

module energy_criterion_oct_m
  use convergence_criterion_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use unit_oct_m

  implicit none

  private
  public ::                  &
    energy_criterion_t

  type, extends(convergence_criterion_t) :: energy_criterion_t
    private
  contains
    final :: energy_criterion_finalize
  end type energy_criterion_t

  interface energy_criterion_t
    procedure energy_criterion_constructor
  end interface energy_criterion_t


contains

  ! ---------------------------------------------------------
  function energy_criterion_constructor(tol_abs, tol_rel, unit) result(crit)
    FLOAT,                       intent(in) :: tol_abs
    FLOAT,                       intent(in) :: tol_rel
    type(unit_t),        target, intent(in) :: unit
    class(energy_criterion_t),  pointer :: crit

    PUSH_SUB(energy_criterion_constructor)

    SAFE_ALLOCATE(crit)

    crit%tol_abs = tol_abs
    crit%tol_rel = tol_rel
    crit%unit => unit
    crit%label = 'energy'

    POP_SUB(energy_criterion_constructor)
  end function energy_criterion_constructor

  ! ---------------------------------------------------------
  subroutine energy_criterion_finalize(this)
    type(energy_criterion_t),   intent(inout) :: this

    PUSH_SUB(energy_criterion_finalize)

    call convergence_criterion_end(this)

    POP_SUB(energy_criterion_finalize)
  end subroutine energy_criterion_finalize

end module energy_criterion_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
