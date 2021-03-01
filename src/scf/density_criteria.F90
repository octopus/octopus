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

module density_criteria_oct_m
  use convergence_criteria_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                  &
    density_criteria_t

  type, extends(convergence_criteria_t) :: density_criteria_t
    private

  contains
    final     :: density_criteria_end
  end type density_criteria_t

  interface density_criteria_t
    procedure density_criteria_constructor
  end interface density_criteria_t


contains

  ! ---------------------------------------------------------
  function density_criteria_constructor(tol_abs, tol_rel) result(crit)
    FLOAT,                       intent(in) :: tol_abs
    FLOAT,                       intent(in) :: tol_rel
    class(density_criteria_t),      pointer :: crit

    PUSH_SUB(density_criteria_constructor)

    SAFE_ALLOCATE(crit)

    crit%tol_abs = tol_abs
    crit%tol_rel = tol_rel
    crit%label = 'dens'

    POP_SUB(density_criteria_constructor)
  end function density_criteria_constructor

  ! ---------------------------------------------------------
  subroutine density_criteria_end(this)
    type(density_criteria_t),   intent(inout) :: this

    PUSH_SUB(density_criteria_end)

    call convergence_criteria_end(this)

    POP_SUB(density_criteria_end)
  end subroutine density_criteria_end


end module density_criteria_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
