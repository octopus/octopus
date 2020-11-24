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
  use quantity_oct_m

  implicit none

  private
  public ::                  &
    density_criteria_t

  type, extends(convergence_criteria_t) :: density_criteria_t
    private

  contains
    procedure :: write_info   => criteria_write_info
    final     :: density_criteria_end
  end type density_criteria_t

  interface density_criteria_t
    procedure density_criteria_constructor
  end interface density_criteria_t


contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function density_criteria_constructor(tol_abs, tol_rel) result(crit)
    FLOAT,                       intent(in) :: tol_abs
    FLOAT,                       intent(in) :: tol_rel
    class(density_criteria_t),      pointer :: crit

    PUSH_SUB(density_criteria_constructor)

    SAFE_ALLOCATE(crit)

    crit%tol_abs = tol_abs
    crit%tol_rel = tol_rel
    crit%quantity = DENSITY

    POP_SUB(density_criteria_constructor)
  end function density_criteria_constructor

  ! ---------------------------------------------------------
  subroutine criteria_write_info(this, iunit)
    class(density_criteria_t),  intent(inout) :: this
    integer,                    intent(in)    :: iunit

    PUSH_SUB(criteria_write_info)

    write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'abs_dens = ', this%val_abs, ' (', this%tol_abs, ')'
    write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'rel_dens = ', this%val_rel, ' (', this%tol_rel, ')'
     
    POP_SUB(criteria_write_info)
  end subroutine criteria_write_info

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
