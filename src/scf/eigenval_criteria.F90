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
  use quantity_oct_m
  use unit_oct_m
  use unit_system_oct_m 

  implicit none

  private
  public ::                  &
    eigenval_criteria_t

  type, extends(convergence_criteria_t) :: eigenval_criteria_t
    private

  contains
    procedure :: write_info   => criteria_write_info
  end type eigenval_criteria_t

  interface eigenval_criteria_t
    procedure eigenval_criteria_constructor
  end interface eigenval_criteria_t


contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function eigenval_criteria_constructor(tol, absolute) result(crit)
    FLOAT,                       intent(in) :: tol
    logical,                     intent(in) :: absolute
    class(eigenval_criteria_t),  pointer :: crit

    PUSH_SUB(eigenval_criteria_constructor)

    SAFE_ALLOCATE(crit)

    crit%tol = tol
    crit%absolute = absolute
    crit%quantity = SUMEIGENVAL

    POP_SUB(eigenval_criteria_constructor)
  end function eigenval_criteria_constructor

  ! ---------------------------------------------------------
  subroutine criteria_write_info(this, iunit)
    class(eigenval_criteria_t),  intent(inout) :: this
    integer,                        intent(in)    :: iunit

    PUSH_SUB(criteria_write_info)

    if(this%absolute) then
      write(iunit, '(6x, a, es15.8,a,es15.8,4a)') 'abs_evsum = ', this%val, &
            ' (', units_from_atomic(units_out%energy, this%tol), ')', &
            ' [',  trim(units_abbrev(units_out%energy)), ']'
    else
      write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'rel_evsum = ', this%val, ' (', this%tol, ')'
    end if
     
    POP_SUB(criteria_write_info)
  end subroutine criteria_write_info

end module eigenval_criteria_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
