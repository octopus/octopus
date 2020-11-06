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

module convergence_criteria_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use unit_oct_m
  use unit_system_oct_m 

  implicit none

  private
  public ::                  &
    convergence_criteria_t,  &
    criteria_list_t,         &
    criteria_iterator_t

  type :: convergence_criteria_t
    private

    integer, public :: type    !< Type of the quantity
    FLOAT,   public :: tol     !< Tolerance of the convergence criteria
    logical, public :: absolute
    
    integer, public :: quantity  !< The quantity that needs to be converged
  
    FLOAT, public   :: val    !< Current value of the criteria
    FLOAT, pointer  :: value_diff
    FLOAT, pointer  :: norm

  contains
    procedure :: is_converged => convergence_criteria_is_converged
    procedure :: set_pointers => criteria_set_quantity_pointers
    procedure :: write_info   => criteria_write_info
  end type convergence_criteria_t

  !> These classes extend the list and list iterator to make a criteria list
  type, extends(linked_list_t) :: criteria_list_t
    private
  contains
    procedure :: add => criteria_list_add_node
  end type criteria_list_t

  type, extends(linked_list_iterator_t) :: criteria_iterator_t
    private
  contains
    procedure :: get_next => criteria_iterator_get_next
  end type criteria_iterator_t

  interface convergence_criteria_t
    procedure convergence_criteria_constructor
  end interface convergence_criteria_t


contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function convergence_criteria_constructor(tol, quantity, absolute) result(crit)
    FLOAT,                       intent(in) :: tol
    integer,                     intent(in) :: quantity
    logical,                     intent(in) :: absolute
    class(convergence_criteria_t),  pointer :: crit

    PUSH_SUB(convergence_criteria_constructor)

    SAFE_ALLOCATE(crit)

    crit%tol = tol
    crit%absolute = absolute
    crit%quantity = quantity

    POP_SUB(convergence_criteria_constructor)
  end function convergence_criteria_constructor

  ! ---------------------------------------------------------
  subroutine criteria_list_add_node(this, criteria)
    class(criteria_list_t)                :: this
    class(convergence_criteria_t), target :: criteria

    PUSH_SUB(criteria_list_add_node)

    call this%add_ptr(criteria)

    POP_SUB(criteria_list_add_node)
  end subroutine criteria_list_add_node

  ! ---------------------------------------------------------
  function criteria_iterator_get_next(this) result(criteria)
    class(criteria_iterator_t),     intent(inout) :: this
    class(convergence_criteria_t),  pointer       :: criteria

    PUSH_SUB(criteria_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (convergence_criteria_t)
      criteria => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(criteria_iterator_get_next)
  end function criteria_iterator_get_next

  ! ---------------------------------------------------------
  !> Is the convergence reached ?
  logical function convergence_criteria_is_converged(this)
    class(convergence_criteria_t),  intent(inout) :: this

    PUSH_SUB(convergence_criteria_is_converged)

    this%val = abs(this%value_diff)
    if(.not. this%absolute) then
      ASSERT(associated(this%norm))
      if(abs(this%norm) <= M_EPSILON) then
        this%val = M_HUGE
      else
        this%val = this%val / abs(this%norm)
      end if
    end if

    convergence_criteria_is_converged = this%val < this%tol

    POP_SUB(convergence_criteria_is_converged)
  end function convergence_criteria_is_converged

  ! ---------------------------------------------------------
  !> Setting pointers to the in, out and norm values
  subroutine criteria_set_quantity_pointers(this, value_diff, value_norm)
    class(convergence_criteria_t),  intent(inout) :: this
    FLOAT, target,                  intent(in)    :: value_diff
    FLOAT, target,                  intent(in)    :: value_norm
  
    PUSH_SUB(criteria_set_quantity_pointers)
  
    this%value_diff => value_diff
    this%norm => value_norm
  
    POP_SUB(criteria_set_quantity_pointers)
  end subroutine criteria_set_quantity_pointers

  ! ---------------------------------------------------------
  subroutine criteria_write_info(this, iunit)
    class(convergence_criteria_t),  intent(inout) :: this
    integer,                        intent(in)    :: iunit

    PUSH_SUB(criteria_write_info)

    select case(this%quantity)
    case(ENERGY)
      if(this%absolute) then
         write(iunit, '(6x, a, es15.8,a,es15.8,4a)') 'abs_en = ', this%val, &
              ' (', units_from_atomic(units_out%energy, this%tol), ')', &
              ' [',  trim(units_abbrev(units_out%energy)), ']'
      else
        write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'rel_en = ', this%val, ' (', this%tol, ')'
      end if
    case(DENSITY)
      if(this%absolute) then
        write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'abs_dens = ', this%val, ' (', this%tol, ')'
      else
        write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'rel_dens = ', this%val, ' (', this%tol, ')'
      end if
    case(EIGENVAL)
      if(this%absolute) then
        write(iunit, '(6x, a, es15.8,a,es15.8,4a)') 'abs_evsum = ', this%val, &
              ' (', units_from_atomic(units_out%energy, this%tol), ')', &
              ' [',  trim(units_abbrev(units_out%energy)), ']'
      else
        write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'rel_evsum = ', this%val, ' (', this%tol, ')'
      end if
    case default
      ASSERT(.false.)
    end select 
     
    POP_SUB(criteria_write_info)
  end subroutine criteria_write_info

end module convergence_criteria_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
