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

  implicit none

  private
  public ::                  &
    convergence_criteria_t,  &
    criteria_list_t,         &
    criteria_iterator_t,     &
    convergence_criteria_end

  type, abstract :: convergence_criteria_t
    private

    integer, public :: type    !< Type of the quantity
    FLOAT,   public :: tol_abs !< Tolerance of the convergence criteria
    FLOAT,   public :: tol_rel !< Tolerance of the convergence criteria
    
    integer, public :: quantity  !< The quantity that needs to be converged
  
    FLOAT, public   :: val_abs    !< Current value of the criteria
    FLOAT, public   :: val_rel  
    FLOAT, pointer  :: value_diff
    FLOAT, pointer  :: norm

  contains
    procedure :: is_converged => convergence_criteria_is_converged
    procedure :: set_pointers => criteria_set_quantity_pointers
    procedure(criteria_write_info), deferred :: write_info
  end type convergence_criteria_t

  abstract interface
    ! ---------------------------------------------------------
    subroutine criteria_write_info(this, iunit)
      import convergence_criteria_t
      class(convergence_criteria_t),  intent(inout) :: this
      integer,                        intent(in)    :: iunit
    end subroutine criteria_write_info
  end interface

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


contains
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
  subroutine convergence_criteria_is_converged(this, is_converged)
    class(convergence_criteria_t),  intent(inout) :: this
    logical,                        intent(out)   :: is_converged

    PUSH_SUB(convergence_criteria_is_converged)
   
    this%val_abs = abs(this%value_diff)

    ASSERT(associated(this%norm))
    if(abs(this%norm) <= M_EPSILON) then
      this%val_rel = M_HUGE
    else
      this%val_rel = this%val_abs / abs(this%norm)
    end if

    if(this%tol_abs <= M_ZERO) then
      is_converged = .true.
    else
      is_converged = this%val_abs < this%tol_abs
    end if
    
    if(this%tol_rel > M_ZERO) then
      is_converged = is_converged .and. this%val_rel < this%tol_rel
    end if


    POP_SUB(convergence_criteria_is_converged)
  end subroutine convergence_criteria_is_converged

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
  subroutine convergence_criteria_end(this)
    class(convergence_criteria_t),  intent(inout) :: this

    PUSH_SUB(convergence_criteria_end)

    nullify(this%value_diff)
    nullify(this%norm)

    POP_SUB(convergence_criteria_end)

  end subroutine convergence_criteria_end

end module convergence_criteria_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
