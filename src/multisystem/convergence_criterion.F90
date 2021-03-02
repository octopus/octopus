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

module convergence_criterion_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m
  use unit_oct_m

  implicit none

  private
  public ::                    &
    convergence_criterion_t,    &
    convergence_criterion_end,  &
    criterion_list_t,           &
    criterion_iterator_t


  type, abstract :: convergence_criterion_t
    private
    FLOAT,   public :: tol_abs !< Tolerance of the convergence criterion
    FLOAT,   public :: tol_rel !< Tolerance of the convergence criterion
    character(len=:), allocatable, public :: label
    type(unit_t), pointer, public :: unit => null()

    FLOAT, public   :: val_abs    !< Current value of the criterion
    FLOAT, public   :: val_rel  
    FLOAT, pointer  :: value_diff
    FLOAT, pointer  :: norm

  contains
    procedure :: is_converged => convergence_criterion_is_converged
    procedure :: set_pointers => criterion_set_quantity_pointers
    procedure :: write_info => criterion_write_info
  end type convergence_criterion_t

  !> These classes extend the list and list iterator to make a criterion list
  type, extends(linked_list_t) :: criterion_list_t
    private
  contains
    procedure :: add => criterion_list_add_node
  end type criterion_list_t

  type, extends(linked_list_iterator_t) :: criterion_iterator_t
    private
  contains
    procedure :: get_next => criterion_iterator_get_next
  end type criterion_iterator_t


contains

  ! ---------------------------------------------------------
  subroutine criterion_write_info(this, iunit)
    class(convergence_criterion_t), intent(inout) :: this
    integer,                       intent(in)    :: iunit

    PUSH_SUB(criterion_write_info)

    if (associated(this%unit)) then
      write(iunit, '(6x,a,a,a,es15.8,a,es15.8,4a)') 'abs_', this%label, ' = ', &
        units_from_atomic(this%unit, this%val_abs), &
        ' (', units_from_atomic(this%unit, this%tol_abs), ')', &
        ' [',  trim(units_abbrev(this%unit)), ']'
    else
      write(iunit, '(6x,a,a,a,es15.8,a,es15.8,a)') 'abs_', this%label, ' = ', this%val_abs, ' (', this%tol_abs, ')'
    end if
    write(iunit, '(6x,a,a,a,es15.8,a,es15.8,a)') 'rel_', this%label, ' = ', this%val_rel, ' (', this%tol_rel, ')'
     
    POP_SUB(criterion_write_info)
  end subroutine criterion_write_info

  ! ---------------------------------------------------------
  !> Is the convergence reached ?
  subroutine convergence_criterion_is_converged(this, is_converged)
    class(convergence_criterion_t),  intent(inout) :: this
    logical,                        intent(out)   :: is_converged

    PUSH_SUB(convergence_criterion_is_converged)
   
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


    POP_SUB(convergence_criterion_is_converged)
  end subroutine convergence_criterion_is_converged

  ! ---------------------------------------------------------
  !> Setting pointers to the in, out and norm values
  subroutine criterion_set_quantity_pointers(this, value_diff, value_norm)
    class(convergence_criterion_t),  intent(inout) :: this
    FLOAT, target,                  intent(in)    :: value_diff
    FLOAT, target,                  intent(in)    :: value_norm
  
    PUSH_SUB(criterion_set_quantity_pointers)
  
    this%value_diff => value_diff
    this%norm => value_norm
  
    POP_SUB(criterion_set_quantity_pointers)
  end subroutine criterion_set_quantity_pointers

  ! ---------------------------------------------------------
  subroutine convergence_criterion_end(this)
    class(convergence_criterion_t),  intent(inout) :: this

    PUSH_SUB(convergence_criterion_end)

    nullify(this%value_diff)
    nullify(this%norm)
    nullify(this%unit)
    if (allocated(this%label)) then
      ! No safe_dealloate here, because it was not allocated with safe_allocate.
      deallocate(this%label)
    end if

    POP_SUB(convergence_criterion_end)
  end subroutine convergence_criterion_end

  ! ---------------------------------------------------------
  subroutine criterion_list_add_node(this, criterion)
    class(criterion_list_t)                 :: this
    class(convergence_criterion_t), target :: criterion

    PUSH_SUB(criterion_list_add_node)

    call this%add_ptr(criterion)

    POP_SUB(criterion_list_add_node)
  end subroutine criterion_list_add_node

  ! ---------------------------------------------------------
  function criterion_iterator_get_next(this) result(criterion)
    class(criterion_iterator_t),      intent(inout) :: this
    class(convergence_criterion_t),  pointer       :: criterion

    PUSH_SUB(criterion_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (convergence_criterion_t)
      criterion => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(criterion_iterator_get_next)
  end function criterion_iterator_get_next

end module convergence_criterion_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
