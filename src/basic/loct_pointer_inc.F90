!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

! SAFE ALLOCATE in this file causes the PGI compiler to crash:
! Lowering Error: bad ast optype in expression [ast=1125,asttype=12,datatype=0]
! PGF90-F-0000-Internal compiler error. Errors in Lowering       1 (../../../src/basic/loct_pointer_inc.F90: 33)
! ---------------------------------------------------------
subroutine SUBNAME(loct_pointer_copy_1)(pout, pin)
  TYPE, pointer, intent(out) :: pout(:)
  TYPE, pointer, intent(in)  :: pin(:)

  integer :: nl1, nu1
  integer :: i1

  PUSH_SUB(SUBNAME(loct_pointer_copy_1))

  if(associated(pin)) then
    nl1 = lbound(pin, 1)
    nu1 = ubound(pin, 1)
    allocate(pout(nl1:nu1))
    do i1 = nl1, nu1
      pout(i1) = pin(i1)
    end do
  else
    nullify(pout)
  end if

  POP_SUB(SUBNAME(loct_pointer_copy_1))
end subroutine SUBNAME(loct_pointer_copy_1)

! ---------------------------------------------------------
subroutine SUBNAME(loct_pointer_copy_2)(pout, pin)
  TYPE, pointer, intent(out) :: pout(:, :)
  TYPE, pointer, intent(in)  :: pin(:, :)

  integer :: nl1, nu1, nl2, nu2
  integer :: i1, i2

  PUSH_SUB(SUBNAME(loct_pointer_copy_2))

  if(associated(pin)) then
    nl1 = lbound(pin, 1)
    nu1 = ubound(pin, 1)
    nl2 = lbound(pin, 2)
    nu2 = ubound(pin, 2)
    allocate(pout(nl1:nu1, nl2:nu2))
    do i2 = nl2, nu2
      do i1 = nl1, nu1
        pout(i1, i2) = pin(i1, i2)
      end do
    end do
  else
    nullify(pout)
  end if

  POP_SUB(SUBNAME(loct_pointer_copy_2))
end subroutine SUBNAME(loct_pointer_copy_2)


! ---------------------------------------------------------
subroutine SUBNAME(loct_pointer_copy_3)(pout, pin)
  TYPE, pointer, intent(out) :: pout(:, :, :)
  TYPE, pointer, intent(in)  :: pin(:, :, :)

  integer :: nl1, nu1, nl2, nu2, nl3, nu3
  integer :: i1, i2, i3

  PUSH_SUB(SUBNAME(loct_pointer_copy_3))

  if(associated(pin)) then
    nl1 = lbound(pin, 1)
    nu1 = ubound(pin, 1)
    nl2 = lbound(pin, 2)
    nu2 = ubound(pin, 2)
    nl3 = lbound(pin, 3)
    nu3 = ubound(pin, 3)
    allocate(pout(nl1:nu1, nl2:nu2, nl3:nu3))
    do i3 = nl3, nu3
      do i2 = nl2, nu2
        do i1 = nl1, nu1
          pout(i1, i2, i3) = pin(i1, i2, i3)
        end do
      end do
    end do
  else
    nullify(pout)
  end if

  POP_SUB(SUBNAME(loct_pointer_copy_3))
end subroutine SUBNAME(loct_pointer_copy_3)


! ---------------------------------------------------------
subroutine SUBNAME(loct_pointer_copy_4)(pout, pin)
  TYPE, pointer, intent(out) :: pout(:, :, :, :)
  TYPE, pointer, intent(in)  :: pin(:, :, :, :)

  integer :: nl1, nu1, nl2, nu2, nl3, nu3, nl4, nu4
  integer :: i1, i2, i3, i4
  
  PUSH_SUB(SUBNAME(loct_pointer_copy_4))

  if(associated(pin)) then
    nl1 = lbound(pin, 1)
    nu1 = ubound(pin, 1)
    nl2 = lbound(pin, 2)
    nu2 = ubound(pin, 2)
    nl3 = lbound(pin, 3)
    nu3 = ubound(pin, 3)
    nl4 = lbound(pin, 4)
    nu4 = ubound(pin, 4)
    allocate(pout(nl1:nu1, nl2:nu2, nl3:nu3, nl4:nu4))
    do i4 = nl4, nu4
      do i3 = nl3, nu3
        do i2 = nl2, nu2
          do i1 = nl1, nu1
            pout(i1, i2, i3, i4) = pin(i1, i2, i3, i4)
          end do
        end do
      end do
    end do
  else
    nullify(pout)
  end if

  POP_SUB(SUBNAME(loct_pointer_copy_4))
end subroutine SUBNAME(loct_pointer_copy_4)

! ---------------------------------------------------------
subroutine SUBNAME(loct_allocatable_copy_1)(pout, pin)
  TYPE, allocatable, intent(out) :: pout(:)
  TYPE, allocatable, intent(in)  :: pin(:)

  integer :: nl1, nu1
  integer :: i1

  PUSH_SUB(SUBNAME(loct_allocatable_copy_1))

  if(allocated(pin)) then
    nl1 = lbound(pin, 1)
    nu1 = ubound(pin, 1)
    allocate(pout(nl1:nu1))
    do i1 = nl1, nu1
      pout(i1) = pin(i1)
    end do
  end if

  POP_SUB(SUBNAME(loct_allocatable_copy_1))
end subroutine SUBNAME(loct_allocatable_copy_1)

! ---------------------------------------------------------
subroutine SUBNAME(loct_allocatable_copy_2)(pout, pin)
  TYPE, allocatable, intent(out) :: pout(:, :)
  TYPE, allocatable, intent(in)  :: pin(:, :)

  integer :: nl1, nu1, nl2, nu2
  integer :: i1, i2

  PUSH_SUB(SUBNAME(loct_allocatable_copy_2))

  if(allocated(pin)) then
    nl1 = lbound(pin, 1)
    nu1 = ubound(pin, 1)
    nl2 = lbound(pin, 2)
    nu2 = ubound(pin, 2)
    allocate(pout(nl1:nu1, nl2:nu2))
    do i2 = nl2, nu2
      do i1 = nl1, nu1
        pout(i1, i2) = pin(i1, i2)
      end do
    end do
  end if

  POP_SUB(SUBNAME(loct_allocatable_copy_2))
end subroutine SUBNAME(loct_allocatable_copy_2)


! ---------------------------------------------------------
subroutine SUBNAME(loct_allocatable_copy_3)(pout, pin)
  TYPE, allocatable, intent(out) :: pout(:, :, :)
  TYPE, allocatable, intent(in)  :: pin(:, :, :)

  integer :: nl1, nu1, nl2, nu2, nl3, nu3
  integer :: i1, i2, i3

  PUSH_SUB(SUBNAME(loct_allocatable_copy_3))

  if(allocated(pin)) then
    nl1 = lbound(pin, 1)
    nu1 = ubound(pin, 1)
    nl2 = lbound(pin, 2)
    nu2 = ubound(pin, 2)
    nl3 = lbound(pin, 3)
    nu3 = ubound(pin, 3)
    allocate(pout(nl1:nu1, nl2:nu2, nl3:nu3))
    do i3 = nl3, nu3
      do i2 = nl2, nu2
        do i1 = nl1, nu1
          pout(i1, i2, i3) = pin(i1, i2, i3)
        end do
      end do
    end do
  end if

  POP_SUB(SUBNAME(loct_allocatable_copy_3))
end subroutine SUBNAME(loct_allocatable_copy_3)


! ---------------------------------------------------------
subroutine SUBNAME(loct_allocatable_copy_4)(pout, pin)
  TYPE, allocatable, intent(out) :: pout(:, :, :, :)
  TYPE, allocatable, intent(in)  :: pin(:, :, :, :)

  integer :: nl1, nu1, nl2, nu2, nl3, nu3, nl4, nu4
  integer :: i1, i2, i3, i4
  
  PUSH_SUB(SUBNAME(loct_allocatable_copy_4))

  if(allocated(pin)) then
    nl1 = lbound(pin, 1)
    nu1 = ubound(pin, 1)
    nl2 = lbound(pin, 2)
    nu2 = ubound(pin, 2)
    nl3 = lbound(pin, 3)
    nu3 = ubound(pin, 3)
    nl4 = lbound(pin, 4)
    nu4 = ubound(pin, 4)
    allocate(pout(nl1:nu1, nl2:nu2, nl3:nu3, nl4:nu4))
    do i4 = nl4, nu4
      do i3 = nl3, nu3
        do i2 = nl2, nu2
          do i1 = nl1, nu1
            pout(i1, i2, i3, i4) = pin(i1, i2, i3, i4)
          end do
        end do
      end do
    end do
  end if

  POP_SUB(SUBNAME(loct_allocatable_copy_4))
end subroutine SUBNAME(loct_allocatable_copy_4)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
