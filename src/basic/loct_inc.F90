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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: loct.F90 4640 2008-10-10 08:56:33Z xavier $

! ---------------------------------------------------------
subroutine SUBNAME(loct_pointer_copy_1)(o, i)
  TYPE, pointer, intent(out) :: o(:)
  TYPE, pointer, intent(in)  :: i(:)

  integer :: nl1, nu1
  integer :: i1

  if(associated(i)) then
    nl1 = lbound(i, 1)
    nu1 = ubound(i, 1)
    allocate(o(nl1:nu1))
    forall (i1 = nl1:nu1) o(i1) = i(i1)
  else
    nullify(o)
  end if

end subroutine SUBNAME(loct_pointer_copy_1)


! ---------------------------------------------------------
subroutine SUBNAME(loct_pointer_copy_2)(o, i)
  TYPE, pointer, intent(out) :: o(:, :)
  TYPE, pointer, intent(in)  :: i(:, :)

  integer :: nl1, nu1, nl2, nu2
  integer :: i1, i2

  if(associated(i)) then
    nl1 = lbound(i, 1)
    nu1 = ubound(i, 1)
    nl2 = lbound(i, 2)
    nu2 = ubound(i, 2)
    allocate(o(nl1:nu1, nl2:nu2))
    forall (i1 = nl1:nu1, i2 = nl2:nu2) o(i1, i2) = i(i1, i2)
  else
    nullify(o)
  end if

end subroutine SUBNAME(loct_pointer_copy_2)


! ---------------------------------------------------------
subroutine SUBNAME(loct_pointer_copy_3)(o, i)
  TYPE, pointer, intent(out) :: o(:, :, :)
  TYPE, pointer, intent(in)  :: i(:, :, :)

  integer :: nl1, nu1, nl2, nu2, nl3, nu3
  integer :: i1, i2, i3

  if(associated(i)) then
    nl1 = lbound(i, 1)
    nu1 = ubound(i, 1)
    nl2 = lbound(i, 2)
    nu2 = ubound(i, 2)
    nl3 = lbound(i, 3)
    nu3 = ubound(i, 3)
    allocate(o(nl1:nu1, nl2:nu2, nl3:nu3))
    forall (i1 = nl1:nu1, i2 = nl2:nu2, i3 = nl3:nu3) o(i1, i2, i3) = i(i1, i2, i3)
  else
    nullify(o)
  end if

end subroutine SUBNAME(loct_pointer_copy_3)


! ---------------------------------------------------------
subroutine SUBNAME(loct_pointer_copy_4)(o, i)
  TYPE, pointer, intent(out) :: o(:, :, :, :)
  TYPE, pointer, intent(in)  :: i(:, :, :, :)

  integer :: nl1, nu1, nl2, nu2, nl3, nu3, nl4, nu4
  integer :: i1, i2, i3, i4
  
  if(associated(i)) then
    nl1 = lbound(i, 1)
    nu1 = ubound(i, 1)
    nl2 = lbound(i, 2)
    nu2 = ubound(i, 2)
    nl3 = lbound(i, 3)
    nu3 = ubound(i, 3)
    nl4 = lbound(i, 4)
    nu4 = ubound(i, 4)
    allocate(o(nl1:nu1, nl2:nu2, nl3:nu3, nl4:nu4))
    forall (i1 = nl1:nu1, i2 = nl2:nu2, i3 = nl3:nu3, i4 = nl4:nu4) o(i1, i2, i3, i4) = i(i1, i2, i3, i4)
  else
    nullify(o)
  end if

end subroutine SUBNAME(loct_pointer_copy_4)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
