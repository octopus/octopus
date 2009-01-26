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

  subroutine SUBNAME(loct_pointer_copy_1)(o, i)
    TYPE, pointer :: o(:), i(:)
    integer :: nl, nu
    if(associated(i)) then
      nl = lbound(i, 1)
      nu = ubound(i, 1)
      allocate(o(nl:nu))
      o = i
    end if
  end subroutine SUBNAME(loct_pointer_copy_1)

  subroutine SUBNAME(loct_pointer_copy_2)(o, i)
    TYPE, pointer :: o(:, :), i(:, :)
    integer :: nl1, nu1, nl2, nu2
    if(associated(i)) then
      nl1 = lbound(i, 1)
      nu1 = ubound(i, 1)
      nl2 = lbound(i, 2)
      nu2 = ubound(i, 2)
      allocate(o(nl1:nu1, nl2:nu2))
      o = i
      end if
  end subroutine SUBNAME(loct_pointer_copy_2)

  subroutine SUBNAME(loct_pointer_copy_3)(o, i)
    TYPE, pointer :: o(:, :, :), i(:, :, :)
    integer :: nl1, nu1, nl2, nu2, nl3, nu3
    if(associated(i)) then
      nl1 = lbound(i, 1)
      nu1 = ubound(i, 1)
      nl2 = lbound(i, 2)
      nu2 = ubound(i, 2)
      nl3 = lbound(i, 3)
      nu3 = ubound(i, 3)
      allocate(o(nl1:nu1, nl2:nu2, nl3:nu3))
      o = i
    end if
  end subroutine SUBNAME(loct_pointer_copy_3)

  subroutine SUBNAME(loct_pointer_copy_4)(o, i)
    TYPE, pointer :: o(:, :, :, :), i(:, :, :, :)
    integer :: nl1, nu1, nl2, nu2, nl3, nu3, nl4, nu4
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
      o = i
    end if
  end subroutine SUBNAME(loct_pointer_copy_4)
