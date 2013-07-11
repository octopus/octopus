## Copyright (C) 2009 X. Andrade
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$
##

#
# Check for Fortran 2003 c_ptr support
# ------------------------------------

AC_DEFUN([ACX_FC_SIZEOF],
[
AC_MSG_CHECKING([for Fortran sizeof intrinsic])

AC_LINK_IFELSE(
    AC_LANG_PROGRAM( [], [[
    implicit none
    integer s
    integer array(10)
    real, allocatable :: alc(:)
    character(len=10) C
    type dtype
      integer :: array(10)
      integer, pointer :: p(:)
    end type
    type (dtype) dobj
    s = sizeof(array)
    s = sizeof(array(3))
    s = sizeof(alc)
    s = sizeof(c(2:5))
    s = sizeof(c)
    s = sizeof(dobj)
    s = sizeof(dobj%array)
    ]]),
    [AC_MSG_RESULT(yes)
     AC_DEFINE(HAVE_FC_SIZEOF, 1,
               [Fortran compiler supports the sizeof intrinsic])],
    [AC_MSG_RESULT(no)]
  )
])
