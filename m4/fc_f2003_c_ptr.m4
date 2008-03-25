## Copyright (C) 2008 T. Burnus
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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id$
##

#
# Check for Fortran 2003 c_ptr support
# ------------------------------------

AC_DEFUN([ACX_FC_F2003_C_PTR],
[
AC_MSG_CHECKING([for Fortran 2003 c_ptr type])

AC_LINK_IFELSE(
    AC_LANG_PROGRAM( [], [
    use iso_c_binding
    implicit none
    type(c_ptr) :: ptr
    ptr = c_null_ptr
    if (c_associated(ptr)) stop 3
    ]),
    [AC_MSG_RESULT(yes)
     AC_DEFINE(F2003_C_PTR, 1,
               [compiler supports c_ptr type of Fortran 2003])],
    [AC_MSG_RESULT(no)])

]
)
