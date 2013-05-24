## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

################################################
# Check whether the compiler accepts very long lines.
# ----------------------------------
AC_DEFUN([ACX_LONG_FORTRAN_LINES],
[AC_MSG_CHECKING([whether the compiler accepts very long lines])
AC_COMPILE_IFELSE( AC_LANG_PROGRAM( [], [
write(*, *) '456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678904567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789001234567890123456789045678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890'
 ]), 
 [acx_long_lines_ok=yes; AC_DEFINE(LONG_LINES, 1, [compiler supports long lines])], [acx_long_lines_ok=no])
AC_SUBST([LONG_LINES], [$acx_long_lines_ok])
AC_MSG_RESULT($acx_long_lines_ok)
])


################################################
# Check whether the compiler accepts preprocessor "# line-number" lines.
# ----------------------------------
AC_DEFUN([ACX_F90_ACCEPTS_LINE_NUMBERS],
[
AC_MSG_CHECKING([whether the compiler accepts "line-number" lines cast by the preprocessor])
AC_COMPILE_IFELSE(
    AC_LANG_PROGRAM( [], [# 1]),
    [acx_f90_accepts_line_numbers_ok=yes
    AC_DEFINE(F90_ACCEPTS_LINE_NUMBERS, 1, [compiler supports line-number lines])],
    [acx_f90_accepts_line_numbers_ok=no])
AC_SUBST(F90_ACCEPTS_LINE_NUMBERS, $acx_f90_accepts_line_numbers_ok)
AC_MSG_RESULT($acx_f90_accepts_line_numbers_ok)
]
)

################################################
# Check for the presence of a given function in Fortran.
# It substitutes AC_CHECK_FUNC, since the latter
# seems to fail with some autotools versions, due to a call to some broken
# version of AC_LANG_FUNC_LINK_TRY.
AC_DEFUN([ACX_FORTRAN_CHECK_FUNC],
[
AC_MSG_CHECKING([for $1])
AC_LANG_PUSH(Fortran)dnl
AC_LINK_IFELSE([AC_LANG_CALL([], [$1])], 
[
acx_fortran_check_func=yes
AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_$1]),1, [Define if the $1 function can be called from Fortran])], 
[
acx_fortran_check_func=no
])dnl
AC_LANG_POP(Fortran)dnl
AC_MSG_RESULT($acx_fortran_check_func)
])


################################################
# AC_LANG_FUNC_LINK_TRY(Fortran)(FUNCTION)
# ----------------------------------
m4_define([AC_LANG_FUNC_LINK_TRY(Fortran)],
[AC_LANG_PROGRAM([], [call [$1]])])

AC_DEFUN([ACX_CHECK_CPP],
[
     for CPP in "$CPP" "$CPP -ansi"; do
         AC_MSG_CHECKING([whether $CPP is usable for Fortran preprocessing])
	 acx_fpp_ok=yes

#      	 AC_EGREP_CPP([hi], AC_LANG_PROGRAM([],[
#@%:@define ADD_I(x) x @%:@@%:@ i
#ADD_I(h)]),
#	   [], [acx_fpp_ok=no; AC_MSG_RESULT([preprocessor does not concatenate tokens])])

      	 AC_EGREP_CPP([hi], AC_LANG_PROGRAM([],[
#define ADD_I(x) x ## i
ADD_I(h)]),
	   [], [acx_fpp_ok=no; AC_MSG_RESULT([preprocessor does not concatenate tokens])])

         # in Fortran this is string concatenation, must not be stripped
         AC_EGREP_CPP([string2], [string1 // string2],
	   [], [acx_fpp_ok=no; AC_MSG_RESULT([preprocessor strips C++ style comment])])

	if test x"$acx_fpp_ok" = xyes; then
          AC_MSG_RESULT([yes])
	  break
	fi
     done

     if test x"$acx_fpp_ok" = xno; then
     	AC_MSG_ERROR([Preprocessor is not usable for Fortran.])
     fi
])

#  for my_cpp in "$FCCPP" "$CPP" "`which cpp` -ansi" "/lib/cpp -ansi"; do


###############################################
# ACX_FORTRAN_LOC checks for the presence of the loc intrinsics
# --------------------------
AC_DEFUN([ACX_FORTRAN_LOC], [
    AC_MSG_CHECKING(for loc)
    AC_LINK_IFELSE(
    [
    AC_LANG_PROGRAM(,
      [	
          implicit none
      
          integer :: ii, jj

          ii = loc(jj)
      ])
    ],
    [
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_FORTRAN_LOC, 1,
                [Define if the compiler provides the loc instrinsic])
    ],
    [
      AC_MSG_RESULT(no)
    ])
])
