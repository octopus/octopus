## Copyright (C) 2008 X. Andrade
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
## $Id: fcflags.m4 4621 2008-10-08 17:21:03Z xavier $
##
################################################

################################################################
# Check whether the compiler accepts the __m128d type
# ----------------------------------
AC_DEFUN([ACX_M128D],
[AC_MSG_CHECKING([whether the compiler accepts the __m128d type])
AC_COMPILE_IFELSE( AC_LANG_PROGRAM( [
#include <emmintrin.h>
], [
__m128d a __attribute__((aligned(16)));
 ]), 
 [AC_DEFINE(HAVE_M128D, 1, [compiler supports the m128d type]) [acx_m128d=yes]], [acx_m128d=no])
AC_MSG_RESULT($acx_m128d)])

#################################################################
# Enables architecture specific code
AC_DEFUN([ACX_ARCH],
[
AC_REQUIRE([AC_CANONICAL_HOST])

vector=no
assembler=no

case "${host}" in
x86_64*)
ACX_M128D
oct_arch=x86_64
vector=$acx_m128d
;;
i?86*)
ACX_M128D
vector=$acx_m128d
oct_arch=x86
;;	
ia64*)
oct_arch=ia64
assembler=yes
;;
sparc*)
oct_arch=sparc
;;
alphaev*)
oct_arch=alpha
;;
mips*)
oct_arch=mips
;;
powerpc64*)
oct_arch=powerpc64
;;
powerpc*)
oct_arch=powerpc
;;
*)
oct_arch=unknown
;;
esac

# We allow explicit disabling of SSE2
ac_enable_sse2=yes
AC_ARG_ENABLE(sse2, AS_HELP_STRING([--disable-sse2], [Disable the use of SSE2 instructions]), 
	[ac_enable_sse2=${enableval}])

if test x"${ac_enable_sse2}" = x"no"; then
vector=disabled
fi

# We allow explicit disabling of assembler
ac_enable_assembler=yes
AC_ARG_ENABLE(assembler, AS_HELP_STRING([--disable-assembler], [Disable the use of assembler code]), 
	[ac_enable_assembler=${enableval}])

if test x"${ac_enable_assembler}" = x"no"; then
assembler=disabled
fi

AC_DEFINE_UNQUOTED(OCT_ARCH, $oct_arch, [The architecture of this system])

AM_CONDITIONAL(COMPILE_AS, test x$assembler = xyes)
AM_CONDITIONAL(COMPILE_VEC, test x$vector = xyes)

if test x$vector = xyes ; then
AC_DEFINE(HAVE_VEC, 1, [Define to 1 if vectorial routines are to be compiled])
fi

if test x$assembler = xyes ; then
AC_DEFINE(HAVE_AS, 1, [Define to 1 if assembler routines are to be compiled])
fi

AC_MSG_NOTICE([ Architecture specific code:
***************************
This is a $oct_arch processor:
vectorial code: $vector 
assembler code: $assembler
***************************])
])
