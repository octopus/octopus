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
AC_LINK_IFELSE( AC_LANG_PROGRAM( [
#include <emmintrin.h>
], [
__m128d a __attribute__((aligned(16)));
 ]), 
 [AC_DEFINE(HAVE_M128D, 1, [compiler supports the m128d type]) [acx_m128d=yes]], [acx_m128d=no])
AC_MSG_RESULT($acx_m128d)])

################################################################
# Check whether the compiler accepts Blue Gene extensions
# ----------------------------------
AC_DEFUN([ACX_BLUE_GENE],
[AC_MSG_CHECKING([for Blue Gene intrinsics])
AC_LINK_IFELSE( AC_LANG_PROGRAM( [
], [[
  double aa, bb;
  double _Complex cc;

  cc = __cmplx(aa, bb);
  cc = __fpneg(cc);
 ]]), 
 [AC_DEFINE(HAVE_BLUE_GENE, 1, [compiler supports Blue Gene intrinsics]) [acx_blue_gene=yes]], [acx_blue_gene=no])
AC_MSG_RESULT($acx_blue_gene)])

#################################################################
# Enables architecture-specific code
AC_DEFUN([ACX_ARCH],
[
AC_REQUIRE([AC_CANONICAL_HOST])

vector=no
assembler=no
blue_gene=no

case "${host}" in
x86_64*)
ACX_M128D
oct_arch=x86_64
vector=$acx_m128d
vector_type="(sse2)"
assembler=no
AC_DEFINE(OCT_ARCH_X86_64, 1, [This an x86_64 system])
;;
i?86*)
ACX_M128D
vector=$acx_m128d
oct_arch=x86
vector_type="(sse2)"
if test x$vector = xyes ; then
# We allow explicit disabling of SSE2
ac_enable_sse2=no
AC_ARG_ENABLE(sse2, AS_HELP_STRING([--enable-sse2], [Enable the use of SSE2 instructions (x86)]), 
	[ac_enable_sse2=${enableval}])

if test x"${ac_enable_sse2}" = x"no"; then
vector=disabled
fi
fi
AC_DEFINE(OCT_ARCH_X86_64, 1, [This an x86 system])
;;	
ia64*)
oct_arch=ia64
assembler=yes
AC_DEFINE(OCT_ARCH_IA64, 1, [This an Itanium system])
;;
sparc*)
oct_arch=sparc
AC_DEFINE(OCT_ARCH_SPARC, 1, [This a Sparc system])
;;
alphaev*)
oct_arch=alpha
AC_DEFINE(OCT_ARCH_ALPHA, 1, [This an Alpha system])
;;
mips*)
oct_arch=mips
AC_DEFINE(OCT_ARCH_MIPS, 1, [This a MIPS system])
;;
powerpc*)
oct_arch=powerpc
AC_DEFINE(OCT_ARCH_POWERPC, 1, [This a PowerPC system])
ACX_BLUE_GENE
blue_gene=$acx_blue_gene
vector=$acx_blue_gene
vector_type="(bg)"
;;
*)
oct_arch=unknown
;;
esac

ac_enable_sse2=yes
AC_ARG_ENABLE(sse2, AS_HELP_STRING([--disable-sse2], [Disable the use of SSE2 instructions (x86_64)]), 
	[ac_enable_sse2=${enableval}])

if test x"${ac_enable_sse2}" = x"no"; then
vector=disabled
fi

# We allow explicit disabling of assembler
ac_enable_assembler=yes
AC_ARG_ENABLE(assembler, AS_HELP_STRING([--disable-assembler], [Disable the use of assembler code (ia64)]), 
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

AC_MSG_NOTICE([ Architecture-specific code:
***************************
This is a $oct_arch processor:
vectorial code: $vector $vector_type
assembler code: $assembler
***************************])
])
