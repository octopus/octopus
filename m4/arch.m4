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
AC_LINK_IFELSE([AC_LANG_PROGRAM( [
#include <emmintrin.h>
], [
__m128d a __attribute__((aligned(16)));
 ])], 
 [AC_DEFINE(HAVE_M128D, 1, [compiler supports the m128d type]) [acx_m128d=yes]], [acx_m128d=no])
AC_MSG_RESULT($acx_m128d)])

################################################################
# Check whether the compiler accepts the __m256d type
# ----------------------------------
AC_DEFUN([ACX_M256D],
[AC_MSG_CHECKING([whether the compiler accepts the __m256d type])
AC_LINK_IFELSE([AC_LANG_PROGRAM( [
#include <immintrin.h>
], [
__m256d a __attribute__((aligned(32)));
 ])], 
 [acx_m256d=yes], [acx_m256d=no])
AC_MSG_RESULT($acx_m256d)])

################################################################
# Check whether the hardware accepts AVX instructions
# ----------------------------------
AC_DEFUN([ACX_AVX],
[AC_MSG_CHECKING([whether AVX instructions can be used])
acx_save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS"
AC_RUN_IFELSE([AC_LANG_PROGRAM( [
#include <immintrin.h>
], [
changequote(,)
__m256d a __attribute__((aligned(32)));
__m256d b __attribute__((aligned(32)));
__m256d c __attribute__((aligned(32)));
double d[4];

a = _mm256_add_pd(b, c);
_mm256_storeu_pd(d, a);
printf("",  *d);
changequote([, ])
 ])], 
 [acx_m256d=yes], [acx_m256d=no], [acx_m256d=no;echo -n "cross-compiling; assuming... "])
CFLAGS="$acx_save_CFLAGS"
AC_MSG_RESULT($acx_m256d)])

################################################################
# Check whether the hardware accepts FMA4 instructions
# ----------------------------------
AC_DEFUN([ACX_FMA4],
[AC_MSG_CHECKING([whether FMA4 instructions can be used])
acx_save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS"
AC_RUN_IFELSE([AC_LANG_PROGRAM( [
#include <x86intrin.h>
], [
__m128d a __attribute__((aligned(16)));
__m128d b __attribute__((aligned(16)));
__m128d c __attribute__((aligned(16)));
__m128d d __attribute__((aligned(16)));
d = (__m128d) _mm_macc_pd(a, b, c);
 ])], 
 [acx_fma4=yes], [acx_fma4=no], [acx_fma4=no;echo -n "cross-compiling; assuming... "])
CFLAGS="$acx_save_CFLAGS"
AC_MSG_RESULT($acx_fma4)])

################################################################
# Check whether the compiler accepts Blue Gene extensions
# ----------------------------------
AC_DEFUN([ACX_BLUE_GENE],
[AC_MSG_CHECKING([for Blue Gene intrinsics])
AC_LINK_IFELSE([AC_LANG_PROGRAM( [
], [[
  double aa, bb;
  double _Complex cc;

  cc = __cmplx(aa, bb);
  cc = __fpneg(cc);
 ]])], 
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
##########################################
x86_64*|*apple-darwin*) #workaround for a bug in autoconf/OS X

oct_arch=x86_64
assembler=no
AC_DEFINE(OCT_ARCH_X86_64, 1, [This is an x86_64 system])

#SSE2
ACX_M128D
vector=$acx_m128d
vector_type="(sse2)"

#FMA4
AC_ARG_ENABLE(fma4, AS_HELP_STRING([--enable-fma4], [Enable the use of FMA4 vectorial instructions (x86_64)]), 
	[ac_enable_fma4=${enableval}])
if test "x$vector" = "xno" ; then
 ac_enable_fma4=no
fi
if test "x$ac_enable_fma4" = "x" ; then
  ACX_FMA4
elif test "x$ac_enable_fma4" = "xyes" ; then
  AC_MSG_NOTICE([FMA4 instruction support enabled])
  acx_fma4=yes
else # no
  AC_MSG_NOTICE([FMA4 instruction support disabled])
  acx_fma4=no
fi
if test "x$acx_fma4" = "xyes" ; then
  AC_DEFINE(HAVE_FMA4, 1, [compiler and hardware supports the FMA4 instructions])
fi

#AVX
AC_ARG_ENABLE(avx, AS_HELP_STRING([--enable-avx], [Enable the use of AVX vectorial instructions (x86_64)]), 
	[ac_enable_avx=${enableval}])
if test "x$vector" = "xno" ; then
 ac_enable_avx=no
fi
if test "x$ac_enable_avx" = "x" ; then
  ACX_M256D
  if test "x$acx_m256d" = "xyes" ; then
    ACX_AVX
  fi
elif test "x$ac_enable_avx" = "xyes" ; then
  AC_MSG_NOTICE([AVX instruction support enabled])
  acx_m256d=yes
else # no
  AC_MSG_NOTICE([AVX instruction support disabled])
  acx_m256d=no
fi
if test "x$acx_m256d" = "xyes" ; then
  AC_DEFINE(HAVE_M256D, 1, [compiler and hardware support the m256d type and AVX instructions])
  vector=$acx_m256d
  vector_type="(avx)"
fi
;;
##########################################
i?86*)
ACX_M128D
vector=$acx_m128d
oct_arch=x86
vector_type="(sse2)"
if test "x$vector" = "xyes" ; then
# We allow explicit disabling of SSE2
ac_enable_vectors=no
AC_ARG_ENABLE(vectors, AS_HELP_STRING([--enable-vectors], [Enable the use of vectorial instructions (x86)]), 
	[ac_enable_vectors=${enableval}])

if test x"${ac_enable_vectors}" = x"no"; then
vector=disabled
fi
fi
AC_DEFINE(OCT_ARCH_X86, 1, [This is an x86 system])
;;	
##########################################
ia64*)
oct_arch=ia64
AC_DEFINE(OCT_ARCH_IA64, 1, [This is an Itanium system])
;;
##########################################
sparc*)
oct_arch=sparc
AC_DEFINE(OCT_ARCH_SPARC, 1, [This is a Sparc system])
;;
##########################################
alphaev*)
oct_arch=alpha
AC_DEFINE(OCT_ARCH_ALPHA, 1, [This is an Alpha system])
;;
##########################################
mips*)
oct_arch=mips
AC_DEFINE(OCT_ARCH_MIPS, 1, [This is a MIPS system])
;;
##########################################
powerpc*)
oct_arch=powerpc
AC_DEFINE(OCT_ARCH_POWERPC, 1, [This is a PowerPC system])
ACX_BLUE_GENE
blue_gene=$acx_blue_gene
vector=$acx_blue_gene
vector_type="(bg)"
;;
##########################################
*)
oct_arch=unknown
;;
esac

ac_enable_vectors=yes
AC_ARG_ENABLE(vectors, AS_HELP_STRING([--disable-vectors], [Disable the use of vectorial instructions (x86_64 and Blue Gene)]), 
	[ac_enable_vectors=${enableval}])

if test x"${ac_enable_vectors}" = x"no"; then
vector=disabled
fi

AC_DEFINE_UNQUOTED(OCT_ARCH, $oct_arch, [The architecture of this system])

AM_CONDITIONAL(COMPILE_VEC, test "x$vector" = "xyes")

if test "x$vector = xyes" ; then
AC_DEFINE(HAVE_VEC, 1, [Define to 1 if vectorial routines are to be compiled])
fi

AC_MSG_NOTICE([Architecture-specific code:
***************************
This is a $oct_arch processor:
vectorial code: $vector $vector_type
***************************])
])
