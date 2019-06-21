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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
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
 [acx_m128d=yes], [acx_m128d=no])
AC_MSG_RESULT($acx_m128d)])

################################################################
# Check whether the hardware accepts SSE2 instructions
# Note: printf in program is to avoid vector intrinsics being optimized out
# ----------------------------------
AC_DEFUN([ACX_SSE2],
[AC_MSG_CHECKING([whether SSE2 instructions can be used])
acx_save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS"
AC_RUN_IFELSE([AC_LANG_PROGRAM([
#include <emmintrin.h>
#include <stdio.h>
], [
changequote(,)
__m128d a __attribute__((aligned(16)));
__m128d b __attribute__((aligned(16)));
__m128d c __attribute__((aligned(16)));
double d[4];

a = _mm_add_pd(b, c);
_mm_storeu_pd(d, a);
printf("",  *d);
changequote([, ])
 ])], 
 [acx_m128d=yes], [acx_m128d=no], [acx_m128d=yes;echo -n "cross-compiling; assuming... "])
# assume yes (rather than no as for FMA4 and AVX) since SSE2 is very common, especially when the compiler has m128d
CFLAGS="$acx_save_CFLAGS"
AC_MSG_RESULT($acx_m128d)])

################################################################
# Routine for checking on SSE2 generally
# ----------------------------------
AC_DEFUN([ACX_SSE2_WRAPPER],
[AC_ARG_ENABLE(sse2, AS_HELP_STRING([--enable-sse2], [Enable the use of SSE2 vectorial instructions]),
	[ac_enable_sse2=${enableval}])
if test "x$vector" = "xno" ; then
 ac_enable_sse2=no
fi
if test "x$ac_enable_sse2" = "x" ; then
  ACX_M128D
  if test "x$acx_m128d" = "xyes" ; then
    ACX_SSE2
  fi
elif test "x$ac_enable_sse2" = "xyes" ; then
  AC_MSG_NOTICE([SSE2 instruction support enabled])
  acx_m128d=yes
else # no
  AC_MSG_NOTICE([SSE2 instruction support disabled])
  acx_m128d=no
fi
if test "x$acx_m128d" = "xyes" ; then
  AC_DEFINE(HAVE_M128D, 1, [compiler and hardware support the m128d type and SSE2 instructions])
  vector=$acx_m128d
  vector_type="(sse2)"
fi])

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
# Note: printf in program is to avoid vector intrinsics being optimized out
# ----------------------------------
AC_DEFUN([ACX_AVX],
[AC_MSG_CHECKING([whether AVX instructions can be used])
acx_save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS"
AC_RUN_IFELSE([AC_LANG_PROGRAM([
#include <immintrin.h>
#include <stdio.h>
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
# Check whether the hardware accepts FMA3 instructions
# ----------------------------------
AC_DEFUN([ACX_FMA3],
[AC_MSG_CHECKING([whether FMA3 instructions can be used])
acx_save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS"
AC_RUN_IFELSE([AC_LANG_PROGRAM( [
#include <x86intrin.h>
#include <stdio.h>
], [
changequote(,)
__m256d a __attribute__((aligned(32)));
__m256d b __attribute__((aligned(32)));
__m256d c __attribute__((aligned(32)));
__m256d d __attribute__((aligned(32)));
double e[4];
d = (__m256d) _mm256_fmadd_pd(a, b, c);
_mm256_storeu_pd(e, d);
printf("",  *e);
changequote([, ])
 ])],
 [acx_fma3=yes], [acx_fma3=no], [acx_fma3=no;echo -n "cross-compiling; assuming... "])
CFLAGS="$acx_save_CFLAGS"
AC_MSG_RESULT($acx_fma3)])

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

################################################################
# Check whether the compiler accepts Blue Gene extensions
# ----------------------------------
AC_DEFUN([ACX_BLUE_GENE_Q],
[AC_MSG_CHECKING([for Blue Gene Q intrinsics])
AC_LINK_IFELSE([AC_LANG_PROGRAM( [
], [[
  vector4double aa, bb, cc;

  cc = vec_madd(aa, bb, cc);
 ]])], 
 [AC_DEFINE(HAVE_BLUE_GENE_Q, 1, [compiler supports Blue Gene Q intrinsics]) [acx_blue_gene_q=yes]], [acx_blue_gene_q=no])
AC_MSG_RESULT($acx_blue_gene_q)])

################################################################
# Check whether the hardware accepts AVX512 instructions
# ----------------------------------
AC_DEFUN([ACX_AVX512],
[AC_MSG_CHECKING([whether AVX512 instructions can be used])
acx_save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS"
AC_RUN_IFELSE([AC_LANG_PROGRAM( [
#include <stdio.h>
#include <immintrin.h>
], [
changequote(,)
__m512d a __attribute__((aligned(64)));
__m512d b __attribute__((aligned(64)));
__m512d c __attribute__((aligned(64)));
__m512d d __attribute__((aligned(64)));
double e[8];
d = (__m512d) _mm512_fmadd_pd(a, b, c);
_mm512_storeu_pd(e, d);
printf("",  *e);
changequote([, ])
 ])],
 [acx_avx512=yes], [acx_avx512=no], [acx_avx512=no;echo -n "cross-compiling; assuming... "])
CFLAGS="$acx_save_CFLAGS"
AC_MSG_RESULT($acx_avx512)])


#################################################################
# Enables architecture-specific code
AC_DEFUN([ACX_ARCH],
[
AC_REQUIRE([AC_CANONICAL_HOST])

blue_gene=no
vector_type=""

ac_enable_vectors=yes
AC_ARG_ENABLE(vectors, AS_HELP_STRING([--enable-vectors], [Enable the use of vectorial instructions]),
	[ac_enable_vectors=${enableval}])

if test x"${ac_enable_vectors}" = x"no"; then
  vector=disabled
else
  vector=yes
fi

case "${host}" in
##########################################
x86_64*|*apple-darwin*) #workaround for a bug in autoconf/OS X
oct_arch=x86_64
;;
##########################################
i?86*)
oct_arch=x86
;;	
##########################################
ia64*)
oct_arch=ia64
# Itanium
;;
##########################################
sparc*)
oct_arch=sparc
;;
##########################################
alphaev*)
oct_arch=alpha
;;
##########################################
mips*)
oct_arch=mips
;;
##########################################
powerpc*)
oct_arch=powerpc
;;
##########################################
*)
oct_arch=unknown
;;
esac

if test x"${vector}" != x"disabled"; then
case "${oct_arch}" in
##########################################
x86_64)

#SSE2
ACX_SSE2_WRAPPER

#FMA3
AC_ARG_ENABLE(fma3, AS_HELP_STRING([--enable-fma3], [Enable the use of FMA3 vectorial instructions (x86_64)]),
	[ac_enable_fma3=${enableval}])
if test "x$vector" = "xno" ; then
 ac_enable_fma3=no
fi
if test "x$ac_enable_fma3" = "x" ; then
  ACX_FMA3
elif test "x$ac_enable_fma3" = "xyes" ; then
  AC_MSG_NOTICE([FMA3 instruction support enabled])
  acx_fma3=yes
else # no
  AC_MSG_NOTICE([FMA3 instruction support disabled])
  acx_fma3=no
fi
if test "x$acx_fma3" = "xyes" ; then
  AC_DEFINE(HAVE_FMA3, 1, [compiler and hardware supports the FMA3 instructions])
fi

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

#AVX512
AC_ARG_ENABLE(avx512, AS_HELP_STRING([--enable-avx512], [Enable the use of AVX512 vectorial instructions (x86_64)]),
	[ac_enable_avx512=${enableval}])
if test "x$vector" = "xno" ; then
 ac_enable_avx512=no
fi
if test "x$ac_enable_avx512" = "x" ; then
  ACX_AVX512
elif test "x$ac_enable_avx512" = "xyes" ; then
  AC_MSG_NOTICE([AVX instruction support enabled])
  acx_avx512=yes
else # no
  AC_MSG_NOTICE([AVX instruction support disabled])
  acx_avx512=no
fi
if test "x$acx_avx512" = "xyes" ; then
  AC_DEFINE(HAVE_M512D, 1, [compiler and hardware support the m512d type and AVX512 instructions])
  vector=$acx_avx512
  vector_type="(avx512)"
fi
;;
##########################################
x86)
ACX_SSE2_WRAPPER
;;
##########################################
powerpc)
ACX_BLUE_GENE
ACX_BLUE_GENE_Q
if test x$acx_blue_gene_q = x"yes"; then
  blue_gene=$acx_blue_gene_q
  vector=$acx_blue_gene_q
  vector_type="(blue gene/q)"
else
  blue_gene=$acx_blue_gene
  vector=$acx_blue_gene
  vector_type="(blue gene/p)"
fi
;;
esac
fi

AC_DEFINE_UNQUOTED(OCT_ARCH, $oct_arch, [The architecture of this system])

if test "x$vector" = "xyes" ; then
AC_DEFINE(HAVE_VEC, 1, [Define to 1 if vectorial routines are to be compiled])
fi

AC_MSG_NOTICE([Architecture-specific code:
***************************
This is a $oct_arch processor:
vectorial code: $vector $vector_type
***************************])
])
