/*
 Copyright (C) 2006 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: operate.c 2146 2006-05-23 17:36:00Z xavier $
*/

#ifndef OCTOPUS_BEAK_H
#define OCTOPUS_BEAK_H

#include <config.h>

/* If __builtin_prefetch is not present (which should have been caught
   by the configure script) one needs to define a dummy preprocessor
   macro. */

#if !defined(HAVE_BUILTIN_PREFETCH)
#define __builtin_prefetch(a, b, c)
#endif

#ifdef SINGLE_PRECISION
typedef float ffloat;
#else
typedef double ffloat;
#endif

typedef struct {
  ffloat re;
  ffloat im;
} comp;

/* These constants have to match the definition in
   src/grid/nl_operator.F90 */

#define OP_FORTRAN 0
#define OP_C       1
#define OP_VEC     2
#define OP_AS      3

#define M_REAL     1
#define M_CMPLX    2

#define MAX_OP_N   100

/* HERE WE DETECT THE PROCESSOR TYPE AND VECTORIAL CAPABILITIES */

/* If this is a x86_64 machine we always can use vectors */
#if defined (__amd64__) || defined(__x86_64__)
#if defined (__SSE2__) && defined(HAVE_EMMINTRIN_H) && defined(FC_USES_MALLOC)
# define OCT_AMD64
# define USE_VECTORS
#endif
#endif

/* Itanium only has single precision (SSE) */
#if defined (__ia64__) || defined(__ia64) || defined(_M_IA64)
# define OCT_ITANIUM
# if defined(SINGLE_PRECISION) && defined(HAVE_XMMINTRIN_H)
#  define USE_VECTORS
# endif
#endif

/* x86 is the most complex case */
#if defined(__i386__) || defined(__i386) || defined(_M_IX86) || defined(_X86_)
# define OCT_X86

# if defined(SINGLE_PRECISION)

/* we only need sse */
#  if defined (__SSE__) || defined(__SSE2__) || defined(HAVE_XMMINTRIN_H)
#   define USE_VECTORS
#  endif

# else /* DOUBLE_PRECISION */

/* we need SSE2 and aligned memory */
#  if defined(__SSE2__) && defined(HAVE_EMMINTRIN_H) && defined(FC_USES_MALLOC)
#   if defined(HAVE_16_BYTES_ALIGNED_MALLOC)
#    define USE_VECTORS
#   else /* not HAVE_16_BYTES_ALIGNED_MALLOC */
#    if defined(HAVE_POSIX_MEMALIGN)
#     define USE_VECTORS
#     define USE_FAKE_MALLOC
#    endif
#   endif
#  endif

# endif
#endif /* x86 */

#endif /* OCTOPUS_BEAK_H */
