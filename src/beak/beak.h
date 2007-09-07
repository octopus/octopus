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
   by the configure script) one needs to define dummy a preprocessor
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

void FC_FUNC_(zoperate_c,ZOPERATE_C)(const int * opnp, 
				     const int * opn, 
				     const ffloat * restrict w, 
				     const int * opi, 
				     const comp * fi, 
				     comp * restrict fo);

#if defined(HAVE_C_SSE2) && defined(HAVE_EMMINTRIN_H) && defined(FC_USES_MALLOC)

#if defined(HAVE_16_BYTES_ALIGNED_MALLOC)

#define USE_VECTORS

#else /* not HAVE_16_BYTES_ALIGNED_MALLOC */

#if defined(HAVE_POSIX_MEMALIGN)
#define USE_VECTORS
#define USE_FAKE_MALLOC
#endif

#endif /* HAVE_16_BYTES_ALIGNED_MALLOC */

#endif /* HAVE_GCC_VECTORS && __SSE2__ && HAVE_EMMINTRIN_H && FC_USES_MALLOC */


#endif
