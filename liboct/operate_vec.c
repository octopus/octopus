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


#include <config.h>

#if defined(HAVE_GCC_VECTORS) && defined(__SSE2__) && defined(HAVE_EMMINTRIN_H)

#if defined(HAVE_16_BYTES_ALIGNED_MALLOC)

#define USE_VECTORS

#else /* DONT HAVE_16_BYTES_ALIGNED_MALLOC */

#if defined(HAVE_POSIX_MEMALIGN) && defined(HAVE_FAKE_MALLOC)
#define USE_VECTORS
#define USE_FAKE_MALLOC
#endif 

#endif /* HAVE_16_BYTES_ALIGNED_MALLOC */
#endif /* defined(HAVE_GCC_VECTORS) && defined(__SSE2__) && defined(HAVE_EMMINTRIN_H) */


#ifdef USE_FAKE_MALLOC
#define _XOPEN_SOURCE 600
#endif 

#include <stdlib.h>

#ifdef USE_FAKE_MALLOC

#include <errno.h>

/* 

Override calls to malloc with calls to posix_memaling. 

Not the most elegant thing in the world, but it appears that most x86
Fortran compilers can't be instructed to align allocated memory to a
16 bytes boundary as required by SSE.

*/

void * malloc (size_t size){
  int err;
  void * ptr;
  err = posix_memalign((void *) &ptr, 16, size);
  if( err == 0 ) return ptr;
  else {
    errno = ENOMEM;
    return NULL;
  }
}

#endif /* USE_FAKE_MALLOC */



#ifdef USE_VECTORS

typedef double v2df __attribute__ ((vector_size (16)));

#include <emmintrin.h>

#define set_vec(vec, a, b) \
(vec) = _mm_loadl_pd((vec), (a)); \
(vec) = _mm_loadh_pd((vec), (b));

void FC_FUNC(zoperate,ZOPERATE)(const int * opnp, 
				const int * opn, 
				const double * restrict w, 
				const int * opi, 
				const v2df * fi, 
				v2df * restrict fo){

  register v2df a  __attribute__ ((__aligned__ (16)));
  register v2df b  __attribute__ ((__aligned__ (16)));
  register v2df c  __attribute__ ((__aligned__ (16)));
  register v2df d  __attribute__ ((__aligned__ (16)));
  register v2df e  __attribute__ ((__aligned__ (16)));
  register v2df f  __attribute__ ((__aligned__ (16)));

  const int n  = opn[0];
  const int np = opnp[0];
  const int * restrict index = opi;
  const int nm2  = n - UNROLL + 1;
  int i, j;

  v2df * restrict vw;

  vw = malloc(n*16);

  for(j = 0; j < n ; j++) {
    set_vec(vw[j], w+j, w+j);
  }

  for(i = 0; i < np ; i++) {
    
    a = vw[0] * fi[index[0]-1];
    b = _mm_setzero_pd();
    c = _mm_setzero_pd();
    d = _mm_setzero_pd();
    e = _mm_setzero_pd();
    f = _mm_setzero_pd();

    for(j = 1; j < nm2; j += UNROLL){
      a += vw[j  ] * fi[index[j+0]-1];
      b += vw[j+1] * fi[index[j+1]-1];
      c += vw[j+2] * fi[index[j+2]-1];
      d += vw[j+3] * fi[index[j+3]-1];
      e += vw[j+4] * fi[index[j+4]-1];
      f += vw[j+5] * fi[index[j+5]-1];
    }
    
    for(; j < n; j++) a += vw[j] * fi[index[j]-1];
    
    index += n;

    fo[i] = a + b + c + d + e + f;
    
  }
  
  free(vw);

}

#endif /* USE_VECTORS */
