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

#include <stdio.h>

#include <config.h>

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

#define _XOPEN_SOURCE 600
#include <stdlib.h>

#if defined(USE_FAKE_MALLOC)
#include <errno.h>

/* 
Override calls to malloc with calls to posix_memalign. 

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


#if defined(USE_VECTORS)

#include <emmintrin.h>

void FC_FUNC_(zoperate_sse,ZOPERATE_SSE)(const int * opnp, 
				const int * opn, 
				const double * restrict w, 
				const int * opi, 
				const __m128d * fi, 
				__m128d * restrict fo){

  /* GCC 3.x requires these variables to be declared static to have the proper alignment */
#if __GNUC__ <= 3
#define register static
#endif

  register __m128d a  __attribute__ ((__aligned__ (16)));
  register __m128d b  __attribute__ ((__aligned__ (16)));
  register __m128d c  __attribute__ ((__aligned__ (16)));
  register __m128d d  __attribute__ ((__aligned__ (16)));
  register __m128d e  __attribute__ ((__aligned__ (16)));
  register __m128d f  __attribute__ ((__aligned__ (16)));

#undef register

  const int n  = opn[0];
  const int np = opnp[0];
  const int * restrict index = opi;
  const int nm2  = n - UNROLL + 1;
  int i, j;

  __m128d * restrict vw;

#pragma omp parallel private(a, b, c, d, e, f, index, i, j, vw)
  {

    vw = malloc(n*16);

    for(j = 0; j < n ; j++) {
      vw[j] =_mm_set1_pd(w[j]);
    }
    

#pragma omp for
    for(i = 0; i < np; i++) {

#ifdef HAVE_C_OMP
      index = opi + n*i;
#endif

      a = _mm_mul_pd(vw[0], fi[index[0]-1]);
      b = _mm_setzero_pd();
      c = _mm_setzero_pd();
      d = _mm_setzero_pd();
      e = _mm_setzero_pd();
      f = _mm_setzero_pd();

      for(j = 1; j < nm2; j += UNROLL){
        a = _mm_add_pd(a, _mm_mul_pd(vw[j  ], fi[index[j+0]-1]));
        b = _mm_add_pd(b, _mm_mul_pd(vw[j+1], fi[index[j+1]-1]));
	c = _mm_add_pd(c, _mm_mul_pd(vw[j+2], fi[index[j+2]-1]));
        d = _mm_add_pd(d, _mm_mul_pd(vw[j+3], fi[index[j+3]-1]));
        e = _mm_add_pd(e, _mm_mul_pd(vw[j+4], fi[index[j+4]-1]));
        f = _mm_add_pd(f, _mm_mul_pd(vw[j+5], fi[index[j+5]-1]));
      }
    
      for(; j < n; j++) a = _mm_add_pd(a, _mm_mul_pd(vw[j], fi[index[j]-1]));

#ifndef HAVE_C_OMP    
      index += n;
#endif
      a = _mm_add_pd(a, b);
      c = _mm_add_pd(c, d);
      e = _mm_add_pd(e, f);
      a = _mm_add_pd(a, c);
      fo[i] = _mm_add_pd(a, e);

    }

    free(vw);

  }

}

void FC_FUNC_(doperate_sse,DOPERATE_SSE)(const int * opnp, 
					 const int * opn, 
					 const double * restrict w, 
					 const int * opi, 
					 const double * fi, 
					 double * restrict fo){

  const int n = opn[0];
  const int np = opnp[0];

  int i, j, nm2;
  const double * restrict mfi;
  __m128d * restrict vw;

  mfi   = fi - 1;

  {

    vw = malloc(n*16);
    
    for(j = 0; j < n ; j++) {
      vw[j] =_mm_set1_pd(w[j]);
    }
    
    for(i = 0; i < (np-2+1); i+=2) {
      const int * restrict index0 = opi + n*i; 
      const int * restrict index1 = opi + n*i+n; 

      register __m128d a  __attribute__ ((__aligned__ (16)));
      register __m128d c  __attribute__ ((__aligned__ (16)));
      
      c = _mm_setr_pd(mfi[index0[0]], mfi[index1[0]]);
      a = _mm_mul_pd(vw[0], c);
      
      for(j=1; j < n; j++) {
	c = _mm_setr_pd(mfi[index0[j]], mfi[index1[j]]);
	a = _mm_add_pd(a, _mm_mul_pd(vw[j], c));
      }

      _mm_storeu_pd(fo+i, a);

      // This sequence is recommended for the Pentium 4 instead of _mm_storeu_pd:
      //      _mm_store_sd(fo+i, a);
      //      _mm_unpackhi_pd(a, a);
      //      _mm_store_sd(fo+i+1, a);
      
    }

    free(vw);
    
    for(; i < np; i++) {
      const int * restrict index; 
      index = opi + n*i;
      fo[i] = 0.0;
      for(j=0; j < n; j++) fo[i] += w[j] * mfi[index[j]];
    }

  }
  
}


#endif /* USE_VECTORS */
