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

#include "beak.h"

#if defined(SINGLE_PRECISION)

#if defined(HAVE_XMMINTRIN_H)
# include <xmmintrin.h>
#endif

void FC_FUNC_(doperate_sse,DOPERATE_SSE)(const int * opnp, 
					 const int * opn, 
					 const float * restrict w, 
					 const int * opi, 
					 const float * fi, 
					 float * restrict fo){

#if !defined(USE_VECTORS)
  fprintf(stderr, "Not available: this is a bug.\n");
  exit(1);
#else

  const int n = opn[0];
  const int np = opnp[0];

  int i, j;
  const float * restrict mfi;
  __m128 vw[MAX_OP_N];
  mfi   = fi - 1;

  {

    for(j = 0; j < n ; j++) vw[j] =_mm_set1_ps(w[j]);
    
    for(i = 0; i < (np-4+1); i+=4) {
      const int * restrict index0 = opi + n*i;
      const int * restrict index1 = opi + n*(i+1);
      const int * restrict index2 = opi + n*(i+2);
      const int * restrict index3 = opi + n*(i+3);

      register __m128 a, c;

      a = _mm_mul_ps(vw[0], _mm_setr_ps(mfi[index0[0]], mfi[index1[0]], mfi[index2[0]], mfi[index3[0]]));

      for(j = 1; j < n; j++) {
	c = _mm_setr_ps(mfi[index0[j]], mfi[index1[j]], mfi[index2[j]], mfi[index3[j]]);
	a = _mm_add_ps(a, _mm_mul_ps(vw[j], c));
      }

      _mm_storeu_ps(fo+i, a);

    }
    
    for(; i < np; i++) {
      const int * restrict index; 
      index = opi + n*i;
      fo[i] = 0.0;
      for(j=0; j < n; j++) fo[i] += w[j] * mfi[index[j]];
    }

  }

#endif

}

void FC_FUNC_(zoperate_sse,ZOPERATE_SSE)(const int * opnp, 
					 const int * opn, 
					 const float * restrict w, 
					 const int * opi, 
					 const float * fi, 
					 float * restrict fo){

  fprintf(stderr, "Not available: this is a bug.\n");
  exit(1);

}


#else /* DOUBLE PRECISION */

#if defined(HAVE_EMMINTRIN_H)
# include <emmintrin.h>
#endif

void FC_FUNC_(doperate_sse,DOPERATE_SSE)(const int * opnp, 
					 const int * opn, 
					 const double * restrict w, 
					 const int * opi, 
					 const double * fi, 
					 double * restrict fo){

#if !defined(USE_VECTORS)
  fprintf(stderr, "Not available: this is a bug.\n");
  exit(1);
#else
  const int n = opn[0];
  const int np = opnp[0];

  int i, j, nm2;
  const double * restrict mfi;
  __m128d vw[MAX_OP_N];

  mfi   = fi - 1;

  {
    for(j = 0; j < n ; j++) vw[j] =_mm_set1_pd(w[j]);
    
    for(i = 0; i < (np-2+1); i+=2) {
      const int * restrict index0 = opi + n*i;
      const int * restrict index1 = opi + n*i+n;

      register __m128d a, c;

      a = _mm_mul_pd(vw[0], _mm_setr_pd(mfi[index0[0]], mfi[index1[0]]) );

      for(j = 1; j < n; j++) {
	c = _mm_setr_pd(mfi[index0[j]], mfi[index1[j]]);
	a = _mm_add_pd(a, _mm_mul_pd(vw[j], c));
      }

      _mm_storeu_pd(fo+i  , a);

    }

    for(; i < np; i++) {
      const int * restrict index; 
      index = opi + n*i;
      fo[i] = 0.0;
      for(j=0; j < n; j++) fo[i] += w[j] * mfi[index[j]];
    }

  }
#endif  
}

#if !defined(USE_VECTORS)

void FC_FUNC_(zoperate_sse,ZOPERATE_SSE)(){
  fprintf(stderr, "Not available: this is a bug.\n");
  exit(1);
}

#else

void FC_FUNC_(zoperate_sse,ZOPERATE_SSE)(const int * opnp, 
				const int * opn, 
				const double * restrict w, 
				const int * opi, 
				const __m128d * fi, 
				__m128d * restrict fo){

  const int n  = opn[0];
  const int np = opnp[0];
  const int * restrict index = opi;
  const int nm2  = n - 6 + 1;
  int i, j;
  __m128d vw[MAX_OP_N];
  register __m128d a, b, c, d, e, f;

#pragma omp parallel private(a, b, c, d, e, f, index, i, j, vw)
  {

    for(j = 0; j < n ; j++) vw[j] =_mm_set1_pd(w[j]);

#pragma omp for
    for(i = 0; i < np; i++) {

      index = opi + n*i;

      a = _mm_mul_pd(vw[0], fi[index[0]-1]);
      b = _mm_setzero_pd();
      c = _mm_setzero_pd();
      d = _mm_setzero_pd();
      e = _mm_setzero_pd();
      f = _mm_setzero_pd();

      for(j = 1; j < nm2; j += 6){
        a = _mm_add_pd(a, _mm_mul_pd(vw[j  ], fi[index[j+0]-1]));
        b = _mm_add_pd(b, _mm_mul_pd(vw[j+1], fi[index[j+1]-1]));
	c = _mm_add_pd(c, _mm_mul_pd(vw[j+2], fi[index[j+2]-1]));
        d = _mm_add_pd(d, _mm_mul_pd(vw[j+3], fi[index[j+3]-1]));
        e = _mm_add_pd(e, _mm_mul_pd(vw[j+4], fi[index[j+4]-1]));
        f = _mm_add_pd(f, _mm_mul_pd(vw[j+5], fi[index[j+5]-1]));
      }
    
      for(; j < n; j++) a = _mm_add_pd(a, _mm_mul_pd(vw[j], fi[index[j]-1]));

      a = _mm_add_pd(a, b);
      a = _mm_add_pd(a, d);
      a = _mm_add_pd(a, f);
      a = _mm_add_pd(a, c);
      fo[i] = _mm_add_pd(a, e);

    }

  }

}
#endif

#endif
