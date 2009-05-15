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
#include "beak.h"
#include <stdio.h>

#if defined(USE_VECTORS)
#include <emmintrin.h>

#if defined(SINGLE_PRECISION)
#define VECSIZE 4
#define ADD   _mm_add_ps
#define MUL   _mm_mul_ps
#define LOAD  _mm_loadu_ps
#define STORE _mm_storeu_ps
#else 
#define VECSIZE 2
#define ADD   _mm_add_pd
#define MUL   _mm_mul_pd
#define LOAD  _mm_loadu_pd
#define STORE _mm_storeu_pd
#endif

#include <assert.h>

void FC_FUNC_(doperate_ri_vec,DOPERATE_RI_VEC)(const int * opn, 
					       const ffloat * restrict w, 
					       const int * opnri,
					       const int * opri,
					       const int * rimap_inv,
					       const int * rimap_inv_max,
					       const ffloat * fi, 
					       ffloat * restrict fo){
  
  const int n = opn[0];
  const int nri = opnri[0];
  
  int l, i, j;
  const int * restrict index;
  int indexj;

#ifdef SINGLE_PRECISION
  __m128 vw[MAX_OP_N] __attribute__((aligned(16)));
  for(j = 0; j < n ; j++) vw[j] =_mm_set1_ps(w[j]);
#else
  __m128d vw[MAX_OP_N] __attribute__((aligned(16)));
  for(j = 0; j < n ; j++) vw[j] =_mm_set1_pd(w[j]);
#endif

  assert(((long long) vw)%16 == 0);

  for (l = 0; l < nri ; l++) {
    register ffloat a;

    index  = opri + n * l;
    
    i = rimap_inv[l];

    for (; i < rimap_inv_max[l] - 4*VECSIZE + 1; i+=4*VECSIZE){

#ifdef SINGLE_PRECISION
      register __m128 a0, a1, a2, a3;
      a0 = a1 = a2 = a3 = _mm_setzero_ps();
#else
      register __m128d a0, a1, a2, a3;
      a0 = a1 = a2 = a3 = _mm_setzero_pd();
#endif

      for(j = 0; j < n; j++){
	indexj = index[j] + i;
	a0 = ADD (a0, MUL (vw[j], LOAD (fi + indexj            ) ));
	a1 = ADD (a1, MUL (vw[j], LOAD (fi + indexj + 1*VECSIZE) ));
	a2 = ADD (a2, MUL (vw[j], LOAD (fi + indexj + 2*VECSIZE) ));
	a3 = ADD (a3, MUL (vw[j], LOAD (fi + indexj + 3*VECSIZE) ));
      }

      STORE (fo + i            , a0);
      STORE (fo + i + 1*VECSIZE, a1);
      STORE (fo + i + 2*VECSIZE, a2);
      STORE (fo + i + 3*VECSIZE, a3);

    }

    for (; i < rimap_inv_max[l]; i++){
      a = 0.0;
      for(j = 0; j < n; j++) a += w[j]*(fi + index[j])[i];
      fo[i] = a;
    }

  }

}

#ifdef SINGLE_PRECISION

void FC_FUNC_(zoperate_ri_vec,ZOPERATE_RI_VEC)(const int * opn, 
					       const ffloat * restrict w, 
					       const int * opnri,
					       const int * opri,
					       const int * rimap_inv,
					       const int * rimap_inv_max,
					       const __m64 * fi, 
					       __m64 * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j;
  const int * restrict index;
  const __m64 * ffi[MAX_OP_N];
  __m128 vw[MAX_OP_N] __attribute__((aligned(16)));

  for(j = 0; j < n ; j++) vw[j] =_mm_set1_ps(w[j]);

  for (l = 0; l < nri ; l++) {

    index = opri + n * l;

    for(j = 0; j < n ; j++) ffi[j] = fi + index[j];

    for (i = rimap_inv[l]; i < (rimap_inv_max[l] - 8 + 1) ; i+=8){
      register __m128 a0, a1, a2, a3;
      a0 = a1 = a2 = a3 = _mm_setzero_ps();
      
      for(j = 0; j < n; j++) {
	a0 = _mm_add_ps(a0, _mm_mul_ps(vw[j], _mm_loadu_ps((float *) (ffi[j]+i+0) )));
	a1 = _mm_add_ps(a1, _mm_mul_ps(vw[j], _mm_loadu_ps((float *) (ffi[j]+i+2) )));
	a2 = _mm_add_ps(a2, _mm_mul_ps(vw[j], _mm_loadu_ps((float *) (ffi[j]+i+4) )));
	a3 = _mm_add_ps(a3, _mm_mul_ps(vw[j], _mm_loadu_ps((float *) (ffi[j]+i+6) )));
      }

      _mm_storeu_ps((float *) (fo+i  ), a0);
      _mm_storeu_ps((float *) (fo+i+2), a1);
      _mm_storeu_ps((float *) (fo+i+4), a2);
      _mm_storeu_ps((float *) (fo+i+6), a3);

    }

    for (; i < rimap_inv_max[l]; i++){
      register __m128 a0, a1;
      
      a0 = _mm_setzero_ps();
      for(j = 0; j < n; j++) {
	a1 = _mm_setzero_ps();
	a1 = _mm_loadl_pi(a1, ffi[j] + i);
	a0 = _mm_add_ps(a0, _mm_mul_ps(vw[j], a1));
      }
      _mm_storel_pi(fo+i, a0);
    }

  }

}

#else /* DOUBLE PRECISION */

void FC_FUNC_(zoperate_ri_vec,ZOPERATE_RI_VEC)(const int * opn, 
					       const ffloat * restrict w, 
					       const int * opnri,
					       const int * opri,
					       const int * rimap_inv,
					       const int * rimap_inv_max,
					       const __m128d * fi, 
					       __m128d * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j, aligned;
  const int * restrict index;
  const __m128d * ffi[MAX_OP_N];

  __m128d vw[MAX_OP_N] __attribute__((aligned(16)));

  for(j = 0; j < n ; j++) vw[j] =_mm_set1_pd(w[j]);

  /* check whether we got aligned vectors or not */
  aligned = 1;
  aligned = aligned && (((long long) fi)%16 == 0);
  aligned = aligned && (((long long) fo)%16 == 0);

  if(aligned){

    for (l = 0; l < nri ; l++) {

      index = opri + n * l;

      register __m128d a0, a1, a2, a3;
      i = rimap_inv[l];

      a0 = _mm_setzero_pd();
      for(j = 0; j < n; j++){
	ffi[j] = fi + index[j];
	a0 = _mm_add_pd(a0, _mm_mul_pd(vw[j], ffi[j][i]));
      }
      fo[i++] = a0;

      for (; i < (rimap_inv_max[l] - 4 + 1) ; i+=4){

	a0 = a1 = a2 = a3 = _mm_setzero_pd();
      
	for(j = 0; j < n; j++) {
	  a0 = _mm_add_pd(a0, _mm_mul_pd(vw[j], ffi[j][i+0]));
	  a1 = _mm_add_pd(a1, _mm_mul_pd(vw[j], ffi[j][i+1]));
	  a2 = _mm_add_pd(a2, _mm_mul_pd(vw[j], ffi[j][i+2]));
	  a3 = _mm_add_pd(a3, _mm_mul_pd(vw[j], ffi[j][i+3]));
	}
	fo[i  ] = a0;
	fo[i+1] = a1;
	fo[i+2] = a2;
	fo[i+3] = a3;
      }

      for (; i < rimap_inv_max[l]; i++){
	a0 = _mm_setzero_pd();
	for(j = 0; j < n; j++) {
	  a0 = _mm_add_pd(a0, _mm_mul_pd(vw[j],ffi[j][i]));
	}
	fo[i] = a0;
      }

    }

  } else { 
    /* not aligned */
    
    for (l = 0; l < nri ; l++) {

      index = opri + n * l;
      register __m128d a0, a1, a2, a3;;
      i = rimap_inv[l];

      a0 = _mm_setzero_pd();      
      for(j = 0; j < n; j++) {
	ffi[j] = fi + index[j];
	a0 = _mm_add_pd(a0, _mm_mul_pd(vw[j], _mm_loadu_pd((double *)(ffi[j] + i))));
      }
      _mm_storeu_pd(fo + i, a0);
      i++;

      for (; i < (rimap_inv_max[l] - 4 + 1) ; i+=4){

	a0 = a1 = a2 = a3 = _mm_setzero_pd();
      
	for(j = 0; j < n; j++) {
	  a0 = _mm_add_pd(a0, _mm_mul_pd(vw[j], _mm_loadu_pd((double *) (ffi[j] + i + 0))));
	  a1 = _mm_add_pd(a1, _mm_mul_pd(vw[j], _mm_loadu_pd((double *) (ffi[j] + i + 1))));
	  a2 = _mm_add_pd(a2, _mm_mul_pd(vw[j], _mm_loadu_pd((double *) (ffi[j] + i + 2))));
	  a3 = _mm_add_pd(a3, _mm_mul_pd(vw[j], _mm_loadu_pd((double *) (ffi[j] + i + 3))));
	}
	_mm_storeu_pd((double *) (fo + i    ), a0);
	_mm_storeu_pd((double *) (fo + i + 1), a1);
	_mm_storeu_pd((double *) (fo + i + 2), a2);
	_mm_storeu_pd((double *) (fo + i + 3), a3);
      }

      for (; i < rimap_inv_max[l]; i++){
      
	a0 = _mm_setzero_pd();
	for(j = 0; j < n; j++) {
	  a0 = _mm_add_pd(a0, _mm_mul_pd(vw[j], _mm_loadu_pd(ffi[j] + i)));
	}
	_mm_storeu_pd(fo + i, a0);
      }

    }
  }

}

#endif

#else

#include <stdlib.h>

void FC_FUNC_(doperate_ri_vec,DOPERATE_RI_VEC)(){
  fprintf(stderr, "Not available: this is a bug.\n");
  exit(1);
}

void FC_FUNC_(zoperate_ri_vec,ZOPERATE_RI_VEC)(){
  fprintf(stderr, "Not available: this is a bug.\n");
  exit(1);
}

#endif
