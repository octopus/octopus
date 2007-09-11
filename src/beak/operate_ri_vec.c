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

#if defined(USE_VECTORS) && !defined(SINGLE_PRECISION)
#include <emmintrin.h>

void FC_FUNC_(doperate_ri_vec,DOPERATE_RI_VEC)(const int * opnp, 
					       const int * opn, 
					       const ffloat * restrict w, 
					       const int * opnri,
					       const int * opri,
					       const int * rimap_inv,
					       const ffloat * fi, 
					       ffloat * restrict fo){
  
  const int n = opn[0];
  const int nri = opnri[0];
  
  int l, i, j;
  const int * restrict index;
  __m128d vw[MAX_OP_N];
  register __m128d a0, a1, a2, a3;
  const ffloat * ffi[MAX_OP_N];

  for(j = 0; j < n ; j++) vw[j] =_mm_set1_pd(w[j]);

  i = 0;
  for (l = 0; l < nri ; l++) {
    index  = opri + n * l;

    for(j = 0; j < n ; j++) ffi[j] = fi + index[j];

    for (; i < rimap_inv[l] - 8 + 1; i+=8){

      a0 = a1 = a2 = a3 = _mm_setzero_pd();

      for(j = 0; j < n; j++){
	a0 = _mm_add_pd(a0, _mm_mul_pd(vw[j], _mm_loadu_pd(ffi[j] + i + 0)));
	a1 = _mm_add_pd(a1, _mm_mul_pd(vw[j], _mm_loadu_pd(ffi[j] + i + 2)));
	a2 = _mm_add_pd(a2, _mm_mul_pd(vw[j], _mm_loadu_pd(ffi[j] + i + 4)));
	a3 = _mm_add_pd(a3, _mm_mul_pd(vw[j], _mm_loadu_pd(ffi[j] + i + 6)));
      }

      _mm_storeu_pd(fo + i    , a0);
      _mm_storeu_pd(fo + i + 2, a1);
      _mm_storeu_pd(fo + i + 4, a2);
      _mm_storeu_pd(fo + i + 6, a3);

    }

    for (; i < rimap_inv[l]; i++){
      register double a = 0.0;
      for(j = 0; j < n; j++) a += w[j] * ffi[j][i];
      fo[i] = a;
    }

  }

}

void FC_FUNC_(zoperate_ri_vec,ZOPERATE_RI_VEC)(const int * opnp, 
					       const int * opn, 
					       const ffloat * restrict w, 
					       const int * opnri,
					       const int * opri,
					       const int * rimap_inv,
					       const __m128d * fi, 
					       __m128d * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j;
  const int * restrict index;
  const __m128d * ffi[MAX_OP_N];
  __m128d vw[MAX_OP_N];
  register __m128d a0, a1, a2, a3;

  for(j = 0; j < n ; j++) vw[j] =_mm_set1_pd(w[j]);

  i = 0;
  for (l = 0; l < nri ; l++) {
    index = opri + n * l;

    for(j = 0; j < n ; j++) ffi[j] = fi + index[j];

    for (; i < (rimap_inv[l] - 4 + 1) ; i+=4){
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

    for (; i < rimap_inv[l]; i++){
      
      a0 = _mm_setzero_pd();
      for(j = 0; j < n; j++) {
	a0 = _mm_add_pd(a0, _mm_mul_pd(vw[j],ffi[j][i]));
      }
      fo[i] = a0;
    }

  }

}
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
