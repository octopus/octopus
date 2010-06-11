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

#ifndef HAVE_VEC
#error Internal error, compiling vector code without vector support
#endif

#include <assert.h>

#define VEC_SSE2
#include "vectors.h"

void FC_FUNC_(doperate_ri_vec, DOPERATE_RI_VEC)(const int * opn, 
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

  for (l = 0; l < nri ; l++) {
    register ffloat a;

    index  = opri + n * l;
    
    i = rimap_inv[l];

    for (; i < rimap_inv_max[l] - 4*VEC_SIZE + 1; i+=4*VEC_SIZE){

      register VEC_TYPE a0, a1, a2, a3;
      a0 = a1 = a2 = a3 = VEC_ZERO;

      for(j = 0; j < n; j++){
	register VEC_TYPE wj = VEC_SCAL(w[j]);
	indexj = index[j] + i;

	a0 = VEC_FMA(wj, VEC_LDU(fi + indexj             ), a0);
	a1 = VEC_FMA(wj, VEC_LDU(fi + indexj + 1*VEC_SIZE), a1);
	a2 = VEC_FMA(wj, VEC_LDU(fi + indexj + 2*VEC_SIZE), a2);
	a3 = VEC_FMA(wj, VEC_LDU(fi + indexj + 3*VEC_SIZE), a3);
      }

      VEC_STU(fo + i             , a0);
      VEC_STU(fo + i + 1*VEC_SIZE, a1);
      VEC_STU(fo + i + 2*VEC_SIZE, a2);
      VEC_STU(fo + i + 3*VEC_SIZE, a3);

    }

    for (; i < rimap_inv_max[l]; i++){
      a = 0.0;
      for(j = 0; j < n; j++) a += w[j]*(fi + index[j])[i];
      fo[i] = a;
    }

  }

}

void FC_FUNC_(zoperate_ri_vec,ZOPERATE_RI_VEC)(const int * opn, 
					       const ffloat * restrict w, 
					       const int * opnri,
					       const int * opri,
					       const int * rimap_inv,
					       const int * rimap_inv_max,
					       const double * fi, 
					       double * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j, aligned;
  const int * restrict index;
  const double * ffi[MAX_OP_N];

  /* check whether we got aligned vectors or not */
  aligned = 1;
  aligned = aligned && (((long long) fi)%16 == 0);
  aligned = aligned && (((long long) fo)%16 == 0);

  if(aligned){
   
    for (l = 0; l < nri ; l++) {
      register VEC_TYPE a0, a1, a2, a3;;

      index = opri + n * l;
      i = rimap_inv[l];

      a0 = VEC_ZERO;
      for(j = 0; j < n; j++) {
	ffi[j] = fi + index[j]*2;
	a0 = VEC_FMA(VEC_SCAL(w[j]), VEC_LD(ffi[j] + i*2), a0);
      }
      VEC_ST(fo + i*2, a0);
      i++;

      for (; i < (rimap_inv_max[l] - 4 + 1) ; i+=4){

	a0 = a1 = a2 = a3 = _mm_setzero_pd();
      
	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  a0 = VEC_FMA(wj, VEC_LD(ffi[j] + i*2    ), a0);
	  a1 = VEC_FMA(wj, VEC_LD(ffi[j] + i*2 + 2), a1);
	  a2 = VEC_FMA(wj, VEC_LD(ffi[j] + i*2 + 4), a2);
	  a3 = VEC_FMA(wj, VEC_LD(ffi[j] + i*2 + 6), a3);
	}
	VEC_ST(fo + i*2    , a0);
	VEC_ST(fo + i*2 + 2, a1);
	VEC_ST(fo + i*2 + 4, a2);
	VEC_ST(fo + i*2 + 6, a3);
      }

      for (; i < rimap_inv_max[l]; i++){
      
	a0 = VEC_ZERO;
	for(j = 0; j < n; j++) {
	  a0 = VEC_FMA(VEC_SCAL(w[j]), VEC_LD(ffi[j] + i*2), a0);
	}
	VEC_ST(fo + i*2, a0);
      }

    }

  } else { 
    /* not aligned */
    
    for (l = 0; l < nri ; l++) {
      register VEC_TYPE a0, a1, a2, a3;;

      index = opri + n * l;
      i = rimap_inv[l];

      a0 = VEC_ZERO;
      for(j = 0; j < n; j++) {
	ffi[j] = fi + index[j]*2;
	a0 = VEC_FMA(VEC_SCAL(w[j]), VEC_LDU(ffi[j] + i*2), a0);
      }
      VEC_STU(fo + i*2, a0);
      i++;

      for (; i < (rimap_inv_max[l] - 4 + 1) ; i+=4){

	a0 = a1 = a2 = a3 = _mm_setzero_pd();
      
	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  a0 = VEC_FMA(wj, VEC_LDU(ffi[j] + i*2    ), a0);
	  a1 = VEC_FMA(wj, VEC_LDU(ffi[j] + i*2 + 2), a1);
	  a2 = VEC_FMA(wj, VEC_LDU(ffi[j] + i*2 + 4), a2);
	  a3 = VEC_FMA(wj, VEC_LDU(ffi[j] + i*2 + 6), a3);
	}
	VEC_STU(fo + i*2    , a0);
	VEC_STU(fo + i*2 + 2, a1);
	VEC_STU(fo + i*2 + 4, a2);
	VEC_STU(fo + i*2 + 6, a3);
      }

      for (; i < rimap_inv_max[l]; i++){
      
	a0 = VEC_ZERO;
	for(j = 0; j < n; j++) {
	  a0 = VEC_FMA(VEC_SCAL(w[j]), VEC_LDU(ffi[j] + i*2), a0);
	}
	VEC_STU(fo + i*2, a0);
      }

    }
  }

}

