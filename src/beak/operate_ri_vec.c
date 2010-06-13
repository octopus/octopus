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
#include <stdio.h>
#include <assert.h>

#include "beak.h"
#include "vectors.h"

void FC_FUNC_(doperate_ri_vec, DOPERATE_RI_VEC)(const int * opn, 
						const ffloat * restrict w, 
						const int * opnri,
						const int * opri,
						const int * rimap_inv,
						const int * rimap_inv_max,
						const ffloat * restrict fi, 
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

    for (; i < rimap_inv_max[l] - DEPTH*VEC_SIZE + 1; i += DEPTH*VEC_SIZE){
      register VEC_TYPE a0, a1, a2, a3;
#if DEPTH > 4
      register VEC_TYPE a4, a5, a6, a7;
#endif
      a0 = a1 = a2 = a3 = VEC_ZERO;
#if DEPTH > 4
      a4 = a5 = a6 = a7 = VEC_ZERO;
#endif

      for(j = 0; j < n; j++){
	register VEC_TYPE wj = VEC_SCAL(w[j]);
	indexj = index[j] + i;

	a0 = VEC_FMA(wj, VEC_LDU(fi + indexj             ), a0);
	a1 = VEC_FMA(wj, VEC_LDU(fi + indexj + 1*VEC_SIZE), a1);
	a2 = VEC_FMA(wj, VEC_LDU(fi + indexj + 2*VEC_SIZE), a2);
	a3 = VEC_FMA(wj, VEC_LDU(fi + indexj + 3*VEC_SIZE), a3);
#if DEPTH > 4
	a4 = VEC_FMA(wj, VEC_LDU(fi + indexj + 4*VEC_SIZE), a4);
	a5 = VEC_FMA(wj, VEC_LDU(fi + indexj + 5*VEC_SIZE), a5);
	a6 = VEC_FMA(wj, VEC_LDU(fi + indexj + 6*VEC_SIZE), a6);
	a7 = VEC_FMA(wj, VEC_LDU(fi + indexj + 7*VEC_SIZE), a7);
#endif
      }

      VEC_STU(fo + i             , a0);
      VEC_STU(fo + i + 1*VEC_SIZE, a1);
      VEC_STU(fo + i + 2*VEC_SIZE, a2);
      VEC_STU(fo + i + 3*VEC_SIZE, a3);
#if DEPTH > 4
      VEC_STU(fo + i + 4*VEC_SIZE, a4);
      VEC_STU(fo + i + 5*VEC_SIZE, a5);
      VEC_STU(fo + i + 6*VEC_SIZE, a6);
      VEC_STU(fo + i + 7*VEC_SIZE, a7);
#endif
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
					       const double * restrict fi, 
					       double * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j, aligned;
  const int * restrict index;

  /* check whether we got aligned vectors or not */
  aligned = 1;
  aligned = aligned && (((long long) fi)%16 == 0);
  aligned = aligned && (((long long) fo)%16 == 0);

  if(aligned){
   
    for (l = 0; l < nri ; l++) {
      index = opri + n * l;
      i = rimap_inv[l];

      for (; i < (rimap_inv_max[l] - DEPTH*VEC_SIZE/2 + 1) ; i += DEPTH*VEC_SIZE/2){
	register VEC_TYPE a0, a1, a2, a3;
#if DEPTH > 4
	register VEC_TYPE a4, a5, a6, a7;
#endif

	a0 = a1 = a2 = a3 = VEC_ZERO;
#if DEPTH > 4
      	a4 = a5 = a6 = a7 = VEC_ZERO;
#endif
	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  int indexj = (index[j] + i)*2;
	  a0 = VEC_FMA(wj, VEC_LD(fi + indexj             ), a0);
	  a1 = VEC_FMA(wj, VEC_LD(fi + indexj + 1*VEC_SIZE), a1);
	  a2 = VEC_FMA(wj, VEC_LD(fi + indexj + 2*VEC_SIZE), a2);
	  a3 = VEC_FMA(wj, VEC_LD(fi + indexj + 3*VEC_SIZE), a3);
#if DEPTH > 4
	  a4 = VEC_FMA(wj, VEC_LD(fi + indexj + 4*VEC_SIZE), a4);
	  a5 = VEC_FMA(wj, VEC_LD(fi + indexj + 5*VEC_SIZE), a5);
	  a6 = VEC_FMA(wj, VEC_LD(fi + indexj + 6*VEC_SIZE), a6);
	  a7 = VEC_FMA(wj, VEC_LD(fi + indexj + 7*VEC_SIZE), a7);
#endif
	}
	VEC_ST(fo + i*2    , a0);
	VEC_ST(fo + i*2 + 1*VEC_SIZE, a1);
	VEC_ST(fo + i*2 + 2*VEC_SIZE, a2);
	VEC_ST(fo + i*2 + 3*VEC_SIZE, a3);
#if DEPTH > 4
	VEC_ST(fo + i*2 + 4*VEC_SIZE, a4);
	VEC_ST(fo + i*2 + 5*VEC_SIZE, a5);
	VEC_ST(fo + i*2 + 6*VEC_SIZE, a6);
	VEC_ST(fo + i*2 + 7*VEC_SIZE, a7);
#endif
      }

      for (; i < rimap_inv_max[l]; i++){
	register VEC_TYPE a0 = VEC_ZERO;
#if VEC_SIZE < 2
	register VEC_TYPE a1 = VEC_ZERO;
#endif
	for(j = 0; j < n; j++) {
	  a0 = VEC_FMA(VEC_SCAL(w[j]), VEC_LD(fi + (index[j] + i)*2), a0);
#if VEC_SIZE < 2
	  a1 = VEC_FMA(VEC_SCAL(w[j]), VEC_LD(fi + (index[j] + i)*2 + 1), a1);
#endif
	}
	VEC_ST(fo + i*2, a0);
#if VEC_SIZE < 2
	VEC_ST(fo + i*2 + 1, a1);
#endif
      }

    }

  } else { 
    /* not aligned */

    for (l = 0; l < nri ; l++) {
      index = opri + n * l;
      i = rimap_inv[l];

      for (; i < (rimap_inv_max[l] - DEPTH*VEC_SIZE/2 + 1) ; i += DEPTH*VEC_SIZE/2){
	register VEC_TYPE a0, a1, a2, a3;
#if DEPTH > 4
	register VEC_TYPE a4, a5, a6, a7;
#endif

	a0 = a1 = a2 = a3 = VEC_ZERO;
#if DEPTH > 4
      	a4 = a5 = a6 = a7 = VEC_ZERO;
#endif
	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  int indexj = (index[j] + i)*2;
	  a0 = VEC_FMA(wj, VEC_LDU(fi + indexj             ), a0);
	  a1 = VEC_FMA(wj, VEC_LDU(fi + indexj + 1*VEC_SIZE), a1);
	  a2 = VEC_FMA(wj, VEC_LDU(fi + indexj + 2*VEC_SIZE), a2);
	  a3 = VEC_FMA(wj, VEC_LDU(fi + indexj + 3*VEC_SIZE), a3);
#if DEPTH > 4
	  a4 = VEC_FMA(wj, VEC_LDU(fi + indexj + 4*VEC_SIZE), a4);
	  a5 = VEC_FMA(wj, VEC_LDU(fi + indexj + 5*VEC_SIZE), a5);
	  a6 = VEC_FMA(wj, VEC_LDU(fi + indexj + 6*VEC_SIZE), a6);
	  a7 = VEC_FMA(wj, VEC_LDU(fi + indexj + 7*VEC_SIZE), a7);
#endif
	}
	VEC_STU(fo + i*2    , a0);
	VEC_STU(fo + i*2 + 1*VEC_SIZE, a1);
	VEC_STU(fo + i*2 + 2*VEC_SIZE, a2);
	VEC_STU(fo + i*2 + 3*VEC_SIZE, a3);
#if DEPTH > 4
	VEC_STU(fo + i*2 + 4*VEC_SIZE, a4);
	VEC_STU(fo + i*2 + 5*VEC_SIZE, a5);
	VEC_STU(fo + i*2 + 6*VEC_SIZE, a6);
	VEC_STU(fo + i*2 + 7*VEC_SIZE, a7);
#endif
      }

      for (; i < rimap_inv_max[l]; i++){
	register VEC_TYPE a0 = VEC_ZERO;
#if VEC_SIZE < 2
	register VEC_TYPE a1 = VEC_ZERO;
#endif
	for(j = 0; j < n; j++) {
	  a0 = VEC_FMA(VEC_SCAL(w[j]), VEC_LDU(fi + (index[j] + i)*2), a0);
#if VEC_SIZE < 2
	  a1 = VEC_FMA(VEC_SCAL(w[j]), VEC_LDU(fi + (index[j] + i)*2 + 1), a1);
#endif
	}
	VEC_STU(fo + i*2, a0);
#if VEC_SIZE < 2
	VEC_STU(fo + i*2 + VEC_SIZE, a1);
#endif
      }

    }

  }

}

