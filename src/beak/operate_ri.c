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

#include <assert.h>

void FC_FUNC_(doperate_ri,DOPERATE_RI)(const int * opn, 
				       const ffloat * restrict w, 
				       const int * opnri,
				       const int * opri,
				       const int * rimap_inv,
				       const int * rimap_inv_max,
				       const ffloat * fi, 
				       ffloat * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j, indexj;
  const int * restrict index;
  register ffloat a0, a1, a2, a3;
  register ffloat a4, a5, a6, a7;
  register ffloat wj;
  assert(MAX_OP_N >= n);

  for (l = 0; l < nri ; l++) {

    index  = opri + n*l;

    i = rimap_inv[l];

    for (; i < rimap_inv_max[l] - 8 + 1; i+=8){
      a0 = a1 = a2 = a3 = 0.0;
      a4 = a5 = a6 = a7 = 0.0;
      for(j = 0; j < n; j++){
	wj = w[j];
	indexj = index[j];
	a0 += wj*(fi + indexj)[i + 0];
	a1 += wj*(fi + indexj)[i + 1];
	a2 += wj*(fi + indexj)[i + 2];
	a3 += wj*(fi + indexj)[i + 3];
	
	a4 += wj*(fi + indexj)[i + 4];
	a5 += wj*(fi + indexj)[i + 5];
	a6 += wj*(fi + indexj)[i + 6];
	a7 += wj*(fi + indexj)[i + 7];
      }
      fo[i  ] = a0;
      fo[i+1] = a1;
      fo[i+2] = a2;
      fo[i+3] = a3;
      
      fo[i+4] = a4;
      fo[i+5] = a5;
      fo[i+6] = a6;
      fo[i+7] = a7;

    }
    
    for (; i < rimap_inv_max[l]; i++){
      a0 = 0.0;
      a1 = 0.0;
      a2 = 0.0;
      j = 0;
      for(; j < n - 5; j += 6){
	a0 += w[j    ]*(fi + index[j    ])[i];
	a1 += w[j + 1]*(fi + index[j + 1])[i];
	a2 += w[j + 2]*(fi + index[j + 2])[i];
	a0 += w[j + 3]*(fi + index[j + 3])[i];
	a1 += w[j + 4]*(fi + index[j + 4])[i];
	a2 += w[j + 5]*(fi + index[j + 5])[i];
      }
      a1 += a2;
      for(; j < n; j++) a0 += w[j]*(fi + index[j])[i];
      fo[i] = a0 + a1;
    }
    
  }

}

void FC_FUNC_(zoperate_ri,ZOPERATE_RI)(const int * opn, 
				       const ffloat * restrict w, 
				       const int * opnri,
				       const int * opri,
				       const int * rimap_inv,
				       const int * rimap_inv_max,
				       const comp * fi, 
				       comp * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j;
  const int * restrict index;
  const int * restrict index1;
  const comp * ffi[MAX_OP_N];
  register ffloat a0, a1, a2, a3;
  register ffloat a4, a5, a6, a7;

  for (l = 0; l < nri ; l++) {

    index = opri + n*l;
    index1= opri + n*(l+1);

    i = rimap_inv[l];

    a0 = 0.0;
    a1 = 0.0;
    for(j = 0; j < n ; j++) {
      ffi[j] = fi + index[j];
      a0 += w[j]*ffi[j][i].re;
      a1 += w[j]*ffi[j][i].im;
      /* prefetch */
#if defined(OCT_ITANIUM)
      __builtin_prefetch(fi + rimap_inv[l+1] + index1[j], 0, 3);
#endif
    }
    fo[i].re = a0;
    fo[i].im = a1;
    i++;

    for (; i < (rimap_inv_max[l] - 4 + 1) ; i+=4){
      a0 = a1 = a2 = a3 = 0.0;
      a4 = a5 = a6 = a7 = 0.0;
      for(j = 0; j < n; j++) {
	a0 += w[j]*ffi[j][i+0].re;
	a1 += w[j]*ffi[j][i+0].im;
	a2 += w[j]*ffi[j][i+1].re;
	a3 += w[j]*ffi[j][i+1].im;
	a4 += w[j]*ffi[j][i+2].re;
	a5 += w[j]*ffi[j][i+2].im;
	a6 += w[j]*ffi[j][i+3].re;
	a7 += w[j]*ffi[j][i+3].im;
      }
      fo[i  ].re = a0;
      fo[i  ].im = a1;
      fo[i+1].re = a2;
      fo[i+1].im = a3;
      fo[i+2].re = a4;
      fo[i+2].im = a5;
      fo[i+3].re = a6;
      fo[i+3].im = a7;

    }

    for (; i < rimap_inv_max[l]; i++){
      
      a0 = 0.0;
      a1 = 0.0;

      for(j = 0; j < n; j++) {
	a0 += w[j]*ffi[j][i].re;
	a1 += w[j]*ffi[j][i].im;
      }

      fo[i].re = a0;
      fo[i].im = a1;
    }

  }

}


void FC_FUNC_(dgauss_seidel,DGAUSS_SEIDEL)(const int * opn, 
				           const ffloat * restrict w, 
				           const int * opnri,
				           const int * opri,
				           const int * rimap_inv,
				           const int * rimap_inv_max,
					   ffloat * restrict factor,
				           ffloat * pot, 
				           ffloat * restrict rho){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j;
  const int * index;
  register ffloat a0;
  const ffloat * ffi[MAX_OP_N];
  register const ffloat fac = factor[0];
  assert(MAX_OP_N >= n);

  for (l = 0; l < nri ; l++) {

    index  = opri + n*l;

    i = rimap_inv[l];

    a0 = 0.0;
    for(j = 0; j < n ; j++) {
      ffi[j] = pot + index[j];
      a0 += w[j]*ffi[j][i];
    }
    pot[i] += fac*(a0 - rho[i]);
    i++;

    for (; i < rimap_inv_max[l]; i++){
      a0 = 0.0;
      for(j = 0; j < n; j++) a0 += w[j] * ffi[j][i];
      pot[i] += fac*(a0 - rho[i]);
    }

  }

}

