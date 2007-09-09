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

#if defined(HAVE_IA64INTRIN_H)
#include <ia64intrin.h>
#endif

void FC_FUNC_(doperate_ri,DOPERATE_RI)(const int * opnp, 
				       const int * opn, 
				       const ffloat * restrict w, 
				       const int * opnri,
				       const int * opri,
				       const int * rimap_inv,
				       const ffloat * fi, 
				       ffloat * restrict fo){

  const int n = opn[0];
  const int np = opnp[0];
  const int nri = opnri[0];

  int l, i, j, nm2;
  const int * restrict index;
  const int * restrict index1;
  register ffloat a0, a1, a2, a3;
  register ffloat a4, a5, a6, a7;
  const ffloat * restrict ffi[30];

  i = 0;
  for (l = 0; l < nri ; l++) {
    index  = opri + n * l;
    index1 = opri + n * (l+1);

    for(j = 0; j < n ; j++) {
      ffi[j] = fi + index[j];
      /* prefecth */
#if defined(OCT_ITANIUM) && defined(HAVE_IA64INTRIN_H)
      __lfetch(3, fi + rimap_inv[l] + index1[j]);
#endif
    }

    for (; i < rimap_inv[l] - 8 + 1; i+=8){
      a0 = a1 = a2 = a3 = 0.0;
      a4 = a5 = a6 = a7 = 0.0;
      for(j = 0; j < n; j++){

	a0 += w[j] * ffi[j][i+0];
	a1 += w[j] * ffi[j][i+1];
	a2 += w[j] * ffi[j][i+2];
	a3 += w[j] * ffi[j][i+3];
	
	a4 += w[j] * ffi[j][i+4];
	a5 += w[j] * ffi[j][i+5];
	a6 += w[j] * ffi[j][i+6];
	a7 += w[j] * ffi[j][i+7];
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
    
    for (; i < rimap_inv[l]; i++){
      a0 = 0.0;
      for(j = 0; j < n; j++) a0 += w[j] * ffi[j][i];
      fo[i] = a0;
    }

  }

}

void FC_FUNC_(zoperate_ri,ZOPERATE_RI)(const int * opnp, 
				       const int * opn, 
				       const ffloat * restrict w, 
				       const int * opnri,
				       const int * opri,
				       const int * rimap_inv,
				       const comp * fi, 
				       comp * restrict fo){

  const int n = opn[0];
  const int np = opnp[0];
  const int nri = opnri[0];

  int l, i, j, nm2;
  const int * restrict index;
  const int * restrict index1;
  const comp * restrict ffi[30];
  comp a0, a1, a2, a3;

  i = 0;
  for (l = 0; l < nri ; l++) {
    index = opri + n * l;
    index1= opri + n * (l+1);

    for(j = 0; j < n ; j++) {
      ffi[j] = fi + index[j];
      /* prefetch */
#if defined(OCT_ITANIUM) && defined(HAVE_IA64INTRIN_H)
      __lfetch(3, fi + rimap_inv[l] + index1[j]);
#endif
    }

    for (; i < (rimap_inv[l] - 4 + 1) ; i+=4){
      a0.re = a1.re = a2.re = a3.re = 0.0;
      a0.im = a1.im = a2.im = a3.im = 0.0;
      
      for(j = 0; j < n; j++) {
	a0.re += w[j] * ffi[j][i+0].re;
	a0.im += w[j] * ffi[j][i+0].im;
	a1.re += w[j] * ffi[j][i+1].re;
	a1.im += w[j] * ffi[j][i+1].im;
	a2.re += w[j] * ffi[j][i+2].re;
	a2.im += w[j] * ffi[j][i+2].im;
	a3.re += w[j] * ffi[j][i+3].re;
	a3.im += w[j] * ffi[j][i+3].im;
      }
      fo[i  ].re = a0.re;
      fo[i  ].im = a0.im;
      fo[i+1].re = a1.re;
      fo[i+1].im = a1.im;
      fo[i+2].re = a2.re;
      fo[i+2].im = a2.im;
      fo[i+3].re = a3.re;
      fo[i+3].im = a3.im;

    }

    for (; i < rimap_inv[l]; i++){
      
      a0.re = 0.0;
      a0.im = 0.0;

      for(j = 0; j < n; j++) {
	a0.re += w[j] * ffi[j][i].re;
	a0.im += w[j] * ffi[j][i].im;
      }

      fo[i].re = a0.re;
      fo[i].im = a0.im;
    }

  }

}
