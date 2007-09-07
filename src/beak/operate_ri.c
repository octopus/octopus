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
  const ffloat * restrict mfi;
  register ffloat a;
  const ffloat * restrict ffi[30];

  i = 0;
  for (l = 0; l < nri ; l++) {
    index = opri + n * l;
    for(j = 0; j < n ; j++) ffi[j] = fi + index[j];
    for (; i < rimap_inv[l]; i++){
      a = 0.0;
      for(j = 0; j < n; j++) a += w[j] * ffi[j][i];
      fo[i] = a;
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
  const comp * restrict mfi;
  comp a;

  mfi   = fi - 1;

  i = 0;
  for (l = 0; l < nri ; l++) {
    index = opri + n * l;
    for (; i < rimap_inv[l]; i++){
      mfi = fi + i;
      
      a.re = 0.0;
      a.im = 0.0;
      
      for(j = 0; j < n; j++) {
	a.re += w[j] * mfi[index[j]].re;
	a.im += w[j] * mfi[index[j]].im;
      }
      fo[i].re = a.re;
      fo[i].im = a.im;
    }

  }

}
