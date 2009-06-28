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

#ifndef HAVE_BLUE_GENE
#error Internal error, compiling blue gene code without compiler support
#endif

#include <assert.h>

void FC_FUNC_(zoperate_bg,ZOPERATE_BG)(const int * opn, 
				       const ffloat * restrict w, 
				       const int * opnri,
				       const int * opri,
				       const int * rimap_inv,
				       const int * rimap_inv_max,
				       const double * restrict fi, 
				       double * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j, indexj;
  const int * restrict index;
  register double _Complex a0, a1, a2, a3, a4, a5, a6, a7;

  __alignx(16, fi);
  __alignx(16, fo);
  __alignx(16, w);

  for (l = 0; l < nri ; l++) {

    index = opri + n*l;

    i = rimap_inv[l];
    for (; i < (rimap_inv_max[l] - 8 + 1) ; i += 8){

      a0 = a1 = a2 = a3 = a4 = a5 = a6 = a7 = __cmplx(0.0, 0.0);
      for(j = 0; j < n; j++) {
	indexj = 2*(index[j] + i);
	a0 = __fxcpmadd(a0, __lfpd(fi + indexj +  0), w[j]);
	a1 = __fxcpmadd(a1, __lfpd(fi + indexj +  2), w[j]);
	a2 = __fxcpmadd(a2, __lfpd(fi + indexj +  4), w[j]);
	a3 = __fxcpmadd(a3, __lfpd(fi + indexj +  6), w[j]);
	a4 = __fxcpmadd(a4, __lfpd(fi + indexj +  8), w[j]);
	a5 = __fxcpmadd(a5, __lfpd(fi + indexj + 10), w[j]);
	a6 = __fxcpmadd(a6, __lfpd(fi + indexj + 12), w[j]);
	a7 = __fxcpmadd(a7, __lfpd(fi + indexj + 14), w[j]);
      }
      __stfpd(fo + 2*i +  0, a0);
      __stfpd(fo + 2*i +  2, a1);
      __stfpd(fo + 2*i +  4, a2);
      __stfpd(fo + 2*i +  6, a3);
      __stfpd(fo + 2*i +  8, a4);
      __stfpd(fo + 2*i + 10, a5);
      __stfpd(fo + 2*i + 12, a6);
      __stfpd(fo + 2*i + 14, a7);
    }

    for (; i < rimap_inv_max[l]; i++){
      a0 = __cmplx(0.0, 0.0);
      for(j = 0; j < n; j++) {
	a0 = __fxcpmadd(a0, __lfpd(fi + 2*index[j] + 2*i), w[j]);
      }
      __stfpd(fo + 2*i, a0);
    }

  }

}
