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

void FC_FUNC_(doperate_olu_c,DOPERATE_OLU_C)(const int * opnp, 
					 const int * opn, 
					 const ffloat * restrict w, 
					 const int * opi, 
					 const ffloat * fi, 
					 ffloat * restrict fo){

  const int n = opn[0];
  const int np = opnp[0];

  int i, j, nm2;
  const ffloat * restrict mfi = fi - 1;

  for(i = 0; i < (np - 8 + 1); i += 8) {
    const int * restrict index0 = opi + n*i    ; 
    const int * restrict index1 = opi + n*i+n  ;
    const int * restrict index2 = opi + n*i+n*2;
    const int * restrict index3 = opi + n*i+n*3;
    const int * restrict index4 = opi + n*i+n*4; 
    const int * restrict index5 = opi + n*i+n*5; 
    const int * restrict index6 = opi + n*i+n*6;
    const int * restrict index7 = opi + n*i+n*7;

    ffloat register a0 = w[0] * mfi[index0[0]];
    ffloat register a1 = w[0] * mfi[index1[0]];
    ffloat register a2 = w[0] * mfi[index2[0]];
    ffloat register a3 = w[0] * mfi[index3[0]];
    ffloat register a4 = w[0] * mfi[index4[0]];
    ffloat register a5 = w[0] * mfi[index5[0]];
    ffloat register a6 = w[0] * mfi[index6[0]];
    ffloat register a7 = w[0] * mfi[index7[0]];

    for(j=1; j < n; j++) {
      a0 += w[j] * mfi[index0[j]];
      a1 += w[j] * mfi[index1[j]];
      a2 += w[j] * mfi[index2[j]];
      a3 += w[j] * mfi[index3[j]];
      a4 += w[j] * mfi[index4[j]];
      a5 += w[j] * mfi[index5[j]];
      a6 += w[j] * mfi[index6[j]];
      a7 += w[j] * mfi[index7[j]];
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

  for(; i < np; i++) {
    const int * restrict index; 
    index = opi + n*i;
    fo[i] = 0.0;
    for(j=0; j < n; j++) fo[i] += w[j] * mfi[index[j]];
  }

}
