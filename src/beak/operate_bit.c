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

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

void FC_FUNC_(generate_ribit,GENERATE_RIBIT)
     (const int *nri, const int * opn, const int *opri, int * opri2){
  
  int ii, jj, kk, pos, ibit, bound;
  int * bm;
  pos = 0;
  
  for(ii = 0; ii < nri[0]; ii++){

    // generate the bitmap
    ibit = 1;
    bm = opri2 + pos++;
    *bm = 0;
    bound = 32;

    for(jj = 0; jj < opn[0]; jj++){
      if(jj == bound){
	bm = opri2 + pos++;
	*bm = 0;
	bound += 32;
	ibit = 1;
      }
    
      if (ii == 0) *bm |= ibit;
      else if (opri[opn[0]*ii + jj] != opri[opn[0]*(ii - 1) + jj] ) *bm |= ibit;
      if(*bm & ibit) opri2[pos++] = opri[opn[0]*ii + jj];

      ibit <<= 1;
    }
  
  }
  /*
  jj = open("coeff",  O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH );
  write(jj, opri, nri[0]*opn[0]*4);
  close(jj);

  jj = open("ccoeff",  O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH );
  write(jj, opri2, (pos + 1)*4);
  close(jj);

  printf("size %f %f\n", (pos + 1.0)*4.0/1024.0, nri[0]*opn[0]*4.0/1024.0);
  */
}

void FC_FUNC_(doperate_bit,DOPERATE_BIT)(const int * opn, 
					 const ffloat * restrict w, 
					 const int * opnri,
					 const int * opri2,
					 const int * rimap_inv,
					 const int * rimap_inv_max,
					 const ffloat * fi, 
					 ffloat * restrict fo){
  
  const int n = opn[0];
  const int nri = opnri[0];
  int l, i, j, pos, bm, bound;
  register ffloat a0, a1, a2, a3;
  register ffloat a4, a5, a6, a7;
  const ffloat * ffi[MAX_OP_N];

  assert(MAX_OP_N >= n);

  pos = 0;

  for (l = 0; l < nri ; l++) {

    i = rimap_inv[l];

    bm = opri2[pos++];
    bound = 32;
    for(j = 0; j < n; j++){
      if(j == bound) {
	bm = opri2[pos++];
	bound += 32;
      }
      if(bm & 1) ffi[j] = fi + opri2[pos++];
      bm >>= 1;
    }

    for (; i < rimap_inv_max[l] - 8 + 1; i+=8){
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
  
    for (; i < rimap_inv_max[l]; i++){
      a0 = 0.0;
      for(j = 0; j < n; j++) a0 += w[j] * ffi[j][i];
      fo[i] = a0;
    }

  }

}

void FC_FUNC_(zoperate_bit,ZOPERATE_BIT)(const int * opn, 
					 const ffloat * restrict w, 
					 const int * opnri,
					 const int * opri2,
					 const int * rimap_inv,
					 const int * rimap_inv_max,
					 const comp * fi, 
					 comp * restrict fo){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j, pos, bm, bound;
  const comp * ffi[MAX_OP_N];
  register ffloat a0, a1, a2, a3;
  register ffloat a4, a5, a6, a7;

  pos = 0;

  for (l = 0; l < nri ; l++) {

    bm = opri2[pos++];
    bound = 32;

    for(j = 0; j < n; j+=2){
      if(j == bound) {
	bm = opri2[pos++];
	bound += 32;
      }
      if(bm & 1) ffi[j] = fi + opri2[pos++];
      if(bm & 2) ffi[j + 1] = fi + opri2[pos++];
      bm >>= 2;
    }
 
    for (i = rimap_inv[l]; i < (rimap_inv_max[l] - 4 + 1) ; i+=4){
      a0 = a1 = a2 = a3 = 0.0;
      a4 = a5 = a6 = a7 = 0.0;
    
      for(j = 0; j < n; j++) {
	a0 += w[j] * ffi[j][i+0].re;
	a1 += w[j] * ffi[j][i+0].im;
	a2 += w[j] * ffi[j][i+1].re;
	a3 += w[j] * ffi[j][i+1].im;
	a4 += w[j] * ffi[j][i+2].re;
	a5 += w[j] * ffi[j][i+2].im;
	a6 += w[j] * ffi[j][i+3].re;
	a7 += w[j] * ffi[j][i+3].im;
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
	a0 += w[j] * ffi[j][i].re;
	a1 += w[j] * ffi[j][i].im;
      }

      fo[i].re = a0;
      fo[i].im = a1;
    }

  }

}
