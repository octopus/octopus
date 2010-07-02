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

 $Id: operate_ri_vec.c 2146 2006-05-23 17:36:00Z xavier $
*/

#ifdef ALIGNED
#define LOAD VEC_LD
#define STORE VEC_ST
#else
#define LOAD VEC_LDU
#define STORE VEC_STU
#endif

{
  const int n = opn[0];
  const int nri = opnri[0];
  int unroll = DEPTH*VEC_SIZE >> ldf;
  if(unroll == 0) unroll = 1;

  int l, i, j;
  const int * restrict index;
   
  for (l = 0; l < nri ; l++) {
    index = opri + n * l;
    i = rimap_inv[l];

    for (; i < (rimap_inv_max[l] - unroll + 1) ; i += unroll){
      int k;
      for(k = 0; k < (1<<ldf); k += DEPTH*VEC_SIZE){
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
	  int indexj = (index[j] + i)<<ldf;
	  a0 = VEC_FMA(wj, LOAD(fi + indexj              + k), a0);
	  a1 = VEC_FMA(wj, LOAD(fi + indexj + 1*VEC_SIZE + k), a1);
	  a2 = VEC_FMA(wj, LOAD(fi + indexj + 2*VEC_SIZE + k), a2);
	  a3 = VEC_FMA(wj, LOAD(fi + indexj + 3*VEC_SIZE + k), a3);
#if DEPTH > 4
	  a4 = VEC_FMA(wj, LOAD(fi + indexj + 4*VEC_SIZE + k), a4);
	  a5 = VEC_FMA(wj, LOAD(fi + indexj + 5*VEC_SIZE + k), a5);
	  a6 = VEC_FMA(wj, LOAD(fi + indexj + 6*VEC_SIZE + k), a6);
	  a7 = VEC_FMA(wj, LOAD(fi + indexj + 7*VEC_SIZE + k), a7);
#endif
	}
	STORE(fo + (i<<ldf)              + k, a0);
	STORE(fo + (i<<ldf) + 1*VEC_SIZE + k, a1);
	STORE(fo + (i<<ldf) + 2*VEC_SIZE + k, a2);
	STORE(fo + (i<<ldf) + 3*VEC_SIZE + k, a3);
#if DEPTH > 4
	STORE(fo + (i<<ldf) + 4*VEC_SIZE + k, a4);
	STORE(fo + (i<<ldf) + 5*VEC_SIZE + k, a5);
	STORE(fo + (i<<ldf) + 6*VEC_SIZE + k, a6);
	STORE(fo + (i<<ldf) + 7*VEC_SIZE + k, a7);
#endif
      }
    }

    for (; i < rimap_inv_max[l]; i++){
      int k;
      for(k = 0; k < (1<<ldf); k++){
	double a = 0.0;
	  for(j = 0; j < n; j++) a += w[j]*fi[((index[j] + i)<<ldf) + k];
	  fo[(i<<ldf) + k] = a;
      }
    }

  } /* l */

}

#undef LOAD
#undef STORE
