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
  const int unroll16 = max1(16*VEC_SIZE >> ldf);
  const int unroll8  = max1(8*VEC_SIZE >> ldf);
  const int unroll4  = max1(4*VEC_SIZE >> ldf);
  const int unroll2  = max1(2*VEC_SIZE >> ldf);
  const int unroll1  = max1(1*VEC_SIZE >> ldf);

  int l, i, j;
  const int * restrict index;
   
  for (l = 0; l < nri ; l++) {
    index = opri + n * l;
    i = rimap_inv[l];

#if DEPTH >= 16
    for (; i < (rimap_inv_max[l] - unroll16 + 1) ; i += unroll16){
      int k;
      for(k = 0; k < (1<<ldf); k += 16*VEC_SIZE){
	register VEC_TYPE a0, a1, a2, a3;
	register VEC_TYPE a4, a5, a6, a7;
	register VEC_TYPE a8, a9, aa, ab;
	register VEC_TYPE ac, ad, ae, af;

	a0 = a1 = a2 = a3 = VEC_ZERO;
	a4 = a5 = a6 = a7 = VEC_ZERO;
	a8 = a9 = aa = ab = VEC_ZERO;
	ac = ad = ae = af = VEC_ZERO;

	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  int indexj = (index[j] + i)<<ldf;
	  a0 = VEC_FMA(wj, LOAD(fi + indexj              + k), a0);
	  a1 = VEC_FMA(wj, LOAD(fi + indexj + 1*VEC_SIZE + k), a1);
	  a2 = VEC_FMA(wj, LOAD(fi + indexj + 2*VEC_SIZE + k), a2);
	  a3 = VEC_FMA(wj, LOAD(fi + indexj + 3*VEC_SIZE + k), a3);
	  a4 = VEC_FMA(wj, LOAD(fi + indexj + 4*VEC_SIZE + k), a4);
	  a5 = VEC_FMA(wj, LOAD(fi + indexj + 5*VEC_SIZE + k), a5);
	  a6 = VEC_FMA(wj, LOAD(fi + indexj + 6*VEC_SIZE + k), a6);
	  a7 = VEC_FMA(wj, LOAD(fi + indexj + 7*VEC_SIZE + k), a7);;
	  a8 = VEC_FMA(wj, LOAD(fi + indexj + 8*VEC_SIZE + k), a8);
	  a9 = VEC_FMA(wj, LOAD(fi + indexj + 9*VEC_SIZE + k), a9);
	  aa = VEC_FMA(wj, LOAD(fi + indexj + 10*VEC_SIZE + k), aa);
	  ab = VEC_FMA(wj, LOAD(fi + indexj + 11*VEC_SIZE + k), ab);
	  ac = VEC_FMA(wj, LOAD(fi + indexj + 12*VEC_SIZE + k), ac);
	  ad = VEC_FMA(wj, LOAD(fi + indexj + 13*VEC_SIZE + k), ad);
	  ae = VEC_FMA(wj, LOAD(fi + indexj + 14*VEC_SIZE + k), ae);
	  af = VEC_FMA(wj, LOAD(fi + indexj + 15*VEC_SIZE + k), af);
	}
	STORE(fo + (i<<ldf)              + k, a0);
	STORE(fo + (i<<ldf) + 1*VEC_SIZE + k, a1);
	STORE(fo + (i<<ldf) + 2*VEC_SIZE + k, a2);
	STORE(fo + (i<<ldf) + 3*VEC_SIZE + k, a3);
	STORE(fo + (i<<ldf) + 4*VEC_SIZE + k, a4);
	STORE(fo + (i<<ldf) + 5*VEC_SIZE + k, a5);
	STORE(fo + (i<<ldf) + 6*VEC_SIZE + k, a6);
	STORE(fo + (i<<ldf) + 7*VEC_SIZE + k, a7);
	STORE(fo + (i<<ldf) + 8*VEC_SIZE + k, a8);
	STORE(fo + (i<<ldf) + 9*VEC_SIZE + k, a9);
	STORE(fo + (i<<ldf) + 10*VEC_SIZE + k, aa);
	STORE(fo + (i<<ldf) + 11*VEC_SIZE + k, ab);
	STORE(fo + (i<<ldf) + 12*VEC_SIZE + k, ac);
	STORE(fo + (i<<ldf) + 13*VEC_SIZE + k, ad);
	STORE(fo + (i<<ldf) + 14*VEC_SIZE + k, ae);
	STORE(fo + (i<<ldf) + 15*VEC_SIZE + k, af);
      }
    }
#endif

#if DEPTH >= 8
    for (; i < (rimap_inv_max[l] - unroll8 + 1) ; i += unroll8){
      int k;
      for(k = 0; k < (1<<ldf); k += 8*VEC_SIZE){
	register VEC_TYPE a0, a1, a2, a3;
	register VEC_TYPE a4, a5, a6, a7;

	a0 = a1 = a2 = a3 = VEC_ZERO;
	a4 = a5 = a6 = a7 = VEC_ZERO;

	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  int indexj = (index[j] + i)<<ldf;
	  a0 = VEC_FMA(wj, LOAD(fi + indexj              + k), a0);
	  a1 = VEC_FMA(wj, LOAD(fi + indexj + 1*VEC_SIZE + k), a1);
	  a2 = VEC_FMA(wj, LOAD(fi + indexj + 2*VEC_SIZE + k), a2);
	  a3 = VEC_FMA(wj, LOAD(fi + indexj + 3*VEC_SIZE + k), a3);
	  a4 = VEC_FMA(wj, LOAD(fi + indexj + 4*VEC_SIZE + k), a4);
	  a5 = VEC_FMA(wj, LOAD(fi + indexj + 5*VEC_SIZE + k), a5);
	  a6 = VEC_FMA(wj, LOAD(fi + indexj + 6*VEC_SIZE + k), a6);
	  a7 = VEC_FMA(wj, LOAD(fi + indexj + 7*VEC_SIZE + k), a7);;
	}
	STORE(fo + (i<<ldf)              + k, a0);
	STORE(fo + (i<<ldf) + 1*VEC_SIZE + k, a1);
	STORE(fo + (i<<ldf) + 2*VEC_SIZE + k, a2);
	STORE(fo + (i<<ldf) + 3*VEC_SIZE + k, a3);
	STORE(fo + (i<<ldf) + 4*VEC_SIZE + k, a4);
	STORE(fo + (i<<ldf) + 5*VEC_SIZE + k, a5);
	STORE(fo + (i<<ldf) + 6*VEC_SIZE + k, a6);
	STORE(fo + (i<<ldf) + 7*VEC_SIZE + k, a7);
      }
    }
#endif

#if DEPTH >= 4
    for (; i < (rimap_inv_max[l] - unroll4 + 1) ; i += unroll4){
      int k;
      for(k = 0; k < (1<<ldf); k += 4*VEC_SIZE){
	register VEC_TYPE a0, a1, a2, a3;

	a0 = a1 = a2 = a3 = VEC_ZERO;

	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  int indexj = (index[j] + i)<<ldf;
	  a0 = VEC_FMA(wj, LOAD(fi + indexj              + k), a0);
	  a1 = VEC_FMA(wj, LOAD(fi + indexj + 1*VEC_SIZE + k), a1);
	  a2 = VEC_FMA(wj, LOAD(fi + indexj + 2*VEC_SIZE + k), a2);
	  a3 = VEC_FMA(wj, LOAD(fi + indexj + 3*VEC_SIZE + k), a3);
	}
	STORE(fo + (i<<ldf)              + k, a0);
	STORE(fo + (i<<ldf) + 1*VEC_SIZE + k, a1);
	STORE(fo + (i<<ldf) + 2*VEC_SIZE + k, a2);
	STORE(fo + (i<<ldf) + 3*VEC_SIZE + k, a3);
      }
    }
#endif

#if DEPTH >= 2
    for (; i < (rimap_inv_max[l] - unroll2 + 1) ; i += unroll2){
      int k;
      for(k = 0; k < (1<<ldf); k += 2*VEC_SIZE){
	register VEC_TYPE a0, a1;

	a0 = a1 = VEC_ZERO;

	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  int indexj = (index[j] + i)<<ldf;
	  a0 = VEC_FMA(wj, LOAD(fi + indexj              + k), a0);
	  a1 = VEC_FMA(wj, LOAD(fi + indexj + 1*VEC_SIZE + k), a1);
	}
	STORE(fo + (i<<ldf)              + k, a0);
	STORE(fo + (i<<ldf) + 1*VEC_SIZE + k, a1);
      }
    }
#endif

#if DEPTH >= 1
    for (; i < (rimap_inv_max[l] - unroll1 + 1) ; i += unroll1){
      int k;
      for(k = 0; k < (1<<ldf); k += VEC_SIZE){
	register VEC_TYPE a0 = VEC_ZERO;
	for(j = 0; j < n; j++) {
	  register VEC_TYPE wj = VEC_SCAL(w[j]);
	  int indexj = (index[j] + i)<<ldf;
	  a0 = VEC_FMA(wj, LOAD(fi + indexj + k), a0);
	}
	STORE(fo + (i<<ldf) + k, a0);
      }
    }
#endif

#if VEC_SIZE > 1
    for (; i < rimap_inv_max[l]; i++){
      int k;
      for(k = 0; k < (1<<ldf); k++){
	double a = 0.0;
	for(j = 0; j < n; j++) a += w[j]*fi[((index[j] + i)<<ldf) + k];
	fo[(i<<ldf) + k] = a;
      }
    }
#endif

  } /* l */

}

#undef LOAD
#undef STORE
