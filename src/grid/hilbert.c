/*
 Copyright (C) 2013 X. Andrade

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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#include <assert.h>

/* This code is based on the implementation of the Hilbert curves
   presented in:

   J. Skilling, Programming the Hilbert curve, AIP Conf. Proc. 707, 381 (2004); http://dx.doi.org/10.1063/1.1751381

*/

void Int4toTranspose(const int dim, const int h, int bits, int * x){
  /* the code uses some funny way of storing the bits */

  int idir, ibit, ifbit;

  for(idir = 0; idir < dim; idir++) x[idir] = 0;

  ifbit = 0;
  for(ibit = 0; ibit < bits; ibit++){
    for(idir = dim - 1; idir >= 0; idir--){
      x[idir] += (((h>>ifbit)&1)<<ibit);
      ifbit++;
    }
  }

}

void InttoTranspose(const int dim, const long long int h, int bits, int * x){
  /* the code uses some funny way of storing the bits */

  int idir, ibit, ifbit;

  for(idir = 0; idir < dim; idir++) x[idir] = 0;

  ifbit = 0;
  for(ibit = 0; ibit < bits; ibit++){
    for(idir = dim - 1; idir >= 0; idir--){
      x[idir] += (((h>>ifbit)&1)<<ibit);
      ifbit++;
    }
  }

}


void TransposetoAxes(int* X, int b, int n ){ /* position, #bits, dimension */

  int N = 2 << (b-1), P, Q, t;
  int i;

  /* Gray decode by H ^ (H/2) */
  t = X[n-1] >> 1;
  for(i = n - 1; i > 0; i--) X[i] ^= X[i-1];
  X[0] ^= t;

  /* Undo excess work */
  for( Q = 2; Q != N; Q <<= 1 ) {
    P = Q - 1;
    for( i = n-1; i >= 0 ; i-- ){
      if( X[i] & Q ) {
	X[0] ^= P; /* invert */
      } else{
	t = (X[0]^X[i]) & P; X[0] ^= t; X[i] ^= t; /* exchange */
      }
    }
  }

}

void TransposetoInt(const int dim, long long int * h, int bits, int * x){
  /* the code uses some funny way of storing the bits */

  int idir, ibit, ifbit;

  *h = 0;

  ifbit = 0;
  for(ibit = 0; ibit < bits; ibit++){
    for(idir = dim - 1; idir >= 0; idir--){
      *h += (((x[idir]>>ibit)&1)<<ifbit);
      ifbit++;
    }
  }
}

void TransposetoInt4(const int dim, int * h, int bits, int * x){
  /* the code uses some funny way of storing the bits */

  int idir, ibit, ifbit;

  *h = 0;

  ifbit = 0;
  for(ibit = 0; ibit < bits; ibit++){
    for(idir = dim - 1; idir >= 0; idir--){
      *h += (((x[idir]>>ibit)&1)<<ifbit);
      ifbit++;
    }
  }
}

void AxestoTranspose(int* X, int b, int n){ /* position, #bits, dimension */
  int M = 1 << (b-1), P, Q, t;
  int i;
  /* Inverse undo */
  for( Q = M; Q > 1; Q >>= 1 ) {
    P = Q - 1;
    for( i = 0; i < n; i++ ) {
      if( X[i] & Q ) X[0] ^= P; /* invert */
      else { /* exchange */
        t = (X[0]^X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
    }
  }

  /* Gray encode */
  for( i = 1; i < n; i++ ) X[i] ^= X[i-1];
  t = 0;
  for( Q = M; Q > 1; Q >>= 1 ) if( X[n-1] & Q ) t ^= Q-1;
  for( i = 0; i < n; i++ ) X[i] ^= t;
}

/* 64 bit integer versions */
void FC_FUNC_(hilbert_index_to_point, HILBERT_INDEX_TO_POINT)(const int * dim, const int * nbits, const long long int * index, int * point){
  InttoTranspose(*dim, *index, *nbits, point);
  TransposetoAxes(point, *nbits, *dim);
}

void FC_FUNC_(hilbert_point_to_index, HILBERT_POINT_TO_INDEX)(const int * dim, const int * nbits, long long int * index, const int * point){
  /* copy point to not overwrite original data */
  int point_copy[*dim];
  for (int i = 0; i < *dim; i++) point_copy[i] = point[i];
  AxestoTranspose(point_copy, *nbits, *dim);
  TransposetoInt(*dim, index, *nbits, point_copy);
}

/* 32 bit integer versions */
void FC_FUNC_(hilbert_index_to_point_int, HILBERT_INDEX_TO_POINT_INT)(const int * dim, const int * nbits, const int * index, int * point){
  Int4toTranspose(*dim, *index, *nbits, point);
  TransposetoAxes(point, *nbits, *dim);
}

void FC_FUNC_(hilbert_point_to_index_int, HILBERT_POINT_TO_INDEX_INT)(const int * dim, const int * nbits, int * index, const int * point){
  /* copy point to not overwrite original data */
  int point_copy[*dim];
  for (int i = 0; i < *dim; i++) point_copy[i] = point[i];
  AxestoTranspose(point_copy, *nbits, *dim);
  TransposetoInt4(*dim, index, *nbits, point_copy);
}
