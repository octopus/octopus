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

 $Id$
*/

#include <config.h>

#include <assert.h>

/* This code is based on the implementation of the Hilbert curves
   presented in:

   J. Skilling, Programming the Hilbert curve, AIP Conf. Proc. 707, 381 (2004); http://dx.doi.org/10.1063/1.1751381

*/

void InttoTranspose(const int dim, const long long int h, int * x){
  /* the code uses some funny way of storing the bits */
  
  int idir, ibit, ifbit;

  for(idir = 0; idir < dim; idir++) x[idir] = 0;
  
  ifbit = 0;
  for(ibit = 0; ibit < 21; ibit++){
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

void FC_FUNC_(hilbert_index_to_point, HILBERT_INDEX_TO_POINT)(const int * dim, const int * nbits, const long long int * index, int * point){
  
  InttoTranspose(*dim, *index, point);
  TransposetoAxes(point, *nbits, *dim);

}

