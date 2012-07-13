/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

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

 $Id$
*/


#include <config.h>

/**
 * Optimizes the order of the FFT grid.
 * The best FFT grid dimensions are given by 2^a*3^b*5^c*7^d*11^e*13^f
 * where a,b,c,d are arbitrary and e,f are 0 or 1.
 * (http://www.fftw.org/doc/Complex-DFTs.html)
 * par is the parity: the result must satisfy n % 2 == par, provided par >= 0.
 */
void fft_optimize(int *n, int par)
{
  if(*n <= 2) return;

  for(;; (*n)++){
    int i, n2;

    if((par >= 0) && (*n % 2 != par)) continue;
    
    /* For debugging:                 */
    /* printf("%i has factors ", *n); */

    n2 = *n;
    for(i = 2; i <= n2; i++){
      if(n2 % i == 0){
        /* For debugging:    */
	/* printf("%i ", i); */
	if(i > 13) break;
	n2 = n2 / i;
	if(i != 11 && i != 13) i--;
      }
    }
    /* For debugging: */
    /* printf("\n");  */
    if(n2 == 1) return;
  }
}

void FC_FUNC_(oct_fft_optimize, OCT_FFT_OPTIMIZE)
  (int *n, int *par)
{
  fft_optimize(n, *par);
}
