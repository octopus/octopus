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

#define CACHELINE 8 //in doubles

inline void doperate_fallback(const int np, const int n,
		       const double * restrict w, 
		       const int * opi, 
		       const double * restrict fi, 
		       double * restrict fo){
  int i,j;
  const int * restrict index;
  
  index = opi;

  __builtin_prefetch (w, 0, 3);
  __builtin_prefetch (w + CACHELINE, 0, 3);
#if CACHELINE < 16
  __builtin_prefetch (w + 2*CACHELINE, 0, 3);
#endif

  for(i = 0; i < np ; i++) {
    
    fo[i] = 0.0;
    for(j = 0; j < n; j++) fo[i] += w[j] * fi[index[j]-1];
    index += n;
    
  }

}

void FC_FUNC_(doperate,DOPERATE)(const int * opnp, 
	       const int * opn, 
	       const double * restrict w, 
	       const int * opi, 
	       const double * restrict fi, 
	       double * restrict fo){

  const int n = opn[0];
  const int np = opnp[0];

  int i, j;
  const int * restrict index;
  register double a, b, c;

  if( __builtin_expect((n-1)%3,0) ){
    doperate_fallback(np, n, w, opi, fi, fo);
    return;
  }

  index = opi;

  /* this instructs the compiler to brings the array w into the cache
     and keep it there */
  __builtin_prefetch (w, 0, 3);
  __builtin_prefetch (w + CACHELINE, 0, 3);
#if CACHELINE < 16
  __builtin_prefetch (w + 2*CACHELINE, 0, 3);
#endif

  for(i = 0; i < np; i++) {

    a = 0.0 ; b = 0.0; c = 0.0;

    for(j = 0; j < n-1; j += 3) {

      a += w[j]   * fi[index[j  ]-1];
      b += w[j+1] * fi[index[j+1]-1];
      c += w[j+2] * fi[index[j+2]-1];

    }

    fo[i] = a + b + c + w[n-1] * fi[index[n-1]-1];

    index += n;

  }

}

