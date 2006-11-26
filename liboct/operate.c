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

#define CACHELINE 64 
#define UNROLL    6

/* If the __builtin_expect and __builtin_prefetch are not present (which
   should have be caught by the configure script, one needs  to define
   dummy preprocessor macros */
#if !defined(HAVE_BUILTIN_EXPECT)
#define __builtin_expect(a, b) (a)
#endif
#if !defined(HAVE_BUILTIN_PREFETCH)
#define __builtin_prefetch(a, b, c)
#endif

#include "operate_vec.c"

void FC_FUNC(doperate,DOPERATE)(const int * opnp, 
				const int * opn, 
				const double * restrict w, 
				const int * opi, 
				const double * fi, 
				double * restrict fo){

  const int n = opn[0];
  const int np = opnp[0];

  int i, j, nm2;
  const int * restrict index;
  const double * restrict mfi;
  register double a;

  /* this instructs the compiler to brings the array w into the cache
     and keep it there */
  __builtin_prefetch (w, 0, 3);
  __builtin_prefetch (w + CACHELINE/8, 0, 3);
#if CACHELINE < 128
  __builtin_prefetch (w + 2*CACHELINE/8, 0, 3);
#endif

  index = opi;
  mfi   = fi - 1;
  nm2   = n - UNROLL + 1;

  for(i = 0; i < np; i++) {
    a = w[0] * mfi[*index++];

    for(j = 1; j < nm2; j += UNROLL) {
      a += w[j    ] * mfi[*index++];
      a += w[j + 1] * mfi[*index++];
      a += w[j + 2] * mfi[*index++];
      a += w[j + 3] * mfi[*index++];
      a += w[j + 4] * mfi[*index++];
      a += w[j + 5] * mfi[*index++];
    }

    for(; j < n; j++)
      a += w[j] * mfi[*index++];

    fo[i] = a;
  }

}

#ifndef USE_VECTORS

typedef struct {
  double re;
  double im;
} comp;


void FC_FUNC(zoperate,ZOPERATE)(const int * opnp, 
				const int * opn, 
				const double * restrict w, 
				const int * opi, 
				const comp * fi, 
				comp * restrict fo){

  const int n = opn[0];
  const int np = opnp[0];


  int i,j, nm2;
  const int * restrict index;
  const comp * restrict mfi;
  register comp a;

  __builtin_prefetch (w, 0, 3);
  __builtin_prefetch (w + CACHELINE/8, 0, 3);
#if CACHELINE < 128
  __builtin_prefetch (w + 2*CACHELINE/8, 0, 3);
#endif

  index = opi;
  mfi   = fi - 1;
  nm2   = n - UNROLL + 1;

  for(i = 0; i < np ; i++) {
    a.re = w[0] * mfi[*index  ].re;
    a.im = w[0] * mfi[*index++].im;

    for(j = 1; j < nm2; j += UNROLL) {
      a.re += w[j    ] * mfi[*index  ].re;
      a.im += w[j    ] * mfi[*index++].im;
      a.re += w[j + 1] * mfi[*index  ].re;
      a.im += w[j + 1] * mfi[*index++].im;
      a.re += w[j + 2] * mfi[*index  ].re;
      a.im += w[j + 2] * mfi[*index++].im;
      a.re += w[j + 3] * mfi[*index  ].re;
      a.im += w[j + 3] * mfi[*index++].im;
      a.re += w[j + 4] * mfi[*index  ].re;
      a.im += w[j + 4] * mfi[*index++].im;
      a.re += w[j + 5] * mfi[*index  ].re;
      a.im += w[j + 5] * mfi[*index++].im;
    }
    for(; j < n; j++){
      a.re += w[j] * mfi[*index  ].re;
      a.im += w[j] * mfi[*index++].im;
    }

    fo[i].re = a.re;
    fo[i].im = a.im;
  }
}
#endif
