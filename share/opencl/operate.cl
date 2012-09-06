/*
 Copyright (C) 2010 X. Andrade, N. Suberviola

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

 $Id: operate.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#include <cl_global.h>

#ifndef STENCIL_SIZE
#error Internal error: STENCIL_SIZE not declared
#endif

__kernel void operate(const int nn,
		      const int nri,
		      __global int const * restrict ri,
		      __global int const * imin,
		      __global int const * imax,
		      __constant double * restrict weights,
		      __global double const * restrict fi, const int ldfi,
		      __global double * fo, const int ldfo){
  
  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int l = get_global_id(1);
  
  if(l >= nri) return;
  
  const int imaxl = imax[l];
  for(int i = imin[l]; i < imaxl; i++){
    double a0 = (double) (0.0);
    for(int j = 0; j < nn; j++){
      a0 += weights[j]*fi[((i + ri[nn*l + j])<<ldfi) + ist];
    }
    fo[(i<<ldfo) + ist] = a0;
  }
  
}

__kernel void operate_map(const int nn,
		 	  const int np,
			  __global int const * restrict ri,
			  __global int const * restrict map,
			  __constant double * restrict weights,
			  __global double const * restrict fi, const int ldfi,
			  __global double * restrict fo, const int ldfo,
			  __local int * indexl){

  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int ip  = get_global_id(1);
  const int lip = get_local_id(1);

  __local int * index = indexl + nn*lip;

  if(ip < np){
    for(int j = ist; j < nn; j += nst){
      index[j] = ri[map[ip] + j];
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE);
  
  if(ip < np) {

#if STENCIL_SIZE > 27

    double a0 = (double) (0.0);
    double a1 = (double) (0.0);
    
    for(int j = 0; j < nn - 2 + 1; j += 2){
      a0 += weights[j    ]*fi[((index[j    ] + ip)<<ldfi) + ist];
      a1 += weights[j + 1]*fi[((index[j + 1] + ip)<<ldfi) + ist];
    }
    
    // if nn is odd, we still have to do the last iteration
    if(nn & 1) a0 += weights[nn - 1]*fi[((index[nn - 1] + ip)<<ldfi) + ist];
    
    fo[(ip<<ldfo) + ist] = a0 + a1;
    
#else

    double a0 = (double) (0.0);

    a0 += weights[ 0]*fi[((index[ 0] + ip)<<ldfi) + ist];
#if STENCIL_SIZE > 1
    a0 += weights[ 1]*fi[((index[ 1] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 2
    a0 += weights[ 2]*fi[((index[ 2] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 3
    a0 += weights[ 3]*fi[((index[ 3] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 4
    a0 += weights[ 4]*fi[((index[ 4] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 5
    a0 += weights[ 5]*fi[((index[ 5] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 6
    a0 += weights[ 6]*fi[((index[ 6] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 7
    a0 += weights[ 7]*fi[((index[ 7] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 8
    a0 += weights[ 8]*fi[((index[ 8] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 9
    a0 += weights[ 9]*fi[((index[ 9] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 10
    a0 += weights[10]*fi[((index[10] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 11
    a0 += weights[11]*fi[((index[11] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 12
    a0 += weights[12]*fi[((index[12] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 13
    a0 += weights[13]*fi[((index[13] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 14
    a0 += weights[14]*fi[((index[14] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 15
    a0 += weights[15]*fi[((index[15] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 16
    a0 += weights[16]*fi[((index[16] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 17
    a0 += weights[17]*fi[((index[17] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 18
    a0 += weights[18]*fi[((index[18] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 19
    a0 += weights[19]*fi[((index[19] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 20
    a0 += weights[20]*fi[((index[20] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 21
    a0 += weights[21]*fi[((index[21] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 22
    a0 += weights[22]*fi[((index[22] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 23
    a0 += weights[23]*fi[((index[23] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 24
    a0 += weights[24]*fi[((index[24] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 25
    a0 += weights[25]*fi[((index[25] + ip)<<ldfi) + ist];
#endif
#if STENCIL_SIZE > 26
    a0 += weights[26]*fi[((index[26] + ip)<<ldfi) + ist];
#endif

    fo[(ip<<ldfo) + ist] = a0;

#endif

  }

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

