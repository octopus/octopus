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

#ifdef EXT_KHR_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif EXT_AMD_FP64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
/*#pragma OPENCL EXTENSION cl_amd_printf:enable*/
#endif

#if   VECSIZE ==  8
typedef double8 vectype;
#elif VECSIZE ==  4
typedef double4 vectype;
#elif VECSIZE ==  2
typedef double2 vectype;
#elif VECSIZE ==  1
typedef double vectype;
#else
#error VECSIZE not defined
#endif

__kernel void operate(const int nn,
		       const int nri,
		       __global int const * restrict ri,
		       __global int const * imin,
		       __global int const * imax,
		       __constant vectype * restrict weights,
		       __global vectype const * restrict fi, const int ldfi,
		       __global vectype * fo, const int ldfo){

  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int l = get_global_id(1);

  if(l >= nri) return;

  const int imaxl = imax[l];
  for(int i = imin[l]; i < imaxl; i++){
    vectype a0 = (vectype) (0.0);
    for(int j = 0; j < nn; j++){
      a0 += weights[j]*fi[((i + ri[nn*l + j])<<ldfi) + ist];
    }
    fo[(i<<ldfo) + ist] = a0;
  }
  
}

__kernel void operate4(const int nn,
		       const int n1,
		       const int n4,
		       __global int const * restrict ri,
		       __global int const * restrict map1,
		       __global int const * restrict map4,
		       __constant vectype * restrict weights,
		       __global vectype const * restrict fi, const int ldfi,
		       __global vectype * restrict fo, const int ldfo){

  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int k = get_global_id(1);

  if(k < n1) {
    
    const int2 idx = vload2(k, map1);
    // idx.s0 is the point
    // idx.s1 is the map position

    vectype a0 = (vectype) (0.0);
    vectype a1 = (vectype) (0.0);
    
    for(int j = 0; j < nn - 2 + 1; j+=2){
      a0 += weights[j    ]*fi[((idx.s0 + ri[idx.s1 + j]    )<<ldfi) + ist];
      a1 += weights[j + 1]*fi[((idx.s0 + ri[idx.s1 + j + 1])<<ldfi) + ist];
    }
    
    if(nn & 1) a0 += weights[nn - 1]*fi[((idx.s0 + ri[idx.s1 + nn - 1])<<ldfi) + ist];
    
    fo[(idx.s0<<ldfo) + ist] = a0 + a1;
    
  } else if (k - n1 < n4) {
    
    const int2 idx = vload2(k - n1, map4);
    
    vectype a0 = (vectype) (0.0);
    vectype a1 = (vectype) (0.0);
    vectype a2 = (vectype) (0.0);
    vectype a3 = (vectype) (0.0);
    
    for(int j = 0; j < nn; j++){
      a0 += weights[j]*fi[((idx.s0 + 0 + ri[idx.s1 + j])<<ldfi) + ist];
      a1 += weights[j]*fi[((idx.s0 + 1 + ri[idx.s1 + j])<<ldfi) + ist];
      a2 += weights[j]*fi[((idx.s0 + 2 + ri[idx.s1 + j])<<ldfi) + ist];
      a3 += weights[j]*fi[((idx.s0 + 3 + ri[idx.s1 + j])<<ldfi) + ist];
    }
    
    fo[((idx.s0 + 0)<<ldfo) + ist] = a0;
    fo[((idx.s0 + 1)<<ldfo) + ist] = a1;
    fo[((idx.s0 + 2)<<ldfo) + ist] = a2;
    fo[((idx.s0 + 3)<<ldfo) + ist] = a3;
  }
}

__kernel void operate_map(const int nn,
		 	  const int np,
			  __global int const * restrict ri,
			  __global int const * restrict map,
			  __constant vectype * restrict weights,
			  __global vectype const * restrict fi, const int ldfi,
			  __global vectype * restrict fo, const int ldfo,
			  __local int * indexl){

  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int ip  = get_global_id(1);
  const int lip = get_local_id(1);

  __local int * index = indexl + nn*lip;

  for(int j = ist; j < nn; j += nst){
    if(ip < np) index[j] = ri[map[ip] + j];
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  if(ip >= np) return;

  vectype a0 = (vectype) (0.0);
  vectype a1 = (vectype) (0.0);

  for(int j = 0; j < nn - 2 + 1; j += 2){
    a0 += weights[j    ]*fi[((index[j    ] + ip)<<ldfi) + ist];
    a1 += weights[j + 1]*fi[((index[j + 1] + ip)<<ldfi) + ist];
  }

  // if nn is odd, we still have to do the last iteration
  if(nn & 1) a0 += weights[nn - 1]*fi[((index[nn - 1] + ip)<<ldfi) + ist];

  fo[(ip<<ldfo) + ist] = a0 + a1;

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

