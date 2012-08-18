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

__kernel void operate4(const int nn,
		       const int n1,
		       const int n1n4,
		       __global int const * restrict ri,
		       __global int const * restrict map_split,
		       __constant double * restrict weights,
		       __global double const * restrict fi, const int ldfi,
		       __global double * restrict fo, const int ldfo,
		       __local int * indexl){

  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int k = get_global_id(1);
  const int lk = get_local_id(1);
  
  __local int * restrict index = indexl + nn*lk;

  const int2 idx = vload2(k, map_split);
  // idx.s0 is the point index
  // idx.s1 is the map position

  if(k < n1n4){
    int j = ist<<2;
    for(; j < nn - 4 + 1; j += (nst<<2)){
      vstore4(vload4(0, ri + idx.s1 + j), 0, index + j);
      }
    for(; j < nn; j+= nst){
      index[j] = ri[idx.s1 + j];
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE);
  
  if(k < n1) {

    double a0 = (double) (0.0);
    double a1 = (double) (0.0);
    
    for(int j = 0; j < nn - 2 + 1; j+=2){
      int2 rij = vload2(0, index + j);
      a0 += weights[j    ]*fi[((idx.s0 + rij.s0)<<ldfi) + ist];
      a1 += weights[j + 1]*fi[((idx.s0 + rij.s1)<<ldfi) + ist];
    }
    
    if(nn & 1) a0 += weights[nn - 1]*fi[((idx.s0 + index[nn - 1])<<ldfi) + ist];
    
    fo[(idx.s0<<ldfo) + ist] = a0 + a1;
    
  } else if (k < n1n4) {
    
    double a0 = (double) (0.0);
    double a1 = (double) (0.0);
    double a2 = (double) (0.0);
    double a3 = (double) (0.0);
    
    for(int j = 0; j < nn; j++){
      a0 += weights[j]*fi[((idx.s0 + 0 + index[j])<<ldfi) + ist];
      a1 += weights[j]*fi[((idx.s0 + 1 + index[j])<<ldfi) + ist];
      a2 += weights[j]*fi[((idx.s0 + 2 + index[j])<<ldfi) + ist];
      a3 += weights[j]*fi[((idx.s0 + 3 + index[j])<<ldfi) + ist];
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
			  __constant double * restrict weights,
			  __global double const * restrict fi, const int ldfi,
			  __global double * restrict fo, const int ldfo,
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

  double a0 = (double) (0.0);
  double a1 = (double) (0.0);

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

