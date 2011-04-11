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
#endif

__kernel void operate(const int nn,
		       const int nri,
		       __global const int * ri,
		       __global const int * imin,
		       __global const int * imax,
		       __constant double * weights,
		       __global const double * fi, const int ldfi,
		       __global double * fo, const int ldfo,
		       __local int * indexl,
		       __local double * weightsl){

  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int l = get_global_id(1);

  // copy the indices and the weights to local memory
  __local int * index = indexl + nn*get_local_id(1);
  for(int j = ist; j < nn; j+=nst){
    index[j] = ri[nn*l + j]*ldfi;
    weightsl[j] = weights[j];
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  if(l >= nri) return;

  const int imaxl = imax[l];
  for(int i = imin[l]; i < imaxl; i++){
    double a0 = 0.0;
    for(int j = 0; j < nn; j++){
      a0 += weightsl[j]*fi[ldfi*i + index[j] + ist];
    }
    fo[ldfo*i + ist] = a0;
  }
  
}

__kernel void operate_map(const int nn,
			  const int np,
			  __global const int * ri,
			  __global const int * map,
			  __constant double * weights,
			  __global const double * fi, const int ldfi,
			  __global double * fo, const int ldfo,
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

  double a0 = 0.0;
  double a1 = 0.0;

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

