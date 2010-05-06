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

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void doperate(const int nn,
		       const int nri,
		       __global const int * ri,
		       __global const int * imin,
		       __global const int * imax,
		       __constant const double * weights,
		       __global const double * fi, const int ldfi,
		       __global double * fo, const int ldfo,
		       __local int * indexl,
		       __local double * weightsl){

  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int l = get_global_id(1);

  // copy the indexes and the weights to local memory
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

__kernel void zoperate(const int nn,
		       const int nri,
		       __global const int * ri,
		       __global const int * imin,
		       __global const int * imax,
		       __constant const double * weights,
		       __global const double2 * fi, const int ldfi,
		       __global double2 * fo, const int ldfo,
		       __local  int * indexl,
		       __local double * weightsl){

  const int ist = get_global_id(0);
  const int nst = get_global_size(0);
  const int l = get_global_id(1);

  // copy the indexes and the weights to local memory
  __local int * index = indexl + nn*get_local_id(1);
  for(int j = ist; j < nn; j+=nst){
    index[j] = ri[nn*l + j]*ldfi;
    weightsl[j] = weights[j];
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  if(l >= nri) return;

  const int imaxl = imax[l];
  for(int i = imin[l]; i < imaxl; i++){
    double2 a0 = 0.0;
    for(int j = 0; j < nn; j++){
      a0 += weightsl[j]*fi[ldfi*i + index[j] + ist];
    }
    fo[ldfo*i + ist] = a0;
  }
  
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

