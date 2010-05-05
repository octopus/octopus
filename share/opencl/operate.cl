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
		       __global double * fo, const int ldfo){

  int ist = get_global_id(0);
  int l = get_global_id(1);
  __global int * index = ri + nn*l;
   
  if(l >= nri) return;

  for(int i = imin[l]; i < imax[l]; i++){
    double a0 = 0.0;
    for(int j = 0; j < nn; j++){
      a0 += weights[j]*fi[ldfi*(index[j] + i) + ist];
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
		       __global double2 * fo, const int ldfo){

  int ist = get_global_id(0);
  int l = get_global_id(1);
  __global int * index = ri + nn*l;

  if(l >= nri) return;

  for(int i = imin[l]; i < imax[l]; i++){
    double2 a0 = 0.0;
    for(int j = 0; j < nn; j++){
      a0 += weights[j]*fi[ldfi*(index[j] + i) + ist];
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

