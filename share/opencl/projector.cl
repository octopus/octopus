/*
 Copyright (C) 2010 X. Andrade

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

 $Id: projector.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void projector_bra(const int npoints,
			    const int nprojs,
			    __global const int * map,
			    __global const double * scal,
			    __global const double * matrix, const int ldmatrix,
			    __global const double * psi, const int ldpsi,
			    __global double * projection, const int ldprojection,
			    __local double * lprojection
			    ){
  
  int ist = get_global_id(0);
  int ipj = get_global_id(1);
  int k = get_global_id(2);
  int nk = get_global_size(2);


  
  double aa = 0.0;
  for(int ip = k; ip < npoints; ip += nk){
      if(ipj < nprojs) aa += matrix[ip + ldmatrix*ipj]*psi[ldpsi*(map[ip] - 1) + ist];
  }
  aa *= scal[ipj];
  
  // this can be improved by doing a parallel reduction
  if(ipj < nprojs && k == 0) lprojection[ist + ldprojection*ipj] = aa;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int k2 = 1; k2 < nk; k2++){
    if(ipj < nprojs && k2 == k) lprojection[ist + ldprojection*ipj] += aa;
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if(ipj < nprojs && k == 0) projection[ist + ldprojection*ipj] = lprojection[ist + ldprojection*ipj];

}

__kernel void projector_ket(const int npoints,
			    const int nprojs,
			    __global const int * map,
			    __global const double * matrix, const int ldmatrix,
			    __global const double * projection, const int ldprojection,
			    __global double * psi, const int ldpsi
			    ){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1);

  if(ip >= npoints) return;

  double aa = 0.0;
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa += matrix[ip + ldmatrix*ipj]*projection[ist + ldprojection*ipj];
  }
  psi[ldpsi*(map[ip] - 1) + ist] += aa;


}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

