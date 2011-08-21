/*
 Copyright (C) 2011 X. Andrade

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

 $Id: pack.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#ifdef EXT_KHR_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif EXT_AMD_FP64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

__kernel void get_points(const int sp,
			 const int ep,
			 const int offset,
			 __global double const * restrict psi, const int ldpsi,
			 __global double * restrict points, const int ldpoints){
  const int ist = get_global_id(0);
  const int ip  = get_global_id(1);
  const int sip = ip + sp - 1;

  points[ldpoints*ip + ist + offset] = psi[ldpsi*sip + ist];

}

__kernel void set_points(const int sp,
			 const int ep,
			 const int offset,
			 __global double const * restrict points, const int ldpoints,
			 __global double * restrict psi, const int ldpsi){
  const int ist = get_global_id(0);
  const int ip  = get_global_id(1);
  const int sip = ip + sp - 1;

  psi[ldpsi*sip + ist] = points[ldpoints*ip + ist + offset];

}
/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
