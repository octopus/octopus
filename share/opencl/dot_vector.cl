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

 $Id: vpsi.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#ifdef EXT_KHR_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif EXT_AMD_FP64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

__kernel void ddot_vector(const int np,
			 const __global double * xx, const int ldxx,
			 const __global double * yy, const int ldyy,
			 __global double * dot){
  
  int ist = get_global_id(0);
  double tmp;

  tmp = 0.0;
  for(int ip = 0; ip < np; ip++){
    tmp += xx[(ip<<ldxx) + ist]*yy[(ip<<ldyy) + ist];
  }
  dot[ist] = tmp;
}

__kernel void zdot_vector(const int np,
			 const __global double2 * xx, const int ldxx,
			 const __global double2 * yy, const int ldyy,
			 __global double2 * dot){
  
  int ist = get_global_id(0);
  double2 tmp, a1, a2;

  tmp = (double2) 0.0;
  for(int ip = 0; ip < np; ip++){
    a1 = xx[(ip<<ldxx) + ist];
    a2 = yy[(ip<<ldyy) + ist];
    tmp += (double2)(a1.x*a2.x + a1.y*a2.y, a1.x*a2.y - a1.y*a2.y);
  }
  dot[ist] = tmp;
}

__kernel void ddot_matrix(const int np,
			  __global double const * restrict xx, const int ldxx,
			  __global double const * restrict yy, const int ldyy,
			  __global double * restrict dot, const int lddot){
  
  int ist = get_global_id(0);
  int jst = get_global_id(1);
  double tmp;

  tmp = 0.0;
  for(int ip = 0; ip < np; ip++){
    tmp += xx[(ip<<ldxx) + ist]*yy[(ip<<ldyy) + jst];
  }
  dot[ist + lddot*jst] = tmp;
}

__kernel void zdot_matrix(const int np,
			  __global double2 const * restrict xx, const int ldxx,
			  __global double2 const * restrict yy, const int ldyy,
			  __global double2 * restrict dot, const int lddot){
  
  int ist = get_global_id(0);
  int jst = get_global_id(1);
  double2 tmp, a1, a2;

  tmp = (double2) 0.0;
  for(int ip = 0; ip < np; ip++){
    a1 = xx[(ip<<ldxx) + ist];
    a2 = yy[(ip<<ldyy) + jst];
    tmp += (double2)(a1.x*a2.x + a1.y*a2.y, a1.x*a2.y - a1.y*a2.y);
  }
  dot[ist + lddot*jst] = tmp;
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
