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

__kernel void vpsi(const int offset, 
		   const int np,
		   __global double const * restrict vv, 
		   __global double const * restrict psi, const int ldpsi,
		   __global double * restrict vpsi, const int ldvpsi){

  const int ist = get_global_id(0);
  const int ip  = get_global_id(1);

  if(ip < np){
    vpsi[(ip<<ldvpsi) + ist] += vv[offset + ip]*psi[(ip<<ldpsi) + ist];
  }

}

__kernel void vpsi_spinors(const int np,
			   __global double const * restrict vv, const int ldvv,
			   __global double2 const * restrict psi, const int ldpsi,
			   __global double2 * restrict vpsi, const int ldvpsi){
  const int ist = 2*get_global_id(0);
  const int ip = get_global_id(1);
  
  if(ip < np){

    const double vi1 = vv[         ip];
    const double vi2 = vv[ldvv   + ip];
    const double vi3 = vv[2*ldvv + ip];
    const double vi4 = vv[3*ldvv + ip];
    
    const double2 psi1 = psi[ip*ldpsi + ist];
    const double2 psi2 = psi[ip*ldpsi + ist + 1];
    
    vpsi[ip*ldvpsi + ist] += vi1*psi1 + (double2)(vi3*psi2.x - vi4*psi2.y, vi3*psi2.y + vi4*psi2.x);
    vpsi[ip*ldvpsi + ist + 1] += vi2*psi2 + (double2)(vi3*psi1.x + vi4*psi1.y, vi3*psi1.y - vi4*psi1.x);

  }
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
