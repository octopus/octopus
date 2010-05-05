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

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void dvpsi(const int offset, 
		    const __global double * vv, 
		    const __global double * psi, const int ldpsi,
		    __global double * vpsi, const int ldvpsi){

  int ist = get_global_id(0);
  int ip = get_global_id(1);

  double vi = vv[offset + ip];
  vpsi[ip*ldvpsi + ist] = vi*psi[ip*ldpsi + ist];

}

__kernel void zvpsi(const int offset, 
		    const __global double * vv, 
		    const __global double2 * psi, const int ldpsi,
		    __global double2 * vpsi, const int ldvpsi){
  int ist = get_global_id(0);
  int ip = get_global_id(1);

  double vi = vv[offset + ip];
  vpsi[ip*ldvpsi + ist] = vi*psi[ip*ldpsi + ist];

}

__kernel void zvpsi_spinors(const __global double * vv, const int ldvv,
			    const __global double2 * psi, const int ldpsi,
			    __global double2 * vpsi, const int ldvpsi){
  int ist = 2*get_global_id(0);
  int ip = get_global_id(1);

  double vi1 = vv[         ip];
  double vi2 = vv[ldvv   + ip];
  double vi3 = vv[2*ldvv + ip];
  double vi4 = vv[3*ldvv + ip];

  double2 psi1 = psi[ip*ldpsi + ist];
  double2 psi2 = psi[ip*ldpsi + ist + 1];
  
  vpsi[ip*ldvpsi + ist] = vi1*psi1 + (double2)(vi3*psi2.x - vi4*psi2.y, vi3*psi2.y + vi4*psi2.x);
  vpsi[ip*ldvpsi + ist + 1] = vi2*psi2 + (double2)(vi3*psi1.x + vi4*psi1.y, vi3*psi1.y - vi4*psi1.x);
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
