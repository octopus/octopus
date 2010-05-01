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

__kernel void dvpsi(const int npsi, const int offset, const __global double * vv, 
		    const __global double * psi, const int ldpsi,
		    __global double * vpsi, const int ldvpsi){
  int i = get_global_id(0);

  double vi = vv[offset + i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*ldvpsi + i] = vi*psi[j*ldpsi + i];
  }

}

__kernel void zvpsi(const int npsi, const int offset, const __global double * vv, 
		    const __global double2 * psi, const int ldpsi,
		    __global double2 * vpsi, const int ldvpsi){
  int i = get_global_id(0);

  double vi = vv[offset + i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*ldvpsi + i] = vi*psi[j*ldpsi + i];
  }

}

__kernel void zvpsi_spinors(const int npsi, 
			    const __global double * vv, const int ldvv,
			    const __global double2 * psi, const int ldpsi,
			    __global double2 * vpsi, const int ldvpsi){
  int i = get_global_id(0);

  double vi1 = vv[         i];
  double vi2 = vv[ldvv   + i];
  double vi3 = vv[2*ldvv + i];
  double vi4 = vv[3*ldvv + i];
  for(int j = 0; j < 2*npsi; j+=2){
    double2 psi1 = psi[      j*ldpsi + i];
    double2 psi2 = psi[(j + 1)*ldpsi + i];

    vpsi[      j*ldvpsi + i] = vi1*psi1 + (double2)(vi3*psi2.x - vi4*psi2.y, vi3*psi2.y + vi4*psi2.x);
    vpsi[(j + 1)*ldvpsi + i] = vi2*psi2 + (double2)(vi3*psi1.x + vi4*psi1.y, vi3*psi1.y - vi4*psi1.x);
  }

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
