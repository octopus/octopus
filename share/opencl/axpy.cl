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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <cl_global.h>
#include <cl_complex.h>
#include <cl_rtype.h>

__kernel void X(axpy)(const int np,
		      const rtype aa,
		      const __global rtype * restrict xx, const int ldxx,
		      __global rtype * restrict yy, const int ldyy){
  int ist = get_global_id(0);
  int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);
  
  if(ip < np) yy[(ip<<ldyy) + ist] += MUL(aa, xx[(ip<<ldxx) + ist]);

}

__kernel void X(axpy_vec)(const int np,
			  const __constant rtype * restrict aa, 
			  const __global rtype * restrict xx, const int ldxx,
			  __global rtype * restrict yy, const int ldyy){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);
  
  if(ip < np) yy[(ip<<ldyy) + ist] += MUL(aa[ist], xx[(ip<<ldxx) + ist]);

}

__kernel void X(scal_vec)(const int np,
			  const __constant rtype * restrict aa, 
			  __global rtype * restrict xx, const int ldxx){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);
  
  if(ip < np) xx[(ip<<ldxx) + ist] = MUL(aa[ist], xx[(ip<<ldxx) + ist]);

}

__kernel void X(xpay_vec)(const int np,
			  const __constant rtype * restrict aa, 
			  const __global rtype * restrict xx, const int ldxx,
			  __global rtype * restrict yy, const int ldyy){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);
  
  if(ip < np) yy[(ip<<ldyy) + ist] = xx[(ip<<ldxx) + ist] + MUL(aa[ist], yy[(ip<<ldyy) + ist]);

}

/* The X(batch_axpy_function) kernels should be called on a global grid of (np, ndim, 1) */

__kernel void dbatch_axpy_function(
    const int np,  //< number of mesh points
    const int nst, //< number of states
    const int ndim, //< number of spin components
    __global const double* __restrict xx_buffer, const int ldxx, //< batch of states
    __global const double* __restrict aa_buffer,                 //< buffer of weights
    __global double* __restrict psi_buffer, const int ldpsi      //< single state accumulating the result
) {

  int ip   = get_global_id(0);
  int idim = get_global_id(1);

  double tmp = 0.0;

  if(ip   >= np) return;
  if(idim >= ndim) return;

  for(int ist=0; ist<nst; ist++) {
     tmp += aa_buffer[ist] * xx_buffer[ist + (ip<<ldxx)];
  }
  psi_buffer[ip + (idim<<ldpsi)] += tmp;
}

__kernel void zbatch_axpy_function(
    const int np,   // number of mesh points
    const int nst,  // number of states
    const int ndim, // number of spin components
    __global const double2* __restrict xx_buffer, const int ldxx, // batch of states
    __global const double2* __restrict aa_buffer,                 // buffer of weights
    __global double2* __restrict psi_buffer, const int ldpsi)     // single state accumulating the result
{
  int ip   = get_global_id(0);
  int idim = get_global_id(1);

  double2 tmp = double2(0.0, 0.0);

  if(ip   >= np) return;
  if(idim >= ndim) return;

  for(int ist=0; ist<nst; ist++) {
     tmp += aa_buffer[ist] * xx_buffer[ist + (ip<<ldxx)];
  }
  psi_buffer[ip + (idim<<ldpsi)] += tmp;
}


/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
