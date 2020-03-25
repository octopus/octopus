/*
 Copyright (C) 2012-2020 X. Andrade, M. Lueders

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
