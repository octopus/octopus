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

#include <cl_global.h>

__kernel void projector_bra(const int nmat,
			    __global int const * restrict offsets,
			    __global double const * restrict matrix,
			    __global int const * restrict map,
			    __global double const * restrict scal,
			    __global double const * restrict psi, const int ldpsi,
			    __global double * restrict projection, const int ldprojection,
			    __local  double * restrict lmatrix
			    ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);
  const int itot = get_local_id(0) + get_local_size(0)*get_local_id(1);
  const int ntot = get_local_size(0)*get_local_size(1);

#define BSIZE 128

  __local int lmap[BSIZE];

  for(int ii = itot; ii < 5; ii += ntot) lmap[ii] = offsets[5*imat + ii];
  barrier(CLK_LOCAL_MEM_FENCE);

  const int npoints       = lmap[0];
  const int nprojs        = lmap[1];
  const int matrix_offset = lmap[2];
  const int map_offset    = lmap[3];
  const int scal_offset   = lmap[4];

  const int nppj = npoints*ipj;

  double aa = 0.0;
  
  for(int sp = 0; sp < npoints; sp += BSIZE){
    int size = min(BSIZE, npoints - sp);

    for(int ip = itot; ip < size; ip += ntot) lmap[ip] = (map[map_offset + sp + ip] - 1)<<ldpsi;
    for(int ip = get_local_id(0); ip < size; ip += get_local_size(0)) {
      if(ipj < nprojs) lmatrix[ip + ipj*BSIZE] = matrix[matrix_offset + nppj + sp + ip];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int ip = 0; ip < size; ip++){
      if(ipj < nprojs) aa += lmatrix[ip + ipj*BSIZE]*psi[lmap[ip] + ist];
    }
  }
  
  if(ipj < nprojs) projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa;

#undef BSIZE

}

__kernel void projector_ket(const int nmat,
			    const int imat_offset,
			    __global int const * restrict offsets,
			    __global double const * restrict matrix,
			    __global int const * restrict map,
			    __global double const * restrict projection, const int ldprojection,
			    __global double * restrict psi, const int ldpsi
			    ){
  
  const int ist = get_global_id(0);
  const int ip = get_global_id(1);
  const int imat = get_global_id(2) + imat_offset;

  const int npoints       = offsets[5*imat + 0];
  const int nprojs        = offsets[5*imat + 1];
  const int matrix_offset = offsets[5*imat + 2];
  const int map_offset    = offsets[5*imat + 3];
  const int scal_offset   = offsets[5*imat + 4];

  if(ip >= npoints) return;

  double aa = 0.0;
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa += matrix[matrix_offset + ip + npoints*ipj]*projection[ist + ((scal_offset + ipj)<<ldprojection)];
  }

  psi[((map[map_offset + ip] - 1)<<ldpsi) + ist] += aa;

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

