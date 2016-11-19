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

 $Id$
*/

#include <cl_global.h>
#include <cl_complex.h>

#define OFFSET_SIZE 6 /* defined in src/hamiltonian/hamiltonian_base.F90 */

__kernel void projector_bra(const int nmat,
			    __global int const * restrict offsets,
			    __global double const * restrict matrix,
			    __global int const * restrict map,
			    __global double const * restrict scal,
			    __global double const * restrict psi, const int ldpsi,
			    __global double * restrict projection, const int ldprojection
			    ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

#ifdef SHARED_MEM
  __local int loff[OFFSET_SIZE];

  for(int ii = get_local_id(0); ii < OFFSET_SIZE; ii += get_local_size(0)) loff[ii] = offsets[OFFSET_SIZE*imat + ii];

  barrier(CLK_LOCAL_MEM_FENCE);

  const int npoints       = loff[0];
  const int nprojs        = loff[1];
  const int matrix_offset = loff[2];
  const int map_offset    = loff[3];
  const int scal_offset   = loff[4];

#else

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

#endif
  
  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

  double aa = 0.0;
  for(int ip = 0; ip < npoints; ip++){
    aa += matrix[matrix_offset + ip + nppj]*psi[((map[map_offset + ip] - 1)<<ldpsi) + ist];
  }
  projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa;

}

__kernel void projector_bra_phase(const int nmat,
				  __global int const * restrict offsets,
				  __global double const * restrict matrix,
				  __global int const * restrict map,
				  __global double const * restrict scal,
				  __global double2 const * restrict psi, const int ldpsi,
				  __global double2 * restrict projection, const int ldprojection,
				  __global double2 const * restrict phases, const int phases_offset
				  ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

#ifdef SHARED_MEM
  __local int loff[OFFSET_SIZE];

  for(int ii = get_local_id(0); ii < OFFSET_SIZE; ii += get_local_size(0)) loff[ii] = offsets[OFFSET_SIZE*imat + ii];

  barrier(CLK_LOCAL_MEM_FENCE);

  const int npoints       = loff[0];
  const int nprojs        = loff[1];
  const int matrix_offset = loff[2];
  const int map_offset    = loff[3];
  const int scal_offset   = loff[4];
  
#else

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];
  
#endif


  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

  double2 aa = 0.0;
  for(int ip = 0; ip < npoints; ip++){
    double2 phasepsi = complex_mul(phases[phases_offset + map_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa += matrix[matrix_offset + ip + nppj]*phasepsi;
  }
  projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa;

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

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ip >= npoints) return;

  double aa = 0.0;
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa += matrix[matrix_offset + ip + npoints*ipj]*projection[ist + ((scal_offset + ipj)<<ldprojection)];
  }

  psi[((map[map_offset + ip] - 1)<<ldpsi) + ist] += aa;

}

__kernel void projector_ket_phase(const int nmat,
				  const int imat_offset,
				  __global int const * restrict offsets,
				  __global double const * restrict matrix,
				  __global int const * restrict map,
				  __global double2 const * restrict projection, const int ldprojection,
				  __global double2 * restrict psi, const int ldpsi,
				  __global double2 const * restrict phases, const int phases_offset
				  ){
  
  const int ist = get_global_id(0);
  const int ip = get_global_id(1);
  const int imat = get_global_id(2) + imat_offset;

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ip >= npoints) return;

  double2 aa = 0.0;
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa += matrix[matrix_offset + ip + npoints*ipj]*projection[ist + ((scal_offset + ipj)<<ldprojection)];
  }

  psi[((map[map_offset + ip] - 1)<<ldpsi) + ist] += complex_mul(complex_conj(phases[phases_offset + map_offset + ip]), aa);

}


__kernel void projector_mix(const int nmat,
			    __global int const * restrict offsets,
			    __global double const * restrict mix,
			    __global double * const restrict projection, const int ldprojection,
			    __global double * restrict mixprojection
			    ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];
  const int mix_offset    = offsets[OFFSET_SIZE*imat + 5];

  if(mix_offset == -1 || ipj >= nprojs) return;
  
  double aa = 0.0;
  for(int jpj = 0; jpj < nprojs; jpj++){
    aa += mix[mix_offset + nprojs*ipj + jpj]*projection[ist + ((scal_offset + jpj)<<ldprojection)];
  }

  mixprojection[ist + ((scal_offset + ipj)<<ldprojection)] = aa;

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

