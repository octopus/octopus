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

#define OFFSET_SIZE 6 /* defined in src/hamiltonian/hamiltonian_base.F90 */

// Define the CUDA warpReduce function

#ifdef CUDA
// shuffle instructions wrappers
#if __CUDACC_VER_MAJOR__ >= 9

#define MASK_ALL_WARP 0xFFFFFFFF

#define warpShflDown(var, delta)   __shfl_down_sync (MASK_ALL_WARP, var, delta)

#else

#define warpShflDown(var, delta)   __shfl_down (var, delta)

#endif


__device__ inline  double warpReduce(double val)
{
#pragma unroll
  for (int offset = warpSize/2; offset > 0; offset /= 2){
    val += warpShflDown(val, offset);
  }
  return val;
}

__device__ inline  double2 warpReduce2(double2 val)
{
#pragma unroll
  for (int offset = warpSize/2; offset > 0; offset /= 2){
    val.x += warpShflDown(val.x, offset);
    val.y += warpShflDown(val.y, offset);
  }
  return val;
}

#endif


__kernel void projector_bra(const int nmat,
			    __global int const * restrict offsets,
			    __global double const * restrict matrix,
			    __global int const * restrict map,
			    __global double const * restrict scal,
			    __global double const * restrict psi, const int ldpsi,
			    __global double * restrict projection, const int ldprojection
			    ){
  
#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size;
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

#ifdef CUDA
  const int slice = npoints%my_warp_size==0 ? npoints/my_warp_size : npoints/my_warp_size+1;
  const int start = slice * get_local_id(0) ;
  const int end   = min( slice*(get_local_id(0)+1), npoints );
  const int step = 1;
#else
  const int start = 0;
  const int end = npoints;
  const int step = 1;
#endif

  double aa = 0.0;
  for(int ip = start; ip < end; ip+=step){
    aa += matrix[matrix_offset + ip + nppj]*psi[((map[map_offset + ip] - 1)<<ldpsi) + ist];
  }

#ifdef CUDA
  aa = warpReduce(aa);
  if(get_local_id(0) == 0)
#endif
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
  

#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size;  // the kernel is to be called for (at least) all ist<nst_linear.
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

#ifdef CUDA
  const int slice = npoints%my_warp_size==0 ? npoints/my_warp_size : npoints/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( slice*((get_local_id(0)%my_warp_size)+1), npoints );
  const int step  = 1;
#else
  const int start = 0;
  const int end = npoints;
  const int step = 1;
#endif

  double2 aa = 0.0;
  for(int ip = start; ip < end; ip+=step){
    double2 phasepsi = complex_mul(phases[phases_offset + map_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa += matrix[matrix_offset + ip + nppj]*phasepsi;
  }

#ifdef CUDA
  aa = warpReduce2(aa);
  if(get_local_id(0)%my_warp_size==0) 
#endif
    projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa;

}

__kernel void projector_bra_phase_spiral(const int nmat,
				  __global int const * restrict offsets,
				  __global double const * restrict matrix,
				  __global int const * restrict map,
				  __global double const * restrict scal,
				  __global double2 const * restrict psi, const int ldpsi,
				  __global double2 * restrict projection, const int ldprojection,
				  __global double2 const * restrict phases, const int phases_offset
				  ){
  

#ifdef CUDA
  const int my_warp_size = warpSize;
#else
  const int my_warp_size=1;
#endif

  const int ist = get_global_id(0)/my_warp_size;
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

#ifdef CUDA
  const int slice = npoints%my_warp_size==0 ? npoints/my_warp_size : npoints/my_warp_size+1;
  const int start = slice * ( get_local_id(0)%my_warp_size ) ;
  const int end   = min( slice*((get_local_id(0)%my_warp_size)+1), npoints );
  const int step  = 1;
#else
  const int start = 0;
  const int end = npoints;
  const int step = 1;
#endif

  double2 aa = 0.0;
  for(int ip = start; ip < end; ip+=step){
    double2 phasepsi = complex_mul(phases[phases_offset + map_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa += matrix[matrix_offset + ip + nppj]*phasepsi;
  }

#ifdef CUDA
  aa = warpReduce2(aa);
  if(get_local_id(0)%my_warp_size==0) 
#endif
    projection[ist + ((scal_offset + ipj)<<ldprojection)] = scal[scal_offset + ipj]*aa;

}


__kernel void projector_commutator_bra(const int nmat,
				       __global int const * restrict offsets,
				       __global double const * restrict matrix,
				       __global int const * restrict map,
				       __global double const * restrict scal,
				       __global double const * restrict position,
				       __global double const * restrict psi, const int ldpsi,
				       __global double * restrict projection, const int ldprojection
				       ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

  double aa0 = 0.0;
  double aa1 = 0.0;
  double aa2 = 0.0;
  double aa3 = 0.0;
  
  for(int ip = 0; ip < npoints; ip++){
    aa0 += matrix[matrix_offset + ip + nppj]*psi[((map[map_offset + ip] - 1)<<ldpsi) + ist];
    aa1 += position[(map_offset + ip)*3 + 0]*matrix[matrix_offset + ip + nppj]*psi[((map[map_offset + ip] - 1)<<ldpsi) + ist];
    aa2 += position[(map_offset + ip)*3 + 1]*matrix[matrix_offset + ip + nppj]*psi[((map[map_offset + ip] - 1)<<ldpsi) + ist];
    aa3 += position[(map_offset + ip)*3 + 2]*matrix[matrix_offset + ip + nppj]*psi[((map[map_offset + ip] - 1)<<ldpsi) + ist];
  }

  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0] = scal[scal_offset + ipj]*aa0;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1] = scal[scal_offset + ipj]*aa1;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2] = scal[scal_offset + ipj]*aa2;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3] = scal[scal_offset + ipj]*aa3;
  
}


__kernel void projector_commutator_bra_phase(const int nmat,
					     __global int const * restrict offsets,
					     __global double const * restrict matrix,
					     __global int const * restrict map,
					     __global double const * restrict scal,
					     __global double const * restrict position,
					     __global double2 const * restrict psi, const int ldpsi,
					     __global double2 * restrict projection, const int ldprojection,
					     __global double2 const * restrict phases, const int phases_offset
					     ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);
  const int imat = get_global_id(2);

  const int npoints       = offsets[OFFSET_SIZE*imat + 0];
  const int nprojs        = offsets[OFFSET_SIZE*imat + 1];
  const int matrix_offset = offsets[OFFSET_SIZE*imat + 2];
  const int map_offset    = offsets[OFFSET_SIZE*imat + 3];
  const int scal_offset   = offsets[OFFSET_SIZE*imat + 4];

  if(ipj >= nprojs) return;

  const int nppj = npoints*ipj;

  double2 aa0 = 0.0;
  double2 aa1 = 0.0;
  double2 aa2 = 0.0;
  double2 aa3 = 0.0;
  
  for(int ip = 0; ip < npoints; ip++){
    double2 phasepsi = complex_mul(phases[phases_offset + map_offset + ip], psi[((map[map_offset + ip] - 1)<<ldpsi) + ist]);
    aa0 += matrix[matrix_offset + ip + nppj]*phasepsi;
    aa1 += position[(map_offset + ip)*3 + 0]*matrix[matrix_offset + ip + nppj]*phasepsi;
    aa2 += position[(map_offset + ip)*3 + 1]*matrix[matrix_offset + ip + nppj]*phasepsi;
    aa3 += position[(map_offset + ip)*3 + 2]*matrix[matrix_offset + ip + nppj]*phasepsi;
  }

  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0] = scal[scal_offset + ipj]*aa0;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1] = scal[scal_offset + ipj]*aa1;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2] = scal[scal_offset + ipj]*aa2;
  projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3] = scal[scal_offset + ipj]*aa3;
  
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


__kernel void projector_commutator_ket(const int nmat,
				       const int imat_offset,
				       __global int const * restrict offsets,
				       __global double const * restrict matrix,
				       __global int const * restrict map,
				       __global double const * restrict position,
				       __global double const * restrict projection, const int ldprojection,
				       __global double * restrict cpsi1,
				       __global double * restrict cpsi2,
				       __global double * restrict cpsi3, const int ldpsi
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

  double aa0 = 0.0;
  double aa1 = 0.0;
  double aa2 = 0.0;
  double aa3 = 0.0;
    
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa0 += matrix[matrix_offset + ip + npoints*ipj]*projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0];
    aa1 += matrix[matrix_offset + ip + npoints*ipj]*projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1];
    aa2 += matrix[matrix_offset + ip + npoints*ipj]*projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2];
    aa3 += matrix[matrix_offset + ip + npoints*ipj]*projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3];
  }

  cpsi1[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 0]*aa0 - aa1;
  cpsi2[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 1]*aa0 - aa2;
  cpsi3[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 2]*aa0 - aa3;

}


__kernel void projector_commutator_ket_phase(const int nmat,
					     const int imat_offset,
					     __global int const * restrict offsets,
					     __global double const * restrict matrix,
					     __global int const * restrict map,
					     __global double const * restrict position,
					     __global double2 const * restrict projection, const int ldprojection,
					     __global double2 * restrict cpsi1,
					     __global double2 * restrict cpsi2,
					     __global double2 * restrict cpsi3, const int ldpsi,
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

  double2 aa0 = 0.0;
  double2 aa1 = 0.0;
  double2 aa2 = 0.0;
  double2 aa3 = 0.0;
    
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa0 += matrix[matrix_offset + ip + npoints*ipj]*projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 0];
    aa1 += matrix[matrix_offset + ip + npoints*ipj]*projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 1];
    aa2 += matrix[matrix_offset + ip + npoints*ipj]*projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 2];
    aa3 += matrix[matrix_offset + ip + npoints*ipj]*projection[4*(ist + ((scal_offset + ipj)<<ldprojection)) + 3];
  }

  double2 phase = complex_conj(phases[phases_offset + map_offset + ip]);
  aa0 = complex_mul(phase, aa0);
  aa1 = complex_mul(phase, aa1);
  aa2 = complex_mul(phase, aa2);
  aa3 = complex_mul(phase, aa3);
  
  cpsi1[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 0]*aa0 - aa1;
  cpsi2[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 1]*aa0 - aa2;
  cpsi3[((map[map_offset + ip] - 1)<<ldpsi) + ist] += position[(map_offset + ip)*3 + 2]*aa0 - aa3;

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

