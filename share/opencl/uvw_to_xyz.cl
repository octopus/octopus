/*
 Copyright (C) 2020 S. Ohlmann

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

#ifndef DIMENSION
#error DIMENSION not defined!
#elif DIMENSION > 3
#error Transformation not defined for more than 3 dimensions!
#endif

__kernel void uvw_to_xyz(const int np,
           __global double const * restrict matrix,
           __global double * gradu, const int ldgradu,
           __global double * gradx, const int ldgradx
#if DIMENSION > 1
           ,__global double * gradv, const int ldgradv
           ,__global double * grady, const int ldgrady
#endif
#if DIMENSION > 2
           ,__global double * gradw, const int ldgradw
           ,__global double * gradz, const int ldgradz
#endif
           ){
  int ist = get_global_id(0);
  int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);
  double tmp[DIMENSION];
  for(int idim = 0; idim < DIMENSION; ++idim) {
    tmp[idim] = 0.0;
  }
  
  if(ip < np) {
#if DIMENSION == 1
    gradx[(ip<<ldgradx) + ist] = matrix[0] * gradu[(ip<<ldgradu) + ist];
#elif DIMENSION == 2
    tmp[0] = matrix[0]*gradu[(ip<<ldgradu) + ist] + matrix[2]*gradv[(ip<<ldgradv) + ist];
    tmp[1] = matrix[1]*gradu[(ip<<ldgradu) + ist] + matrix[3]*gradv[(ip<<ldgradv) + ist];
    gradx[(ip<<ldgradx) + ist] = tmp[0];
    grady[(ip<<ldgrady) + ist] = tmp[1];
#elif DIMENSION == 3
    tmp[0] = matrix[0]*gradu[(ip<<ldgradu) + ist] + matrix[3]*gradv[(ip<<ldgradv) + ist] + matrix[6]*gradw[(ip<<ldgradw) + ist];
    tmp[1] = matrix[1]*gradu[(ip<<ldgradu) + ist] + matrix[4]*gradv[(ip<<ldgradv) + ist] + matrix[7]*gradw[(ip<<ldgradw) + ist];
    tmp[2] = matrix[2]*gradu[(ip<<ldgradu) + ist] + matrix[5]*gradv[(ip<<ldgradv) + ist] + matrix[8]*gradw[(ip<<ldgradw) + ist];
    gradx[(ip<<ldgradx) + ist] = tmp[0];
    grady[(ip<<ldgrady) + ist] = tmp[1];
    gradz[(ip<<ldgradz) + ist] = tmp[2];
#endif
  }
}

/*
__kernel void zrotate_states(const int nst,
           const int np,
           __global double2 const * restrict uu, const int lduu,
           __global double2 const * restrict psi, const int ldpsi,
           __global double2 * restrict upsi, const int ldupsi){
  const int ist = get_global_id(0);
  const int ip  = get_global_id(1);
  
  double2 a0 = (double2) (0.0);

  for(int ist2 = 0; ist2 < nst; ist2++){
    double2 xx = uu[ist*lduu + ist2];
    double2 yy = psi[ip*ldpsi + ist2];
    a0 += (double2) (xx.s0*yy.s0 - xx.s1*yy.s1, xx.s0*yy.s1 + xx.s1*yy.s0);
  }
  
  upsi[ip*ldupsi + ist] = a0;

}
*/

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
