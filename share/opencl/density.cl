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

__kernel void density_real(const int nst,
			   const int np,
			   const int offset,
			   __constant double * restrict weights,
			   const __global double * restrict psi, const int ldpsi,
			   __global double * restrict density){
  
  int ip  = get_global_id(0);
  if(ip >= np) return;

  double dd = 0.0;

  for(int ist = 0; ist < nst; ist++){
    double ff = psi[(ip<<ldpsi) + ist];
    dd += weights[ist]*ff*ff;
  }

  density[offset + ip] += dd;

}

__kernel void density_complex(const int nst,
			      const int np,
			      const int offset,
			      __constant double * restrict weights,
			      const __global double2 * restrict psi, const int ldpsi,
			      __global double * restrict density){
  
  int ip  = get_global_id(0);
  if(ip >= np) return;

  double dd = 0.0;

  for(int ist = 0; ist < nst; ist ++){
    double2 ff = psi[(ip<<ldpsi) + ist];
    ff = ff*ff;
    dd += weights[ist]*(ff.x + ff.y);
  }

  density[offset + ip] += dd;

}

__kernel void density_spinors(const int nst,
			      const int np,
            const int pnp,
			      __constant double * restrict weights,
			      const __global double2 * restrict psi, const int ldpsi,
			      __global double * restrict density){
  
  int ip  = get_global_id(0);
  if(ip >= np) return;

  double dd1 = 0.0;
  double dd2 = 0.0;
  double dd3 = 0.0;
  double dd4 = 0.0;

  for(int ist = 0; ist < nst; ist ++){

    double2 psi1 = psi[(ip<<ldpsi) + 2*ist + 0];
    double2 psi2 = psi[(ip<<ldpsi) + 2*ist + 1];
    double2 ff1 = complex_mul(psi1,complex_conj(psi1));
    double2 ff2 = complex_mul(psi2,complex_conj(psi2));
    double2 ff3 = complex_mul(psi1,complex_conj(psi2));

    dd1 += weights[ist]*ff1.x;
    dd2 += weights[ist]*ff2.x;
    dd3 += weights[ist]*ff3.x;
    dd4 += weights[ist]*ff3.y;
  }

  density[ip + 0*pnp] += dd1;
  density[ip + 1*pnp] += dd2;
  density[ip + 2*pnp] += dd3;
  density[ip + 3*pnp] += dd4;

}



__kernel void current_accumulate(const int nst,
				 const int np,
				 __constant double * restrict weights,
				 const __global double2 * restrict psi, const int ldpsi,
				 const __global double2 * restrict gpsi1, 
				 const __global double2 * restrict gpsi2, 
				 const __global double2 * restrict gpsi3, const int ldgpsi,
				 __global double * restrict current){
  
  int ip  = get_global_id(0);
  if(ip >= np) return;

  double dd1 = 0.0;
  double dd2 = 0.0;
  double dd3 = 0.0;

  for(int ist = 0; ist < nst; ist ++){
    double2 ff = complex_conj(psi[(ip<<ldpsi) + ist]);
    double2 cc1 = complex_mul(ff, gpsi1[(ip<<ldgpsi) + ist]);
    double2 cc2 = complex_mul(ff, gpsi2[(ip<<ldgpsi) + ist]);
    double2 cc3 = complex_mul(ff, gpsi3[(ip<<ldgpsi) + ist]);
    dd1 += weights[ist]*cc1.y;
    dd2 += weights[ist]*cc2.y;
    dd3 += weights[ist]*cc3.y;
  }

  current[ip*3 + 0] = dd1;
  current[ip*3 + 1] = dd2;
  current[ip*3 + 2] = dd3;

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
