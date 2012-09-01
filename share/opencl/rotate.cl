/*
 Copyright (C) 2011 X. Andrade

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

 $Id: pack.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#include <cl_global.h>

__kernel void drotate_states(const int nst,
			     const int np,
			     __global double const * restrict uu, const int lduu,
			     __global double const * restrict psi, const int ldpsi,
			     __global double * restrict upsi, const int ldupsi){
  const int ist = get_global_id(0);
  const int ip  = get_global_id(1);
  
  double a0 = 0.0;

  for(int ist2 = 0; ist2 < nst; ist2++){
    a0 += uu[ist*lduu + ist2]*psi[ip*ldpsi + ist2];
  }
  
  upsi[ip*ldupsi + ist] = a0;

}

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

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
