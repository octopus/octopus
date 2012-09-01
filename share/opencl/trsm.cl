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

inline double2 complex_mul(const double2 a, const double2 b){
  return (double2) (a.x*b.x - a.y*b.y, a.y*b.x + a.x*b.y);
}

inline double2 complex_conj_mul(const double2 a, const double2 b){
  return (double2) (a.x*b.x + a.y*b.y, a.y*b.x - a.x*b.y);
}

inline double2 complex_div(const double2 a, const double2 b){
  double2 c = b*b;
  return complex_conj_mul(a, b)/(c.x + c.y);
}

__kernel void dtrsm(const int nst,
		    __global double const * restrict ss, const int ldss,
		    __global double * restrict psi, const int ldpsi){

  const int ip = get_global_id(0);
  
  for(int ist = 0; ist < nst; ist++){
    double a0 = 0.0;
    for(int jst = 0; jst < ist; jst++){
      a0 -= ss[ist*ldss + jst]*psi[ip*ldpsi + jst];
    }
    psi[ip*ldpsi + ist] = (psi[ip*ldpsi + ist] + a0)/ss[ist*ldss + ist];
  }

}

__kernel void ztrsm(const int nst,
		    __global double2 const * restrict ss, const int ldss,
		    __global double2 * restrict psi, const int ldpsi){

  const int ip = get_global_id(0);
  
  for(int ist = 0; ist < nst; ist++){
    double2 a0 = (double2) (0.0);
    for(int jst = 0; jst < ist; jst++){
      a0 -= complex_mul(ss[ist*ldss + jst], psi[ip*ldpsi + jst]);
    }
    psi[ip*ldpsi + ist] = complex_div(psi[ip*ldpsi + ist] + a0, ss[ist*ldss + ist]);
  }

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
