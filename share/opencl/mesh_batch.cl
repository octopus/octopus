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


__kernel void ddot_matrix(const int np,
        __global double const * restrict xx, const int ldxx,
        __global double const * restrict yy, const int ldyy,
        __global double * restrict dot, const int lddot){
  
  int ist = get_global_id(0);
  int jst = get_global_id(1);
  double tmp;

  if(ist >= lddot) return;

  tmp = 0.0;
  for(int ip = 0; ip < np; ip++){
    tmp += xx[(ip<<ldxx) + ist]*yy[(ip<<ldyy) + jst];
  }
  dot[ist + lddot*jst] = tmp;
}

__kernel void zdot_matrix(const int np,
        __global double2 const * restrict xx, const int ldxx,
        __global double2 const * restrict yy, const int ldyy,
        __global double2 * restrict dot, const int lddot){
  
  int ist = get_global_id(0);
  int jst = get_global_id(1);

  if(ist >= lddot) return;
     
  double2 tmp = (double2) (0.0);
  for(int ip = 0; ip < np; ip++){
    double2 a1 = xx[(ip<<ldxx) + ist];
    double2 a2 = yy[(ip<<ldyy) + jst];
    tmp += complex_mul(complex_conj(a1), a2);
  }
  dot[ist + lddot*jst] = tmp;
}


__kernel void zdot_matrix_spinors(const int np,
          const int nst_xx,
          const int nst_yy,
          __global double2 const * restrict xx, const int ldxx,
          __global double2 const * restrict yy, const int ldyy,
          __global double2 * restrict dot, const int lddot){
  
  int ist = get_global_id(0);
  int jst = get_global_id(1);

  if(ist >= nst_xx || jst >= nst_yy) return;

  double2 tmp1 = (double2) (0.0);
  double2 tmp2 = (double2) (0.0);

  for(int ip = 0; ip < np; ip++){
    double2 a1 = xx[(ip<<ldxx) + 2*ist];
    double2 a2 = yy[(ip<<ldyy) + 2*jst];
    double2 b1 = xx[(ip<<ldxx) + 2*ist+1];
    double2 b2 = yy[(ip<<ldyy) + 2*jst+1];

    tmp1 += complex_mul(complex_conj(a1), a2);
    tmp2 += complex_mul(complex_conj(b1), b2);
  }
  dot[ist + lddot*jst] = tmp1 + tmp2;
}




/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
