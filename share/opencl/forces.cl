/*
 Copyright (C) 2013 X. Andrade

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

 $Id: pack.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#include <cl_global.h>
#include <cl_complex.h>
#include <cl_rtype.h>

__kernel void X(density_gradient)(const int idir,
				  const int nst,
				  const int np,
				  __constant double * weights,
				  const __global rtype * grad_psi,
				  const __global rtype * psi, const int ldpsi,
				  __global double * grad_density){
  
  int ip    = get_global_id(0);

  if(ip >= np) return;

  double dd = 0.0;

  for(int ist = 0; ist < nst; ist ++){
    rtype ff = psi[(ip<<ldpsi) + ist];
    rtype gff = grad_psi[(ip<<ldpsi) + ist];

    dd += weights[ist]*REAL(MUL(CONJ(ff), gff));
  }

  grad_density[np*idir + ip] = dd;

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
