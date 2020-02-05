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
#include <cl_complex.h>
#include <cl_rtype.h>

__kernel void X(curl)(const int np,
           __global rtype * restrict ff, const int ldff,
           __global const rtype * restrict gradx, const int ldgradx,
           __global const rtype * restrict grady, const int ldgrady,
           __global const rtype * restrict gradz, const int ldgradz
           ){
  int ist = get_global_id(0);
  int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);

  /* there can be several vectors in one batch -> idim is the spatial dimension in each vector */
  int idim = ist % 3;
  int offset = ist / 3;
  /* get the right pointers for the current component of the curl */
  const rtype *grad1, *grad2;
  const int *ldgrad1, *ldgrad2;
  if(idim == 0) {
    grad1 = grady;
    grad2 = gradz;
    ldgrad1 = &ldgrady;
    ldgrad2 = &ldgradz;
  } else if (idim == 1) {
    grad1 = gradz;
    grad2 = gradx;
    ldgrad1 = &ldgradz;
    ldgrad2 = &ldgradx;
  } else {
    grad1 = gradx;
    grad2 = grady;
    ldgrad1 = &ldgradx;
    ldgrad2 = &ldgrady;
  }
  
  /* now subtract the right component: the modulo is needed to get the right
   * anticyclic indices */
  if(ip < np) {
    ff[(ip<<ldff) + ist] = grad1[(ip<<*ldgrad1) + ((idim+2)%3) + offset] - grad2[(ip<<*ldgrad2) + ((idim+1)%3) + offset];
  }
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
/* vim: set filetype=c: */
