/*
 Copyright (C) 2012 X. Andrade

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

#define MCM_POINT 3
#define MCM_COUNT 4

__kernel void dmesh_to_cube(const int nmap,
			    const int stridex,
			    const int stridey,
			    const int stridez,
			    const int centerx,
			    const int centery,
			    const int centerz,
			    __global int const * restrict map,
			    __global double const * restrict mesh_function,
			    __global double * restrict cube_function){
  
  const int imap = get_global_id(0);
  
  if(imap >= nmap) return;

  const int ix = centerx + map[5*imap + 0] - 1;
  const int iy = centery + map[5*imap + 1] - 1;
  const int iz = centerz + map[5*imap + 2] - 1;
  const int ip = map[5*imap + MCM_POINT] - 1;
  const int count = map[5*imap + MCM_COUNT];

  for(int ii = 0; ii < count; ii++){
    cube_function[ix*stridex + iy*stridey + (iz + ii)*stridez] = mesh_function[ip + ii];
  }

}

__kernel void dcube_to_mesh(const int nmap,
			    const int stridex,
			    const int stridey,
			    const int stridez,
			    const int centerx,
			    const int centery,
			    const int centerz,
			    __global int const * restrict map,
			    __global double const * restrict cube_function,
			    __global double * restrict mesh_function){
  
  const int imap = get_global_id(0);
  
  if(imap >= nmap) return;

  const int ix = centerx + map[5*imap + 0] - 1;
  const int iy = centery + map[5*imap + 1] - 1;
  const int iz = centerz + map[5*imap + 2] - 1;
  const int ip = map[5*imap + MCM_POINT] - 1;
  const int count = map[5*imap + MCM_COUNT];

  for(int ii = 0; ii < count; ii++){
    mesh_function[ip + ii] = cube_function[ix*stridex + iy*stridey + (iz + ii)*stridez];
  }

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
