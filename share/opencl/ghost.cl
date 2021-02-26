/*
 Copyright (C) 2020 S.Ohlmann

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

__kernel void ghost_reorder(const int np, const int offset,
                            const __global int * restrict map,
                            const __global double * restrict xx, const int ldxx,
                            __global double * restrict yy, const int ldyy){

  int ist = get_global_id(0);
  int ip = get_global_id(1) + get_global_size(1)*get_global_id(2);

  if(ip < np) yy[((offset + ip)<<ldyy) + ist] = xx[((map[ip]-1)<<ldxx) + ist];
  //if(ip < np) printf("%d %d %d %d %d\n", ist, ip, map[ip]-1, np, offset);
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
