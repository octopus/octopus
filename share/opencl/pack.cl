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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: pack.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void dpack(const int npsi,
		    const int ipsi,
		    const __global double * src, 
		    __global double * dest){
  int ip = get_global_id(0);

  dest[ipsi + ip*npsi] = src[ip];
}

__kernel void zpack(const int npsi,
		    const int ipsi,
		    const __global double2 * src, 
		    __global double2 * dest){
  int ip = get_global_id(0);

  dest[ipsi + ip*npsi] = src[ip];
}

__kernel void dunpack(const int npsi,
		    const int ipsi,
		    const __global double * src, 
		    __global double * dest){
  int i = get_global_id(0);

  dest[i] = src[ipsi + i*npsi];
}

__kernel void zunpack(const int npsi,
		    const int ipsi,
		    const __global double2 * src, 
		    __global double2 * dest){
  int i = get_global_id(0);

  dest[i] = src[ipsi + i*npsi];
}


/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
