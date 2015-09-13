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

 $Id$
*/

#include <cl_global.h>

__kernel void subarray_gather(const __global int * blength,
			      const __global int * offsets,
			      const __global int * dest,
			      const __global double * array, const int ldarray,
			      __global double * subarray, const int ldsubarray){

  const int ist = get_global_id(0);
  const int ii0 = get_local_id(1);
  const int dii = get_local_size(1);
  const int ibl = get_global_id(2);

  for(int ii = ii0; ii < blength[ibl]; ii += dii){
    subarray[((dest[ibl] + ii)<<ldsubarray) + ist] = array[((offsets[ibl] - 1 + ii)<<ldarray) + ist];
  }

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
