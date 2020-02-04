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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <cl_global.h>
#include <cl_complex.h>


__kernel void complex_conj_combine(const int nst, __global const double2 * __restrict src, __global double2 * __restrict dest){
  int ist = get_global_id(0);
  
  if(ist < nst) dest[ist] = complex_conj( src[2*ist] + src[2*nst + 2*ist+1] );

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
