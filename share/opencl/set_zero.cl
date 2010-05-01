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

 $Id: set_zero.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void set_zero(__global double * aa){
  aa[get_global_id(0)] = 0.0;
}

__kernel void dset_zero_part(const int nst, const int start, const int end, __global double * aa, const int ldaa){
  int i = get_global_id(0);

  if(i < end - start + 1){
    for(int j = 0; j < nst; j++){
      aa[j*ldaa + start + i] = 0.0;
    }
  }
}

__kernel void zset_zero_part(const int nst, const int start, const int end, __global double2 * aa, const int ldaa){
  int i = get_global_id(0);

  if(i < end - start + 1){
    for(int j = 0; j < nst; j++){
      aa[j*ldaa + start + i] = 0.0;
    }
  }
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
