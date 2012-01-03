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

 $Id: copy.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#ifdef EXT_KHR_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif EXT_AMD_FP64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

__kernel void copy(const __global double * restrict xx, const int ldxx,
		   __global double * restrict yy, const int ldyy){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1);
  
  yy[(ip<<ldyy) + ist] = xx[(ip<<ldxx) + ist];

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
