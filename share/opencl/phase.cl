/*
 Copyright (C) 2011 X. Andrade

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

 $Id: phase.cl $
*/

#ifdef EXT_KHR_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define HAVE_SINCOS
#elif EXT_AMD_FP64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

__kernel void phase(const int offset, 
		    const __global double * phase, 
		    __global double2 * psi, const int ldpsi){

  const int ist = get_global_id(0);
  const int ip  = get_global_id(1);

#ifdef HAVE_SINCOS
  double cc;
  double ss = sincos(phase[offset + ip], &cc);
#else
#warning "Using single-precision sincos functions."
  float fcc;
  double ss = (double) sincos((float) phase[offset + ip], &fcc);
  double cc = (double) fcc;
#endif

  double2 ff = psi[(ip<<ldpsi) + ist];

  psi[(ip<<ldpsi) + ist] = (double2)(cc*ff.x + ss*ff.y, cc*ff.y - ss*ff.x);

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
