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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <cl_global.h>
#include <cl_complex.h>



__kernel void phase_spiral_apply(
				const int ip_start, 
        const int ip_end,
				__global double2 * psi, const int ldpsi,
				const __global double2 * restrict phase_spiral,
        const __global int * restrict spin_label)
{ 
  
  // This routine expects phase_spiral(np+1:np_part, 1:2) and spin_label(1:nst)

  const int ist = get_global_id(0);
  const int ip0  = get_global_id(1);
  const int ip = ip0 + ip_start;

  if(ip >= ip_end) return;

  const int offset_psi = spin[ist];            // 0 or 1
  const int offset_phase = offset_psi * (ip_end-ip_start);  // assuming phase_spiral(1:np,1:2) in Fortran notation

  psi[(ip<<ldpsi) + 2*ist + offset_psi] = complex_mul(phase_spiral[offset_phase + ip-ip_start], psi[(ip<<ldpsi) + 2*ist + offset_psi]);
  
}
/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
