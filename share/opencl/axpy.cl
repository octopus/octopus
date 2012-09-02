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

 $Id: vpsi.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#include <cl_global.h>
#include <cl_complex.h>
#include <cl_rtype.h>

__kernel void X(axpy)(const rtype aa,
		      const __global rtype * restrict xx, const int ldxx,
		      __global rtype * restrict yy, const int ldyy){
  int ist = get_global_id(0);
  int ip = get_global_id(1);

  yy[ip*ldyy + ist] += MUL(aa, xx[ip*ldxx + ist]);

}

__kernel void X(axpy_vec)(const __constant rtype * restrict aa, 
			  const __global rtype * restrict xx, const int ldxx,
			  __global rtype * restrict yy, const int ldyy){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1);
  
  yy[(ip<<ldyy) + ist] += MUL(aa[ist], xx[(ip<<ldxx) + ist]);

}

__kernel void X(scal_vec)(const __constant rtype * restrict aa, 
			  __global rtype * restrict xx, const int ldxx){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1);
  
  xx[(ip<<ldxx) + ist] = MUL(aa[ist], xx[(ip<<ldxx) + ist]);

}


/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
