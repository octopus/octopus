/*
 Copyright (C) 2013 X. Andrade

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

 $Id: boundaries.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#include <cl_global.h>

__kernel void boundaries_periodic(const int nper, 
				  __global const int * __restrict per_points,
				  __global double * __restrict ff,
				  const int ldff){
  const int ist  = get_global_id(0);
  const int iper = get_global_id(1);
  
  if(iper >= nper) return;
  
  const int ip_bnd = per_points[iper*2    ] - 1;
  const int ip_inn = per_points[iper*2 + 1] - 1;

  ff[(ip_bnd<<ldff) + ist] = ff[(ip_inn<<ldff) + ist];

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
