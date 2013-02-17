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

 $Id: rmmdiis.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#include <cl_global.h>
#include <cl_complex.h>
#include <cl_rtype.h>

__kernel void X(rmmdiis_iter_3)(const int np,
				__global const rtype * restrict resb1,
				__global const rtype * restrict resb2,
				__global const rtype * restrict resb3,
				__global const rtype * restrict psib1,
				__global const rtype * restrict psib2,
				__global const rtype * restrict psib3,
				const int ld,
				__global rtype * restrict sum){
  
  const int ist = get_global_id(0);
  
  rtype s1 = 0.0;
  rtype s2 = 0.0;
  rtype s3 = 0.0;
  rtype s4 = 0.0;
  rtype s5 = 0.0;
  rtype s6 = 0.0;
  rtype s7 = 0.0;
  rtype s8 = 0.0;
  rtype s9 = 0.0;
  rtype sa = 0.0;

  for(int ip = 0; ip < np; ip++){
    rtype r1 = resb1[(ip<<ld) + ist];
    rtype r2 = resb2[(ip<<ld) + ist];
    rtype r3 = resb3[(ip<<ld) + ist];
    rtype p1 = psib1[(ip<<ld) + ist];
    rtype p2 = psib2[(ip<<ld) + ist];
    rtype p3 = psib3[(ip<<ld) + ist];
    
    s1 += MUL(CONJ(r2), r1);
    s2 += MUL(CONJ(p2), p1);
    s3 += MUL(CONJ(r2), r2);
    s4 += MUL(CONJ(p2), p2);
    s5 += MUL(CONJ(r3), r1);
    s6 += MUL(CONJ(p3), p1);
    s7 += MUL(CONJ(r3), r2);
    s8 += MUL(CONJ(p3), p2);
    s9 += MUL(CONJ(r3), r3);
    sa += MUL(CONJ(p3), p3);

  }
  
  sum[(0<<ld) + ist] = s1;
  sum[(1<<ld) + ist] = s2;  
  sum[(2<<ld) + ist] = s3;
  sum[(3<<ld) + ist] = s4; 
  sum[(4<<ld) + ist] = s5;
  sum[(5<<ld) + ist] = s6;
  sum[(6<<ld) + ist] = s7;
  sum[(7<<ld) + ist] = s8;
  sum[(8<<ld) + ist] = s9;
  sum[(9<<ld) + ist] = sa;

}


/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
