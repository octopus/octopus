/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

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
*/

#include <math.h>

/* optimizes the order of the fft
	 p is the maximum prime allowed in n */
void fft_optimize(int *n, int p, int par)
{
	if(*n <= 2) return;

  for(;; (*n)++){
		int i, n2;

		if((par > 0) && (*n % 2 != par)) continue;

		n2 = *n;
		for(i = 2; i<=n2; i++){
			if(n2 % i == 0){
				if(i > p) break;
				n2 = n2 / i; 
				i--; 
			}
		}
		if(i > n2) return;
	}
}
