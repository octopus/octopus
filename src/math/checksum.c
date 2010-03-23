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

 $Id: checksum.c 4577 2008-09-29 20:18:55Z xavier $
*/

#include <config.h>
typedef long long checksum_t;

/* This function implements a very stupid checksum function. The only
   important thing is that it produces different results for arrays
   with the same numbers but in different orders. For more serious
   applications a better function must be used. */

void FC_FUNC_(checksum_calculate, CHECKSUM_CALCULATE)(const int * algorithm, const int * narray, const unsigned int * array, checksum_t * sum){
  int i;
  checksum_t mult;
  *sum = 0;
  mult = 1;
  for(i = 0; i < *narray; i++){
    *sum += mult*array[i];
    mult++;
  }
  *sum %= (1<<31);
}
