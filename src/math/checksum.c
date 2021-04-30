/*
 Copyright (C) 2010 X. Andrade
 Copyright (C) 2021 S. Ohlmann

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

#include <config.h>
typedef long long checksum_t;

/* This function implements a very stupid checksum function. The only
   important thing is that it produces different results for arrays
   with the same numbers but in different orders. For more serious
   applications a better function must be used. */

void FC_FUNC_(checksum_calculate, CHECKSUM_CALCULATE)(const int * algorithm, const checksum_t * narray, const checksum_t * array, checksum_t * sum){
  int i;
  checksum_t mult;
  *sum = 0;
  mult = 1;
  for(i = 0; i < *narray; i++){
    *sum += mult*array[i];
    mult++;
  }
}
