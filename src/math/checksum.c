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

void FC_FUNC_(checksum_calculate, CHECKSUM_CALCULATE)(const int * algorithm, const int * narray, const unsigned int * array, checksum_t * sum){
  int i;

  *sum = 0;
  for(i = 0; i < *narray; i++) *sum += array[i];

}

int FC_FUNC_(checksum_compare_int, CHECKSUM_COMPARE_INT)(const int * algorithm, const checksum_t * sum1, const checksum_t * sum2){
  return (sum1 == sum2);
}
