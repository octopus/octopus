/*
 Copyright (C) 2018 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <sys/stat.h>
#include <algorithm>
#include <cassert>

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include "fortran_types.h"

#include "set.hpp"

extern "C" void FC_FUNC_(pseudo_set_init, PSEUDO_SET_INIT)
  (pseudopotential::set ** pseudo_set, STR_F_TYPE const filename_f, fint * ierr STR_ARG1){

  *ierr = 0;
  
  char * filename_c;
  TO_C_STR1(filename_f, filename_c);

  *pseudo_set = new pseudopotential::set(filename_c);
  assert(*pseudo_set);

  free(filename_c);
}

extern "C" void FC_FUNC_(pseudo_set_end, PSEUDO_SET_END)(pseudopotential::set ** pseudo_set){
  delete *pseudo_set;
}
