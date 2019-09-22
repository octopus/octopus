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

#include "element.hpp"

extern "C" void FC_FUNC_(pseudo_set_init_low, PSEUDO_SET_INIT_LOW)
  (pseudopotential::set ** pseudo_set, STR_F_TYPE const filename_f, const fint * automatic, fint * ierr STR_ARG1){

  *ierr = 0;
  
  char * filename_c;
  TO_C_STR1(filename_f, filename_c);

  *pseudo_set = new pseudopotential::set(filename_c, bool(*automatic));
  assert(*pseudo_set);

  free(filename_c);
}

extern "C" void FC_FUNC_(pseudo_set_nullify, PSEUDO_SET_NULLIFY)(pseudopotential::set ** pseudo_set){
  *pseudo_set = NULL;
}

extern "C" void FC_FUNC_(pseudo_set_end, PSEUDO_SET_END)(pseudopotential::set ** pseudo_set){
  if(*pseudo_set){
    delete *pseudo_set;
    *pseudo_set = NULL;
  }
}

extern "C" fint FC_FUNC_(pseudo_set_has_low, PSEUDO_SET_HAS_LOW)
  (const pseudopotential::set ** pseudo_set, pseudopotential::element **el){
  return (*pseudo_set)->has(**el);
}

extern "C" void FC_FUNC_(pseudo_set_file_path_low, PSEUDO_SET_FILE_PATH_LOW)
  (const pseudopotential::set ** pseudo_set, pseudopotential::element **el, STR_F_TYPE const path_f STR_ARG1){
  std::string path = (*pseudo_set)->file_path(**el);
  TO_F_STR1(path.c_str(), path_f);
}

extern "C" fint FC_FUNC_(pseudo_set_lmax, PSEUDO_SET_LMAX)
  (const pseudopotential::set ** pseudo_set, pseudopotential::element **el){
  return (*pseudo_set)->lmax(**el);
}

extern "C" fint FC_FUNC_(pseudo_set_llocal, PSEUDO_SET_LLOCAL)
  (const pseudopotential::set ** pseudo_set, pseudopotential::element **el){
  return (*pseudo_set)->llocal(**el);
}

extern "C" double FC_FUNC_(pseudo_set_spacing, PSEUDO_SET_SPACING)
  (const pseudopotential::set ** pseudo_set, pseudopotential::element **el, double * etol){
  return (*pseudo_set)->spacing(**el, *etol);
}

extern "C" double FC_FUNC_(pseudo_set_radius, PSEUDO_SET_RADIUS)
  (const pseudopotential::set ** pseudo_set, pseudopotential::element **el){
  return (*pseudo_set)->radius(**el);
}
