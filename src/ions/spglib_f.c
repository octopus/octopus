/*
 Copyright (C) 2009 X. Andrade

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

 $Id$
*/

#include <config.h>
#include <spglib.h>
#include "string_f.h" /* Fortran <-> c string compatibility issues */

int FC_FUNC_(spglib_get_multiplicity, SPGLIB_GET_MULTIPLICITY)
     (const double lattice[3][3], const double position[][3],
      const int types[], const int * num_atom, const double * symprec){

  return spg_get_multiplicity(lattice, position, types, *num_atom, *symprec);
}

int FC_FUNC_(spglib_get_symmetry, SPGLIB_GET_SYMMETRY)
     (int rotation[][3][3], double translation[][3], const int * max_size, const double lattice[3][3],
      const double position[][3], const int types[], const int * num_atom, const double * symprec){

  return spg_get_symmetry(rotation, translation, *max_size, lattice, position, types, *num_atom, *symprec);
}

int FC_FUNC_(spglib_get_international, SPGLIB_GET_INTERNATIONAL)(STR_F_TYPE symbol, const double lattice[3][3], const double position[][3],
								 const int types[], const int * num_atom, const double * symprec STR_ARG1){
  char symbol_c[11];
  int space_group = spg_get_international(symbol_c, lattice, position, types, *num_atom, *symprec);
  TO_F_STR1(symbol_c, symbol);
  return space_group;
}

int FC_FUNC_(spglib_get_schoenflies, SPGLIB_GET_SCHOENFLIES)(STR_F_TYPE symbol, const double lattice[3][3], const double position[][3],
							     const int types[], const int * num_atom, const double * symprec STR_ARG1){
  char symbol_c[10];
  int space_group = spg_get_schoenflies(symbol_c, lattice, position, types, *num_atom, *symprec);
  TO_F_STR1(symbol_c, symbol);
  return space_group;
}
