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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: spglib_f.c 2146 2006-05-23 17:36:00Z xavier $
*/

#include <config.h>
#include <spglib.h>

int FC_FUNC_(spglib_get_max_multiplicity, SPGLIB_GET_MAX_MULTIPLICITY)
     (const double lattice[3][3], const double position[][3],
      const int types[], const int * num_atom, const double * symprec){

  return spg_get_max_multiplicity(lattice, position, types, *num_atom, *symprec);
}

int FC_FUNC_(spglib_get_symmetry, SPGLIB_GET_SYMMETRY)
     (int rotation[][3][3], double translation[][3], const int * max_size, const double lattice[3][3],
      const double position[][3], const int types[], const int * num_atom, const double * symprec){

  return spg_get_symmetry(rotation, translation, *max_size, lattice, position, types, *num_atom, *symprec);
}

void FC_FUNC_(spglib_show_symmetry, SPGLIB_SHOW_SYMMETRY)(const double lattice[3][3], const double position[][3],
							  const int types[], const int * num_atom, const double * symprec){
  spg_show_symmetry(lattice, position, types, *num_atom, *symprec);
}

int FC_FUNC_(spglib_get_group_number, SPGLIB_GET_GROUP_NUMBER)(const double lattice[3][3], const double position[][3],
							       const int types[], const int * num_atom, const double * symprec){
  char symbol[11];
  return spg_get_international(symbol, lattice, position, types, *num_atom, *symprec);
}
