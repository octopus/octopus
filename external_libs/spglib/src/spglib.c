/* spglib.c */
/* Copyright (C) 2008 Atsushi Togo */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#ifdef DEBUG
#include "debug.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spglib.h"
#include "bravais.h"
#include "mathfunc.h"
#include "primitive.h"
#include "cell.h"
#include "symmetry.h"
#include "pointgroup.h"
#include "spacegroup.h"
#include "spacegroup_database.h"

/*
  ------------------------------------------------------------------

  lattice: Lattice vectors (in Cartesian)

  [ [ a_x, b_x, c_x ],
    [ a_y, b_y, c_y ],
    [ a_z, b_z, c_z ] ]

  position: Atomic positions (in fractional coordinates)
  
  [ [ x1_a, x1_b, x1_c ], 
    [ x2_a, x2_b, x2_c ], 
    [ x3_a, x3_b, x3_c ], 
    ...                   ]

  types: Atom types, i.e., species identified by number

  [ type_1, type_2, type_3, ... ]

  rotation: Rotation matricies of symmetry operations

    each rotation is:
    [ [ r_11, r_12, r_13 ],
      [ r_21, r_22, r_23 ],
      [ r_31, r_32, r_33 ] ]

  translation: Translation vectors of symmetry operations

    each translation is:
    [ t_1, t_2, t_3 ]

  symprec: Tolerance of atomic positions (in fractional coordinate)
           in finding symmetry operations

  ------------------------------------------------------------------

  Definitio of the operation:
    r : rotation     3x3 matrix
    t : translation  vector

    x_new = r * x + t:
    [ x_new_1 ]   [ r_11 r_12 r_13 ]   [ x_1 ]   [ t_1 ]
    [ x_new_2 ] = [ r_21 r_22 r_23 ] * [ x_2 ] + [ t_2 ]
    [ x_new_3 ]   [ r_31 r_32 r_33 ]   [ x_3 ]   [ t_3 ]

  ------------------------------------------------------------------
 */




/* Return symmetry operations, which are rotations and translations */
/* in fractional coordinates. The first indices of rotation and transtion */
/* correspond each other. */
int spg_get_symmetry(int rotation[][3][3], double translation[][3],
		     const int max_size, const double lattice[3][3],
		     const double position[][3], const int types[],
		     const int num_atom, const double symprec)
{
  /* max_size is used for allocating memory space for returning symmetry operations. */

  int i, j, size;
  Symmetry symmetry;
  Bravais bravais;
  Cell cell;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);
  bravais = brv_get_brv_lattice(cell.lattice, symprec);
  symmetry = sym_get_operation(&bravais, &cell, symprec);

  if (symmetry.size > max_size) {
    fprintf(stderr, "spglib: Indicated max size(=%d) is less than number ", max_size);
    fprintf(stderr, "spglib: of symmetry operations(=%d).\n", symmetry.size);
    return 0;
  }

  for (i = 0; i < symmetry.size; i++) {
    mat_copy_matrix_i3(rotation[i], symmetry.rot[i]);
    for (j = 0; j < 3; j++)
      translation[i][j] = symmetry.trans[i][j];
  }

  size = symmetry.size;

  cel_delete_cell(&cell);
  sym_delete_symmetry(&symmetry);

  return size;
}

/* Retrun bravais lattice estimated from lattice vectors only. */
/* If crystal structure is virtural structure, */
/* internal atomic positions are returied to determine space group */
/* and thus bravais lattice. */
void spg_get_bravais_lattice(double bravais_lattice[3][3],
			     const double lattice[3][3],
                             const double symprec)
{
  Bravais bravais;

  bravais = brv_get_brv_lattice(lattice, symprec);
  mat_copy_matrix_d3(bravais_lattice, bravais.lattice);
}

void spg_get_smallest_lattice(double smallest_lattice[3][3],
			      const double lattice[3][3],
			      const double symprec)
{
  brv_smallest_lattice_vector(smallest_lattice, lattice, symprec);
}

/* Return exact number of symmetry operations. */
/* This is slower than spg_get_max_multiplicity. */
int spg_get_multiplicity(const double lattice[3][3], const double position[][3],
			 const int types[], const int num_atom,
			 const double symprec)
{
  Symmetry symmetry;
  Bravais bravais;
  Cell cell;
  int size;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);
  bravais = brv_get_brv_lattice(cell.lattice, symprec);
  symmetry = sym_get_operation(&bravais, &cell, symprec);

  size = symmetry.size;

  cel_delete_cell(&cell);
  sym_delete_symmetry(&symmetry);

  return size;
}

/* Return the possiblly largest number of symmetry operations. */
/* This is faster than spg_check_symmetry. */
int spg_get_max_multiplicity(const double lattice[3][3], const double position[][3],
			     const int types[], const int num_atom,
			     const double symprec)
{
  Cell cell;
  int multiplicity;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);
  /* 96 is a magic number, which is the twice of the number of rotations */
  /* in the highest point symmetry Oh. */
  multiplicity = sym_get_multiplicity(&cell, symprec) * 96;
  cel_delete_cell(&cell);
  
  return multiplicity;
}

/* Find a primitive cell in the input cell. */
/* lattice, position, and types are overwritten. num_atom is returned. */
int spg_find_primitive(double lattice[3][3], double position[][3],
                       int types[], const int num_atom, const double symprec)
{
  int i, j;
  Cell cell, primitive;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);

  /* find primitive cell */
  if (sym_get_multiplicity(&cell, symprec) > 1) {

    primitive = prm_get_primitive(&cell, symprec);
    mat_copy_matrix_d3(lattice, primitive.lattice);

    for (i=0; i<primitive.size; i++) {
      types[i] = primitive.types[i];
            
      for (j=0; j<3; j++)
	position[i][j] = primitive.position[i][j];
    }
        
    cel_delete_cell(&primitive);
  }

  cel_delete_cell(&cell);
    
  return primitive.size;
}

/* Print-out space and point groups. This may be useful for */
/* testing, tasting, or debuging. */
void spg_show_symmetry(const double lattice[3][3], const double position[][3],
		       const int types[], const int num_atom, const double symprec)
{
  Cell cell;
  Spacegroup spacegroup;
  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);

  spacegroup = tbl_get_spacegroup(&cell, symprec);

  if (spacegroup.number) {
    printf("Space group No.%d\n", spacegroup.number);
    printf(" International: %s%s\n", spacegroup.bravais_symbol,
	   spacegroup.international);
    printf(" International(long): %s%s\n", spacegroup.bravais_symbol,
	   spacegroup.international_long);
    printf(" Schoenflies: %s\n", spacegroup.schoenflies);
    printf(" Multiplicity: %d\n", spacegroup.multi);
    printf("Point group\n");
    printf(" International: %s\n", spacegroup.pointgroup.international);
    printf(" Schoenflies: %s\n", spacegroup.pointgroup.schoenflies);
  }
    
  cel_delete_cell(&cell);
}

/* Return space group information in international table and number. */
int spg_get_international(char symbol[21], const double lattice[3][3],
			  const double position[][3],
			  const int types[], const int num_atom,
			  const double symprec)
{
  Cell cell;
  Spacegroup spacegroup;
  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);

  spacegroup = tbl_get_spacegroup(&cell, symprec);
  strcpy(symbol, spacegroup.bravais_symbol);
  strcpy(&symbol[1], spacegroup.international);
  
  cel_delete_cell(&cell);

  return spacegroup.number;
}
    
/* Return space group information in schoenflies and number. */
int spg_get_schoenflies(char symbol[10], const double lattice[3][3],
                        const double position[][3],
                        const int types[], const int num_atom,
			const double symprec)
{
  Cell cell;
  Spacegroup spacegroup;
  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);

  spacegroup = tbl_get_spacegroup(&cell, symprec);
  strcpy(symbol, spacegroup.schoenflies);

  return spacegroup.number;
}
    

