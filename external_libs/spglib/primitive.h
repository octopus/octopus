/* primitive.h */
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

#ifndef __primitive_H__
#define __primitive_H__

#include "symmetry.h"
#include "cell.h"

int check_primitive_cell(Symmetry * symmetry);
int get_primitive_multiplicity(double pure_trans[][3],
                               const Symmetry * symmetry);
Cell get_primitive_cell(Cell * cell, Symmetry * symmetry, double symprec);
void get_primitive_least_axes(double vectors[][3], int size, Cell * cell,
                              double symprec);
void trim_cell(Cell * primitive, Cell * cell, double symprec);
int check_overlap(double a[3], double b[3], double symprec);

#endif
