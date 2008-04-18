/* spglib.h */
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

#ifndef __spglib_H__
#define __spglib_H__

int spg_get_symmery(int rotation[][3][3], double translation[][3],
                    int max_size, double lattice[3][3],
                    double position[][3], int types[], int num_atom,
                    double symprec);
void spg_get_bravais_lattice(double bravais_lattice[3][3], double lattice[3][3],
                             double position[][3], int types[], int num_atom,
                             double symprec);
int spg_find_primitive(double lattice[3][3], double position[][3],
                       int types[], int num_atom, double symprec);
void spg_showsymmetry(double lattice[3][3], double position[][3],
                       int types[], int num_atom, double symprec);
int spg_get_international(char symbol[21], double lattice[3][3],
                          double position[][3], int types[],
                          int num_atom, double symprec);
int spg_get_schoenflies(char symbol[10], double lattice[3][3],
                        double position[][3], int types[], int num_atom,
                        double symprec);
void spg_show_symmetry_operations(double lattice[3][3], double position[][3],
                                  int types[], int num_atom, double symprec);

#endif
