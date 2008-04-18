/* cell.h */
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

#ifndef __cell_H__
#define __cell_H__

typedef struct {
    int size;
    double lattice[3][3];
    int *types;
    double (*position)[3];
} Cell;

Cell new_cell(int size);
void delete_cell(Cell * cell);
Cell copy_cell(Cell * cell_orig);
void get_cell_position(double position[][3], Cell * cell);
void get_cell_types(int types[], Cell * cell);
void set_cell(Cell * cell, double lattice[3][3], double position[][3],
              int types[]);

#endif
