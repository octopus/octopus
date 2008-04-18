/* cell.c */
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

#include <stdlib.h>
#include <stdio.h>
#include "cell.h"

Cell new_cell(int size)
{
    int i;
    Cell cell;
    cell.size = size;
    if ((cell.types = (int *) malloc(sizeof(int) * size)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    if ((cell.position =
         (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    return cell;
}

Cell copy_cell(Cell * cell_orig)
{
    Cell cell;

    cell = new_cell(cell_orig->size);
    set_cell(&cell, cell_orig->lattice, cell_orig->position,
             cell_orig->types);

    return cell;
}

void delete_cell(Cell * cell)
{
    free(cell->position);
    free(cell->types);
}

void get_cell_position(double position[][3], Cell * cell)
{
    int i;
    for (i = 0; i < cell->size; i++) {
        position[i][0] = cell->position[i][0];
        position[i][1] = cell->position[i][1];
        position[i][2] = cell->position[i][2];
    }
}

void get_cell_types(int types[], Cell * cell)
{
    int i;
    for (i = 0; i < cell->size; i++)
        types[i] = *(cell->types + i);
}

void set_cell(Cell * cell, double lattice[3][3], double position[][3],
              int types[])
{
    int i, j;
    copy_matrix_d3(cell->lattice, lattice);
    for (i = 0; i < cell->size; i++) {
        for (j = 0; j < 3; j++)
            cell->position[i][j] = position[i][j];
        cell->types[i] = types[i];
    }
}

