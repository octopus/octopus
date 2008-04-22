/* primitive.c */
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

/* #define DEBUG */

#include <stdio.h>
#include <stdlib.h>
#include "primitive.h"
#include "mathfunc.h"
#include "cell.h"
#include "symmetry.h"

int check_primitive_cell(Symmetry * symmetry)
{
    double pure_trans[symmetry->size][3];
    int multi;
    multi = get_primitive_multiplicity(pure_trans, symmetry);
    if (multi > 0)
        return multi;
    else {
        fprintf(stderr, "No identity operation.\n");
        exit(1);
    }
}

Cell get_primitive_cell(Cell * cell, Symmetry * symmetry, double symprec)
{
    int i, j, multi;
    double pure_trans[symmetry->size][3], new_lattice[3][3];
    Cell primitive;

    debug_print("*** get_primitive_cell ***\n");

    multi = get_primitive_multiplicity(pure_trans, symmetry);
    primitive = new_cell(cell->size / multi);

    if (multi > 1) {
        double vectors[multi + 2][3];

        for (i = 0; i < multi - 1; i++)
            copy_vector_d3(vectors[i], pure_trans[i + 1]);

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                if (i == j)
                    vectors[i+multi-1][j] = 1;
                else
                    vectors[i+multi-1][j] = 0;

        for (i = 0; i < multi + 2; i++)
            debug_print("%d: %f %f %f\n", i + 1, vectors[i][0],
                        vectors[i][1], vectors[i][2]);

        /* vectors are reduced axis about orignal lattice. */
        get_primitive_least_axes(vectors, multi + 2, cell, symprec);


        for (i = 0; i < 3; i++) {
            debug_print("found axis %d: %f %f %f\n", i + 1, vectors[i][0],
                        vectors[i][1], vectors[i][2]);
            multiply_matrix_vector_d3(new_lattice[i], cell->lattice,
                                      vectors[i]);
        }

        transpose_matrix_d3(new_lattice, new_lattice);

        debug_print("Original cell lattice.\n");
        debug_print_matrix_d3(cell->lattice);

        debug_print
            ("Found primitive lattice after choosing least axes.\n");
        debug_print_matrix_d3(new_lattice);

        get_smallest_primitive(primitive.lattice, new_lattice, symprec);

        debug_print
            ("Found primitive lattice after choosing smallest axes.\n");
        debug_print_matrix_d3(primitive.lattice);
    }

    else {                      /* multi == 1: No primitive found */
        printf("Primitive cell was not found.\n");
        copy_matrix_d3(primitive.lattice, cell->lattice);
    }

    debug_print("Number of atoms in primitive cell: %d\n", primitive.size);
    debug_print("Found primitive lattice\n");
    debug_print_matrix_d3(primitive.lattice);
    trim_cell(&primitive, cell, symprec);

    return primitive;
}

void trim_cell(Cell * primitive, Cell * cell, double symprec)
{
    int i, j, k, count, ratio;
    int table[cell->size][cell->size], check_table[cell->size];
    double axis_inv[3][3], tmp_matrix[3][3], position[cell->size][3];

    debug_print("*** trim_cell ***\n");

    primitive->size = Nint(Dabs(get_determinant_d3(primitive->lattice)) /
                           Dabs(get_determinant_d3(cell->lattice)) *
                           cell->size);

    ratio = cell->size / primitive->size;

    inverse_matrix_d3(tmp_matrix, primitive->lattice, symprec * symprec * symprec);
    multiply_matrix_d3(axis_inv, tmp_matrix, cell->lattice);

    /* Send atoms into the primitive cell */
    debug_print("Reduced new position in new axes\n");
    for (i = 0; i < cell->size; i++) {
        multiply_matrix_vector_d3(position[i], axis_inv,
                                  cell->position[i]);
        for (j = 0; j < 3; j++)
            position[i][j] = position[i][j] - Nint(position[i][j]);
        debug_print("%d: %f %f %f\n", i + 1, position[i][0],
                    position[i][1], position[i][2]);
    }

    /* Create overlapping table */
    for (i = 0; i < cell->size; i++)
        for (j = 0; j < cell->size; j++)
            table[i][j] = 0;

    for (i = 0; i < cell->size; i++) {
        count = 0;
        for (j = 0; j < cell->size; j++) {
            if (check_overlap(position[i], position[j], symprec)) {
                table[i][count] = j;
                count++;
            }
        }
        if (count != ratio) {
            fprintf(stderr, "Bug: Primitive cell could not found.\n");
            exit(0);
        }
    }

#ifdef DEBUG
    debug_print("Overlaping table\n");
    for (i = 0; i < cell->size; i++) {
        debug_print("%2d: ", count);
        for (j = 0; j < cell->size; j++)
            debug_print("%2d ", table[i][j]);
        debug_print("\n");
    }
#endif

    /* Copy positions. Positions of overlapped atoms are averaged. */
    for (i = 0; i < cell->size; i++)
        check_table[i] = 0;
    count = 0;

    for (i = 0; i < cell->size; i++)

        if (!check_table[i]) {
            primitive->types[count] = cell->types[i];

            for (j = 0; j < 3; j++)
                primitive->position[count][j] = 0;

            for (j = 0; j < ratio; j++) {	/* overlap atoms */

                for (k = 0; k < 3; k++)

                    /* boundary treatment */
                    if (Dabs(position[table[i][0]][k] -
                             position[table[i][j]][k]) > 0.5)

                        if (position[table[i][j]][k] < 0)
                            primitive->position[count][k]
                                = primitive->position[count][k] +
                                position[table[i][j]][k] + 1;

                        else
                            primitive->position[count][k]
                                = primitive->position[count][k] +
                                position[table[i][j]][k] - 1;

                    else
                        primitive->position[count][k]
                            = primitive->position[count][k] +
                            position[table[i][j]][k];
                check_table[table[i][j]] = 1;
            }

            for (j = 0; j < 3; j++) {	/* take average and reduce */

                primitive->position[count][j] =
                    primitive->position[count][j] / ratio;

                primitive->position[count][j] =
                    primitive->position[count][j] -
                    Nint(primitive->position[count][j] - symprec);
            }
            count++;
        }

    if (count != primitive->size) {
        fprintf(stderr, "Bug: Primitive cell could not found.\n");
        exit(0);
    }

    debug_print("Trimed position\n");
    debug_print_vectors_with_label(primitive->position, primitive->types,
                                   primitive->size);

}

int check_overlap(double a[3], double b[3], double symprec)
{
    if ((Dabs(a[0] - b[0]) < symprec
         || Dabs(Dabs(a[0] - b[0]) - 1) < symprec)
        && (Dabs(a[1] - b[1]) < symprec
            || Dabs(Dabs(a[1] - b[1]) - 1) < symprec)
        && (Dabs(a[2] - b[2]) < symprec
            || Dabs(Dabs(a[2] - b[2]) - 1) < symprec))
        return 1;
    return 0;
}

void get_primitive_least_axes(double vectors[][3], int size, Cell * cell,
                              double symprec)
{
    int i, j, k;
    double min_volume, initial_volume, volume, min_vectors[3][3],
        tmp_lattice[3][3];

    debug_print("*** get_primitive_least_axes ***\n");

    initial_volume = Dabs(get_determinant_d3(cell->lattice));
    min_volume = initial_volume;
    debug_print("initial volume: %f\n", initial_volume);

    /* check volumes of all possible lattices, find smallest volume */
    for (i = 0; i < size; i++)
        for (j = i + 1; j < size; j++)
            for (k = j + 1; k < size; k++) {
                multiply_matrix_vector_d3(tmp_lattice[0], cell->lattice,
                                          vectors[i]);
                multiply_matrix_vector_d3(tmp_lattice[1], cell->lattice,
                                          vectors[j]);
                multiply_matrix_vector_d3(tmp_lattice[2], cell->lattice,
                                          vectors[k]);
                volume = Dabs(get_determinant_d3(tmp_lattice));
                if (volume < min_volume - symprec && volume > symprec) {
                    min_volume = volume;
                    copy_vector_d3(min_vectors[0], vectors[i]);
                    copy_vector_d3(min_vectors[1], vectors[j]);
                    copy_vector_d3(min_vectors[2], vectors[k]);
                }
            }

    if (min_volume < initial_volume - symprec) {
        debug_print("minimum volume: %f\n", min_volume);
        for (i = 0; i < 3; i++)
            copy_vector_d3(vectors[i], min_vectors[i]);
    } else {
        /* This should not happen. */
        fprintf(stderr, "No primitive least axes are found.");
        exit(1);
    }
}


int get_primitive_multiplicity(double pure_trans[][3],
                               const Symmetry * symmetry)
{
    int i, j, multi = 0;

    debug_print("*** get_primitive_multiplicity ***\n");

    for (i = 0; i < symmetry->size; i++)
        if (symmetry->rot[i][0][0] == 1 &&
            !symmetry->rot[i][0][1] &&
            !symmetry->rot[i][0][2] &&
            !symmetry->rot[i][1][0] &&
            symmetry->rot[i][1][1] == 1 &&
            !symmetry->rot[i][1][2] &&
            !symmetry->rot[i][2][0] &&
            !symmetry->rot[i][2][1] && symmetry->rot[i][2][2] == 1) {
            for (j = 0; j < 3; j++)
                pure_trans[multi][j] = symmetry->trans[i][j];
            multi++;
            debug_print("%d: %f %f %f\n", i + 1, pure_trans[i][0],
                        pure_trans[i][1], pure_trans[i][2]);
        }
    return multi;
}
