/* symmetry.c */
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
#include "cell.h"
#include "symmetry.h"
#include "mathfunc.h"
#include "bravais.h"

Symmetry get_symmetry_operation(Cell * cell, LatticeSymmetry * lattice_sym,
                                double symprec)
{
    double test_trans[3], tmp_symmetry[3][3];
    double tmp_vector[3],
        trans[cell->size * lattice_sym->size][3];
    int rot[cell->size * lattice_sym->size][3][3];
    double tmp_val;
    int i, j, k, l, m, count, num_sym = 0;
    Symmetry symmetry;

    for (i = 0; i < lattice_sym->size; i++)

        for (j = 0; j < cell->size; j++) {	/* test translation */

            cast_matrix_3i_to_3d(tmp_symmetry, lattice_sym->rot[i]);
            multiply_matrix_vector_d3(tmp_vector, tmp_symmetry,
                                      cell->position[0]);

            for (k = 0; k < 3; k++)
                test_trans[k] = cell->position[j][k] - tmp_vector[k];

            count = 0;

            for (k = 0; k < cell->size; k++) {	/* test nonsymmorphic operation for an atom */
                multiply_matrix_vector_d3(tmp_vector, tmp_symmetry,
                                          cell->position[k]);

                for (l = 0; l < cell->size; l++) {	/* check overlap of atom_k and atom_l */

                    for (m = 0; m < 3; m++) {	/* pos_l = S*pos_k + test_translation ?  */

                        tmp_val =
                            cell->position[l][m] - tmp_vector[m] -
                            test_trans[m];
                        tmp_val = Dabs(tmp_val - Nint(tmp_val));

                        if (tmp_val > symprec)
                            break;
                    }

                    if (tmp_val < symprec
                        && cell->types[k] == cell->types[l]) {
                        count++;
                        break;
                    }
                }               /* loop l */

                if (count < k + 1)
                    break;
            }                   /* loop k */

            if (count == cell->size) {	/* all atoms OK ? */
                copy_matrix_i3(rot[num_sym], lattice_sym->rot[i]);
                for (k = 0; k < 3; k++)
                    trans[num_sym][k] =
                        test_trans[k] - (double) Nint(test_trans[k] -
                                                      symprec);
                num_sym++;
            }
        }                       /* loop j */


    /* New a symmetry object */
    symmetry = new_symmetry(num_sym);

    for (i = 0; i < num_sym; i++) {
        copy_matrix_i3(symmetry.rot[i], rot[i]);

        for (j = 0; j < 3; j++)
            symmetry.trans[i][j] = trans[i][j];
    }

#ifdef DEBUG
    debug_print("*** get_symmetry_operation ***\n");

    for (i = 0; i < symmetry.size; i++) {

        debug_print("--- %d ---\n", i + 1);
        debug_print_matrix_i3(symmetry.rot[i]);
        debug_print("%f %f %f\n", symmetry.trans[i][0],
                    symmetry.trans[i][1], symmetry.trans[i][2]);
    }
#endif



    return symmetry;
}


Symmetry new_symmetry(int size)
{
    Symmetry symmetry;
    int i;
    symmetry.size = size;
    if ((symmetry.rot =
         (int (*)[3][3]) malloc(sizeof(int[3][3]) * size)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    if ((symmetry.trans =
         (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
        fprintf(stderr, "Memory could not be allocated.");
        exit(1);
    }
    return symmetry;
}

void delete_symmetry(Symmetry * symmetry)
{
    free(symmetry->rot);
    free(symmetry->trans);
}
