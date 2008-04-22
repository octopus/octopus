/* bravais.c */
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
#include "bravais.h"
#include "mathfunc.h"

LatticeSymmetry get_symmetry_candidate(Bravais * bravais,
                                       double lattice[3][3],
                                       double symprec)
{
    int i, j;
    double coordinate[3][3], tmp_matrix_d3[3][3], check_val;
    LatticeSymmetry lattice_sym;

    debug_print("*** get_symmetry_candidate ***\n");

    /* Obtain ratio of lattice and bravais lattice. B^-1*P */
    inverse_matrix_d3(tmp_matrix_d3, bravais->lattice, symprec * symprec * symprec);
    multiply_matrix_d3(coordinate, tmp_matrix_d3, lattice);
    debug_print("Ratio of lattice and bravais lattice. B^-1*P.\n");
    debug_print_matrix_d3(coordinate);

    /* check coordinate: elements of coordinate have to be integer or half intenger */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            check_val = coordinate[i][j] * bravais->factor;
    if (Dabs(check_val - (double) Nint(check_val)) > symprec) {
        fprintf(stderr, "Bug\n");
    }

    bravais->lattice_symmetry = get_bravais_symmetry(bravais->holohedry);

    for (i = 0; i < bravais->lattice_symmetry.size; i++) {

        cast_matrix_3i_to_3d(tmp_matrix_d3,
                             bravais->lattice_symmetry.rot[i]);

        if (!get_similar_matrix_d3(tmp_matrix_d3, tmp_matrix_d3,
                                   coordinate, symprec)) {

            debug_print_matrix_d3(tmp_matrix_d3);
            fprintf(stderr, "BUG: Determinant is zero.\n");
            exit(1);
        }

        cast_matrix_3d_to_3i(lattice_sym.rot[i], tmp_matrix_d3);
        lattice_sym.size = bravais->lattice_symmetry.size;
    }

    return lattice_sym;
}

LatticeSymmetry get_bravais_symmetry(int holohedry)
{
    int i, j, k;
    LatticeSymmetry lattice_sym;
    int identity[3][3] = {
        {1, 0, 0},              /* order 1 */
        {0, 1, 0},
        {0, 0, 1}
    };

    int inversion[3][3] = {
        {-1, 0, 0},             /* order 2 */
        {0, -1, 0},
        {0, 0, -1}
    };

    int generator6[3][3] = {
        {1, -1, 0},             /* order 6 */
        {1, 0, 0},
        {0, 0, 1}
    };

    int generator3[3][3] = {
        {0, 1, 0},              /* order 3 */
        {0, 0, 1},
        {1, 0, 0}
    };

    int generator2xy[3][3] = {
        {0, 1, 0},              /* order 2 */
        {1, 0, 0},
        {0, 0, 1}
    };

    int generator2y[3][3] = {
        {-1, 0, 0},             /* order 2 */
        {0, 1, 0},
        {0, 0, -1}
    };

    int generator2z[3][3] = {
        {-1, 0, 0},             /* order 2 */
        {0, -1, 0},
        {0, 0, 1}
    };

    debug_print("*** get_bravais_symmetry ***\n");

    /* all clear */
    for (i = 0; i < 48; i++)
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                lattice_sym.rot[i][j][k] = 0;
    lattice_sym.size = 0;

    /* indentity: this is seed. */
    copy_matrix_i3(lattice_sym.rot[0], identity);

    /* inversion */
    generate_point_symmetry(lattice_sym.rot, inversion, 1, 2);

    switch (holohedry) {
    case CUBIC:
        generate_point_symmetry(lattice_sym.rot, generator2y, 2, 2);
        generate_point_symmetry(lattice_sym.rot, generator2z, 4, 2);
        generate_point_symmetry(lattice_sym.rot, generator2xy, 8, 2);
        generate_point_symmetry(lattice_sym.rot, generator3, 16, 3);
        lattice_sym.size = 48;
        break;
    case HEXA:
        generate_point_symmetry(lattice_sym.rot, generator2xy, 2, 2);
        generate_point_symmetry(lattice_sym.rot, generator6, 4, 6);
        lattice_sym.size = 24;
        break;
    case TRIGO:
        generate_point_symmetry(lattice_sym.rot, generator2xy, 2, 2);
        generate_point_symmetry(lattice_sym.rot, generator3, 4, 3);
        lattice_sym.size = 12;
        break;
    case TETRA:
        generate_point_symmetry(lattice_sym.rot, generator2y, 2, 2);
        generate_point_symmetry(lattice_sym.rot, generator2z, 4, 2);
        generate_point_symmetry(lattice_sym.rot, generator2xy, 8, 2);
        lattice_sym.size = 16;
        break;
    case ORTHO:
        generate_point_symmetry(lattice_sym.rot, generator2y, 2, 2);
        generate_point_symmetry(lattice_sym.rot, generator2z, 4, 2);
        lattice_sym.size = 8;
        break;
    case MONOCLI:
        generate_point_symmetry(lattice_sym.rot, generator2y, 2, 2);
        lattice_sym.size = 4;
        break;
    case TRICLI:
        lattice_sym.size = 2;
        break;
    default:
        fprintf(stderr, "BUG: no Bravais lattice found.\n");
    }

    return lattice_sym;
}

void generate_point_symmetry(int point_symmetry[][3][3],
                             int generator[3][3], int n_sym, int n_gen)
{
    static int i, j, count;  /* this is declared static to avoid
				inlining, which triggers a bug in
				gcc-4.3 */
    int tmp_matrix[3][3];

    /* n_sym is number of symmetry operations, which was already counted. */
    /* n_gen is order, number of operations of the generator class. */
    for (i = 0; i < n_sym; i++) {
        for (j = 0; j < n_gen - 1; j++) {	/* "-1" comes from E (identity) in a class */
            count = i * (n_gen - 1) + j + n_sym;
            multiply_matrix_i3(tmp_matrix, generator,
                               point_symmetry[count - n_sym]);
            copy_matrix_i3(point_symmetry[count], tmp_matrix);
            debug_print("--- %d ---\n", count + 1);
            debug_print_matrix_i3(tmp_matrix);
        }
    }
}

Bravais get_bravais_lattice(double lattice_orig[3][3], double symprec)
{
    Bravais bravais;
    double min_lattice[3][3];
    int i, holohedry, found;

    get_smallest_primitive(min_lattice, lattice_orig, symprec);

    debug_print("*** get_bravais_lattice ***\n");
    debug_print("Original lattice\n");
    debug_print_matrix_d3(lattice_orig);

    /* loop: cubic, hexa, trigo, tetra, ortho, monocli, tricli */
    for (holohedry = CUBIC; holohedry > 0; holohedry--) {
        bravais.holohedry = holohedry;

        for (i = 0; i < 3; i++) {	/* look for principal axis */

            found = get_bravais_lattice_in_loop(&bravais, min_lattice, i,
                                                symprec);
            if (found)
                break;
        }

        if (found)
            break;
    }
    
    debug_print("Bravais lattice\n");
    debug_print_matrix_d3(bravais.lattice);
    debug_print("holohedry: %d, centering: %d, factor: %d, volume: %f\n",
                bravais.holohedry, bravais.centering, bravais.factor,
                get_determinant_d3(bravais.lattice));

    set_bravais_lattice(&bravais, symprec);

    debug_print("--- return to get_bravais_lattice ---\n");
    debug_print("centering: %d, factor: %d, holohedry: %d\n",
                bravais.centering, bravais.factor, bravais.holohedry);
    debug_print("Bravais lattice after post-checking\n");
    debug_print_matrix_d3(bravais.lattice);

    return bravais;
}

int get_bravais_lattice_in_loop(Bravais *bravais, double min_lattice[3][3],
                                int axis, double symprec)
{
    int i, found, holohedry, orthogonal, centering, factor, c0, c1;
    double lattice[3][3];
    int combination[3][2] = { {1, 2}, {2, 0}, {0, 1} };	/* combination of axes for iteration */
    int angle_90[3], edge_equal[3];
    
    found = 0;
    orthogonal = 0;
    holohedry = bravais->holohedry;
    check_angle90_equaledge(angle_90, edge_equal, min_lattice, symprec);

    debug_print("angle_90: %d %d %d\n", angle_90[0], angle_90[1],
                angle_90[2]);
    debug_print("equal:    %d %d %d\n", edge_equal[0], edge_equal[1],
                edge_equal[2]);

    if (holohedry == CUBIC || holohedry == TETRA || holohedry == ORTHO)
        orthogonal = 1;

    c0 = combination[axis][0];
    c1 = combination[axis][1];

    /* cubic, tetragonal, orthorhombic */
    if (angle_90[0] && angle_90[1] && angle_90[2] && orthogonal) {

        copy_matrix_d3(lattice, min_lattice);
        found = check_holohedry(lattice, holohedry, symprec);

        if (found) {
            factor = 1;
            centering = NO_CENTER;
            debug_print("ortho, centering:%d\n", centering);
        }
    }

    /* hexagonal */
    if (!found && holohedry == HEXA && angle_90[c0] && angle_90[c1]
        && edge_equal[axis]) {
        if (bravais_hexa(axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            factor = 1;
            centering = NO_CENTER;
            debug_print("hexa, centering:%d\n", centering);
        }
    }

    /* orthogonal1: A-, B-, C- or body-centered, a and b are orthogonal. */
    if (!found && orthogonal && angle_90[axis]) {
        centering =
            bravais_ortho1(axis, c0, c1, min_lattice, lattice,
                           symprec);
        if (centering)
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("ortho1, centering:%d\n", centering);
            factor = 2;
        }
    }

    /* orthogonal2: A-, B-, or C-centeredm, a=b, a and c, b and c are orthogonal. */
    if (!found && holohedry == ORTHO && angle_90[c0]
        && angle_90[c1] && edge_equal[axis]) {
        factor = 2;
        centering = axis + 1;  /* C-centered */
        for (i = 0; i < 3; i++) {	/* Right-handed ? */
            lattice[i][c0] =
                min_lattice[i][c0] - min_lattice[i][c1];
            lattice[i][c1] =
                min_lattice[i][c0] + min_lattice[i][c1];
            lattice[i][axis] = min_lattice[i][axis];
        }
        found = check_holohedry(lattice, holohedry, symprec);
        debug_print("ortho2, centering:%d\n", centering);
    }

    /* orthogonal3: body- or face-centered, a, b is not axes of Bravais lattice. */
    if (!found && orthogonal) {
        centering =
            bravais_ortho3(axis, c0, c1, min_lattice, lattice,
                           symprec);
        if (centering)
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            factor = 2;
            debug_print("ortho3, centering:%d\n", centering);
        }
    }

    /* orthogonal4: body-centered, a=b=c, a,b,c are not axes of Bravais lattice */
    if (!found && orthogonal && edge_equal[0] && edge_equal[1]
        && edge_equal[2]) {
        if (bravais_ortho4
            (axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            centering = BODY;
            factor = 2;
            debug_print("ortho4, centering:%d\n", centering);
        }
    }

    /* orthogonal5: face-centered, a=b, a,b,c are not axes of Bravais lattice. */
    if (!found && orthogonal && edge_equal[axis]) {
        if (bravais_ortho5
            (axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            centering = FACE;
            factor = 2;
            debug_print("ortho5, centering:%d\n", centering);
        }
    }

    /* orthogonal6: face-centered, a!=b!=c, a,b,c are not axes of Bravais lattice. */
    if (!found && orthogonal) {
        if (bravais_ortho6
            (axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            centering = FACE;
            factor = 2;
            debug_print("ortho6, centering:%d\n", centering);
        }
    }

    /* rhombohedral1: a=b=c, |cos(angle_a)| = |cos(angle_b)| = |cos(angle_c)| */
    if (!found && holohedry == TRIGO && edge_equal[0]
        && edge_equal[1] && edge_equal[2]) {
        if (bravais_rhombo1
            (axis, combination, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("rhombo1\n");
            factor = 1;
            centering = NO_CENTER;
        }
    }

    /* rhombohedral2: a=b, c is trigonal axis. */
    if (!found && holohedry == TRIGO && edge_equal[axis]) {
        if (bravais_rhombo2
            (axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("rhombo2\n");
            factor = 1;
            centering = NO_CENTER;
        }
    }

    /* rhombohedral3: a=b, c is in plane normal to trigonal axis. */
    if (!found && holohedry == TRIGO && edge_equal[axis]) {
        if (bravais_rhombo3
            (axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("rhombo3\n");
            factor = 1;
            centering = NO_CENTER;
        }
    }

    /* rhombohedral4: two vectors are in plane normal to traigonal axis. */
    if (!found && holohedry == TRIGO && edge_equal[axis]) {
        if (bravais_rhombo4
            (axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("rhombo4\n");
            factor = 1;
            centering = NO_CENTER;
        }
    }

    /* monoclinic1 */
    if (!found && holohedry == MONOCLI && angle_90[c0]
        && angle_90[c1]) {
        for (i = 0; i < 3; i++) {
            lattice[i][0] = min_lattice[i][c1];
            lattice[i][1] = min_lattice[i][axis];
            lattice[i][2] = min_lattice[i][c0];
        }
        found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("monocli1\n");
            factor = 1;
            centering = NO_CENTER;
        }
    }

    /* monoclinic2: one-face-centered, two vecs face center. */
    if (!found && holohedry == MONOCLI && edge_equal[axis]) {
        if (bravais_monocli2
            (axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("monocli2\n");
            factor = 2;
            centering = C_FACE;
        }
    }

    /* monoclinic3: one-face-centered, one vec face center. */
    if (!found && holohedry == MONOCLI) {
        if (bravais_monocli3
            (axis, c0, c1, min_lattice, lattice, symprec))
            found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("monocli3\n");
            factor = 2;
            centering = C_FACE;
        }
    }

    /* triclinic */
    if (!found && holohedry == TRICLI) {
        copy_matrix_d3(lattice, min_lattice);
        found = check_holohedry(lattice, holohedry, symprec);
        if (found) {
            debug_print("triclinic\n");
            factor = 1;
            centering = NO_CENTER;
        }
    }
        
    if (found) {
        copy_matrix_d3(bravais->lattice, lattice);
        bravais->centering = centering;
        bravais->factor = factor;
    }

    return found;

}

void set_bravais_lattice(Bravais *bravais, double symprec)
{
    double metric[3][3], lattice[3][3];
    int angle_90[3], edge_equal[3];
    int i, centering, holohedry, found = 0;
    
    debug_print("*** set_bravais_lattice ***\n");

    copy_matrix_d3(lattice, bravais->lattice);
    centering = bravais->centering;
    holohedry = bravais->holohedry;
    get_metric(metric, lattice);

    check_angle90_equaledge(angle_90, edge_equal, lattice, symprec);

    debug_print("angle_90: %d %d %d\n", angle_90[0], angle_90[1],
                angle_90[2]);
    debug_print("equal:    %d %d %d\n", edge_equal[0], edge_equal[1],
                edge_equal[2]);
    debug_print("centering: %d\n", centering);

    /* orthorhombic, tetragonal, cubic */
    if (angle_90[0] && angle_90[1] && angle_90[2]) {

        /* cubic */
        if (edge_equal[0] && edge_equal[1] && edge_equal[2] &&
            holohedry==CUBIC) {
            /* non-, face-, and body-centered */
            if (centering == NO_CENTER || centering == BODY
                || centering == FACE) {
                holohedry = CUBIC;
                found = 1;
                copy_matrix_d3(bravais->lattice, lattice);
            }
        }

        /* tetragonal */
        if (!found && (holohedry == TETRA || holohedry == CUBIC) && 
            (edge_equal[0] || edge_equal[1] || edge_equal[2])) {
            if (centering == NO_CENTER || centering == BODY) {	/* non- and body-centered */
                holohedry = TETRA;
                found = 1;
                if (edge_equal[0])
                    for (i = 0; i < 3; i++) {
                        bravais->lattice[i][0] = lattice[i][1];
                        bravais->lattice[i][1] = lattice[i][2];
                        bravais->lattice[i][2] = lattice[i][0];
                    }
                if (edge_equal[1])
                    for (i = 0; i < 3; i++) {
                        bravais->lattice[i][0] = lattice[i][2];
                        bravais->lattice[i][1] = lattice[i][0];
                        bravais->lattice[i][2] = lattice[i][1];
                    }
                if (edge_equal[2])
                    copy_matrix_d3(bravais->lattice, lattice);
            }
        }

        /* orthorhombic */
        if (!found) {
            holohedry = ORTHO;
            found = 1;
            copy_matrix_d3(bravais->lattice, lattice);
        }
    }

    /* hexagonal */
    if (!found && angle_90[0] && angle_90[1] && edge_equal[2] &&
        Dabs(2 * metric[1][0] + metric[1][1]) < symprec) {
        holohedry = HEXA;
        found = 1;
        copy_matrix_d3(bravais->lattice, lattice);
    }

    /* rhombohedral */
    if (!found && edge_equal[0] && edge_equal[1] && edge_equal[2] &&
        Dabs(metric[1][0] - metric[2][1]) < symprec &&
        Dabs(metric[1][0] - metric[2][0]) < symprec) {
        holohedry = TRIGO;
        found = 1;
        copy_matrix_d3(bravais->lattice, lattice);
    }

    /* monoclinic */
    if (!found && ( (angle_90[0] && angle_90[1]) ||
                    (angle_90[1] && angle_90[2]) ||
		    (angle_90[2] && angle_90[0]) )) {
        if (angle_90[0] && angle_90[1])
            for (i = 0; i < 3; i++) {
                bravais->lattice[i][0] = lattice[i][0];
                bravais->lattice[i][1] = -lattice[i][2];
                bravais->lattice[i][2] = lattice[i][1];
            }
        if (angle_90[1] && angle_90[2])
            for (i = 0; i < 3; i++) {
                bravais->lattice[i][0] = lattice[i][1];
                bravais->lattice[i][1] = -lattice[i][0];
                bravais->lattice[i][2] = lattice[i][2];
            }
        if (angle_90[2] && angle_90[0])
            copy_matrix_d3(bravais->lattice, lattice);
        holohedry = MONOCLI;
        found = 1;
    }

    /* triclinic */
    if (!found) {
        holohedry = TRICLI;
        found = 1;
        copy_matrix_d3(bravais->lattice, lattice);
    }

    if (get_determinant_d3(bravais->lattice) < symprec)
        for (i = 0; i < 3; i++) {
            bravais->lattice[i][0] = -lattice[i][1];
            bravais->lattice[i][1] = -lattice[i][0];
            bravais->lattice[i][2] = -lattice[i][2];
        }

    bravais->holohedry = holohedry;
}


int bravais_hexa(int axis, int c0, int c1, double min_lattice[3][3],
                 double lattice[3][3], double symprec)
{
    int j;
    double ratio, metric[3][3];
    get_metric(metric, min_lattice);

    ratio = metric[c1][c0] / metric[c0][c0];	/* ab/aa, |b|cos/|a|=cos */
    if (Dabs(ratio + 0.5) < symprec) {	/* theta = 120 deg */
        for (j = 0; j < 3; j++) {
            lattice[j][0] = min_lattice[j][c0];
            lattice[j][1] = min_lattice[j][c1];
            lattice[j][2] = min_lattice[j][axis];
        }
        return 1;
    }
    if (Dabs(ratio - 0.5) < symprec) {	/* theta = 60 deg */
        for (j = 0; j < 3; j++) {
            lattice[j][0] = min_lattice[j][c0];
            lattice[j][1] = min_lattice[j][c1] - min_lattice[j][c0];
            lattice[j][2] = min_lattice[j][axis];
        }
        return 1;
    }
    return 0;
}

int bravais_ortho1(int axis, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec)
{
    int j, centering = NO_CENTER;
    double ratio_a, ratio_b, metric[3][3];
    get_metric(metric, min_lattice);

    ratio_a = metric[axis][c0] / metric[c0][c0];	/* ca/aa = |c|cos/|a| */
    ratio_b = metric[axis][c1] / metric[c1][c1];	/* cb/bb = |c|cos/|b| */

    for (j = 0; j < 3; j++) {
        lattice[j][c0] = min_lattice[j][c0];
        lattice[j][c1] = min_lattice[j][c1];
    }

    if (Dabs(Dabs(ratio_a) - 0.5) < symprec &&	/* |a|=2*|c||cos|; B centered */
        Dabs(ratio_b) < symprec)	/* theta (b and c) = 90 deg */
        centering = c1 + 1;     /* B-centered */

    if (Dabs(ratio_a) < symprec &&	/* theta (a and c) = 90 deg */
        Dabs(Dabs(ratio_b) - 0.5) < symprec)	/* |b|=2*|c||cos|; A centered */
        centering = c0 + 1;     /* A-centered */

    if (Dabs(Dabs(ratio_a) - 0.5) < symprec &&	/* |a|=2*|c||cos|; B centered */
        Dabs(Dabs(ratio_b) - 0.5) < symprec)	/* |b|=2*|c||cos|; A centered */
        centering = BODY;       /* body-centered */

    if (centering) {
        for (j = 0; j < 3; j++)
            lattice[j][axis] =
                (min_lattice[j][axis] - ratio_a * min_lattice[j][c0]
                 - ratio_b * min_lattice[j][c1]) * 2;
    }

    return centering;
}


int bravais_ortho3(int axis, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec)
{
    int j, centering = NO_CENTER;
    double vector_a[3], vector_b[3], ratio_a, ratio_b, metric[3][3];
    get_metric(metric, min_lattice);

    ratio_a = metric[axis][c0] / metric[axis][axis];	/* ca/cc = |a|cos/|c| */
    ratio_b = metric[axis][c1] / metric[axis][axis];	/* cb/cc = |b|cos/|c| */
    if (Dabs(Dabs(ratio_a) - 0.5) < symprec &&	/* |c|=2*|a||cos| */
        Dabs(Dabs(ratio_b) - 0.5) < symprec) {	/* |c|=2*|b||cos| */

        for (j = 0; j < 3; j++) {	/* projection to a-b surface of Bravais lattice */
            vector_a[j] = min_lattice[j][c0] - ratio_a * min_lattice[j][axis];
            vector_b[j] = min_lattice[j][c1] - ratio_b * min_lattice[j][axis];
        }

        /* |a|=|b|: body-centered */
        if (Dabs(inner_product_d3(vector_a, vector_a) -
                 inner_product_d3(vector_b, vector_b)) < symprec) {
            centering = BODY;   /* body-centered */
            for (j = 0; j < 3; j++) {	/* Right-handed ? */
                lattice[j][c0] = vector_a[j] - vector_b[j];
                lattice[j][c1] = vector_a[j] + vector_b[j];
                lattice[j][axis] = min_lattice[j][axis];
            }
        } else
            /* a and b is orthogonal: face-centered */
        if (Dabs(inner_product_d3(vector_a, vector_b)) < symprec) {
            centering = FACE;   /* face-centered */
            for (j = 0; j < 3; j++) {
                lattice[j][c0] = 2 * vector_a[j];
                lattice[j][c1] = 2 * vector_b[j];
                lattice[j][axis] = min_lattice[j][axis];
            }
        }
    }
    return centering;
}

int bravais_ortho4(int axis, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec)
{
    int j;
    double vector_a[3], vector_b[3], ratio_a, ratio_b, tmp_matrix_d3[3][3];

    for (j = 0; j < 3; j++) {
        vector_a[j] = min_lattice[j][c0] - min_lattice[j][c1];	/* Right handed ? */
        vector_b[j] = min_lattice[j][c0] + min_lattice[j][c1];
    }
    /* project on to an another axis */
    transpose_matrix_d3(tmp_matrix_d3, min_lattice);
    ratio_a = inner_product_d3(tmp_matrix_d3[axis], vector_a) / inner_product_d3(vector_a, vector_a);	/* ac/aa = |c|cos/|a| */
    ratio_b = inner_product_d3(tmp_matrix_d3[axis], vector_b) / inner_product_d3(vector_b, vector_b);	/* bc/bb = |c|cos/|b| */
    /* one of the vectors a and b can be an axis of Bravais lattice */
    if (Dabs(Dabs(ratio_a) - 0.5) < symprec) {	/* |a|=2*|c||cos| */
        for (j = 0; j < 3; j++) {
            lattice[j][c0] = vector_a[j];
            vector_a[j] = min_lattice[j][axis] - ratio_a * vector_a[j];
            vector_b[j] = 0.5 * vector_b[j];
            lattice[j][c1] = vector_a[j] + vector_b[j];
            lattice[j][axis] = vector_a[j] - vector_b[j];
        }
        return 1;
    }
    if (Dabs(Dabs(ratio_b) - 0.5) < symprec) {	/* |a|=2*|c||cos| */
        for (j = 0; j < 3; j++) {
            lattice[j][c1] = vector_b[j];
            vector_b[j] = min_lattice[j][axis] - ratio_b * vector_b[j];
            vector_a[j] = 0.5 * vector_a[j];
            lattice[j][c0] = vector_b[j] + vector_a[j];
            lattice[j][axis] = vector_b[j] - vector_a[j];
        }
        return 1;
    }
    return 0;
}

int bravais_ortho5(int axis, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec)
{
    int j;
    double vector_a[3], vector_b[3], ratio_a, ratio_b, tmp_matrix_d3[3][3];

    for (j = 0; j < 3; j++) {
        vector_a[j] = min_lattice[j][c0] - min_lattice[j][c1];	/* Right handed ? */
        vector_b[j] = min_lattice[j][c0] + min_lattice[j][c1];
    }
    /* project on to an another axis */
    transpose_matrix_d3(tmp_matrix_d3, min_lattice);

    ratio_a = inner_product_d3(tmp_matrix_d3[axis], vector_a) / inner_product_d3(vector_a, vector_a);	/* ac/aa = |c|cos/|a| */

    ratio_b = inner_product_d3(tmp_matrix_d3[axis], vector_b) / inner_product_d3(vector_b, vector_b);	/* bc/bb = |c|cos/|b| */

    if (Dabs(Dabs(ratio_a) - 0.5) < symprec &&	/* |a|=2*|c||cos| */
        (Dabs(ratio_b) < symprec ||  Dabs(ratio_a) < symprec ) &&  /* b and c are orthogonal or  a and c are orthogonal */
        Dabs(Dabs(ratio_b) - 0.5) < symprec) {	/* |b|=2*|c||cos| */
        for (j = 0; j < 3; j++) {
            lattice[j][axis] = (min_lattice[j][axis] - ratio_a * vector_a[j]
                             - ratio_b * vector_b[j]) * 2;
            lattice[j][c0] = vector_a[j];
            lattice[j][c1] = vector_b[j];
        }
        return 1;
    }
    return 0;
}

int bravais_ortho6(int axis, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec)
{
    int j;
    double vector_a[3], vector_b[3], metric[3][3];

    get_metric(metric, min_lattice);

    for (j = 0; j < 3; j++) {
        vector_a[j] = min_lattice[j][c0] - min_lattice[j][c1];	/* Right handed ? */
        vector_b[j] = min_lattice[j][c0] + min_lattice[j][c1];
    }

    if (Dabs(metric[axis][axis] - inner_product_d3(vector_a, vector_a)) <
        symprec)
        return bravais_ortho6_ext(axis, c0, c1, min_lattice, lattice,
                                  vector_a, vector_b, symprec);

    if (Dabs(metric[axis][axis] - inner_product_d3(vector_b, vector_b)) <
        symprec)
        return bravais_ortho6_ext(axis, c0, c1, min_lattice, lattice,
                                  vector_b, vector_a, symprec);

    return 0;
}

int bravais_ortho6_ext(int axis, int c0, int c1, double min_lattice[3][3],
                       double lattice[3][3], double vector_a[3],
                       double vector_b[3], double symprec)
{
    int j;
    double ratio_a, ratio_b, tmp_matrix_d3[3][3];

    for (j = 0; j < 3; j++) {
        lattice[j][c0] = vector_a[j] - min_lattice[j][axis];	/* Right handed ? */
        lattice[j][c1] = vector_a[j] + min_lattice[j][axis];
    }
    transpose_matrix_d3(tmp_matrix_d3, lattice);

    ratio_a = inner_product_d3(tmp_matrix_d3[c0], vector_b) /
        inner_product_d3(tmp_matrix_d3[c0], tmp_matrix_d3[c0]);

    ratio_b = inner_product_d3(tmp_matrix_d3[c1], vector_b) /
        inner_product_d3(tmp_matrix_d3[c1], tmp_matrix_d3[c1]);

    if (Dabs(Dabs(ratio_a) - 0.5) < symprec &&
        Dabs(Dabs(ratio_b) - 0.5) < symprec) {
        for (j = 0; j < 3; j++) {
            lattice[j][axis] =
                vector_b[j] - ratio_a * lattice[j][c0] -
                ratio_b * lattice[j][c1];
        }
        return 1;
    }
    return 0;
}

int bravais_rhombo1(int axis, int combination[3][2], double min_lattice[3][3],
                    double lattice[3][3], double symprec)
{
    int j, k, c0, c1, sign[3];
    double metric[3][3];

    c0 = combination[axis][0];
    c1 = combination[axis][1];

    get_metric(metric, min_lattice);

    if (Dabs(Dabs(metric[0][1]) - Dabs(metric[0][2])) < symprec &&
        Dabs(Dabs(metric[0][1]) - Dabs(metric[1][2])) < symprec) {
        copy_matrix_d3(lattice, min_lattice);

        for (j = 0; j < 3; j++) {
            sign[j] = 1;
            if (metric[combination[j][0]][combination[j][1]] < 0)
                sign[j] = -1;
        }

        for (j = 0; j < 3; j++) {
            if (sign[0] + sign[1] + sign[2] == -1) {
                if (sign[j] == 1)
                    for (k = 0; k < 3; k++)
                        lattice[k][j] = -lattice[k][j];
            } else if (sign[0] + sign[1] + sign[2] == 1)
                if (sign[j] == -1)
                    for (k = 0; k < 3; k++)
                        lattice[k][j] = -lattice[k][j];
        }
        return 1;
    }
    return 0;
}

int bravais_rhombo2(int axis, int c0, int c1, double min_lattice[3][3],
                    double lattice[3][3], double symprec)
{
    int j;
    double ratio_a, ratio_b, vector_a[3], vector_b[3], metric[3][3];

    get_metric(metric, min_lattice);

    for (j = 0; j < 3; j++) {
        vector_a[j] = min_lattice[j][c0];
        vector_b[j] = min_lattice[j][c1];
    }

    ratio_a = metric[axis][c0] / metric[axis][axis];
    ratio_b = metric[axis][c1] / metric[axis][axis];

    if (Dabs(Dabs(ratio_a) - 1.0 / 3) < symprec &&	/* |c|=3*|a|cos */
        Dabs(Dabs(ratio_b) - 1.0 / 3) < symprec) {	/* |c|=3*|b|cos */

        /* Projection onto plane normal c axis */
        for (j = 0; j < 3; j++) {
            vector_a[j] = vector_a[j] - ratio_a * min_lattice[j][axis];
            vector_b[j] = vector_b[j] - ratio_b * min_lattice[j][axis];
        }

        if (Dabs(Dabs(2 * inner_product_d3(vector_a, vector_b)) -
                 inner_product_d3(vector_a, vector_a)) < symprec) {

            /* angle between vec_a and vec_b is goint to 120 deg. */
            if (inner_product_d3(vector_a, vector_b) > 0)
                for (j = 0; j < 3; j++)
                    vector_b[j] = -vector_b[j];

            for (j = 0; j < 3; j++) {
                lattice[j][0] = min_lattice[j][axis] / 3 + vector_a[j];
                lattice[j][1] = min_lattice[j][axis] / 3 + vector_b[j];
                lattice[j][2] =
                    min_lattice[j][axis] / 3 - vector_a[j] - vector_b[j];
            }
            return 1;
        }
    }
    return 0;
}

int bravais_rhombo3(int axis, int c0, int c1, double min_lattice[3][3],
                    double lattice[3][3], double symprec)
{
    int j;
    double ratio_a, ratio_b, vector_a[3], vector_b[3], tmp_matrix_d3[3][3];

    for (j = 0; j < 3; j++) {
        vector_a[j] = min_lattice[j][c0] - min_lattice[j][c1];
        vector_b[j] = min_lattice[j][c0] + min_lattice[j][c1];
    }

    transpose_matrix_d3(tmp_matrix_d3, min_lattice);

    ratio_a = inner_product_d3(tmp_matrix_d3[axis], vector_a) /	/* ac/cc */
        inner_product_d3(tmp_matrix_d3[axis], tmp_matrix_d3[axis]);	/* |a|cos/|c| */

    ratio_b = inner_product_d3(tmp_matrix_d3[axis], vector_b) /	/* bc/cc */
        inner_product_d3(tmp_matrix_d3[axis], tmp_matrix_d3[axis]);	/* |b|cos/|c| */

    /* |vec_a| = |c| and |c|=2*|a||cos(angle)| (60,120 degs) */
    if (Dabs(inner_product_d3(tmp_matrix_d3[axis], tmp_matrix_d3[axis]) -
             inner_product_d3(vector_a, vector_a)) < symprec &&
        Dabs(Dabs(ratio_a) - 0.5) < symprec) {

        for (j = 0; j < 3; j++) {
            lattice[j][c0] = min_lattice[j][c0];
            lattice[j][c1] = min_lattice[j][c1];
            lattice[j][axis] =
                min_lattice[j][c1] + 2 * ratio_a * min_lattice[j][axis];
        }
        return 1;
    }

    /* |vec_b| = |c| (original b faces opposite) and |c|=2*|b||cos(angle)|| */
    if (Dabs(inner_product_d3(tmp_matrix_d3[axis], tmp_matrix_d3[axis]) -
             inner_product_d3(vector_b, vector_b)) < symprec &&
        Dabs(Dabs(ratio_b) - 0.5) < symprec) {

        for (j = 0; j < 3; j++) {
            lattice[j][c0] = min_lattice[j][c0];
            lattice[j][c1] = -min_lattice[j][c1];
            lattice[j][axis] =
                -min_lattice[j][c1] + 2 * ratio_b * min_lattice[j][axis];
        }
        return 1;
    }
    return 0;
}

int bravais_rhombo4(int axis, int c0, int c1, double min_lattice[3][3],
                    double lattice[3][3], double symprec)
{
    int j;
    double ratio_a, ratio_b, vector_a[3], vector_b[3], metric[3][3];
    get_metric(metric, min_lattice);

    for (j = 0; j < 3; j++) {
        vector_a[j] = min_lattice[j][c0];
        vector_b[j] = min_lattice[j][c1];
    }

    /* 2ab-aa=0, 2*|b|cos=|a|, theta=60,120 deg */
    if (Dabs(Dabs(2 * metric[c0][c1]) - metric[c0][c0]) < symprec) {
        ratio_a = metric[axis][c0] / metric[c0][c0];
        ratio_b = metric[axis][c1] / metric[c1][c1];

        if (Dabs(Dabs(ratio_a) - 0.5) < symprec
            && ( Dabs(Dabs(ratio_b) - 0.5) < symprec || Dabs(ratio_a) < symprec )
            && ( Dabs(Dabs(ratio_b) - 0.5) < symprec || Dabs(Dabs(ratio_a) - 0.5) < symprec)
            && Dabs(ratio_b) < symprec) {

            for (j = 0; j < 3; j++) {
                lattice[j][axis] = min_lattice[j][axis];

                if (Dabs(Dabs(ratio_a) - 0.5) < symprec)
                    lattice[j][c0] = min_lattice[j][axis] -
                        /* 120 deg *//* 120 deg *//* 120 deg */ 2 *
                        ratio_a * vector_a[j];

                if (Dabs(Dabs(ratio_b) - 0.5) < symprec)
                    lattice[j][c1] =
                        min_lattice[j][axis] - 2 * ratio_b * vector_b[j];

                if (Dabs(ratio_a) < symprec) {

                    /* angle between a and b is 60 deg. */
                    if (inner_product_d3(vector_a, vector_b) > 0)
                        lattice[j][c0] =
                            min_lattice[j][axis] - 2 * ratio_b * vector_b[j] -
                            vector_a[j];
                    else        /* 120 deg */
                        lattice[j][c0] =
                            min_lattice[j][axis] - 2 * ratio_b * vector_b[j] +
                            vector_a[j];
		}

                if (Dabs(ratio_b) < symprec){

                    /* angle between a and b is 60 deg. */
                    if (inner_product_d3(vector_a, vector_b) > 0)
                        lattice[j][c1] =
                            min_lattice[j][axis] - 2 * ratio_a * vector_a[j] -
                            vector_b[j];
                    else        /* 120 deg */
                        lattice[j][c1] =
                            min_lattice[j][axis] - 2 * ratio_a * vector_a[j] +
                            vector_b[j];
		}
            }
            return 1;
        }
    }
    return 0;
}

int bravais_monocli2(int axis, int c0, int c1, double min_lattice[3][3],
                     double lattice[3][3], double symprec)
{
    int j;
    double vectors[4][3];

    for (j = 0; j < 3; j++)     /* case 1 */
        vectors[0][j] = min_lattice[j][axis];

    if (bravais_monocli2_ext
        (axis, c0, c1, min_lattice, lattice, vectors[0], symprec))
        return 1;

    /* case 2 */
    for (j = 0; j < 3; j++) {
        vectors[0][j] = min_lattice[j][axis] - min_lattice[j][c0];
        vectors[1][j] = min_lattice[j][axis] + min_lattice[j][c0];
        vectors[2][j] = min_lattice[j][axis] - min_lattice[j][c1];
        vectors[3][j] = min_lattice[j][axis] + min_lattice[j][c1];
    }
    /* find shortest vector */
    qsort(vectors, 4, sizeof(vectors[0]), vector_compare);

    for (j = 0; j < 4; j++)     /* iterate from shorter to longer */
        if (bravais_monocli2_ext
            (axis, c0, c1, min_lattice, lattice, vectors[j], symprec))
            return 1;

    return 0;
}

int bravais_monocli2_ext(int axis, int c0, int c1, double min_lattice[3][3],
                         double lattice[3][3], double vector[3],
                         double symprec)
{
    int j;
    double vector_a[3], vector_b[3];

    for (j = 0; j < 3; j++) {
        vector_a[j] = min_lattice[j][c0] - min_lattice[j][c1];
        vector_b[j] = min_lattice[j][c0] + min_lattice[j][c1];
    }

    if (Dabs(inner_product_d3(vector, vector_a)) < symprec) {
        for (j = 0; j < 3; j++) {
            lattice[j][0] = vector[j];
            lattice[j][1] = vector_a[j];
            lattice[j][2] = vector_b[j];
        }
        return 1;
    }

    if (Dabs(inner_product_d3(vector, vector_b)) < symprec) {
        for (j = 0; j < 3; j++) {
            lattice[j][0] = vector[j];
            lattice[j][1] = vector_b[j];
            lattice[j][2] = vector_a[j];
        }
        return 1;
    }
    return 0;
}

int bravais_monocli3(int axis, int c0, int c1, double min_lattice[3][3],
                     double lattice[3][3], double symprec)
{
    int j;
    double vectors[4][3], vector_a[3], vector_b[3], vector_c[3], ratio_a,
        ratio_b, metric[3][3];

    get_metric(metric, min_lattice);

    ratio_a = metric[axis][c0] / metric[c0][c0];	/* ca/aa=|c|cos/|a| */
    ratio_b = metric[axis][c1] / metric[c1][c1];	/* cb/bb=|c|cos/|b| */

    if (Dabs(Dabs(ratio_b) - 0.5) > symprec
        && Dabs(Dabs(ratio_a) - 0.5) > symprec)
        return 0;               /* not one-face-centered */

    if (Dabs(Dabs(ratio_b) - 0.5) < symprec) {
        for (j = 0; j < 3; j++) {
            vector_a[j] = min_lattice[j][c0];
            vector_b[j] = min_lattice[j][c1];
            vector_c[j] = min_lattice[j][axis];
        }
    }

    if (Dabs(Dabs(ratio_a) - 0.5) < symprec) {	/* swap a and b */
        for (j = 0; j < 3; j++) {
            vector_a[j] = min_lattice[j][c1];
            vector_b[j] = -min_lattice[j][c0];
            vector_c[j] = min_lattice[j][axis];
        }
        ratio_b = ratio_a;
    }

    /* case 1 */
    if (bravais_monocli3_ext
        (lattice, vector_a, vector_b, vector_c, symprec))
        return 1;

    for (j = 0; j < 3; j++) {   /* case 2 */
        vectors[0][j] = vector_a[j] - vector_c[j];
        vectors[1][j] = vector_a[j] + vector_c[j];
        vectors[2][j] =
            vector_a[j] - vector_c[j] + 2 * ratio_b * vector_b[j];
        vectors[3][j] =
            vector_a[j] + vector_c[j] - 2 * ratio_b * vector_b[j];
    }

    /* find shortest vector */
    for (j = 0; j < 4; j++)
        debug_print("%f %f %f\n", vectors[j][0], vectors[j][1],
                    vectors[j][2]);

    qsort(vectors, 4, sizeof(vectors[0]), vector_compare);

    for (j = 0; j < 4; j++)
        debug_print("%f %f %f\n", vectors[j][0], vectors[j][1],
                    vectors[j][2]);

    for (j = 0; j < 4; j++)     /* iterate from shorter to longer */
        if (bravais_monocli3_ext
            (lattice, vectors[j], vector_b, vector_c, symprec))
            return 1;

    return 0;
}

int bravais_monocli3_ext(double lattice[3][3], double vector_a[3],
                         double vector_b[3], double vector_c[3],
                         double symprec)
{
    int j;

    if (Dabs(inner_product_d3(vector_a, vector_b)) < symprec) {
        for (j = 0; j < 3; j++) {
            lattice[j][0] = vector_a[j];
            lattice[j][1] = vector_b[j];
            lattice[j][2] = 2 * (vector_c[j] -
                                 inner_product_d3(vector_c, vector_b) /
                                 inner_product_d3(vector_b,
                                                  vector_b) * vector_b[j]);
        }
        return 1;
    }
    return 0;
}

int check_holohedry(double lattice[3][3], int holohedry, double symprec)
{
    double metric[3][3];
    int angle_90[3] = { 0, 0, 0 }, edge_equal[3] = {0, 0, 0};
    int orthogonal = 0, all_edge_equal = 0;

    get_metric(metric, lattice);
    check_angle90_equaledge(angle_90, edge_equal, metric, symprec);

    if (angle_90[0] && angle_90[1] && angle_90[2])
        orthogonal = 1;
    if (edge_equal[0] && edge_equal[1] && edge_equal[2])
        all_edge_equal = 1;

    if (holohedry == TRICLI)
        return 1;               /* triclinic */

    if (holohedry == MONOCLI && angle_90[0] && angle_90[2])
        return 1;               /* monoclinic */

    if (holohedry == ORTHO && orthogonal)
        return 1;               /* orthogonalrhombic */

    if (holohedry == TETRA && orthogonal
        && (edge_equal[0] || edge_equal[1] || edge_equal[2]))
        return 1;               /* tetragonal */

    if (holohedry == TRIGO && all_edge_equal
        && (Dabs(metric[0][1] - metric[1][2]) < symprec)
        && (Dabs(metric[0][1] - metric[0][2]) < symprec))
        return 1;               /* trigonal */

    if (holohedry == HEXA && edge_equal[2] && angle_90[0] && angle_90[1]
        && Dabs(2 * metric[0][1] + metric[0][0]) < symprec)
        return 1;               /* hexagonal */

    if (holohedry == CUBIC && orthogonal && all_edge_equal)
        return 1;               /* cubic */

    return 0;
}

void check_angle90_equaledge(int angle_90[3], int edge_equal[3],
                             double lattice[3][3], double symprec)
{
    int i, c0, c1;
    int combination[3][2] = { {1, 2}, {2, 0}, {0, 1} };	/* combination of axes for iteration */
    double metric[3][3];

    get_metric(metric, lattice);
    for (i = 0; i < 3; i++) {

        c0 = combination[i][0];
        c1 = combination[i][1];
        angle_90[i] = 0;
        edge_equal[i] = 0;

        if (Dabs(metric[c0][c1]) < symprec)	/* orthogonal */
            angle_90[i] = 1;

        if (Dabs(metric[c0][c0] - metric[c1][c1]) < symprec)	/* Equal edges */
            edge_equal[i] = 1;
    }

}

void get_smallest_primitive(double lattice_new[3][3], double lattice[3][3],
                            double symprec)
{
    int i, j, c0, c1, c2, flag = 1;
    double metric[3][3], tmp_matrix[3][3], projection[3][3];
    int combination[3][2] = { {1, 2}, {2, 0}, {0, 1} };	/* combination of axes for iteration */
    get_metric(metric, lattice);
    copy_matrix_d3(lattice_new, lattice);

    /*   Loop over all clear */
    while (flag) {
        flag = 0;

        for (i = 0; i < 3; i++) {
            c0 = combination[i][0];
            c1 = combination[i][1];

            /* Angle shoule be larger than 90 => Scalar product is negative */
            if (metric[c0][c1] > symprec) {
                lattice_new[0][c0] = -lattice_new[0][c0];
                lattice_new[1][c0] = -lattice_new[1][c0];
                lattice_new[2][c0] = -lattice_new[2][c0];
                get_metric(metric, lattice_new);
            }

            /* Replace the first (a) or second (b) vector by the sum vector (a+b) */
            for (j = 0; j < 2; j++) {
                c2 = combination[i][j];	/* take a or b */

                if (metric[c0][c0] + 2 * metric[c0][c1] + metric[c1][c1]
                    < metric[c2][c2] - symprec) {	/* |a+b| < |a| or |b| */

                    lattice_new[0][c2] =
                        lattice_new[0][c0] + lattice_new[0][c1];
                    lattice_new[1][c2] =
                        lattice_new[1][c0] + lattice_new[1][c1];
                    lattice_new[2][c2] =
                        lattice_new[2][c0] + lattice_new[2][c1];
                    get_metric(metric, lattice_new);
                    flag = 1;
                    break;
                }
            }
        }
    }

    debug_print("New lattice before forcing right handed orientation\n");
    debug_print_matrix_d3(lattice_new);

    /* Choose first vector as most overwrap with original first vector. */
    transpose_matrix_d3(tmp_matrix, lattice_new);
    multiply_matrix_d3(projection, tmp_matrix, lattice);	/* project */

    debug_print("projection a: %f %f %f\n", projection[0][0],
                projection[1][0], projection[2][0]);
    debug_print("projection b: %f %f %f\n", projection[0][1],
                projection[1][1], projection[2][1]);
    debug_print("projection c: %f %f %f\n", projection[0][2],
                projection[1][2], projection[2][2]);

    /* Choose first axis (well projected one) */
    i = 0;
    
    if (Dabs(projection[1][0]) - Dabs(projection[0][0]) > symprec)
        i = 1;

    if (Dabs(projection[2][0]) - Dabs(projection[i][0]) > symprec)
        i = 2;


    /* Swap axes */
    copy_matrix_d3(tmp_matrix, lattice_new);

    for (j = 0; j < 3; j++) {   
        lattice_new[j][0] = tmp_matrix[j][i];
        lattice_new[j][i] = tmp_matrix[j][0];
    }

    /* Flip first axis */
    if (projection[i][0] < -symprec)
        for (j = 0; j < 3; j++)
            lattice_new[j][0] = -lattice_new[j][0];



    /* Choose second axis (better projected one) */
    i = 1;

    if (Dabs(projection[2][0]) - Dabs(projection[1][0]) > symprec)
        i = 2;

    /* Swap axes */
    copy_matrix_d3(tmp_matrix, lattice_new);

    for (j = 0; j < 3; j++) {   
        lattice_new[j][1] = tmp_matrix[j][i];
        lattice_new[j][i] = tmp_matrix[j][1];
    }
    
    /* Flip second axis */
    if (projection[i][0] < -symprec)
        for (j = 0; j < 3; j++)
            lattice_new[j][1] = -lattice_new[j][1];
    

    /*   Right-handed orientation */
    if (get_determinant_d3(lattice_new) < -symprec*symprec*symprec) {

        /* Flip third axis */
        for (j = 0; j < 3; j++)
            lattice_new[j][2] = -lattice_new[j][2];
    }        
}
