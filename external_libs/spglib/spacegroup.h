/* spacegroup.h */
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

#ifndef __spacegroup_H__
#define __spacegroup_H__

#include "bravais.h"
#include "cell.h"
#include "symmetry.h"

typedef struct {
    char international[10];
    char schoenflies[10];
} Pointgroup;

typedef struct {
    int number;
    char international_long[40];
    char international[20];
    char schoenflies[10];
    char bravais_symbol[2];
    int multi;
    Pointgroup pointgroup;
} Spacegroup;

enum {
    IDENTITY = 8,
    PURE_TRANSLATION = 7,
    SECONDARY_2_AXIS = 4,
    TERTIARY_2_AXIS = 21,
    PRIMARY_2_AXIS = 9,
    PRIMARY_2_1_AXIS = 20,
    PRIMARY_3_AXIS = 10,
    PRIMARY_3_3_1_3_2_AXIS = 110,
    PRIMARY_3_1_AXIS = 22,
    PRIMARY_3_2_AXIS = 23,
    PRIMARY_4_AXIS = 12,
    PRIMARY_4_2_AXIS = 25,
    PRIMARY_4_1_4_3_AXIS = 124,
    PRIMARY_4_1_AXIS = 24,
    PRIMARY_4_3_AXIS = 26,
    PRIMARY_6_AXIS = 14,
    PRIMARY_6_3_AXIS = 29,
    PRIMARY_6_2_AXIS = 28,
    PRIMARY_6_4_AXIS = 30,
    PRIMARY_6_1_AXIS = 27,
    PRIMARY_6_5_AXIS = 31,
    INVERSION = 5,
    PRIMARY_M3_AXIS = 3,
    PRIMARY_M4_AXIS = 2,
    PRIMARY_M6_AXIS = 1,
    M_PLANE = 15,
    PRIMARY_N = 18,
    SECONDARY_A_B = 16,
    SECONDARY_C = 17,
    SECONDARY_N = 118,
    TERTIARY_M = 115,
    TERTIARY_C = 19,
    A_B_C_PLANE = 116,
    TERTIARY_D = 117,
    SECONDARY_M = 215,
    SECONDARY_C_TRIGO = 216,    /* ? */
    PRIMARY_M = 315,
    TERTIARY_M_HEXA = 317,
    TERTIARY_C_HEXA = 418,
    TERTIARY_N = 218,
    N_PLANE = 318,
    D_PLANE = 217,
    G_PLANE = 119
};

Spacegroup get_spacegroup(Cell * cell, double symprec);
Bravais get_primitive_bravais(Cell *cell, double symprec);
void get_correct_bravais(Bravais *bravais, Cell *cell,
                         int holohedry, double symprec);
int get_spacegroup_number(Bravais * bravais, Cell * primitive,
                          Symmetry * primitive_sym, double symprec);
int get_rotation_class(Bravais * bravais, int order, int rot[3][3],
                       double trans[3], double symprec);
int get_improper_rotation_class_2axis(Bravais * bravais, int rot[3][3],
                                      double raw_trans[3], double symprec);
int get_proper_rotation_class(Bravais * bravais, int order, int rot[3][3],
                              double raw_trans[3], double symprec);
int get_proper_rotation_class_6axis(int rot[3][3], int trans_order,
                                    int nonsymmorphic);
int get_proper_rotation_class_4axis(int centering, int rot[3][3],
                                    int trans_order, int nonsymmorphic[3]);
int get_proper_rotation_class_3axis(int holohedry, int rot[3][3],
                                    int trans_order, int nonsymmorphic);
int get_proper_rotation_class_2axis(Bravais * bravais, int rot[3][3],
                                    int trans_order);
int get_translation_order(int centering, int order, double raw_trans[3],
                          int rot[3][3], double symprec);
int get_class_order(int rot[3][3]);
double get_abs_diff_two_vectors(double a[3], double b[3]);
Symmetry get_conventional_symmetry(Bravais * bravais, Cell * primitive,
                                   Symmetry * primitive_sym,
                                   double symprec);
void get_order_times_translation(double trans[3], int order, int rot[3][3],
                                 double raw_trans[3]);
void print_rotation_class(int rot_class);

#endif
