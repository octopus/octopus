/* bravais.h */
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

#ifndef __bravasis_H__
#define __bravasis_H__

enum holohedry { TRICLI = 1, MONOCLI, ORTHO, TETRA, TRIGO, HEXA, CUBIC };
enum centering { NO_CENTER = 0, BODY = -1, FACE = -3, A_FACE = 1, B_FACE =
        2, C_FACE = 3 };

typedef struct {
    int rot[48][3][3];
    int size;
} LatticeSymmetry;

typedef struct {
    int holohedry;
    int centering;
    int factor;
    double lattice[3][3];
    LatticeSymmetry lattice_symmetry;
} Bravais;


LatticeSymmetry get_symmetry_candidate(Bravais * bravais,
                                       double lattice[3][3],
                                       double symprec);
LatticeSymmetry get_bravais_symmetry(int holohedry);
void generate_point_symmetry(int point_symmetry[][3][3],
                             int generator[3][3], int n_sym, int n_gen);
Bravais get_bravais_lattice(double lattice_orig[3][3], double symprec);
void set_bravais_lattice(Bravais *bravais, double symprec);
int get_bravais_lattice_in_loop(Bravais *bravais, double min_lattice[3][3],
                                int axis, double symprec);
int bravais_hexa(int i, int c0, int c1, double min_lattice[3][3],
                 double lattice[3][3], double symprec);
int bravais_ortho1(int i, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec);
int bravais_ortho3(int i, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec);
int bravais_ortho4(int i, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec);
int bravais_ortho5(int i, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec);
int bravais_ortho6(int i, int c0, int c1, double min_lattice[3][3],
                   double lattice[3][3], double symprec);
int bravais_ortho6_ext(int i, int c0, int c1, double min_lattice[3][3],
                       double lattice[3][3], double vector_a[3],
                       double vector_b[3], double symprec);
int bravais_rhombo1(int i, int combination[3][2], double min_lattice[3][3],
                    double lattice[3][3], double symprec);
int bravais_rhombo2(int i, int c0, int c1, double min_lattice[3][3],
                    double lattice[3][3], double symprec);
int bravais_rhombo3(int i, int c0, int c1, double min_lattice[3][3],
                    double lattice[3][3], double symprec);
int bravais_rhombo4(int i, int c0, int c1, double min_lattice[3][3],
                    double lattice[3][3], double symprec);
int bravais_monocli2(int i, int c0, int c1, double min_lattice[3][3],
                     double lattice[3][3], double symprec);
int bravais_monocli2_ext(int i, int c0, int c1, double min_lattice[3][3],
                         double lattice[3][3], double vector[3],
                         double symprec);
int bravais_monocli3(int i, int c0, int c1, double min_lattice[3][3],
                     double lattice[3][3], double symprec);
int bravais_monocli3_ext(double lattice[3][3], double vector_a[3],
                         double vector_b[3], double vector_c[3],
                         double symprec);
int check_holohedry(double lattice[3][3], int holohedry, double symprec);
void check_angle90_equaledge(int angle_90[3], int edge_equal[3],
                             double lattice[3][3], double symprec);
void get_smallest_primitive(double lattice_new[3][3], double lattice[3][3],
                            double symprec);

#endif
