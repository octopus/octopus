/* test.c */
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

#define DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bravais.h"
#include "mathfunc.h"
#include "primitive.h"
#include "cell.h"
#include "symmetry.h"
#include "pointgroup.h"
#include "spacegroup.h"
#include "spglib.h"

int test_bravais_lattice(double symprec);
int test_get_bravais_symmetry();
int test_inverse_matrix_d3(double a[3][3]);
int test_get_symmetry_candidate();
int test_get_symmetry_operation();
Cell initial_setup(void);
int test_get_primitive_cell();
int test_cell();
Symmetry test_get_symmetry(Cell * cell);
int test_get_pointgroup_old();
void test_get_spacegroup();
int test_spg_get_symmetry();
int test_spg_check_symmetry();
void test_spg_show_symmetry();
void test_spg_find_primitive();
void test_spg_show_symmetry_operations();

int main()
{
    test_spg_show_symmetry_operations();
}

void test_spg_show_symmetry_operations()
{
    double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}};
    double position[][3] =
        {
            {0,0,0},
            {0.5,0.5,0.25},
            {0.3,0.3,0},
            {0.7,0.7,0},
            {0.2,0.8,0.25},
            {0.8,0.2,0.25},
            {0,0,0.5},
            {0.5,0.5,0.75},
            {0.3,0.3,0.5},
            {0.7,0.7,0.5},
            {0.2,0.8,0.75},
            {0.8,0.2,0.75}
        };
    int types[] = {1,1,2,2,2,2,1,1,2,2,2,2};
    int num_atom = 12;

    int size;

    spg_show_symmetry_operations(lattice, position, types, num_atom, 1e-5);
    
}

void test_spg_find_primitive()
{
    double lattice[3][3] = { {4, 0, 0}, {0, 4, 0}, {0, 0, 4} };
    double position[][3] = {
        {0, 0, 0},
        {0.5, 0.5, 0.5}
    };
    int types[] = { 1, 1 };
    int i, num_atom = 2;

    num_atom = spg_find_primitive(lattice, position, types, num_atom, 1e-5);
    
    printf("Lattice\n");
    print_matrix_d3(lattice);
    for (i=0; i<num_atom; i++)
        printf("%d: %f %f %f\n", types[i], position[i][0], position[i][1],
               position[i][2]);

}

void test_spg_show_symmetry()
{
/*     double lattice[3][3] = { {4, 0, 0}, {0, 4, 0}, {0, 0, 4} }; */
/*     double position[][3] = { */
/*         {0, 0, 0}, */
/*         {0.5, 0.5, 0.5} */
/*     }; */
/*     int types[] = { 1, 2 }; */
/*     int num_atom = 2; */

    double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}};
    double position[][3] =
        {
            {0,0,0},
            {0.5,0.5,0.25},
            {0.3,0.3,0},
            {0.7,0.7,0},
            {0.2,0.8,0.25},
            {0.8,0.2,0.25},
            {0,0,0.5},
            {0.5,0.5,0.75},
            {0.3,0.3,0.5},
            {0.7,0.7,0.5},
            {0.2,0.8,0.75},
            {0.8,0.2,0.75}
        };
    int types[] = {1,1,2,2,2,2,1,1,2,2,2,2};
    int num_atom = 12;

    int size;

    spg_show_symmetry(lattice, position, types, num_atom, 1e-5);
}

int test_spg_check_symmetry()
{
    double lattice[3][3] = { {4, 0, 0}, {0, 4, 0}, {0, 0, 4} };
    double position[][3] = {
        {0, 0, 0},
        {0.5, 0.5, 0.5}
    };
    int types[] = { 1, 2 };
    int num_atom = 2;
    int size;

    size = spg_check_symmetry(lattice, position, types, num_atom, 1e-5);
    printf("Number of symmetry operations: %d\n", size);
}


int test_spg_get_symmetry()
{
    double lattice[3][3] = { {4, 0, 0}, {0, 4, 0}, {0, 0, 4} };
    double position[][3] = {
        {0, 0, 0},
        {0.5, 0.5, 0.5}
    };
    int types[] = { 1, 2 };
    int num_atom = 2;
    int max_size = 50;
    int i, j, size;
    int rotation[max_size][3][3];
    double translation[max_size][3];

    size =
        spg_get_symmetry(rotation, translation, max_size, lattice, position,
                        types, num_atom, 1e-5);

    for (i = 0; i < size; i++) {
        printf("--- %d ---\n", i + 1);
        print_matrix_i3(rotation[i]);
        printf("%f %f %f\n", translation[i][0], translation[i][1],
               translation[i][2]);
    }

}

void test_get_spacegroup()
{
    Cell cell;
    cell = initial_setup();

    get_spacegroup(&cell, 1e-5);
    delete_cell(&cell);
}


int test_get_spacegroup_old()
{
    Cell cell, primitive, tmpcell;
    int i;
    Symmetry symmetry;
    Bravais bravais;
    double symprec = 1e-8;
    char pointgroup[6];
    LatticeSymmetry lattice_sym;
    cell = initial_setup();

    /* get symmetry information of input cell */
    bravais = get_bravais_lattice(cell.lattice, symprec);
    lattice_sym = get_symmetry_candidate(&bravais, cell.lattice, symprec);
    symmetry = get_symmetry_operation(&cell, &lattice_sym, symprec);

    /* SEE HERE ! */
    printf("Supercell symmetry if cell is supercell\n");
    for (i = 0; i < symmetry.size; i++) {
        printf("--- %d ---\n", i + 1);
        print_matrix_i3(symmetry.rot[i]);
        printf("%f %f %f\n", symmetry.trans[i][0], symmetry.trans[i][1],
               symmetry.trans[i][2]);
    }

    /* find primitive cell */
    primitive = get_primitive_cell(&cell, &symmetry, symprec);

    delete_symmetry(&symmetry);
    delete_cell(&cell);

    /* get symmetry information of primitive cell */
    bravais = get_bravais_lattice(primitive.lattice, symprec);
    lattice_sym =
        get_symmetry_candidate(&bravais, primitive.lattice, symprec);
    symmetry =
        get_symmetry_operation(&primitive, &lattice_sym, symprec);
    bravais.holohedry = check_pointgroup(bravais.holohedry, &symmetry);
    printf("%s\n", pointgroup);
    get_spacegroup_number(&bravais, &primitive, &symmetry, symprec);

    delete_symmetry(&symmetry);
    delete_cell(&primitive);
}

int test_check_pointgroup()
{
    Symmetry symmetry;
    Bravais bravais;
    char pointgroup[6];
    Cell cell;
    double symprec = 1e-5;
    cell = initial_setup();
    bravais = get_bravais_lattice(cell.lattice, symprec);
    symmetry = test_get_symmetry(&cell);
    bravais.holohedry = check_pointgroup(bravais.holohedry, &symmetry);
    printf("%s\n", pointgroup);
}

Cell initial_setup(void)
{
    Cell cell;
/*     double lattice[3][3] = { {0, 0.5, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0} }; */
/*     double lattice[3][3] = { {1, -0.5, 0}, {0, sqrt(3)/2, 0}, {0, 0, 3} }; */
/*     double position[][3] = { */
/*         {0, 0, 0} */
/*         {0.5, 0.5, 0.5} */
/*     }; */
/*     int types[] = { 1}; */

/*   double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}}; */
/*   double position[][3] = */
/*     { */
/*       {0,0,0}, */
/*       {0.5,0.5,0.5}, */
/*       {0.3,0.3,0}, */
/*       {0.7,0.7,0}, */
/*       {0.2,0.8,0.5}, */
/*       {0.8,0.2,0.5}, */
/*     }; */
/*   int types[] = {1,1,2,2,2,2}; */
/*     cell = new_cell(6); */

/*   double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}}; */
/*   double position[][3] = */
/*     { */
/*       {0,0,0}, */
/*       {0.5,0.5,0.25}, */
/*       {0.3,0.3,0}, */
/*       {0.7,0.7,0}, */
/*       {0.2,0.8,0.25}, */
/*       {0.8,0.2,0.25}, */
/*       {0,0,0.5}, */
/*       {0.5,0.5,0.75}, */
/*       {0.3,0.3,0.5}, */
/*       {0.7,0.7,0.5}, */
/*       {0.2,0.8,0.75}, */
/*       {0.8,0.2,0.75} */
/*     }; */
/*   int types[] = {1,1,2,2,2,2,1,1,2,2,2,2}; */
/*     cell = new_cell(12); */

    double lattice[3][3] = { {4, 0, 0}, {0, 4, 0}, {0, 0, 4} };
    double position[][3] = {
        {0, 0, 0},
        {0.5, 0.5, 0.5}
    };
    int types[] = { 1, 2 };
    cell = new_cell(2);

    set_cell(&cell, lattice, position, types);
    return cell;
}

int test_cell()
{
    int i, j;
    Cell cell;
    cell = initial_setup();
    double position[cell.size][3];
    int types[cell.size];
    get_cell_position(position, &cell);
    get_cell_types(types, &cell);
    for (i = 0; i < 6; i++)
        printf("%d: %f %f %f\n", types[i], position[i][0], position[i][1],
               position[i][2]);
    delete_cell(&cell);
}

int test_get_primitive_cell()
{
    Cell cell, primitive;
    int i;
    Symmetry symmetry;
    double symprec = 1e-5;
    cell = initial_setup();
    symmetry = test_get_symmetry(&cell);

    primitive = get_primitive_cell(&cell, &symmetry, symprec);
    delete_symmetry(&symmetry);

    symmetry = test_get_symmetry(&primitive);
    delete_symmetry(&symmetry);

    delete_cell(&primitive);
    delete_cell(&cell);
}

Symmetry test_get_symmetry(Cell * cell)
{
    int i, point_symmetry[48][3][3];
    double symprec = 1e-5;
    int num_point_sym;
    Bravais bravais;
    Symmetry symmetry;
    LatticeSymmetry lattice_sym;

    bravais = get_bravais_lattice(cell->lattice, symprec);

    debug_print("Original lattice\n");
    debug_print_matrix_d3(cell->lattice);
    debug_print("Bravais lattice\n");
    debug_print_matrix_d3(bravais.lattice);
    if (!num_point_sym) {
        fprintf(stderr, "BUG: Symmetry candidate could not be found.");
        exit(1);
    }
    lattice_sym = get_symmetry_candidate(&bravais, cell->lattice, symprec);
    symmetry = get_symmetry_operation(cell, &lattice_sym, symprec);
    for (i = 0; i < symmetry.size; i++) {
        printf("----- %i ----\n", i + 1);
        print_matrix_i3(symmetry.rot[i]);
        printf("%f %f %f\n", symmetry.trans[i][0], symmetry.trans[i][1],
               symmetry.trans[i][2]);
    }
    return symmetry;
}

int test_get_symmetry_operation()
{
    Cell cell, primitive;
    cell = initial_setup();
    int i, point_symmetry[48][3][3];
    double symprec = 1e-5;
    int num_point_sym;
    Bravais bravais;
    Symmetry symmetry;
    LatticeSymmetry lattice_sym;

    bravais = get_bravais_lattice(cell.lattice, symprec);
    debug_print_matrix_d3(bravais.lattice);
    lattice_sym = get_symmetry_candidate(&bravais, cell.lattice, symprec);
    symmetry = get_symmetry_operation(&cell, &lattice_sym, symprec);
    primitive = get_primitive_cell(&cell, &symmetry, symprec);

    delete_cell(&cell);
    delete_symmetry(&symmetry);
}

int test_inverse_matrix_d3(double a[3][3])
{
    double m[3][3];
    double b[3][3];
    debug_print("original\n");
    debug_print_matrix_d3(a);
    if (!inverse_matrix_d3(m, a, 1e-10)) {
        debug_print("no inverse matrix\n");
        return 0;
    }
    debug_print("inverse\n");
    debug_print_matrix_d3(m);
    multiply_matrix_d3(b, m, a);
    debug_print("multiply\n");
    debug_print_matrix_d3(b);
}

int test_get_bravais_symmetry()
{
    LatticeSymmetry lattice_sym;
    int i;
    for (i = 7; i > 0; i--)
        lattice_sym = get_bravais_symmetry(i);
}

int test_bravais_lattice(double symprec)
{
    /* ortho */
/*   double lattice[3][3] = {{1,0,0},{0,2,0},{0,0,3}}; */
    /* hexagonal */
/*   double lattice[3][3] = {{1,0.5,0},{0,sqrt(3)/2,0},{0,0,10}}; */
    /* tetragonal body centered */
/*   double lattice[3][3] = {{2,0,1},{0,2,1},{0,0,5}}; */
    /* ortho 1 C-centered */
/*     double lattice[3][3] = {{1,0.5,0},{0,0.5,0},{0,0,2}}; */
    /* ortho 1 B-centered */
/*     double lattice[3][3] = {{0.5,0,0},{0,1,0},{0.5,0,1}}; */
    /* ortho 1 A-centered */
/*     double lattice[3][3] = {{1,0,0},{0,0.5,0},{0,0.5,1}}; */
    /* ortho 1 body-centered */
/*     double lattice[3][3] = {{1,0,0},{0,2,0},{0.5,0.5,1}}; */
    /* ortho 2 C-centered */
/*     double lattice[3][3] = {{2,-2,0},{1,1,0},{0,0,1}}; */
    /* ortho 3 body-centered */
/*     double lattice[3][3] = {{2,2,0},{-3,3,0},{1,1,2}}; */
    /* ortho 3 face-centered */
/*     double lattice[3][3] = {{1,0,0},{0,2,0},{1,1,2}}; */
    /* ortho 3 body-centered */
/*     double lattice[3][3] = {{1,0,0},{0,2,0},{1,1,2}}; */
    /* ortho 4 body-centered */
/*     double lattice[3][3] = {{1,-1,1},{2,2,2},{3,3,-3}}; */
    /* ortho 5 face-centered */
/*     double lattice[3][3] = {{1,-1,0},{0,0,2},{3,3,3}}; */
    /* ortho 6 face-centered */
/*     double lattice[3][3] = {{0,-1,1},{2,0,2},{3,-3,0}}; */
    /* rhombo 1 */
/*     double lattice[3][3] = {{1,-0.5,-0.5},{0,sqrt(3)/2,-sqrt(3)/2},{1,1,1}}; */
    /* rhombo 2 */
/*     double lattice[3][3] = {{1,-0.5,0},{0,sqrt(3)/2,0},{1,1,3}}; */
    /* rhombo 3 */
/*     double lattice[3][3] = {{1,-0.5,0},{0,sqrt(3)/2,-sqrt(3)},{1,1,0}}; */
    /* rhombo 4 */
/*     double lattice[3][3] = {{1,-1.5,-1.5},{0,sqrt(3)/2,-sqrt(3)/2},{1,0,0}}; */
    /* monocli 1 */
/*     double lattice[3][3] = {{1,0,0.1},{0,2,0},{0,0,3}}; */
    /* monocli 2 */
/*     double lattice[3][3] = {{1,0.05,-0.05},{0,1,1},{0,1.5,-1.5}}; */
/*     double lattice[3][3] = {{1.05,0.05,-0.05},{1,1,1},{1.5,1.5,-1.5}}; */
/*     double lattice[3][3] = {{0.95,0.05,-0.05},{1,1,1},{-1.5,1.5,-1.5}}; */
/*   double lattice[3][3] = {{0.95,0.05,-0.05},{-1,1,1},{-1.5,1.5,-1.5}}; */
/*   double lattice[3][3] = {{1.05,0.05,-0.05},{-1,1,1},{1.5,1.5,-1.5}}; */
    /* monocli 3 case1 */
/*     double lattice[3][3] = {{1,0,0.05},{0,2,1},{0,0,1.5}}; */
    /* monocli 3 case1 */
/*     double lattice[3][3] = {{0,1,0.05},{2,0,1},{0,0,-1.5}}; */
    /* monocli 3 case2 */
/*     double lattice[3][3] = {{1.05,0,0.05},{1,2,1},{1.5,0,1.5}}; */
    double lattice[3][3] = { {0.95, 0, 0.05}, {-1, 2, 1}, {-1.5, 0, 1.5} };
/*     double lattice[3][3] = {{0.95,0,0.05},{1,2,1},{-1.5,0,1.5}}; */

    Bravais bravais;
    bravais = get_bravais_lattice(lattice, symprec);
    debug_print_matrix_d3(bravais.lattice);
}
