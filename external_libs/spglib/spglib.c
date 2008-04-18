/* spglib.c */
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ispglib.h"
#include "bravais.h"
#include "mathfunc.h"
#include "primitive.h"
#include "cell.h"
#include "symmetry.h"
#include "pointgroup.h"
#include "spacegroup.h"
#include "spacegroup_database.h"

int spg_get_symmetry(int rotation[][3][3], double translation[][3],
                    int max_size, double lattice[3][3],
                    double position[][3], int types[], int num_atom,
                    double symprec)
{
    int i, j, size;
    Symmetry symmetry;
    Bravais bravais;
    LatticeSymmetry lattice_sym;
    Cell cell;

    cell = new_cell(num_atom);
    set_cell(&cell, lattice, position, types);

    /* get symmetry information of input cell */
    bravais = get_bravais_lattice(cell.lattice, symprec);
    lattice_sym = get_symmetry_candidate(&bravais, cell.lattice, symprec);
    symmetry = get_symmetry_operation(&cell, &lattice_sym, symprec);

    /* get correct Bravais lattice if the cell is primitive */
    if (check_primitive_cell(&symmetry) == 1) {

        bravais = get_primitive_bravais(&cell, symprec);
        lattice_sym = get_symmetry_candidate(&bravais, cell.lattice, symprec);
        delete_symmetry(&symmetry);
        symmetry = get_symmetry_operation(&cell, &lattice_sym, symprec);

    }

    if (symmetry.size > max_size) {
        fprintf(stderr, "Indicated max size(=%d) is less than number ",
                max_size);
        fprintf(stderr, "of symmetry operations(=%d).\n", symmetry.size);
        exit(1);
    }

    for (i = 0; i < symmetry.size; i++) {
        copy_matrix_i3(rotation[i], symmetry.rot[i]);
        for (j = 0; j < 3; j++)
            translation[i][j] = symmetry.trans[i][j];
    }

    size = symmetry.size;

    delete_cell(&cell);
    delete_symmetry(&symmetry);

    return size;
}

void spg_get_bravais_lattice(double bravais_lattice[3][3], double lattice[3][3],
                             double position[][3], int types[], int num_atom,
                             double symprec)
{
    int holohedry;
    Bravais bravais;
    Cell cell;
    Symmetry symmetry;
    LatticeSymmetry lattice_sym;

    cell = new_cell(num_atom);
    set_cell(&cell, lattice, position, types);

    bravais = get_primitive_bravais(&cell, symprec);

    copy_matrix_d3(bravais_lattice, bravais.lattice);

    delete_cell(&cell);
}

int spg_check_symmetry(double lattice[3][3], double position[][3],
                      int types[], int num_atom, double symprec)
{
    Symmetry symmetry;
    Bravais bravais;
    LatticeSymmetry lattice_sym;
    Cell cell;
    int size;

    cell = new_cell(num_atom);
    set_cell(&cell, lattice, position, types);

    /* get symmetry information of input cell */
    bravais = get_bravais_lattice(cell.lattice, symprec);
    lattice_sym = get_symmetry_candidate(&bravais, cell.lattice, symprec);
    symmetry = get_symmetry_operation(&cell, &lattice_sym, symprec);

    size = symmetry.size;

    delete_cell(&cell);
    delete_symmetry(&symmetry);

    return size;
}

/* lattice, position, and types are overwritten. num_atom is returned. */
int spg_find_primitive(double lattice[3][3], double position[][3],
                       int types[], int num_atom, double symprec)
{
    int i, j;
    Cell cell, primitive;
    Symmetry symmetry;
    Bravais bravais;
    LatticeSymmetry lattice_sym;

    cell = new_cell(num_atom);
    set_cell(&cell, lattice, position, types);

    /* get symmetry information of input cell */
    bravais = get_bravais_lattice(cell.lattice, symprec);
    lattice_sym = get_symmetry_candidate(&bravais, cell.lattice, symprec);
    symmetry = get_symmetry_operation(&cell, &lattice_sym, symprec);

    /* find primitive cell */
    if (check_primitive_cell(&symmetry) > 1) {

        primitive = get_primitive_cell(&cell, &symmetry, symprec);
        copy_matrix_d3(lattice, primitive.lattice);
        num_atom = primitive.size;

        for (i=0; i<num_atom; i++) {
            types[i] = primitive.types[i];
            
            for (j=0; j<3; j++)
                position[i][j] = primitive.position[i][j];
        }
        
        delete_cell(&primitive);
    }

    delete_cell(&cell);
    delete_symmetry(&symmetry);
    
    return num_atom;
}


void spg_show_symmetry(double lattice[3][3], double position[][3],
                      int types[], int num_atom, double symprec)
{
    Cell cell;
    Spacegroup spacegroup;
    cell = new_cell(num_atom);
    set_cell(&cell, lattice, position, types);

    spacegroup = get_spacegroup(&cell, symprec);

    if (spacegroup.number) {
        printf("Space group No.%d\n", spacegroup.number);
        printf(" International: %s%s\n", spacegroup.bravais_symbol,
               spacegroup.international);
        printf(" International(long): %s%s\n", spacegroup.bravais_symbol,
               spacegroup.international_long);
        printf(" Schoenflies: %s\n", spacegroup.schoenflies);
        printf(" Multiplicity: %d\n", spacegroup.multi);
        printf("Point group\n");
        printf(" International: %s\n", spacegroup.pointgroup.international);
        printf(" Schoenflies: %s\n", spacegroup.pointgroup.schoenflies);
    }
    
    delete_cell(&cell);
}

int spg_get_international(char symbol[21], double lattice[3][3],
                           double position[][3],
                           int types[], int num_atom, double symprec)
{
    Cell cell;
    Spacegroup spacegroup;
    cell = new_cell(num_atom);
    set_cell(&cell, lattice, position, types);

    spacegroup = get_spacegroup(&cell, symprec);
    strcpy(symbol, spacegroup.bravais_symbol);
    strcpy(&symbol[1], spacegroup.international);

    return spacegroup.number;
}
    
int spg_get_schoenflies(char symbol[10], double lattice[3][3],
                        double position[][3],
                        int types[], int num_atom, double symprec)
{
    Cell cell;
    Spacegroup spacegroup;
    cell = new_cell(num_atom);
    set_cell(&cell, lattice, position, types);

    spacegroup = get_spacegroup(&cell, symprec);
    strcpy(symbol, spacegroup.schoenflies);

    return spacegroup.number;
}
    


void spg_show_symmetry_operations(double lattice[3][3], double position[][3],
                                  int types[], int num_atom, double symprec)
{
    int i, j, order;
    Cell cell, primitive;
    Symmetry symmetry, conv_symmetry;
    Bravais bravais;
    LatticeSymmetry lattice_sym;

    cell = new_cell(num_atom);
    set_cell(&cell, lattice, position, types);

    /* get symmetry information of input cell */
    bravais = get_bravais_lattice(cell.lattice, symprec);
    lattice_sym = get_symmetry_candidate(&bravais, cell.lattice, symprec);
    symmetry = get_symmetry_operation(&cell, &lattice_sym, symprec);

    /* find primitive cell */
    if (check_primitive_cell(&symmetry) > 1) {

        printf("Cell is not primitive cell. The operations below are \n");
        printf("determined by primitive cell or ");
        printf("conventional unit cell.\n");

        primitive = get_primitive_cell(&cell, &symmetry, symprec);

        spg_show_symmetry_operations(primitive.lattice, primitive.position,
                                     primitive.types, primitive.size,
                                     symprec);

        delete_symmetry(&symmetry);
        delete_cell(&primitive);
        delete_cell(&cell);
        return;
    }


    /* get correct Bravais lattice if the cell is primitive */
    bravais = get_primitive_bravais(&cell, symprec);
    lattice_sym = get_symmetry_candidate(&bravais, cell.lattice, symprec);
    delete_symmetry(&symmetry);
    symmetry = get_symmetry_operation(&cell, &lattice_sym, symprec);

    conv_symmetry = get_conventional_symmetry(&bravais, &cell, &symmetry,
                                              symprec);
    

    printf("\n");
    printf("          Lattice                 Bravais lattice\n");
    printf("     a       b       c           a       b       c\n");
    for (i = 0; i<3; i++)
        printf("[%7.2f %7.2f %7.2f]   [%7.2f %7.2f %7.2f]\n",
               lattice[i][0], lattice[i][1], lattice[i][2],
               bravais.lattice[i][0], bravais.lattice[i][1],
               bravais.lattice[i][2]);
    printf("\n");

    for (i = 0; i < symmetry.size; i++) {

        order = get_class_order(conv_symmetry.rot[i]);

        printf("----------------------\n");

        printf("%3d: ", i + 1);

        print_rotation_class(
            get_rotation_class(&bravais, order, conv_symmetry.rot[i],
                               conv_symmetry.trans[i], symprec));

        printf("----------------------\n");

        for (j = 0; j<3; j++)
            printf("[%2d %2d %2d]\n", symmetry.rot[i][j][0], symmetry.rot[i][j][1],
                   symmetry.rot[i][j][2]);

        printf("\n");

        printf("[%5.2f,%5.2f,%5.2f]\n", symmetry.trans[i][0], symmetry.trans[i][1],
               symmetry.trans[i][2]);
        
        printf("\n");
    }

    delete_symmetry(&conv_symmetry);
    delete_symmetry(&symmetry);
    delete_cell(&cell);
}

