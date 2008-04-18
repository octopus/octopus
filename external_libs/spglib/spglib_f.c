#include "spglib.h"
#include <string.h>

void spg_check_symmetry_(int *size, double lattice[3][3], double position[][3],
                         int types[], int *num_atom, double *symprec);
void spg_get_symmetry_(int *nsym, int rot[][3][3], double trans[][3], int *size,
                         double lattice[3][3], double position[][3], int types[],
                       int *num_atom, double *symprec);
void spg_get_bravais_lattice_(double bravais_lattice[3][3], double lattice[3][3],
                             double position[][3], int types[], int *num_atom,
                              double *symprec);
void spg_get_international_(int *spacegroup, char symbol[21], double lattice[3][3],
                           double position[][3], int types[], int *num_atom,
                           double *symprec);
void spg_get_schoenflies_(int *spacegroup, char symbol[10], double lattice[3][3],
                         double position[][3], int types[], int *num_atom,
                         double *symprec);
void spg_find_primitive_(double lattice[3][3], double position[][3],
                        int types[], int *num_atom, double *symprec);





void spg_check_symmetry_(int *size, double lattice[3][3], double position[][3],
                         int types[], int *num_atom, double *symprec)
{
    *size = spg_check_symmetry(lattice, position, types, *num_atom, *symprec);
}

void spg_get_symmetry_(int *nsym, int rot[][3][3], double trans[][3], int *size,
                         double lattice[3][3], double position[][3], int types[],
                         int *num_atom, double *symprec)
{
    *nsym = spg_get_symmetry(rot, trans, *size, lattice, position,
                             types, *num_atom, *symprec);
}

void spg_get_bravais_lattice_(double bravais_lattice[3][3], double lattice[3][3],
                             double position[][3], int types[], int *num_atom,
                             double *symprec)
{
    spg_get_bravais_lattice(bravais_lattice, lattice, position, types,
                            *num_atom, *symprec);
}

void spg_get_international_(int *spacegroup, char symbol[21], double lattice[3][3],
                           double position[][3], int types[], int *num_atom,
                           double *symprec)
{
    char symbol_c[21];
    int i, length;

    *spacegroup = spg_get_international(symbol_c, lattice, position, types,
                                        *num_atom, *symprec);
    length = strlen(symbol_c);
    strncpy(symbol, symbol_c, length);

    for (i=length; i<21; i++)
        symbol[i] = ' ';
}


void spg_get_schoenflies_(int *spacegroup, char symbol[10], double lattice[3][3],
                         double position[][3], int types[], int *num_atom,
                         double *symprec)
{
    char symbol_c[10];
    int i, length;
    
    *spacegroup = spg_get_schoenflies(symbol_c, lattice, position, types,
                                      *num_atom, *symprec);
    length = strlen(symbol_c);
    strncpy(symbol, symbol_c, length);

    for (i=length; i<10; i++)
        symbol[i] = ' ';
}

void spg_find_primitive_(double lattice[3][3], double position[][3],
                       int types[], int *num_atom, double *symprec)
{

    *num_atom = spg_find_primitive(lattice, position, types, *num_atom,
                                    *symprec);
}
