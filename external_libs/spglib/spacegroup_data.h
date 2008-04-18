/* spacegroup_data.h */
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

#ifndef __spacegroup_data_H__
#define __spacegroup_data_H__

#include "bravais.h"
#include "symmetry.h"

typedef struct {
    int class_table[32];
    int spacegroup;
} SpacegroupData;

int get_spacegroup_data(Symmetry * symmetry, Bravais * bravais,
                        int rot_class[], double symprec);
int get_spacegroup_data_special_case(Symmetry * symmetry,
                                     Bravais * bravais, int rot_class[],
                                     double symprec);
int check_class_table(int a[32], int b[32]);
void get_class_table(int class_talbe[32], int rot_class[], int size);
SpacegroupData get_spacegroup_data_table(int num);


#endif
