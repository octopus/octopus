/* mathfunc.c */
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
#include "mathfunc.h"

void get_metric(double metric[3][3], double lattice[3][3])
{
    double lattice_t[3][3];
    transpose_matrix_d3(lattice_t, lattice);
    multiply_matrix_d3(metric, lattice_t, lattice);
}

double get_determinant_d3(double a[3][3])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        + a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

double get_determinant_i3(int a[3][3])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        + a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

void copy_matrix_d3(double a[3][3], double b[3][3])
{
    a[0][0] = b[0][0];
    a[0][1] = b[0][1];
    a[0][2] = b[0][2];
    a[1][0] = b[1][0];
    a[1][1] = b[1][1];
    a[1][2] = b[1][2];
    a[2][0] = b[2][0];
    a[2][1] = b[2][1];
    a[2][2] = b[2][2];
}

void copy_matrix_i3(int a[3][3], int b[3][3])
{
    a[0][0] = b[0][0];
    a[0][1] = b[0][1];
    a[0][2] = b[0][2];
    a[1][0] = b[1][0];
    a[1][1] = b[1][1];
    a[1][2] = b[1][2];
    a[2][0] = b[2][0];
    a[2][1] = b[2][1];
    a[2][2] = b[2][2];
}

void copy_vector_d3(double a[3], double b[3])
{
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
}

int check_identity_matrix_i3(int a[3][3], int b[3][3])
{
    if (a[0][0] - b[0][0] ||
        a[0][1] - b[0][1] ||
        a[0][2] - b[0][2] ||
        a[1][0] - b[1][0] ||
        a[1][1] - b[1][1] ||
        a[1][2] - b[1][2] ||
        a[2][0] - b[2][0] || a[2][1] - b[2][1] || a[2][2] - b[2][2])
        return 0;
    else
        return 1;
}

/* m=axb */
void multiply_matrix_d3(double m[3][3], double a[3][3], double b[3][3])
{
    int i, j;                   /* a_ij */
    double c[3][3];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            c[i][j] =
                a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    copy_matrix_d3(m, c);
}

void multiply_matrix_i3(int m[3][3], int a[3][3], int b[3][3])
{
    int i, j;                   /* a_ij */
    int c[3][3];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            c[i][j] =
                a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    copy_matrix_i3(m, c);
}

/* m=axv */
void multiply_matrix_vector_d3(double v[3], double a[3][3], double b[3])
{
    int i;
    double c[3];
    for (i = 0; i < 3; i++)
        c[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2];
    for (i = 0; i < 3; i++)
        v[i] = c[i];
}

void multiply_matrix_vector_id3(double v[3], int a[3][3], double b[3])
{
    int i;
    double c[3];
    for (i = 0; i < 3; i++)
        c[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2];
    for (i = 0; i < 3; i++)
        v[i] = c[i];
}

void cast_matrix_3i_to_3d(double m[3][3], int a[3][3])
{
    m[0][0] = (double) a[0][0];
    m[0][1] = (double) a[0][1];
    m[0][2] = (double) a[0][2];
    m[1][0] = (double) a[1][0];
    m[1][1] = (double) a[1][1];
    m[1][2] = (double) a[1][2];
    m[2][0] = (double) a[2][0];
    m[2][1] = (double) a[2][1];
    m[2][2] = (double) a[2][2];
}

void cast_matrix_3d_to_3i(int m[3][3], double a[3][3])
{
    m[0][0] = Nint(a[0][0]);
    m[0][1] = Nint(a[0][1]);
    m[0][2] = Nint(a[0][2]);
    m[1][0] = Nint(a[1][0]);
    m[1][1] = Nint(a[1][1]);
    m[1][2] = Nint(a[1][2]);
    m[2][0] = Nint(a[2][0]);
    m[2][1] = Nint(a[2][1]);
    m[2][2] = Nint(a[2][2]);
}

/* m^-1 */
/* ruby code for auto generating */
/* 3.times {|i| 3.times {|j| */
/*       puts "m[#{j}][#{i}]=(a[#{(i+1)%3}][#{(j+1)%3}]*a[#{(i+2)%3}][#{(j+2)%3}] */
/*	 -a[#{(i+1)%3}][#{(j+2)%3}]*a[#{(i+2)%3}][#{(j+1)%3}])/det;" */
/* }} */
int inverse_matrix_d3(double m[3][3], double a[3][3], double precision)
{
    double det;
    double c[3][3];
    det = get_determinant_d3(a);
    if (Dabs(det) < precision) {
        fprintf(stderr, "No inverse matrix\n");
        return 0;
    }

    c[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    c[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    c[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    c[0][1] = (a[2][1] * a[0][2] - a[2][2] * a[0][1]) / det;
    c[1][1] = (a[2][2] * a[0][0] - a[2][0] * a[0][2]) / det;
    c[2][1] = (a[2][0] * a[0][1] - a[2][1] * a[0][0]) / det;
    c[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    c[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
    c[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
    copy_matrix_d3(m, c);
    return 1;
}

int inverse_matrix_i3(int m[3][3], int a[3][3])
{
    int det;
    int c[3][3];
    det = get_determinant_i3(a);
    if (!det)
        return 0;

    m[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    m[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    m[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    m[0][1] = (a[2][1] * a[0][2] - a[2][2] * a[0][1]) / det;
    m[1][1] = (a[2][2] * a[0][0] - a[2][0] * a[0][2]) / det;
    m[2][1] = (a[2][0] * a[0][1] - a[2][1] * a[0][0]) / det;
    m[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    m[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
    m[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
    copy_matrix_i3(m, c);
    return 1;
}

/* m = b^-1 a b */
int get_similar_matrix_d3(double m[3][3], double a[3][3], double b[3][3],
                          double precision)
{
    double c[3][3];
    if (!inverse_matrix_d3(c, b, precision)) {
        fprintf(stderr, "No similar matrix due to 0 determinant.\n");
        return 0;
    }
    multiply_matrix_d3(m, a, b);
    multiply_matrix_d3(m, c, m);
    return 1;
}

void print_matrix_d3(double a[3][3])
{
    int i;
    for (i = 0; i < 3; i++) {
        printf("%f %f %f\n", a[i][0], a[i][1], a[i][2]);
    }
}

void print_matrix_i3(int a[3][3])
{
    int i;
    for (i = 0; i < 3; i++) {
        printf("%d %d %d\n", a[i][0], a[i][1], a[i][2]);
    }
}

void print_vectors_d3(double a[][3], int size)
{
    int i;
    for (i = 0; i < size; i++) {
        printf("%d: %f %f %f\n", i + 1, a[i][0], a[i][1], a[i][2]);
    }
}

void print_vectors_with_label(double a[][3], int b[], int size)
{
    int i;
    for (i = 0; i < size; i++) {
        printf("%d: %f %f %f\n", b[i], a[i][0], a[i][1], a[i][2]);
    }
}

void transpose_matrix_d3(double a[3][3], double b[3][3])
{
    double c[3][3];
    c[0][0] = b[0][0];
    c[0][1] = b[1][0];
    c[0][2] = b[2][0];
    c[1][0] = b[0][1];
    c[1][1] = b[1][1];
    c[1][2] = b[2][1];
    c[2][0] = b[0][2];
    c[2][1] = b[1][2];
    c[2][2] = b[2][2];
    copy_matrix_d3(a, c);
}

void transpose_matrix_i3(int a[3][3], int b[3][3])
{
    int c[3][3];
    c[0][0] = b[0][0];
    c[0][1] = b[1][0];
    c[0][2] = b[2][0];
    c[1][0] = b[0][1];
    c[1][1] = b[1][1];
    c[1][2] = b[2][1];
    c[2][0] = b[0][2];
    c[2][1] = b[1][2];
    c[2][2] = b[2][2];
    copy_matrix_i3(a, c);
}

double inner_product_d3(double a[3], double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Dabs(double a)
{
    if (a < 0)
        return -a;
    else
        return a;
}

int Nint(double a)
{
    if (a < 0)
        return (int) (a - 0.5);
    else
        return (int) (a + 0.5);
}

int vector_compare(const void *_v0, const void *_v1)
{
    double *v0 = (double *) _v0;
    double *v1 = (double *) _v1;

    if (inner_product_d3(v0, v0) < inner_product_d3(v1, v1))
        return -1;
    else
        return 1;
}
