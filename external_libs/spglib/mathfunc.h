/* mathfunc.h */
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

#ifndef __mathfunc_H__
#define __mathfunc_H__

#ifdef DEBUG
#define debug_print(...) printf(__VA_ARGS__)
#define debug_print_matrix_d3( a ) print_matrix_d3( a )
#define debug_print_matrix_i3( a ) print_matrix_i3( a )
#define debug_print_vectors_d3(...) print_vectors_d3(__VA_ARGS__)
#define debug_print_vectors_with_label(...) print_vectors_with_label(__VA_ARGS__)
#else
#define debug_print(...)
#define debug_print_matrix_d3( a )
#define debug_print_matrix_i3( a )
#define debug_print_vectors_d3(...)
#define debug_print_vectors_with_label(...)
#endif

void get_metric(double metric[3][3], double lattice[3][3]);
double get_determinant_d3(double a[3][3]);
double get_determinant_i3(int a[3][3]);
void copy_matrix_d3(double a[3][3], double b[3][3]);
void copy_matrix_i3(int a[3][3], int b[3][3]);
void copy_vector_d3(double a[3], double b[3]);
int check_identity_matrix_i3(int a[3][3], int b[3][3]);
void multiply_matrix_d3(double m[3][3], double a[3][3], double b[3][3]);
void multiply_matrix_i3(int m[3][3], int a[3][3], int b[3][3]);
void multiply_matrix_vector_d3(double v[3], double a[3][3], double b[3]);
void cast_matrix_3i_to_3d(double m[3][3], int a[3][3]);
void cast_matrix_3d_to_3i(int m[3][3], double a[3][3]);
int inverse_matrix_d3(double m[3][3], double a[3][3], double precision);
int inverse_matrix_i3(int m[3][3], int a[3][3]);
int get_similar_matrix_d3(double m[3][3], double a[3][3], double s[3][3],
                          double precision);
void print_matrix_d3(double a[3][3]);
void print_matrix_i3(int a[3][3]);
void print_vectors_d3(double a[][3], int size);
void print_vectors_with_label(double a[][3], int b[], int size);
void transpose_matrix_d3(double a[3][3], double b[3][3]);
void transpose_matrix_i3(int a[3][3], int b[3][3]);
double inner_product_d3(double a[3], double b[3]);
double Dabs(double a);
int Nint(double a);
int vector_compare(const void *_v0, const void *_v1);

#endif
