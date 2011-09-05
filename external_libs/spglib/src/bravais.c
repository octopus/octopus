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

#include <stdio.h>
#include <stdlib.h>
#include "bravais.h"
#include "debug.h"
#include "mathfunc.h"

static int get_brv_cubic(Bravais *bravais, const double min_lattice[3][3],
			 const double symprec);
static int get_brv_tetra(Bravais *bravais, const double min_lattice[3][3], const double symprec);
static int get_brv_hexa(Bravais *bravais, const double min_lattice[3][3], const double symprec);
static int get_brv_rhombo(Bravais *bravais, const double min_lattice[3][3], const double symprec);
static int get_brv_ortho(Bravais *bravais, const double min_lattice[3][3], const double symprec);
static int get_brv_monocli(Bravais *bravais, const double min_lattice[3][3], const double symprec);
static int brv_exhaustive_search(double lattice[3][3], const double min_lattice[3][3],
				 int (*check_bravais)(const double lattice[3][3], const double symprec),
				 const int relative_axes[][3], const int num_axes, const Centering centering,
				 const double symprec);
static int brv_cubic_I_center(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_cubic_F_center(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_tetra_primitive(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_tetra_one(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_tetra_two(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_tetra_three(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static Centering brv_ortho_base_center(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static Centering get_base_center(const double brv_lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_ortho_I_center(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_ortho_F_center(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_rhombo_two(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static int brv_rhombo_three(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static void set_brv_monocli(Bravais *bravais, const double symprec);
static int brv_monocli_primitive(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static Centering brv_monocli_base_center(double lattice[3][3], const double min_lattice[3][3], const double symprec);
static void check_angle90(int angle_90[3], const double lattice[3][3], const double symprec);
static void check_equal_edge(int edge_equal[3], const double lattice[3][3], const double symprec);
static int check_cubic(const double lattice[3][3], const double symprec);
static int check_hexa(const double lattice[3][3], const double symprec);
static int check_tetra(const double lattice[3][3], const double symprec);
static int check_ortho(const double lattice[3][3], const double symprec);
static int check_monocli(const double lattice[3][3], const double symprec);
static int check_rhombo(const double lattice[3][3], const double symprec);
static void get_right_hand_lattice(double lattice[3][3], const double symprec);
static void get_projection(double projection[3][3], const double min_lattice[3][3], const double lattice[3][3]);

/* math */
static void get_metric(double metric[3][3], const double lattice[3][3]);

/**********************/
/**********************/
/** Public functions **/
/**********************/
/**********************/
Bravais brv_get_brv_lattice(const double lattice_orig[3][3], const double symprec)
{
  Bravais bravais;
  double min_lattice[3][3];
  int i;
  Holohedry holohedries[] = {
    CUBIC,
    HEXA,
    RHOMB,
    TETRA,
    ORTHO,
    MONOCLI,
    TRICLI
  };

  brv_smallest_lattice_vector(min_lattice, lattice_orig, symprec);

#ifdef DEBUG
  double metric[3][3];
  int angle_90[3], edge_equal[3];
  debug_print("*** brv_get_brv_lattice ***\n");
  debug_print("Original lattice\n");
  debug_print_matrix_d3(lattice_orig);
  debug_print("Minimum lattice\n");
  debug_print_matrix_d3(min_lattice);
  debug_print("Metric tensor of minimum lattice\n");
  get_metric(metric, min_lattice);
  debug_print_matrix_d3(metric);
  check_angle90(angle_90, min_lattice, symprec);
  check_equal_edge(edge_equal, min_lattice, symprec);
  printf("angle_90: %d %d %d\n", angle_90[0], angle_90[1], angle_90[2]);
  printf("equal:    %d %d %d\n", edge_equal[0], edge_equal[1], edge_equal[2]);
#endif

  for (i = 0; i < 7; i++) {
    bravais.holohedry = holohedries[i];
    if (brv_get_brv_lattice_in_loop(&bravais, min_lattice, symprec))
      break;
  }

  debug_print("Bravais lattice\n");
  debug_print_holohedry(&bravais);
  debug_print_matrix_d3(bravais.lattice);
  debug_print_matrix_d3(bravais.lattice);

  return bravais;
}

/* Note: bravais is overwritten. */
int brv_get_brv_lattice_in_loop(Bravais *bravais, const double min_lattice[3][3],
                                const double symprec)
{
  switch (bravais->holohedry) {

  case CUBIC:
    if (get_brv_cubic(bravais, min_lattice, symprec))
      goto ok;
    break;
    
  case TETRA:
    if (get_brv_tetra(bravais, min_lattice, symprec))
      goto ok;
    break;
  
  case ORTHO:
    if (get_brv_ortho(bravais, min_lattice, symprec))
      goto ok;
    break;

  case HEXA:
  case TRIGO:
    if (get_brv_hexa(bravais, min_lattice, symprec))
      goto ok;
    break;

  case RHOMB:
    if (get_brv_rhombo(bravais, min_lattice, symprec))
      goto ok;
    break;

  case MONOCLI:
    if (get_brv_monocli(bravais, min_lattice, symprec))
      goto ok;
    break;

  case TRICLI:
  default:
    debug_print("triclinic\n");
    mat_copy_matrix_d3(bravais->lattice, min_lattice);
    bravais->centering = NO_CENTER;
    goto ok;
  }

  return 0;
  
 ok:
  /* Flip if determinant is minus. */
  get_right_hand_lattice(bravais->lattice, symprec);
  return 1;
}

void brv_smallest_lattice_vector(double min_lattice[3][3], const double lattice[3][3],
				 const double symprec)
{
  int i, j, c0, c1, c2, ok_length = 1;
  double metric[3][3], tmp_matrix[3][3], projection[3][3];

  get_metric(metric, lattice);
  mat_copy_matrix_d3(min_lattice, lattice);

  debug_print("*** brv_smallest_lattice_vector ***\n");

  /*   Loop over being all clear */
  while (ok_length) {
    ok_length = 0;

    for (i = 0; i < 3; i++) {
      c0 = (i + 1) % 3;
      c1 = (i + 2) % 3;

      /* Angle shoule be larger than 90 => Scalar product is negative */
      if (metric[c0][c1] > symprec) {
	min_lattice[0][c0] = -min_lattice[0][c0];
	min_lattice[1][c0] = -min_lattice[1][c0];
	min_lattice[2][c0] = -min_lattice[2][c0];
	get_metric(metric, min_lattice);
      }

      /* Replace the first (a) or second (b) vector by the sum vector (a+b) */
      for (j = 0; j < 2; j++) {
	c2 = (i + j + 1) % 3; /* take a or b */

	if (metric[c0][c0] + 2 * metric[c0][c1] + metric[c1][c1]
	    < metric[c2][c2] - symprec) {	/* |a+b| < |a| or |b| */

	  min_lattice[0][c2] =
	    min_lattice[0][c0] + min_lattice[0][c1];
	  min_lattice[1][c2] =
	    min_lattice[1][c0] + min_lattice[1][c1];
	  min_lattice[2][c2] =
	    min_lattice[2][c0] + min_lattice[2][c1];
	  get_metric(metric, min_lattice);

	  ok_length = 1;
	  break;
	}
      }

    }
  }


  debug_print("New lattice before forcing right handed orientation\n");
  debug_print_matrix_d3(min_lattice);

  /* Choose first vector as most overwrapping with the original first vector. */
  get_projection(projection, min_lattice, lattice);

  debug_print("Projection\n");
  debug_print("%f %f %f\n", projection[0][0], projection[1][0], projection[2][0]);
  debug_print("%f %f %f\n", projection[0][1], projection[1][1], projection[2][1]);
  debug_print("%f %f %f\n", projection[0][2], projection[1][2], projection[2][2]);

  /* Choose first axis (well projected one) */
  i = 0;
    
  if (mat_Dabs(projection[1][0]) - mat_Dabs(projection[0][0]) > symprec)
    i = 1;

  if (mat_Dabs(projection[2][0]) - mat_Dabs(projection[i][0]) > symprec)
    i = 2;


  /* Swap axes */
  mat_copy_matrix_d3(tmp_matrix, min_lattice);

  for (j = 0; j < 3; j++) {
    min_lattice[j][0] = tmp_matrix[j][i];
    min_lattice[j][i] = tmp_matrix[j][0];
  }

  /* Flip first axis */
  if (projection[i][0] < -symprec)
    for (j = 0; j < 3; j++)
      min_lattice[j][0] = -min_lattice[j][0];


  /* Choose second axis (better projected one) */
  i = 1;

  if (mat_Dabs(projection[2][0]) - mat_Dabs(projection[1][0]) > symprec) {
    i = 2;

    /* Swap axes */
    mat_copy_matrix_d3(tmp_matrix, min_lattice);

    for (j = 0; j < 3; j++) {   
      min_lattice[j][1] = tmp_matrix[j][i];
      min_lattice[j][i] = tmp_matrix[j][1];
    }
  }    

  /* Flip second axis */
  if (projection[i][0] < -symprec)
    for (j = 0; j < 3; j++)
      min_lattice[j][1] = -min_lattice[j][1];
    
  /*   Right-handed orientation */
  if (mat_get_determinant_d3(min_lattice) < -symprec*symprec*symprec) {

    /* Flip third axis */
    for (j = 0; j < 3; j++)
      min_lattice[j][2] = -min_lattice[j][2];
  }        


}

static void get_projection(double projection[3][3], const double min_lattice[3][3], const double lattice[3][3])
{
  double tmp_matrix[3][3];
  
  mat_transpose_matrix_d3(tmp_matrix, min_lattice);
  mat_multiply_matrix_d3(projection, tmp_matrix, lattice);
}

/***********/
/*  Cubic  */
/***********/
static int get_brv_cubic(Bravais *bravais, const double min_lattice[3][3],
			 const double symprec)
{
  int edge_equal[3];

  check_equal_edge(edge_equal, min_lattice, symprec);
  if (!(edge_equal[0] && edge_equal[1] && edge_equal[2]))
    return 0;

  /* Cubic-P */
  if (check_cubic(min_lattice, symprec)) {
    mat_copy_matrix_d3(bravais->lattice, min_lattice);
    bravais->centering = NO_CENTER;
    debug_print("cubic, no-centering\n");
    return 1;
  }

  /* Cubic-I */
  if (brv_cubic_I_center(bravais->lattice, min_lattice, symprec)) {
    bravais->centering = BODY;
    debug_print("cubic, I-center\n");
    return 1;
  }

  /* Cubic-F */
  if (brv_cubic_F_center(bravais->lattice, min_lattice, symprec)) {
    bravais->centering = FACE;
    debug_print("cubic, F-center\n");
    return 1;
  }

  return 0;
}
  
static int check_cubic(const double lattice[3][3], const double symprec)
{
  int angle_90[3], edge_equal[3];

  check_angle90(angle_90, lattice, symprec);
  check_equal_edge(edge_equal, lattice, symprec);
 
  if (angle_90[0] && angle_90[1] && angle_90[2] && 
      edge_equal[0] && edge_equal[1] && edge_equal[2]) {

    return 1;
  }

  return 0;
}

static int brv_cubic_F_center(double lattice[3][3], const double min_lattice[3][3],
			      const double symprec)
{
  const int relative_axes[22][3] = {
    {-1, 1, 1},
    { 1,-1, 1},
    { 1, 1,-1},
    { 1, 1, 1},
    { 0, 1, 1}, /* 5 */
    { 1, 0, 1},
    { 1, 1, 0},
    { 0, 1,-1},
    {-1, 0, 1},
    { 1,-1, 0}, /* 10 */
    { 2, 1, 1},
    { 1, 2, 1},
    { 1, 1, 2},
    { 2, 1,-1},
    {-1, 2, 1}, /* 15 */
    { 1,-1, 2},
    { 2,-1,-1},
    {-1, 2,-1},
    {-1,-1, 2},
    { 2,-1, 1}, /* 20 */
    { 1, 2,-1},
    {-1, 1, 2},
  };

  return brv_exhaustive_search(lattice, min_lattice, check_cubic, relative_axes,
			       22, FACE, symprec);
}

static int brv_cubic_I_center(double lattice[3][3], const double min_lattice[3][3],
			      const double symprec)
{
  const int relative_axes[6][3] = {
    { 0, 1, 1},
    { 1, 0, 1},
    { 1, 1, 0},
    { 0, 1,-1},
    {-1, 0, 1},
    { 1,-1, 0},
  };

  return brv_exhaustive_search(lattice, min_lattice, check_cubic, relative_axes,
			       6, BODY, symprec);
}

/****************/
/*  Tetragonal  */
/****************/
static int get_brv_tetra(Bravais *bravais, const double min_lattice[3][3],
			 const double symprec)
{
  int angle_90[3], edge_equal[3];

  check_angle90(angle_90, min_lattice, symprec);
  check_equal_edge(edge_equal, min_lattice, symprec);

  /* Tetragonal-P */
  if ((angle_90[0] && angle_90[1] && angle_90[2]) &&
      (edge_equal[0] || edge_equal[1] || edge_equal[2])) {
    if (brv_tetra_primitive(bravais->lattice, min_lattice, symprec)) {
      bravais->centering = NO_CENTER;
      debug_print("tetra, no-centering\n");
      return 1;
    }
  }

  /* Tetragonal-I */
  /* There are three patterns. */
  /* One or two or three primitive axes orient to the body center. */

  /* One */
  if ((angle_90[0] && edge_equal[0]) || (angle_90[1] && edge_equal[1])
      || (angle_90[2] && edge_equal[2])) {
    if (brv_tetra_one(bravais->lattice, min_lattice, symprec)) {
      debug_print("tetra1, I-center\n");
      bravais->centering = BODY;
      return 1;
    }
  }

  /* Three */
  /* All thses axes orient to the body center. */
  if (edge_equal[0] && edge_equal[1] && edge_equal[2]) {
    if (brv_tetra_three(bravais->lattice, min_lattice, symprec)) {
      bravais->centering = BODY;
      debug_print("tetra three, I-center\n");
      return 1;
    }
  }

  /* Two */
  if (edge_equal[0] || edge_equal[1] || edge_equal[2]) {
    if (brv_tetra_two(bravais->lattice, min_lattice, symprec)) {
      bravais->centering = BODY;
      debug_print("tetra two, I-center\n");
      return 1;
    }
  }
  return 0;
}

static int check_tetra(const double lattice[3][3], const double symprec)
{
  int angle_90[3], edge_equal[3];

  check_angle90(angle_90, lattice, symprec);
  check_equal_edge(edge_equal, lattice, symprec);
 
  if (angle_90[0] && angle_90[1] && angle_90[2] && edge_equal[2]) {
    debug_print("*** check_tetra ***\n");
    debug_print_matrix_d3(lattice);
    return 1;
  }

  return 0;
}

static int brv_tetra_primitive(double lattice[3][3], const double min_lattice[3][3],
			       const double symprec)
{
  const int relative_axes[3][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1}
  };

  return brv_exhaustive_search(lattice, min_lattice, check_tetra, relative_axes,
			       3, NO_CENTER, symprec);
}

static int brv_tetra_one(double lattice[3][3], const double min_lattice[3][3],
			 const double symprec)
{
  const int relative_axes[15][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 1, 1, 2},
    {-1, 1, 2},
    { 1,-1, 2},
    {-1,-1, 2},
    { 1, 2, 1},
    {-1, 2, 1},
    { 1, 2,-1},
    {-1, 2,-1},
    { 2, 1, 1},
    { 2,-1, 1},
    { 2, 1,-1},
    { 2,-1,-1}
  };

  return brv_exhaustive_search(lattice, min_lattice, check_tetra, relative_axes,
			       15, BODY, symprec);
}

static int brv_tetra_two(double lattice[3][3], const double min_lattice[3][3],
			 const double symprec)
{
  const int relative_axes[13][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 1, 1, 0},
    { 1,-1, 0},
    { 0, 1, 1},
    { 0, 1,-1},
    { 1, 0, 1},
    {-1, 0, 1},
    {-1, 1, 1},
    { 1,-1, 1},
    { 1, 1,-1},
    { 1, 1, 1},
  };

  return brv_exhaustive_search(lattice, min_lattice, check_tetra, relative_axes,
			       13, BODY, symprec);
}

static int brv_tetra_three(double lattice[3][3], const double min_lattice[3][3],
			   const double symprec)
{
  const int relative_axes[6][3] = {
    { 0, 1, 1},
    { 1, 0, 1},
    { 1, 1, 0},
    { 0, 1,-1},
    {-1, 0, 1},
    { 1,-1, 0},
  };

  return brv_exhaustive_search(lattice, min_lattice, check_tetra, relative_axes,
			       6, BODY, symprec);
}

/******************/
/*  Orthorhombic  */
/******************/
static int get_brv_ortho(Bravais *bravais, const double min_lattice[3][3],
			 const double symprec)
{
  Centering centering;

  /* orthorhombic-P */
  if (check_ortho(min_lattice, symprec)) {
    mat_copy_matrix_d3(bravais->lattice, min_lattice);
    bravais->centering = NO_CENTER;
    goto end;
  }

  /* orthorhombic-C (or A,B) */
  centering = brv_ortho_base_center(bravais->lattice, min_lattice, symprec);
  if (centering) {
    bravais->centering = centering;
    goto end;
  }

  /* orthorhombic-I */
  if (brv_ortho_I_center(bravais->lattice, min_lattice, symprec)) {
    bravais->centering = BODY;
    goto end;
  }

  /* orthorhombic-F */
  if (brv_ortho_F_center(bravais->lattice, min_lattice, symprec)) {
    bravais->centering = FACE;
    goto end;
  }

  /* Not found */
  return 0;

  /* Found */
 end:
  return 1;
}

static int check_ortho(const double lattice[3][3], const double symprec)
{
  int angle_90[3];

  check_angle90(angle_90, lattice, symprec);
 
  if (angle_90[0] && angle_90[1] && angle_90[2]) {
    return 1;
  }
  
  return 0;
}

static Centering brv_ortho_base_center(double lattice[3][3],
				       const double min_lattice[3][3],
				       const double symprec)
{
  const int relative_axes_one[15][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 1, 0, 2},
    {-1, 0, 2},
    { 0, 1, 2},
    { 0,-1, 2},
    { 0, 2, 1},
    { 0, 2,-1},
    { 1, 2, 0},
    {-1, 2, 0},
    { 2, 1, 0},
    { 2,-1, 0},
    { 2, 0, 1},
    { 2, 0,-1}
  };

  const int relative_axes_two[9][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 0, 1, 1},
    { 0, 1,-1},
    { 1, 0, 1},
    {-1, 0, 1},
    { 1, 1, 0},
    { 1,-1, 0},
  };


  /* One axis orients to the base center. */
  if (brv_exhaustive_search(lattice, min_lattice, check_ortho, relative_axes_one,
			    15, BASE, symprec))
    return get_base_center(lattice, min_lattice, symprec);


  /* Two axes orient to the base center. */
  if (brv_exhaustive_search(lattice, min_lattice, check_ortho, relative_axes_two,
			    9, BASE, symprec))
    return get_base_center(lattice, min_lattice, symprec);


  return NO_CENTER;
}

static int brv_ortho_I_center(double lattice[3][3], const double min_lattice[3][3],
			      const double symprec)
{
  /* Basically same as Tetra-I */

  const int relative_axes_one[15][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 1, 1, 2},
    {-1, 1, 2}, /*  5 */
    { 1,-1, 2},
    {-1,-1, 2},
    { 1, 2, 1},
    {-1, 2, 1},
    { 1, 2,-1}, /* 10 */
    {-1, 2,-1},
    { 2, 1, 1},
    { 2,-1, 1},
    { 2, 1,-1},
    { 2,-1,-1}  /* 15 */
  };

  const int relative_axes_two[13][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 1, 1, 0},
    { 1,-1, 0}, /*  5 */
    { 0, 1, 1},
    { 0, 1,-1},
    { 1, 0, 1},
    {-1, 0, 1},
    {-1, 1, 1}, /* 10 */
    { 1,-1, 1},
    { 1, 1,-1},
    { 1, 1, 1},
  };

  const int relative_axes_three[6][3] = {
    { 0, 1, 1},
    { 1, 0, 1},
    { 1, 1, 0},
    { 0, 1,-1},
    {-1, 0, 1}, /*  5 */
    { 1,-1, 0},
  };


  /* One axis orients to the I center. */
  if (brv_exhaustive_search(lattice, min_lattice, check_ortho, relative_axes_one,
			    15, BODY, symprec))
    return 1;

  /* Two axes orient to the I center. */
  if (brv_exhaustive_search(lattice, min_lattice, check_ortho, relative_axes_two,
			    13, BODY, symprec))
    return 1;

  /* Three axes orient to the I center. */
  if (brv_exhaustive_search(lattice, min_lattice, check_ortho, relative_axes_three,
			    6, BODY, symprec))
    return 1;

  return 0;
}

static int brv_ortho_F_center(double lattice[3][3], const double min_lattice[3][3],
			      const double symprec)
{
  /* At least two axes orient to the face center. */
  const int relative_axes_two[21][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 1, 1, 0},
    { 1,-1, 0}, /* 5 */
    { 0, 1, 1},
    { 0, 1,-1},
    { 1, 0, 1},
    {-1, 0, 1},
    { 1, 0, 2}, /* 10 */
    {-1, 0, 2},
    { 0, 1, 2},
    { 0,-1, 2},
    { 0, 2, 1},
    { 0, 2,-1}, /* 15 */
    { 1, 2, 0},
    {-1, 2, 0},
    { 2, 1, 0},
    { 2,-1, 0},
    { 2, 0, 1}, /* 20 */
    { 2, 0,-1}
  };

  const int relative_axes_three[22][3] = {
    {-1, 1, 1},
    { 1,-1, 1},
    { 1, 1,-1},
    { 1, 1, 1},
    { 0, 1, 1}, /* 5 */
    { 1, 0, 1},
    { 1, 1, 0},
    { 0, 1,-1},
    {-1, 0, 1},
    { 1,-1, 0}, /* 10 */
    { 2, 1, 1},
    { 1, 2, 1},
    { 1, 1, 2},
    { 2, 1,-1},
    {-1, 2, 1}, /* 15 */
    { 1,-1, 2},
    { 2,-1,-1},
    {-1, 2,-1},
    {-1,-1, 2},
    { 2,-1, 1}, /* 20 */
    { 1, 2,-1},
    {-1, 1, 2},
  };


  /* Two axes orient to the F center. */
  if (brv_exhaustive_search(lattice, min_lattice, check_ortho, relative_axes_two,
			    21, FACE, symprec))
    return 1;

  /* Three axes orient to the F center. */
  if (brv_exhaustive_search(lattice, min_lattice, check_ortho, relative_axes_three,
			    22, FACE, symprec))
    return 1;

  return 0;
}



/************************************/
/*  Hexagonal and Trigonal systems  */
/************************************/
static int get_brv_hexa(Bravais *bravais, const double min_lattice[3][3],
			const double symprec)
{
  const int relative_axes[5][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    {-1, 0, 0},
    { 0,-1, 0},
  };

  if (brv_exhaustive_search(bravais->lattice, min_lattice, check_hexa,
			    relative_axes, 5, NO_CENTER, symprec)) {
    bravais->centering = NO_CENTER;
    return 1;
  }

  return 0;
}

static int check_hexa(const double lattice[3][3], const double symprec)
{
  int angle_90[3], edge_equal[3];
  double ratio, metric[3][3];

  check_angle90(angle_90, lattice, symprec);
  check_equal_edge(edge_equal, lattice, symprec);

  if (angle_90[0] && angle_90[1] && edge_equal[2]) {
    get_metric(metric, lattice);
    ratio = metric[0][1] / metric[0][0];

    if ( mat_Dabs(ratio + 0.5) < symprec ) {
      return 1;
    }
  }

  return 0;
}

/*************************/
/*  Rhombohedral system  */
/*************************/
static int get_brv_rhombo(Bravais *bravais, const double min_lattice[3][3],
			  const double symprec)
{
  int edge_equal[3];
  check_equal_edge(edge_equal, min_lattice, symprec);

  /* One or two or three axes orient to the rhombochedral lattice points. */

  /* Three */
  if (edge_equal[0] && edge_equal[1] && edge_equal[2]) {
    if (brv_rhombo_three(bravais->lattice, min_lattice, symprec)) {
      bravais->centering = NO_CENTER;
      return 1;
    }
  }

  /* Two of three are in the base plane or Two are the rhombo axes. */
  if (edge_equal[0] || edge_equal[1] || edge_equal[2]) {
    if (brv_rhombo_two(bravais->lattice, min_lattice, symprec)) {
      bravais->centering = NO_CENTER;
      return 1;
    }
  }

  return 0;
}

static int check_rhombo(const double lattice[3][3], const double symprec)
{
  double metric[3][3];
  int edge_equal[3];

  get_metric(metric, lattice);
  check_equal_edge(edge_equal, lattice, symprec);

  if (edge_equal[0] && edge_equal[1] && edge_equal[2] &&
      (mat_Dabs((metric[0][1] - metric[1][2]) / metric[0][1]) < symprec * symprec) &&
      (mat_Dabs((metric[0][1] - metric[0][2]) / metric[0][1]) < symprec * symprec)) {
    return 1;
  }

  return 0;
}

static int brv_rhombo_two(double lattice[3][3], const double min_lattice[3][3], const double symprec)
{
  const int relative_axes_rhombo[9][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 0,-1, 0},
    { 0, 0,-1}, /*  5 */
    {-1, 1, 1},
    { 1,-1, 1},
    { 1, 1,-1},
    { 1, 1, 1}
  };

  const int relative_axes_base[13][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    {-1, 1, 1},
    { 1,-1, 1}, /*  5 */
    { 1, 1,-1},
    { 1, 1, 1},
    { 0, 1, 1},
    { 1, 0, 1},
    { 1, 1, 0}, /* 10 */
    { 0, 1,-1},
    {-1, 0, 1},
    { 1,-1, 0}
  };

  if (brv_exhaustive_search(lattice, min_lattice, check_rhombo,
			    relative_axes_rhombo, 9, NO_CENTER, symprec)) {
    return 1;
  }

  if (brv_exhaustive_search(lattice, min_lattice, check_rhombo, relative_axes_base,
			    13, NO_CENTER, symprec)) {
    return 1;
  }

  return 0;
}

static int brv_rhombo_three(double lattice[3][3], const double min_lattice[3][3], 
			    const double symprec)
{
  const int relative_axes[5][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 0,-1, 0},
    { 0, 0,-1},
  };

  return brv_exhaustive_search(lattice, min_lattice, check_rhombo, relative_axes,
			       5, NO_CENTER, symprec);
}

/***********************/
/*  Monoclinic system  */
/***********************/
static int get_brv_monocli(Bravais *bravais, const double min_lattice[3][3],
			   const double symprec)
{
  Centering centering;

  /* Monoclinic definition: */
  /*  second setting, i.e., alpha=gamma=90 deg. beta is not 90 deg. */

  /* Monoclinic-P */
  if (brv_monocli_primitive(bravais->lattice, min_lattice, symprec)) {
    bravais->centering = NO_CENTER;
    debug_print("Monoclinic-P\n");
    goto end;
  }

  /* Monoclinic base center*/
  centering = brv_monocli_base_center(bravais->lattice, min_lattice, symprec);
  if (centering) {
    bravais->centering = centering;
    debug_print("Monoclinic base center\n");
    set_brv_monocli(bravais, symprec);
    goto end;
  }

  return 0;

 end:
  return 1;
}


static int check_monocli(const double lattice[3][3], const double symprec)
{
  int angle_90[3];

  check_angle90(angle_90, lattice, symprec);

  if (angle_90[0] && angle_90[2]) {
    return 1;
  }

  /* Not found */
  return 0;
}

static void set_brv_monocli(Bravais *bravais, const double symprec)
{
  int i, angle_90[3];
  Centering centering;
  double tmp_lattice[3][3], lattice[3][3];

  mat_copy_matrix_d3(lattice, bravais->lattice);
  centering = bravais->centering;
  mat_copy_matrix_d3(tmp_lattice, lattice);
  check_angle90(angle_90, lattice, symprec);

  if (angle_90[1] && angle_90[2]) {
    for (i = 0; i < 3; i++) {
      lattice[i][0] =  tmp_lattice[i][1];
      lattice[i][1] = -tmp_lattice[i][0];
      lattice[i][2] =  tmp_lattice[i][2];
    }

    if (centering == B_FACE)
      centering = A_FACE;

  } else {
    if (angle_90[0] && angle_90[1]) {
      for (i = 0; i < 3; i++) {
	lattice[i][0] =  tmp_lattice[i][0];
	lattice[i][1] = -tmp_lattice[i][2];
	lattice[i][2] =  tmp_lattice[i][1];
      }

      if (centering == B_FACE)
	centering = C_FACE;
    }
  }

  mat_copy_matrix_d3(tmp_lattice, lattice);

  if (centering == A_FACE) {
    for (i = 0; i < 3; i++) {
      lattice[i][0] = -tmp_lattice[i][2];
      lattice[i][1] =  tmp_lattice[i][1];
      lattice[i][2] =  tmp_lattice[i][0];
    }
    centering = C_FACE;
  }

  if (centering == BODY) {

    for (i = 0; i < 3; i++) {
      lattice[i][0] = tmp_lattice[i][0] + tmp_lattice[i][2];
      lattice[i][1] = tmp_lattice[i][1];
      lattice[i][2] = tmp_lattice[i][2];
    }
    centering = C_FACE;
  }

  mat_copy_matrix_d3(bravais->lattice, lattice);
  bravais->centering = centering;
}

static int brv_monocli_primitive(double lattice[3][3],
				 const double min_lattice[3][3],
				 const double symprec)
{
  const int relative_axes[3][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1}
  };

  return brv_exhaustive_search(lattice, min_lattice, check_monocli, relative_axes,
			       3, NO_CENTER, symprec);
}

static Centering brv_monocli_base_center(double lattice[3][3],
					 const double min_lattice[3][3],
					 const double symprec)
{
  int found;
  Centering centering = NO_CENTER;
  
  int relative_axes_one[15][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 1, 0, 2},
    {-1, 0, 2}, /*  5 */
    { 0, 1, 2},
    { 0,-1, 2},
    { 0, 2, 1},
    { 0, 2,-1},
    { 1, 2, 0}, /* 10 */
    {-1, 2, 0},
    { 2, 1, 0},
    { 2,-1, 0},
    { 2, 0, 1},
    { 2, 0,-1}  /* 15 */
  };

  int relative_axes_two_three[25][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 0, 1, 1},
    { 0, 1,-1}, /*  5 */
    { 1, 0, 1},
    {-1, 0, 1},
    { 1, 1, 0},
    { 1,-1, 0},
    { 1, 0, 2}, /* 10 */
    {-1, 0, 2},
    { 0, 1, 2},
    { 0,-1, 2},
    { 0, 2, 1},
    { 0, 2,-1}, /* 15 */
    { 1, 2, 0},
    {-1, 2, 0},
    { 2, 1, 0},
    { 2,-1, 0},
    { 2, 0, 1}, /* 20 */
    { 2, 0,-1},
    {-1, 1, 1},
    { 1,-1, 1},
    { 1, 1,-1},
    { 1, 1, 1}  /* 25 */
  };

  /* One axis orients to the base center. */
  found = brv_exhaustive_search(lattice, min_lattice, check_monocli, relative_axes_one,
				15, BASE, symprec);

  /* Two or three axes orient to the base center. */
  if (!found) {
   found = brv_exhaustive_search(lattice, min_lattice, check_monocli,
				 relative_axes_two_three, 25, BASE, symprec);
  }

  if (found) {
    brv_smallest_lattice_vector(lattice, lattice, symprec);
    centering = get_base_center(lattice, min_lattice, symprec);
    debug_print("Monocli centering: %d\n",  centering);
    return centering;
  }

  return NO_CENTER;
}


/*******************************/
/* The other private functions */
/*******************************/
static void check_angle90(int angle_90[3], const double lattice[3][3],
			  const double symprec)
{
  int i, c0, c1;
  double metric[3][3];

  get_metric(metric, lattice);
  for (i = 0; i < 3; i++) {

    c0 = (i + 1) % 3;
    c1 = (i + 2) % 3;

    if (metric[c0][c1] * metric[c1][c0] / metric[c0][c0] / metric[c1][c1]
	< symprec * 2) {	/* orthogonal */
      angle_90[i] = 1;
    }
    else {
      angle_90[i] = 0;
    }
  }
}

static void check_equal_edge(int edge_equal[3], const double lattice[3][3],
			     const double symprec)
{
  int i, c0, c1;
  double metric[3][3];

  get_metric(metric, lattice);
  for (i = 0; i < 3; i++) {

    c0 = (i + 1) % 3;
    c1 = (i + 2) % 3;

    if (mat_Dabs((metric[c0][c0] - metric[c1][c1]) / metric[c0][c0]) < symprec * 2 ) {	/* Equal edges */
      edge_equal[i] = 1;
    }
    else {
      edge_equal[i] = 0;
    }
  }
}

static void get_metric(double metric[3][3], const double lattice[3][3])
{
  double lattice_t[3][3];
  mat_transpose_matrix_d3(lattice_t, lattice);
  mat_multiply_matrix_d3(metric, lattice_t, lattice);
}

static int brv_exhaustive_search(double lattice[3][3], const double min_lattice[3][3],
				 int (*check_bravais)(const double lattice[3][3], const double symprec),
				 const int relative_axes[][3], const int num_axes,
				 const Centering centering, const double symprec)
{
  int i, j, k, l, factor, coordinate[3][3];
  double tmp_matrix[3][3];

  switch (centering) {
  case NO_CENTER:
    factor = 1;
    break;
  case BODY:
  case BASE:
  case A_FACE:
  case B_FACE:
  case C_FACE:
    factor = 2;
    break;
  case FACE:
    factor = 4;
    break;
  }

  for (i = 0; i < num_axes; i++) {
    for (j = 0; j < num_axes; j ++) {
      for (k = 0; k < num_axes; k++) {

	for (l = 0; l < 3; l++) {
	  coordinate[l][0] = relative_axes[i][l];
	  coordinate[l][1] = relative_axes[j][l];
	  coordinate[l][2] = relative_axes[k][l];
	}
	    
	if (mat_Dabs(mat_get_determinant_i3(coordinate)) == factor) { 

	  mat_cast_matrix_3i_to_3d(tmp_matrix, coordinate);
	  mat_multiply_matrix_d3(lattice, min_lattice, tmp_matrix);

	  if ((*check_bravais)(lattice, symprec) &&
	      mat_Dabs(mat_get_determinant_d3(lattice)) > symprec)
	    return 1;
	}
      }
    }
  }

  return 0;
}

static Centering get_base_center(const double brv_lattice[3][3],
				 const double min_lattice[3][3],
				 const double symprec) {
  int i;
  Centering centering = NO_CENTER;
  double tmp_matrix[3][3], axis[3][3];

  if (mat_inverse_matrix_d3(tmp_matrix, brv_lattice, symprec)) {
    mat_multiply_matrix_d3(axis, tmp_matrix, min_lattice);
  } else {
    fprintf(stderr, "spglib BUG: %s in %s\n", __FUNCTION__, __FILE__);
    return NO_CENTER;
  }

  debug_print_matrix_d3(min_lattice);
  debug_print_matrix_d3(brv_lattice);
  debug_print("%f %f %f\n", axis[0][0], axis[0][1], axis[0][2]);
  debug_print("%f %f %f\n", axis[1][0], axis[1][1], axis[1][2]);
  debug_print("%f %f %f\n", axis[2][0], axis[2][1], axis[2][2]);

  /* C center */
  for (i = 0; i < 3; i++) {
    if ((mat_Dabs(mat_Dabs(axis[0][i]) - 0.5) < symprec) &&
	(mat_Dabs(mat_Dabs(axis[1][i]) - 0.5) < symprec) &&
	(!(mat_Dabs(mat_Dabs(axis[2][i]) - 0.5) < symprec))) {
      centering = C_FACE;
      goto end;
    }
  }

  /* A center */
  for (i = 0; i < 3; i++) {
    if ((!(mat_Dabs(mat_Dabs(axis[0][i]) - 0.5) < symprec)) &&
	(mat_Dabs(mat_Dabs(axis[1][i]) - 0.5) < symprec) &&
	(mat_Dabs(mat_Dabs(axis[2][i]) - 0.5) < symprec)) {
      centering = A_FACE;
      goto end;
    }
  }

  /* B center */
  for (i = 0; i < 3; i++) {
    if ((mat_Dabs(mat_Dabs(axis[0][i]) - 0.5) < symprec) &&
	(!(mat_Dabs(mat_Dabs(axis[1][i]) - 0.5) < symprec)) &&
	(mat_Dabs(mat_Dabs(axis[2][i]) - 0.5) < symprec)) {
      centering = B_FACE;
      goto end;
    }
  }

  /* body center */
  for (i = 0; i < 3; i++) {
    if ((mat_Dabs(mat_Dabs(axis[0][i]) - 0.5) < symprec) &&
	(mat_Dabs(mat_Dabs(axis[1][i]) - 0.5) < symprec) &&
	(mat_Dabs(mat_Dabs(axis[2][i]) - 0.5) < symprec)) {
      centering = BODY;
      goto end;
    }
  }

  /* This should not happen. */
  fprintf(stderr, "spglib BUG: %s in %s\n", __FUNCTION__, __FILE__);
  return NO_CENTER;

 end:
  return centering;
}

static void get_right_hand_lattice(double lattice[3][3], const double symprec)
{
  int i;

  if (mat_get_determinant_d3(lattice) < symprec) {
    for (i = 0; i < 3; i++) {
      lattice[i][0] = -lattice[i][0];
      lattice[i][1] = -lattice[i][1];
      lattice[i][2] = -lattice[i][2];
    }
  }
}
