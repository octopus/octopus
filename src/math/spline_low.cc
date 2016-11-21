/*
 Copyright (C) 2016 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira, X. Andrade

// Spline routines produced at the Lawrence Livermore National
// Laboratory.  Written by Xavier Andrade (xavier@llnl.gov), Erik
// Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_spline.h>

#include "string_f.h"

#include <fortran_types.h>

#include <assert.h>
#include <iostream>
#include <vector>

#define SPLINE_FLAT_BC 0.0       /* Flat boundary condition (y'=0) */
#define SPLINE_NATURAL_BC 1.e31  /* Natural boundary condition (Y"=0) */

#define GSL    1
#define NATIVE 2

static int library;

static void spline(const double *x, const double *y, int n, double yp1, double ypn, double *y2)
{

  int i,k;
  double p,qn,sig,un,*u = new double[n];

  if ( yp1 >= 1.e30 )
  {
    y2[0] = 0.0;
    u[0] = 0.0;
  }
  else
  {
    y2[0] = -0.5;
    assert ( x[1] - x[0] > 0.0 );
    u[0] = ( 3.0 / (x[1]-x[0]) ) * ( (y[1]-y[0]) / (x[1]-x[0]) - yp1 );
  }

  for ( i = 1; i < n-1; i++ )
  {
    assert ( x[i+1] > x[i] );
    sig = ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] );
    p = sig * y2[i-1] + 2.0;
    y2[i] = ( sig - 1.0 ) / p;
    u[i] = ( 6.0 * ( ( y[i+1] - y[i] ) / ( x[i+1] - x[i] ) -
                     ( y[i] - y[i-1] ) / ( x[i] - x[i-1] ) ) /
             ( x[i+1] - x[i-1] ) - sig * u[i-1] ) / p;
  }

  if ( ypn >= 1.e30 )
  {
    qn = 0.0;
    un = 0.0;
  }
  else
  {
    qn = 0.5;
    un = ( 3.0 / (x[n-1]-x[n-2]) ) * 
         ( ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]) );
  }

  y2[n-1] = ( un - qn * u[n-2] ) / ( qn * y2[n-2] + 1.0 );

  for ( k = n-2; k >= 0; k-- )
  {
    y2[k] = y2[k] * y2[k+1] + u[k];
  }

  delete [] u;
}

static void splint (const double *xa, const double *ya, const double *y2a, int n, double x, double *y)
{
  int k,khi,klo;
  double a,b,h;

  klo = 0;
  khi = n-1;

  while ( khi - klo > 1 )
  {
    k = ( khi + klo ) / 2;
    if ( xa[k] > x )
      khi = k;
    else
      klo = k;
  }

  //ewd DEBUG
  if (khi > n-1) {
    std::cout << "ERROR.SPLINT:  khi = " << khi << ", n = " << n << std::endl;
    return;
  }
  if (klo > n-1) {
    std::cout << "ERROR.SPLINT:  klo = " << klo << ", n = " << n << std::endl;
    return;
  }
  
  h = xa[khi] - xa[klo];
  assert ( h > 0.0 );

  a = ( xa[khi] - x ) / h;
  b = ( x - xa[klo] ) / h;

  *y = a * ya[klo] + b * ya[khi] + h * h * (1.0/6.0) *
       ( (a*a*a-a) * y2a[klo] + (b*b*b-b) * y2a[khi] );

}

static void splintd (const double *xa, const double *ya, const double *y2a,
 int n, double x, double *y, double *dy)
{
  int k,khi,klo;
  double a,b,h;

  klo = 0;
  khi = n-1;

  while ( khi - klo > 1 )
  {
    k = ( khi + klo ) / 2;
    if ( xa[k] > x )
      khi = k;
    else
      klo = k;
  }

  h = xa[khi] - xa[klo];
  assert ( h > 0.0 );

  a = ( xa[khi] - x ) / h;
  b = ( x - xa[klo] ) / h;

  *y = a * ya[klo] + b * ya[khi] + h * h * (1.0/6.0) *
       ( (a*a*a-a) * y2a[klo] + (b*b*b-b) * y2a[khi] );

  *dy = ( ya[khi] - ya[klo] ) / h +
        h * ( ( (1.0/6.0) - 0.5 * a * a ) * y2a[klo] +
              ( 0.5 * b * b - (1.0/6.0) ) * y2a[khi] );
}

class Spline {

 public:

  Spline(){
  }

  const int size() const {
    return x_.size();
  }
  
  void fit(const double *x, const double *y, const int n, const double yp1, const double ypn){
    x_.resize(n);
    y_.resize(n);
    y2_.resize(n);
    
    for(int ii = 0; ii < n; ii++){
      x_[ii] = x[ii];
      y_[ii] = y[ii];
    }
    spline(x, y, n, yp1, ypn, &y2_[0]);
  }

  void get_x(double * x) const {
    for(int ii = 0; ii < size(); ii++) x[ii] = x_[ii];
  }

  void get_y(double * y) const {
    for(int ii = 0; ii < size(); ii++) y[ii] = y_[ii];
  }

  double value(const double & x) const {
    double y;
    splint(&x_[0], &y_[0], &y2_[0], x_.size(), x, &y);
    return y;
  }
  
  double derivative(const double & x) const {
    double y, dy;
    splintd(&x_[0], &y_[0], &y2_[0], x_.size(), x, &y, &dy);
    return dy;
  }
  
  void derivative(const double & x, double & y, double & dy) const {
    splintd(&x_[0], &y_[0], &y2_[0], x_.size(), x, &y, &dy);
  }

  double integral() const {
    double sum = 0.5*(x_[1] - x_[0])*y_[0];
    for(int ii = 1; ii < size() - 2; ii++){
      double dx = 0.5*(x_[ii + 1] - x_[ii - 1]);
      sum += dx*y_[ii];
    }
    sum += 0.5*(x_[size() - 1] - x_[size() - 2])*y_[size() - 1];    
    return sum;
  }
  
 private :

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> y2_;
  
};

/* Interpolation */
extern "C" void FC_FUNC_(oct_spline_end, OCT_SPLINE_END)
     (void **spl, void **acc)
{
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    delete spline;
  } else {
    gsl_spline_free((gsl_spline *)(*spl));
    gsl_interp_accel_free((gsl_interp_accel *)(*acc));
  }
}

extern "C" void FC_FUNC_(oct_spline_fit, OCT_SPLINE_FIT)
  (const fint *nrc, const double *x, const double *y, void **spl, void **acc, const fint * lib){
  library = *lib;
  if(*lib == NATIVE){
    Spline * spline = new Spline();
    spline->fit(x, y, *nrc, SPLINE_NATURAL_BC, SPLINE_FLAT_BC);    
    *spl = (void *) spline;
  } else {
    /* the GSL headers actually specify size_t instead of const int for nrc */
    *acc = (void *)gsl_interp_accel_alloc();
    *spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
    gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
    fflush(stdout);
  }
}

extern "C" double FC_FUNC_(oct_spline_eval, OCT_SPLINE_EVAL)
     (const double *x, const void **spl, void **acc)
{
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    return spline->value(*x);
  } else {
    /* the GSL headers specify double x instead of const double x */
    return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
  }
}


template <typename Type, int stride>
void oct_spline_eval_array(fint nn, Type * xf, const void **spl, void **acc){
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    for(fint ii = 0; ii < nn; ii++) xf[ii] = spline->value(xf[stride*ii]);
  } else {
    for(fint ii = 0; ii < nn; ii++){
      xf[stride*ii] = Type(gsl_spline_eval((gsl_spline *)(*spl), xf[stride*ii], (gsl_interp_accel *)(*acc)));
    }
  }
}

extern "C" void FC_FUNC_(oct_spline_eval_array, OCT_SPLINE_EVAL_ARRAY)
     (const fint * nn, double *xf, const void **spl, void **acc)
{
  oct_spline_eval_array<double, 1>(*nn, xf, spl, acc);
}

extern "C" void FC_FUNC_(oct_spline_eval_array4, OCT_SPLINE_EVAL_ARRAY4)
     (const fint * nn, float *xf, const void **spl, void **acc)
{
  oct_spline_eval_array<float, 1>(*nn, xf, spl, acc);
}

/* use a stride of 2 to store into just the real part of a Fortran complex array */
extern "C" void FC_FUNC_(oct_spline_eval_arrayz, OCT_SPLINE_EVAL_ARRAYZ)
  (const fint * nn, double *xf, const void **spl, void **acc)
{
  oct_spline_eval_array<double, 2>(*nn, xf, spl, acc);
}

/* use a stride of 2 to store into just the real part of a Fortran complex array */
extern "C" void FC_FUNC_(oct_spline_eval_arrayc, OCT_SPLINE_EVAL_ARRAYC)
  (const fint * nn, float *xf, const void **spl, void **acc)
{
  oct_spline_eval_array<float, 2>(*nn, xf, spl, acc);
}

/* This function returns the number of points with which a spline
	 was constructed (the size component of the gsl_spline struct). */
extern "C" fint FC_FUNC_(oct_spline_npoints, OCT_SPLINE_NPOINTS)
  (const void **spl, void **acc)
{
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    return spline->size();
  } else {
    return (fint)((gsl_spline *)(*spl))->size;
  }
}

/* This function places in the x array the x values of a given spline spl*/ 
extern "C" void FC_FUNC_(oct_spline_x, OCT_SPLINE_X)
  (const void **spl, void **acc, double *x)
{
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    spline->get_x(x);
  } else {
    int size, i;
	
    size = (int)((gsl_spline *)(*spl))->size;
    for(i=0; i<size; i++)
      x[i] = ((gsl_spline *)(*spl))->x[i];
  }
}

/* This function places in the y array the y values of a given spline spl*/ 
extern "C" void FC_FUNC_(oct_spline_y, OCT_SPLINE_Y)
  (const void **spl, void **acc, double *y)
{
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    spline->get_y(y);  
  } else {
   int size, i;
	
   size = (int)((gsl_spline *)(*spl))->size;
   for(i=0; i<size; i++)
     y[i] = ((gsl_spline *)(*spl))->y[i];
  }
}

/* Returns the integral of the spline stored in spl, between a and b */
extern "C" double FC_FUNC_(oct_spline_eval_integ, OCT_SPLINE_EVAL_INTEG)
     (const void **spl, const double *a, const double *b, void **acc)
{
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    const int nn = spline->size()*2;
    double dx = fabs(*b - *a)/(nn - 1);
    
    double integral = (spline->value(*a) + spline->value(*b))*0.5;
    
    double x = *a + dx;
    for(int ii = 0; ii < nn - 2; ii++){
      integral += spline->value(x);
      x += dx;
    }
    
    integral *= dx;
    return integral;
  
  } else {
    /* the GSL headers specify double a, double b */
    return gsl_spline_eval_integ((gsl_spline *)(*spl), *a, *b, (gsl_interp_accel *)(* acc));
  }
}

extern "C" double FC_FUNC_(oct_spline_eval_integ_full, OCT_SPLINE_EVAL_INTEG_FULL)
     (const void **spl, void **acc)
{
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    return spline->integral();
  } else {
    /* the GSL headers specify double a, double b */
    const int size = (int)((gsl_spline *)(*spl))->size;
    const double a = ((gsl_spline *)(*spl))->x[0];
    const double b = ((gsl_spline *)(*spl))->x[size - 1];
    return gsl_spline_eval_integ((gsl_spline *)(*spl), a, b, (gsl_interp_accel *)(* acc));
  }
}

/* Performs the derivative of a spline */
extern "C" double FC_FUNC_(oct_spline_eval_der, OCT_SPLINE_EVAL_DER)
     (const double *x, const void **spl, void **acc)
{
  if(library == NATIVE){
    Spline * spline = (Spline *)*spl;
    return spline->derivative(*x);  
  } else {
    /* the GSL headers specify double x */
    return gsl_spline_eval_deriv((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
  }
}

/* Performs the second derivative of a spline */
extern "C" double FC_FUNC_(oct_spline_eval_der2, OCT_SPLINE_EVAL_DER2)
     (const double *x, const void **spl, void **acc)
{
  /* the GSL headers specify double x */
  return gsl_spline_eval_deriv2((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}
