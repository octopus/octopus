/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id$
*/

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_sf_legendre.h>

/* normalization constants for the spherical harmonics */
static double sph_cnsts[9] = {
	0.282094791773878, /* l = 0 */
	0.488602511902920, 0.488602511902920, 0.488602511902920, /* l = 1 */
	0.182091405098680, 0.364182810197360, 0.630783130505040, 0.364182810197360, 0.182091405098680
};

/* Computes real spherical harmonics ylm in the direction of vector r:
	 ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
	 ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
	 with (theta,phi) the polar angles of r, c a positive normalization
	 constant and plm associated legendre polynomials. */

double FC_FUNC_(oct_ylm, OCT_YLM)
     (const double *x, const double *y, const double *z, const int *l, const int *m)
{
  double r, r2, rr, rx, ry, rz, cosphi, sinphi, cosm, sinm, phase;
  int i;

  if(l[0] == 0) return sph_cnsts[0];

  r2 = x[0]*x[0] + y[0]*y[0] + z[0]*z[0];

  /* if r=0, direction is undefined => make ylm=0 except for l=0 */
  if(r2 < 1.0e-15) return 0.0;

  switch(l[0]){
  case 1:
    rr = 1.0/sqrt(r2);
    switch(m[0]){
    case -1: return -sph_cnsts[1]*rr*y[0];
    case  0: return  sph_cnsts[2]*rr*z[0];
    case  1: return -sph_cnsts[3]*rr*x[0];
    }
  case 2:
    switch(m[0]){
    case -2: return  sph_cnsts[4]*6.0*x[0]*y[0]/r2;
    case -1: return -sph_cnsts[5]*3.0*y[0]*z[0]/r2;
    case  0: return  sph_cnsts[6]*0.5*(3.0*z[0]*z[0]/r2 - 1.0);
    case  1: return -sph_cnsts[7]*3.0*x[0]*z[0]/r2;
    case  2: return  sph_cnsts[8]*3.0*(x[0]*x[0] - y[0]*y[0])/r2;
    }
  }

  /* get phase */
  rr = 1.0/sqrt(r2);
  
  rx = x[0]*rr;
  ry = y[0]*rr;
  rz = z[0]*rr;

  r = hypot(rx, ry);
  if(r < 1e-20) r = 1e-20; /* one never knows... */

  cosphi = rx/r;
  sinphi = ry/r;

  /* compute sin(mphi) and cos(mphi) by adding cos/sin */
  cosm = 1.; sinm=0.;
  for(i = 0; i < abs(m[0]); i++){
    double a = cosm, b = sinm;
    cosm = a*cosphi - b*sinphi;
    sinm = a*sinphi + b*cosphi;
  }
  phase = m[0] < 0 ? sinm : cosm;
  phase = m[0] == 0 ? phase : sqrt(2.0)*phase;
	
  /* adding small number (~= 10^-308) to avoid floating invalids */
  rz = rz + DBL_MIN;

  r = gsl_sf_legendre_sphPlm(l[0], abs(m[0]), rz);

  /* I am not sure whether we are including the Condon-Shortley factor (-1)^m */
  return r*phase;
}

