#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_spline.h>

#include "config.h"

// Fortran does not have the asinh intrinsic, 
// so we use the one from libm.a
double PROTO(asinh)(double *x)
{
  return asinh(*x);
}

// complementary error function
// we use the one in gsl
double PROTO(erfc)(double *x)
{
	return gsl_sf_erfc(*x);
}

// error function
// we use the one in gsl
double PROTO(erf)(double *x)
{
	return gsl_sf_erf(*x);
}

// normalization constants for the spherical harmonics
static double sph_cnsts[9] = {
	0.282094791773878, // l = 0
	0.488602511902920, 0.488602511902920, 0.488602511902920, // l = 1
	0.182091405098680, 0.364182810197360, 0.630783130505040, 0.364182810197360, 0.182091405098680
};

// Computes real spherical harmonics ylm in the direction of vector r:
//    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
//    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
// with (theta,phi) the polar angles of r, c a positive normalization
// constant and plm associated legendre polynomials.
double PROTO(ylm)(double *x, double *y, double *z, int *l, int *m)
{
	double r, rx, ry, rz, cosphi, sinphi, cosm, sinm, phase;
	int i;

	if(*l == 0)
		return sph_cnsts[0];

	r = sqrt((*x)*(*x) + (*y)*(*y) + (*z)*(*z));

	// if r=0, direction is undefined => make ylm=0 except for l=0
	if(r == 0.) return 0.;

	rx = (*x)/r; ry = (*y)/r; rz = (*z)/r;

	switch(*l){
	case 1:
		switch(*m){
		case -1: return -sph_cnsts[1]*ry;
		case  0: return  sph_cnsts[2]*rz;
		case  1: return -sph_cnsts[3]*rx;
		}
	case 2:
		switch(*m){
		case -2: return  sph_cnsts[4]*6*rx*ry;
		case -1: return -sph_cnsts[5]*3*ry*rz;
		case  0: return  sph_cnsts[6]*0.5*(3*rz*rz - 1);
		case  1: return -sph_cnsts[7]*3*rx*rz;
		case  2: return  sph_cnsts[8]*3*(rx*rx - ry*ry);
		}
	}

	// get phase
	r = sqrt(rx*rx + ry*ry);
	if(fabs(r) < 1e-20) // one never knows...
		r = 1e-20;
	cosphi = rx/r; sinphi = ry/r;
		
	// compute sin(mphi) and cos(mphi) by adding cos/sin
	cosm = 1.; sinm=0.;
	for(i=0; i<abs(*m); i++){
		double a = cosm, b = sinm;
		cosm = a*cosphi - b*sinphi;
		sinm = a*sinphi + b*cosphi;
	}
	phase = *m<0 ? sinm : cosm;

	r = gsl_sf_legendre_sphPlm(*l, abs(*m), rz);

	return r*phase;
}

/* Interface to the GSL interpolation functions */
void PROTO(spline_end)(void **spl, void **acc)
{
	gsl_spline_free((gsl_spline *)(*spl));
	gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

void PROTO(spline_fit)(int *nrc, double *x, double *y, void **spl, void **acc)
{
	*acc = (void *)gsl_interp_accel_alloc();
	*spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
	gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
	fflush(stdout);
}

double PROTO(spline_eval)(double *x, void **spl, void **acc)
{
	return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}


