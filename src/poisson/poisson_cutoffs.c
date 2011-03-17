#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <assert.h>
#include <config.h>

#define M_EULER_MASCHERONI 0.5772156649015328606065120

/*
This file contains the definition of functions used by Fortran module
poisson_cutoff_m. There are four functions to be defined:

c_poisson_cutoff_3d_1d_finite: For the cutoff of 3D/0D cases, in which however
                               we want the region of the cutoff to be a cylinder.
c_poisson_cutoff_2d_0d:

c_poisson_cutoff_2d_1d:

*/



/************************************************************************/
/************************************************************************/
/* C_POISSON_CUTOFF_3D_1D_FINITE */
/************************************************************************/
/************************************************************************/

struct parameters
{ 
  double kx; 
  double kr; 
  double x0; 
  double r0; 
  int index; 
}; 

/* the index is the type of integral. This is necessary because depending of the value of k the
   the integral needs to be performed in different ways */
#define PFC_F_0   0
#define PFC_F_KX  1
#define PFC_F_KR  2

static const double tol = 1e-7;
static double zero = 1e-100;

static double f(double w, void *p)
{
  struct parameters *params = (struct parameters *)p;
  const double  kx = (params->kx);
  double  kr = (params->kr);
  double  x0 = (params->x0);
  double  r0 = (params->r0);
  int index  = (params->index);
  double result;

  if (w == 0) w = zero;

  if (fabs(kx - w) > tol){ 
    result = sin((kx - w)*x0)/(kx - w)*
      sqrt(pow(w, 2.0*index)*pow(r0, 2.0*index))/(kr*kr + w*w)
      *gsl_sf_bessel_Kn(index, r0*fabs(w));
  }
  else{
    result = x0*sqrt(pow(w, 2.0*index))/(kr*kr + w*w)*gsl_sf_bessel_Kn(index, r0*fabs(w));
  }
  return result;
}


static double f_kx_zero(double w, void *p)
{
  struct parameters *params = (struct parameters *)p;
  double  kr = (params->kr);
  double  x0 = (params->x0);
  double result;
  
  if (w == 0.0) w = zero;
  
  result = 2.0*M_PI*w*gsl_sf_bessel_J0(kr*w)*
    (log(x0 + sqrt(w*w + x0*x0)) - log(-x0 + sqrt(w*w + x0*x0)));
  return result;
}


static double f_kr_zero(double w, void *p)
{
  struct parameters *params = (struct parameters *)p;
  double  kx = (params->kx);
  double  r0 = (params->r0);
  double result;

  result = 4.0*M_PI*cos(kx*w)*(sqrt(r0*r0 + w*w) - fabs(w));
  return result;
}


static double int_aux(int index, double kx, double kr, double x0, double r0, int which_int)
{
  /*variables*/
  double result, error;
  double xinf = 100.0/r0;
  struct parameters params;

  /*pass the given parameters to program*/
  const size_t wrk_size = 5000;
  /* I believe that 1e-6 is over-killing, a substantial gain in time
     is obtained by setting this threshold to 1e-3 */
  /* const double epsilon_abs = 1e-6; */
  /* const double epsilon_rel = 1e-6; */
  const double epsilon_abs = 1e-3;
  const double epsilon_rel = 1e-3;

  /*creating the workspace*/
  gsl_integration_workspace * ws = gsl_integration_workspace_alloc (wrk_size);

  /*create the function*/
  gsl_function F;

  assert(which_int == PFC_F_0 || which_int == PFC_F_KR || which_int == PFC_F_KX);

  params.kx = kx;
  params.kr = kr;
  params.x0 = x0;
  params.r0 = r0;
  params.index = index;

  F.params = &params;
      
  /*select the integration method*/
  switch(which_int){
  case PFC_F_0:
    F.function = &f;
    gsl_integration_qag(&F, -xinf, xinf, epsilon_abs, epsilon_rel, 
			wrk_size, 3, ws, &result, &error);
    break;
  case PFC_F_KX:
    F.function = &f_kx_zero;
    gsl_integration_qag (&F, 0.0, r0, epsilon_abs, epsilon_rel, 
			 wrk_size, 3, ws, &result, &error);
    break;
  case PFC_F_KR:
    F.function = &f_kr_zero;
    gsl_integration_qag (&F, 0.0, x0, epsilon_abs, epsilon_rel, 
			 wrk_size, 3, ws, &result, &error);
    break;
  }

  /*free the integration workspace*/
  gsl_integration_workspace_free (ws);
  
  /*return*/
  return result;
}


double poisson_finite_cylinder(double kx, double kr, double x0,double r0) {
  double result;

  if ((kx>=tol) & (kr>=tol)) {
    result = 4.0*kr*r0*gsl_sf_bessel_Jn(1, kr*r0)*int_aux(0, kx, kr, x0, r0, PFC_F_0)
      - 4.0*gsl_sf_bessel_Jn(0, kr*r0)*int_aux(1, kx, kr, x0, r0, PFC_F_0)
      + 4.0*M_PI/(kx*kx + kr*kr)*(1.0 + exp(-kr*x0)*(kx/kr*sin(kx*x0) - cos(kx*x0)));
  }
  else if ((kx < tol) & (kr >= tol)){
    result = int_aux(0, kx, kr, x0, r0, PFC_F_KX);
  }
  else if ((kx >= tol) & (kr < tol)){
    result = int_aux(0, kx, kr, x0, r0, PFC_F_KR);
  }
  else 
    result = -2.0*M_PI*( log(r0/(x0 + sqrt(r0*r0 + x0*x0)))
			 *r0*r0 + x0*(x0 - sqrt(r0*r0 + x0*x0)));
  
  /* the 1/(4pi) factor is due to the factor on the poissonsolver3d (end)*/
  return result/(4.0*M_PI);
}

/* --------------------- Interface to Fortran ---------------------- */
double FC_FUNC_(c_poisson_cutoff_3d_1d_finite, C_POISSON_CUTOFF_3D_1D_FINITE)
  (double *gx, double *gperp, double *xsize, double *rsize)
{
  return poisson_finite_cylinder(*gx, *gperp, *xsize, *rsize);
}





/************************************************************************/
/************************************************************************/
/* C_POISSON_CUTOFF_2D_0D */
/************************************************************************/
/************************************************************************/
static double bessel_J0(double w, void *p)
{
  return gsl_sf_bessel_J0(w);
}

/* --------------------- Interface to Fortran ---------------------- */
double FC_FUNC_(c_poisson_cutoff_2d_0d, C_POISSON_CUTOFF_2D_0D)
  (double *x, double *y)
{
  double result, error;
  const size_t wrk_size = 500;
  const double epsilon_abs = 1e-3;
  const double epsilon_rel = 1e-3;

  gsl_integration_workspace * ws = gsl_integration_workspace_alloc (wrk_size);
  gsl_function F;

  F.function = &bessel_J0;
  gsl_integration_qag(&F, *x, *y, epsilon_abs, epsilon_rel, 
  			wrk_size, 3, ws, &result, &error);
  gsl_integration_workspace_free(ws);

  return result;
}





/************************************************************************/
/************************************************************************/
/* INTCOSLOG */
/************************************************************************/
/************************************************************************/
/* Int_0^mu dy cos(a*y) log(abs(b*y)) = 
   (1/a) * (log(|b*mu|)*sin(a*mu) - Si(a*mu) )    */
double FC_FUNC(intcoslog, INTCOSLOG)(double *mu, double *a, double *b)
{
  if(fabs(*a)>0.0){
    return (1.0/(*a)) * (log((*b)*(*mu))*sin((*a)*(*mu)) - gsl_sf_Si((*a)*(*mu)) );
  }else{
    return (*mu)*(log((*mu)*(*b)) - 1.0);
  }
}





/************************************************************************/
/************************************************************************/
/* C_POISSON_CUTOFF_2D_1D */
/************************************************************************/
/************************************************************************/
struct parameters_2d_1d
{ 
  double gx;
  double gy;
  double rc; 
}; 

static double cutoff_2d_1d(double w, void *p)
{
  double k0arg, k0;
  struct parameters_2d_1d *params = (struct parameters_2d_1d *)p;
  double  gx = (params->gx);
  double  gy = (params->gy);

  k0arg = fabs(w*gx);
  if(k0arg < 0.05){
    k0 = - (log(k0arg/2) + M_EULER_MASCHERONI);
  }
  else if(k0arg < 50.0 ){
    k0 = gsl_sf_bessel_K0(k0arg);
  }else{
    k0 = sqrt(2.0*M_PI/k0arg)*exp(-k0arg);
  }
  return 4.0*cos(w*gy)*k0;
}


/* --------------------- Interface to Fortran ---------------------- */
double FC_FUNC_(c_poisson_cutoff_2d_1d, C_POISSON_CUTOFF_2D_1D)
     (double *gy, double *gx, double *rc)
{
  double result, error, res, b;
  struct parameters_2d_1d params;
  const size_t wrk_size = 5000;
  const double epsilon_abs = 1e-3;
  const double epsilon_rel = 1e-3;

  double mu;
  gsl_integration_workspace * ws = gsl_integration_workspace_alloc (wrk_size);
  gsl_function F;

  mu = 0.1/(*gx);
  b  = fabs((*gx))/2.0;

  params.gx = *gx;
  params.gy = *gy;
  params.rc = *rc;

  F.function = &cutoff_2d_1d;
  F.params = &params;

  res = -4.0 * FC_FUNC(intcoslog, INTCOSLOG)(&mu, gy, &b);

  if( fabs(*gy) > 0.0) {
    res = res - (4.0 * M_EULER_MASCHERONI / (*gy)) * sin( (*gy)*mu );
  }else{
    res = res - 4.0 * M_EULER_MASCHERONI * mu;
  }
  
  gsl_integration_qag(&F, mu, (*rc), epsilon_abs, epsilon_rel, 
  			wrk_size, 3, ws, &result, &error);
  res = res + result;

  gsl_integration_workspace_free (ws);
  return res;
}




/************************************************************************/
/************************************************************************/
/* C_POISSON_CUTOFF_1D_0D */
/************************************************************************/
/************************************************************************/

struct parameters_1d_0d
{
  double g;
  double a;
};

static double cutoff_1d_0d(double w, void *p)
{
  struct parameters_1d_0d *params = (struct parameters_1d_0d *)p;
  double g = (params->g);
  double a = (params->a);
  return 2.0*(cos(w)/sqrt(w*w+a*a*g*g));
}



/* --------------------- Interface to Fortran ---------------------- */
double FC_FUNC_(c_poisson_cutoff_1d_0d, C_POISSON_CUTOFF_1D_0D)
     (double *g, double *a, double *rc)
{
  double result, error;
  struct parameters_1d_0d params;
  const size_t wrk_size = 5000;
  const double epsilon_abs = 1e-3;
  const double epsilon_rel = 1e-3;
  int status;
  gsl_integration_workspace * ws;
  gsl_function F;

  if((*g) <= 0.0){return 2.0*asinh((*rc)/(*a));};
  if((*g)*(*rc) > 100.0*M_PI){return 2.0*gsl_sf_bessel_K0((*a)*(*g));};

  gsl_set_error_handler_off();

  ws = gsl_integration_workspace_alloc (wrk_size);

  params.g = *g;
  params.a = *a;

  F.function = &cutoff_1d_0d;
  F.params = &params;

  status = gsl_integration_qag(&F, 0.0, (*rc)*(*g), epsilon_abs, epsilon_rel, 
    wrk_size, 3, ws, &result, &error);

  gsl_integration_workspace_free (ws);

  if(status){
    return 0.0;
  }else{
    return result;
  }

}
