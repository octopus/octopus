#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <assert.h>

struct parameters { double kx; double kr; double x0; double r0; int index; }; 

/* the index is the tipe of integral this is necessary cause dependent of the value of k the */
/* integral needs to be performed in several ways*/
#define PFC_F_0   0
#define PFC_F_KX  1
#define PFC_F_KR  2

static const double tol = 1e-7;
static double zero = 1e-100;

static double f(double w, void *p) {
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

static double f_kx_zero(double w, void *p) {
  struct parameters *params = (struct parameters *)p;
  double  kr = (params->kr);
  double  x0 = (params->x0);
  double result;
  
  if (w == 0.0) w = zero;
  
  result = 2.0*M_PI*w*gsl_sf_bessel_J0(kr*w)*
    (log(x0 + sqrt(w*w + x0*x0)) - log(-x0 + sqrt(w*w + x0*x0)));
  return result;
}


static double f_kr_zero(double w, void *p) {
  struct parameters *params = (struct parameters *)p;
  double  kx = (params->kx);
  double  r0 = (params->r0);
  double result;

  result = 4.0*M_PI*cos(kx*w)*(sqrt(r0*r0 + w*w) - fabs(w));
  return result;
}


static double int_aux(int index, double kx, double kr, double x0, double r0, int which_int) {

  assert(which_int == PFC_F_0 || which_int == PFC_F_KR || which_int == PFC_F_KX);

  /*variables*/
  double result, error;
  double xinf = 100.0/r0;
  struct parameters params = {kx, kr, x0, r0, index};

  /*pass the given parameters to program*/
  const size_t wrk_size = 5000;
  /* I believe that 1e-6 is over-killing, a substantial gain in time
     is obtained by setting this threshold to 1e-3 */
  /* const double epsilon_abs = 1e-6; */
  /* const double epsilon_rel = 1e-6; */
  const double epsilon_abs = 1e-3;
  const double epsilon_rel = 1e-3;

  /*creating the workspace*/
  gsl_integration_workspace * ws
    = gsl_integration_workspace_alloc (wrk_size);

  /*create the function*/
  gsl_function F;
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
