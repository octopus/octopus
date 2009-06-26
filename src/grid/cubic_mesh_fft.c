#include <config.h>
#include <stdlib.h>
#include <assert.h>

#include "cubic_mesh.h"

#include <fftw3.h>

void FC_FUNC_(cubic_mesh_init_fourier, CUBIC_MESH_INIT_FOURIER)(cubic_mesh_t ** this, const cubic_mesh_t ** real){

  int ii;

  this[0] = (cubic_mesh_t *) malloc(sizeof(cubic_mesh_t));

  this[0]->type = TYPE_CMPLX;

  this[0]->nx = real[0]->nx;
  this[0]->ny = real[0]->ny;
  this[0]->nz = real[0]->nz/2 + 1;

  this[0]->npoints = this[0]->nx*this[0]->ny*this[0]->nz;
  this[0]->nvalues = 2*this[0]->npoints;

  this[0]->ox = real[0]->ox;
  this[0]->oy = real[0]->oy;
  this[0]->oz = real[0]->oz;

  this[0]->np = real[0]->np;
  this[0]->np_part = real[0]->np_part;
  this[0]->lxyz = real[0]->lxyz;

  this[0]->ff = (double *) malloc(sizeof(double)*this[0]->nvalues);

  for(ii = 0; ii < this[0]->nvalues; ii++) this[0]->ff[ii] = 0.0;
  printf("Init fourier %d %d %d\n", this[0]->nx, this[0]->ny, this[0]->nz);

}

void FC_FUNC_(cubic_mesh_fourier_transform, CUBIC_MESH_FOURIER_TRANSFORM)(const cubic_mesh_t ** real, cubic_mesh_t ** fourier){

  fftw_plan plan;

  plan = fftw_plan_dft_r2c_3d(real[0]->nx, real[0]->ny, real[0]->nz,
			      real[0]->ff, (fftw_complex *) fourier[0]->ff, FFTW_ESTIMATE);

  fftw_execute(plan);

  fftw_destroy_plan(plan);

}

void FC_FUNC_(cubic_mesh_fourier_trans_back, CUBIC_MESH_FOURIER_TRANSFORM_BACK)(cubic_mesh_t ** fourier, const cubic_mesh_t ** real){
  
  fftw_plan plan;
  int i;
  plan = fftw_plan_dft_c2r_3d(real[0]->nx, real[0]->ny, real[0]->nz,
			      (fftw_complex *) fourier[0]->ff, real[0]->ff, FFTW_ESTIMATE);

  fftw_execute(plan);
  
  for(i = 0; i < real[0] -> nvalues; i++) real[0]->ff[i] *= 1.0/(real[0]->nx*((double)real[0]->ny)*real[0]->nz);
  
  fftw_destroy_plan(plan);
  
}

