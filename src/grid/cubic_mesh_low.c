#include <config.h>
#include <stdlib.h>
#include <assert.h>
#include "cubic_mesh.h"

void FC_FUNC_(cubic_mesh_init_c, CUBIC_MESH_INIT_C)(cubic_mesh_t ** this,
						    const int * nx, const int * ny, const int * nz,
						    const int * ox, const int * oy, const int * oz,
						    const int * np, const int * np_part, const int * lxyz){

  int ii;

  this[0] = (cubic_mesh_t *) malloc(sizeof(cubic_mesh_t));

  this[0]->type = TYPE_FLOAT;

  this[0]->nx = nx[0];
  this[0]->ny = ny[0];
  this[0]->nz = nz[0];

  this[0]->npoints = this[0]->nx*this[0]->ny*this[0]->nz;
  this[0]->nvalues = this[0]->npoints;

  this[0]->ox = ox[0];
  this[0]->oy = oy[0];
  this[0]->oz = oz[0];

  this[0]->np = np[0];
  this[0]->np_part = np_part[0];
  this[0]->lxyz = lxyz;

  this[0]->ff = (double *) malloc(sizeof(double)*this[0]->nvalues);

  for(ii = 0; ii < this[0]->nvalues; ii++) this[0]->ff[ii] = 0.0;

}

void FC_FUNC_(dcubic_mesh_from_mesh, DCUBIC_MESH_FROM_MESH)(cubic_mesh_t ** this, const double * func){
  const int nx = this[0]->nx;
  const int ny = this[0]->ny;
  const int nz = this[0]->nz;

  int ip, ix, iy, iz, ii;


  for(ip = 0; ip < this[0]->np; ip++){
    ix = this[0]->lxyz[ip                     ] - this[0]->ox;
    iy = this[0]->lxyz[ip +   this[0]->np_part] - this[0]->oy;
    iz = this[0]->lxyz[ip + 2*this[0]->np_part] - this[0]->oz;

    //    printf("%d %d %d - %d %d %d\n", ix, iy, iz, nx, ny, nz);

    assert(ix >= 0);
    assert(iy >= 0);
    assert(iz >= 0);

    assert(ix < nx);
    assert(iy < ny);
    assert(iz < nz);

    ii = index_c(nx, ny, nz, ix, iy, iz);

    assert(ii >= 0);
    assert(ii < this[0]->npoints);

    this[0]->ff[ii] = func[ip];
  }
}

void FC_FUNC_(dcubic_mesh_to_mesh, DCUBIC_MESH_TO_MESH)(const cubic_mesh_t ** this, double * func){
  const int nx = this[0]->nx;
  const int ny = this[0]->ny;
  const int nz = this[0]->nz;

  int ix, iy, iz, ip, ii;

  for(ip = 0; ip < this[0]->np; ip++){
    ix = this[0]->lxyz[ip                     ] - this[0]->ox;
    iy = this[0]->lxyz[ip +   this[0]->np_part] - this[0]->oy;
    iz = this[0]->lxyz[ip + 2*this[0]->np_part] - this[0]->oz;

    ii = index_c(nx, ny, nz, ix, iy, iz);

    assert(ii < this[0]->npoints);
    assert(ii >= 0);

    func[ip] = this[0]->ff[ii];
  }
}

void FC_FUNC_(cubic_mesh_end, CUBIC_MESH_END)(cubic_mesh_t ** this){
  free(this[0]->ff);
  free(this[0]);
}

void FC_FUNC_(cubic_mesh_dimensions, CUBIC_MESH_DIMENSIONS)
     (const cubic_mesh_t ** this, int * nx, int * ny, int * nz){
  nx[0] = this[0]->nx;
  ny[0] = this[0]->ny;
  nz[0] = this[0]->nz;
}

void FC_FUNC_(dcubic_mesh_set_point, DCUBIC_MESH_SET_POINT)
     (const cubic_mesh_t ** this, const int * ix, const int * iy, const int * iz, const double * vv){
  int ii;

  assert(ix[0] > 0 && ix[0] <= this[0]->nx);
  assert(iy[0] > 0 && iy[0] <= this[0]->ny);
  assert(iz[0] > 0 && iz[0] <= this[0]->nz);
  assert(this[0]->type == TYPE_FLOAT);

  ii = index_c(this[0]->nx, this[0]->ny, this[0]->nz, ix[0] - 1, iy[0] - 1, iz[0] - 1);

  this[0]->ff[ii] = vv[0];
}

void FC_FUNC_(zcubic_mesh_set_point, ZCUBIC_MESH_SET_POINT)
     (const cubic_mesh_t ** this, const int * ix, const int * iy, const int * iz, const double * vv){
  int ii;

  assert(ix[0] > 0 && ix[0] <= this[0]->nx);
  assert(iy[0] > 0 && iy[0] <= this[0]->ny);
  assert(iz[0] > 0 && iz[0] <= this[0]->nz);
  assert(this[0]->type == TYPE_CMPLX);

  ii = index_c(this[0]->nx, this[0]->ny, this[0]->nz, ix[0] - 1, iy[0] - 1, iz[0] - 1);

  this[0]->ff[2*ii    ] = vv[0];
  this[0]->ff[2*ii + 1] = vv[1];
}

void FC_FUNC_(cubic_mesh_multiply, DCUBIC_MESH_MULTIPLY)
     (const cubic_mesh_t ** aa, cubic_mesh_t ** bb){
  
  int ia, ib;
  
  printf("aa %d %d %d\n", aa[0]->nx, aa[0]->ny, aa[0]->nz);
  printf("bb %d %d %d\n", bb[0]->nx, bb[0]->ny, bb[0]->nz);

  ib = 0;
  for(ia = 0; ia < aa[0]-> npoints; ia++) {
    //    printf("%d %d %d %d\n", ia, ib, aa[0]-> npoints, bb[0]-> npoints);
    assert(ia < aa[0] -> nvalues);
    assert(ib + 1 < bb[0] -> nvalues);

    bb[0]->ff[ib    ] *= aa[0]->ff[ia];
    bb[0]->ff[ib + 1] *= aa[0]->ff[ia];
    ib+=2;
    //    printf("kernel %d %f \n", ia, aa[0]->ff[ia]);
  }
}
