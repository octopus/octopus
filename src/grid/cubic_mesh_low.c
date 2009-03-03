#include <config.h>
#include <stdlib.h>
#include <assert.h>

typedef struct {
  const int * lxyz;
  double * ff;
  int nx;
  int ny;
  int nz;
  int nvalues;
  int npoints;
  int np_part;
  int np;
  int ox;
  int oy;
  int oz;
} cubic_mesh_t;

#define index_c(nx, ny, nz, ix, iy, iz) (ix) + (nx)*((iy) + (ny)*(iz))

void FC_FUNC_(cubic_mesh_init_c, CUBIC_MESH_INIT_C)(cubic_mesh_t ** this,
						    const int * nx, const int * ny, const int * nz,
						    const int * ox, const int * oy, const int * oz,
						    const int * np, const int * np_part, const int * lxyz){

  int ii;

  this[0] = (cubic_mesh_t *) malloc(sizeof(cubic_mesh_t));

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
