#ifndef CUBIC_MESH
#define CUBIC_MESH

#include <config.h>

#define TYPE_CMPLX 1
#define TYPE_FLOAT 2

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
  int type;
} cubic_mesh_t;

inline static int index_c(const int nx, const int ny, const int nz, const int ix, const int iy, const int iz){
  return ix + nx*(iy + ny*iz);
}

#endif
