#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <fftw3.h>
#include <omp.h>
#include <fstream>
 
#define DOUBLE

#ifdef DOUBLE
#define Complex  cufftDoubleComplex
#define Real double
#define Transform CUFFT_Z2Z
#define TransformExec cufftExecZ2Z
#else
#define Complex  cufftComplex
#define Real float
#define Transform CUFFT_C2C
#define TransformExec cufftExecC2C
#endif

#define TILE_DIM  4

using namespace std;

class dim3 {
 public:
 int x;
 int y;
 int z;

 dim3(int x, int y, int z) {
   this->x = x;
   this->y = y;
   this->z = z;
 }
};

void transpose_c(fftw_complex *idata, fftw_complex *odata,
        int width, int height, int gridDimx,
        int gridDimy, int blockDimx, int blockDimy)
{
  #pragma omp parallel 
  {
  #pragma omp for
  for (int by=0; by<gridDimy; by++)
    for (int bx=0; bx<gridDimx; bx++) {
      int yIndex = by * TILE_DIM;
      for (int ty=0; ty< blockDimy; ty++) {
          int xIndex = bx * TILE_DIM;
          int index_in  = xIndex + width * yIndex;
          int index_out = yIndex + height * xIndex;
          for (int tx=0; tx< blockDimx; tx++) {
                if (xIndex < width && yIndex < height) {
                        odata[index_out][0] = idata[index_in][0];
                        odata[index_out][1] = idata[index_in][1];
                }
                xIndex++;
                index_in++;
                index_out += height;
          }
          yIndex++;
      }
    }
  }
}

void transpose_spread_c(fftw_complex *idata, fftw_complex *odata,
        int width, int height, int bign_h, int gridDimx,
        int gridDimy, int blockDimx, int blockDimy)
{
  #pragma omp parallel
  {
  #pragma omp for
  for (int by=0; by<gridDimy; by++)
    for (int bx=0; bx<gridDimx; bx++) {
      int yIndex = by * TILE_DIM;
      for (int ty=0; ty< blockDimy; ty++) {
          int xIndex = bx * TILE_DIM;
          int index_in = xIndex + width * yIndex;
          int index_out = yIndex + height * xIndex;
          for (int tx=0; tx< blockDimx; tx++) {
          	int div = (index_out) / bign_h;
          	int mod = (index_out) % bign_h;
          	int index_out1 = div * (bign_h << 1) + mod+bign_h;
                if (xIndex < width && yIndex < height) {
                        odata[index_out1][0] = idata[index_in][0];
                        odata[index_out1][1] = idata[index_in][1];
                        odata[index_out1-bign_h][0] = 0.0;
                        odata[index_out1-bign_h][1] = 0.0;
                }
                xIndex++;
                index_in++;
                index_out += height;
          }
          yIndex++;
      }
    }
  }
}

void transpose_spread_i_c(fftw_complex *idata, fftw_complex *odata,
        int width, int height, int bign_h, int gridDimx,
        int gridDimy, int blockDimx, int blockDimy)
{
  #pragma omp parallel 
  {
  #pragma omp for
  for (int by=0; by<gridDimy; by++)
    for (int bx=0; bx<gridDimx; bx++) {
      int yIndex = by * TILE_DIM;
      for (int ty=0; ty< blockDimy; ty++) {
          int xIndex = bx * TILE_DIM;
          int index_in  = xIndex + width * yIndex;
          int index_out = yIndex + height * xIndex;
          for (int tx=0; tx< blockDimx; tx++) {
          	int div = (index_in) / bign_h;
          	int mod = (index_in) % bign_h;
          	int index_in1 = div * (bign_h << 1) + mod;
                if (xIndex < width && yIndex < height) {
                        odata[index_out][0] = idata[index_in1][0];
                        odata[index_out][1] = idata[index_in1][1];
                }
                xIndex++;
                index_out += height;
                index_in++;
          }
          yIndex++;
      }
    }
  }
}

void spread_c(fftw_complex *src, int spitch, fftw_complex *dst, int dpitch, int gridDimx, 
	int gridDimy, int blockDimx)
{

 fftw_complex zero;
 zero[0]=0.0;
 zero[1]=0.0;

 #pragma omp parallel 
 {
 #pragma omp for
 for (int by=0; by<gridDimy; by++) {
    int bid = by * gridDimx -1;
    for (int bx=0; bx<gridDimx; bx++) {
       bid++;
       for (int tid=0; tid< blockDimx; tid+=2) {
       	  fftw_complex res;
          res[0] = (tid >= spitch) ? src[(bid * spitch + tid-spitch)/2][0] : zero[0];
          res[1] = (tid >= spitch) ? src[(bid * spitch + tid-spitch)/2][1] : zero[1];
       	  if( tid < dpitch) { 
       	     dst[(bid * dpitch + tid)/2][0] = res[0];
       	     dst[(bid * dpitch + tid)/2][1] = res[1];
	  }
       }
    }
 }
 }
}

void spread_i_c(fftw_complex *src, int spitch, fftw_complex *dst, int dpitch, int gridDimx,
        int gridDimy, int blockDimx)
{
 #pragma omp parallel 
 {
 #pragma omp for
 for (int by=0; by<gridDimy; by++)
    for (int bx=0; bx<gridDimx; bx++) {
       int bid = by * gridDimx + bx;
       for (int tid=0; tid< blockDimx; tid+=2) {
   	  fftw_complex res;
          res[0] = src[(bid * dpitch + tid)/2][0];
          res[1] = src[(bid * dpitch + tid)/2][1];
   	  if( tid < dpitch) {
		dst[(bid * spitch + tid)/2][0] = res[0];
		dst[(bid * spitch + tid)/2][1] = res[1];
          }
       }
    }
 }
}

void multiply_kernel_c(int nx, int ny, int nz, fftw_complex *d_data, Real *d_kernel, Real scal,
			int gridDimx, int gridDimy, int blockDimx) {

 int td = blockDimx;

 #pragma omp parallel 
 {
 #pragma omp for
 for (int by=0; by<gridDimy; by++)
    for (int bx=0; bx<gridDimx; bx++) {

        for (int tj=0; tj< td; tj++) {

 		int blockData = (nx*ny*nz)/(gridDimx*gridDimy);

 		int jj = (by*gridDimx + bx)*blockData;

 		for (int k=0; k<blockData/td; k++) {
     			d_data[jj + tj+ k*td][0] *= d_kernel[jj + tj+ k*td]*scal;
     			d_data[jj + tj+ k*td][1] *= d_kernel[jj + tj+ k*td]*scal;
 		}
        }
    }
 }
}

extern "C" void fftw_3d_psolver_general_plan_(int *N,
                 fftw_plan *plan, int *geo, fftw_complex *in, fftw_complex *out) {

 int NX = N[0];
 int NY = N[1];
 int NZ = N[2];

 //int geo1 = geo[0];
 int geo2 = geo[1];
 int geo3 = geo[2];

 int n1d[3]= {1, 1, 1};

 int ysize = NY/2 + geo2 * NY/2;
 int zsize = NZ/2 + geo3 * NZ/2;

 fftw_plan_with_nthreads(omp_get_max_threads());

 n1d[0] = NX;
 plan[0] = fftw_plan_many_dft_r2c(1, n1d,ysize*zsize, (Real*)in,
              NULL, 1, NX, out,
              NULL, 1, NX/2+1,FFTW_MEASURE);
 plan[1] = fftw_plan_many_dft_r2c(1, n1d,ysize*zsize, (Real*)out,
              NULL, 1, NX, in,
              NULL, 1, NX/2+1,FFTW_MEASURE);

 plan[2] = fftw_plan_many_dft_c2r(1, n1d,ysize*zsize, in,
              NULL, 1, NX/2+1, (Real*)out,
              NULL, 1, NX, FFTW_MEASURE);
 plan[3] = fftw_plan_many_dft_c2r(1, n1d,ysize*zsize, out,
              NULL, 1, NX/2+1, (Real*)in,
              NULL, 1, NX, FFTW_MEASURE);

 n1d[0] = NY;
 plan[4] = fftw_plan_many_dft(1, n1d,(NX/2+1)*zsize, in,
              NULL, 1, NY, out,
              NULL, 1, NY, FFTW_FORWARD, FFTW_MEASURE);
 plan[5] = fftw_plan_many_dft(1, n1d,(NX/2+1)*zsize, out,
              NULL, 1, NY, in,
              NULL, 1, NY, FFTW_FORWARD, FFTW_MEASURE);

 plan[6] = fftw_plan_many_dft(1, n1d,(NX/2+1)*zsize, in,
              NULL, 1, NY, out,
              NULL, 1, NY, FFTW_BACKWARD, FFTW_MEASURE);
 plan[7] = fftw_plan_many_dft(1, n1d,(NX/2+1)*zsize, out,
              NULL, 1, NY, in,
              NULL, 1, NY, FFTW_BACKWARD, FFTW_MEASURE);

 n1d[0] = NZ;
 plan[8] = fftw_plan_many_dft(1, n1d,(NX/2+1)*NY, in,
              NULL, 1, NZ, out,
              NULL, 1, NZ, FFTW_FORWARD, FFTW_MEASURE);
 plan[9] = fftw_plan_many_dft(1, n1d,(NX/2+1)*NY, out,
              NULL, 1, NZ, in,
              NULL, 1, NZ, FFTW_FORWARD, FFTW_MEASURE);

 plan[10] = fftw_plan_many_dft(1, n1d,(NX/2+1)*NY, in,
              NULL, 1, NZ, out,
              NULL, 1, NZ, FFTW_BACKWARD, FFTW_MEASURE);
 plan[11] = fftw_plan_many_dft(1, n1d,(NX/2+1)*NY, out,
              NULL, 1, NZ, in,
              NULL, 1, NZ, FFTW_BACKWARD, FFTW_MEASURE);
}

extern "C" void fftw_3d_psolver_general_(int *N, fftw_plan *plan, int *geo,
          fftw_complex *d_data,  fftw_complex *d_data2, Real *d_kernel,
          Real *scal_p) {

 fftw_init_threads();


 int NX = N[0];
 int NY = N[1];
 int NZ = N[2];

 Real scal = *scal_p;

 int geo1 = geo[0];
 int geo2 = geo[1];
 int geo3 = geo[2];

 int ysize=NY/2+geo2*NY/2;
 int zsize=NZ/2+geo3*NZ/2;

 // transpose kernel parameters
 dim3 grid((NX/2+1+TILE_DIM-1)/TILE_DIM,(ysize*zsize+TILE_DIM-1)/TILE_DIM,1);
 dim3 threads(TILE_DIM,TILE_DIM,1);

 // spread kernel parameters
 dim3 nblocks(zsize,ysize,1);

 // multiply kernel paramters
 int nThreads = NX/2+1;
 dim3 nBlocks(NZ,NY,1);

 fftw_complex *dst = d_data;
 fftw_complex *src = d_data2;

 // X transform 
   if (geo1==0) {
     src = d_data;
     dst = d_data2;
     spread_c(src, NX/2, dst, NX, nblocks.x, nblocks.y, NX);
     fftw_execute(plan[1]);
   }
   else {
     fftw_execute(plan[0]);
   }

   if (geo2==0) {
     transpose_spread_c(src, dst,NX/2+1,ysize*zsize,NY/2, grid.x, grid.y, 
		threads.x, threads.y);
   } else {
     transpose_c(src, dst,NX/2+1,ysize*zsize, grid.x, grid.y,                      
                threads.x, threads.y);
   }


   if (geo1==0) 
   	fftw_execute(plan[5]);
   else
        fftw_execute(plan[4]);

   grid.x = (NY+TILE_DIM-1)/TILE_DIM;
   grid.y = ((NX/2+1)*zsize+TILE_DIM-1)/TILE_DIM;

   if (geo3==0) {
     transpose_spread_c(src, dst,NY,(NX/2+1)*NZ/2,NZ/2, grid.x, grid.y,
                threads.x, threads.y);
   } else {
     transpose_c(src, dst,NY,(NX/2+1)*NZ, grid.x, grid.y,
                threads.x, threads.y);
   }

   if (geo1==0)
        fftw_execute(plan[9]);
   else
        fftw_execute(plan[8]);

   // multiply with kernel

   multiply_kernel_c(NX/2+1,NY,NZ,src,d_kernel,scal, nBlocks.x, nBlocks.y, nThreads);

   // inverse transform

  // Z transform, on entire cube 

   if (geo1==0)
        fftw_execute(plan[10]);
   else
        fftw_execute(plan[11]);
 
   
   grid.x = (zsize*(NX/2+1)+TILE_DIM-1)/TILE_DIM;
   grid.y = (NY+TILE_DIM-1)/TILE_DIM;

   if (geo3==0) {
     transpose_spread_i_c(dst,src,NZ/2*(NX/2+1),NY,NZ/2, grid.x, grid.y,
                threads.x, threads.y);
   } else {
     transpose_c(dst, src,NZ*(NX/2+1),NY, grid.x, grid.y,
                threads.x, threads.y);
   }

   // Y transform

   if (geo1==0)
        fftw_execute(plan[6]);
   else
        fftw_execute(plan[7]);

   grid.x = (ysize*zsize+TILE_DIM-1)/TILE_DIM;
   grid.y = (NX/2+1+TILE_DIM-1)/TILE_DIM;

   if (geo2==0) {
      transpose_spread_i_c(dst, src,ysize*zsize,NX/2+1, NY/2, grid.x, grid.y,
                threads.x, threads.y);
   } else
      transpose_c(dst, src,ysize*zsize,NX/2+1, grid.x, grid.y,
                threads.x, threads.y);

   // X transform

   nblocks.x=zsize;
   nblocks.y=ysize;
   if (geo1==0) {
        fftw_execute(plan[2]);
      spread_i_c(dst, NX/2, src, NX, nblocks.x, nblocks.y, NX/2);
   }
   else
        fftw_execute(plan[3]);

}
