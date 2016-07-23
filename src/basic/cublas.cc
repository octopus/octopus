/*
 Copyright (C) 2016 X. Andrade

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

 $Id$
*/

#include <config.h>

#include <cstdio>
#include <iostream>

#include <fortran_types.h>

#ifdef HAVE_CUDA
#include <cuda.h>
#include <cublas_v2.h>
#else

#include <complex>

typedef int cublasHandle_t;
typedef int CUdeviceptr;
typedef std::complex<double> cuDoubleComplex;
#endif

#define CUBLAS_SAFE_CALL(x)                                       \
  do {                                                            \
    cublasStatus_t safe_call_result = x;                          \
    if(safe_call_result != CUBLAS_STATUS_SUCCESS) {               \
      std::cerr << "\nerror: " #x " failed with error\n";	  \
      exit(1);                                                    \
    }                                                             \
  } while(0)

extern "C" void FC_FUNC_(cublas_init, CUBLAS_INIT)(cublasHandle_t ** handle){
#ifdef HAVE_CUDA
  *handle = new cublasHandle_t;
  CUBLAS_SAFE_CALL(cublasCreate(*handle));
  CUBLAS_SAFE_CALL(cublasSetPointerMode(**handle, CUBLAS_POINTER_MODE_DEVICE));
#endif
}

extern "C" void FC_FUNC_(cublas_end, CUBLAS_END)(cublasHandle_t ** handle){
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasDestroy(**handle));
  delete *handle;
#endif
}

extern "C" void FC_FUNC_(cuda_blas_ddot, CUDA_BLAS_DDOT)
  (cublasHandle_t ** handle, fint8 *n, CUdeviceptr ** x, fint8 * offx, fint8 * incx,
   CUdeviceptr ** y, fint8 * offy, fint8 * incy, CUdeviceptr ** result, fint8 * off_result){
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasDdot(**handle, *n, ((double *) **x) + *offx, *incx,
			      ((double *) **y) + *offy, *incy, ((double *)**result) + *off_result));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_zdotc, CUDA_BLAS_ZDOTC)
  (cublasHandle_t ** handle, fint8 *n, CUdeviceptr ** x, fint8 * offx, fint8 * incx,
   CUdeviceptr ** y, fint8 * offy, fint8 * incy, CUdeviceptr ** result, fint8 * off_result){
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasZdotc(**handle, *n, ((cuDoubleComplex *) **x) + *offx, *incx,
			      ((cuDoubleComplex *) **y) + *offy, *incy, ((cuDoubleComplex *)**result) + *off_result));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_dgemm, CUDA_BLAS_DGEMM)
  (cublasHandle_t ** handle, fint * transa, fint * transb, fint * m, fint * n, fint * k,
   const double *alpha, CUdeviceptr ** A, fint * lda, CUdeviceptr ** B, fint * ldb, 
   const double *beta, CUdeviceptr ** C, fint * ldc){
  
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasDgemm(**handle, (cublasOperation_t) *transa, (cublasOperation_t) *transb, *m, *n, *k,
			       alpha, (double *) **A, *lda, (double *) **B, *ldb,
			       beta, (double *) **C, *ldc));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_zgemm, CUDA_BLAS_ZGEMM)
  (cublasHandle_t ** handle, fint * transa, fint * transb, fint * m, fint * n, fint * k,
   const cuDoubleComplex * alpha, CUdeviceptr ** A, fint * lda, CUdeviceptr ** B, fint * ldb, 
   const cuDoubleComplex * beta, CUdeviceptr ** C, fint * ldc){
  
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasZgemm(**handle, (cublasOperation_t) *transa, (cublasOperation_t) *transb, *m, *n, *k,
			       alpha, (cuDoubleComplex *) **A, *lda, (cuDoubleComplex *) **B, *ldb,
			       beta, (cuDoubleComplex *) **C, *ldc));
#endif
}