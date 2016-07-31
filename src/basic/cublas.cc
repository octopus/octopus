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

extern "C" void FC_FUNC_(cuda_blas_dnrm2, CUDA_BLAS_DNRM2)
  (cublasHandle_t ** handle, fint8 *n, CUdeviceptr ** x, fint8 * offx, fint8 * incx,
   CUdeviceptr ** result, fint8 * off_result){
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasDnrm2(**handle, *n, ((double *) **x) + *offx, *incx, ((double *)**result) + *off_result));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_znrm2, CUDA_BLAS_ZNRM2)
  (cublasHandle_t ** handle, fint8 *n, CUdeviceptr ** x, fint8 * offx, fint8 * incx,
   CUdeviceptr ** result, fint8 * off_result){
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasDznrm2(**handle, *n, ((cuDoubleComplex *) **x) + *offx, *incx, ((double *)**result) + *off_result));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_dgemm, CUDA_BLAS_DGEMM)
  (cublasHandle_t ** handle, fint * transa, fint * transb, fint8 * m, fint8 * n, fint8 * k,
   CUdeviceptr ** alpha, CUdeviceptr ** A, fint8 * lda, CUdeviceptr ** B, fint8 * ldb, 
   CUdeviceptr ** beta, CUdeviceptr ** C, fint8 * ldc){
  
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasDgemm(**handle, (cublasOperation_t) *transa, (cublasOperation_t) *transb, *m, *n, *k,
			       (double *) ** alpha, (double *) **A, *lda, (double *) **B, *ldb,
			       (double *) ** beta, (double *) **C, *ldc));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_zgemm, CUDA_BLAS_ZGEMM)
  (cublasHandle_t ** handle, fint * transa, fint * transb, fint8 * m, fint8 * n, fint8 * k,
   CUdeviceptr ** alpha, CUdeviceptr ** A, fint8 * lda, CUdeviceptr ** B, fint8 * ldb, 
   CUdeviceptr ** beta, CUdeviceptr ** C, fint8 * ldc){
  
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasZgemm(**handle, (cublasOperation_t) *transa, (cublasOperation_t) *transb, *m, *n, *k,
			       (cuDoubleComplex *) **alpha, (cuDoubleComplex *) **A, *lda, (cuDoubleComplex *) **B, *ldb,
			       (cuDoubleComplex *) **beta, (cuDoubleComplex *) **C, *ldc));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_dsyrk, CUDA_BLAS_DSYRK)
  (cublasHandle_t ** handle, fint * uplo, fint * trans, fint8 * n, fint8 * k,
   CUdeviceptr **alpha, CUdeviceptr ** A, fint8 * lda,
   CUdeviceptr **beta, CUdeviceptr ** C, fint8 * ldc){
  
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasDsyrk(**handle, (cublasFillMode_t) *uplo, (cublasOperation_t) *trans, *n, *k,
			       (double *) **alpha, (double *) **A, *lda,
			       (double *) **beta, (double *) **C, *ldc));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_zherk, CUDA_BLAS_ZHERK)
  (cublasHandle_t ** handle, fint * uplo, fint * trans, fint8 * n, fint8 * k,
   CUdeviceptr **alpha, CUdeviceptr ** A, fint8 * lda,
   CUdeviceptr **beta, CUdeviceptr ** C, fint8 * ldc){
  
#ifdef HAVE_CUDA
  CUBLAS_SAFE_CALL(cublasZherk(**handle, (cublasFillMode_t) *uplo, (cublasOperation_t) *trans, *n, *k,
			       (double *) **alpha, (cuDoubleComplex *) **A, *lda,
			       (double *) **beta, (cuDoubleComplex *) **C, *ldc));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_dtrsm, CUDA_BLAS_DTRSM)
  (cublasHandle_t ** handle, fint * side, fint * uplo, fint * trans, fint * diag,
   fint8 * m, fint8 * n,
   CUdeviceptr **alpha, CUdeviceptr ** A, fint8 * lda,
   CUdeviceptr ** B, fint8 * ldb){
  
#ifdef HAVE_CUDA  
  CUBLAS_SAFE_CALL(cublasDtrsm(**handle, (cublasSideMode_t) *side, (cublasFillMode_t) *uplo,
			       (cublasOperation_t) *trans, (cublasDiagType_t) *diag,
			       *m, *n, (double *) **alpha, (double *) **A, *lda, (double *) **B, *ldb));
#endif
}

extern "C" void FC_FUNC_(cuda_blas_ztrsm, CUDA_BLAS_ZTRSM)
  (cublasHandle_t ** handle, fint * side, fint * uplo, fint * trans, fint * diag,
   fint8 * m, fint8 * n,
   CUdeviceptr **alpha, CUdeviceptr ** A, fint8 * lda,
   CUdeviceptr ** B, fint8 * ldb){
  
#ifdef HAVE_CUDA  
  CUBLAS_SAFE_CALL(cublasZtrsm(**handle, (cublasSideMode_t) *side, (cublasFillMode_t) *uplo,
			       (cublasOperation_t) *trans, (cublasDiagType_t) *diag,
			       *m, *n, (cuDoubleComplex *) **alpha, (cuDoubleComplex *) **A, *lda, (cuDoubleComplex *) **B, *ldb));
#endif
}
