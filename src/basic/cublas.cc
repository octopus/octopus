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
typedef int cublasHandle_t;
typedef int CUdeviceptr;
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
