/*
 Copyright (C) 2008 X. Andrade

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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: fft_cuda.c 2146 2006-05-23 17:36:00Z xavier $
*/

#include <config.h>

#if HAVE_FFT == cuda

#include <cufft.h>

typedef struct {
  int nx, ny, nz;
  cufftHandle * r2c;
  cufftHandle * c2r;
  cufftHandle * c2c;
} fft_t;

void FC_FUNC_(fft_all_init, FFT_ALL_INIT)
(){
}
    
void FC_FUNC_(fft_all_end, FFT_ALL_END)
(){
}
    
void FC_FUNC_(fft_init, FFT_INIT)
     (int n[3], int * is_real, fft_t * fft){
  
  cufftResult result;

  fft->nx = n[1];
  fft->ny = n[2];
  fft->nz = n[3];

  fft->r2c = 0;
  fft->c2r = 0;
  fft->c2c = 0;

  if(*is_real){
    result = cufftPlan3d(fft->r2c, fft->nx, fft->ny, fft->nz, CUFFT_R2C);
    result = cufftPlan3d(fft->c2r, fft->nx, fft->ny, fft->nz, CUFFT_C2R);
  } else {
    result = cufftPlan3d(fft->c2c, fft->nx, fft->ny, fft->nz, CUFFT_C2C);
  }

}

void FC_FUNC_(fft_copy, FFT_COPY)
     (fft_t * fft_i, fft_t * fft_o){
}

void FC_FUNC_(fft_end, FFT_END)
     (fft_t * fft){
  if(fft->r2c) cufftDestroy(fft->r2c);
  if(fft->c2r) cufftDestroy(fft->c2r);
  if(fft->c2c) cufftDestroy(fft->c2c);
}

void FC_FUNC_(fft_getdim_real, FFT_GETDIM_REAL)
(const fft_t * fft, int d[3]){
}

void FC_FUNC_(fft_getdim_complex, FFT_GETDIM_COMPLEX)
(const fft_t * fft, int d[3]){
}

void FC_FUNC_(dfft_forward, DFFT_FORWARD)
(const fft_t * fft, const double * r, double *c){
}

void FC_FUNC_(dfft_backward, DFFT_BACKWARD)
(const fft_t * fft, const double * c, double *r){
}
    
void FC_FUNC_(zfft_forward, ZFFT_FORWARD)
(const fft_t * fft, const double * in, double * out){
}
    
void FC_FUNC_(zfft_backward, ZFFT_BACKWARD)
(const fft_t * fft, const double * in, double * out){
}

int FC_FUNC_(pad_feq, PAD_FEQ)
(const * i, const * n, const * mode){
}

#endif
