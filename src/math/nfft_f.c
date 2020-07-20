/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira

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

*/

#include <config.h>


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#ifdef HAVE_NFFT
#include "nfft3util.h"
#include "nfft3.h"
#else
typedef int nfft_plan;
#endif

// NFFT FUNCTIONS
void FC_FUNC(oct_nfft_init_1d,OCT_NFFT_INIT_1D)
   (nfft_plan *plan, int *N1, int *M){
#ifdef HAVE_NFFT
    nfft_init_1d(plan, *N1, *M);
#endif
}

void FC_FUNC(oct_nfft_init_2d,OCT_NFFT_INIT_2D)
   (nfft_plan *plan, int *N1, int *N2, int *M){
#ifdef HAVE_NFFT
    nfft_init_2d(plan, *N1, *N2, *M);
#endif
}

void FC_FUNC(oct_nfft_init_3d,OCT_NFFT_INIT_3D)
   (nfft_plan *plan, int *N1, int *N2, int *N3, int *M){
#ifdef HAVE_NFFT
    nfft_init_3d(plan, *N1, *N2, *N3, *M);
#endif
}

void FC_FUNC(oct_nfft_init_guru,OCT_NFFT_INIT_GURU) 
  (nfft_plan *ths, int *d, int *N, int *M, int *n, int *m, unsigned *nfft_flags, unsigned *fftw_flags){
#ifdef HAVE_NFFT
  nfft_init_guru (ths, *d, N, *M, n, *m, *nfft_flags, *fftw_flags);
#endif
}

void 	FC_FUNC(oct_nfft_check, OCT_NFFT_CHECK)(nfft_plan *ths){
#ifdef HAVE_NFFT
  nfft_check (ths);
#endif
}

void 	FC_FUNC(oct_nfft_finalize,OCT_NFFT_FINALIZE)(nfft_plan *ths){
#ifdef HAVE_NFFT
  nfft_finalize (ths);
#endif
}

void 	FC_FUNC(oct_nfft_trafo,OCT_NFFT_TRAFO)(nfft_plan *ths){
#ifdef HAVE_NFFT
  nfft_trafo (ths);
#endif
}  
  
void 	FC_FUNC(oct_nfft_adjoint,OCT_NFFT_ADJOINT)(nfft_plan *ths){
#ifdef HAVE_NFFT
  nfft_adjoint (ths);
#endif
}
  
void FC_FUNC(oct_nfft_precompute_one_psi_1d, OCT_NFFT_PRECOMPUTE_ONE_PSI_1D)(nfft_plan *plan, int *m, double *X1){

#ifdef HAVE_NFFT
   int ii;
   int M = *m;
 
   for (ii=0;ii< M;ii++) plan->x[ii] =  X1[ii];
  
   if(plan->nfft_flags & PRE_ONE_PSI) nfft_precompute_one_psi(plan);
#endif 
}

void FC_FUNC(oct_nfft_precompute_one_psi_2d, OCT_NFFT_PRECOMPUTE_ONE_PSI_2D)(nfft_plan *plan, int *M, double* X1, double* X2){

#ifdef HAVE_NFFT
	int ii;
   int jj;

   for (ii=0; ii< M[0]; ii++){
     for (jj=0; jj< M[1]; jj++){
       plan->x[2*(M[1] * ii + jj) + 0] =  X1[ii];
       plan->x[2*(M[1] * ii + jj) + 1] =  X2[jj];
     }
    }

 
   if(plan->nfft_flags & PRE_ONE_PSI)nfft_precompute_one_psi(plan);
#endif
}

void FC_FUNC(oct_nfft_precompute_one_psi_3d, OCT_NFFT_PRECOMPUTE_ONE_PSI_3D)
     (nfft_plan *plan, int *M, double* X1, double* X2, double* X3){
#ifdef HAVE_NFFT
   int ii,jj,kk;

   for (ii=0;ii< M[0];ii++){
     for (jj=0;jj< M[1];jj++){
       for (kk=0;kk< M[2];kk++){
         plan->x[3*(M[1]*M[2]*ii + M[2]*jj + kk) + 0] =  X1[ii];
         plan->x[3*(M[1]*M[2]*ii + M[2]*jj + kk) + 1] =  X2[jj];
         plan->x[3*(M[1]*M[2]*ii + M[2]*jj + kk) + 2] =  X3[kk];
       }
     }
   }
  
   if(plan->nfft_flags & PRE_ONE_PSI) nfft_precompute_one_psi(plan);
#endif
}

// Type dependent functions 


// ********** COMPLEX ************
void FC_FUNC(zoct_set_f, ZOCT_SET_F)
    (nfft_plan *plan, int *M, int *DIM, double complex *VAL, int *IX, int *IY, int *IZ){

#ifdef HAVE_NFFT
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;
  double complex val = *VAL;
  
   switch (dim){
     case 1:
       plan->f[ix-1] = val;
     break;
     case 2:
     plan->f[(ix-1)*M[1] + (iy-1)] = val;
     break;
     case 3:
       plan->f[(ix-1)*M[1]*M[2] + (iy-1)*M[2] + (iz-1)] = val;
     break;
     }
#endif
}

void FC_FUNC(zoct_get_f, ZOCT_GET_F)
    (nfft_plan *plan, int *M, int *DIM, double complex *val, int *IX, int *IY, int *IZ){

#ifdef HAVE_NFFT
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;
   

   switch (dim){
     case 1:
       *val = plan->f[ix-1];
     break;
     case 2:
       *val = plan->f[(ix-1)*M[1] + (iy-1)];
     break;
     case 3:
       *val = plan->f[(ix-1)*M[1]*M[2] + (iy-1)*M[2] + (iz-1)];
     break;
   }
#endif
}


void FC_FUNC(zoct_set_f_hat, ZOCT_SET_F_HAT)
    (nfft_plan *plan, int *DIM, double complex *VAL, int *IX, int *IY, int *IZ){

#ifdef HAVE_NFFT
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;
  double complex val = *VAL;

   switch (dim){
     case 1:
       plan->f_hat[ix-1] = val;
     break;
     case 2:
       plan->f_hat[(ix-1)*plan->N[1] + (iy-1)] = val;
     break;
     case 3:
       plan->f_hat[(ix-1)*plan->N[1]*plan->N[2] + (iy-1)*plan->N[2] + (iz-1)] = val;
     break;
   }
#endif
}

void FC_FUNC(zoct_get_f_hat, ZOCT_GET_F_HAT)
    (nfft_plan *plan, int *DIM, double complex *val, int *IX, int *IY, int *IZ){

#ifdef HAVE_NFFT
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;


   switch (dim){
     case 1:
       *val = plan->f_hat[ix-1];
     break;
     case 2:
       *val = plan->f_hat[(ix-1)*plan->N[1] + (iy-1)];
     break;
     case 3:
       *val = plan->f_hat[(ix-1)*plan->N[1]*plan->N[2] + (iy-1)*plan->N[2] + (iz-1)];
     break;
   }
#endif
}

// ********** DOUBLE ************
void FC_FUNC(doct_set_f, DOCT_SET_F)
    (nfft_plan *plan, int *M, int *DIM, double *VAL, int *IX, int *IY, int *IZ){

#ifdef HAVE_NFFT
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;
  double val = *VAL;
  

   switch (dim){
     case 1:
       plan->f[ix-1] = val;
     break;
     case 2:
       plan->f[(ix-1)*M[1] + (iy-1)] = val;
     break;
     case 3:
       plan->f[(ix-1)*M[1]*M[2] + (iy-1)*M[2] + (iz-1)] = val;
     break;
   }
#endif
}

void FC_FUNC(doct_get_f, DOCT_GET_F)
    (nfft_plan *plan, int *M, int *DIM, double *val, int *IX, int *IY, int *IZ){

#ifdef HAVE_NFFT
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;
   

   switch (dim){
     case 1:
       *val = plan->f[ix-1];
     break;
     case 2:
       *val = plan->f[(ix-1)*M[1] + (iy-1)];
     break;
     case 3:
       *val = plan->f[(ix-1)*M[1]*M[2] + (iy-1)*M[2] + (iz-1)];
     break;
   }
#endif
}


void FC_FUNC(doct_set_f_hat, DOCT_SET_F_HAT)
    (nfft_plan *plan, int *DIM, double *VAL, int *IX, int *IY, int *IZ){

#ifdef HAVE_NFFT
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;
  double  val = *VAL;

   switch (dim){
     case 1:
       plan->f_hat[ix-1] = val;
     break;
     case 2:
       plan->f_hat[(ix-1)*plan->N[1] + (iy-1)] = val;
     break;
     case 3:
       plan->f_hat[(ix-1)*plan->N[1]*plan->N[2] + (iy-1)*plan->N[2] + (iz-1)] = val;
     break;
   }
#endif
}

void FC_FUNC(doct_get_f_hat, DOCT_GET_F_HAT)
    (nfft_plan *plan, int *DIM, double *val, int *IX, int *IY, int *IZ){
	
#ifdef HAVE_NFFT
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;


   switch (dim){
     case 1:
       *val = plan->f_hat[ix-1];
     break;
     case 2:
       *val = plan->f_hat[(ix-1)*plan->N[1] + (iy-1)];
     break;
     case 3:
       *val = plan->f_hat[(ix-1)*plan->N[1]*plan->N[2] + (iy-1)*plan->N[2] + (iz-1)];
     break;
   }
#endif
}


