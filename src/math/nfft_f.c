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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: oct_gsl_f.c 7086 2010-11-15 05:45:08Z acastro $
*/

#include <config.h>

#if defined(HAVE_NFFT) 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"

#include "string_f.h"


// NFFT FUNCTIONS
void FC_FUNC(oct_nfft_init_1d,OCT_NFFT_INIT_1D)
   (nfft_plan *plan, int *N1, int *M)
{
    nfft_init_1d(plan, *N1, *M);  
}

void FC_FUNC(oct_nfft_init_2d,OCT_NFFT_INIT_2D)
   (nfft_plan *plan, int *N1, int *N2, int *M)
{
    nfft_init_2d(plan, *N1, *N2, *M);  
}

void FC_FUNC(oct_nfft_init_3d,OCT_NFFT_INIT_3D)
   (nfft_plan *plan, int *N1, int *N2, int *N3, int *M)
{
    nfft_init_3d(plan, *N1, *N2, *N3, *M);  
}

void FC_FUNC(oct_nfft_init_guru,OCT_NFFT_INIT_GURU) 
  (nfft_plan *ths, int *d, int *N, int *M, int *n, int *m, unsigned *nfft_flags, unsigned *fftw_flags)
{
  nfft_init_guru (ths, *d, N, *M, n, *m, *nfft_flags, *fftw_flags);
}

void 	FC_FUNC(oct_nfft_check, OCT_NFFT_CHECK) 
  (nfft_plan *ths)
{
  nfft_check (ths);
}

void 	FC_FUNC(oct_nfft_finalize,OCT_NFFT_FINALIZE) 
  (nfft_plan *ths)	
{
  nfft_finalize (ths);
}

void 	FC_FUNC(oct_nfft_trafo,OCT_NFFT_TRAFO)
   (nfft_plan *ths)
{
  nfft_trafo (ths);
}  
  
void 	FC_FUNC(oct_nfft_adjoint,OCT_NFFT_ADJOINT)
   (nfft_plan *ths)
{
  nfft_adjoint (ths);
}
  
void FC_FUNC(oct_nfft_precompute_one_psi_1d, OCT_NFFT_PRECOMPUTE_ONE_PSI_1D)
     (nfft_plan *plan, int *m, double *X1)
{

   int ii;
   int M = *m;
 
   for (ii=0;ii< M;ii++){
     plan->x[ii] =  X1[ii];
   }
  
   if(plan->nfft_flags & PRE_ONE_PSI)
     nfft_precompute_one_psi(plan);
 
}

void FC_FUNC(oct_nfft_precompute_one_psi_2d, OCT_NFFT_PRECOMPUTE_ONE_PSI_2D)
     (nfft_plan *plan, int *m, double* X1, double* X2)
{
   int ii;
   int jj;
   int M = *m;
   
   for (ii=0; ii< M; ii++){
     for (jj=0; jj< M; jj++){
       plan->x[2*(M * ii + jj) + 0] =  X1[ii];
       plan->x[2*(M * ii + jj) + 1] =  X2[jj];
     }
   }


   if(plan->nfft_flags & PRE_ONE_PSI)
     nfft_precompute_one_psi(plan);
}

void FC_FUNC(oct_nfft_precompute_one_psi_3d, OCT_NFFT_PRECOMPUTE_ONE_PSI_3D)
     (nfft_plan *plan, int *m, double* X1, double* X2, double* X3)
{
   int ii,jj,kk;
   int M = *m;

   for (ii=0;ii< M;ii++){
     for (jj=0;jj< M;jj++){
       for (kk=0;kk< M;kk++){
         plan->x[3*(M*M*ii+M*jj+kk) + 0] =  X1[ii];
         plan->x[3*(M*M*ii+M*jj+kk) + 1] =  X2[jj];
         plan->x[3*(M*M*ii+M*jj+kk) + 2] =  X3[kk];
       }
     }
   }
  
   if(plan->nfft_flags & PRE_ONE_PSI)
     nfft_precompute_one_psi(plan);
}


void FC_FUNC(oct_set_f, OCT_SET_F)
    (nfft_plan *plan, int *m, int *DIM, double complex *VAL, int *IX, int *IY, int *IZ)
{
  int M = *m;
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
       plan->f[(ix-1)*M + (iy-1)] = val;
     break;
     case 3:
       plan->f[(ix-1)*M*M + (iy-1)*M + (iz-1)] = val;
     break;
   }

}

void FC_FUNC(oct_get_f, OCT_GET_F)
    (nfft_plan *plan, int *m, int *DIM, double complex *val, int *IX, int *IY, int *IZ)
{
  int M = *m;
  int dim = *DIM;
  int ix = *IX;
  int iy = *IY;
  int iz = *IZ;
   

   switch (dim){
     case 1:
       *val = plan->f[ix-1];
     break;
     case 2:
       *val = plan->f[(ix-1)*M + (iy-1)];
     break;
     case 3:
       *val = plan->f[(ix-1)*M*M + (iy-1)*M + (iz-1)];
     break;
   }

}


void FC_FUNC(oct_set_f_hat, OCT_SET_F_HAT)
    (nfft_plan *plan, int *m, int *DIM, double complex *VAL, int *IX, int *IY, int *IZ)
{
  int M = *m;
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

}

void FC_FUNC(oct_get_f_hat, OCT_GET_F_HAT)
    (nfft_plan *plan, int *m, int *DIM, double complex *val, int *IX, int *IY, int *IZ)
{
  int M = *m;
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
       *val = plan->f_hat[(ix-1)*plan->N[1]*plan->N[2] + (iy-1)*plan->N[1] + (iz-1)];
     break;
   }

}

#endif 
