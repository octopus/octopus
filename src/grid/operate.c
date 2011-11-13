/*
 Copyright (C) 2006 X. Andrade

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

 $Id: operate.c 2146 2006-05-23 17:36:00Z xavier $
*/


#include <config.h>

#include <assert.h>
#include <vectors.h>

void FC_FUNC_(operate_ri_vec,OPERATE_RI_VEC)(const int * opn, 
					       const double * restrict w, 
					       const int * opnri,
					       const int * opri,
					       const int * rimap_inv,
					       const int * rimap_inv_max,
					       const double * restrict fi,
					       const int * ldfp,
					       double * restrict fo){
  const int ldf = ldfp[0];

  /* check whether we got aligned vectors or not */
  int aligned = 1;
  aligned = aligned && (((long long) fi)%(8*VEC_SIZE) == 0);
  aligned = aligned && (((long long) fo)%(8*VEC_SIZE) == 0);
  aligned = aligned && ((1<<ldf)%VEC_SIZE == 0);

  if(aligned){
#define ALIGNED
#include "operate_inc.c"
#undef ALIGNED

  } else {
    /* not aligned */
#include "operate_inc.c"
  }

}

void FC_FUNC_(dgauss_seidel,DGAUSS_SEIDEL)(const int * opn, 
				           const double * restrict w, 
				           const int * opnri,
				           const int * opri,
				           const int * rimap_inv,
				           const int * rimap_inv_max,
					   double * restrict factor,
				           double * pot, 
				           double * restrict rho){

  const int n = opn[0];
  const int nri = opnri[0];

  int l, i, j;
  const int * index;
  register double a0;
  register const double fac = factor[0];

  for (l = 0; l < nri ; l++) {
    index  = opri + n*l;
    i = rimap_inv[l];

    for (; i < rimap_inv_max[l]; i++){
      a0 = 0.0;
      for(j = 0; j < n; j++) a0 += w[j]*(pot + index[j])[i];
      pot[i] += fac*(a0 - rho[i]);
    }

  }

}

