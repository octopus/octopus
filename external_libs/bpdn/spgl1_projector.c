/*
 Copyright (C) 2011 X. Andrade
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include <config.h>
#include <stdlib.h>
#include <math.h>
#include "oneProjectorCore.h"

void FC_FUNC_(spgl1_projector, SPGL1_PROJECTOR)
     (double * restrict xPtr, double * restrict bPtr, 
      double * restrict tau, int* restrict nn){
  int ii;
  double * signs;

  /* we need to do some preprocessing, copy the input into the output
     vector and remove the sign */
  signs = (double *) malloc(sizeof(double)*(*nn));
  for(ii = 0; ii < *nn; ii++){
    double bb = bPtr[ii];
    signs[ii] = (bPtr[ii] >= 0.0)?1.0:-1.0;
    xPtr[ii] = fabs(bb);
    bPtr[ii] = fabs(bb);
  }

  projectI(xPtr, bPtr, *tau, *nn);

  /* now we put the sign back */
  for(ii = 0; ii < *nn; ii++){
    xPtr[ii] = xPtr[ii]*signs[ii];
    bPtr[ii] = bPtr[ii]*signs[ii];
  }
  
  free(signs);
}
