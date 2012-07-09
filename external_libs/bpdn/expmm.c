/*
!! Copyright (C) 2012 X. Andrade
!! 
!! This program is a free translation to Fortran of part of the spgl1
!! library, version 1.0.
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! $Id$
*/

#include <config.h>
#include <math.h>

void FC_FUNC(expmm, EXPMM)(const int * restrict nx, const int * restrict ny, 
			   const double * restrict x, double * restrict y,
			   const double * restrict dx, const double * restrict dy,
			   const int * restrict component){
  int ix, iy;
  double dexpr, dexpi;
  double aexpr, aexpi;
  double aexp[2];
  double tt;

  for(iy = 0; iy < *ny; iy++){
    dexpr = cos((*dx)*(*dy)*iy);
    dexpi = sin((*dx)*(*dy)*iy);
    aexp[0] = 1.0;
    aexp[1] = 0.0;
    tt = 0.0;

    for(ix = 0; ix < *nx; ix++){
      tt += aexp[*component]*x[ix];

      aexpr = aexp[0];
      aexpi = aexp[1];
      aexp[0] = dexpr*aexpr - dexpi*aexpi;
      aexp[1] = dexpr*aexpi + dexpi*aexpr;
    }
    y[iy] = -tt;
  }
  
}
