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
  double dexpr0, dexpi0;
  double aexpr0, aexpi0;
  double dexpr1, dexpi1;
  double aexpr1, aexpi1;
  double aexp0[2];
  double aexp1[2];
  double xx, tt0, tt1;

  for(iy = 0; iy < (*ny) - 2 + 1; iy += 2){
    dexpr0 = cos((*dx)*(*dy)*iy);
    dexpr1 = cos((*dx)*(*dy)*(iy + 1));
    dexpi0 = sin((*dx)*(*dy)*iy);
    dexpi1 = sin((*dx)*(*dy)*(iy + 1));

    aexp0[0] = 1.0;
    aexp1[0] = 1.0;
    aexp0[1] = 0.0;
    aexp1[1] = 0.0;

    tt0 = 0.0;
    tt1 = 0.0;

    for(ix = 0; ix < *nx; ix++){
      xx = x[ix];

      tt0 += aexp0[*component]*xx;
      tt1 += aexp1[*component]*xx;

      aexpr0 = aexp0[0];
      aexpi0 = aexp0[1];
      aexpr1 = aexp1[0];
      aexpi1 = aexp1[1];

      aexp0[0] = dexpr0*aexpr0 - dexpi0*aexpi0;
      aexp0[1] = dexpr0*aexpi0 + dexpi0*aexpr0;
      aexp1[0] = dexpr1*aexpr1 - dexpi1*aexpi1;
      aexp1[1] = dexpr1*aexpi1 + dexpi1*aexpr1;
    }
    y[iy    ] = -tt0;
    y[iy + 1] = -tt1;
  }

  for(; iy < *ny; iy++){
    dexpr0 = cos((*dx)*(*dy)*iy);
    dexpi0 = sin((*dx)*(*dy)*iy);
    aexp0[0] = 1.0;
    aexp0[1] = 0.0;
    tt0 = 0.0;

    for(ix = 0; ix < *nx; ix++){
      tt0 += aexp0[*component]*x[ix];

      aexpr0 = aexp0[0];
      aexpi0 = aexp0[1];
      aexp0[0] = dexpr0*aexpr0 - dexpi0*aexpi0;
      aexp0[1] = dexpr0*aexpi0 + dexpi0*aexpr0;
    }
    y[iy] = -tt0;
  }
  
}
