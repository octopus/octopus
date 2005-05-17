!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#define R_TREAL     1

#define X(x)        d ## x

#define R_TYPE      FLOAT
#define R_MPITYPE   MPI_FLOAT
#define R_TOTYPE(x) real(x, PRECISION)

#define R_ABS(x)    abs(x)
#define R_CONJ(x)   (x)
#define R_REAL(x)   (x)
#define R_AIMAG(x)  (M_ZERO)
