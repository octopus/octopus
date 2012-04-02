!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#define R_TCOMPLEX 1

#define R_TYPE      CMPLX
#define R_SINGLE    complex(4)
#define R_DOUBLE    complex(8)
#define R_MPITYPE   MPI_CMPLX
#define R_TYPE_VAL  TYPE_CMPLX
#define R_TOTYPE(x) cmplx(x, M_ZERO, REAL_PRECISION)
#define R_TOPREC(x) cmplx(real(x), aimag(x), REAL_PRECISION)

#define R_ABS(x)    abs(x)
#define R_CONJ(x)   conjg(x)
#define R_REAL(x)   real(x, REAL_PRECISION)
#define R_AIMAG(x)  aimag(x)

#define R_SIZEOF    16
#define R_ADD       2
#define R_MUL       6

#define X(x)        z ## x
#define pX(x)        pz ## x
#define aX(x,y)        x ## z ## y

#if defined(DISABLE_DEBUG)
#define TS(x)       x
#else
#define TS(x)       TSZ_ ## x
#endif


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
