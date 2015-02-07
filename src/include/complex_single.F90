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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#define R_TCOMPLEX 1
#define SINGLE_PRECISION 1

#define R_TYPE      complex(4)
#define R_BASE      real(4)
#define R_SINGLE    complex(4)
#define R_DOUBLE    complex(8)
#define R_MPITYPE   MPI_COMPLEX
#define R_TYPE_VAL  TYPE_CMPLX_SINGLE
#define R_TYPE_CL   'RTYPE_COMPLEX_SINGLE'
#define R_TYPE_IOBINARY TYPE_DOUBLE_COMPLEX
#define R_TOTYPE(x) cmplx(x, 0.0_4, 4)
#define R_TOPREC(x) cmplx(real(x), aimag(x), 4)

#define R_ABS(x)    abs(x)
#define R_CONJ(x)   conjg(x)
#define R_REAL(x)   real(x, 4)
#define R_AIMAG(x)  aimag(x)

#define R_SIZEOF    8
#define R_ADD       2
#define R_MUL       6

#define X(x)        c ## x
#define pX(x)       pc ## x
#define aX(x,y)     x ## c ## y

#if defined(DISABLE_DEBUG)
#define TS(x)       x
#else
#define TS(x)       TSZ_ ## x
#endif


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
