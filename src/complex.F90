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

#define R_TCOMPLEX 1
#define R_FUNC(x) z ## x
#define X(x) z ## x
#define R_TYPE complex(r8)
#define R_MPITYPE MPI_DOUBLE_COMPLEX
#define R_ABS(x) abs(x)
#define R_CONJ(x) conjg(x)
#define R_REAL(x) real(x, r8)
#define R_AIMAG(x) aimag(x)
#define R_DOT zdotc
#define R_NRM2 dznrm2
#define R_TOTYPE(x) cmplx(x, 0.0_r8, r8)
