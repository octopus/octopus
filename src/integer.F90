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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#define R_TINTEGER  1

#define X(x)        i ## x

#define R_TYPE      integer
#define R_MPITYPE   MPI_INTEGER
#define R_TOTYPE(x) (x)

#define R_ABS(x)    abs(x)
#define R_CONJ(x)   (x)
#define R_REAL(x)   (x)
#define R_AIMAG(x)  (M_ZERO)

#define TS(x)       TSI_ ## x
