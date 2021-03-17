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

#define R_TINTEGER8  1

#define X(x)        l ## x

#define R_TYPE      integer(8)
#define R_BASE      integer(8)
#define R_TYPE_VAL  TYPE_INTEGER8
#define R_MPITYPE   MPI_LONG_LONG
#define R_TYPE_IOBINARY TYPE_INT_64
#define R_TOTYPE(x) (x)

#define R_SIZEOF    8

#if defined(DISABLE_DEBUG)
#define TS(x)       x
#else
#define TS(x)       TSL_ ## x
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
