!! Copyright (C) 2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "config_F90.h"

#define NP      gr%m%np
#define NP_PART gr%m%np_part

#define NDIM    gr%sb%dim
#define LAP     f_der%der_discr%lapl

#define __STRING(x)     #x

#if !defined(NDEBUG) && defined(LONG_LINES)
#  define ASSERT(expr) \
  if(.not.(expr)) call assert_die (__STRING(expr), __FILE__, __LINE__)
#else
#  define ASSERT(expr)
!!#  define IN inout
#endif

#define DOUBLE real(8)
#define SINGLE real(4)

#if defined(SINGLE_PRECISION)
#  define PRECISION 4
#  define FLOAT     real(4)
#  define MPI_FLOAT MPI_REAL
#  define CMPLX     complex(4)
#  define MPI_CMPLX MPI_COMPLEX
#  define PREC(x)   s ## x
#  define CNST(x)   x ## _4
#else
#  define PRECISION 8
#  define FLOAT     real(8)
#  define MPI_FLOAT MPI_DOUBLE_PRECISION
#  define CMPLX     complex(8)
#  define MPI_CMPLX MPI_DOUBLE_COMPLEX
#  define PREC(x)   d ## x
#  define CNST(x)   x ## _8
#endif

! what do you wish for dinner, dear?
#if defined(COMPLEX_WFNS)
#  include "complex.F90"
#else
#  include "real.F90"
#endif
