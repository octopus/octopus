/*
 Copyright (C) 2006 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: beak.h 2146 2006-05-23 17:36:00Z xavier $
*/

#ifndef OCTOPUS_BEAK_H
#define OCTOPUS_BEAK_H

#include <config.h>

/* If __builtin_prefetch is not present (which should have been caught
   by the configure script) one needs to define a dummy preprocessor
   macro. */

#if !defined(HAVE_BUILTIN_PREFETCH)
#define __builtin_prefetch(a, b, c)
#endif

#ifdef SINGLE_PRECISION
typedef float ffloat;
#else
typedef double ffloat;
#endif

typedef struct {
  ffloat re;
  ffloat im;
} comp;

/* These constants have to match the definition in
   src/grid/nl_operator.F90 */

#define OP_FORTRAN 0
#define OP_VEC     1

#define TYPE_REAL     1
#define TYPE_CMPLX    2

#define MAX_OP_N   400

#endif /* OCTOPUS_BEAK_H */
