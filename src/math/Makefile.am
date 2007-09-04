## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id: Makefile.am 2995 2007-06-13 17:49:22Z xavier $

noinst_LIBRARIES = libmath.a

libmath_a_SOURCES = \
	blas.F90 dablas.c linalg_basic.F90 linalg_adv.F90 gauss_legendre.F90 math.F90 \
	ode_solver.F90 root_solver.F90 spline.F90 sparskit.F90 fft.F90

noinst_HEADERS = \
	linalg_basic_blas.F90 linalg_adv_lapack.F90 math_cg_inc.F90 math_inc.F90 \
	ode_solver_inc.F90 root_solver_inc.F90 sparskit_inc.F90 \
	fftw2.F90 fftw3.F90 linalg_basic_mp.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libmath_a_SOURCES)

linalg_adv.o: linalg_adv_lapack.F90
linalg_basic.o: linalg_basic_blas.F90 linalg_basic_mp.F90
sparskit.o: sparskit_inc.F90
math.o: math_cg_inc.F90 math_inc.F90
