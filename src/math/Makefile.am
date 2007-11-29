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
	blas.F90 dablas.c ylm.c oct_gsl_f.c loct_math.F90 lalg_basic.F90 lalg_adv.F90 gauss_legendre.F90 math.F90 \
	ode_solver.F90 root_solver.F90 spline.F90 sparskit.F90 fft.F90

noinst_HEADERS = \
	lalg_basic_blas.F90 lalg_adv_lapack.F90 math_cg_inc.F90 math_inc.F90 \
	ode_solver_inc.F90 root_solver_inc.F90 sparskit_inc.F90 \
	fftw2.F90 fftw3.F90 lalg_basic_mp.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libmath_a_SOURCES)

blas.o : blas.F90 
dablas.o : dablas.c 
lalg_basic.o : lalg_basic.F90 lalg_basic_mp.F90 lalg_basic_blas.F90 blas.o 
lalg_adv.o : lalg_adv.F90 
gauss_legendre.o : gauss_legendre.F90 
math.o : math.F90 math_inc.F90 math_cg_inc.F90 lalg_basic.o loct_math.o
ode_solver.o : ode_solver.F90 ode_solver_inc.F90 
root_solver.o : root_solver.F90 root_solver_inc.F90 ode_solver.o lalg_adv.o 
spline.o : spline.F90 loct_math.o
sparskit.o : sparskit.F90 
fft.o : fft.F90 fftw3.F90 fftw2.F90 lalg_basic.o loct_math.o
oct_gsl_f.o : oct_gsl_f.c
ylm.o : ylm.c