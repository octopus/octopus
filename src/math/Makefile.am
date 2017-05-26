## Process this file with automake to produce Makefile.in

## Copyright (C) 2002-2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##

noinst_LIBRARIES = libmath.a

libmath_a_SOURCES =                      \
	blas.F90                         \
	checksum.c                       \
	checksum_interface.F90           \
	dablas.c                         \
	fft.F90                          \
	fftw.F90                         \
	spline_filter.F90                \
	gauss_legendre.F90               \
	loct_math.F90                    \
	lalg_basic.F90                   \
	math.F90                         \
	metis.F90                        \
	metis_f.c                        \
	lalg_adv.F90                     \
	lapack.F90                       \
	oct_gsl_f.c                      \
	minimizer_low.c                  \
	ode_solver.F90                   \
	permutations.F90                 \
	pfft.F90			 \
	root_solver.F90                  \
	solvers.F90                      \
	splines.F90                      \
	sparskit.F90                     \
	symm_op.F90                      \
	tdfunction.F90                   \
	ylm.c                            \
	young.F90                        \
	scalapack.F90                    \
	nfft.F90    			 \
	nfft_f.c                         \
	pblas.F90                        \
	sgfft.F90                        \
	matrix.F90                       \
	matrix_ops.F90                   \
	compressed_sensing.F90           \
	fftw_low.c                       \
	minimizer.F90                    \
        pnfft.F90                        \
	qshep.F90                        \
	spline_low.cc                    \
	quickrnd.F90                     \
	accel_blas.F90                   \
	cufft.cc 

noinst_HEADERS =                         \
	fft_inc.F90                      \
	lalg_basic_blas_inc.F90          \
	lalg_adv_lapack_inc.F90          \
	math_inc.F90                     \
	ode_solver_inc.F90               \
	root_solver_inc.F90              \
	solvers_inc.F90                  \
	sparskit_inc.F90                 \
	symm_op_inc.F90                  \
	matrix_inc.F90                   \
	matrix_ops_inc.F90               \
	nfft_inc.F90                     \
	pnfft_inc.F90                    \
	accel_blas_inc.F90

include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libmath_a_SOURCES)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

blas.o : blas.F90 
checksum.o : checksum.c 
checksum_interface.o : checksum_interface.F90 
dablas.o : dablas.c 
fft.o : fft.F90 fft_inc.F90 pnfft.o pfft.o nfft.o loct_math.o lalg_basic.o fftw.o 
fftw.o : fftw.F90 
spline_filter.o : spline_filter.F90 loct_math.o splines.o 
gauss_legendre.o : gauss_legendre.F90 
loct_math.o : loct_math.F90 
lalg_basic.o : lalg_basic.F90 lalg_basic_blas_inc.F90 blas.o 
math.o : math.F90 math_inc.F90 loct_math.o lalg_basic.o 
metis.o : metis.F90 
metis_f.o : metis_f.c 
lalg_adv.o : lalg_adv.F90 lalg_adv_lapack_inc.F90 scalapack.o 
lapack.o : lapack.F90 
oct_gsl_f.o : oct_gsl_f.c 
minimizer_low.o : minimizer_low.c 
ode_solver.o : ode_solver.F90 ode_solver_inc.F90 
permutations.o : permutations.F90 math.o 
pfft.o : pfft.F90 math.o fftw.o 
root_solver.o : root_solver.F90 root_solver_inc.F90 ode_solver.o lalg_adv.o
solvers.o : solvers.F90 solvers_inc.F90 loct_math.o lalg_basic.o blas.o
splines.o : splines.F90 loct_math.o
sparskit.o : sparskit.F90
symm_op.o : symm_op.F90 symm_op_inc.F90
tdfunction.o : tdfunction.F90 splines.o math.o loct_math.o fft.o
ylm.o : ylm.c
young.o : young.F90 math.o
scalapack.o : scalapack.F90
nfft.o : nfft.F90
nfft_f.o : nfft_f.c
pblas.o : pblas.F90
sgfft.o : sgfft.F90
matrix.o : matrix.F90 matrix_inc.F90
matrix_ops.o : matrix_ops.F90 matrix_ops_inc.F90 scalapack.o matrix.o
compressed_sensing.o : compressed_sensing.F90 blas.o
fftw_low.o : fftw_low.c
minimizer.o : minimizer.F90 lalg_basic.o
pnfft.o : pnfft.F90 pfft.o math.o loct_math.o pfft.o
qshep.o : qshep.F90
spline_low.o : spline_low.cc
quickrnd.o : quickrnd.F90
accel_blas.o : accel_blas.F90 accel_blas_inc.F90
cufft.o : cufft.cc
