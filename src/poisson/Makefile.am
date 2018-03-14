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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##

noinst_LIBRARIES = libpoisson.a

# Fortran sources - keep in alphabetical order
libpoisson_f_srcs =                      \
	poisson_cg.F90                   \
	poisson_corrections.F90          \
	poisson_cutoff.F90               \
	poisson.F90                      \
	poisson_fft.F90                  \
	poisson_fmm.F90                  \
	poisson_isf.F90                  \
	poisson_libisf.F90               \
	poisson_multigrid.F90            \
	poisson_no.F90                   \
	scaling_function.F90

# C/C++ sources - keep in alphabetical order
libpoisson_c_srcs =                      \
	poisson_cutoffs.c

libpoisson_a_SOURCES = $(libpoisson_c_srcs) $(libpoisson_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	lazy_8_inc.F90                   \
	poisson_inc.F90                  \
	poisson_direct_sm_inc.F90        \
	poisson_init_direct_inc.F90      \
	solver_1d_solve_inc.F90


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	makedepf90 $(libpoisson_a_SOURCES)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

poisson_cutoffs.o : poisson_cutoffs.c 
poisson_cg.o : poisson_cg.F90 poisson_corrections.o 
poisson_corrections.o : poisson_corrections.F90 
poisson_cutoff.o : poisson_cutoff.F90 
poisson.o : poisson.F90 solver_1d_solve_inc.F90 poisson_inc.F90 poisson_init_direct_inc.F90 poisson_no.o poisson_multigrid.o poisson_libisf.o poisson_fmm.o poisson_fft.o poisson_isf.o poisson_corrections.o poisson_cg.o poisson_direct_sm_inc.F90 
poisson_fft.o : poisson_fft.F90 poisson_cutoff.o 
poisson_fmm.o : poisson_fmm.F90 
poisson_isf.o : poisson_isf.F90 scaling_function.o 
poisson_libisf.o : poisson_libisf.F90 
poisson_multigrid.o : poisson_multigrid.F90 poisson_corrections.o 
poisson_no.o : poisson_no.F90 
scaling_function.o : scaling_function.F90 lazy_8_inc.F90
