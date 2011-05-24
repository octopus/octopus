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

noinst_LIBRARIES = libpoisson.a

libpoisson_a_SOURCES = \
	poisson_cutoffs.c poisson_cutoff.F90 poisson_isf.F90 poisson_fft.F90 \
	poisson_corrections.F90 poisson_cg.F90 poisson_multigrid.F90 poisson.F90\
	poisson_sete.F90 scaling_function.F90

noinst_HEADERS = \
	solver_1d_inc.F90    \
	solver_2d_inc.F90    \
	solver_3d_inc.F90    \
	poisson_fmm_inc.F90  \
	lazy_8_inc.F90

include $(srcdir)/../common-rules.make

depend:
	makedepf90 $(libpoisson_a_SOURCES)

poisson_sete.o : poisson_sete.F90
poisson_cutoffs.o : poisson_cutoffs.c
poisson_cutoff.o : poisson_cutoff.F90 
poisson_isf.o : poisson_isf.F90 
poisson_fft.o : poisson_fft.F90 poisson_cutoff.o 
poisson_corrections.o : poisson_corrections.F90 
poisson_cg.o : poisson_cg.F90 poisson_corrections.o 
poisson_multigrid.o : poisson_multigrid.F90 poisson_corrections.o 
poisson.o : poisson.F90 poisson_fmm_inc.F90 solver_3d_inc.F90 solver_2d_inc.F90 solver_1d_inc.F90 poisson_multigrid.o poisson_fft.o poisson_isf.o poisson_cg.o poisson_corrections.o poisson_sete.o
scaling_function.o : scaling_function.F90 lazy_8_inc.F90
poisson_isf.o : poisson_isf.F90 scaling_function.o
