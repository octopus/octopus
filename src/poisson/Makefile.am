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
	solver_isf.F90 solver_fft.F90 \
	corrections.F90 conj_grad.F90 multigrid.F90 poisson.F90

noinst_HEADERS = \
	solver_1D.F90 solver_2D.F90 solver_3D.F90

AM_CPPFLAGS = -I$(srcdir)/../include -I$(top_builddir)/src/include

AM_FCFLAGS = \
	@F90_MODULE_FLAG@$(top_builddir)/src/basic    \
	@F90_MODULE_FLAG@$(top_builddir)/src/math     \
	@F90_MODULE_FLAG@$(top_builddir)/src/species  \
	@F90_MODULE_FLAG@$(top_builddir)/src/grid

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libspecies_a_SRCS)

cube_functions.o:  cube_functions_inc.F90
derivatives.o:  derivatives_inc.F90
double_grid.o:  double_grid_apply.F90
functions.o: functions_inc.F90
mesh.o:  mesh_init.F90
mesh_function.o:  mesh_function_inc.F90
out.o:  out_inc.F90
par_vec.o:  par_vec_inc.F90
submesh.o:  submesh_inc.F90
