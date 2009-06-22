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

noinst_LIBRARIES = libgrid.a

libgrid_a_SOURCES = \
	curv_gygi.F90 curv_briggs.F90 curv_modine.F90 curvilinear.F90 \
	index.F90 \
	par_vec.F90 mesh.F90 scalar_mesh.F90 nl_operator.F90 \
	submesh.F90 stencil_star.F90 stencil_variational.F90 stencil_cube.F90 \
	stencil.F90        \
	stencil_starplus.F90 derivatives.F90 cube_function.F90 fourier_space.F90 mesh_function.F90 \
	basins.F90         \
	io_function.F90    \
	multigrid.F90      \
	gridhier.F90       \
	double_grid.F90    \
	grid.F90           \
	batch.F90          \
	mesh_init.F90      \
	metis.c            \
	ob_interface.F90   \
	zoltan_low.c       \
	zoltan.F90         \
	hypercube.F90      \
        modelmb_1part.F90  \
	cubic_mesh.F90     \
	cubic_mesh_low.c   \
	cubic_mesh_fft.c   \
	transfer_table.F90 \
	fourier_kernel.F90 \
	mesh_batch.F90

noinst_HEADERS = \
	par_vec_inc.F90 derivatives_inc.F90 cube_function_inc.F90 fourier_space_inc.F90 \
	mesh_function_inc.F90 double_grid_apply.F90 \
	submesh_inc.F90 io_function_inc.F90 nl_operator_inc.F90 \
	batch_inc.F90       \
	multigrid_inc.F90   \
	gridhier_inc.F90    \
	cubic_mesh.h        \
	mesh_batch_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libgrid_a_SOURCES)

curv_gygi.o : curv_gygi.F90  
curv_briggs.o : curv_briggs.F90  
curv_modine.o : curv_modine.F90  
curvilinear.o : curvilinear.F90  curv_modine.o curv_gygi.o curv_briggs.o 
index.o : index.F90 hypercube.o
par_vec.o : par_vec.F90 par_vec_inc.F90 index.o stencil.o
mesh.o : mesh.F90 par_vec.o index.o curvilinear.o
mesh_init.o : mesh_init.F90 mesh.o par_vec.o zoltan.o stencil_star.o stencil.o
modelmb_1part.o : modelmb_1part.F90 mesh.o hypercube.o
scalar_mesh.o : scalar_mesh.F90
nl_operator.o : nl_operator.F90 nl_operator_inc.F90 par_vec.o mesh.o index.o batch.o stencil.o
submesh.o : submesh.F90 submesh_inc.F90  mesh.o batch.o
stencil_star.o : stencil_star.F90 nl_operator.o stencil.o
stencil_variational.o : stencil_variational.F90 nl_operator.o stencil.o
stencil_cube.o : stencil_cube.F90 stencil.o
stencil_starplus.o : stencil_starplus.F90 stencil.o
derivatives.o : derivatives.F90 derivatives_inc.F90 stencil_variational.o stencil_starplus.o stencil_star.o stencil_cube.o  nl_operator.o mesh.o batch.o
cube_function.o : cube_function.F90 cube_function_inc.F90  mesh.o 
fourier_space.o : fourier_space.F90 fourier_space_inc.F90 cube_function.o mesh.o 
mesh_function.o : mesh_function.F90 mesh_function_inc.F90 par_vec.o mesh.o cube_function.o modelmb_1part.o
basins.o : basins.F90 par_vec.o index.o mesh.o 
io_function.o : io_function.F90 io_function_inc.F90 par_vec.o mesh.o  cube_function.o 
multigrid.o : multigrid.F90 multigrid_inc.F90 mesh.o curvilinear.o derivatives.o mesh_init.o transfer_table.o
double_grid.o : double_grid.F90 double_grid_apply.F90 submesh.o mesh_function.o mesh.o curvilinear.o 
grid.o : grid.F90  multigrid.o mesh.o double_grid.o curvilinear.o ob_interface.o mesh_init.o
gridhier.o : gridhier.F90 gridhier_inc.F90 multigrid.o
ob_interface.o : ob_interface.F90 derivatives.o mesh.o
zoltan_low.o : zoltan_low.c
zoltan.o : zoltan.F90
batch.o : batch.F90 batch_inc.F90
stencil.o : stencil.F90
hypercube.o : hypercube.F90
cubic_mesh.o : cubic_mesh.F90 mesh.o
cubic_mesh_low.o : cubic_mesh_low.c cubic_mesh.h
cubic_mesh_fft.o : cubic_mesh_fft.c cubic_mesh.h
transfer_table.o : transfer_table.F90
mesh_batch.o : mesh_batch.F90 mesh_batch_inc.F90 mesh.o
fourier_kernel.F90 : fourier_kernel.o cubic_mesh.o
