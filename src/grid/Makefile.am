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

libgrid_a_SOURCES =       \
	curv_gygi.F90           \
	curv_briggs.F90         \
	curv_modine.F90         \
	curvilinear.F90         \
	index.F90               \
	par_vec.F90             \
	mesh.F90                \
	nl_operator.F90         \
	submesh.F90             \
	stencil_star.F90        \
	stencil_variational.F90 \
	stencil_cube.F90        \
	stencil.F90             \
	stencil_starplus.F90    \
	derivatives.F90         \
	cube_function.F90       \
	fourier_space.F90       \
	mesh_function.F90       \
	basins.F90              \
	io_function.F90         \
	multigrid.F90           \
	ob_grid.F90             \
	gridhier.F90            \
	double_grid.F90         \
	grid.F90                \
	batch.F90               \
	mesh_init.F90           \
	mesh_partition.F90      \
	metis.c                 \
	ob_interface.F90        \
	partition.F90           \
	partitioner.F90         \
	zoltan_low.c            \
	zoltan.F90              \
	hypercube.F90           \
	modelmb_1part.F90       \
	transfer_table.F90      \
	mesh_batch.F90          \
	subarray.F90            \
	symmetrizer.F90         \
	boundaries.F90          \
	operate.c               \
	mesh_cube_map.F90       \
	fourier_shell.F90


noinst_HEADERS =                  \
	par_vec_inc.F90           \
	derivatives_inc.F90       \
	cube_function_inc.F90     \
	fourier_space_inc.F90     \
	mesh_function_inc.F90     \
	double_grid_apply_inc.F90 \
	submesh_inc.F90           \
	io_function_inc.F90       \
	nl_operator_inc.F90       \
	batch_inc.F90             \
	multigrid_inc.F90         \
	gridhier_inc.F90          \
	mesh_batch_inc.F90        \
	subarray_inc.F90          \
	symmetrizer_inc.F90       \
	boundaries_inc.F90        \
	operate_inc.c

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libgrid_a_SOURCES)

curv_gygi.o : curv_gygi.F90  
curv_briggs.o : curv_briggs.F90  
curv_modine.o : curv_modine.F90  
curvilinear.o : curvilinear.F90  curv_modine.o curv_gygi.o curv_briggs.o 
index.o : index.F90 hypercube.o
par_vec.o : par_vec.F90 par_vec_inc.F90 index.o stencil.o subarray.o
mesh.o : mesh.F90 par_vec.o index.o curvilinear.o mesh_cube_map.o
ob_grid.o : ob_grid.F90 mesh.o ob_interface.o
mesh_init.o : mesh_init.F90 mesh.o par_vec.o stencil.o ob_grid.o partition.o partitioner.o mesh_partition.o mesh_cube_map.o
mesh_partition.o : mesh_partition.F90 mesh.o zoltan.o stencil_star.o stencil.o partition.o partitioner.o 
modelmb_1part.o : modelmb_1part.F90 mesh.o hypercube.o
nl_operator.o : nl_operator.F90 nl_operator_inc.F90 par_vec.o mesh.o index.o batch.o stencil.o boundaries.o
submesh.o : submesh.F90 submesh_inc.F90  mesh.o batch.o
stencil_star.o : stencil_star.F90 nl_operator.o stencil.o
stencil_variational.o : stencil_variational.F90 nl_operator.o stencil.o
stencil_cube.o : stencil_cube.F90 stencil.o
stencil_starplus.o : stencil_starplus.F90 stencil.o
derivatives.o : derivatives.F90 derivatives_inc.F90 stencil_variational.o stencil_starplus.o stencil_star.o stencil_cube.o  nl_operator.o mesh.o batch.o transfer_table.o mesh_function.o boundaries.o
cube_function.o : cube_function.F90 cube_function_inc.F90  mesh.o mesh_cube_map.o par_vec.o
fourier_space.o : fourier_space.F90 fourier_space_inc.F90 cube_function.o mesh.o 
mesh_function.o : mesh_function.F90 mesh_function_inc.F90 par_vec.o mesh.o cube_function.o modelmb_1part.o
basins.o : basins.F90 par_vec.o index.o mesh.o 
io_function.o : io_function.F90 io_function_inc.F90 par_vec.o mesh.o mesh_function.o cube_function.o 
multigrid.o : multigrid.F90 multigrid_inc.F90 mesh.o curvilinear.o derivatives.o mesh_init.o transfer_table.o boundaries.o
double_grid.o : double_grid.F90 double_grid_apply_inc.F90 submesh.o mesh_function.o mesh.o curvilinear.o 
grid.o : grid.F90  multigrid.o mesh.o double_grid.o curvilinear.o ob_interface.o mesh_init.o nl_operator.o
gridhier.o : gridhier.F90 gridhier_inc.F90 multigrid.o
ob_interface.o : ob_interface.F90 derivatives.o mesh.o
zoltan_low.o : zoltan_low.c
zoltan.o : zoltan.F90
batch.o : batch.F90 batch_inc.F90
stencil.o : stencil.F90
hypercube.o : hypercube.F90
transfer_table.o : transfer_table.F90
mesh_batch.o : mesh_batch.F90 mesh_batch_inc.F90 mesh.o mesh_function.o
symmetrizer.o : symmetrizer.F90 symmetrizer_inc.F90 mesh.o
partition.o : partition.F90 mesh.o
partitioner.o : partitioner.F90 partition.o
subarray.o : subarray.F90 subarray_inc.F90 batch.o
boundaries.o : boundaries.F90 boundaries_inc.F90 subarray.o batch.o par_vec.o mesh.o
operate.o : operate.c operate_inc.c
mesh_cube_map.o : mesh_cube_map.F90 index.o
fourier_shell.o : fourier_shell.F90 cube_function.o mesh.o
