## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
## Copyright (C) 2012 M. Oliveira
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

noinst_LIBRARIES = libgrid.a

# Fortran sources - keep in alphabetical order
libgrid_f_srcs =                         \
	basins.F90                       \
	batch.F90                        \
	batch_ops.F90                    \
	boundaries.F90                   \
	cube.F90                         \
	cube_function.F90                \
	curv_briggs.F90                  \
	curv_gygi.F90                    \
	curvilinear.F90                  \
	curv_modine.F90                  \
	derivatives.F90                  \
	double_grid.F90                  \
	fourier_shell.F90                \
	fourier_space.F90                \
	grid.F90                         \
	hypercube.F90                    \
	index.F90                        \
	io_function.F90                  \
	mesh_batch.F90                   \
	mesh_cube_map.F90                \
	mesh_cube_parallel_map.F90       \
	mesh.F90                         \
	mesh_function.F90                \
	mesh_init.F90                    \
	mesh_interpolation.F90           \
	mesh_partition.F90               \
	modelmb_1part.F90                \
	multigrid.F90                    \
	nl_operator.F90                  \
	operate_f.F90                    \
	partition.F90                    \
	partition_transfer.F90           \
	par_vec.F90                      \
	restart.F90                      \
	sphericalmesh.F90                \
	stencil_cube.F90                 \
	stencil.F90                      \
	stencil_star.F90                 \
	stencil_stargeneral.F90          \
	stencil_starplus.F90             \
	stencil_variational.F90          \
	subarray.F90                     \
	submesh.F90                      \
	symmetrizer.F90                  \
	transfer_table.F90               \
	volume.F90                       \
	vtk.F90

# C/C++ sources - keep in alphabetical order
libgrid_c_srcs =                         \
	hilbert.c                        \
	operate.c

libgrid_a_SOURCES = $(libgrid_c_srcs) $(libgrid_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	batch_inc.F90                    \
	batch_ops_inc.F90                \
	boundaries_inc.F90               \
	cube_function_inc.F90            \
	derivatives_inc.F90              \
	double_grid_apply_inc.F90        \
	fourier_space_inc.F90            \
	io_function_inc.F90              \
	mesh_batch_inc.F90               \
	mesh_function_inc.F90            \
	mesh_interpolation_inc.F90       \
	multigrid_inc.F90                \
	nl_operator_inc.F90              \
	operate_inc.c                    \
	partition_transfer_inc.F90       \
	par_vec_inc.F90                  \
	restart_inc.F90                  \
	subarray_inc.F90                 \
	submesh_inc.F90                  \
	symmetrizer_inc.F90              \
	vtk_inc.F90

include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libgrid_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

basins.o : basins.F90 par_vec.o mesh.o index.o 
batch.o : batch.F90 batch_inc.F90 
batch_ops.o : batch_ops.F90 batch_ops_inc.F90 batch.o 
boundaries.o : boundaries.F90 boundaries_inc.F90 subarray.o par_vec.o mesh.o batch.o 
cube.o : cube.F90 index.o 
cube_function.o : cube_function.F90 cube_function_inc.F90 par_vec.o partition_transfer.o mesh_cube_parallel_map.o mesh_cube_map.o mesh.o index.o cube.o 
curv_briggs.o : curv_briggs.F90 
curv_gygi.o : curv_gygi.F90 
curvilinear.o : curvilinear.F90 curv_modine.o curv_gygi.o curv_briggs.o 
curv_modine.o : curv_modine.F90 
derivatives.o : derivatives.F90 derivatives_inc.F90 transfer_table.o stencil_variational.o stencil_stargeneral.o stencil_starplus.o stencil_star.o stencil_cube.o par_vec.o nl_operator.o mesh_function.o mesh.o boundaries.o batch.o 
double_grid.o : double_grid.F90 double_grid_apply_inc.F90 submesh.o par_vec.o mesh_function.o mesh.o curvilinear.o 
fourier_shell.o : fourier_shell.F90 mesh.o cube.o cube_function.o 
fourier_space.o : fourier_space.F90 fourier_space_inc.F90 mesh.o cube_function.o cube.o 
grid.o : grid.F90 stencil_cube.o stencil.o nl_operator.o multigrid.o mesh_init.o mesh.o index.o double_grid.o derivatives.o curvilinear.o cube.o 
hypercube.o : hypercube.F90 
index.o : index.F90 hypercube.o 
io_function.o : io_function.F90 io_function_inc.F90 vtk.o par_vec.o mesh_function.o mesh.o index.o cube.o cube_function.o 
mesh_batch.o : mesh_batch.F90 mesh_batch_inc.F90 partition.o par_vec.o mesh_function.o mesh.o index.o batch.o 
mesh_cube_map.o : mesh_cube_map.F90 index.o 
mesh_cube_parallel_map.o : mesh_cube_parallel_map.F90 par_vec.o partition_transfer.o partition.o mesh_cube_map.o mesh.o index.o cube.o 
mesh.o : mesh.F90 partition.o par_vec.o mesh_cube_map.o index.o hypercube.o curvilinear.o 
mesh_function.o : mesh_function.F90 mesh_function_inc.F90 par_vec.o modelmb_1part.o mesh.o index.o hypercube.o cube_function.o batch.o 
mesh_init.o : mesh_init.F90 subarray.o stencil.o restart.o partition.o par_vec.o mesh_partition.o mesh_cube_map.o mesh.o index.o hypercube.o curvilinear.o cube.o 
mesh_interpolation.o : mesh_interpolation.F90 mesh_interpolation_inc.F90 par_vec.o mesh.o 
mesh_partition.o : mesh_partition.F90 stencil_star.o stencil.o restart.o partition.o mesh.o index.o hypercube.o curvilinear.o 
modelmb_1part.o : modelmb_1part.F90 mesh.o hypercube.o 
multigrid.o : multigrid.F90 multigrid_inc.F90 transfer_table.o stencil.o par_vec.o mesh_init.o mesh.o derivatives.o curvilinear.o boundaries.o 
nl_operator.o : nl_operator.F90 nl_operator_inc.F90 stencil.o partition.o par_vec.o operate_f.o mesh.o index.o boundaries.o batch.o 
operate_f.o : operate_f.F90 
partition.o : partition.F90 
partition_transfer.o : partition_transfer.F90 partition_transfer_inc.F90 subarray.o par_vec.o 
par_vec.o : par_vec.F90 subarray.o stencil.o partition.o index.o 
restart.o : restart.F90 restart_inc.F90 par_vec.o mesh_batch.o mesh.o io_function.o index.o batch.o 
sphericalmesh.o : sphericalmesh.F90 mesh_interpolation.o mesh.o 
stencil_cube.o : stencil_cube.F90 stencil.o 
stencil.o : stencil.F90 
stencil_star.o : stencil_star.F90 stencil.o nl_operator.o 
stencil_stargeneral.o : stencil_stargeneral.F90 stencil.o 
stencil_starplus.o : stencil_starplus.F90 stencil.o 
stencil_variational.o : stencil_variational.F90 nl_operator.o 
subarray.o : subarray.F90 subarray_inc.F90 batch.o 
submesh.o : submesh.F90 submesh_inc.F90 par_vec.o mesh.o batch.o 
symmetrizer.o : symmetrizer.F90 symmetrizer_inc.F90 par_vec.o mesh.o 
transfer_table.o : transfer_table.F90 
volume.o : volume.F90 
vtk.o : vtk.F90 vtk_inc.F90 mesh.o index.o cube.o cube_function.o
