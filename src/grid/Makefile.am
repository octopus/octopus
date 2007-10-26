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
	simul_box.F90 curv_gygi.F90 curv_briggs.F90 curv_modine.F90 curvlinear.F90 \
	mesh_lib.F90 par_vec.F90 mesh.F90 scalar_mesh.F90 nl_operator.F90 \
	submesh.F90 stencil_star.F90 stencil_variational.F90 stencil_cube.F90 \
	stencil_starplus.F90 derivatives.F90 cube_function.F90 mesh_function.F90 \
	functions.F90 basins.F90 out.F90 multigrid.F90 double_grid.F90 grid.F90

noinst_HEADERS = \
	mesh_init.F90 par_vec_inc.F90 derivatives_inc.F90 cube_function_inc.F90 \
	mesh_function_inc.F90 functions_inc.F90 double_grid_apply.F90 \
	submesh_inc.F90 out_inc.F90 nl_operator_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libgrid_a_SOURCES)

simul_box.o : simul_box.F90 
curv_gygi.o : curv_gygi.F90 simul_box.o 
curv_briggs.o : curv_briggs.F90 simul_box.o 
curv_modine.o : curv_modine.F90 simul_box.o 
curvlinear.o : curvlinear.F90 simul_box.o curv_modine.o curv_gygi.o curv_briggs.o 
mesh_lib.o : mesh_lib.F90 
par_vec.o : par_vec.F90 par_vec_inc.F90 mesh_lib.o 
mesh.o : mesh.F90 mesh_init.F90 simul_box.o par_vec.o mesh_lib.o curvlinear.o 
scalar_mesh.o : scalar_mesh.F90 
nl_operator.o : nl_operator.F90 nl_operator_inc.F90 par_vec.o mesh.o mesh_lib.o 
submesh.o : submesh.F90 submesh_inc.F90 simul_box.o mesh.o 
stencil_star.o : stencil_star.F90 nl_operator.o 
stencil_variational.o : stencil_variational.F90 nl_operator.o 
stencil_cube.o : stencil_cube.F90 
stencil_starplus.o : stencil_starplus.F90 
derivatives.o : derivatives.F90 derivatives_inc.F90 stencil_variational.o stencil_starplus.o stencil_star.o stencil_cube.o simul_box.o nl_operator.o mesh.o 
cube_function.o : cube_function.F90 cube_function_inc.F90 simul_box.o mesh.o 
mesh_function.o : mesh_function.F90 mesh_function_inc.F90 par_vec.o mesh.o 
functions.o : functions.F90 functions_inc.F90 simul_box.o mesh_function.o mesh.o derivatives.o cube_function.o 
basins.o : basins.F90 par_vec.o mesh_lib.o mesh.o 
out.o : out.F90 out_inc.F90 par_vec.o simul_box.o mesh.o functions.o cube_function.o 
multigrid.o : multigrid.F90 mesh.o functions.o curvlinear.o 
double_grid.o : double_grid.F90 double_grid_apply.F90 submesh.o simul_box.o mesh_function.o functions.o mesh.o curvlinear.o 
grid.o : grid.F90 simul_box.o multigrid.o mesh.o functions.o double_grid.o curvlinear.o 
