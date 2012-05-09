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

noinst_LIBRARIES = libsystem.a

libsystem_a_SOURCES =         \
	energy.F90            \
	preconditioners.F90   \
	linear_solver.F90     \
	eigen_cg.F90          \
	eigen_lobpcg.F90      \
	eigen_rmmdiis.F90     \
	eigen.F90             \
	subspace.F90          \
	xc_oep.F90            \
	xc_ks_inversion.F90   \
	v_ks.F90              \
	output_me.F90         \
	output.F90            \
	system.F90            \
	pert.F90              \
	exponential.F90       \
	forces.F90            \
	elf.F90               \
	magnetic.F90

noinst_HEADERS =                        \
	energy_inc.F90                  \
	xc_kli_inc.F90                  \
	xc_oep_x_inc.F90                \
	xc_oep_sic_inc.F90              \
	xc_oep_inc.F90                  \
	preconditioners_inc.F90         \
	pert_inc.F90                    \
	linear_solver_inc.F90           \
	eigen_cg_inc.F90                \
	eigen_plan_inc.F90              \
	eigen_evolution_inc.F90         \
	eigen_lobpcg_inc.F90            \
	eigen_mg_inc.F90                \
	eigen_rmmdiis_inc.F90           \
	subspace_inc.F90                \
	output_h_inc.F90                \
	output_linear_response_inc.F90  \
	output_me_inc.F90               \
	output_states_inc.F90           \
	output_etsf_inc.F90             \
	output_berkeleygw_inc.F90       \
	forces_inc.F90                  \
  eigen_arpack_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libsystem_a_SOURCES)

preconditioners.o : preconditioners.F90 preconditioners_inc.F90 
subspace.o : preconditioners.o subspace_inc.F90
linear_solver.o : linear_solver.F90 linear_solver_inc.F90 preconditioners.o 
eigen_cg.o: eigen_cg.F90 eigen_cg_inc.F90 preconditioners.o
eigen_lobpcg.o: eigen_cg.F90 eigen_lobpcg_inc.F90 preconditioners.o
eigen.o: eigen.F90 eigen_evolution_inc.F90 eigen_plan_inc.F90 eigen_mg_inc.F90 eigen_arpack_inc.F90 eigen_cg.o subspace.o eigen_lobpcg.o eigen_rmmdiis.o exponential.o
eigen_rmmdiis.o: eigen_rmmdiis.F90 eigen_rmmdiis_inc.F90 subspace.o
xc_oep.o : xc_oep.F90 xc_oep_inc.F90 xc_oep_sic_inc.F90 xc_oep_x_inc.F90 xc_kli_inc.F90 linear_solver.o 
xc_ks_inversion.o : eigen.o
v_ks.o : v_ks.F90 xc_oep.o energy.o magnetic.o xc_ks_inversion.o
system.o : system.F90 v_ks.o energy.o elf.o output.o
pert.o : pert.F90 pert_inc.F90
output.o : output.F90 output_states_inc.F90 output_h_inc.F90 output_linear_response_inc.F90 output_etsf_inc.F90 output_berkeleygw_inc.F90 output_me.o elf.o magnetic.o
output_me.o: output_me_inc.F90
exponential.o : exponential.F90
energy.o : energy.F90 energy_inc.F90
forces.o : forces.F90 forces_inc.F90
elf.o : elf.F90
magnetic.o : magnetic.F90
