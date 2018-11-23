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

noinst_LIBRARIES = libsystem.a

# Fortran sources - keep in alphabetical order
libsystem_f_srcs =                       \
	current.F90                      \
	dos.F90                          \
	eigen_cg.F90                     \
	eigen_lobpcg.F90                 \
	eigen_rmmdiis.F90                \
	eigensolver.F90                  \
	elf.F90                          \
	energy_calc.F90                  \
	exponential.F90                  \
	forces.F90                       \
	linear_solver.F90                \
	magnetic.F90                     \
	output.F90                       \
	output_me.F90                    \
	partial_charges.F90              \
	pert.F90                         \
	preconditioners.F90              \
	stress.F90                       \
	subspace.F90                     \
	system.F90                       \
	v_ks.F90                         \
	xc_ks_inversion.F90              \
	xc_oep.F90

# C/C++ sources - keep in alphabetical order
libsystem_c_srcs =

libsystem_a_SOURCES = $(libsystem_c_srcs) $(libsystem_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	eigen_cg_inc.F90                 \
	eigen_evolution_inc.F90          \
	eigen_lobpcg_inc.F90             \
	eigen_plan_inc.F90               \
	eigen_rmmdiis_inc.F90            \
	energy_calc_inc.F90              \
	forces_inc.F90                   \
	linear_solver_inc.F90            \
	output_berkeleygw_inc.F90        \
	output_etsf_inc.F90              \
	output_h_inc.F90                 \
	output_linear_response_inc.F90   \
	output_me_inc.F90                \
	output_modelmb_inc.F90           \
	output_states_inc.F90            \
	pert_inc.F90                     \
	preconditioners_inc.F90          \
	subspace_inc.F90                 \
	xc_kli_inc.F90                   \
	xc_kli_pauli_inc.F90             \
	xc_oep_inc.F90                   \
	xc_oep_sic_inc.F90               \
	xc_oep_x_inc.F90


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libsystem_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

current.o : current.F90 
dos.o : dos.F90
eigen_cg.o : eigen_cg.F90 eigen_cg_inc.F90 preconditioners.o 
eigen_lobpcg.o : eigen_lobpcg.F90 eigen_lobpcg_inc.F90 preconditioners.o 
eigen_rmmdiis.o : eigen_rmmdiis.F90 eigen_rmmdiis_inc.F90 subspace.o preconditioners.o 
eigensolver.o : eigensolver.F90 eigen_evolution_inc.F90 eigen_plan_inc.F90 subspace.o preconditioners.o exponential.o energy_calc.o eigen_rmmdiis.o eigen_lobpcg.o eigen_cg.o 
elf.o : elf.F90 
energy_calc.o : energy_calc.F90 energy_calc_inc.F90 
exponential.o : exponential.F90 
forces.o : forces.F90 forces_inc.F90 v_ks.o
linear_solver.o : linear_solver.F90 linear_solver_inc.F90 preconditioners.o 
magnetic.o : magnetic.F90 
output.o : output.F90 output_modelmb_inc.F90 output_linear_response_inc.F90 output_h_inc.F90 output_states_inc.F90 output_etsf_inc.F90 v_ks.o output_me.o magnetic.o elf.o current.o dos.o 
output_me.o : output_me.F90 output_me_inc.F90 
partial_charges.o : partial_charges.F90 
pert.o : pert.F90 pert_inc.F90 
preconditioners.o : preconditioners.F90 preconditioners_inc.F90 
stress.o : stress.F90 v_ks.o 
subspace.o : subspace.F90 subspace_inc.F90 preconditioners.o 
system.o : system.F90 v_ks.o output.o energy_calc.o elf.o 
v_ks.o : v_ks.F90 xc_oep.o xc_ks_inversion.o magnetic.o energy_calc.o current.o 
xc_ks_inversion.o : xc_ks_inversion.F90 eigensolver.o 
xc_oep.o : xc_oep.F90 xc_oep_inc.F90 xc_oep_sic_inc.F90 xc_oep_x_inc.F90 xc_kli_inc.F90 xc_kli_pauli_inc.F90 linear_solver.o
