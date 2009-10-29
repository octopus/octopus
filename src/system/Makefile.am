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
	xc_oep.F90            \
	v_ks.F90              \
	output_me.F90         \
	h_sys_output.F90      \
	system.F90            \
	xc_oep_kernel.F90     \
	pert.F90              \
	exponential_split.F90 \
	exponential.F90       \
	forces.F90            \
	elf.F90               \
	magnetic.F90

noinst_HEADERS =                        \
	energy_inc.F90                  \
	xc_kli_inc.F90                  \
	xc_oep_x_inc.F90                \
	xc_oep_sic_inc.F90              \
	xc_oep_kernel_inc.F90           \
	xc_oep_inc.F90                  \
	preconditioners_inc.F90         \
	pert_inc.F90                    \
	linear_solver_inc.F90           \
	output_h_inc.F90                \
	output_linear_response_inc.F90  \
	output_me_inc.F90               \
	output_states_inc.F90           \
	output_etsf_inc.F90             \
	forces_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libsystem_a_SOURCES)

preconditioners.o : preconditioners.F90 preconditioners_inc.F90 
linear_solver.o : linear_solver.F90 linear_solver_inc.F90 preconditioners.o 
xc_oep.o : xc_oep.F90 xc_oep_inc.F90 xc_oep_sic_inc.F90 xc_oep_x_inc.F90 xc_kli_inc.F90 linear_solver.o 
v_ks.o : v_ks.F90 xc_oep.o energy.o magnetic.o
system.o : system.F90 v_ks.o energy.o elf.o
xc_oep_kernel.o : xc_oep_kernel.F90 xc_oep_kernel_inc.F90 xc_oep.o system.o 
pert.o : pert.F90 pert_inc.F90
h_sys_output.o : h_sys_output.F90 output_states_inc.F90 output_h_inc.F90 output_linear_response_inc.F90 output_etsf_inc.F90 output_me.o elf.o magnetic.o
output_me.o: output_me_inc.F90
exponential_split.o : exponential_split.F90
exponential.o : exponential.F90 exponential_split.o
energy.o : energy.F90 energy_inc.F90
forces.o : forces.F90 forces_inc.F90
elf.o : elf.F90
magnetic.o : magnetic.F90
