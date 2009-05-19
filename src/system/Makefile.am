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
	xc_OEP.F90            \
	v_ks.F90              \
	output_me.F90         \
	output_h_sys.F90      \
	system.F90            \
	xc_OEP_kernel.F90     \
	pert.F90              \
	exponential_split.F90 \
	exponential.F90

noinst_HEADERS =                     \
	energy_inc.F90               \
	xc_KLI_inc.F90               \
	xc_OEP_x_inc.F90             \
	xc_OEP_SIC_inc.F90           \
	xc_OEP_kernel_inc.F90        \
	xc_OEP_inc.F90               \
	preconditioners_inc.F90      \
	pert_inc.F90                 \
	linear_solver_inc.F90        \
	output_h.F90                 \
	output_linear_response.F90   \
	output_me_inc.F90            \
	output_states.F90            \
	output_etsf.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libsystem_a_SOURCES)

preconditioners.o : preconditioners.F90 preconditioners_inc.F90 
linear_solver.o : linear_solver.F90 linear_solver_inc.F90 preconditioners.o 
xc_OEP.o : xc_OEP.F90 xc_OEP_inc.F90 xc_OEP_SIC_inc.F90 xc_OEP_x_inc.F90 xc_KLI_inc.F90 linear_solver.o 
v_ks.o : v_ks.F90 xc_OEP.o energy.o
system.o : system.F90 v_ks.o energy.o
xc_OEP_kernel.o : xc_OEP_kernel.F90 xc_OEP_kernel_inc.F90 xc_OEP.o system.o 
pert.o : pert.F90 pert_inc.F90
output_h_sys.o : output_h_sys.F90 output_states.F90 output_h.F90 output_linear_response.F90 output_etsf.F90 output_me.o
output_me.o: output_me_inc.F90
exponential_split.o : exponential_split.F90
exponential.o : exponential.F90 exponential_split.o
energy.o : energy.F90 energy_inc.F90
