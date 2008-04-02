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

noinst_LIBRARIES = libh_sys.a

libh_sys_a_SOURCES = \
	lasers.F90 specie_pot.F90 \
	rkb_projector.F90 kb_projector.F90 hgh_projector.F90 projector.F90 \
	epot.F90 h.F90 preconditioners.F90 linear_solver.F90 \
	scf_tol.F90 xc_OEP.F90 v_ks.F90 \
	output_h_sys.F90 \
	system.F90 xc_OEP_kernel.F90 \
	pert.F90                     \
	gauge_field.F90

noinst_HEADERS = \
	xc_KLI_inc.F90 xc_OEP_x_inc.F90 xc_OEP_SIC_inc.F90 xc_OEP_kernel_inc.F90 \
	kb_projector_inc.F90 hgh_projector_inc.F90 projector_inc.F90 epot_inc.F90 \
	xc_OEP_inc.F90 h_inc.F90 preconditioners_inc.F90 pert_inc.F90 \
	linear_solver_inc.F90 output_h.F90 output_linear_response.F90 output_states.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libh_sys_a_SOURCES)

lasers.o : lasers.F90 
specie_pot.o : specie_pot.F90 
rkb_projector.o : rkb_projector.F90 specie_pot.o 
kb_projector.o : kb_projector.F90 kb_projector_inc.F90 specie_pot.o 
hgh_projector.o : hgh_projector.F90 hgh_projector_inc.F90 specie_pot.o 
projector.o : projector.F90 projector_inc.F90 rkb_projector.o kb_projector.o hgh_projector.o specie_pot.o 
epot.o : epot.F90 epot_inc.F90 projector.o lasers.o specie_pot.o gauge_field.o
h.o : h.F90 h_inc.F90 lasers.o projector.o epot.o gauge_field.o
preconditioners.o : preconditioners.F90 preconditioners_inc.F90 h.o 
linear_solver.o : linear_solver.F90 linear_solver_inc.F90 preconditioners.o h.o 
scf_tol.o : scf_tol.F90 
xc_OEP.o : xc_OEP.F90 xc_OEP_inc.F90 xc_OEP_SIC_inc.F90 xc_OEP_x_inc.F90 xc_KLI_inc.F90 scf_tol.o linear_solver.o h.o 
v_ks.o : v_ks.F90 xc_OEP.o h.o 
system.o : system.F90 v_ks.o h.o 
xc_OEP_kernel.o : xc_OEP_kernel.F90 xc_OEP_kernel_inc.F90 xc_OEP.o system.o h.o 
pert.o : pert.F90 pert_inc.F90 specie_pot.o projector.o h.o epot.o 
output_h_sys.o : output_h_sys.F90 output_states.F90 output_h.F90 output_linear_response.F90 h.o
gauge_field.o : gauge_field.F90
