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
	tdfunction.F90 lasers.F90 specie_pot.F90 \
	rkb_projector.F90 kb_projector.F90 hgh_projector.F90 projector.F90 \
	epot.F90 h.F90 preconditioners.F90 linear_solver.F90 \
	scf_tol.F90 xc_OEP.F90 v_ks.F90 systm.F90 xc_OEP_kernel.F90

noinst_HEADERS = \
	xc_KLI_inc.F90 xc_OEP_x_inc.F90 xc_OEP_SIC_inc.F90 \
	kb_projector_inc.F90 hgh_projector_inc.F90 projector_inc.F90 epot_inc.F90 \
	xc_OEP_inc.F90 h_inc.F90 preconditioners_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libh_sys_a_SOURCES)

epot.o: epot_inc.F90
h.o: h_inc.F90
hgh_projector.o: hgh_projector_inc.F90
kb_projector.o: kb_projector_inc.F90
linear_solver.o: linear_solver_inc.F90
projector.o:  projector_inc.F90
preconditioners.o: preconditioners_inc.F90
xc_OEP.o: xc_KLI_inc.F90 xc_OEP_x_inc.F90 xc_OEP_SIC_inc.F90 xc_OEP_inc.F90
xc_OEP_kernel.o: xc_OEP_kernel_inc.F90
