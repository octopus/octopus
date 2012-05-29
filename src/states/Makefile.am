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

noinst_LIBRARIES = libstates.a

libstates_a_SOURCES = \
	physics_op.F90 \
	smear.F90 states.F90 \
	states_block.F90 states_dim.F90 states_calc.F90 \
	excited_states.F90 linear_response.F90 restart.F90 \
	modelmb_density_matrix.F90 modelmb_particles.F90 \
	ob_green.F90                 \
	modelmb_exchange_syms.F90    \
	density.F90                  \
	states_io.F90

noinst_HEADERS = \
	physics_op_inc.F90 \
	excited_states_inc.F90         \
	states_calc_inc.F90 \
	modelmb_density_matrix_inc.F90 \
	states_block_inc.F90 linear_response_inc.F90 restart_inc.F90 \
	states_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libstates_a_SOURCES)

ob_green.o: ob_green.F90
ksinversion.o : ksinversion.F90 modelmb_particles.o states.o
modelmb_density_matrix.o : modelmb_density_matrix.F90 modelmb_density_matrix_inc.F90 modelmb_particles.o states.o
modelmb_exchange_syms.o : modelmb_exchange_syms.F90 modelmb_density_matrix.o modelmb_particles.o states.o
physics_op.o : physics_op_inc.F90
states.o : states.F90 states_inc.F90 states_dim.o ob_green.o modelmb_particles.o smear.o
states_dim.o : states_dim.F90
states_calc.o : states_calc.F90 states_calc_inc.F90 states.o states_block.o smear.o
states_block.o : states_block.F90 states_block_inc.F90 states.o
excited_states.o : excited_states.F90 excited_states_inc.F90 states.o states_calc.o
linear_response.o : linear_response.F90 linear_response_inc.F90 states.o states_calc.o smear.o
restart.o : restart.F90 restart_inc.F90 states.o linear_response.o states_calc.o states_io.o smear.o
density.o : states.o
states_io.o : states.o
