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
	smear.F90 states.F90 magnetic.F90 elf.F90 states_block.F90 \
	states_dim.F90  \
	states_calc.F90 \
	excited_states.F90 mix.F90 linear_response.F90 restart.F90 \
	density_matrix.F90 modelmb_particles.F90 ksinversion.F90

noinst_HEADERS = \
	excited_states_inc.F90 mix_inc.F90 states_kpoints.F90 \
	states_calc_inc.F90 \
	states_block_inc.F90 linear_response_inc.F90 restart_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libstates_a_SOURCES)

ksinversion.o : ksinversion.F90 modelmb_particles.o states.o
density_matrix.o : density_matrix.F90 modelmb_particles.o states.o
states.o : states.F90 states_dim.o
states_dim.o : states_dim.F90 states_kpoints.F90
states_calc.o : states_calc.F90 states_calc_inc.F90 states.o states_block.o
magnetic.o : magnetic.F90 states.o 
elf.o : elf.F90 states.o 
states_block.o : states_block.F90 states_block_inc.F90 states.o
excited_states.o : excited_states.F90 excited_states_inc.F90 states.o states_calc.o
mix.o : mix.F90 mix_inc.F90 
linear_response.o : linear_response.F90 linear_response_inc.F90 states.o states_calc.o
restart.o : restart.F90 restart_inc.F90 states.o linear_response.o states_calc.o
