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

noinst_LIBRARIES = libstates.a

# Fortran sources - keep in alphabetical order
libstates_f_srcs =                         \
	density.F90                      \
	excited_states.F90               \
	linear_response.F90              \
	modelmb_density_matrix.F90       \
	modelmb_exchange_syms.F90        \
	modelmb_particles.F90            \
	opt_control_state.F90            \
	physics_op.F90                   \
	smear.F90                        \
	states_block.F90                 \
	states_calc.F90                  \
	states_dim.F90                   \
	states.F90                       \
	states_group.F90                 \
	states_io.F90                    \
	states_parallel.F90              \
	states_restart.F90

# C/C++ sources - keep in alphabetical order
libstates_c_srcs =

libstates_a_SOURCES = $(libstates_c_srcs) $(libstates_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	density_inc.F90                  \
	excited_states_inc.F90           \
	linear_response_inc.F90          \
	modelmb_density_matrix_inc.F90   \
	modelmb_exchange_syms_inc.F90    \
	physics_op_inc.F90               \
	states_block_inc.F90             \
	states_calc_inc.F90              \
	states_inc.F90                   \
	states_parallel_inc.F90


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libstates_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

density.o : density.F90 density_inc.F90 states_dim.o states.o smear.o 
excited_states.o : excited_states.F90 excited_states_inc.F90 states_dim.o states_calc.o states.o 
linear_response.o : linear_response.F90 linear_response_inc.F90 states_dim.o states_calc.o states.o smear.o 
modelmb_density_matrix.o : modelmb_density_matrix.F90 modelmb_density_matrix_inc.F90 states.o modelmb_particles.o 
modelmb_exchange_syms.o : modelmb_exchange_syms.F90 modelmb_exchange_syms_inc.F90 states.o modelmb_density_matrix.o modelmb_particles.o 
modelmb_particles.o : modelmb_particles.F90 
opt_control_state.o : opt_control_state.F90 states.o 
physics_op.o : physics_op.F90 physics_op_inc.F90 
smear.o : smear.F90 
states_block.o : states_block.F90 states_block_inc.F90 states.o 
states_calc.o : states_calc.F90 states_calc_inc.F90 states_parallel.o states_dim.o states.o smear.o physics_op.o 
states_dim.o : states_dim.F90 
states.o : states.F90 states_inc.F90 states_dim.o states_group.o smear.o modelmb_particles.o 
states_group.o : states_group.F90 states_dim.o 
states_io.o : states_io.F90 states_dim.o states.o smear.o 
states_parallel.o : states_parallel.F90 states_parallel_inc.F90 states.o 
states_restart.o : states_restart.F90 states_io.o states_dim.o states_calc.o states.o smear.o linear_response.o
