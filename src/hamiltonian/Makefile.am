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

noinst_LIBRARIES = libhamiltonian.a

# Fortran sources - keep in alphabetical order
libhamiltonian_f_srcs =                  \
	berry.F90                        \
	born_charges.F90                 \
	boundary_op.F90                  \
	energy.F90                       \
	epot.F90                         \
	functionals.F90                  \
	gauge_field.F90                  \
	hamiltonian_base.F90             \
	hamiltonian.F90                  \
	hgh_projector.F90                \
	hirshfeld.F90                    \
	ion_interaction.F90              \
	kb_projector.F90                 \
	kick.F90                         \
	lasers.F90                       \
	libvdwxc.F90                     \
	oct_exchange.F90                 \
	pcm.F90                          \
	pcm_eom.F90                      \
	projector.F90                    \
	projector_matrix.F90             \
	rkb_projector.F90                \
	scdm.F90                         \
	scf_tol.F90                      \
	scissor.F90                      \
	species_pot.F90                  \
	vdw_ts.F90                       \
	xc.F90

# C/C++ sources - keep in alphabetical order
libhamiltonian_c_srcs =                  \
	vdw_ts_low.c

libhamiltonian_a_SOURCES = $(libhamiltonian_c_srcs) $(libhamiltonian_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
        functionals_list.F90             \
        fxc_inc.F90                      \
	hamiltonian_base_inc.F90         \
	hamiltonian_inc.F90              \
	hgh_projector_inc.F90            \
	kb_projector_inc.F90             \
	kxc_inc.F90                      \
	lasers_inc.F90                   \
	oct_exchange_inc.F90             \
	projector_inc.F90                \
	scdm_inc.F90                     \
	scissor_inc.F90                  \
	vxc_inc.F90

include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libhamiltonian_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

berry.o : berry.F90 
born_charges.o : born_charges.F90 
boundary_op.o : boundary_op.F90 
energy.o : energy.F90 
epot.o : epot.F90 xc.o species_pot.o projector.o lasers.o kick.o ion_interaction.o gauge_field.o 
functionals.o : functionals.F90 libvdwxc.o functionals_list.F90
gauge_field.o : gauge_field.F90 projector.o 
hamiltonian_base.o : hamiltonian_base.F90 hamiltonian_base_inc.F90 projector_matrix.o projector.o hgh_projector.o kb_projector.o epot.o 
hamiltonian.o : hamiltonian.F90 hamiltonian_inc.F90 functionals.o xc.o scissor.o scdm.o pcm.o projector.o oct_exchange.o lasers.o gauge_field.o epot.o hamiltonian_base.o energy.o boundary_op.o 
hgh_projector.o : hgh_projector.F90 hgh_projector_inc.F90 
hirshfeld.o : hirshfeld.F90 species_pot.o 
ion_interaction.o : ion_interaction.F90 
kb_projector.o : kb_projector.F90 kb_projector_inc.F90 
kick.o : kick.F90 pcm.o
lasers.o : lasers.F90 lasers_inc.F90 
libvdwxc.o : libvdwxc.F90 
oct_exchange.o : oct_exchange.F90 oct_exchange_inc.F90 xc.o 
pcm.o : pcm.F90 pcm_eom.o
pcm_eom.o : pcm_eom.F90 
projector.o : projector.F90 projector_inc.F90 rkb_projector.o kb_projector.o hgh_projector.o 
projector_matrix.o : projector_matrix.F90 
rkb_projector.o : rkb_projector.F90
scdm.o : scdm.F90 scdm_inc.F90 xc.o
scf_tol.o : scf_tol.F90
scissor.o : scissor.F90 scissor_inc.F90
species_pot.o : species_pot.F90
vdw_ts.o : vdw_ts.F90 hirshfeld.o
xc.o : xc.F90 kxc_inc.F90 fxc_inc.F90 vxc_inc.F90 functionals.o libvdwxc.o
