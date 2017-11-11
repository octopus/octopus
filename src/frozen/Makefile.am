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

noinst_LIBRARIES = libfrozen.a

# Fortran sources - keep in alphabetical order
libfrozen_f_srcs =                       \
	base_config.F90                  \
	base_density.F90                 \
	base_functional.F90              \
	base_geometry.F90                \
	base_hamiltonian.F90             \
	base_handle.F90                  \
	base_model.F90                   \
	base_potential.F90               \
	base_states.F90                  \
	base_system.F90                  \
	base_term.F90                    \
	basis.F90                        \
	config_dict.F90                  \
	domain.F90                       \
	fio_config.F90                   \
	fio_curvilinear.F90              \
	fio_density.F90                  \
	fio_external.F90                 \
	fio_geometry.F90                 \
	fio_grid.F90                     \
	fio_hamiltonian.F90              \
	fio_handle.F90                   \
	fio_index.F90                    \
	fio_mesh.F90                     \
	fio_model.F90                    \
	fio_simulation.F90               \
	fio_simul_box.F90                \
	fio_states.F90                   \
	fio_system.F90                   \
	frozen_config.F90                \
	frozen_density.F90               \
	frozen_external.F90              \
	frozen_hamiltonian.F90           \
	frozen_handle.F90                \
	frozen_model.F90                 \
	frozen_states.F90                \
	frozen_system.F90                \
	functional.F90                   \
	geo_build.F90                    \
	geo_intrf.F90                    \
	grid_intrf.F90                   \
	interface_xc.F90                 \
	intrpl.F90                       \
	live_geometry.F90                \
	live_handle.F90                  \
	live_model.F90                   \
	live_system.F90                  \
	simulation.F90                   \
	ssys_config.F90                  \
	ssys_density.F90                 \
	ssys_external.F90                \
	ssys_hamiltonian.F90             \
	ssys_handle.F90                  \
	ssys_ionic.F90                   \
	ssys_model.F90                   \
	ssys_states.F90                  \
	ssys_system.F90                  \
	ssys_tnadd.F90                   \
	storage.F90                      \
	strng.F90

# C/C++ sources - keep in alphabetical order
libfrozen_c_srcs =

libfrozen_a_SOURCES = $(libfrozen_c_srcs) $(libfrozen_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =		\
	kfunctionals_list.F90   \
	util.h			\
	template.h		\
	tsingle_inc.F90		\
	tpair_inc.F90		\
	tdarr_inc.F90		\
	tlist_inc.F90		\
	thash_inc.F90		\
	tdict_inc.F90		\
	iterator_inc.F90

include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libfrozen_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

base_config.o : base_config.F90 storage.o simulation.o intrpl.o base_handle.o base_hamiltonian.o 
base_density.o : base_density.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 storage.o simulation.o config_dict.o 
base_functional.o : base_functional.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 storage.o simulation.o functional.o config_dict.o base_system.o base_states.o base_density.o 
base_geometry.o : base_geometry.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 geo_intrf.o geo_build.o config_dict.o 
base_hamiltonian.o : base_hamiltonian.F90 iterator_inc.F90 tdict_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 strng.o storage.o simulation.o config_dict.o base_term.o base_system.o base_potential.o base_functional.o 
base_handle.o : base_handle.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 simulation.o config_dict.o base_model.o 
base_model.o : base_model.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 simulation.o config_dict.o base_system.o base_hamiltonian.o base_geometry.o 
base_potential.o : base_potential.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 storage.o simulation.o config_dict.o base_system.o base_states.o base_density.o 
base_states.o : base_states.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 simulation.o config_dict.o base_density.o 
base_system.o : base_system.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 simulation.o config_dict.o base_states.o base_geometry.o base_density.o 
base_term.o : base_term.F90 iterator_inc.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 template.h util.h tlist_inc.F90 config_dict.o base_system.o 
basis.o : basis.F90 
config_dict.o : config_dict.F90 tsingle_inc.F90 tdarr_inc.F90 tlist_inc.F90 tpair_inc.F90 thash_inc.F90 util.h template.h tdict_inc.F90 strng.o 
domain.o : domain.F90 template.h util.h tlist_inc.F90 grid_intrf.o basis.o 
fio_config.o : fio_config.F90 intrpl.o fio_handle.o base_hamiltonian.o base_config.o 
fio_curvilinear.o : fio_curvilinear.F90 
fio_density.o : fio_density.F90 storage.o simulation.o intrpl.o base_density.o 
fio_external.o : fio_external.F90 storage.o simulation.o intrpl.o base_system.o base_potential.o base_geometry.o 
fio_geometry.o : fio_geometry.F90 geo_intrf.o base_geometry.o 
fio_grid.o : fio_grid.F90 fio_simul_box.o fio_mesh.o fio_curvilinear.o 
fio_hamiltonian.o : fio_hamiltonian.F90 fio_external.o base_potential.o base_hamiltonian.o 
fio_handle.o : fio_handle.F90 simulation.o fio_simulation.o fio_model.o base_model.o base_handle.o 
fio_index.o : fio_index.F90 
fio_mesh.o : fio_mesh.F90 fio_index.o 
fio_model.o : fio_model.F90 fio_system.o fio_hamiltonian.o base_system.o base_model.o base_hamiltonian.o 
fio_simulation.o : fio_simulation.F90 simulation.o grid_intrf.o fio_grid.o 
fio_simul_box.o : fio_simul_box.F90 
fio_states.o : fio_states.F90 fio_density.o base_states.o base_density.o 
fio_system.o : fio_system.F90 fio_states.o fio_geometry.o base_system.o base_states.o base_geometry.o 
frozen_config.o : frozen_config.F90 storage.o frozen_handle.o base_hamiltonian.o base_config.o 
frozen_density.o : frozen_density.F90 simulation.o fio_density.o basis.o base_density.o 
frozen_external.o : frozen_external.F90 simulation.o fio_external.o basis.o base_potential.o 
frozen_hamiltonian.o : frozen_hamiltonian.F90 simulation.o frozen_external.o base_system.o base_potential.o base_hamiltonian.o 
frozen_handle.o : frozen_handle.F90 frozen_model.o fio_handle.o base_model.o base_handle.o 
frozen_model.o : frozen_model.F90 frozen_system.o frozen_hamiltonian.o base_system.o base_model.o base_hamiltonian.o 
frozen_states.o : frozen_states.F90 frozen_density.o base_states.o base_density.o 
frozen_system.o : frozen_system.F90 frozen_states.o base_system.o base_states.o 
functional.o : functional.F90 storage.o simulation.o interface_xc.o 
geo_build.o : geo_build.F90 tsingle_inc.F90 tdarr_inc.F90 tpair_inc.F90 thash_inc.F90 tdict_inc.F90 template.h util.h tlist_inc.F90 basis.o strng.o 
geo_intrf.o : geo_intrf.F90 
grid_intrf.o : grid_intrf.F90 
interface_xc.o : interface_xc.F90 
intrpl.o : intrpl.F90 storage.o simulation.o domain.o 
live_geometry.o : live_geometry.F90 geo_intrf.o base_geometry.o 
live_handle.o : live_handle.F90 live_model.o base_model.o base_handle.o 
live_model.o : live_model.F90 live_system.o base_system.o base_model.o 
live_system.o : live_system.F90 live_geometry.o base_system.o base_geometry.o 
simulation.o : simulation.F90 grid_intrf.o domain.o 
ssys_config.o : ssys_config.F90 storage.o ssys_handle.o functional.o frozen_handle.o frozen_config.o fio_handle.o fio_config.o base_handle.o base_hamiltonian.o base_config.o 
ssys_density.o : ssys_density.F90 base_density.o 
ssys_external.o : ssys_external.F90 storage.o base_potential.o base_density.o 
ssys_hamiltonian.o : ssys_hamiltonian.F90 base_hamiltonian.o 
ssys_handle.o : ssys_handle.F90 live_handle.o frozen_handle.o base_handle.o 
ssys_ionic.o : ssys_ionic.F90 base_term.o base_system.o base_geometry.o 
ssys_model.o : ssys_model.F90 base_model.o 
ssys_states.o : ssys_states.F90 base_states.o 
ssys_system.o : ssys_system.F90 base_system.o 
ssys_tnadd.o : ssys_tnadd.F90 base_hamiltonian.o base_functional.o 
storage.o : storage.F90 simulation.o 
strng.o : strng.F90
