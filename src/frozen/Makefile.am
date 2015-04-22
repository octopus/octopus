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

libfrozen_a_SOURCES =		\
	basis.F90		\
	strng.F90		\
	interface_xc.F90	\
	config_dict.F90		\
	domain.F90		\
	simulation.F90		\
	functional.F90		\
	intrpl.F90		\
	storage.F90		\
	base_geom.F90		\
	base_density.F90	\
	base_states.F90		\
	base_system.F90		\
	base_term.F90		\
	base_potential.F90	\
	base_functional.F90	\
	base_hamiltonian.F90	\
	base_model.F90		\
	base_handle.F90		\
	base_config.F90		\
	fio_simul_box.F90	\
	fio_curvilinear.F90	\
	fio_mesh.F90		\
	fio_grid.F90		\
	fio_simulation.F90	\
	fio_density.F90		\
	fio_states.F90		\
	fio_system.F90		\
	fio_external.F90	\
	fio_hamiltonian.F90	\
	fio_model.F90		\
	fio_handle.F90		\
	fio_config.F90		\
	frozen_density.F90	\
	frozen_states.F90	\
	frozen_system.F90	\
	frozen_external.F90	\
	frozen_hamiltonian.F90	\
	frozen_model.F90	\
	frozen_handle.F90	\
	frozen_config.F90	\
	live_hamiltonian.F90	\
	live_handle.F90		\
	live_config.F90		\
	ssys_density.F90	\
	ssys_states.F90		\
	ssys_system.F90		\
	ssys_ionic.F90		\
	ssys_functional.F90	\
	ssys_external.F90	\
	ssys_tnadd.F90		\
	ssys_hamiltonian.F90	\
	ssys_model.F90		\
	ssys_handle.F90		\
	ssys_config.F90

noinst_HEADERS =		\
	kfunctionals_list.F90   \
	util.h			\
	template.h		\
	tsingle.F90		\
	tdarr.F90		\
	tlist.F90		\
	thash.F90		\
	tdict.F90		\
	iterator_code.F90	\
	intrpl_inc.F90		

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libfrozen_a_SOURCES)

basis.o 	: basis.F90
strng.o 	: strng.F90
interface_xc.o	: interface_xc.F90
config_dict.o	: config_dict.F90 strng.o 
domain.o 	: domain.F90 config_dict.o
simulation.o	: simulation.F90 domain.o config_dict.o basis.o
functional.o 	: functional.F90 simulation.o interface_xc.o
intrpl.o 	: intrpl.F90 simulation.o domain.o basis.o
storage.o 	: storage.F90 intrpl.o simulation.o basis.o

base_geom.o 		: base_geom.F90 config_dict.o strng.o basis.o
base_density.o 		: base_density.F90 simulation.o storage.o config_dict.o strng.o 
base_states.o 		: base_states.F90 base_density.o simulation.o config_dict.o
base_system.o 		: base_system.F90 base_states.o base_geom.o simulation.o config_dict.o
base_term.o 		: base_term.F90 base_system.o config_dict.o
base_potential.o	: base_potential.F90 base_system.o simulation.o storage.o config_dict.o
base_functional.o	: base_functional.F90 base_system.o base_states.o base_density.o simulation.o storage.o functional.o config_dict.o
base_hamiltonian.o 	: base_hamiltonian.F90 base_functional.o base_potential.o base_term.o base_system.o simulation.o storage.o config_dict.o strng.o
base_model.o 		: base_model.F90 base_hamiltonian.o base_system.o simulation.o config_dict.o
base_handle.o 		: base_handle.F90 base_model.o intrpl.o config_dict.o
base_config.o 		: base_config.F90 base_handle.o base_hamiltonian.o intrpl.o 

fio_simul_box.o 	: fio_simul_box.F90 base_geom.o
fio_curvilinear.o 	: fio_curvilinear.F90 fio_simul_box.o base_geom.o
fio_mesh.o 		: fio_mesh.F90 fio_simul_box.o
fio_grid.o 		: fio_grid.F90 fio_curvilinear.o fio_mesh.o fio_simul_box.o
fio_simulation.o 	: fio_simulation.F90 fio_grid.o simulation.o
fio_density.o 		: fio_density.F90 fio_simulation.o base_density.o
fio_states.o 		: fio_states.F90 fio_density.o base_states.o
fio_system.o 		: fio_system.F90 fio_states.o base_system.o
fio_external.o 		: fio_external.F90 fio_system.o fio_simulation.o base_potential.o
fio_hamiltonian.o       : fio_hamiltonian.F90 fio_external.o fio_system.o base_hamiltonian.o
fio_model.o             : fio_model.F90 fio_hamiltonian.o fio_system.o fio_simulation.o base_model.o
fio_handle.o 		: fio_handle.F90 fio_model.o base_handle.o
fio_config.o 		: fio_config.F90 fio_handle.o base_hamiltonian.o base_config.o intrpl.o

frozen_density.o 	: frozen_density.F90 fio_density.o base_density.o simulation.o intrpl.o basis.o
frozen_states.o 	: frozen_states.F90 frozen_density.o fio_states.o fio_density.o base_states.o
frozen_system.o		: frozen_system.F90 frozen_states.o fio_system.o base_system.o
frozen_external.o 	: frozen_external.F90 fio_external.o base_potential.o simulation.o intrpl.o basis.o
frozen_hamiltonian.o	: frozen_hamiltonian.F90 frozen_external.o frozen_system.o fio_hamiltonian.o fio_external.o base_hamiltonian.o
frozen_model.o		: frozen_model.F90 frozen_hamiltonian.o frozen_system.o fio_model.o base_model.o
frozen_handle.o 	: frozen_handle.F90 frozen_model.o fio_handle.o base_handle.o
frozen_config.o 	: frozen_config.F90 frozen_handle.o base_hamiltonian.o base_config.o

live_hamiltonian.o	: live_hamiltonian.F90 base_hamiltonian.o simulation.o
live_handle.o		: live_handle.F90 base_handle.o
live_config.o		: live_config.F90 live_handle.o base_hamiltonian.o base_config.o

ssys_density.o		: ssys_density.F90 base_density.o
ssys_states.o		: ssys_states.F90 ssys_density.o base_states.o
ssys_system.o		: ssys_system.F90 ssys_states.o base_system.o
ssys_ionic.o		: ssys_ionic.F90 base_term.o
ssys_functional.o	: ssys_functional.F90 base_functional.o
ssys_external.o		: ssys_external.F90 base_potential.o
ssys_tnadd.o		: ssys_tnadd.F90 ssys_functional.o ssys_system.o base_hamiltonian.o
ssys_hamiltonian.o	: ssys_hamiltonian.F90 ssys_tnadd.o ssys_external.o ssys_ionic.o ssys_system.o base_hamiltonian.o simulation.o
ssys_model.o		: ssys_model.F90 ssys_hamiltonian.o ssys_system.o base_model.o
ssys_handle.o		: ssys_handle.F90 ssys_model.o live_handle.o frozen_handle.o base_handle.o
ssys_config.o		: ssys_config.F90 ssys_handle.o live_config.o frozen_config.o fio_config.F90 base_config.o functional.o
