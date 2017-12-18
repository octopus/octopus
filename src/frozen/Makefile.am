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
	refcount.F90            \
	message.F90		\
	envelope.F90		\
	channel.F90		\
	msgbus.F90		\
	geo_intrf.F90		\
	geo_build.F90		\
	grid_intrf.F90		\
	interface_xc.F90	\
	domain.F90		\
	simulation.F90		\
	smlt_intrf.F90          \
	storage.F90		\
	dnst.F90		\
	dnst_intrf.F90		\
	functional.F90		\
	intrpl.F90		\
	base_geometry.F90	\
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
	fio_index.F90		\
	fio_mesh.F90		\
	fio_grid.F90		\
	fio_density.F90		\
	fio_external.F90	\
	fio_hamiltonian.F90	\
	fio_model.F90		\
	fio_handle.F90		\
	fio_config.F90		\
	frozen_density.F90	\
	frozen_external.F90	\
	frozen_hamiltonian.F90	\
	frozen_model.F90	\
	frozen_handle.F90	\
	frozen_config.F90	\
	live_handle.F90		\
	ssys_ionic.F90		\
	ssys_external.F90	\
	ssys_tnadd.F90		\
	ssys_handle.F90

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
	tbase_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libfrozen_a_SOURCES)

basis.o 	: basis.F90
strng.o		: strng.F90
refcount.o      : refcount.F90
message.o	: message.F90 refcount.o
envelope.o	: envelope.F90 message.o refcount.o
channel.o	: channel.F90 envelope.o tlist_inc.F90 template.h util.h
msgbus.o	: msgbus.F90 channel.o envelope.o message.o tlist_inc.F90 template.h util.h
geo_intrf.o	: geo_intrf.F90
geo_build.o	: geo_build.F90 										\
	strng.o basis.o												\
	tdarr_inc.F90 tdict_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 			\
	template.h util.h
grid_intrf.o    : grid_intrf.F90
interface_xc.o	: interface_xc.F90
domain.o 	: domain.F90 											\
	geo_intrf.o grid_intrf.o basis.o									\
	tlist_inc.F90 template.h util.h
simulation.o	: simulation.F90 domain.o grid_intrf.o refcount.o
smlt_intrf.o    : smlt_intrf.F90 simulation.o
storage.o 	: storage.F90 simulation.o
dnst.o 		: dnst.F90 storage.o simulation.o
dnst_intrf.o	: dnst_intrf.F90 dnst.o
functional.o 	: functional.F90 storage.o simulation.o interface_xc.o
intrpl.o 	: intrpl.F90 storage.o simulation.o domain.o

base_geometry.o 	: base_geometry.F90 									\
	strng.o geo_build.o geo_intrf.o										\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90	\
	template.h util.h
base_density.o 		: base_density.F90 									\
	dnst_intrf.o dnst.o storage.o simulation.o msgbus.o strng.o						\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 	\
	template.h util.h
base_states.o 		: base_states.F90 									\
	base_density.o simulation.o strng.o									\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 	\
	template.h util.h
base_system.o 		: base_system.F90 									\
	base_states.o base_density.o base_geometry.o simulation.o strng.o					\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 	\
	template.h util.h
base_term.o 		: base_term.F90 									\
	base_system.o strng.o											\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 	\
	template.h util.h
base_potential.o	: base_potential.F90 									\
	base_system.o base_states.o base_density.o simulation.o storage.o strng.o				\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 	\
	template.h util.h
base_functional.o	: base_functional.F90 									\
	base_system.o base_states.o base_density.o simulation.o storage.o functional.o strng.o			\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 	\
	template.h util.h
base_hamiltonian.o 	: base_hamiltonian.F90 									\
	base_functional.o base_potential.o base_term.o base_system.o						\
	simulation.o storage.o strng.o										\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 						\
	tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 								\
	template.h util.h
base_model.o 		: base_model.F90 									\
	base_hamiltonian.o base_system.o base_states.o base_density.o base_geometry.o simulation.o strng.o	\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 	\
	template.h util.h
base_handle.o 		: base_handle.F90 									\
	base_model.o base_hamiltonian.o base_system.o base_states.o base_density.o base_geometry.o		\
	smlt_intrf.o simulation.o strng.o									\
	tbase_inc.F90 tdict_inc.F90 tdarr_inc.F90 thash_inc.F90 tlist_inc.F90 tpair_inc.F90 tsingle_inc.F90 	\
	template.h util.h
base_config.o 		: base_config.F90 									\
	base_handle.o base_hamiltonian.o intrpl.o 

fio_simul_box.o 	: fio_simul_box.F90
fio_curvilinear.o 	: fio_curvilinear.F90
fio_index.o 		: fio_index.F90
fio_mesh.o 		: fio_mesh.F90 fio_index.o 
fio_grid.o 		: fio_grid.F90 fio_mesh.o fio_curvilinear.o fio_simul_box.o
fio_density.o 		: fio_density.F90									\
	base_density.o intrpl.o simulation.o dnst_intrf.o dnst.o storage.o basis.o
fio_external.o 		: fio_external.F90 base_potential.o base_geometry.o intrpl.o storage.o
fio_hamiltonian.o       : fio_hamiltonian.F90 fio_external.o base_hamiltonian.o base_potential.o
fio_model.o             : fio_model.F90 									\
	fio_hamiltonian.o fio_density.o 									\
	base_model.o base_hamiltonian.o base_density.o
fio_handle.o 		: fio_handle.F90									\
	fio_model.o fio_density.o base_handle.o base_model.o base_density.o simulation.o
fio_config.o 		: fio_config.F90 fio_handle.o base_hamiltonian.o base_config.o intrpl.o

frozen_density.o 	: frozen_density.F90 fio_density.o base_density.o simulation.o
frozen_external.o 	: frozen_external.F90 fio_external.o base_potential.o simulation.o basis.o
frozen_hamiltonian.o	: frozen_hamiltonian.F90 frozen_external.o		 				\
	fio_hamiltonian.o fio_external.o base_hamiltonian.o simulation.o
frozen_model.o		: frozen_model.F90 frozen_hamiltonian.o frozen_density.o 				\
	base_model.o base_hamiltonian.o base_density.o
frozen_handle.o 	: frozen_handle.F90 frozen_model.o fio_handle.o 					\
	base_handle.o base_model.o base_geometry.o 
frozen_config.o 	: frozen_config.F90 frozen_handle.o fio_config.o fio_handle.o				\
	base_hamiltonian.o base_config.o storage.o

live_handle.o		: live_handle.F90 base_handle.o base_geometry.o

ssys_ionic.o		: ssys_ionic.F90 base_term.o base_geometry.o
ssys_external.o		: ssys_external.F90 base_potential.o storage.o
ssys_tnadd.o		: ssys_tnadd.F90 base_hamiltonian.o
ssys_handle.o		: ssys_handle.F90 live_handle.o frozen_handle.o base_handle.o

