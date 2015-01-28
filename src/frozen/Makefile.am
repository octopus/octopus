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
	imesh.F90		\
	igrid.F90		\
	interface_xc.F90	\
	lda.F90			\
	gga.F90			\
	functional.F90		\
	basis.F90		\
	strng.F90		\
	domain.F90		#\
	simulation.F90		\
	intrpl.F90		\
	storage.F90		\
	bgeom.F90		\
	bdnst.F90		\
	bstts.F90		\
	bsyst.F90		\
	bterm.F90		\
	bpotn.F90		\
	bepot.F90		\
	bhmlt.F90		\
	bmodl.F90		\
	bhndl.F90		\
	bcnfg.F90		\
	fio_geom.F90		\
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
	frozen_geom.F90		\
	frozen_simulation.F90	\
	frozen_density.F90	\
	frozen_states.F90	\
	frozen_system.F90	\
	frozen_external.F90	\
	frozen_hamiltonian.F90	\
	frozen_model.F90	\
	frozen_handle.F90	\
	frozen_config.F90	\
	ssys_grid.F90		\
	ssys_states.F90		\
	ssys_system.F90		\
	ssys_hamiltonian.F90	\
	ssys_model.F90		\
	ssys_handle.F90		\
	ssys_config.F90		#\
	frozen_kinetic.F90       \
	frozen_tnadd.F90         \
	frozen_ionic.F90         \
	frozen_simul_box.F90     \
	solver.F90               \
	hartree.F90              \
	frozen_parse.F90         

noinst_HEADERS =		\
	util.h			\
	template.h		\
	tsingle.F90		\
	tpair.F90		\
	tlist.F90		\
	thash.F90		\
	tdict.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libfrozen_a_SOURCES)


imesh.o : imesh.F90
igrid.o : igrid.F90
interface_xc.o : interface_xc.F90
lda.o : lda.F90 interface_xc.o
gga.o : gga.F90 interface_xc.o
functional.o : functional.F90 gga.o lda.o interface_xc.o

basis.o : basis.F90
strng.o : strng.F90
domain.o : domain.F90
# simulation.o : simulation.F90 domain.o igrid.o
# intrpl.o : intrpl.F90 simulation.o domain.o #basis.o
# storage.o : storage.F90 intrpl.o simulation.o
# bgeom.o : bgeom.F90 strng.o basis.o 
# bdnst.o : bdnst.F90 simulation.o storage.o strng.o 
# bstts.o : bstts.F90 bdnst.o simulation.o
# bsyst.o : bsyst.F90 bstts.o simulation.o bgeom.o 
# bterm.o : bterm.F90 bsyst.o
# bpotn.o : bpotn.F90 bsyst.o storage.o 
# bepot.o : bepot.F90 bpotn.o 
# bhmlt.o : bhmlt.F90 bpotn.o bterm.o bsyst.o 
# bmodl.o : bmodl.F90 bhmlt.o bsyst.o simulation.o
# bhndl.o : bhndl.F90 bmodl.o
# bcnfg.o : bcnfg.F90 bhmlt.o intrpl.o

# fio_geom.o : fio_geom.F90 bgeom.o
# fio_simul_box.o : fio_simul_box.F90 fio_geom.o
# fio_curvilinear.o : fio_curvilinear.F90 fio_simul_box.o
# fio_mesh.o : fio_mesh.F90 fio_simul_box.o imesh.o
# fio_grid.o : fio_grid.F90 fio_curvilinear.o fio_mesh.o igrid.o
# fio_simulation.o : fio_simulation.F90 fio_grid.o simulation.o
# fio_density.o : fio_density.F90 fio_simulation.o bdnst.o
# fio_states.o : fio_states.F90 fio_density.o bstts.o
# fio_system.o : fio_system.F90 fio_states.o bsyst.o
# fio_external.o : fio_external.F90 fio_simulation.o bepot.o
# fio_hamiltonian.o : fio_hamiltonian.F90 fio_external.o bhmlt.o
# fio_model.o : fio_model.F90 fio_hamiltonian.o fio_system.o fio_simulation.o bmodl.o
# fio_handle.o : fio_handle.F90 fio_model.o bhndl.o
# fio_config.o : fio_config.F90 bcnfg.o 

# frozen_geom.o : frozen_geom.F90 bgeom.o
# frozen_simulation.o : fio_simulation.F90 frozen_geom.o simulation.o
# frozen_density.o : frozen_density.F90 bdnst.o
# frozen_states.o : frozen_states.F90 frozen_density.o fio_states.o bstts.o
# frozen_system.o : frozen_system.F90 frozen_states.o fio_system.o bsyst.o
# frozen_external.o : frozen_external.F90 bepot.o #frozen_system.o fio_external.o 
# #frozen_ionic.o : frozen_ionic.F90 frozen_density.o frozen_grid.o
# #frozen_kinetic.o : frozen_kinetic.F90 frozen_system.o functional.o #frozen_grid.o
# frozen_hamiltonian.o : frozen_hamiltonian.F90 frozen_external.o frozen_system.o fio_hamiltonian.o bhmlt.o
# frozen_model.o : frozen_model.F90 frozen_hamiltonian.o frozen_system.o fio_model.o bmodl.o
# frozen_handle.o : frozen_handle.F90 frozen_model.o fio_handle.o bhndl.o
# frozen_config.o : frozen_config.F90 bcnfg.o

# #solver.o : solver.F90
# #hartree.o : hartree.F90 solver.o

# ssys_grid.o : ssys_grid.F90 igrid.o
# ssys_density.o : ssys_density.F90 bdnst.o
# ssys_states.o : ssys_states.F90 bstts.o #ssys_density.o
# ssys_system.o : ssys_system.F90 bsyst.o #ssys_states.o sigma_system.o
# ssys_hamiltonian.o : ssys_hamiltonian.F90 bhmlt.o #ssys_system.o sigma_hamiltonian.o 
# ssys_model.o : ssys_model.F90 bmodl.o #ssys_hamiltonian.o sigma_calc.o
# ssys_handle.o : ssys_handle.F90 bhndl.o
# ssys_config.o : ssys_config.F90 ssys_model.o frozen_config.o frozen_handle.o fio_config.F90 bcnfg.o functional.o

