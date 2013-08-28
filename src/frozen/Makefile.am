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

libfrozen_a_SOURCES =            \
	iatom.F90                \
	igrid.F90                \
	basis.F90                \
	domain.F90               \
	simulation.F90           \
	interpolation.F90        \
	storage.F90              \
	interface_xc.F90         \
	lda.F90                  \
	gga.F90                  \
	functional.F90           \
	fio_geometry.F90         \
	fio_simul_box.F90        \
	fio_mesh.F90             \
	fio_grid.F90             \
	fio_simulation.F90       #\
	fio_density.F90          #\
	fio_external.F90         \
	fio_config.F90           \
	frozen_geometry.F90      \
	frozen_grid.F90          \
	frozen_simulation.F90    \
	frozen_density.F90       \
	frozen_external.F90      \
	frozen_ionic.F90         \
	frozen_tnadd.F90         \
	frozen_config.F90        #\
	sigma_simulation.F90     \
	sigma_states.F90         \
	sigma_external.F90       \
	sigma_hartree.F90        \
	sigma_ionic.F90          \
	sigma_tnadd.F90          \
	sigma.F90                #\
	sigma_system.F90         \
	sigma_hamiltonian.F90    \
	frozen_simul_box.F90     \
	fio_wrap.F90             #\
	tgrid.F90                \
	interface_geometry.F90   \
	frozen_geometry.F90      \
	frozen_wrap.F90          \
	solver.F90               \
	hartree.F90              \
	sigma_geometry.F90       \
	sigma.F90                #\
	frozen_tnadd.F90         #\
	frozen_parse.F90         

noinst_HEADERS =		\
	template.h		\
	tpair.F90		\
	tlist.F90		\
	thash.F90		\
	tgeo.F90		\
	tbase_density.F90	\
	tstates.F90		\
	tsystem.F90		\
	tterm.F90		\
	tionic_term.F90		\
	tpotential.F90		\
	texternal_potential.F90	\
	thamiltonian.F90	\
	tcalc.F90		\
	twrap.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libfrozen_a_SOURCES)

iatom.o : iatom.F90
igrid.o : igrid.F90
basis.o : basis.F90
domain.o : domain.F90
simulation.o : simulation.F90 domain.o igrid.o
interpolation.o : interpolation.F90 simulation.o domain.o basis.o
storage.o : storage.F90 interpolation.o simulation.o

interface_xc.o : interface_xc.F90 simulation.o
lda.o : lda.F90 interface_xc.o
gga.o : gga.F90 interface_xc.o
functional.o : functional.F90 gga.o lda.o

fio_geometry.o : fio_geometry.F90
fio_simul_box.o : fio_simul_box.F90 fio_geometry.o
fio_mesh.o : fio_mesh.F90 fio_simul_box.o
fio_grid.o : fio_grid.F90 fio_mesh.o fio_simul_box.o fio_geometry.o
fio_simulation.o : fio_simulation.F90 fio_grid.o fio_geometry.o
#fio_density.o : fio_density.F90 fio_simulation.o fio_geometry.o iatom.o
#fio_external.o : fio_external.F90 fio_density.o fio_simulation.o storage.o
#fio_config.o : fio_config.F90 fio_external.o fio_density.o igeometry.o basis.o

#frozen_geometry.o : frozen_geometry.F90 igeometry.o
#frozen_simul_box.o : frozen_simul_box.F90
#frozen_grid.o : frozen_grid.F90
#frozen_simulation.o : frozen_simulation.F90 fio_config.o fio_simulation.o simulation.o #domain.o basis.o
#frozen_density.o : frozen_density.F90 frozen_geometry.o frozen_simulation.o fio_config.o iatom.o interpolation.o basis.o
#frozen_external.o : frozen_external.F90 frozen_density.o frozen_simulation.o fio_config.o fio_external.o storage.o
#frozen_hartree.o : frozen_hartree.F90 frozen_simulation.o fio_config.o storage.o 
#frozen_ionic.o : frozen_ionic.F90 frozen_density.o frozen_simulation.o fio_config.o storage.o
#frozen_tnadd.o : frozen_tnadd.F90 frozen_density.o frozen_simulation.o functional.o storage.o simulation.o
#frozen_config.o : frozen_config.F90 frozen_external.o frozen_ionic.o frozen_simulation.o fio_config.o igeometry.o basis.o

#solver.o : solver.F90
#hartree.o : hartree.F90 solver.o
#sigma_geometry.o : sigma_geometry.F90

#sigma_simulation.o : sigma_simulation.F90 frozen_config.o frozen_simulation.o simulation.o #domain.o basis.o
#sigma_states.o : sigma_states.F90 sigma_simulation.o frozen_config.o #base_states.o sub_states.o
#sigma_external.o : sigma_external.F90 sigma_simulation.o #base_external.o sub_external.o
#sigma_hartree.o : sigma_hartree.F90 sigma_simulation.o #base_hartree.o sub_hartree.o
#sigma_ionic.o : sigma_ionic.F90 sigma_simulation.o #sub_ionic.o
#sigma_tnadd.o : sigma_tnadd.F90 sigma_simulation.o #sigma_system.o sigma_states.o #tnadd.o
#sigma.o : sigma.F90 sigma_states.o sigma_simulation.o frozen_config.o #sigma_hamiltonian.o sigma_system.o base.o sub.o

