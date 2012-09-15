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

include $(srcdir)/../common-rules.make

bin_PROGRAMS = octopus oct-test

octopus_SOURCES = \
	invert_ks.F90 memory.F90 \
	one_shot.F90 unocc.F90 static_pol.F90 geom_opt.F90 phonons_fd.F90 \
	casida.F90 pulpo.F90 gcm.F90 run.F90 main.F90

octopus_LDADD = $(all_LIBS)

oct_test_SOURCES = test.F90
oct_test_LDADD = \
        $(top_builddir)/src/system/libsystem.a           \
        $(top_builddir)/src/hamiltonian/libhamiltonian.a \
        $(top_builddir)/src/xc/libxc.a                   \
        $(top_builddir)/src/poisson/libpoisson.a         \
        $(top_builddir)/src/states/libstates.a           \
        $(top_builddir)/src/grid/libgrid.a               \
        $(top_builddir)/src/ions/libions.a               \
        $(top_builddir)/src/species/libspecies.a         \
        $(top_builddir)/src/math/libmath.a               \
        $(top_builddir)/src/basic/libbasic.a \
        $(other_LIBS)

noinst_HEADERS = \
	casida_inc.F90

CLEANFILES += octopus oct-test octopus_mpi oct-test_mpi octopus_single oct-test_single

depend:
	@makedepf90 $(octopus_SOURCES) $(oct_test_sources)

unocc.o : unocc.F90
static_pol.o : static_pol.F90 
geom_opt.o : geom_opt.F90 
phonons_fd.o : phonons_fd.F90 
casida.o : casida.F90 casida_inc.F90 
pulpo.o : pulpo.F90 
gcm.o : gcm.F90
run.o : run.F90 unocc.o static_pol.o pulpo.o phonons_fd.o geom_opt.o casida.o memory.o gcm.o invert_ks.o
main.o : main.F90 run.o 
memory.o : memory.F90
