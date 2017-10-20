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

include $(srcdir)/../common-rules.make

bin_PROGRAMS = octopus oct-test octopus_mpi

octopus_SOURCES =                       \
	invert_ks.F90                   \
	memory.F90                      \
	unocc.F90                       \
	static_pol.F90                  \
	geom_opt.F90                    \
	phonons_fd.F90                  \
	casida.F90                      \
	pulpo.F90                       \
	run.F90                         \
	main.F90

octopus_LDADD = $(all_LIBS)

oct_test_SOURCES = test.F90
oct_test_LDADD = \
        $(top_builddir)/src/system/libsystem.a           \
        $(top_builddir)/src/hamiltonian/libhamiltonian.a \
        $(top_builddir)/src/basis_set/libbasis_set.a     \
        $(top_builddir)/src/poisson/libpoisson.a         \
        $(top_builddir)/src/states/libstates.a           \
        $(top_builddir)/src/frozen/libfrozen.a           \
        $(top_builddir)/src/grid/libgrid.a               \
        $(top_builddir)/src/ions/libions.a               \
        $(top_builddir)/src/species/libspecies.a         \
        $(top_builddir)/src/math/libmath.a               \
        $(top_builddir)/src/basic/libbasic.a \
        $(other_LIBS)

noinst_HEADERS = \
	casida_inc.F90

octopus_mpi_SOURCES = octopus_mpi.F90

octopus_mpi_LDADD = \
	$(top_builddir)/src/basic/libbasic.a \
	$(other_LIBS)

CLEANFILES += octopus oct-test

# debug files from MacOS gcc with -g
clean-local :
	rm -rf *.dSYM


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(octopus_SOURCES) $(oct_test_sources)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

invert_ks.o : invert_ks.F90 
memory.o : memory.F90 
unocc.o : unocc.F90 
static_pol.o : static_pol.F90 
geom_opt.o : geom_opt.F90 
phonons_fd.o : phonons_fd.F90 
casida.o : casida.F90 casida_inc.F90 
pulpo.o : pulpo.F90 
run.o : run.F90 unocc.o static_pol.o pulpo.o phonons_fd.o memory.o invert_ks.o geom_opt.o casida.o 
main.o : main.F90 run.o 
