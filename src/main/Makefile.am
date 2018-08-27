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

bin_PROGRAMS = octopus

octopus_SOURCES =                       \
	invert_ks.F90                   \
	memory.F90                      \
	unocc.F90                       \
	static_pol.F90                  \
	geom_opt.F90                    \
	phonons_fd.F90                  \
	casida.F90                      \
	pulpo.F90                       \
	test.F90                        \
	run.F90                         \
	main.F90

octopus_LDADD = $(all_LIBS)

noinst_HEADERS = \
	casida_inc.F90

CLEANFILES += octopus

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
	@makedepf90 $(octopus_SOURCES)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

invert_ks.o : invert_ks.F90 
memory.o : memory.F90 
unocc.o : unocc.F90 
static_pol.o : static_pol.F90 
geom_opt.o : geom_opt.F90 
phonons_fd.o : phonons_fd.F90 
casida.o : casida.F90 casida_inc.F90 
pulpo.o : pulpo.F90 
run.o : run.F90 unocc.o static_pol.o pulpo.o phonons_fd.o memory.o invert_ks.o geom_opt.o casida.o test.o
test.o : test.F90
main.o : main.F90 run.o 
