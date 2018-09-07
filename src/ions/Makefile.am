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

noinst_LIBRARIES = libions.a

# Fortran sources - keep in alphabetical order
libions_f_srcs =                         \
	atom.F90                         \
	box.F90                          \
	box_union.F90                    \
	geometry.F90                     \
	ion_dynamics.F90                 \
	kpoints.F90                      \
	periodic_copy.F90                \
	read_coords.F90                  \
	simul_box.F90                    \
	spglib_f08.f90                   \
	symmetries.F90                   \
	vibrations.F90                   \
	xyz_adjust.F90

# C/C++ sources - keep in alphabetical order
libions_c_srcs =                         \
	symmetries_finite.c

libions_a_SOURCES = $(libions_c_srcs) $(libions_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libions_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

atom.o : atom.F90 
box.o : box.F90 
box_union.o : box_union.F90 box.o 
geometry.o : geometry.F90 read_coords.o atom.o 
ion_dynamics.o : ion_dynamics.F90 simul_box.o read_coords.o geometry.o 
kpoints.o : kpoints.F90 symmetries.o geometry.o 
periodic_copy.o : periodic_copy.F90 simul_box.o geometry.o 
read_coords.o : read_coords.F90 
simul_box.o : simul_box.F90 symmetries.o kpoints.o geometry.o atom.o 
spglib_f08.o : spglib_f08.f90
symmetries.o : symmetries.F90 geometry.o spglib_f08.o
vibrations.o : vibrations.F90 simul_box.o geometry.o 
xyz_adjust.o : xyz_adjust.F90 geometry.o
