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

noinst_LIBRARIES = libions.a

libions_a_SOURCES =         \
	xyz_file.F90        \
	geometry.F90        \
	xyz_adjust.F90      \
	simul_box.F90       \
	solids.F90          \
	vibrations.F90      \
	crystal.F90         \
	ion_dynamics.F90    \
	ion_interaction.F90 \
	symmetries.F90      \
	spglib_f.c

noinst_HEADERS =

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libions_a_SOURCES)

xyz_file.o : xyz_file.F90 
geometry.o : geometry.F90 xyz_file.o
xyz_adjust.o : xyz_adjust.F90 geometry.o 
solids.o: solids.F90 simul_box.o
simul_box.o : simul_box.F90 geometry.o
vibrations.o : vibrations.F90 simul_box.o geometry.o
crystal.o : crystal.F90 
ion_dynamics.o : ion_dynamics.F90 simul_box.o geometry.o xyz_file.o
ion_interaction.o : ion_interaction.F90 solids.o simul_box.o geometry.o
spglib.o : spglib.F90
spglib_f.o : spglib_f.c
