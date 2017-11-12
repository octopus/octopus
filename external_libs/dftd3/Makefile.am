## Process this file with automake to produce Makefile.in

## Copyright (C) 2017 Xavier Andrade
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

noinst_LIBRARIES = libdftd3.a

libdftd3_a_SOURCES =     \
	api.f90          \
	common.f90       \
	core.f90         \
	pars.f90         \
	sizes.f90

api.o: api.f90 sizes.o core.o common.o
core.o: core.f90 pars.o sizes.o common.o
pars.o: pars.f90 common.o sizes.o

# how to compile a f90 file :)
.f90.o:
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $<

CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc*
