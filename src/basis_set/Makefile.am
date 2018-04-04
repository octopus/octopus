## Process this file with automake to produce Makefile.in

## Copyright (C) 2017 N. Tancogne-Dejean
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

noinst_LIBRARIES = libbasis_set.a

# Fortran sources - keep in alphabetical order
libbasis_set_a_SOURCES =             \
		atomic_orbital.F90   \
		orbitalset.F90       \
		orbitalset_utils.F90

# Include files - keep this in alphabetical order
noinst_HEADERS =                        \
		atomic_orbital_inc.F90  \
		orbitalset_inc.F90      \
		orbitalset_utils_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libbasis_set_a_SOURCES)

atomic_orbital.o : atomic_orbital.F90 atomic_orbital_inc.F90 orbitalset.o
orbitalset.o : orbitalset.F90 orbitalset_inc.F90
orbitalset_utils.o : orbitalset_utils.F90 orbitalset_utils_inc.F90 atomic_orbital.o orbitalset.o
