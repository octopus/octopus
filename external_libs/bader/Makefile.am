## Copyright (C) 2015 X. Andrade
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id: Makefile.am 7479 2011-02-18 11:33:08Z xavier $


noinst_LIBRARIES = libbader.a

libbader_a_SOURCES =           \
	bader_mod.f90          \
	charge_mod.f90         \
	chgcar_mod.f90         \
        command_arg.f90        \
	cube_mod.f90           \
        io_mod.f90             \
        ions_mod.f90           \
	kind_mod.f90           \
	main.f90               \
	matrix_mod.f90         \
	multipole_mod.f90      \
	options_mod.f90        \
	results_mod.f90        \
	voronoi_mod.f90        

chgcar_mod.o: matrix_mod.o charge_mod.o
io_mod.o: chgcar_mod.o cube_mod.o
bader_mod.o: kind_mod.o matrix_mod.o options_mod.o ions_mod.o charge_mod.o io_mod.o
main.o: voronoi_mod.o multipole_mod.o
voronoi_mod.o: io_mod.o
cube_mod.o: kind_mod.o charge_mod.o
results_mod.o: chgcar_mod.o bader_mod.o
multipole_mod.o: io_mod.o bader_mod.o
matrix_mod.o: kind_mod.o
charge_mod.o: matrix_mod.o
options_mod.o: kind_mod.o
ions_mod.o: kind_mod.o

CLEANFILES = *~ *.bak *.mod *.MOD *.il *.d *.pc* ifc*
