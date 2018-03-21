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

noinst_LIBRARIES = libbasic.a

# Fortran sources - keep in alphabetical order
libbasic_f_srcs =                        \
	accel.F90                        \
	blacs.F90                        \
	blacs_proc_grid.F90              \
	calc_mode_par.F90                \
	command_line.F90                 \
	comm.F90                         \
	debug.F90                        \
	distributed.F90                  \
	global.F90                       \
	hardware.F90                     \
	ialist.F90                       \
	iihash.F90                       \
	io_binary_f.F90                  \
	io_csv_f.F90                     \
	io.F90                           \
	json.F90                         \
	json_parser.F90                  \
	kinds.F90                        \
	loct.F90                         \
	loct_pointer.F90                 \
	lookup.F90                       \
	messages.F90                     \
	mpi_debug.F90                    \
	mpi.F90                          \
	mpi_lib.F90                      \
	multicomm.F90                    \
	openscad.F90                     \
	parser.F90                       \
	path.F90                         \
	polyhedron.F90                   \
	profiling.F90                    \
	sort.F90                         \
	space.F90                        \
	string.F90                       \
	test_parameters.F90              \
	topology.F90                     \
	types.F90                        \
	unit.F90                         \
	unit_system.F90                  \
	utils.F90                        \
	varinfo.F90                      \
	write_iter.F90

# C/C++ sources - keep in alphabetical order
libbasic_c_srcs =                        \
	cublas.cc                        \
	cuda.cc                          \
	getopt_f.c                       \
	io_binary.c                      \
	io_csv.c                         \
	oct_f.c                          \
	oct_gdlib_f.c                    \
	recipes.c                        \
	signals.c                        \
	sort_low.cc                      \
	varia.c                          \
	varinfo_low.c                    \
	write_iter_low.cc

libbasic_a_SOURCES = $(libbasic_c_srcs) $(libbasic_f_srcs)

# Include files - keep this in alphabetical order
noinst_HEADERS =             \
	accel_inc.F90        \
	comm_inc.F90         \
	io_binary.h          \
        io_binary_f_inc.F90  \
	loct_pointer_inc.F90 \
	mpi_debug_inc.F90    \
	mpi_lib_inc.F90      \
	sort_inc.F90         \
	varia.h


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libbasic_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

accel.o : accel.F90 accel_inc.F90 unit_system.o profiling.o parser.o types.o mpi.o messages.o loct.o io.o global.o 
blacs.o : blacs.F90 
blacs_proc_grid.o : blacs_proc_grid.F90 profiling.o messages.o mpi.o blacs.o global.o 
calc_mode_par.o : calc_mode_par.F90 multicomm.o messages.o global.o 
command_line.o : command_line.F90 
comm.o : comm.F90 comm_inc.F90 profiling.o types.o mpi.o messages.o global.o 
debug.o : debug.F90 parser.o global.o 
distributed.o : distributed.F90 varinfo.o profiling.o parser.o multicomm.o mpi.o messages.o loct.o io.o global.o 
global.o : global.F90 varinfo.o mpi.o loct.o hardware.o 
hardware.o : hardware.F90 
ialist.o : ialist.F90 profiling.o messages.o global.o 
iihash.o : iihash.F90 profiling.o messages.o global.o ialist.o 
io_binary_f.o : io_binary_f.F90 io_binary_f_inc.F90 io_binary.h profiling.o mpi.o messages.o global.o 
io_csv_f.o : io_csv_f.F90 io_binary.h messages.o global.o 
io.o : io.F90 parser.o mpi.o messages.o loct.o global.o debug.o 
json.o : json.F90 profiling.o messages.o kinds.o global.o 
json_parser.o : json_parser.F90 profiling.o mpi.o messages.o kinds.o json.o io.o global.o 
kinds.o : kinds.F90 
loct.o : loct.F90 
loct_pointer.o : loct_pointer.F90 loct_pointer_inc.F90 profiling.o messages.o global.o 
lookup.o : lookup.F90 profiling.o messages.o loct_pointer.o global.o 
messages.o : messages.F90 varinfo.o unit.o string.o parser.o mpi.o loct.o debug.o global.o 
mpi_debug.o : mpi_debug.F90 mpi_debug_inc.F90 mpi.o messages.o loct.o global.o 
mpi.o : mpi.F90 blacs.o 
mpi_lib.o : mpi_lib.F90 mpi_lib_inc.F90 profiling.o mpi_debug.o mpi.o messages.o global.o 
multicomm.o : multicomm.F90 varinfo.o utils.o profiling.o parser.o mpi.o messages.o loct.o io.o global.o 
openscad.o : openscad.F90 profiling.o polyhedron.o messages.o io.o global.o 
parser.o : parser.F90 varinfo.o unit.o mpi.o loct.o global.o 
path.o : path.F90 profiling.o messages.o loct.o global.o 
polyhedron.o : polyhedron.F90 profiling.o io.o messages.o global.o 
profiling.o : profiling.F90 varinfo.o types.o string.o sort.o parser.o mpi.o messages.o loct.o io.o global.o 
sort.o : sort.F90 sort_inc.F90 messages.o global.o 
space.o : space.F90 parser.o messages.o json.o global.o 
string.o : string.F90 loct.o 
test_parameters.o : test_parameters.F90 
topology.o : topology.F90 utils.o profiling.o parser.o mpi.o messages.o loct.o io.o global.o 
types.o : types.F90 
unit.o : unit.F90 global.o 
unit_system.o : unit_system.F90 varinfo.o unit.o parser.o messages.o io.o global.o 
utils.o : utils.F90 string.o unit_system.o unit.o mpi.o messages.o loct.o io.o global.o 
varinfo.o : varinfo.F90 string.o 
write_iter.o : write_iter.F90
