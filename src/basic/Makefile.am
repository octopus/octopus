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

noinst_LIBRARIES = libbasic.a

libbasic_a_SOURCES =     \
	io_csv.c         \
	loct.F90         \
	string.F90       \
	varinfo.F90      \
	mpi.F90          \
	global.F90       \
	parser.F90       \
	datasets.F90     \
	messages.F90     \
	io.F90           \
	io_binary.c      \
	json.F90         \
	json_parser.F90  \
	mpi_lib.F90      \
	c_pointer.F90    \
	unit.F90         \
	unit_system.F90  \
	mpi_debug.F90    \
	profiling.F90    \
	utils.F90        \
	topology.F90     \
	multicomm.F90    \
	getopt_f.c       \
	command_line.F90 \
	distributed.F90  \
	hardware.F90     \
	ialist.F90       \
	iihash.F90       \
	calc_mode.F90    \
	recipes.c        \
	varia.c          \
	oct_f.c          \
	oct_gdlib_f.c    \
	write_iter.F90   \
	write_iter_low.c \
	varinfo_low.c    \
	signals.c        \
	io_binary_f.F90  \
	io_csv_f.F90     \
	lookup.F90       \
	types.F90        \
	opencl.F90       \
	space.F90        \
	blacs.F90        \
	blacs_proc_grid.F90  \
	comm.F90             \
	octcl_kernel.F90

if COMPILE_PAPI
libbasic_a_SOURCES += papi.c papi_f.F90
endif

noinst_HEADERS =               \
	mpi_debug_inc.F90      \
	mpi_lib_inc.F90        \
	loct_inc.F90           \
	varia.h                \
	io_binary.h            \
	io_csv.h               \
	comm_inc.F90           \
	opencl_inc.F90

include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libbasic_a_SOURCES)

types.o : types.F90
blacs.o : blacs.F90
mpi.o : mpi.F90 blacs.o
hardware.o : hardware.F90
c_pointer.o : c_pointer.F90
varinfo.o : varinfo.F90 string.o c_pointer.o
string.o : string.F90 c_pointer.o
loct.o : loct.F90 loct_inc.F90 c_pointer.o
global.o : global.F90 mpi.o loct.o hardware.o varinfo.o Makefile
unit.o : unit.F90 global.o
parser.o : parser.F90 mpi.o loct.o global.o c_pointer.o unit.o
datasets.o : datasets.F90 parser.o c_pointer.o
messages.o : messages.F90 varinfo.o string.o mpi.o parser.o loct.o global.o datasets.o unit.o
papi_f.o : papi_f.F90
profiling.o : profiling.F90 mpi.o messages.o loct.o io.o global.o papi_f.o

command_line.o : command_line.F90 messages.o profiling.o
io.o : io.F90 mpi.o messages.o parser.o loct.o global.o datasets.o 
unit_system.o : unit_system.F90 unit.o io.o messages.o parser.o global.o datasets.o  varinfo.o
write_iter.o : write_iter.F90 c_pointer.o
io_binary_f.o : io_binary_f.F90 io_binary.h messages.o
io_csv_f.o : io_csv_f.F90 io_csv.h messages.o
json.o : json.F90 io.o messages.o profiling.o global.o
json_parser.o : json_parser.F90 json.o io.o messages.o profiling.o global.o

opencl.o : opencl.F90 opencl_inc.F90 global.o types.o profiling.o messages.o parser.o mpi.o io.o unit_system.o
octcl_kernel.o : octcl_kernel.F90 opencl.o c_pointer.o global.o
mpi_debug.o : mpi_debug.F90 mpi_debug_inc.F90 mpi.o messages.o loct.o global.o 
mpi_lib.o : mpi_lib.F90 mpi_lib_inc.F90 mpi_debug.o mpi.o messages.o global.o profiling.o
topology.o : topology.F90 utils.o mpi.o messages.o parser.o io.o global.o datasets.o profiling.o
multicomm.o : multicomm.F90 varinfo.o utils.o mpi.o messages.o parser.o io.o global.o datasets.o profiling.o
distributed.o: distributed.F90 multicomm.o profiling.o
calc_mode.o : calc_mode.F90 multicomm.o messages.o
blacs_proc_grid.o : blacs_proc_grid.F90 blacs.o mpi.o profiling.o messages.o global.o
comm.o : comm.F90 comm_inc.F90 mpi.o messages.o global.o profiling.o

ialist.o : ialist.F90 global.o messages.o profiling.o
iihash.o : iihash.F90 ialist.o global.o
lookup.o : lookup.F90 global.o messages.o profiling.o
utils.o : utils.F90 messages.o global.o unit.o unit_system.o string.o mpi.o loct.o

space.o : space.F90 datasets.o parser.o messages.o
