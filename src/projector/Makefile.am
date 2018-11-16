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

noinst_LIBRARIES = libprojector.a

# Fortran sources - keep in alphabetical order
libprojector_f_srcs =                    \
	hgh_projector.F90                \
	kb_projector.F90                 \
	projector.F90                    \
	rkb_projector.F90 

# C/C++ sources - keep in alphabetical order
libprojector_c_srcs =                  

libprojector_a_SOURCES = $(libprojector_c_srcs) $(libprojector_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	hgh_projector_inc.F90            \
	kb_projector_inc.F90             \
	projector_inc.F90 

include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libprojector_f_srcs)

# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

hgh_projector.o : hgh_projector.F90 hgh_projector_inc.F90 
kb_projector.o : kb_projector.F90 kb_projector_inc.F90 
projector.o : projector.F90 projector_inc.F90 rkb_projector.o kb_projector.o hgh_projector.o 
rkb_projector.o : rkb_projector.F90
