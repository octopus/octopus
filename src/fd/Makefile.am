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

noinst_LIBRARIES = libfd.a

libfd_a_SOURCES =            \
	fd_interpolation.F90 \
	fd_basis.F90         \
	fd_geometry.F90      \
	fd_states.F90        \
	fd_epot.F90          \
	fd.F90

noinst_HEADERS =

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libfd_a_SOURCES)

fd_interpolation.o : fd_interpolation.F90
fd_basis.o : fd_basis.F90
fd_geometry.o : fd_geometry.F90 fd_basis.o
fd_states.o : fd_states.F90 fd_basis.o fd_interpolation.o 
fd_epot.o : fd_epot.F90 fd_basis.o fd_interpolation.o 
fd.o : fd.F90 fd_epot.o fd_states.o fd_geometry.o fd_basis.o 
