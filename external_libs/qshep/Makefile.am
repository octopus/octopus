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
## $Id$

noinst_LIBRARIES = libqshep.a

libqshep_a_SOURCES = qshep2d.f90 qshep3d.f90 qshepmod.f90

# how to compile a f90 file :)
.f90.o:
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $<

CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc*
