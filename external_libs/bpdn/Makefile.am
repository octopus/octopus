## Copyright (C) 2011 X. Andrade
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


noinst_LIBRARIES = libbpdn.a

libbpdn_a_SOURCES =            \
	heap.c                 \
	oneProjectorCore.c     \
	spgl1_projector.c      \
	bpdn.f90               \
	expmm.c

noinst_HEADERS =               \
	heap.h                 \
	oneProjectorCore.h

spgl1_projector.o : spgl1_projector.c 
oneProjectorCore.o : oneProjectorCore.c oneProjectorCore.h heap.h
heap.o : heap.c heap.h
bpdn.o : bpdn.f90

CLEANFILES = *~ *.bak *.mod *.MOD *.il *.d *.pc* ifc*
