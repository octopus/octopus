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

noinst_LIBRARIES = libscf.a

libscf_a_SOURCES = \
	lcao.F90 \
        rdmft.F90 \
 	scf.F90 \
	ground_state.F90 \
	ob_lippmann_schwinger.F90

noinst_HEADERS = \
	lcao_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libscf_a_SOURCES)

subspace.o: subspace.F90 subspace_inc.F90
lcao.o: lcao.F90 lcao_inc.F90
rdmft.o: rdmft.F90 
scf.o: scf.F90 lcao.o ob_lippmann_schwinger.o
ground_state.o: ground_state.F90 scf.o rdmft.o lcao.o 
ob_lippmann_schwinger.o: ob_lippmann_schwinger.F90

