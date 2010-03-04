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
	subspace.F90 \
	eigen_cg.F90 \
	eigen_lobpcg.F90 \
	eigen_rmmdiis.F90 \
        eigen.F90 \
	lcao.F90 \
	scf.F90 \
	ground_state.F90 \
	ob_lippmann_schwinger.F90

noinst_HEADERS = \
	eigen_cg_inc.F90 \
	eigen_plan_inc.F90 \
	eigen_evolution_inc.F90 \
	eigen_lobpcg_inc.F90 \
	subspace_inc.F90 \
	lcao_inc.F90 \
	eigen_mg_inc.F90 \
	eigen_rmmdiis_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libscf_a_SOURCES)

subspace.o: subspace.F90 subspace_inc.F90
eigen_cg.o: eigen_cg.F90 eigen_cg_inc.F90
eigen_lobpcg.o: eigen_cg.F90 eigen_lobpcg_inc.F90
eigen.o: eigen.F90 eigen_evolution_inc.F90 eigen_plan_inc.F90 eigen_mg_inc.F90 eigen_cg.o subspace.o eigen_lobpcg.o eigen_rmmdiis.o
lcao.o: lcao.F90 lcao_inc.F90 
scf.o: scf.F90 lcao.o eigen.o ob_lippmann_schwinger.o
ground_state.o: ground_state.F90 scf.o lcao.o 
eigen_rmmdiis.o: eigen_rmmdiis.F90 eigen_rmmdiis_inc.F90 subspace.o
ob_lippmann_schwinger.o: ob_lippmann_schwinger.F90 eigen.o

