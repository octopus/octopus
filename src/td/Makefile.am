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

noinst_LIBRARIES = libtd.a

libtd_a_SOURCES = \
	cpmd.F90 pes.F90 spectrum.F90 propagator.F90 td_calc.F90 td_write.F90 td.F90 \
	ob_propagator.F90 ob_src.F90 ob_mem.F90 ob_terms.F90 tdpsf.F90 \
	pes_rc.F90

noinst_HEADERS = \
	pes_mask_inc.F90             \
	pes_mask_out_inc.F90         \
	td_init_inc.F90              \
	cpmd_inc.F90                 \
	propagator_qoct_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libtd_a_SOURCES)

pes.o: pes.F90 pes_mask_inc.F90 pes_mask_out_inc.F90 pes_rc.o tdpsf.o
pes_rc.o : pes_rc.F90
tdpsf.o: tdpsf.F90
spectrum.o: spectrum.F90
propagator.o: propagator.F90 ob_propagator.o propagator_qoct_inc.F90
ob_propagator.o: ob_propagator.F90 ob_mem.o ob_src.o ob_terms.o
ob_mem.o: ob_mem.F90 ob_terms.o
ob_src.o: ob_src.F90 ob_mem.o ob_terms.o
td_write.o: td_write.F90 spectrum.o td_calc.o
td.o: td.F90 td_init_inc.F90 spectrum.o pes.o td_write.o propagator.o cpmd.o
cpmd.o: cpmd.F90 cpmd_inc.F90
