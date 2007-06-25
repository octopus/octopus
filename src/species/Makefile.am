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

noinst_LIBRARIES = libspecies.a

libspecies_a_SOURCES = \
	logrid.F90 atomic.F90 ps_in_grid.F90 ps_psf_file.F90 ps_psf.F90  \
	ps_cpi_file.F90 ps_cpi.F90 ps_fhi_file.F90 ps_fhi.F90 \
	ps_upf.F90 ps_hgh.F90 ps.F90 specie.F90 \
	xyz_file.F90 geometry.F90 xyz_adjust.F90

noinst_HEADERS =

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libspecies_a_SOURCES)
