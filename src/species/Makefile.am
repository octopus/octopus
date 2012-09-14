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

libspecies_a_SOURCES =     \
	atomic.F90         \
	logrid.F90         \
	periodic_table.F90 \
	ps.F90             \
	ps_cpi.F90         \
	ps_cpi_file.F90    \
	ps_fhi.F90         \
	ps_fhi_file.F90    \
	ps_hgh.F90         \
	ps_in_grid.F90     \
	ps_psf.F90         \
	ps_psf_file.F90    \
	ps_upf.F90         \
	species.F90

noinst_HEADERS = ps_pspio_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libspecies_a_SOURCES)

logrid.o : logrid.F90 
atomic.o : atomic.F90 logrid.o periodic_table.o
ps_in_grid.o : ps_in_grid.F90 logrid.o atomic.o 
ps_psf_file.o : ps_psf_file.F90 ps_in_grid.o 
ps_psf.o : ps_psf.F90 ps_psf_file.o ps_in_grid.o logrid.o atomic.o 
ps_cpi_file.o : ps_cpi_file.F90 ps_in_grid.o 
ps_cpi.o : ps_cpi.F90 ps_in_grid.o ps_cpi_file.o logrid.o atomic.o 
ps_fhi_file.o : ps_fhi_file.F90 
ps_fhi.o : ps_fhi.F90 ps_in_grid.o ps_fhi_file.o ps_cpi.o ps_cpi_file.o atomic.o 
ps_upf.o : ps_upf.F90 atomic.o ps_in_grid.o
ps_hgh.o : ps_hgh.F90 logrid.o atomic.o 
ps.o : ps.F90 ps_pspio_inc.F90 ps_upf.o ps_psf.o ps_in_grid.o ps_hgh.o ps_fhi.o ps_cpi.o logrid.o atomic.o 
species.o : species.F90 ps.o 
periodic_table.o : periodic_table.F90
