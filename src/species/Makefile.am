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

noinst_LIBRARIES = libspecies.a

# Fortran sources - keep in alphabetical order
libspecies_f_srcs =                      \
	atomic.F90                       \
	element.F90                      \
	logrid.F90                       \
	ps_cpi.F90                       \
	ps_cpi_file.F90                  \
	ps.F90                           \
	ps_fhi.F90                       \
	ps_fhi_file.F90                  \
	ps_hgh.F90                       \
	ps_in_grid.F90                   \
	ps_psf.F90                       \
	ps_psf_file.F90                  \
	ps_xml.F90                       \
	pseudo.F90                       \
	share_directory.F90              \
	species.F90

# C/C++ sources - keep in alphabetical order
libspecies_c_srcs =                      \
	element_low.cc                   \
	pseudo_low.cc                    \
	share_directory_low.cc           \
	spline.cc

libspecies_a_SOURCES = $(libspecies_c_srcs) $(libspecies_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	anygrid.hpp                      \
	base.hpp                         \
	element.hpp                      \
	detect_format.hpp                \
	ps_pspio_inc.F90                 \
	psml.hpp                         \
	qso.hpp                          \
	share_directory.hpp              \
	spline.h                         \
	upf.hpp                          \
	upf1.hpp                         \
	upf2.hpp


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libspecies_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

atomic.o : atomic.F90 logrid.o 
element.o : element.F90 
logrid.o : logrid.F90 
ps_cpi.o : ps_cpi.F90 ps_in_grid.o ps_cpi_file.o logrid.o atomic.o 
ps_cpi_file.o : ps_cpi_file.F90 ps_in_grid.o 
ps.o : ps.F90 ps_pspio_inc.F90 ps_psf.o ps_in_grid.o ps_xml.o ps_hgh.o ps_fhi.o ps_cpi.o logrid.o atomic.o 
ps_fhi.o : ps_fhi.F90 ps_in_grid.o ps_fhi_file.o ps_cpi.o ps_cpi_file.o atomic.o 
ps_fhi_file.o : ps_fhi_file.F90 
ps_hgh.o : ps_hgh.F90 logrid.o atomic.o 
ps_in_grid.o : ps_in_grid.F90 logrid.o atomic.o 
ps_psf.o : ps_psf.F90 ps_psf_file.o ps_in_grid.o logrid.o atomic.o 
ps_psf_file.o : ps_psf_file.F90 ps_in_grid.o 
ps_xml.o : ps_xml.F90 ps_in_grid.o atomic.o pseudo.o
species.o : species.F90 ps.o logrid.o element.o share_directory.o

