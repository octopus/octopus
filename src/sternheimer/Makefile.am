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

noinst_LIBRARIES = libsternheimer.a

# Fortran sources - keep in alphabetical order
libsternheimer_f_srcs =                  \
	em_resp_calc.F90                 \
	em_resp.F90                      \
	kdotp_calc.F90                   \
	kdotp.F90                        \
	phonons_lr.F90                   \
	sternheimer.F90                  \
	vdw.F90

# C/C++ sources - keep in alphabetical order
libsternheimer_c_srcs =

libsternheimer_a_SOURCES = $(libsternheimer_c_srcs) $(libsternheimer_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	em_resp_calc_inc.F90             \
	em_resp_inc.F90                  \
	kdotp_calc_inc.F90               \
	phonons_lr_inc.F90               \
	sternheimer_inc.F90 


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libsternheimer_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

em_resp_calc.o : em_resp_calc.F90 em_resp_calc_inc.F90 sternheimer.o 
em_resp.o : em_resp.F90 em_resp_inc.F90 sternheimer.o kdotp_calc.o kdotp.o em_resp_calc.o 
kdotp_calc.o : kdotp_calc.F90 kdotp_calc_inc.F90 sternheimer.o 
kdotp.o : kdotp.F90 sternheimer.o kdotp_calc.o 
phonons_lr.o : phonons_lr.F90 phonons_lr_inc.F90 sternheimer.o kdotp_calc.o kdotp.o 
sternheimer.o : sternheimer.F90 sternheimer_inc.F90 
vdw.o : vdw.F90 sternheimer.o em_resp_calc.o em_resp.o
