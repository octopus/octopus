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

noinst_LIBRARIES = libsternheimer.a

libsternheimer_a_SOURCES =   \
	sternheimer.F90      \
	em_resp_calc.F90     \
	em_resp.F90          \
	vdw.F90              \
	phonons_lr.F90       \
	kdotp_calc.F90       \
	kdotp.F90

noinst_HEADERS =                 \
	sternheimer_inc.F90      \
	em_resp_calc_inc.F90     \
        kdotp_calc_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libsternheimer_a_SOURCES)

sternheimer.o : sternheimer.F90 sternheimer_inc.F90
em_resp_calc.o : em_resp_calc.F90 em_resp_calc_inc.F90 sternheimer.o
em_resp.o : em_resp.F90 sternheimer.o em_resp_calc.o kdotp.o kdotp_calc.o
kdotp_calc.o : kdotp_calc.F90 kdotp_calc_inc.F90 sternheimer.o
kdotp.o : kdotp.F90 kdotp_calc.o sternheimer.o
vdw.o : vdw.F90 em_resp.o sternheimer.o em_resp_calc.o
phonons_lr.o : phonons_lr.F90 sternheimer.o kdotp.o kdotp_calc.o
