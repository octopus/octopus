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

noinst_LIBRARIES = libopt_control.a

libopt_control_a_SOURCES = \
	opt_control_constants.F90 filters.F90 parameters.F90 iter.F90 \
	target.F90 propagation.F90 opt_control.F90

noinst_HEADERS = \
	read.F90 defstates.F90 propagation.F90 finalcheck.F90 

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libopt_control_a_SOURCES)


opt_control_constants.o : opt_control_constants.F90 
filters.o : filters.F90 
parameters.o : parameters.F90 filters.o 
iter.o : iter.F90 parameters.o opt_control_constants.o 
target.o : target.F90 opt_control_constants.o 
propagation.o : propagation.F90 parameters.o target.o opt_control_constants.o 
opt_control.o : opt_control.F90 finalcheck.F90 defstates.F90 read.F90 target.o iter.o parameters.o \
  propagation.o opt_control_constants.o filters.o 
