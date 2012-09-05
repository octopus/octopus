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
	opt_control_global.F90 opt_control_state.F90 filters.F90 controlfunction.F90 opt_control_iter.F90 \
	target.F90 propagation.F90 initst.F90 opt_control.F90

noinst_HEADERS = \
	controlfunction_trans_inc.F90 check_input_inc.F90 finalcheck_inc.F90 opt_control_c_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libopt_control_a_SOURCES)


opt_control_global.o : opt_control_global.F90 
filters.o : filters.F90 
controlfunction.o : controlfunction.F90 controlfunction_trans_inc.F90 filters.o 
opt_control_iter.o : opt_control_iter.F90 controlfunction.o opt_control_global.o 
target.o : target.F90 opt_control_global.o 
initst.o: initst.F90 
propagation.o : propagation.F90 controlfunction.o target.o opt_control_global.o 
opt_control.o : opt_control.F90 finalcheck_inc.F90 opt_control_c_inc.F90 initst.o check_input_inc.F90 target.o \
  opt_control_iter.o controlfunction.o propagation.o opt_control_global.o filters.o 
