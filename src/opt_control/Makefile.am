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

noinst_LIBRARIES = libopt_control.a

# Fortran sources - keep in alphabetical order
libopt_control_f_srcs =                  \
	controlfunction.F90              \
	filter.F90                       \
	initst.F90                       \
	opt_control.F90                  \
	opt_control_global.F90           \
	opt_control_iter.F90             \
	propagation.F90                  \
	target.F90

# C/C++ sources - keep in alphabetical order
libopt_control_c_srcs =

libopt_control_a_SOURCES = $(libopt_control_c_srcs) $(libopt_control_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	check_input_inc.F90              \
	controlfunction_trans_inc.F90    \
	finalcheck_inc.F90               \
	opt_control_c_inc.F90            \
	target_classical_inc.F90         \
	target_density_inc.F90           \
	target_excited_inc.F90           \
	target_exclude_inc.F90           \
	target_groundstate_inc.F90       \
	target_gstransformation_inc.F90  \
	target_hhg_inc.F90               \
	target_local_inc.F90             \
	target_spin_inc.F90              \
	target_tdlocal_inc.F90           \
	target_userdefined_inc.F90       \
	target_velocity_inc.F90 


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libopt_control_a_SOURCES)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

controlfunction.o : controlfunction.F90 controlfunction_trans_inc.F90 filter.o 
filter.o : filter.F90 
initst.o : initst.F90 
opt_control.o : opt_control.F90 finalcheck_inc.F90 check_input_inc.F90 opt_control_c_inc.F90 target.o propagation.o opt_control_iter.o opt_control_global.o initst.o filter.o controlfunction.o 
opt_control_global.o : opt_control_global.F90 
opt_control_iter.o : opt_control_iter.F90 controlfunction.o 
propagation.o : propagation.F90 target.o controlfunction.o 
target.o : target.F90 target_spin_inc.F90 target_classical_inc.F90 target_tdlocal_inc.F90 target_local_inc.F90 target_userdefined_inc.F90 target_exclude_inc.F90 target_gstransformation_inc.F90 target_excited_inc.F90 target_groundstate_inc.F90 target_hhg_inc.F90 target_velocity_inc.F90 target_density_inc.F90 opt_control_global.o
