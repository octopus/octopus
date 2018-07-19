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

noinst_LIBRARIES = libtd.a

# Fortran sources - keep in alphabetical order
libtd_f_srcs =                           \
	pes.F90                          \
	pes_flux.F90                     \
	pes_mask.F90                     \
	pes_out.F90                      \
	pes_spm.F90                      \
	potential_interpolation.F90      \
	propagator_base.F90              \
	propagator_cn.F90                \
	propagator_etrs.F90              \
	propagator_expmid.F90            \
	propagator.F90                   \
	propagator_magnus.F90            \
	propagator_qoct.F90              \
	propagator_rk.F90                \
	spectrum.F90                     \
	td_calc.F90                      \
	td.F90                           \
	td_write.F90

# C/C++ sources - keep in alphabetical order
libtd_c_srcs =

libtd_a_SOURCES = $(libtd_c_srcs) $(libtd_f_srcs)

# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	pes_mask_out_inc.F90             \
	pes_flux_out_inc.F90   


include $(srcdir)/../common-rules.make


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).

depend:
	@makedepf90 $(libtd_f_srcs)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

pes.o : pes.F90 pes_flux.o pes_spm.o pes_mask.o 
pes_flux.o : pes_flux.F90 pes_flux_out_inc.F90 
pes_mask.o : pes_mask.F90 pes_mask_out_inc.F90 
pes_out.o : pes_out.F90 
pes_spm.o : pes_spm.F90 
potential_interpolation.o : potential_interpolation.F90 
propagator_base.o : propagator_base.F90 potential_interpolation.o 
propagator_cn.o : propagator_cn.F90 propagator_base.o potential_interpolation.o 
propagator_etrs.o : propagator_etrs.F90 propagator_base.o potential_interpolation.o 
propagator_expmid.o : propagator_expmid.F90 propagator_base.o potential_interpolation.o 
propagator.o : propagator.F90 propagator_rk.o propagator_qoct.o propagator_magnus.o propagator_expmid.o propagator_etrs.o propagator_cn.o propagator_base.o potential_interpolation.o 
propagator_magnus.o : propagator_magnus.F90 propagator_rk.o propagator_base.o potential_interpolation.o
propagator_qoct.o : propagator_qoct.F90 propagator_base.o potential_interpolation.o 
propagator_rk.o : propagator_rk.F90 propagator_base.o potential_interpolation.o 
spectrum.o : spectrum.F90 
td_calc.o : td_calc.F90 
td.o : td.F90 td_write.o propagator_base.o propagator.o potential_interpolation.o pes.o 
td_write.o : td_write.F90 td_calc.o
