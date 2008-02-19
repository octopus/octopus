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
	cpmd.F90 pes.F90 spectrum.F90 \
	td_exp_split.F90 td_exp.F90 td_rti.F90 td_write.F90 td.F90 \
	td_transport.F90 td_trans_rti.F90 td_trans_intf.F90 td_trans_mem.F90

noinst_HEADERS = \
	pes_rc.F90 pes_mask.F90 td_init.F90 td_calc.F90 cpmd_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libtd_a_SOURCES)

pes.o : pes.F90 pes_mask.F90 pes_rc.F90 
spectrum.o : spectrum.F90 
td_exp_split.o : td_exp_split.F90 
td_exp.o : td_exp.F90 td_exp_split.o 
td_rti.o : td_rti.F90 td_trans_rti.o td_exp_split.o td_exp.o 
td_write.o : td_write.F90 td_calc.F90 spectrum.o 
td.o : td.F90 td_init.F90 spectrum.o pes.o td_write.o td_rti.o cpmd.o
td_transport.o : td_transport.F90 td.o td_rti.o td_trans_rti.o 
td_trans_rti.o : td_trans_rti.F90 td_trans_intf.o td_trans_mem.o td_exp.o 
td_trans_intf.o : td_trans_intf.F90 
td_trans_mem.o : td_trans_mem.F90 td_trans_intf.o
cpmd.o : cpmd.F90 cpmd_inc.F90
