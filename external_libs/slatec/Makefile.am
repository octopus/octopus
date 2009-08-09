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

noinst_LIBRARIES = libslatec.a


libslatec_a_SOURCES = \
d1mach.f90 dcgs.f90 dchkw.f90 \
ds2y.f90 dsilus.f90 dslucs.f90 dslui2.f90 dslui.f90 dsmv.f90\
fdump.f90 i1mach.f90 isdcgs.f90 j4save.f90 qs2i1d.f90 xercnt.f90\
xerhlt.f90 xermsg.f90 xerprn.f90 xersve.f90 xgetua.f90 

depend:
	makedepf90 $(libslatec_a_SOURCES)

d1mach.o : d1mach.f90
dcgs.o : dcgs.f90
dchkw.o : dchkw.f90
ds2y.o : ds2y.f90
dsilus.o : dsilus.f90
dslucs.o : dslucs.f90
dslui2.o : dslui2.f90
dslui.o : dslui.f90
dsmv.o : dsmv.f90
fdump.o : fdump.f90
i1mach.o : i1mach.f90
isdcgs.o : isdcgs.f90
j4save.o : j4save.f90
qs2i1d.o : qs2i1d.f90
xercnt.o : xercnt.f90
xerhlt.o : xerhlt.f90
xermsg.o : xermsg.f90
xermsg.o : xerprn.f90
xersve.o : xersve.f90
xgetua.o : xgetua.f90
