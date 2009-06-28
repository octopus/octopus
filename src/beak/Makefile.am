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

noinst_LIBRARIES = libbeak.a

libbeak_a_SOURCES =   \
	beak_utils.c  \
	operate_ri.c  \
	operate_bit.c

if COMPILE_VEC
libbeak_a_SOURCES += operate_ri_vec.c
endif

if COMPILE_AS
 libbeak_a_SOURCES += operate_as.S
endif

if COMPILE_BG
 libbeak_a_SOURCES += operate_bg.c
endif

noinst_HEADERS =           \
	beak.h             \
	operate_ia64.S

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libbeak_a_SOURCES)

beak_utils.o : beak_utils.c beak.h 
operate_ri.o : operate_ri.c beak.h 
operate_ri_vec.o : operate_ri_vec.c beak.h 
operate_as.S : operate_ia64.S
operate_bit.o : operate_bit.c
operate_bg.o : operate_bg.c beak.h
