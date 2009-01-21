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
## $Id: Makefile.am 2593 2006-11-20 10:58:15Z marques $

if COMPILE_NEWUOA
  noinst_LIBRARIES = libnewuoa.a
endif

##libnewuoa_a_SOURCES = bigden.f90 biglag.f90 newuob.f90 trsapp.f90 update.f90 newuoa.f90
libnewuoa_a_SOURCES = newuoa_mod.f90

# how to compile a f90 file :)
.f90.o:
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $<

CLEANFILES = *~ *.bak
