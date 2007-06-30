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

noinst_LIBRARIES = libpoisson_isf.a

libpoisson_isf_a_SOURCES =   \
	Build_Kernel.f90     \
	fft3d.f90            \
	PSolver_Kernel.f90   \
	scaling_function.f90 \
	smooth.f90           \
        gequad.f90           \
	pfft3d.F90           \
	pconv.F90

nodist_libpoisson_isf_a_SOURCES = config_F90.h

# How to compile an f90 file :)
.f90.o:
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $<

# How to compile an F90 file :)
.F90.o:
	@@FCCPP@ @CPPFLAGS@ -I. $< > $*_oct.f90
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

scaling_function.o: lazy_100.inc  lazy_14.inc  \
	lazy_16.inc  lazy_20.inc  lazy_24.inc  \
	lazy_30.inc  lazy_40.inc  lazy_50.inc  \
	lazy_60.inc  lazy_8.inc

pconv.o: config_F90.h

pfft3d.o: config_F90.h

config_F90.h: ../../config.h
	@GREP@ '#define' $^ > $@

CLEANFILES = config_F90.h *~ *.bak *.mod *.il *.d *.pc* ifc*

noinst_HEADERS=lazy_100.inc  lazy_14.inc  \
	lazy_16.inc  lazy_20.inc  lazy_24.inc  \
	lazy_30.inc  lazy_40.inc  lazy_50.inc  \
	lazy_60.inc  lazy_8.inc
