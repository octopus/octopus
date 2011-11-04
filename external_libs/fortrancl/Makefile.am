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

noinst_LIBRARIES = libfortrancl.a

libfortrancl_a_SOURCES =       \
	cl_types.F90           \
	cl_constants.F90       \
	cl_functions.F90       \
	cl_device.F90          \
	cl_platform.F90        \
	cl_context.F90         \
	cl_command_queue.F90   \
	cl.F90

if COMPILE_OPENCL
 libfortrancl_a_SOURCES +=           \
	cl_low.c                     \
	cl_device_low.c              \
	cl_platform_low.c            \
	cl_context_low.c             \
	cl_command_queue_low.c
endif

noinst_HEADERS =               \
	string_f.h

cl_types.o     : cl_types.F90 config_F90.h
cl_constants.o : cl_constants.F90 config_F90.h
cl_functions.o : cl_functions.F90 cl_types.F90 config_F90.h
cl_platform.o  : cl_platform.F90 cl_types.F90 config_F90.h
cl_device.o    : cl_device.F90 cl_types.F90 config_F90.h
cl_context.o   : cl_context.F90 cl_types.F90 config_F90.h
cl_command_queue.o : cl_command_queue.F90 cl_types.F90 config_F90.h
cl.o : cl.F90 cl_constants.F90 cl_types.o cl_functions.o cl_device.o cl_platform.o cl_context.o cl_command_queue.o

config_F90.h: $(top_builddir)/config.h
	grep '^#define' $^ > config_F90.h

CLEANFILES = *~ *.bak *.mod *.MOD *.il *.d *.pc* ifc* config_F90.h

