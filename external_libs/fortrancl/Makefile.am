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

if COMPILE_OPENCL
  noinst_LIBRARIES = libfortrancl.a
endif

AM_CFLAGS = @CL_CFLAGS@

libfortrancl_a_SOURCES =       \
	cl_types.f90           \
	cl_buffer.f90          \
	cl_command_queue.f90   \
	cl_constants.f90       \
	cl_context.f90         \
	cl_device.f90          \
	cl_kernel.f90          \
	cl_platform.f90        \
	cl_program.f90         \
	cl_event.f90           \
	cl.f90                 \
	clamdblas.f90          \
	clamdfft.f90

if COMPILE_OPENCL
 libfortrancl_a_SOURCES +=           \
	cl_buffer_low.c              \
	cl_command_queue_low.c       \
	cl_context_low.c             \
	cl_device_low.c              \
	cl_kernel_low.c              \
	cl_platform_low.c            \
	cl_program_low.c             \
	cl_event_low.f90             \
	clamdblas_low.c              \
	clamdfft_low.c
endif

noinst_HEADERS =               \
	string_f.h             \
	localcl.h

cl_types.o     : cl_types.f90
cl_constants.o : cl_constants.f90
cl_platform.o  : cl_platform.f90 cl_types.o
cl_program.o   : cl_program.f90 cl_types.o
cl_kernel.o    : cl_kernel.f90 cl_types.o
cl_device.o    : cl_device.f90 cl_types.o
cl_context.o   : cl_context.f90 cl_types.o
cl_command_queue.o : cl_command_queue.f90 cl_types.o
cl_buffer.o    : cl_buffer.f90  cl_types.o
cl_event.o     : cl_event.f90  cl_types.o
cl.o : cl.f90 cl_constants.o cl_types.o cl_device.o cl_kernel.o cl_program.o cl_platform.o cl_context.o cl_command_queue.o cl_buffer.o
clamdblas.o    : clamdblas.f90 cl.o
clamdfft.o     : clamdfft.f90 cl.o

CLEANFILES = *~ *.bak *.mod *.MOD *.il *.d *.pc* ifc*

