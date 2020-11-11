## Process this file with automake to produce Makefile.in

## Copyright (C) 2002-2016 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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
## $Id$

if COMPILE_OPENCL
  noinst_LTLIBRARIES = libfortrancl.la
endif

AM_CFLAGS = @CFLAGS_CLBLAS@ @CFLAGS_CLFFT@ @CL_CFLAGS@

libfortrancl_la_SOURCES =       \
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
	clblas.f90          \
	clfft.f90

if COMPILE_OPENCL
 libfortrancl_la_SOURCES +=           \
	utils.c                      \
	cl_buffer_low.c              \
	cl_command_queue_low.c       \
	cl_context_low.c             \
	cl_device_low.c              \
	cl_kernel_low.c              \
	cl_platform_low.c            \
	cl_program_low.c             \
	cl_event_low.c               \
	clblas_low.c              \
	clfft_low.c
endif

noinst_HEADERS =               \
	string_f.h             \
	localcl.h

cl_types.lo     : cl_types.f90
cl_constants.lo : cl_constants.f90
cl_platform.lo  : cl_platform.f90 cl_types.lo
cl_program.lo   : cl_program.f90 cl_types.lo
cl_kernel.lo    : cl_kernel.f90 cl_types.lo
cl_device.lo    : cl_device.f90 cl_types.lo
cl_context.lo   : cl_context.f90 cl_types.lo
cl_command_queue.lo : cl_command_queue.f90 cl_types.lo
cl_buffer.lo    : cl_buffer.f90  cl_types.lo
cl_event.lo     : cl_event.f90  cl_types.lo
cl.lo           : cl.f90 cl_constants.lo cl_types.lo cl_device.lo cl_kernel.lo cl_program.lo cl_platform.lo cl_context.lo cl_command_queue.lo cl_buffer.lo cl_event.lo
clblas.lo       : clblas.f90 cl.lo
clfft.lo        : clfft.f90 cl.lo

CLEANFILES = *~ *.bak *.mod *.MOD *.il *.d *.pc* ifc*

