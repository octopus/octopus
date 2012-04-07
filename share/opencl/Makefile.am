## Process this file with automake to produce Makefile.in

## Copyright (C) 2010 X. Andrade
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
## $Id: Makefile.am 2789 2007-03-24 12:47:37Z lorenzen $

opencldir = $(pkgdatadir)/opencl

dist_opencl_DATA =              \
	vpsi.cl                 \
	set_zero.cl             \
	operate.cl              \
	axpy.cl                 \
	projector.cl            \
	pack.cl                 \
	copy.cl                 \
	subarray.cl             \
	density.cl              \
	phase.cl                \
	mesh_batch.cl           \
	rotate.cl               \
	points.cl               \
	trsm.cl                 \
	mul.cl                  \
	mesh_to_cube.cl

# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_opencl_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

check-full: check

clean-local:
	if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_opencl_DATA); \
	fi

.PHONY: check check-full
