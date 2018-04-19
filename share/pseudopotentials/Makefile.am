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

SUBDIRS = PSF HGH quantum-simulation.org pseudo-dojo.org

ppdir = $(pkgdatadir)/pseudopotentials

dist_pp_DATA =                            \
	standard.set                      \
	sg15.set                          \
	hgh_lda.set                       \
	hscv_lda.set                      \
	hscv_pbe.set                      \
	pseudodojo_lda.set                \
	pseudodojo_lda_stringent.set      \
	pseudodojo_pbe.set                \
	pseudodojo_pbe_stringent.set      \
	pseudodojo_pbesol.set             \
	pseudodojo_pbesol_stringent.set   \
	elements.dat

# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_pp_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

clean-local:
	if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_pp_DATA); \
	fi

.PHONY: check

