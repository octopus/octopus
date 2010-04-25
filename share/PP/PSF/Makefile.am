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
## $Id$

psfdir = $(pkgdatadir)/PP/PSF

dist_psf_DATA = \
	Ar.psf \
	As.psf \
	Au.psf \
	Cd.psf \
	C.psf  \
	Cl.psf \
	Ga.psf \
	H.psf  \
	Li.psf \
	Na.psf \
	Ne.psf \
	N.psf  \
	O.psf  \
	Se.psf \
	Si.psf \
	S.psf  \
	Ti.psf \
	Xe.psf


# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_psf_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

check-full: check

clean-local:
	@if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_psf_DATA); \
	fi

.PHONY: check check-full
