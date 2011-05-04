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

SUBDIRS = PP recipes util opencl

dta1dir = $(pkgdatadir)
dist_dta1_DATA = logo filter_mask.data gequad.data

dta2dir = $(pkgdatadir)
dist_dta2_DATA = variables varinfo

CLEANFILES = variables varinfo varinfo_orig

varinfo: RUN_SCRIPT

variables: RUN_SCRIPT

RUN_SCRIPT:
	$(top_srcdir)/build/mk_varinfo.pl -s $(top_srcdir) -b $(top_builddir)

check:
	@for f in $(dist_dta1_DATA) $(dist_dta2_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

check-full: check

clean-local:
	if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_dta1_DATA); \
	fi

.PHONY: check check-full
