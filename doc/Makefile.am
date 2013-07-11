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
## $Id$

SUBDIRS = pdf html

man_MANS = octopus.1 oct-atomic_occupations.1 oct-casida_spectrum.1 oct-center-geom.1 oct-propagation_spectrum.1 \
        oct-display_partitions.1 oct-harmonic-spectrum.1 oct-vibrational_spectrum.1 oct-rotatory_strength.1 \
	oct-run_periodic_table.1 oct-run_regression_test.pl.1 oct-run_testsuite.sh.1 oct-xyz-anim.1

EXTRA_DIST = \
	$(man_MANS:.1=.pod) \
	update.el           \
	img/cosine.eps      \
	img/cover.ps        \
	img/gaussian.eps    \
	img/octosex.eps     \
	img/octosex.jpg     \
	img/ramp.eps        \
	ChangeLog.2001.gz   \
	ChangeLog.2002.gz   \
	ChangeLog.2003.gz   \
	ChangeLog.2004.gz   \
	ChangeLog.2005.gz   \
	ChangeLog.2006.gz   \
	ChangeLog.2007.gz   \
	ChangeLog.2008.gz

CLEANFILES = *~ *.bak variables.texi $(man_MANS)

.pod.1:
	pod2man -c "octopus documentation" -s 1 $< > $@
