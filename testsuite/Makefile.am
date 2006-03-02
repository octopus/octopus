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

bin_SCRIPTS = oct-run_regression_test.pl oct-run_testsuite

sharedir = $(pkgdatadir)/testsuite
dist_share_DATA = finite_systems_*/*.test periodic_systems_*/*.test

EXTRA_DIST = oct-run_regression_test.pl oct-run_testsuite oct-run_regression_test.pl.pod

rlog=oct-regression_test.out

check: regression_test

regression_test:
	@./oct-run_testsuite -l -s

clean:
	-rm -rf /tmp/tmp.* $(rlog)

.PHONY: regression_test check clean 
