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

man_MANS = oct-run_regression_test.pl.1

sharedir = $(pkgdatadir)/testsuite
dist_share_DATA = finite_systems_*/*.test periodic_systems_*/*.test

EXTRA_DIST = oct-run_regression_test.pl oct-run_testsuite oct-run_regression_test.pl.pod

rlog=oct-regression_test.out

check: regression_test
regression_test:
	@echo; find . -name *.test -print -exec ./oct-run_regression_test.pl -e $(PWD)/../src/octopus$(EXEEXT) -f \{\}  \; | tee $(rlog)
	@echo -e "\033[32;1mPassed: `grep OK $(rlog) | wc -l` / `grep '\[ ' $(rlog) | wc -l` \033[0m"
	@if [ `grep FAIL $(rlog) | wc -l` != "0"  ] ; \
	then echo -e "\033[31;1mFailed: `grep FAIL $(rlog) | wc -l` / `grep '\[ ' $(rlog) | wc -l` \033[0m" ; \
	else echo -e "\033[32;1m\nEverything seems to be OK \n\033[0m" ; fi

clean:
	-rm -rf /tmp/tmp.* $(rlog)

oct-run_regression_test.pl.1: oct-run_regression_test.pl.pod
	pod2man -c "Octopus Documentation" -s 1 oct-run_regression_test.pl.pod > oct-run_regression_test.pl.1


.PHONY: regression_test check clean 
