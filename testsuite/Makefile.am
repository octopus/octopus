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

man_MANS = run_regression.pl.1
bin_SCRIPTS = run_regression.pl
sharedir = $(pkgdatadir)/testsuite
dist_share_DATA = finite_systems_*/*.test

rlog=regression_test.out

regression_test:
	@echo; find . -name *.test -print -exec ./run_regression.pl -f \{\}  \; | tee $(rlog)
	@echo -e "\033[32;1mPassed: `grep OK $(rlog) | wc -l` / `grep '\[ ' $(rlog) | wc -l` \033[0m"
	@echo -e "\033[31;1mFailed: `grep FAIL $(rlog) | wc -l` / `grep '\[ ' $(rlog) | wc -l` \033[0m"

clean:
	-rm -rf /tmp/tmp.* $(rlog)

run_regression.pl.1: run_regression.pl.pod
	pod2man -c "Octopus Documentation" -s 1 run_regression.pl.pod > run_regression.pl.1


.PHONY: regression clean 
