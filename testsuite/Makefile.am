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

SUBDIRS = \
	finite_systems_1d finite_systems_2d finite_systems_3d \
	periodic_systems_1d 

bin_SCRIPTS = oct-run_regression_test.pl oct-run_testsuite

EXTRA_DIST = oct-run_regression_test.pl oct-run_testsuite.in

installcheck-local:
	@$(bindir)/oct-run_testsuite -c -s | tee rlog
	@echo -e "\033[32;1mPassed: `grep OK rlog | wc -l` / `grep '\[ ' rlog | wc -l` \033[0m"
	@if [ `grep FAIL rlog | wc -l` != "0"  ] ; \
	 then echo -e "\033[31;1mFailed: `grep FAIL rlog | wc -l` / `grep '\[ ' rlog | wc -l` \033[0m" ; \
	 else echo -e "\033[32;1m\nEverything seems to be OK \n\033[0m" ; fi
	@rm -f rlog

# This target works only after octopus has been installed because the data in
# $(pkgdatadir)/share is required.
check-full:
	@$(srcdir)/oct-run_testsuite -c -d $(srcdir) -l -f | tee rlog
	@echo -e "\033[32;1mPassed: `grep OK rlog | wc -l` / `grep '\[ ' rlog | wc -l` \033[0m"
	@if [ `grep FAIL rlog | wc -l` != "0"  ] ; \
	 then echo -e "\033[31;1mFailed: `grep FAIL rlog | wc -l` / `grep '\[ ' rlog | wc -l` \033[0m" ; \
	 else echo -e "\033[32;1m\nEverything seems to be OK \n\033[0m" ; fi
	@rm -f rlog

clean-local:
	-rm -rf /tmp/octopus.* 

.PHONY: check-full installcheck-local clean clean-local
