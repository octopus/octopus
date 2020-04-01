## Process this file with automake to produce Makefile.in

## Copyright (C) 2020 S. Ohlmann
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

SUBDIRS = tests

sharedir = $(pkgdatadir)/testsuite/performance

dist_share_DATA = compare_results.py Makefile.performance

bin_SCRIPTS = oct-run_performance_testsuite.sh create_combinations.py

EXTRA_DIST = oct-run_performance_testsuite.sh.in create_combinations.py

# OCTOPUS_SHARE has to be set to an absolute path because
# oct-run_testsuite.sh changes the working directory for every
# testcase.
# If the testsuite should be skipped, e.g. in make distcheck,
# set the environment variable SKIP_CHECK to some value.
check-performance: $(bin_SCRIPTS)
	@base=`basename "$(top_builddir)/share"` && \
	dir=`dirname "$(top_builddir)/share"` && \
	sharedir="`(cd \"$$dir\" 2> /dev/null && pwd || echo \"$$dir\")`/$$base" && \
	OCTOPUS_SHARE=$$sharedir ./oct-run_performance_testsuite.sh -l -d $(srcdir)


.PHONY: check-performance
