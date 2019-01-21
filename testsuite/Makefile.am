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

SUBDIRS = \
	components       \
	finite_systems_1d \
	finite_systems_2d \
	finite_systems_3d \
	functionals      \
	real_time        \
        symmetries       \
	linear_response  \
        lda_u            \
	optimal_control  \
	periodic_systems \
	photo_electron   \
	pseudopotentials \
	tutorials

bin_SCRIPTS = oct-run_regression_test.pl oct-run_testsuite.sh

EXTRA_DIST = oct-run_regression_test.pl oct-run_testsuite.sh.in

# OCTOPUS_SHARE has to be set to an absolute path because
# oct-run_testsuite.sh changes the working directory for every
# testcase.
# If the testsuite should be skipped, e.g. in make distcheck,
# set the environment variable SKIP_CHECK to some value.
check: $(bin_SCRIPTS)
	@CHECKOPTS="" $(MAKE) check-wrapper

check-short: $(bin_SCRIPTS)
	@CHECKOPTS="-g short-run" $(MAKE) check-wrapper

check-long: $(bin_SCRIPTS)
	@CHECKOPTS="-g long-run" $(MAKE) check-wrapper

check-dryrun: $(bin_SCRIPTS)
	@CHECKOPTS="-n" $(MAKE) check-base

check-wrapper: $(bin_SCRIPTS)
	@if test "x$(SKIP_CHECK)" == "x"; then \
	    if test "x$(CHECK_REPORT)" == "x"; then \
	        CHECKREPORTOPTS="" $(MAKE) check-base; \
	    else \
	        CHECKREPORTOPTS="-r $(CHECK_REPORT)" $(MAKE) check-base; \
	    fi \
	else \
            $(MAKE) check-dryrun; \
	fi

check-base: $(bin_SCRIPTS)
	@base=`basename "$(top_builddir)/share"` && \
	dir=`dirname "$(top_builddir)/share"` && \
	sharedir="`(cd \"$$dir\" 2> /dev/null && pwd || echo \"$$dir\")`/$$base" && \
	OCTOPUS_SHARE=$$sharedir ./oct-run_testsuite.sh -c -d $(srcdir) -l $(CHECKOPTS) $(CHECKREPORTOPTS)

# Delete temporary working directories from testsuite runs, if they are older than MINTIME (in minutes)
MINTIME=600
# note: beginning line with '-' means ignore the exit status of the command
clean-local:
	-if test "x$$TEMPDIRPATH" == "x"; then TEMPDIRPATH=/tmp; fi; \
	find $$TEMPDIRPATH -name 'octopus-*' -type d -mmin +$(MINTIME) -exec rm -rf {} \; 2> /dev/null

CLEANFILES = *.log

.PHONY: check check-short check-long clean clean-local check-dryrun check-base check-wrapper
