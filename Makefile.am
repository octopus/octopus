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

ACLOCAL_AMFLAGS = -I m4

if ENABLE_PYTHON
  PYTHON_DIR = python
endif

SUBDIRS =        \
	external_libs liboct_parser liboct libxc src doc \
	share testsuite $(PYTHON_DIR)

EXTRA_DIST =     \
	TODO build/*.pl build/*.m4 build/*.in debian/*

CLEANFILES =     \
	*~ *.bak *.pc*

DISTCLEANFILES = \
	config.log config.status Makefile.in build/octopus.spec


RPM_DIR = $(HOME)/software/rpm
DEB_DIR = deb-work

check-full:
	$(MAKE) -C testsuite check-full

rpm: dist
	mv -f $(PACKAGE)-$(VERSION).tar.gz $(RPM_DIR)/SOURCES/
	rpmbuild --clean -ba build/$(PACKAGE).spec
#	scp $(RPM_DIR)/RPMS/*/octopus*.rpm \
		marques@nautilus.fis.uc.pt:/server/www/tddft.org/programs/octopus/download/
	@echo "-------------------------------------------------------------"
	@echo "  RPMs are built, and are located under:"
	@echo "       $(RPM_DIR)/RPMS/<arch>"
	@echo "       $(RPM_DIR)/SRPMS"
	@echo "-------------------------------------------------------------"

deb: dist
	rm -rf $(DEB_DIR)
	mkdir $(DEB_DIR)
	mv -f $(PACKAGE)-$(VERSION).tar.gz $(DEB_DIR)
	(cd $(DEB_DIR); tar xzf *gz; cd $(PACKAGE)-$(VERSION); debuild -E -ePATH)
	rm -rf $(DEB_DIR)/$(PACKAGE)-$(VERSION)
	scp $(DEB_DIR)/octopus*.deb \
		marques@nautilus.fis.uc.pt:/server/www/tddft.org/programs/octopus/download/
