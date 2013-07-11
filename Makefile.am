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

ACLOCAL_AMFLAGS = -I m4

if ENABLE_PERL
  PERL_DIR = scripts
endif

SUBDIRS =        \
	build external_libs liboct_parser \
	src share doc testsuite $(PERL_DIR)

EXTRA_DIST =     \
	ChangeLog PACKAGING \
	TODO \
	debian/changelog debian/control debian/copyright debian/dirs     \
	debian/info debian/README.Debian debian/rules

DISTCLEANFILES = \
	ChangeLog config.log config.status

DEB_DIR = deb-work

check-full:
	$(MAKE) -C testsuite check-full

rpm: dist
	@topdir=`rpmbuild --showrc | @GREP@ ': *_topdir' | cut -f 2`; \
	mv -f $(PACKAGE)-$(VERSION).tar.gz $$topdir/SOURCES; \
	cp build/$(PACKAGE).spec $$topdir/SPECS; \
	rpmbuild -ba $$topdir/SPECS/$(PACKAGE).spec; \
	echo "-------------------------------------------------------------"; \
	echo "  RPMs are built, and are located under:"; \
	echo "       $$topdir/RPMS/<arch>"; \
	echo "       $$topdir/SRPMS"; \
	echo "-------------------------------------------------------------"
#	scp $(RPM_DIR)/RPMS/*/octopus*.rpm \
		marques@nautilus.fis.uc.pt:/server/www/tddft.org/programs/octopus/download/

deb: dist
	rm -rf $(DEB_DIR)
	mkdir $(DEB_DIR)
	mv -f $(PACKAGE)-$(VERSION).tar.gz $(DEB_DIR)
	(cd $(DEB_DIR); tar xzf *gz; cd $(PACKAGE)-$(VERSION); debuild -E -ePATH)
	rm -rf $(DEB_DIR)/$(PACKAGE)-$(VERSION)
	scp $(DEB_DIR)/octopus*.deb \
		marques@nautilus.fis.uc.pt:/server/www/tddft.org/programs/octopus/download/

.PHONY: ChangeLog

ChangeLog:
	(svn2cl -i -r {20130101}:BASE -o .changelog.tmp \$(top_srcdir) && \
	mv .changelog.tmp ChangeLog) || (rm -f .changelog.tmp; touch ChangeLog; \
	echo 'svn2cl command not found: installing possibly empty changelog.' >&2)

dist-hook: ChangeLog
