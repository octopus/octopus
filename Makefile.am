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

SUBDIRS = liboct src doc share

EXTRA_DIST = TODO autogen.sh build/*.m4 build/*.in samples/inp

CLEANFILES = *~ *.bak *.pc*
DISTCLEANFILES = Makefile.in build/octopus.spec

RPM_SRCDIR = `rpm --showrc |grep "^sourcedir" | sed 's/^.*: *//'`
RPM_ARCH   = `rpm --showrc |grep "^build arch" | sed 's/^.*: *//'`
RPM_RPMS   = `rpm --showrc |grep "^rpmdir" | sed 's/^.*: *//'`
RPM_SRPMS  = `rpm --showrc |grep "^srcrpmdir" | sed 's/^.*: *//'`
RPM_TOPDIR = `rpm --showrc |grep "^TOPDIRdir" | sed 's/^.*: *//'`

rpm: dist
	cp -f $(PACKAGE)-$(VERSION).tar.gz $(RPM_SRCDIR)/$(PACKAGE)-$(VERSION).tar.gz
#	rm -f $(RPM_SRCDIR)/$(PACKAGE)-logo-thumb.gif
#	lynx -source -dump http://www.fftw.org/fftw-logo-thumb.gif > $(RPM_SRCDIR)/fftw-logo-thumb.gif
	rpm --clean -ba build/$(PACKAGE).spec
	rpm -U $(RPM_RPMS)/$(RPM_ARCH)/$(PACKAGE)-$(VERSION)-1.$(RPM_ARCH).rpm
	@echo "-------------------------------------------------------------"
	@echo "  RPMs are built and installed, and are located under:       "
	@echo "       $(RPM_RPMS)/$(RPM_ARCH)"
	@echo "       $(RPM_SRPMS)"
	@echo "-------------------------------------------------------------"
