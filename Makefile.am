## Process this file with automake to produce Makefile.in

SUBDIRS = liboct src doc

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
