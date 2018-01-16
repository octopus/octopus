#aclocaldir = $(datadir)/aclocal
#dist_aclocal_DATA = config/m4/etsfio.m4 \
#	config/m4/libxc.m4 \
#	config/m4/psolver.m4

EXTRA_DIST= config/m4

ACLOCAL_AMFLAGS = -I config/m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = psolver.pc

SUBDIRS = src

if HAVE_PYYAML
SUBDIRS += tests
endif

SUBDIRS += doc

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(pythondir)/report.py
