# $Id$

SUBDIRS = PP recipes util samples

sharedir = $(pkgdatadir)
dist_share_DATA = logo variables varinfo variables.xml

CLEANFILES = variables varinfo variables.xml
DISTCLEANFILES = Makefile.in

variables.xml: ../src/*.F90
	cat ../src/*.F90 | grep '^[[:space:]]*!%' | sed -e's,.*!%,,' | @PYTHON@ ../build/var2xml.py > $@

varinfo: RUN_SCRIPT

variables: RUN_SCRIPT

RUN_SCRIPT:
	(cd ..; build/mk_varinfo.pl)
