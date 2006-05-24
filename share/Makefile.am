# $Id$

SUBDIRS = PP recipes util samples

sharedir = $(pkgdatadir)
dist_share_DATA = logo variables varinfo

CLEANFILES = variables varinfo
DISTCLEANFILES = Makefile.in

varinfo: RUN_SCRIPT

variables: RUN_SCRIPT

RUN_SCRIPT:
	(cd ..; build/mk_varinfo.pl)
