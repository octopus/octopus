# $Id$

SUBDIRS = PP recipes util samples

sharedir = $(pkgdatadir)
dist_share_DATA = logo variables varinfo

CLEANFILES = variables varinfo
DISTCLEANFILES = Makefile.in

varinfo: variables

variables:
	(cd ..; build/mk_varinfo.pl)
