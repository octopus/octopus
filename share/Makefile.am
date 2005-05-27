# $Id$

SUBDIRS = PP recipes util samples

sharedir = $(pkgdatadir)
dist_share_DATA = logo variables varinfo

DISTCLEANFILES = Makefile.in

varinfo:
	(cd ..; build/mk_varinfo.pl)
