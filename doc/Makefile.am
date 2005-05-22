# $Id$

SUBDIRS = pdf

info_TEXINFOS = octopus.texi octopus-faq.texi octopus-devel.texi

man_MANS = octopus.1

CLEANFILES = *~ *.bak

octopus.1: octopus.pod
	pod2man -c "Octopus Documentation" -s 1 octopus.pod > octopus.1