# $Id$

SUBDIRS = pdf html

info_TEXINFOS = octopus.texi

octopus_TEXINFOS = authors.texi copying.texi examples.texi index.texi    \
	installation.texi introduction.texi parser.texi description.texi \
	undocumented.texi external.texi variables.texi

man_MANS = octopus.1

EXTRA_DIST = update.el octopus.pod img/*.eps img/*.jpg

CLEANFILES = *~ *.bak variables.texi

octopus.1: octopus.pod
	pod2man -c "Octopus Documentation" -s 1 octopus.pod > octopus.1

octopus.texi: ../share/varinfo variables.texi

variables.texi: ../share/varinfo
	(cd ../; build/var2html.pl)
	$(MAKE) update

../share/varinfo:
	(cd ../share; $(MAKE) varinfo)

update:
	emacs -batch -l update.el

octopus.info: variables.texi
octopus.ps: variables.texi
octopus.pdf: variables.texi
