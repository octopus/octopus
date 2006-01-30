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


# Checks wether or not the emacs program exiss before attempting to launch it.
update:
	@(which emacs >/dev/null && emacs -batch -l update.el >/dev/null 2>/dev/null) \
	|| echo "WARNING: variables.texi not properly built because emacs is needed."

# If emacs was not launched in the update step, the variables.texi will be illegal,
# and so the makeinfo would crash.
octopus.info:
	@(which emacs >/dev/null && makeinfo octopus.texi) \
	|| echo "WARNING: octopus.info was not regenerated because emacs is needed."

octopus.info: variables.texi
octopus.ps: variables.texi
octopus.pdf: variables.texi
