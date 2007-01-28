# $Id$

SUBDIRS = pdf html

man_MANS = octopus.1 oct-atomic_occupations.1 oct-broad.1 oct-center-geom.1 oct-cross-section.1 \
        oct-display_partitions.1 oct-harmonic-spectrum.1 oct-ph.1 oct-rsf.1 oct-run_periodic_table.1 \
        oct-run_regression_test.pl.1 oct-run_testsuite.1 oct-xyz-anim.1

EXTRA_DIST = update.el *.pod img/*.eps img/*.jpg

CLEANFILES = *~ *.bak variables.texi

#octopus.1: octopus.pod
#	pod2man -c "Octopus Documentation" -s 1 octopus.pod > octopus.1
.pod.1:
	pod2man -c "octopus documentation" -s 1 $< > $@

../share/varinfo:
	(cd ../share; $(MAKE) varinfo)


