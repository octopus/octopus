bin_SCRIPTS = oct-analyze_projections
noinst_HEADERS = oct-analyze_projections.pl

CLEANFILES = oct-analyze_projections

oct-analyze_projections: oct-analyze_projections.pl
	@cat $< | \
	  sed "s#<PERL>#${PERL}#" | \
	  sed "s#<BINDIR>#${bindir}#" > $@
