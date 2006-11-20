bin_SCRIPTS = oct_analyze_projections
noinst_HEADERS = oct_analyze_projections.pl

CLEANFILES = oct_analyze_projections

oct_analyze_projections: oct_analyze_projections.pl
	@cat $< | \
	  sed "s#<PERL>#${PERL}#" | \
	  sed "s#<BINDIR>#${bindir}#" > $@
