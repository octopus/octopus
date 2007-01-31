bin_SCRIPTS = oct-analyze_projections oct-vdW_c6
noinst_HEADERS = oct-analyze_projections.pl oct-vdW_c6.pl

CLEANFILES = oct-analyze_projections oct-vdW_c6

oct-analyze_projections: oct-analyze_projections.pl
	@cat $< | \
	  sed "s#<PERL>#${PERL}#" | \
	  sed "s#<BINDIR>#${bindir}#" > $@

oct-vdW_c6: oct-vdW_c6.pl
	@cat $< | \
	  sed "s#<PERL>#${PERL}#" | \
	  sed "s#<BINDIR>#${bindir}#" > $@
