all:
	tex octopus.texinfo
	makeinfo octopus.texinfo
	texi2html octopus.texinfo

clean:
	@rm -f *.aux *.fn *.log *.tp *.cp *.pg *.vr *.ky *.ps *.toc *~

cleanall:
	@make clean
	@rm -f *.dvi *.info *.html *.pdf
