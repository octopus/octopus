INFODIR=/usr/share/info
INFODIRFILE=$(INFODIR)/dir
INSTALL-INFO=/sbin/install-info

all: dvi html info

dvi:
	texi2dvi octopus.texinfo

html:
	texi2html octopus.texinfo

info:
	makeinfo octopus.texinfo

install-info: info
	gzip octopus.info
	cp octopus.info.gz $(INFODIR)
	$(INSTALL-INFO) $(INFODIR)/octopus.info.gz $(INFODIRFILE)

uninstall-info:
	$(INSTALL-INFO) --delete $(INFODIR)/octopus.info.gz $(INFODIRFILE)
	rm /usr/share/info/octopus.info.gz

clean:
	@rm -f *.aux *.fn *.log *.tp *.cp *.pg *.vr *.ky *.ps *.toc *~

cleanall:
	@make clean
	@rm -f *.dvi *.info *.html *.pdf *.info.gz
