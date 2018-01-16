EXTRA_DIST = \
	PSolver_logo.png \
	extra.css

CLEANFILES = *.aux *.toc *.log *.xml


if BUILD_DEVDOC
insdevdocs = install-devdoc
uninsdevdocs = uninstall-devdoc
else
insdevdocs =
uninsdevdocs =
endif

install-data-hook: $(insdevdocs)

uninstall-local: $(uninsdevdocs)

install-devdoc: Doxygen/html
	test -z "$(docdir)" || /bin/mkdir -p "$(docdir)"
	cp -rp Doxygen/html '$(docdir)'

uninstall-devdoc:
	rm -rf '$(docdir)/html'

doxygen-doc: Doxygen/html

Doxygen/html: Doxygen.stamp web.stamp

Doxygen.stamp: $(top_builddir)/config/Doxyfile
	\rm -rf Doxygen;\
	doxygen $(top_builddir)/config/Doxyfile; \
	touch Doxygen.stamp

web.stamp: Doxygen.stamp
	touch web.stamp

ford-doc: FORD

FORD: FORD.stamp

FORD.stamp: $(top_builddir)/config/ford.md
	\rm -rf FORD
	ford --src_dir $(top_srcdir)/src \
		--include $(top_srcdir)/src --include $(top_builddir)/src \
		--exclude_dir $(top_srcdir)/src/unused \
		--output_dir FORD $(top_srcdir)/config/ford.md

clean-local:
	rm -rf FORD Doxygen lcov* *.aux *.log  *.toc Doxygen *.stamp doxygen.warnings

PHONY = install-devdoc
