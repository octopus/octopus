sharedir = $(pkgdatadir)
dist_share_DATA=../share/variables.xml octopus/ui/oct_main.jpg
CLEANFILES=../share/variables.xml

DESTDIR=/

../share/variables.xml: ../src/*.F90 var2xml.py
	cat ../src/*.F90 | grep '^[[:space:]]*!%' | sed -e's,.*!%,,' | ./var2xml.py > $@

all-local: ../share/variables.xml octopus-gui
	@PYTHON@ ./setup.py build

octopus-gui: octopus-gui.in
	sed -e 's,@pythondir_expanded@,$(pythondir),g' -e's,@datadir_expanded@,$(datadir),g' < $< > $@
	chmod +x $@

var2xml.py: var2xml.py.in
	sed -e 's,@pythondir_expanded@,$(pythondir),g' -e's,@datadir_expanded@,$(datadir),g' < $< > $@
	chmod +x $@

#pythonlibdir=$(shell python -c "from distutils.sysconfig import get_python_lib; print get_python_lib(0,False,'"$(prefix)"')")

install-data-local:
	@PYTHON@ ./setup.py install --install-purelib=$(pythondir) --install-scripts=$(bindir) --root $(DESTDIR)

clean-local:
	@PYTHON@ ./setup.py clean

SUBDIRS = po intl

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = config.rpath mkinstalldirs
