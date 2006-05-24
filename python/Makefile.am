sharedir = $(pkgdatadir)
dist_share_DATA=../share/variables.xml
CLEANFILES=../share/variables.xml

DESTDIR=/

../share/variables.xml: ../src/*.F90
	cat ../src/*.F90 | grep '^[[:space:]]*!%' | sed -e's,.*!%,,' | ./var2xml.py > $@

all-local: ../share/variables.xml
	@PYTHON@ ./setup.py build

install-data-local:
	@PYTHON@ ./setup.py install --root $(DESTDIR)

clean-local:
	@PYTHON@ ./setup.py clean
