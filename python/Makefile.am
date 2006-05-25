sharedir = $(pkgdatadir)
dist_share_DATA=../share/variables.xml
CLEANFILES=../share/variables.xml

DESTDIR=/

../share/variables.xml: ../src/*.F90
	cat ../src/*.F90 | grep '^[[:space:]]*!%' | sed -e's,.*!%,,' | ./var2xml.py > $@

all-local: ../share/variables.xml
	@PYTHON@ ./setup.py build

pythonlibdir=$(shell python -c "from distutils.sysconfig import get_python_lib; print get_python_lib(0,False,'"$(prefix)"')")

install-data-local:
	@PYTHON@ ./setup.py install --install_purelib=$(pythonlibdir) --install-scripts=$(bindir) --root $(DESTDIR)

clean-local:
	@PYTHON@ ./setup.py clean
