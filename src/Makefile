# Generated automatically from Makefile.in by configure.
# Makefile.in generated automatically by automake 1.4-p5 from Makefile.am

# Copyright (C) 1994, 1995-8, 1999, 2001 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.


SHELL = /bin/sh

srcdir = .
top_srcdir = ..
prefix = /home/nfs/marques/work/octopus
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

DESTDIR =

pkgdatadir = $(datadir)/octopus
pkglibdir = $(libdir)/octopus
pkgincludedir = $(includedir)/octopus

top_builddir = ..

ACLOCAL = aclocal
AUTOCONF = autoconf
AUTOMAKE = './automake'
AUTOHEADER = autoheader

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL} $(AM_INSTALL_PROGRAM_FLAGS)
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL_PROGRAM}
transform = s,x,x,

NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
host_alias = i686-pc-linux-gnu
host_triplet = i686-pc-linux-gnu
AS = @AS@
AUTOMAKE = './automake'
CC = gcc
CPP = gcc -E
DLLTOOL = @DLLTOOL@
ECHO = echo
EXEEXT = 
F90 = f90
F90CPP = /lib/cpp
F90CPPFLAGS = -C
F90FLAGS = -YEXT_NAMES=LCS -YEXT_SFX=_
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LIB_FFTW = -L/home/nfs/marques -L/home/nfs/marques/lib -ldfftw -ldrfftw
LIB_GSL = -lgsl -lgslcblas
LIB_LAPACK = -L/home/nfs/marques -L/home/nfs/marques/lib -llapack -lblas
LN_S = ln -s
MAKEINFO = makeinfo
OBJDUMP = @OBJDUMP@
OBJEXT = o
PACKAGE = octopus
RANLIB = ranlib
STRIP = strip
VERSION = 1.0
YACC = bison -y

bin_PROGRAMS = octopus

octopus_SOURCES = global.F90 liboct.F90 io.F90 pulpo.F90 units.F90 math.F90 \
	spline.F90 fft.F90 ps_xc.F90 ps_kb.F90 ps.F90 specie.F90 atom.F90 \
	mesh.F90 hartree.F90 states.F90 xc.F90 systm.F90 lasers.F90 hamiltonian.F90 \
	eigen.F90 mix.F90 lcao.F90 scf.F90 unocc.F90 static_pol.F90 \
	td.F90 run.F90 main.F90


EXTRA_octopus_SOURCES = \
	mesh1D.F90 mesh1D_create.F90 mesh1D_inc.F90 ps1D.F90 hartree1D.F90 \
	mesh3D.F90 mesh3D_create.F90 mesh3D_inc.F90 ps3D.F90 hartree3D.F90 \
	xc_pot.F90 xc_LDA.F90 xc_GGA.F90 xc_KLI.F90 xc_KLI_x.F90 xc_KLI_SIC.F90 \
	hamiltonian_inc.F90 h_external_pot.F90 h_forces.F90 \
	eigen_cg1.F90 eigen_cg2.F90 \
	states_inc.F90 \
	td_init.F90 td_rti.F90 \
	mix_inc.F90


octopus_LDADD = $(LIB_FFTW) $(LIB_GSL) $(LIB_LAPACK) ../liboct/liboct.a

CLEANFILES = *~ *.bak *.mod
DISTCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_HEADER = ../liboct/config.h config.h
CONFIG_CLEAN_FILES = 
bin_PROGRAMS =  octopus$(EXEEXT)
PROGRAMS =  $(bin_PROGRAMS)


DEFS = -DHAVE_CONFIG_H -I. -I$(srcdir) -I../liboct -I.
CPPFLAGS = 
LDFLAGS = 
LIBS = 
octopus_OBJECTS =  global.$(OBJEXT) liboct.$(OBJEXT) io.$(OBJEXT) \
pulpo.$(OBJEXT) units.$(OBJEXT) math.$(OBJEXT) spline.$(OBJEXT) \
fft.$(OBJEXT) ps_xc.$(OBJEXT) ps_kb.$(OBJEXT) ps.$(OBJEXT) \
specie.$(OBJEXT) atom.$(OBJEXT) mesh.$(OBJEXT) hartree.$(OBJEXT) \
states.$(OBJEXT) xc.$(OBJEXT) systm.$(OBJEXT) lasers.$(OBJEXT) \
hamiltonian.$(OBJEXT) eigen.$(OBJEXT) mix.$(OBJEXT) lcao.$(OBJEXT) \
scf.$(OBJEXT) unocc.$(OBJEXT) static_pol.$(OBJEXT) td.$(OBJEXT) \
run.$(OBJEXT) main.$(OBJEXT)
octopus_DEPENDENCIES =  ../liboct/liboct.a
octopus_LDFLAGS = 
F90COMPILE = $(F90) $(AM_F90FLAGS) $(F90FLAGS)
LTF90COMPILE = $(LIBTOOL) --mode=compile $(F90) $(AM_F90FLAGS) $(F90FLAGS)
F90LD = $(F90)
F90LINK = $(LIBTOOL) --mode=link $(F90LD) $(AM_F90FLAGS) $(F90FLAGS) $(LDFLAGS) -o $@
DIST_COMMON =  ./stamp-h2.in Makefile.am Makefile.in config.h.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = gtar
GZIP_ENV = --best
DEP_FILES =  .deps/atom.P .deps/eigen.P .deps/eigen_cg1.P \
.deps/eigen_cg2.P .deps/fft.P .deps/global.P .deps/h_external_pot.P \
.deps/h_forces.P .deps/hamiltonian.P .deps/hamiltonian_inc.P \
.deps/hartree.P .deps/hartree1D.P .deps/hartree3D.P .deps/io.P \
.deps/lasers.P .deps/lcao.P .deps/liboct.P .deps/main.P .deps/math.P \
.deps/mesh.P .deps/mesh1D.P .deps/mesh1D_create.P .deps/mesh1D_inc.P \
.deps/mesh3D.P .deps/mesh3D_create.P .deps/mesh3D_inc.P .deps/mix.P \
.deps/mix_inc.P .deps/ps.P .deps/ps1D.P .deps/ps3D.P .deps/ps_kb.P \
.deps/ps_xc.P .deps/pulpo.P .deps/run.P .deps/scf.P .deps/specie.P \
.deps/spline.P .deps/states.P .deps/states_inc.P .deps/static_pol.P \
.deps/systm.P .deps/td.P .deps/td_init.P .deps/td_rti.P .deps/units.P \
.deps/unocc.P .deps/xc.P .deps/xc_GGA.P .deps/xc_KLI.P \
.deps/xc_KLI_SIC.P .deps/xc_KLI_x.P .deps/xc_LDA.P .deps/xc_pot.P
SOURCES = $(octopus_SOURCES) $(EXTRA_octopus_SOURCES)
OBJECTS = $(octopus_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .F90 .S .c .f90 .lo .o .obj .s
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu src/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status $(BUILT_SOURCES)
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status


config.h: stamp-h2
	@if test ! -f $@; then \
		rm -f stamp-h2; \
		$(MAKE) stamp-h2; \
	else :; fi
stamp-h2: $(srcdir)/config.h.in $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES= CONFIG_HEADERS=src/config.h \
	     $(SHELL) ./config.status
	@echo timestamp > stamp-h2 2> /dev/null
$(srcdir)/config.h.in: $(srcdir)/stamp-h2.in
	@if test ! -f $@; then \
		rm -f $(srcdir)/stamp-h2.in; \
		$(MAKE) $(srcdir)/stamp-h2.in; \
	else :; fi
$(srcdir)/stamp-h2.in: $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOHEADER)
	@echo timestamp > $(srcdir)/stamp-h2.in 2> /dev/null

mostlyclean-hdr:

clean-hdr:

distclean-hdr:
	-rm -f config.h

maintainer-clean-hdr:

mostlyclean-binPROGRAMS:

clean-binPROGRAMS:
	-test -z "$(bin_PROGRAMS)" || rm -f $(bin_PROGRAMS)

distclean-binPROGRAMS:

maintainer-clean-binPROGRAMS:

install-binPROGRAMS: $(bin_PROGRAMS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	@list='$(bin_PROGRAMS)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo " $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`"; \
	    $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	  else :; fi; \
	done

uninstall-binPROGRAMS:
	@$(NORMAL_UNINSTALL)
	list='$(bin_PROGRAMS)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	done

# FIXME: We should only use cygpath when building on Windows,
# and only if it is available.
.c.obj:
	$(COMPILE) -c `cygpath -w $<`

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

mostlyclean-compile:
	-rm -f *.o core *.core
	-rm -f *.$(OBJEXT)

clean-compile:

distclean-compile:
	-rm -f *.tab.c

maintainer-clean-compile:

.s.lo:
	$(LIBTOOL) --mode=compile $(COMPILE) -c $<

.S.lo:
	$(LIBTOOL) --mode=compile $(COMPILE) -c $<

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

distclean-libtool:

maintainer-clean-libtool:

octopus$(EXEEXT): $(octopus_OBJECTS) $(octopus_DEPENDENCIES)
	@rm -f octopus$(EXEEXT)
	$(F90LINK) $(octopus_LDFLAGS) $(octopus_OBJECTS) $(octopus_LDADD) $(LIBS)
.f90.o:
	$(F90COMPILE) -c $<
.F90.f90:
	$(F90CPP) $(F90CPPFLAGS) $< > $@
.f90.obj:
	$(F90COMPILE) -c `cygpath -w $<`
.f90.lo:
	$(LTF90COMPILE) -c $<

tags: TAGS

ID: $(HEADERS) $(SOURCES) $(LISP)
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	here=`pwd` && cd $(srcdir) \
	  && mkid -f$$here/ID $$unique $(LISP)

TAGS:  $(HEADERS) $(SOURCES) config.h.in $(TAGS_DEPENDENCIES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	test -z "$(ETAGS_ARGS)config.h.in$$unique$(LISP)$$tags" \
	  || (cd $(srcdir) && etags $(ETAGS_ARGS) $$tags config.h.in $$unique $(LISP) -o $$here/TAGS)

mostlyclean-tags:

clean-tags:

distclean-tags:
	-rm -f TAGS ID

maintainer-clean-tags:

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = src

distdir: $(DISTFILES)
	here=`cd $(top_builddir) && pwd`; \
	top_distdir=`cd $(top_distdir) && pwd`; \
	distdir=`cd $(distdir) && pwd`; \
	cd $(top_srcdir) \
	  && $(AUTOMAKE) --include-deps --build-dir=$$here --srcdir-name=$(top_srcdir) --output-dir=$$top_distdir --gnu src/Makefile
	@for file in $(DISTFILES); do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    cp -pr $$d/$$file $(distdir)/$$file; \
	  else \
	    test -f $(distdir)/$$file \
	    || ln $$d/$$file $(distdir)/$$file 2> /dev/null \
	    || cp -p $$d/$$file $(distdir)/$$file || :; \
	  fi; \
	done

DEPS_MAGIC := $(shell mkdir .deps > /dev/null 2>&1 || :)

-include $(DEP_FILES)

mostlyclean-depend:

clean-depend:

distclean-depend:
	-rm -rf .deps

maintainer-clean-depend:

%.o: %.c
	@echo '$(COMPILE) -c $<'; \
	$(COMPILE) -Wp,-MD,.deps/$(*F).pp -c $<
	@-cp .deps/$(*F).pp .deps/$(*F).P; \
	tr ' ' '\012' < .deps/$(*F).pp \
	  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
	    >> .deps/$(*F).P; \
	rm .deps/$(*F).pp

%.lo: %.c
	@echo '$(LTCOMPILE) -c $<'; \
	$(LTCOMPILE) -Wp,-MD,.deps/$(*F).pp -c $<
	@-sed -e 's/^\([^:]*\)\.o[ 	]*:/\1.lo \1.o :/' \
	  < .deps/$(*F).pp > .deps/$(*F).P; \
	tr ' ' '\012' < .deps/$(*F).pp \
	  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
	    >> .deps/$(*F).P; \
	rm -f .deps/$(*F).pp
info-am:
info: info-am
dvi-am:
dvi: dvi-am
check-am: all-am
check: check-am
installcheck-am:
installcheck: installcheck-am
all-recursive-am: config.h
	$(MAKE) $(AM_MAKEFLAGS) all-recursive

install-exec-am: install-binPROGRAMS
install-exec: install-exec-am

install-data-am:
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am: uninstall-binPROGRAMS
uninstall: uninstall-am
all-am: Makefile $(PROGRAMS) config.h
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:
	$(mkinstalldirs)  $(DESTDIR)$(bindir)


mostlyclean-generic:

clean-generic:
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*
	-test -z "$(DISTCLEANFILES)" || rm -f $(DISTCLEANFILES)

maintainer-clean-generic:
mostlyclean-am:  mostlyclean-hdr mostlyclean-binPROGRAMS \
		mostlyclean-compile mostlyclean-libtool \
		mostlyclean-tags mostlyclean-depend mostlyclean-generic

mostlyclean: mostlyclean-am

clean-am:  clean-hdr clean-binPROGRAMS clean-compile clean-libtool \
		clean-tags clean-depend clean-generic mostlyclean-am

clean: clean-am

distclean-am:  distclean-hdr distclean-binPROGRAMS distclean-compile \
		distclean-libtool distclean-tags distclean-depend \
		distclean-generic clean-am
	-rm -f libtool

distclean: distclean-am

maintainer-clean-am:  maintainer-clean-hdr maintainer-clean-binPROGRAMS \
		maintainer-clean-compile maintainer-clean-libtool \
		maintainer-clean-tags maintainer-clean-depend \
		maintainer-clean-generic distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-am

.PHONY: mostlyclean-hdr distclean-hdr clean-hdr maintainer-clean-hdr \
mostlyclean-binPROGRAMS distclean-binPROGRAMS clean-binPROGRAMS \
maintainer-clean-binPROGRAMS uninstall-binPROGRAMS install-binPROGRAMS \
mostlyclean-compile distclean-compile clean-compile \
maintainer-clean-compile mostlyclean-libtool distclean-libtool \
clean-libtool maintainer-clean-libtool tags mostlyclean-tags \
distclean-tags clean-tags maintainer-clean-tags distdir \
mostlyclean-depend distclean-depend clean-depend \
maintainer-clean-depend info-am info dvi-am dvi check check-am \
installcheck-am installcheck all-recursive-am install-exec-am \
install-exec install-data-am install-data install-am install \
uninstall-am uninstall all-redirect all-am all installdirs \
mostlyclean-generic distclean-generic clean-generic \
maintainer-clean-generic clean mostlyclean distclean maintainer-clean


# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
