###############################################################################
#          Makefile for the tddft program                                     #
###############################################################################

# This will give the hardware we are running on
ARCH   = $(shell uname -m)
POINTER_SIZE = $(shell bin/pointer_size)
PWD    = $(shell pwd)
export ARCH POINTER_SIZE

include Make.$(ARCH)

PREFIX   ?= $(PWD)
EXEC_DIR ?= $(PREFIX)/bin/$(ARCH)
LIB_DIR  ?= $(PREFIX)/lib/$(ARCH)
SUBDIRS = $(EXEC_DIR) $(LIB_DIR)

LIBOCT = $(LIB_DIR)/liboct.a
TDDFTLIBS:= -L$(LIB_DIR) $(FORLIBS)
export TDDFTLIBS
CPPFLAGS:=$(CPPFLAGS) -D POINTER_SIZE=$(POINTER_SIZE)

# Pseudopotential generation code
ATOM      = $(EXEC_DIR)/atm.x

#
# these names will be used to build the libraries
#
BLASLIB   = libblas.a
LAPACKLIB = liblapack.a
ARPACKLIB = libarpack.a
LASOLIB   = liblaso.a

# don't change anything from now on
all: dirs $(LIBOCT)
	cd src && $(MAKE) all EXEC_DIR=$(EXEC_DIR)

$(LIBOCT):
	cd liboct && $(MAKE) all LIB_DIR=$(LIB_DIR)
prep:
	cd src && $(MAKE) prep DEST_ARCH=alpha CPPFLAGS="$(CPPFLAGS) -D POINTER_SIZE=8 -D IN=inout"
	cd src && $(MAKE) prep DEST_ARCH=linux CPPFLAGS="$(CPPFLAGS) -D POINTER_SIZE=4"

old: dirs
	cd src.old && $(MAKE) all EXEC_DIR=$(EXEC_DIR)

libs: blas lapack arpack laso fftw

blas: dirs
	cd src-lib/LAPACK && $(MAKE) blaslib \
	  BLASLIB=$(LIB_DIR)/$(BLASLIB)

lapack: dirs
	cd src-lib/LAPACK && $(MAKE) lapacklib \
	  LAPACKLIB=$(LIB_DIR)/$(LAPACKLIB)

arpack: dirs
	cd src-lib/ARPACK && $(MAKE) \
	  ARPACKLIB=$(LIB_DIR)/$(ARPACKLIB)

laso: dirs
	cd src-lib/laso && $(MAKE) \
	  LASOLIB=$(LIB_DIR)/$(LASOLIB)

fftw: dirs
	(cd src-lib/fftw && ./configure --enable-type-prefix)
	(cd src-lib/fftw/fftw;$(MAKE) ; mv .libs/*.a $(LIB_DIR); $(MAKE) clean)
	(cd src-lib/fftw/rfftw;$(MAKE); mv .libs/*.a $(LIB_DIR); $(MAKE) clean)

util: dirs
	cd src && $(MAKE) plan EXEC_DIR=$(EXEC_DIR)
#	cd src && $(MAKE) spectrum EXEC_DIR=$(EXEC_DIR)
#	(cd pseudo/atom && $(MAKE) ATOM=$(ATOM))
#	$(FC) $(FFLAGS) -o $(EXEC_DIR)/vpsa2bin.x pseudo/atom/Util/vpsa2bin.f $(UTIL_LIBS)
#	$(FC) $(FFLAGS) -o $(EXEC_DIR)/vpsb2asc.x pseudo/atom/Util/vpsb2asc.f $(UTIL_LIBS)

dirs:
	for subdir in $(SUBDIRS) ; do \
	  test -d $$subdir || mkdir -p $$subdir || exit 1; \
	  chmod 777 $$subdir; \
	done

############################################################################
#     Clean stuff
############################################################################

cleanall: clean cleanlib cleanutil

clean:
	cd src && $(MAKE) clean

cleanold:
	cd src.old && $(MAKE) clean

cleanlib:
	cd src-lib/LAPACK && $(MAKE) cleanlib
	cd src-lib/ARPACK && $(MAKE) clean
	cd src-lib/laso   && $(MAKE) clean
	cd src-lib/fftw   && $(MAKE) clean

cleanutil:
	cd pseudo/atom    && $(MAKE) clean




