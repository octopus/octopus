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

# don't change anything from now on
all: dirs $(LIBOCT)
	cd src && $(MAKE) all EXEC_DIR=$(EXEC_DIR) CPPFLAGS="$(CPPFLAGS) -D THREE_D"

1d: dirs $(LIBOCT)
	cd src && $(MAKE) 1d EXEC_DIR=$(EXEC_DIR) CPPFLAGS="$(CPPFLAGS) -D ONE_D"

$(LIBOCT):
	cd liboct && $(MAKE) all LIB_DIR=$(LIB_DIR)

util: dirs
	cd src && $(MAKE) plan EXEC_DIR=$(EXEC_DIR)
#	cd src && $(MAKE) spectrum EXEC_DIR=$(EXEC_DIR)

dirs:
	for subdir in $(SUBDIRS) ; do \
	  test -d $$subdir || mkdir -p $$subdir || exit 1; \
	  chmod 777 $$subdir; \
	done

############################################################################
#     Clean stuff
############################################################################

cleanall: clean cleanlib

clean:
	cd src && $(MAKE) clean

cleanlib:
	cd liboct && $(MAKE) clean



