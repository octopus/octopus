## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##

# ---------------------------------------------------------------
# Include paths.
# ---------------------------------------------------------------

FCFLAGS_MODS = \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/bpdn      \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/dftd3     \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/spglib-1.9.9/src/

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/include   \
	-I$(top_builddir)/src/include \
        -I$(top_srcdir)/external_libs/spglib-1.9.9/src \
	-I$(top_srcdir)/liboct_parser \
        $(GSL_CFLAGS) $(GD_CFLAGS) \
	@METIS_CFLAGS@ @PARMETIS_CFLAGS@ @CFLAGS_NFFT@ @CFLAGS_FFTW@ @CFLAGS_CUDA@ \
        @CGAL_CPPFLAGS@ \
	-DSHARE_DIR='"$(pkgdatadir)"'

AM_CCASFLAGS = \
	-I$(top_builddir)/

AM_CXXFLAGS = -I$(top_srcdir)/external_libs/rapidxml

# ---------------------------------------------------------------
# Define libraries here.
# ---------------------------------------------------------------

octopus_LIBS =

scalapack_LIBS = @LIBS_ELPA@ @LIBS_SCALAPACK@ @LIBS_BLACS@

core_LIBS = \
	@LIBS_FFTW@  @LIBS_LAPACK@ @LIBS_BLAS@                     \
	$(top_builddir)/liboct_parser/liboct_parser.a \
	@GSL_LIBS@ @LIBS_LIBXC@ @FCEXTRALIBS@

external_LIBS = \
	$(top_builddir)/external_libs/qshep/libqshep.a                  \
	$(top_builddir)/external_libs/spglib-1.9.9/src/libsymspg.a      \
	$(top_builddir)/external_libs/bpdn/libbpdn.a                    \
	$(top_builddir)/external_libs/dftd3/libdftd3.a
# we should not have libyaml here if we used an external one...

FCFLAGS_MODS += @FCFLAGS_LIBXC@ @FCFLAGS_PSPIO@ @FCFLAGS_PSOLVER@ @FCFLAGS_ISF@	\
  @FCFLAGS_FUTILE@ @FCFLAGS_FFTW@ @FCFLAGS_PFFT@ @FCFLAGS_PNFFT@		\
  @FCFLAGS_NETCDF@ @FCFLAGS_ETSF_IO@ @FCFLAGS_LIBVDWXC@ @FCFLAGS_BERKELEYGW@	\
  @FCFLAGS_NLOPT@ @FCFLAGS_LIBFM@ @FCFLAGS_ELPA@ @FCFLAGS_POKE@			\
  @FCFLAGS_LIKWID@

if COMPILE_OPENCL
  external_LIBS += $(top_builddir)/external_libs/fortrancl/libfortrancl.a @LIBS_CLBLAS@ @LIBS_CLFFT@ @CL_LIBS@
  FCFLAGS_MODS += @F90_MODULE_FLAG@$(top_builddir)/external_libs/fortrancl
endif

if COMPILE_METIS
  external_LIBS += $(top_builddir)/external_libs/metis-5.1/libmetis/libmetis.a
  external_LIBS += $(top_builddir)/external_libs/metis-5.1/GKlib/libgk.a
  AM_CPPFLAGS += -I$(top_srcdir)/external_libs/metis-5.1/include/
endif

if COMPILE_LIBYAML
  external_LIBS += $(top_builddir)/external_libs/yaml-0.1.4/src/libyaml.a
endif

# These must be arranged so if LIB1 depends on LIB2, LIB1 must occur before LIB2.
# e.g. ETSF_IO depends on netCDF, ISF depends on LAPACK
outside_LIBS = @LIBS_PSPIO@ @LIBS_POKE@ @LIBS_PSOLVER@ @LIBS_ISF@ @LIBS_FUTILE@	\
  @LIBS_LIBYAML@ @LIBS_NFFT@ @LIBS_PNFFT@ @LIBS_PFFT@ @LIBS_SPARSKIT@		\
  @LIBS_ETSF_IO@ @LIBS_NETCDF@ @LIBS_LIBFM@ @LIBS_LIBVDWXC@ @LIBS_BERKELEYGW@	\
  @LIBS_NLOPT@ @GD_LIBS@ @LIBS_PARMETIS@ @LIBS_METIS@ @LIBS_LIKWID@ @LIBS_CUDA@	\
  @LIBS_MPI@ @CGAL_LDFLAGS@

other_LIBS = $(external_LIBS) $(scalapack_LIBS) $(outside_LIBS) $(core_LIBS) @CXXLIBS@
all_LIBS = $(octopus_LIBS) $(other_LIBS)

# ---------------------------------------------------------------
# How to compile F90 files.
# ---------------------------------------------------------------

SUFFIXES = _oct.f90 .F90 .o .lo

# some definitions for silencing make
cpp_verbose = $(cpp_verbose_@AM_V@)
cpp_verbose_ = $(cpp_verbose_@AM_DEFAULT_V@)
cpp_verbose_0 = @echo "  CPP      $@";
fc_verbose = $(fc_verbose_@AM_V@)
fc_verbose_ = $(fc_verbose_@AM_DEFAULT_V@)
fc_verbose_0 = @echo "  FC       $@";


# Compilation is a two-step process: first we preprocess F90 files
# to generate _oct.f90 files. Then, we compile this _oct.f90 into
# an object file and delete the intermediate file.
.F90.o:
	$(cpp_verbose)@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< | \
	  $(top_srcdir)/build/preprocess.pl - \
	  "@DEBUG@" "@F90_ACCEPTS_LINE_NUMBERS@" > $*_oct.f90
	$(fc_verbose)@FC@ @FCFLAGS@ $(FCFLAGS_MODS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

.F90.lo:
	$(cpp_verbose)@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< | \
	  $(top_srcdir)/build/preprocess.pl - \
	  "@DEBUG@" "@F90_ACCEPTS_LINE_NUMBERS@" > $*_oct.f90
	$(fc_verbose)$(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	  $(LIBTOOLFLAGS) --mode=compile \
	  @FC@ @FCFLAGS@ $(FCFLAGS_MODS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90


# This rule is basically to create a _oct.f90 file by hand for
# debugging purposes. It is identical to the first part of
# the .F90.o rule.
.F90_oct.f90:
	$(cpp_verbose)@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< | \
	  $(top_srcdir)/build/preprocess.pl - \
	  "@DEBUG@" "@F90_ACCEPTS_LINE_NUMBERS@" > $*_oct.f90


# ---------------------------------------------------------------
# Miscellaneous.
# ---------------------------------------------------------------

# ctags.
CTAGS = ctags-exuberant -e

# Cleaning.
CLEANFILES = *~ *.bak *.mod *_oct.f90

# Local Variables:
# mode: Makefile
# coding: utf-8
# End:
