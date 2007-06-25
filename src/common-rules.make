#---------------------------------------------------------------
# include paths
#---------------------------------------------------------------
AM_FCFLAGS = \
	@F90_MODULE_FLAG@$(top_builddir)/src/basic   \
	@F90_MODULE_FLAG@$(top_builddir)/src/math    \
	@F90_MODULE_FLAG@$(top_builddir)/src/species \
	@F90_MODULE_FLAG@$(top_builddir)/src/grid    \
	@F90_MODULE_FLAG@$(top_builddir)/src/poisson \
	@F90_MODULE_FLAG@$(top_builddir)/src/states  \
	@F90_MODULE_FLAG@$(top_builddir)/src/xc      \
	@F90_MODULE_FLAG@$(top_builddir)/src/h_sys   \
	@F90_MODULE_FLAG@$(top_builddir)/src/scf     \
	@F90_MODULE_FLAG@$(top_builddir)/src/td      \
	@F90_MODULE_FLAG@$(top_builddir)/src/opt_control \
	@F90_MODULE_FLAG@$(top_builddir)/src/sternheimer \
	@F90_MODULE_FLAG@$(top_builddir)/libxc/src   \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/qshep

AM_CPPFLAGS = -I$(srcdir)/include -I$(top_builddir)/src/include -I$(top_srcdir)/libstring_f


#---------------------------------------------------------------
# define libraries here
#---------------------------------------------------------------
octopus_LIBS = \
	$(top_builddir)/src/sternheimer/libsternheimer.a \
	$(top_builddir)/src/opt_control/libopt_control.a \
	$(top_builddir)/src/td/libtd.a                   \
	$(top_builddir)/src/scf/libscf.a                 \
	$(top_builddir)/src/h_sys/libh_sys.a             \
	$(top_builddir)/src/xc/libxc.a                   \
	$(top_builddir)/src/states/libstates.a           \
	$(top_builddir)/src/poisson/libpoisson.a         \
	$(top_builddir)/src/grid/libgrid.a               \
	$(top_builddir)/src/species/libspecies.a         \
	$(top_builddir)/src/math/libmath.a               \
	$(top_builddir)/src/basic/libbasic.a

core_LIBS = \
	$(octopus_LIBS) \
	@LIBS_LAPACK@ @LIBS_BLAS@ \
	$(top_builddir)/liboct/liboct.a \
	$(top_builddir)/liboct_parser/liboct_parser.a \
	$(top_builddir)/libxc/src/libxc.a \
	-L$(top_builddir)/libstring_f -lstring_f \
	@GSL_LIBS@ @GD_LIBS@ @FCEXTRALIBS@

external_LIBS = \
	$(top_builddir)/external_libs/expokit/libexpokit.a \
	$(top_builddir)/external_libs/qshep/libqshep.a \
	$(top_builddir)/external_libs/poisson_isf/libpoisson_isf.a

if COMPILE_METIS
  external_LIBS += $(top_builddir)/external_libs/metis-4.0/libmetis.a
endif

if COMPILE_LIBNBC
  external_LIBS += $(top_builddir)/external_libs/libnbc/libnbc.a
endif

all_LIBS = $(core_LIBS) @LIBS_FFT@ @LIBS_TRLAN@ @LIBS_ARPACK@ @LIBS_SPARSKIT@ \
  @LIBS_NETCDF@ $(external_LIBS) @LIBS_MPI@ @LIBS_LAPACK@ @LIBS_BLAS@


#---------------------------------------------------------------
# How to compile F90 files
#---------------------------------------------------------------

# Define empty rule to override the native rule of automake (otherwise the native
# rule will have precedence over the pattern based rules below). Could not figure out
# how to remove .F90 from .SUFFIXES, which would be an alternative route.
.F90.o:

# compilation is a two step process: first we preprocess F90 files to generate _oct.f90 files
%_oct.f90: %.F90
	@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< > $*_oct.f90
	@if [ "@DEBUG@" = "no" ]; then \
		cat $*_oct.f90 | grep -v pop_sub | grep -v push_sub >$*_oct.f91; \
		mv -f $*_oct.f91 $*_oct.f90; \
	fi
	@perl -pi -e 's/\\newline/\n/g; s/\\cardinal/#/g' $*_oct.f90

# then we actually compile _oct.f90 files and cleanup
%.o: %_oct.f90
	@FC@ @FCFLAGS@ @FCFLAGS_NETCDF@ $(AM_FCFLAGS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# ctags
CTAGS = ctags-exuberant -e

# cleaning
CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* *_oct.f90 config_F90.h

#---------------------------------------------------------------
# if it exists include user defined Makefile.local (which is not included in the octopus package)
#---------------------------------------------------------------
-include Makefile.local


# Local Variables:
# mode: Makefile
# coding: utf-8
# End:
