## Copyright (C) 2002-2018 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

AC_PREREQ(2.59)
# Note: the version string is written into ETSF_IO files, so when you update
# this, update the results for test periodic_systems/09-etsf_io.test too.
AC_INIT([Octopus],[arborescens],[octopus-devel@tddft.org],[octopus])
AC_CONFIG_SRCDIR([src/basic/global.F90])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects tar-ustar])
# enable silent rules by default (use --disable-silent-rules for disabling)
AM_SILENT_RULES([yes])

AC_CONFIG_HEADERS([config.h])

# Installation prefix by default
AC_PREFIX_DEFAULT([/usr/local])

# who am i
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P

dnl check whether mpi is enabled
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--enable-mpi(=PATH)], [Parallel version. For ancient versions of MPI, you may need to specify PATH for MPI libs.]))
case $enable_mpi in
  yes) ;;
  no | "") enable_mpi=no ;;
  -* | */* | *.a | *.so | *.so.* | *.o)
    LIBS_MPI="$enable_mpi"
    enable_mpi=yes
    ;;
  *) 
    LIBS_MPI="-l$enable_mpi"
    enable_mpi=yes
    ;;
esac
AM_CONDITIONAL(USE_MPI, test "$enable_mpi" = "yes")

if test x"$enable_mpi" == x"yes"; then
  octopus_default_cc=mpicc
else
  octopus_default_cc=gcc
fi

AC_PROG_CC($octopus_default_cc)
AC_PROG_CC_C99
AM_PROG_CC_C_O

AC_PROG_CXX
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])
AC_CXX_LIBRARY_LDFLAGS

# optimizations extensions supported by the compiler
AC_C_INLINE
AC_C_RESTRICT
ACX_C_BUILTIN_PREFETCH

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([signal.h errno.h])
AC_FUNC_ALLOCA

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_CHECK_SIZEOF(void*)
AC_CHECK_SIZEOF(size_t)
AC_CHECK_TYPES([uint32_t, uint64_t])
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(unsigned long)
AC_CHECK_SIZEOF(unsigned long long)
AX_TLS

# Checks for library functions.
AC_PROG_GCC_TRADITIONAL
AC_FUNC_STAT
AC_CHECK_FUNCS([gettimeofday nanosleep sysconf getpid \
strndup tcgetpgrp scandir alphasort \
sigaction readdir closedir strchr strtod \
perror strcasestr sbrk getopt_long uname ioctl strsignal])

ACX_ARCH

dnl ============================================================================
dnl This part of the configuration file tries to sort out what needs to be built.
dnl In particular, the code can be compiled in debug mode,
dnl and with OpenMP threads.

dnl these are the --enable/--disable possible in octopus
AC_ARG_ENABLE(debug, AS_HELP_STRING([--disable-debug], [debug mode]), [DEBUG=${enableval}])
AC_ARG_ENABLE(openmp, AS_HELP_STRING([--enable-openmp], [Multi-threaded parallel version using OpenMP from Fortran]))
AC_ARG_ENABLE(opencl, AS_HELP_STRING([--enable-opencl], [OpenCL support (experimental)]))

if test -n "${enable_openmp}" -a x"${enable_openmp}" != x"no"; then
  HAVE_OPENMP=1
fi



AC_ARG_ENABLE(nvtx, AS_HELP_STRING([--enable-nvtx], [NVTX profiling support for Cuda version]))
if test x"${enable_nvtx}" == x"yes" ; then
  AC_DEFINE(HAVE_NVTX, 1, [defined if nvtx support is enabled])
fi

AC_ARG_ENABLE(cuda, AS_HELP_STRING([--enable-cuda], [Cuda support (experimental)]))

if test x"${enable_cuda}" == x"yes" ; then
  AC_DEFINE(HAVE_CUDA, 1, [defined if cuda support is enabled])
  ACX_CUDA
fi

if test x"${enable_opencl}" == x"yes" ; then
  AC_MSG_WARN([
  
	Experimental OpenCL support enabled.
  ])
  AX_CHECK_CL
  if test "X$no_cl" == "Xyes"; then
    AC_MSG_FAILURE([Cannot find the OpenCL library.])
  fi
  AC_DEFINE(HAVE_OPENCL, 1, [defined if opencl support is enabled])
  ACX_CLBLAS
  ACX_CLFFT
fi

AM_CONDITIONAL(COMPILE_OPENCL, test x"${enable_opencl}" == x"yes")

dnl check for C MPI
if test x"$enable_mpi" == x"yes"; then
  AC_CHECK_HEADER(mpi.h,,AC_MSG_ERROR([cannot find mpi.h]))
  ACX_MPI([], AC_MSG_ERROR([cannot call MPI from C]))
fi

dnl build/install Perl bits unless --disable-perl is used
AC_ARG_ENABLE([perl],
	      AS_HELP_STRING([--disable-perl], [don't build/install Perl components]),
	      [ac_cv_build_perl=$enableval],
	      [ac_cv_build_perl=yes])
if test "x$ac_cv_build_perl" = "xyes"; then
  AC_PATH_PROG([PERL], [perl], [], [])
  if test -z "$PERL" ; then
     AC_MSG_WARN([perl not found. Disabling Perl parts])
     ac_cv_build_perl="no"
  fi
fi
AM_CONDITIONAL([ENABLE_PERL], [test x"$ac_cv_build_perl" = x"yes"])


dnl try to find out what is the default FORTRAN 90 compiler
dnl note that this has to be done after the AC_ARG_ENABLE(mpi...)
acx_save_fcflags="${FCFLAGS}"
if test x"$enable_mpi" = x"yes"; then
  AC_PROG_FC([mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r], [Fortran 90])
else
  AC_PROG_FC([], Fortran 90)
fi
if test x"$FC" = x; then
  AC_MSG_ERROR([could not find Fortran 90 compiler])
fi
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(f90)
FCFLAGS="${acx_save_fcflags}"

ACX_FCFLAGS

dnl octopus needs the preprocessor. The result goes to FCCPP
ACX_FCCPP

AC_FC_LIBRARY_LDFLAGS

dnl how Fortran mangles function names
AC_FC_WRAPPERS

ACX_FC_INTEGER_SIZE
ACX_CC_FORTRAN_INT

dnl Sometimes Fortran may use the flush function, sometimes not.
ACX_FORTRAN_CHECK_FUNC(flush)

dnl Some compilers define the sizeof intrinsic
ACX_FC_SIZEOF

ACX_FORTRAN_LOC

dnl check whether we have or not a compiler that allows for very long lines...
ACX_LONG_FORTRAN_LINES

dnl check whether the compiler supports compiler_version intrinsic.
ACX_FC_COMPILER_VERSION

dnl check whether the Fortran 90 compiler accepts line markers cast by the preprocessor.
ACX_F90_ACCEPTS_LINE_NUMBERS

dnl check whether the Fortran compiler supports command-line arguments
ACX_FC_COMMAND_LINE_ARGUMENTS

dnl check whether the Fortran compiler supports Fortran 2003's iso_c_binding
ACX_FC_ISO_C_BINDING([], AC_MSG_ERROR([could not find required Fortran 2003 iso_c_binding]))

dnl check flag to include module path
AX_F90_MODULE_FLAG
F90_MODULE_FLAG=$ax_cv_f90_modflag
AC_SUBST([F90_MODULE_FLAG])

if test "${HAVE_OPENMP}"; then
  AX_OPENMP(,AC_MSG_ERROR([Fortran compiler does not support OpenMP]))
  FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
fi


dnl debug mode
if test x"${DEBUG}" = x"no" ; then
  AC_DEFINE(NDEBUG, yes, [octopus compiled without debug mode])
fi
AC_SUBST(DEBUG)

UTILITY_PROGRAMS='$(EXTRA_PROGRAMS)'
AC_SUBST(UTILITY_PROGRAMS)

AC_ARG_WITH(max-dim, [  --with-max-dim=NDIM     maximum number of dimensions Octopus can use;
                          the default is 3 (other values are experimental), must be >= 3],
  [maxdim=$withval], [maxdim=3])

if test $maxdim -lt 3; then
  AC_MSG_ERROR([max-dim must be at least 3])
fi

AC_DEFINE_UNQUOTED(MAX_DIM, $maxdim, [the maximum dimension of the space])

dnl ===================================================================
dnl now we search for libraries

ACX_LIBXC
ACX_LIBVDWXC

dnl check for BLAS
ACX_BLAS([], AC_MSG_ERROR([could not find required BLAS library]))

dnl test zdotc in BLAS
ACX_ZDOTC

dnl check for LAPACK
ACX_LAPACK([], AC_MSG_ERROR([could not find required LAPACK library]))

dnl check for GSL
AC_LANG_PUSH(C)
AX_PATH_GSL

dnl check for GD library
ACX_GDLIB

AC_LANG_POP(C)

dnl check for FFTW
ACX_FFTW

dnl check for PFFT
ACX_PFFT

dnl check for NFFT
AC_LANG_PUSH(C)
ACX_NFFT
AC_LANG_POP(C)

dnl check for PNFFT
ACX_PNFFT

dnl check for NETCDF
ACX_NETCDF

dnl check for ETSF_IO
ACX_ETSF_IO

dnl check for BerkeleyGW
ACX_BERKELEYGW

dnl check for SPARSKIT
ACX_SPARSKIT

dnl check for NLOPT
ACX_NLOPT
  
dnl check for Fortran MPI support
if test x"$enable_mpi" != x"no"; then
  ACX_MPI([], AC_MSG_ERROR([could not compile an MPI test program]))
  ACX_MPI_FC_MODULE
  ACX_MPI2
fi

dnl check for SCALAPACK/BLACS
ACX_BLACS
ACX_SCALAPACK
ACX_ELPA

dnl check for Libfm
ACX_LIBFM

dnl Check for PSOlver library.
dnl PSolver depends on LibYAML and on Futile, so we first check for those one.
AC_LANG_PUSH(C)
ACX_YAML
AC_LANG_POP(C)
ACX_FUTILE
ACX_PSOLVER

dnl check for CGAL library.
AX_LIB_CGAL_CORE

dnl If PSolver was not found, we check for older versions of the library (1.7.x),
dnl which we named ISF
if test "x$acx_psolver_ok" = "xno"; then
  ACX_ISF
fi
AM_CONDITIONAL(COMPILE_LIBYAML, test "x${HAVE_ISF}" == "x1")

dnl check for PSPIO
ACX_PSPIO

dnl check for the Poke library to solve the Poisson equation
dnl http://www.e-cam2020.eu:10080/ESL/poke
ACX_POKE

dnl check for compiled METIS5 and ParMETIS (also checks METIS in the same path)
AC_LANG_PUSH(C)
ACX_PATH_METIS
ACX_PATH_PARMETIS
AC_LANG_POP(C)
AM_CONDITIONAL(COMPILE_METIS, test x${HAVE_COMP_METIS} = x1)

dnl performance analysis using the likwid marker API
ACX_LIKWID

ACX_DFTBPLUS

AC_SUBST(FCEXTRALIBS)

m4_include([src/fdep/fortran_dependencies.m4])
FDEP_F90_GNU_MAKE_DEPS

LT_INIT([disable-shared])

ACX_OCTOPUS_COMPILATION_INFO

#  external_libs/isf/Makefile
#  external_libs/isf/src/Makefile
#  external_libs/isf/flib/Makefile
#  external_libs/isf/wrappers/Makefile

AC_CONFIG_FILES([Makefile
  external_libs/Makefile
  external_libs/dftd3/Makefile
  external_libs/metis-5.1/Makefile
  external_libs/metis-5.1/include/Makefile
  external_libs/metis-5.1/GKlib/Makefile
  external_libs/metis-5.1/libmetis/Makefile
  external_libs/qshep/Makefile
  external_libs/fortrancl/Makefile
  external_libs/bpdn/Makefile
  external_libs/yaml-0.1.4/Makefile
  external_libs/yaml-0.1.4/include/Makefile
  external_libs/yaml-0.1.4/src/Makefile
  external_libs/rapidxml/Makefile
  liboct_parser/Makefile
  src/Makefile
  scripts/Makefile
  doc/Makefile
  doc/html/Makefile
  doc/html/icons/Makefile
  doc/html/js/Makefile
  doc/html/vars/Makefile
  share/Makefile
  share/pseudopotentials/Makefile
  share/pseudopotentials/PSF/Makefile
  share/pseudopotentials/HGH/Makefile
  share/pseudopotentials/HGH/lda/Makefile
  share/pseudopotentials/HGH/lda_sc/Makefile
  share/pseudopotentials/quantum-simulation.org/Makefile
  share/pseudopotentials/quantum-simulation.org/hscv/Makefile
  share/pseudopotentials/quantum-simulation.org/hscv/lda/Makefile
  share/pseudopotentials/quantum-simulation.org/hscv/pbe/Makefile
  share/pseudopotentials/quantum-simulation.org/sg15/Makefile
  share/pseudopotentials/pseudo-dojo.org/Makefile  
  share/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/Makefile
  share/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_stringent/Makefile
  share/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbesol_standard/Makefile
  share/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbesol_stringent/Makefile
  share/pseudopotentials/pseudo-dojo.org/nc-sr-04_pw_standard/Makefile
  share/pseudopotentials/pseudo-dojo.org/nc-sr-04_pw_stringent/Makefile
  share/recipes/Makefile
  share/recipes/da/Makefile
  share/recipes/de/Makefile
  share/recipes/en/Makefile
  share/recipes/es/Makefile
  share/recipes/it/Makefile
  share/recipes/eu/Makefile
  share/util/Makefile
  share/opencl/Makefile
  build/Makefile
  build/octopus.spec
  testsuite/Makefile
  testsuite/performance/Makefile
  testsuite/performance/tests/Makefile
  testsuite/components/Makefile 
  testsuite/errors/Makefile 
  testsuite/finite_systems_1d/Makefile 
  testsuite/finite_systems_2d/Makefile
  testsuite/finite_systems_3d/Makefile
  testsuite/functionals/Makefile
  testsuite/real_time/Makefile
  testsuite/lda_u/Makefile
  testsuite/linear_response/Makefile
  testsuite/maxwell/Makefile
  testsuite/modelmb/Makefile
  testsuite/multisystem/Makefile
  testsuite/optimal_control/Makefile
  testsuite/periodic_systems/Makefile
  testsuite/photo_electron/Makefile
  testsuite/pseudopotentials/Makefile
  testsuite/tutorials/Makefile
  testsuite/symmetries/Makefile])
AC_CONFIG_FILES([testsuite/oct-run_testsuite.sh], [chmod +x testsuite/oct-run_testsuite.sh])
AC_CONFIG_FILES([testsuite/performance/oct-run_performance_testsuite.sh], [chmod +x testsuite/performance/oct-run_performance_testsuite.sh])
AC_CONFIG_FILES([doc/html/varsRightFrame.html])

AC_CONFIG_SUBDIRS([external_libs/spglib-1.9.9])

dnl if test x"$ac_build_libyaml" != x"skip" ; then
dnl    AC_CONFIG_SUBDIRS([external_libs/yaml-0.1.4])
dnl fi

AC_OUTPUT
