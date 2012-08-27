## Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id$

AC_PREREQ(2.59)
AC_INIT([Octopus],[nocturnus],[octopus-devel@tddft.org],[octopus])
AC_CONFIG_SRCDIR([src/basic/global.F90])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

# Installation prefix by default
AC_PREFIX_DEFAULT([/opt/etsf])

# who am i
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_YACC

dnl check whether mpi is enabled
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--enable-mpi(=PATH)], [Parallel version]))
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

if test x"$enable_mpi" == x"yes"; then
  octopus_default_cc=mpicc
else
  octopus_default_cc=gcc
fi

AC_PROG_CC($octopus_default_cc)
AM_PROG_AS

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

# Checks for library functions.
AC_PROG_GCC_TRADITIONAL
AC_FUNC_STAT
AC_CHECK_FUNCS([gettimeofday \
strndup tcgetpgrp scandir alphasort \
sigaction readdir closedir strchr strtod \
perror strcasestr sbrk getopt_long uname ioctl])

ACX_PAPI
ACX_ARCH

dnl ============================================================================
dnl This part of the configuration file tries to sort out what needs to be built.
dnl In particular, the code can be compiled in single-precision, in debug mode,
dnl and with OpenMP threads.

dnl these are the --enable/--disable possible in octopus
AC_ARG_ENABLE(single, AS_HELP_STRING([--enable-single], [octopus in single-precision]))
AC_ARG_ENABLE(debug, AS_HELP_STRING([--disable-debug], [debug mode]), [DEBUG=${enableval}])
AC_ARG_ENABLE(openmp, AS_HELP_STRING([--enable-openmp], [Multi-threaded parallel version using OpenMP]))
AC_ARG_ENABLE(opencl, AS_HELP_STRING([--enable-opencl], [OpenCL support (experimental)]))
AC_ARG_ENABLE(blacs, AS_HELP_STRING([--with-blacs], [Blacs support (experimental)]))
AC_ARG_ENABLE(scalapack, AS_HELP_STRING([--with-scalapack], [ScaLAPACK support (experimental), add before --with-blacs]))
AC_ARG_ENABLE(libfm, AS_HELP_STRING([--with-libfm], [Libfm support (experimental)]))
AC_ARG_ENABLE(libfm, AS_HELP_STRING([--with-pfft], [PFFT support (experimental)]))

if test -n "${enable_single}" -a x"${enable_single}" != x"no"; then
  SINGLE_PRECISION=1
fi

if test -n "${enable_openmp}" -a x"${enable_openmp}" != x"no"; then
  HAVE_OPENMP=1
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
  AC_CHECK_FUNC(clAmdBlasGetVersion, [AC_DEFINE(HAVE_CLAMDBLAS, 1, [defined if clAmdBlas is available])])
  AC_CHECK_FUNC(clAmdFftGetVersion, [AC_DEFINE(HAVE_CLAMDFFT, 1, [defined if clAmdFft is available])])
fi

AM_CONDITIONAL(COMPILE_OPENCL, test x"${enable_opencl}" == x"yes")

dnl workaround for xlf compilers
ac_enable_f90_forall=yes
AC_ARG_ENABLE(f90-forall, AS_HELP_STRING([--disable-f90-forall], [Workaround for buggy xlf compilers]), 
        [ac_enable_f90_forall=${enableval}])
if test x"${ac_enable_f90_forall}" = x"no"; then
  F90_FORALL=no
else
  F90_FORALL=yes
fi
AC_SUBST(F90_FORALL)

dnl check for METIS
if test x"$enable_mpi" != x"no"; then
  ACX_METIS
fi
AM_CONDITIONAL(COMPILE_METIS, test x${HAVE_METIS} = x1)
  
dnl check for NEWUOA
ACX_NEWUOA
AM_CONDITIONAL(COMPILE_NEWUOA, test x${HAVE_NEWUOA} = x1)

dnl check for external Zoltan
ACX_ZOLTAN
AM_CONDITIONAL(COMPILE_ZOLTAN, test x"$compile_zoltan" == x"yes")

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
dnl this is a hack, and should be changed in the future
AC_LANG_PREPROC

dnl how Fortran mangles function names
AC_FC_WRAPPERS
acx_save_libs="${LIBS}"
LIBS="${FCEXTRALIBS} $FCLIBS"

ACX_FC_INTEGER_SIZE
ACX_CC_FORTRAN_INT

dnl Sometimes Fortran may use the flush function, sometimes not.
ACX_FORTRAN_CHECK_FUNC(flush)

dnl Some compilers define the sizeof intrinsic
ACX_FC_SIZEOF

ACX_FORTRAN_LOC

dnl check whether we have or not a compiler that allows for very long lines...
ACX_LONG_FORTRAN_LINES

dnl check whether the Fortran 90 compiler accepts line markers cast by the preprocessor.
ACX_F90_ACCEPTS_LINE_NUMBERS

dnl check whether the Fortran compiler supports command-line arguments
ACX_FC_COMMAND_LINE_ARGUMENTS

dnl check whether the Fortran compiler supports Fortran 2003's c_ptr type
ACX_FC_F2003_C_PTR

dnl check flag to include module path
AX_F90_MODULE_FLAG
F90_MODULE_FLAG=$ax_cv_f90_modflag
AC_SUBST([F90_MODULE_FLAG])

if test "${HAVE_OPENMP}"; then
  AC_DEFINE(HAVE_OPENMP, 1, [enable OpenMP support])
  AX_OPENMP(,AC_MSG_ERROR([Fortran compiler does not support OpenMP]))
  FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
fi

dnl ============================================================================
dnl find what is the extension of the programs supposed to be
MYEXT=""

dnl mpi
AS_IF([test x"$enable_mpi" != x"no"], [MYEXT="${MYEXT}_mpi"])

dnl debug mode
if test "${SINGLE_PRECISION}"; then
  AC_DEFINE(SINGLE_PRECISION, 1, [octopus compiled in single precision])
  MYEXT="${MYEXT}_single"
fi

dnl debug mode
if test x"${DEBUG}" = x"no" ; then
  AC_DEFINE(NDEBUG, yes, [octopus compiled without debug mode])
fi
AC_SUBST(DEBUG)

dnl set the default extension for the programs
EXEEXT="${MYEXT}${EXEEXT}"

dnl can force build of utilties even with mpi
AC_ARG_ENABLE(utils, AS_HELP_STRING([--enable-utils], [Build utilities (may sometimes work even with mpi)]))

dnl we only generate the utilities if single precision and MPI are off
if test \( x"$enable_mpi" = x"no" -a -z "${SINGLE_PRECISION}" \) -o \( x"$enable_utils" = x"yes" \); then
  UTILITY_PROGRAMS='$(EXTRA_PROGRAMS)'
  msg="+ utilities"
else
  UTILITY_PROGRAMS=''
  msg=""
fi
AC_SUBST(UTILITY_PROGRAMS)

AC_ARG_WITH(max-dim, [  --with-max-dim=NDIM     maximum number of dimensions Octopus can use;
                          the default is 3 (other values are experimental), must be >= 3],
  [maxdim=$withval], [maxdim=3])

if test $maxdim -lt 3; then
  AC_MSG_ERROR([max-dim must be at least 3])
fi

AC_DEFINE_UNQUOTED(MAX_DIM, $maxdim, [the maximum dimension of the space])

dnl tell the user what we will be compiling
echo "**********************************************"
echo "***   octopus${EXEEXT} $msg will be generated "
echo "**********************************************"

dnl ===================================================================
dnl now we search for libraries

ACX_LIBXC

dnl check for BLAS
ACX_BLAS([], AC_MSG_ERROR([could not find required BLAS library]))

dnl test zdotc in BLAS
ACX_ZDOTC

dnl check for LAPACK
ACX_LAPACK([], AC_MSG_ERROR([could not find required LAPACK library]))

dnl check for GSL
AC_LANG_PUSH(C)
LIBS="${acx_save_libs}"

AX_PATH_GSL([1],[9], 
  [CFLAGS="$CFLAGS $GSL_CFLAGS"], 
  AC_MSG_ERROR([could not find required GSL library]))

dnl check for GD library
ACX_GDLIB

LIBS="${FCLIBS}"
AC_LANG_POP(C)

dnl check for FFT
ACX_FFT

dnl check for NFFT
AC_LANG_PUSH(C)
ACX_NFFT
AC_LANG_POP(C)

dnl check for NETCDF
ACX_NETCDF

dnl check for ETSF_IO
ACX_ETSF_IO

dnl check for BerkeleyGW
ACX_BERKELEYGW

dnl check for SPARSKIT
ACX_SPARSKIT
  
dnl check for Fortran MPI support
if test x"$enable_mpi" != x"no"; then
  ACX_MPI([], AC_MSG_ERROR([could not compile an MPI test program]))
  ACX_MPI_FC_MODULE
  ACX_MPI2
fi

dnl check for SCALAPACK/BLACS
ACX_BLACS
ACX_SCALAPACK

dnl check for ARPACK
ACX_ARPACK

dnl check for Libfm
ACX_LIBFM

dnl check for PFFT
ACX_PFFT

dnl now I run configure for the internal libraries
AC_CONFIG_SUBDIRS([external_libs/spglibs])

# This is certainly wrong, but I don't know
# how to do it right
LINK="$FC $FCFLAGS $LDFLAGS -o \$@"
AC_SUBST(LINK)
AC_SUBST(FCEXTRALIBS)

ACX_OCTOPUS_COMPILATION_INFO

AC_CONFIG_FILES([Makefile
  external_libs/Makefile
  external_libs/metis-4.0/Makefile
  external_libs/qshep/Makefile
  external_libs/spglib/Makefile
  external_libs/spglib/src/Makefile
  external_libs/newuoa/Makefile
  external_libs/zoltan/Makefile
  external_libs/slatec/Makefile
  external_libs/fortrancl/Makefile
  external_libs/bpdn/Makefile
  liboct_parser/Makefile
  src/Makefile
  src/include/Makefile
  src/basic/Makefile
  src/math/Makefile
  src/species/Makefile
  src/ions/Makefile
  src/grid/Makefile
  src/poisson/Makefile
  src/states/Makefile
  src/xc/Makefile
  src/hamiltonian/Makefile
  src/system/Makefile
  src/scf/Makefile
  src/td/Makefile
  src/opt_control/Makefile
  src/sternheimer/Makefile
  src/main/Makefile
  src/utils/Makefile
  perl/Makefile
  doc/Makefile
  doc/pdf/Makefile
  doc/html/Makefile
  doc/html/icons/Makefile
  doc/html/js/Makefile
  doc/html/vars/Makefile
  share/Makefile
  share/PP/Makefile
  share/PP/PSF/Makefile
  share/PP/HGH/Makefile
  share/PP/UPF/Makefile
  share/recipes/Makefile
  share/recipes/en/Makefile
  share/recipes/es/Makefile
  share/recipes/it/Makefile
  share/util/Makefile
  share/opencl/Makefile
  build/Makefile
  build/octopus.spec
  testsuite/Makefile
  testsuite/components/Makefile 
  testsuite/finite_systems_1d/Makefile 
  testsuite/finite_systems_2d/Makefile
  testsuite/finite_systems_3d/Makefile
  testsuite/functionals/Makefile
  testsuite/real_time/Makefile
  testsuite/linear_response/Makefile
  testsuite/optimal_control/Makefile
  testsuite/open_systems/Makefile
  testsuite/periodic_systems/Makefile
  testsuite/photo_electron/Makefile])
AC_CONFIG_FILES([testsuite/oct-run_testsuite], [chmod +x testsuite/oct-run_testsuite])

AC_OUTPUT
