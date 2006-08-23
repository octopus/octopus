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
AC_INIT([Octopus],[2.1.0rc1],[octopus-devel@tddft.org],[octopus])
AC_CONFIG_SRCDIR([src/global.F90])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

# Installation prefix by default
AC_PREFIX_DEFAULT($HOME)

# who am i
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_YACC

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h float.h locale.h malloc.h stddef.h stdlib.h string.h strings.h sys/ioctl.h sys/time.h termios.h unistd.h signal.h])
AC_FUNC_ALLOCA

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_GETPGRP
AC_PROG_GCC_TRADITIONAL
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_STAT
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([getcwd gettimeofday mkdir pow sqrt strcasecmp strdup \
strerror strncasecmp strndup uname tcgetpgrp scandir alphasort \
sigaction malloc realloc opendir readdir closedir strchr strtod])

# Checks for size of pointer
ACX_POINTER_SIZE

dnl ============================================================================
dnl This part of the configuration file tries to sort out what needs to be built.
dnl In particular, the code can be compile in single precision, mpi, using
dnl complex wave-functions in the ground-state, and in debug mode.

dnl these are the --enable/--disable possible in octopus
AC_ARG_ENABLE(single, AS_HELP_STRING([--enable-single], [octopus in single precision]), [SINGLE_PRECISION=1])
AC_ARG_ENABLE(debug, AS_HELP_STRING([--disable-debug], [debug mode]), [DEBUG=${enableval}])


dnl the MPI has a special treatment
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

dnl build/install python bits unless --disable-python is used
dnl set temporarily to default to not enabling.
AC_ARG_ENABLE([python],
	      AS_HELP_STRING([--disable-python], [don't build/install python components]),
	      [ac_cv_build_python=$enableval],
	      [ac_cv_build_python=yes])
if test $ac_cv_build_python = yes; then
  AC_CONFIG_SUBDIRS(python)
fi
AM_CONDITIONAL([ENABLE_PYTHON], [test $ac_cv_build_python = yes])


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

dnl some functions we may need
AC_CHECK_FUNCS([flush])

dnl check wether we have or not a compiler that allows for very long lines...
ACX_LONG_FORTRAN_LINES

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

dnl we only generate the utilities if single precision and mpi are off
if test x"$enable_mpi" = x"no" -a -z "${SINGLE_PRECISION}"; then
  UTILITY_PROGRAMS='$(EXTRA_PROGRAMS)'
  msg="+ utilities"
else
  UTILITY_PROGRAMS=''
  msg=""
fi
AC_SUBST(UTILITY_PROGRAMS)

dnl tell the user what we will be compiling
echo "**********************************************"
echo "***   octopus${EXEEXT} $msg will be generated "
echo "**********************************************"

dnl ===================================================================
dnl now we search for libraries

dnl check for BLAS
ACX_BLAS([], AC_MSG_ERROR([could not find required blas library]))

dnl check for LAPACK
ACX_LAPACK([], AC_MSG_ERROR([could not find required lapack library]))

dnl check for GSL
AC_LANG_PUSH(C)
LIBS="${acx_save_libs}"

AM_PATH_GSL([1.0], 
  [CFLAGS="$CFLAGS $GSL_CFLAGS"], 
  AC_MSG_ERROR([could not find required gsl library]))

dnl check for GD library
ACX_GDLIB

LIBS="${FCLIBS}"
AC_LANG_POP(C)

dnl check for FFT
ACX_FFT

dnl check for TRLAN
ACX_TRLAN

dnl check for NETCDF
ACX_NETCDF

dnl check for JDQZ
dnl ACX_JDQZ([], AC_MSG_WARN([libjdqz.a library not found. Will not use it]))

dnl check for ARPACK
ACX_ARPACK([], AC_MSG_WARN([ARPACK library not found. Will not use it]))

dnl check for SPARSKIT
ACX_SPARSKIT([], AC_MSG_WARN([SPARSKIT library not found. Will not use it]))

dnl check for METIS
if test x"$enable_mpi" != x"no"; then
  ACX_METIS
fi
AM_CONDITIONAL(COMPILE_METIS, test x${HAVE_METIS} = x1)

dnl check for MPI
if test x"$enable_mpi" != x"no"; then
  ACX_MPI([], AC_MSG_WARN([could not compile an mpi test program]))
fi

# This is certainly wrong, but I don't know
# how to do it right
LINK="$FC $FCFLAGS $LDFLAGS -o \$@"
AC_SUBST(LINK)
AC_SUBST(FCEXTRALIBS)
AC_SUBST(FCDEBUGFLAGS)

AC_CONFIG_FILES([Makefile 
  external_libs/Makefile external_libs/expokit/Makefile external_libs/metis-4.0/Makefile  external_libs/qshep/Makefile
	liboct_parser/Makefile liboct/Makefile libxc/Makefile src/Makefile doc/Makefile
	doc/pdf/Makefile doc/html/Makefile share/Makefile share/PP/Makefile share/PP/TM2/Makefile share/PP/HGH/Makefile
	share/recipes/Makefile share/recipes/en/Makefile share/recipes/es/Makefile 
	share/recipes/it/Makefile share/util/Makefile share/samples/Makefile
	build/octopus.spec testsuite/Makefile testsuite/finite_systems_1d/Makefile 
	testsuite/finite_systems_2d/Makefile testsuite/finite_systems_3d/Makefile 
	testsuite/periodic_systems_1d/Makefile
])

AC_OUTPUT
