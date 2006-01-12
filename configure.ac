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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id$

AC_INIT(Octopus, 2.0rc1, marques@nautilus.fis.uc.pt, octopus)
AC_CONFIG_SRCDIR([src/global.F90])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS(config.h)

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
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_SETVBUF_REVERSED
AC_CHECK_FUNCS([sqrt strcasecmp strdup tcgetpgrp])

# Checks for size of pointer
ACX_POINTER_SIZE

dnl ============================================================================
dnl This part of the configuration file tries to sort out what needs to be built.
dnl In particular, the code can be compile in single precision, mpi, using
dnl complex wave-functions in the ground-state, and in debug mode.

dnl these are the --enable/--disable possible in octopus
AC_ARG_ENABLE(complex, AC_HELP_STRING([--enable-complex], [Ground-state complex wavefunctions]), [COMPLEX_WFNS=1])
AC_ARG_ENABLE(single, AC_HELP_STRING([--enable-single], [octopus in single precision]), [SINGLE_PRECISION=1])
AC_ARG_ENABLE(debug, AC_HELP_STRING([--disable-debug], [debug mode]), [DEBUG=${enableval}])

dnl the MPI has a special treatment
AC_ARG_ENABLE(mpi, AC_HELP_STRING([--enable-mpi(=PATH)], [Paralel version]))
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
AS_IF([test x"$enable_mpi" != x"no"], [MYEXT="${MYEXT}_mpi"])

dnl complex wavefunctions for the ground-state
if test "${COMPLEX_WFNS}"; then
  AC_DEFINE(COMPLEX_WFNS, 1, [octopus compiled for ground-state complex wavefunctions])
  MYEXT="${MYEXT}_cmplx"
  msg="complex wave functions"
else
  msg="real wave functions"
fi

dnl debug mode
if test "${SINGLE_PRECISION}"; then
  AC_DEFINE(SINGLE_PRECISION, 1, [octopus compiled in single precision])
  MYEXT="${MYEXT}_single"
  msg="$msg, single precision"
fi

dnl debug mode
if test x"${DEBUG}" = x"no" ; then
  AC_DEFINE(NDEBUG, yes, [octopus compiled without debug mode])
fi
AC_SUBST(DEBUG)

dnl set the default extension for the programs
EXEEXT="${MYEXT}${EXEEXT}"

dnl we only generate the utilities if single precision, cmplx, and mpi are off
if test x"$enable_mpi" = x"no" -a -z "${COMPLEX_WFNS}" -a -z "${SINGLE_PRECISION}"; then
  UTILITY_PROGRAMS='$(EXTRA_PROGRAMS)'
  msg="$msg + utilities"
else
  UTILITY_PROGRAMS=''
fi
AC_SUBST(UTILITY_PROGRAMS)

dnl tell the user what we will be compiling
echo "**********************************************"
echo "***   octopus${EXEEXT} will be generated using"
echo "***          $msg"
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

dnl check for MPI
if test x"$enable_mpi" != x"no"; then
  ACX_MPI([], AC_MSG_WARN([could not compile an mpi test program]))

	AC_LANG_PUSH(C)
  ACX_METIS([], AC_MSG_WARN([octopus will not be able to use mesh partitioning]))
	AC_LANG_POP(C)
fi

# This is certainly wrong, but I don't know
# how to do it right
LINK="$FC $FCFLAGS $LDFLAGS -o \$@"
AC_SUBST(LINK)
AC_SUBST(FCEXTRALIBS)

AC_CONFIG_FILES([Makefile 
	liboct_parser/Makefile liboct/Makefile libxc/Makefile expokit/Makefile src/Makefile doc/Makefile
	doc/pdf/Makefile doc/html/Makefile share/Makefile share/PP/Makefile share/PP/TM2/Makefile share/PP/HGH/Makefile
	share/recipes/Makefile share/recipes/en/Makefile share/recipes/es/Makefile 
	share/recipes/it/Makefile share/util/Makefile share/samples/Makefile
	build/octopus.spec testsuite/Makefile
])

AC_OUTPUT
