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

AC_INIT(Octopus, 1.3, marques@nautilus.fis.uc.pt, octopus)
AC_CONFIG_SRCDIR([src/global.F90])
AM_INIT_AUTOMAKE

AM_CONFIG_HEADER([liboct/config.h src/config.h])

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

# libtool files
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

dnl ============================================================================
dnl This part of the configuration file tries to sort out what needs to be built.
dnl In particular, the code can be compile in single precision, mpi, using
dnl complex wave-functions in the ground-state, and in debug mode.

dnl these are the --enable possible in octopus
AC_ARG_ENABLE(mpi, AC_HELP_STRING([--enable-mpi(=PATH)], [Paralel version]), MPIDIR=$enableval, MPIDIR=no)
AC_ARG_ENABLE(complex, AC_HELP_STRING([--enable-complex], [Ground-state complex wavefunctions]), [COMPLEX_WFNS=1])
AC_ARG_ENABLE(single, AC_HELP_STRING([--enable-single], [octopus in single precision]), [SINGLE_PRECISION=1])
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug], [octopus in debug mode]), [DEBUG=1])

dnl try to find out what is the default FORTRAN 90 compiler
dnl note that this has to be done after the AC_ARG_ENABLE(mpi...)
if test x"$MPIDIR" != x"no"; then
  AC_PROG_F90([],[MPI])
else
  AC_PROG_F90
fi
AC_LANG_PUSH(Fortran 90)
AC_LANG_PREPROC(Fortran 90)

dnl how Fortran mangles function names
AC_F90_WRAPPERS

dnl ============================================================================
dnl find what is the extension of the programs supposed to be

dnl MPI version
MYEXT=""
AS_IF([test x"$MPIDIR" != x"no"], [MYEXT="${MYEXT}_mpi"])

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
if test "${DEBUG}"; then
  AC_DEFINE(DEBUG, 1, [octopus compiled in debug mode])
  MYEXT="${MYEXT}_debug"
  msg="$msg, debug mode"
fi
AC_SUBST(DEBUG)

dnl set the default extension for the programs
EXEEXT="${MYEXT}${EXEEXT}"

dnl we only generate the utilities if single precision, cmplx, and mpi are off
if test x"$MPIDIR" = x"no" -a -z "${COMPLEX_WFNS}" -a -z "${SINGLE_PRECISION}"; then
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

dnl FFT
ACX_LIB_FFTW

if test -z "${LIB_FFT}"; then
  AC_MSG_WARN([Octopus will be compiled *without* fft support])
else
  AC_DEFINE_UNQUOTED(HAVE_FFT, $fft, [FFT library (fftw2 | fftw3)])
  AC_SUBST(LIB_FFT)
fi

dnl tests for trlan. Default behavior will be *not* using it.   
lib_ok="no"
if test -z "${LIB_TRLAN}"; then
	AC_ARG_WITH(trlan,
  	          AC_HELP_STRING([--with-trlan=DIR], [get it! (www.nersc.gov/research/SIMON/trlan.html)]),
              [TRLANDIR=$withval], [TRLANDIR=""])

  if test "$TRLANDIR" = "no"; then
  	lib_ok="no"
  else
  	save_ldflags="$LDFLAGS"
    AS_IF([test "$TRLANDIR"], [LDFLAGS="${LDFLAGS} -L${TRLANDIR} -L${TRLANDIR}/lib"])

    AC_TRY_COMPILE([],[use trl_info],lib_ok="ok",lib_ok="no")

    LDFLAGS="$save_ldflags"
    LIB_TRLAN="-ltrlan"
    #LIB_TRLAN="${LIB_TRLAN}"
    AS_IF([test "$TRLANDIR"], [LIB_TRLAN="-L${TRLANDIR} -L${TRLANDIR}/lib ${TRLANDIR}/${LIB_TRLAN}"])
  fi
fi

if test "$lib_ok" = "ok"; then
  AC_DEFINE(HAVE_TRLAN, 1, [TRLAN library (www.nersc.gov/research/SIMON/trlan.html)])
  AC_SUBST(LIB_TRLAN)
else
  LIB_TRLAN=""
  AC_MSG_WARN([Will not use trlan])
fi

dnl tests for BLAS
lib_ok="ok"
if test -z "${LIB_BLAS}"; then
	AC_ARG_WITH(blas, 
		AC_HELP_STRING([--with-blas=DIR], [get it! (www.netlib.org/blas)]), 
			[BLASDIR=$withval], [BLASDIR=""])

	if test "$BLASDIR" = "no"; then
		lib_ok="no"
	else
	  save_ldflags="$LDFLAGS"
  	AS_IF([test "$BLASDIR"], [LDFLAGS="${LDFLAGS} -L${BLASDIR} -L${BLASDIR}/lib"])

		ACX_CHECK_FUNC([LIB_BLAS], [ddot])
		ACX_CHECK_LIB([LIB_BLAS], [mkl],  [ddot], [-lmkl], [], [-lmkl])
		ACX_CHECK_LIB([LIB_BLAS], [mkl_def],  [ddot], [-lmkl_def -lguide -lpthread], [], [-lmkl_def -lguide -lpthread])
		ACX_CHECK_LIB([LIB_BLAS], [f77blas], [ddot], [-lf77blas -latlas], [], [-latlas])
		ACX_CHECK_LIB([LIB_BLAS], [blas], [ddot], [-lblas])
		ACX_CHECK_LIB([LIB_BLAS], [dxml], [ddot], [-ldxml])
		ACX_CHECK_LIB([LIB_BLAS], [cxml], [ddot], [-lcxml])
		ACX_CHECK_LIB([LIB_BLAS], [scs], [ddot], [-lscs])
		AS_IF([test -z "${LIB_BLAS}"], AC_MSG_ERROR([couldn't find blas library]))

  	LDFLAGS="$save_ldflags"

		AS_IF([test "$BLASDIR"], [LIB_BLAS="-L${BLASDIR} -L${BLASDIR}/lib ${LIB_BLAS}"])
	fi
fi

if test "$lib_ok" = "ok"; then
	AC_DEFINE(HAVE_BLAS,   1, [Blas   Library (www.netlib.org)])
	AC_SUBST(LIB_BLAS)
fi


lib_ok="ok"
if test -z "${LIB_LAPACK}"; then

	AC_ARG_WITH(lapack, 
		AC_HELP_STRING([--with-lapack=DIR], [get it! (http://www.netlib.org/lapack/)]),
			[LAPACKDIR=$withval], [LAPACKDIR=""])

	if test "$LAPACKDIR" = "no"; then
		lib_ok="no"
	else
  	save_ldflags="$LDFLAGS"
	  AS_IF([test "$LAPACKDIR"], [LDFLAGS="$LDFLAGS -L$LAPACKDIR -L$LAPACKDIR/lib"])

		ACX_CHECK_FUNC([LIB_LAPACK], [dsygv], [], [], [${LIB_BLAS}])
		ACX_CHECK_LIB([LIB_LAPACK], [mkl_lapack], [dsygv], [-lmkl_lapack], [], [${LIB_BLAS}])
		ACX_CHECK_LIB([LIB_LAPACK], [lapack], [dsygv], [-llapack], [], [${LIB_BLAS}])
		ACX_CHECK_LIB([LIB_LAPACK], [dxml], [dsygv], [-ldxml], [], [${LIB_BLAS}])
		ACX_CHECK_LIB([LIB_LAPACK], [cxml], [dsygv], [-lcxml], [], [${LIB_BLAS}])
		ACX_CHECK_LIB([LIB_LAPACK], [scs], [dsygv], [-lscs], [], [${LIB_BLAS}])
		AS_IF([test -z "${LIB_LAPACK}"], [AC_MSG_ERROR([couldn't find lapack library])])

  	LDFLAGS="$save_ldflags"

	  AS_IF([test "$LAPACKDIR"], [LIB_LAPACK="-L${LAPACKDIR} -L${LAPACKDIR}/lib ${LIB_LAPACK}"])
	fi
fi

if test "$lib_ok" = "ok"; then
	AC_DEFINE(HAVE_LAPACK, 1, [Lapack Library (www.netlib.org)])
	AC_SUBST(LIB_LAPACK)
fi

dnl check for GSL
AC_LANG_PUSH(C)
AM_PATH_GSL([1.0])
CFLAGS="$CFLAGS $GSL_CFLAGS"

lib_ok="ok"
if test -z "${LIB_NETCDF}"; then
	AC_ARG_WITH(netcdf, 
		AC_HELP_STRING([--with-netcdf=DIR], [optional (http://www.unidata.ucar.edu/packages/netcdf/)]),
			[NETCDFDIR=$withval], [NETCDFDIR="no"])

	if test "$NETCDFDIR" = "no"; then
		lib_ok="no"
	else
  	save_ldflags="$LDFLAGS"
	  AS_IF([test "$NETCDFDIR"], [LDFLAGS="$LDFLAGS -L$NETCDFDIR -L$NETCDFDIR/lib"])

		ACX_CHECK_FUNC([LIB_NETCDF], [nc_create], [], [], [])
		ACX_CHECK_LIB([LIB_NETCDF], [netcdf], [nc_create], [-lnetcdf], [], [])

  	LDFLAGS="$save_ldflags"

	  AS_IF([test "$NETCDFDIR"], [LIB_NETCDF="-L${NETCDFDIR} -L${NETCDFDIR}/lib ${LIB_NETCDF}"])
	fi
fi
if test "$lib_ok" = "ok"; then
	AC_DEFINE(HAVE_NETCDF, 1, [NetCDF Library (http://www.unidata.ucar.edu/packages/netcdf/)])
	AC_SUBST(LIB_NETCDF)
fi

AC_LANG_POP(C)


if test x"$MPIDIR" != x"no"; then
	if test -z "${LIB_MPI}"; then
    save_ldflags="$LDFLAGS"
   	AS_IF([test "$MPIDIR"], [LDFLAGS="${LDFLAGS} -L${MPIDIR} -L${MPIDIR}/lib"])

  	ACX_CHECK_FUNC([LIB_MPI], [MPI_INIT])
  	ACX_CHECK_LIB([LIB_MPI], [mpi],  [MPI_INIT], [-lmpi])
		AS_IF([test -z "${LIB_MPI}"], AC_MSG_ERROR([couldn't find mpi library]))

  	LDFLAGS="$save_ldflags"
		AS_IF([test "$MPIDIR"], [LIB_MPI="-L${MPIDIR} -L${MPIDIR}/lib ${LIB_MPI}"])
  fi
	
	# let us see if we have a mpi module
	save_ldflags="$LDFLAGS"
	AS_IF([test "$LIB_MPI"], [LDFLAGS="${LDFLAGS} -L${LIB_MPI}"])
	_AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [
use mpi
integer :: ierr
call MPI_Init(ierr)
])], [HAVE_MPI_MOD=1], [HAVE_MPI_MOD=0])

	if test "$HAVE_MPI_MOD" = 1; then
		AC_DEFINE(MPI_MOD, 1, [have mpi module])
	else
		AC_CHECK_HEADERS(mpif.h)
		AC_DEFINE(MPI_H, 1, [have mpi Fortran header file])
	fi

	LDFLAGS="$save_ldflags"
	AC_DEFINE(HAVE_MPI, 1, [MPI library (www.mpich.org)])
  AC_DEFINE(MPI_TD, 1, [soon])
	AC_SUBST(LIB_MPI)
fi

# This is certainly wrong, but I don't know
# how to do it right
LINK="$F90 $F90FLAGS $LDFLAGS -o \$@"
AC_SUBST(LINK)

AC_CONFIG_FILES([Makefile 
	liboct/Makefile src/Makefile doc/Makefile
	share/Makefile share/PP/Makefile share/PP/TM2/Makefile share/PP/HGH/Makefile
	share/recipes/Makefile share/recipes/en/Makefile share/recipes/es/Makefile 
	share/recipes/it/Makefile share/util/Makefile
	build/octopus.spec
])

AC_OUTPUT
