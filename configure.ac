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

AC_INIT(Octopus, 1.1, marques@nautilus.fis.uc.pt, octopus)
AC_CONFIG_SRCDIR([src/global.F90])
AM_INIT_AUTOMAKE

AM_CONFIG_HEADER([liboct/config.h src/config.h])

# Installation prefix by default
AC_PREFIX_DEFAULT($HOME)

# what should I build?
AC_ARG_ENABLE(mpi, AC_HELP_STRING([--enable-mpi(=PATH)], [compile paralel version]), MPIDIR=$enableval, MPIDIR=no)
AC_ARG_ENABLE(complex, AC_HELP_STRING([--enable-complex], [Ground-state complex wavefunctions]), [COMPLEX_WFNS=1])

# set program name
PROGNAME="octopus"
AS_IF([test x"$MPIDIR" != x"no"], [PROGNAME="$PROGNAME-mpi"])
AC_SUBST(PROGNAME)

# complex wavefunctions
if test "${COMPLEX_WFNS}"; then
	AC_DEFINE(COMPLEX_WFNS, 1, [Code compiled for ground-state complex wavefunctions])
	PROGNAME="z$PROGNAME"
  msg="complex wave functions"
else
  msg="real wave functions"
fi


echo "**********************************************"
echo "***   $PROGNAME will be generated"
echo "***       (with $msg)"
echo "**********************************************"

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
AC_FUNC_MALLOC
AC_FUNC_SETVBUF_REVERSED
AC_CHECK_FUNCS([sqrt strcasecmp strdup tcgetpgrp])

# Checks for size of pointer
ACX_POINTER_SIZE

# libtool files
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

# unfortunately we are Fortran 90 guys
if test x"$MPIDIR" != x"no"; then
  AC_PROG_F90([],[MPI])
else
	AC_PROG_F90
fi
AC_LANG_PUSH(Fortran 90)
AC_LANG_PREPROC(Fortran 90)

dnl how Fortran mangles function names
AC_F90_WRAPPERS

dnl tests for FFTW
lib_ok="ok"
if test -z "${LIB_FFTW}"; then
	AC_ARG_WITH(fftw, 
		AC_HELP_STRING([--with-fftw=DIR], [get it! (www.fftw.org)]), 
		[FFTWDIR=$withval], [FFTWDIR=""])

	if test "$FFTWDIR" = "no"; then
		lib_ok="no"
	else
	  save_ldflags="$LDFLAGS"
  	AS_IF([test "$FFTWDIR"], [LDFLAGS="${LDFLAGS} -L${FFTWDIR} -L${FFTWDIR}/lib"])

		ACX_CHECK_FUNC([LIB_FFTW], [fftw3d_f77_create_plan])
		ACX_CHECK_LIB([LIB_FFTW], [fftw],  [fftw3d_f77_create_plan], [-lfftw])
		ACX_CHECK_LIB([LIB_FFTW], [dfftw], [fftw3d_f77_create_plan], [-ldfftw])
		AS_IF([test -z "${LIB_FFTW}"], AC_MSG_ERROR([couldn't find fftw library]))
	
  	ACX_CHECK_FUNC([LIB_RFFTW], [rfftw3d_f77_create_plan], [], [], [${LIB_FFTW}])
   	ACX_CHECK_LIB([LIB_RFFTW], [rfftw], [rfftw3d_f77_create_plan], [-lrfftw], [], ["${LIB_FFTW}"])
   	ACX_CHECK_LIB([LIB_RFFTW], [drfftw], [rfftw3d_f77_create_plan], [-ldrfftw], [], ["${LIB_FFTW}"])
		AS_IF([test -z "${LIB_RFFTW}"], AC_MSG_ERROR([couldn't find rfftw library]))
	
  	LDFLAGS="$save_ldflags"
		LIB_FFTW="${LIB_FFTW} ${LIB_RFFTW}"
		AS_IF([test "$FFTWDIR"], [LIB_FFTW="-L${FFTWDIR} -L${FFTWDIR}/lib ${LIB_FFTW}"])
	fi
fi

if test "$lib_ok" = "ok"; then
	AC_DEFINE(HAVE_FFTW, 1, [FFTW library (www.fftw.org)])
	AC_SUBST(LIB_FFTW)
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
			[NETCDFDIR=$withval], [NETCDF=""])

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
  share/recipes/it/Makefile
	build/octopus.spec
])

AC_OUTPUT
