AC_INIT(Octopus, 0.9, marques@nautilus.fis.uc.pt, octopus)
AC_CONFIG_SRCDIR([src/global.F90])
AM_INIT_AUTOMAKE

AM_CONFIG_HEADER([liboct/config.h src/config.h])

# Installation prefix by default
AC_PREFIX_DEFAULT($HOME)

# who am i
AC_CANONICAL_HOST

AC_ARG_ENABLE(mpi, AC_HELP_STRING([--enable-mpi(=PATH)], [compile paralel version]), 
	MPIDIR=$enableval, MPIDIR=no)
AC_ARG_ENABLE(1d, AC_HELP_STRING([--enable-1d], [Let us work in 1D]), ONE_D=1)
AC_ARG_ENABLE(3d, AC_HELP_STRING([--enable-3d], [or, if you prefer, in 3D]), [THREE_D=1])

# set program name
PROGNAME="octopus"
AS_IF([test x"$MPIDIR" != x"no"], [PROGNAME="$PROGNAME-mpi"])
if test "${ONE_D}"; then
	AC_DEFINE(ONE_D,   1, [Code compiled for 1D systems])
	PROGNAME="$PROGNAME-1d"
else
	AC_DEFINE(THREE_D, 1, [Code compiled for 3D systmes])
fi
AC_SUBST(PROGNAME)

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
AC_CHECK_FUNCS([sqrt strcasecmp strdup])

# Checks for size of pointer
ACX_POINTER_SIZE

# libtool files
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

# unfortunately we are Fortran 90 guys
if test x"$MPIDIR" != x"no"; then
	AC_CHECK_TOOLS(F90, [mpif90])
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
		AS_IF([test "$FFTWDIR"], [LIB_FFTW="-L${FFTWDIR} -L${FFTWDIR}/lib ${LIB_FFTW} ${LIB_RFFTW}"])
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
		ACX_CHECK_LIB([LIB_BLAS], [f77blas], [ddot], [-lf77blas -latlas], [], [-latlas])
		ACX_CHECK_LIB([LIB_BLAS], [blas], [ddot], [-lblas])
		ACX_CHECK_LIB([LIB_BLAS], [dxml], [ddot], [])
		ACX_CHECK_LIB([LIB_BLAS], [cxml], [ddot], [])
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
		ACX_CHECK_LIB([LIB_LAPACK], [lapack], [dsygv], [-llapack], [], [${LIB_BLAS}])
		ACX_CHECK_LIB([LIB_LAPACK], [dxml], [dsygv], [-llapack], [], [${LIB_BLAS}])
		ACX_CHECK_LIB([LIB_LAPACK], [cxml], [dsygv], [-llapack], [], [${LIB_BLAS}])
		AS_IF([test -z "${LIB_LAPACK}"], [AC_MSG_ERROR([couldn't find lapack library])])

  	LDFLAGS="$save_ldflags"

	  AS_IF([test "$LAPACKDIR"], [LIB_LAPACK="-L${LAPACKDIR} -L${LAPACKDIR}/lib ${LIB_LAPACK}"])
	fi
fi

if test "$lib_ok" = "ok"; then
	AC_DEFINE(HAVE_LAPACK, 1, [Lapack Library (www.netlib.org)])
	AC_SUBST(LIB_LAPACK)
fi


dnl now we check for gsl
lib_ok="ok"
if test -z "${LIB_GSL}"; then
	AC_LANG_PUSH(C) 	dnl gsl is to be used from C

	AC_ARG_WITH(gsl,
		AC_HELP_STRING([--with-gsl=DIR], [get it! (http://sources.redhat.com/gsl/)]),
			[GSLDIR=$withval], [GSLDIR=""])

	if test "$GSLDIR" = "no"; then
		lib_ok="no"
	else
	  save_ldflags="$LDFLAGS"
		AS_IF([test "$GSLDIR"], [
    	LDFLAGS="${LDFLAGS} -L${GSLDIR}/lib"
			CPPFLAGS="${CPPFLAGS} -I${GSLDIR}/include"
		])

		ACX_CHECK_FUNC([LIB_GSL], [gsl_complex_sqrt], [-lm], [], [-lm])
	 	ACX_CHECK_LIB([LIB_GSL], [gsl], [gsl_complex_sqrt], [-lgsl -lgslcblas -lm], [], [-lgslcblas -lm])
  	AS_IF([test -z "${LIB_GSL}"], [AC_MSG_ERROR([couldn't find gsl library])])

  	LDFLAGS="$save_ldflags"
	  AS_IF([test "$GSLDIR"], [LIB_GSL="-L${GSLDIR}/lib ${LIB_GSL}"])
	fi

	AC_LANG_POP(C)
fi

if test "$lib_ok" = "ok"; then
 	AC_DEFINE(HAVE_GSL, 1, [Gnu Scientific Library (http://sources.redhat.com/gsl/)])
	AC_SUBST(LIB_GSL)
fi

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
	AC_SUBST(LIB_MPI)
fi

# This is certainly wrong, but I don't know
# how to do it right
LINK="$F90 $F90FLAGS $LDFLAGS -o \$@"
AC_SUBST(LINK)

AC_CONFIG_FILES([Makefile liboct/Makefile src/Makefile doc/Makefile])

AC_OUTPUT