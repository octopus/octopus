AC_INIT(Octopus, 0.9, marques@nautilus.fis.uc.pt, octopus)
AC_CONFIG_SRCDIR([src/global.F90])
AM_INIT_AUTOMAKE

AM_CONFIG_HEADER([liboct/config.h src/config.h])

# Installation prefix by default
AC_PREFIX_DEFAULT($HOME)

# who am i
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_YACC
ACX_PROG_F90

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

dnl FORTRAN 90 checks
ACX_F90_FUNC_MANGLE
AC_LANG_PUSH(Fortran 90)

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
  	if test "$FFTWDIR"; then
    	LDFLAGS="${LDFLAGS} -L${FFTWDIR} -L${FFTWDIR}/lib"
	  fi

  	if test -z "${LIB_FFTW}"; then
  		ACX_CHECK_FUNC_F90(fftw3d_f77_create_plan, [LIB_FFTW=" "])
	  fi
  	if test -z "${LIB_FFTW}"; then
    	ACX_CHECK_LIB_F90(fftw, fftw3d_f77_create_plan, [LIB_FFTW="-lfftw"])
	  fi
  	if test -z "${LIB_FFTW}"; then
    	ACX_CHECK_LIB_F90(dfftw, fftw3d_f77_create_plan, [LIB_FFTW="-ldfftw"])
	  fi
  	if test -z "${LIB_FFTW}"; then
    	AC_MSG_ERROR([couldn't find fftw library])
		fi
	
		dnl and now we check for rfftw
	  if test -z "${LIB_RFFTW}"; then
  		ACX_CHECK_FUNC_F90(rfftw3d_f77_create_plan, [LIB_RFFTW=" "], [], "${LIB_FFTW}")
	  fi
  	if test -z "${LIB_RFFTW}"; then
    	ACX_CHECK_LIB_F90(rfftw, rfftw3d_f77_create_plan, [LIB_RFFTW="-lrfftw"], [], "${LIB_FFTW}")
	  fi
  	if test -z "${LIB_RFFTW}"; then
    	ACX_CHECK_LIB_F90(drfftw, rfftw3d_f77_create_plan, [LIB_RFFTW="-ldrfftw"], [], "${LIB_FFTW}")
	  fi
  	if test -z "${LIB_RFFTW}"; then
    	AC_MSG_ERROR([couldn't find rfftw library])
		fi
	
  	LDFLAGS="$save_ldflags"

	  if test "$FFTWDIR"; then
			LIB_FFTW="-L${FFTWDIR} -L${FFTWDIR}/lib ${LIB_FFTW} ${LIB_RFFTW}"
		fi
	fi
fi

if test "$lib_ok" = "ok"; then
	AC_DEFINE(HAVE_FFTW, 1, [FFTW library (www.fftw.org)])
	AC_SUBST(LIB_FFTW)
fi


dnl tests for BLAS/LAPACK
lib_ok="ok"
if test -z "${LIB_LAPACK}"; then

	AC_ARG_WITH(lapack, 
		AC_HELP_STRING([--with-lapack=DIR], [get it! (http://www.netlib.org/lapack/)]),
		[LAPACKDIR=$withval], [LAPACKDIR=""])

	if test "$LAPACKDIR" = "no"; then
		lib_ok="no"
	else
  	save_ldflags="$LDFLAGS"
	  if test "$LAPACKDIR"; then
  	  LDFLAGS="$LDFLAGS -L$LAPACKDIR -L$LAPACKDIR/lib"
	  fi

  	if test -z "${LIB_LAPACK}"; then
	  	ACX_CHECK_FUNC_F90(dsygv, [LIB_LAPACK=" "])
  	fi
	  if test -z "${LIB_LAPACK}"; then
  	  ACX_CHECK_LIB_F90(lapack, dsygv, [LIB_LAPACK="-llapack -lblas"], [], "-lblas")
	  fi
	  if test -z "${LIB_LAPACK}"; then   dnl for the alphas
  	  ACX_CHECK_LIB_F90(cxml, dsygv, [LIB_LAPACK="-lcxml"])
	  fi
  	if test -z "${LIB_LAPACK}"; then
    	AC_MSG_ERROR([couldn't find lapack/blas library])
		fi

  	LDFLAGS="$save_ldflags"

	  if test "$LAPACKDIR"; then
			LIB_LAPACK="-L${LAPACKDIR} -L${LAPACKDIR}/lib ${LIB_LAPACK}"
		fi
	fi
fi

if test "$lib_ok" = "ok"; then
	AC_DEFINE(HAVE_LAPACK, 1, [Lapack Library (www.netlib.org)])
	AC_DEFINE(HAVE_BLAS,   1, [Blas   Library (www.netlib.org)])
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
  	if test "$GSLDIR"; then
    	LDFLAGS="${LDFLAGS} -L${GSLDIR}/lib"
			CPPFLAGS="${CPPFLAGS} -I${GSLDIR}/include"
  	fi

	  if test -z "${LIB_GSL}"; then
  		AC_CHECK_FUNC(gsl_complex_sqrt, [LIB_GSL=" "])
	  fi
  	if test -z "${LIB_GSL}"; then
    	AC_CHECK_LIB(gsl, gsl_complex_sqrt, [LIB_GSL="-lgsl -lgslcblas"], [], "-lgslcblas")
	  fi
  	if test -z "{$LIB_GSL}"; then
    	AC_MSG_ERROR([couldn't find gsl library])
		fi

  	LDFLAGS="$save_ldflags"

	  if test "$GSLDIR"; then
  		LIB_GSL="-L${GSLDIR}/lib ${LIB_GSL}"
		fi
	fi

	AC_LANG_POP(C)
fi

if test "$lib_ok" = "ok"; then
 	AC_DEFINE(HAVE_GSL, 1, [Gnu Scientific Library (http://sources.redhat.com/gsl/)])
	AC_SUBST(LIB_GSL)
fi

dnl we can work in 1 and 3 dimensions
AC_ARG_ENABLE(1d,
	AC_HELP_STRING([--enable-1d], [Let us work in 1D]), ONE_D=1)
AC_ARG_ENABLE(3d, 
	AC_HELP_STRING([--enable-3d], [or, if you prefer, in 3D]), [THREE_D=1])

if test "${ONE_D}"; then
	AC_DEFINE(ONE_D,   1, [Code compiled for 1D systems])
	PROGNAME="octopus_1d"
else
	AC_DEFINE(THREE_D, 1, [Code compiled for 3D systmes])
	PROGNAME="octopus"
fi
AC_SUBST(PROGNAME)

# This is certainly wrong, but I don't know
# how to do it right
LINK="$F90 $F90FLAGS -o $PROGNAME"
AC_SUBST(LINK)

AC_CONFIG_FILES([Makefile liboct/Makefile src/Makefile doc/Makefile])

AC_OUTPUT