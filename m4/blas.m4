dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_blas.html
dnl
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_LANG([Fortran])
acx_blas_ok=no

dnl Get fortran linker names of BLAS functions to check for.
dnl if not compiling with Fortran, convert the names
m4_if(_AC_LANG, Fortran, [
sgemm=sgemm
dgemm=dgemm
], [
AC_FC_FUNC(sgemm)
AC_FC_FUNC(dgemm)
])

dnl Check if the library was given in the command line
AC_ARG_WITH(blas, [AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
  yes | "") ;;
  no) acx_blas_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_BLAS="$with_blas" ;;
  *) LIBS_BLAS="-l$with_blas" ;;
esac

dnl Backup LIBS 
acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS_BLAS $LIBS $FLIBS"

dnl First, check LIBS_BLAS environment variable
if test $acx_blas_ok = no; then
  if test "x$LIBS_BLAS" != x; then
    AC_MSG_CHECKING([for $sgemm in $LIBS_BLAS])
    AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [])
    if test $acx_blas_ok = no; then
      AC_MSG_RESULT([$acx_blas_ok])
    else
      AC_MSG_RESULT([$acx_blas_ok ($LIBS_BLAS)])
    fi
  fi
fi

dnl BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
  AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
fi

dnl generic BLAS, in default location such as /usr/lib, as given by Ubuntu libblas-dev package
if test $acx_blas_ok = no; then
   AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; LIBS_BLAS="-lblas"], 
   [], [])
fi

dnl BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(atlas, ATL_xerbla,
    [AC_CHECK_LIB(f77blas, $sgemm,
      [AC_CHECK_LIB(cblas, cblas_dgemm,
       [acx_blas_ok=yes LIBS_BLAS="$LIBS_BLAS -lcblas -lf77blas -latlas"],
      [], [-lf77blas -latlas])],
    [], [-latlas])]
  )
fi

dnl atlas BLAS from debian
if test $acx_blas_ok = no; then
   AC_CHECK_LIB(atlas, $sgemm, [acx_blas_ok=yes; LIBS_BLAS="-L/usr/lib/atlas/ -lblas -latlas"], 
   [], [-L/usr/lib/atlas/ -lblas])
fi

dnl BLAS in version 6 of MKL or higher
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(mkl_ia32, $sgemm, [acx_blas_ok=yes; LIBS_BLAS="$LIBS_BLAS -lmkl_ia32 -lguide -lpthread"],
    [], [-lguide -lpthread])
fi

dnl BLAS in mkl_def
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(mkl_def, $sgemm, [acx_blas_ok=yes; LIBS_BLAS="$LIBS_BLAS -lmkl_def -lguide -lpthread"], 
    [], [-lguide -lpthread])
fi

dnl BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(blas, $sgemm,
    [AC_CHECK_LIB(dgemm, $dgemm,
      [AC_CHECK_LIB(sgemm, $sgemm,
        [acx_blas_ok=yes; LIBS_BLAS="$LIBS_BLAS -lsgemm -ldgemm -lblas"],
      [], [-lblas])],
    [], [-lblas])]
  )
fi

dnl BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
  if test "x$GCC" != xyes; then # only works with Sun CC
    AC_CHECK_LIB(sunmath, acosp,
     [AC_CHECK_LIB(sunperf, $sgemm, 
       [LIBS_BLAS="$LIBS_BLAS -xlic_lib=sunperf -lsunmath" acx_blas_ok=yes],[],[-lsunmath])])
  fi
fi

dnl BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
  unset ac_cv_lib_blas_$sgemm
  AC_CHECK_LIB(blas, $sgemm,
    [AC_CHECK_LIB(essl, $sgemm,
      [acx_blas_ok=yes; LIBS_BLAS="$LIBS_BLAS -lessl -lblas"], [], [-lblas])])
fi

dnl Other libraries
for blas in acml goto mkl cxml dxml scs complib.sgimath blas; do
  if test $acx_blas_ok = no; then
    unset ac_cv_lib_$blas_$sgemm
    AC_CHECK_LIB($blas, $sgemm, [acx_blas_ok=yes; LIBS_BLAS="$LIBS_BLAS -l$blas"], [], [])
  fi
done

dnl Generic BLAS library needing g2c
if test $acx_blas_ok = no; then
  unset ac_cv_lib_blas_$sgemm
  AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; LIBS_BLAS="$LIBS_BLAS -lblas -lg2c"], [], [-lg2c])
fi

AC_SUBST(LIBS_BLAS)
LIBS="$acx_blas_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
  AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.])
  $1
else
  $2
fi
])dnl ACX_BLAS
