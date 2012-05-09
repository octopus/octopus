dnl looks for libarpack.a
AC_DEFUN([ACX_ARPACK], [
AC_REQUIRE([ACX_BLAS])
acx_arpack_ok=no

dnl We cannot use ARPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
  acx_arpack_ok=noblas
fi

dnl Backup LIBS 
acx_arpack_save_LIBS="$LIBS"


dnl Check if the library was given in the command line
if test $acx_arpack_ok = no; then
  AC_ARG_WITH(arpack, [AS_HELP_STRING([--with-arpack=<lib>], [use ARPACK library <lib> http://forge.scilab.org/index.php/p/arpack-ng/])])
  case $with_arpack in
    yes | "") ;;
    no) acx_arpack_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_ARPACK="$with_arpack" ;;
    *) LIBS_ARPACK="-l$with_arpack" ;;
  esac
fi



dnl First, check LIBS_ARPACK environment variable
if test $acx_arpack_ok = no; then
  LIBS="$LIBS_ARPACK $LIBS_LAPACK $LIBS_BLAS $acx_arpack_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for arpack library])
  AC_LINK_IFELSE([
	    program main
	    call dsaupd
	    end program main], [acx_arpack_ok=yes], [])
  if test $acx_arpack_ok = no; then
    AC_MSG_RESULT([$acx_arpack_ok])
  else
    AC_MSG_RESULT([$acx_arpack_ok ($LIBS_ARPACK)])
  fi
fi

if test $acx_arpack_ok = no; then
  LIBS="$LIBS_ARPACK -larpack $LIBS_LAPACK $LIBS_BLAS $acx_arpack_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for arpack library with -larpack])
  AC_LINK_IFELSE([
    program main
    call dsaupd
    end program main
], [acx_arpack_ok=yes; LIBS_ARPACK="$LIBS_ARPACK -larpack"], [])
  if test $acx_arpack_ok = no; then
    AC_MSG_RESULT([$acx_arpack_ok])
  else
    AC_MSG_RESULT([$acx_arpack_ok ($LIBS_ARPACK)])
  fi
fi


AC_SUBST(LIBS_ARPACK)
LIBS="$acx_arpack_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_arpack_ok" = xyes; then
  AC_DEFINE(HAVE_ARPACK,1,[Defined if you have ARPACK library.])
  $1
else
    AC_MSG_WARN([Could not find ARPACK library. 
               *** Will compile without ARPACK support])
  $2
fi
])dnl ACX_ARPACK
