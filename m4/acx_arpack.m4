dnl looks for libarpack.a
AC_DEFUN([ACX_ARPACK], [
acx_arpack_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(arpack, [AS_HELP_STRING([--with-arpack=DIR], [http://www.caam.rice.edu/software/ARPACK])])
case $with_arpack in
  yes) ;;
  no | "") acx_arpack_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_ARPACK="$with_arpack" ;;
  *) LIBS_ARPACK="-l$with_arpack" ;;
esac

dnl Backup LIBS 
acx_arpack_save_LIBS="$LIBS"

dnl First, check if it links
if test $acx_arpack_ok = no; then
  LIBS="$LIBS_ARPACK $LIBS_LAPACK $LIBS_BLAS $acx_arpack_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for arpack library])
  AC_LINK_IFELSE([
    program main
    call dsaupd
    end program main
], acx_arpack_ok=yes, [])
  if test $acx_arpack_ok = no; then
    AC_MSG_RESULT([$acx_arpack_ok])
  else
    AC_MSG_RESULT([$acx_arpack_ok ($LIBS_ARPACK)])
  fi
fi

dnl First, check if it links
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

dnl Put the library.
AC_SUBST(LIBS_ARPACK)

dnl Put back LIBS
LIBS="$acx_arpack_save_LIBS"


dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_arpack_ok" = xyes; then
  AC_DEFINE(HAVE_ARPACK,1,[Defined if you have libarpack library.])
  $1
else
  $2
fi

])dnl ACX_ARPACK

