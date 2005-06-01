dnl looks for libjdqz.a
AC_DEFUN([ACX_JDQZ], [
acx_jdqz_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(jdqz, [AC_HELP_STRING([--with-jdqz=DIR], [http://www.math.ruu.nl/people/vorst/jd.html])])
case $with_jdqz in
  yes) ;;
  no | "") acx_jdqz_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_JDQZ="$with_jdqz" ;;
  *) LIBS_JDQZ="-l$with_jdqz" ;;
esac

dnl Backup LIBS 
acx_jdqz_save_LIBS="$LIBS"

dnl First, check if it links
if test $acx_jdqz_ok = no; then
  LIBS="$LIBS_JDQZ $LIBS_LAPACK $LIBS_BLAS $acx_jdqz_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for jdqz library])
  AC_LINK_IFELSE([
    subroutine amul
    end subroutine amul
    subroutine bmul
    end subroutine bmul
    subroutine precon
    end subroutine precon
    program main
    call jdqz
    end program main
], acx_jdqz_ok=yes, [])
  AC_MSG_RESULT($acx_jdqz_ok)
fi

dnl First, check if it links
if test $acx_jdqz_ok = no; then
  LIBS="$LIBS_JDQZ -ljdqz $LIBS_LAPACK $LIBS_BLAS $acx_jdqz_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for jdqz library with -ljdqz])
  AC_LINK_IFELSE([
    subroutine amul
    end subroutine amul
    subroutine bmul
    end subroutine bmul
    subroutine precon
    end subroutine precon
    program main
    call jdqz
    end program main
], [acx_jdqz_ok=yes; LIBS_JDQZ="$LIBS_JDQZ -ljdqz"], [])
  AC_MSG_RESULT($acx_jdqz_ok)
fi

dnl Put the library.
AC_SUBST(LIBS_JDQZ)

dnl Put back LIBS
LIBS="$acx_trlan_save_LIBS"


dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_jdqz_ok" = xyes; then
  AC_DEFINE(HAVE_JDQZ,1,[Defined if you have libjdqz library.])
  $1
else
  if test $acx_jdqz_ok != disable; then
  $2
  fi
fi

])dnl ACX_JDQZ