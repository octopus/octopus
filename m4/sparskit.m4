dnl looks for libskit.a
AC_DEFUN([ACX_SPARSKIT], [
acx_sparskit_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(sparskit, [AS_HELP_STRING([--with-sparskit=DIR],[http://www-users.cs.umn.edu/~saad/software/])])
case $with_sparskit in
  yes | "") ;;
  no ) acx_sparskit_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_SPARSKIT="$with_sparskit" ;;
  *) LIBS_SPARSKIT="-l$with_sparskit" ;;
esac

dnl Backup LIBS 
acx_sparskit_save_LIBS="$LIBS"

dnl First, check if it links
if test $acx_sparskit_ok = no; then
  LIBS="$LIBS_SPARSKIT $LIBS_LAPACK $LIBS_BLAS $acx_sparskit_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for SPARSKIT library])
  AC_LINK_IFELSE([
    subroutine distdot
    end subroutine distdot
    program main
    call bcgstab
    end program main
], acx_sparskit_ok=yes, [])
  if test $acx_sparskit_ok = no; then
    AC_MSG_RESULT([$acx_sparskit_ok])
  else
    AC_MSG_RESULT([$acx_sparskit_ok($LIBS_SPARSKIT)])
  fi
fi

dnl ... check if it links ...
if test $acx_sparskit_ok = no; then
  LIBS="$LIBS_SPARSKIT /usr/lib/libskit.a $LIBS_LAPACK $LIBS_BLAS $acx_sparskit_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for SPARSKIT library with libskit.a])
  AC_LINK_IFELSE([
    subroutine distdot
    end subroutine distdot
    program main
    call bcgstab
    end program main
], [acx_sparskit_ok=yes; LIBS_SPARSKIT="$LIBS_SPARSKIT /usr/lib/libskit.a"], [])
  if test $acx_sparskit_ok = no; then
    AC_MSG_RESULT([$acx_sparskit_ok])
  else
    AC_MSG_RESULT([$acx_sparskit_ok ($LIBS_SPARSKIT)])
  fi
fi

dnl Put the library.
AC_SUBST(LIBS_SPARSKIT)

dnl Put back LIBS
LIBS="$acx_sparskit_save_LIBS"


dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_sparskit_ok" = xyes; then
  AC_DEFINE(HAVE_SPARSKIT,1,[Defined if you have SPARSKIT library.])
  $1
else
  if test $acx_sparskit_ok != disable; then
    AC_MSG_WARN([Could not find SPARSKIT library. 
                *** Will compile without SPARSKIT support])
  fi
  $2
fi

])dnl ACX_SPARSKIT
