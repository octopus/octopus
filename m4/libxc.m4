AC_DEFUN([ACX_LIBXC], [
acx_libxc_ok=no
FCFLAGS_LIBXC=""
LIBS_LIBXC=""

dnl Check if the library was given in the command line
AC_ARG_WITH(libxc-prefix, [AS_HELP_STRING([--with-libxc-prefix=DIR], [Directory where libxc was installed.])])
case $with_libxc_prefix in
  "") LIBS_LIBXC="-lxc"; FCFLAGS_LIBXC="-I/usr/include" ;;
  *) LIBS_LIBXC="-L$with_libxc_prefix/lib -lxc"; FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_prefix/include" ;;
esac

AC_ARG_WITH(libxc-include, [AS_HELP_STRING([--with-libxc-include=DIR], [Directory where libxc Fortran headers were installed.])])
case $with_libxc_include in
  "") ;;
  *)  FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_libxc_save_LIBS="$LIBS"
acx_libxc_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for libxc])
libxc_fcflags="$FCFLAGS_LIBXC"; libxc_libs="$LIBS_LIBXC"
FCFLAGS="$libxc_fcflags $acx_libxc_save_FCFLAGS"
LIBS="$libxc_libs $acx_libxc_save_LIBS"
AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
  use xc_f90_lib_m
  implicit none
  type(xc_f90_pointer_t) :: info
  integer :: i
  i = xc_f90_info_number(info) + XC_GGA_X_LB
]), [acx_libxc_ok=yes; FCFLAGS_LIBXC="$libxc_fcflags"; LIBS_LIBXC="$libxc_libs"], [])
AC_MSG_RESULT([$acx_libxc_ok ($FCFLAGS_LIBXC $LIBS_LIBXC)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libxc_ok" = xyes; then
  AC_DEFINE(HAVE_LIBXC,1,[Defined if you have the LIBXC library.])
  $1
else
  AC_MSG_ERROR([Could not find required libxc library ( >= v 1.1.0).])
  FCFLAGS_LIBXC=""
  LIBS_LIBXC=""
  $2
fi

AC_SUBST(FCFLAGS_LIBXC)
AC_SUBST(LIBS_LIBXC)
FCFLAGS="$acx_libxc_save_FCFLAGS"
LIBS="$acx_libxc_save_LIBS"
])dnl ACX_LIBXC
