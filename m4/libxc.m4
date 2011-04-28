AC_DEFUN([ACX_LIBXC], [
acx_libxc_ok=no
FCFLAGS_LIBXC=""
LIBS_LIBXC=""

dnl Check if the library was given in the command line
AC_ARG_WITH(libxc-prefix, [AS_HELP_STRING([--with-libxc-prefix=DIR], [Directory where libxc was installed.])])
case $with_libxc_prefix in
  no ) acx_libxc_ok=disable ;;
  "") LIBS_LIBXC="-lxc"; FCFLAGS_LIBXC="" ;;
  *) LIBS_LIBXC="-L$with_libxc_prefix/lib -lxc"; FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_prefix/include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_libxc_save_LIBS="$LIBS"
acx_libxc_save_FCFLAGS="$FCFLAGS"

dnl The tests
if test "$acx_libxc_ok" = no; then
  AC_MSG_CHECKING([for libxc])
  libxc_fcflags="$FCFLAGS_LIBXC"; libxc_libs="$LIBS_LIBXC"
  FCFLAGS="$libxc_fcflags $acx_libxc_save_FCFLAGS"
  LIBS="$libxc_libs $acx_libxc_save_LIBS"
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    use xc_f90_lib_m

    type(xc_f90_pointer_t) :: info
    integer :: i
    i = xc_f90_info_number(info)
  ]), [acx_libxc_ok=yes; FCFLAGS_LIBXC="$libxc_fcflags"; LIBS_LIBXC="$libxc_libs"], [])
  AC_MSG_RESULT([$acx_libxc_ok ($LIBS_LIBXC)])
fi

AC_SUBST(FCFLAGS_LIBXC)
AC_SUBST(LIBS_LIBXC)
FCFLAGS="$acx_libxc_save_FCFLAGS"
LIBS="$acx_libxc_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libxc_ok" = xyes; then
  AC_DEFINE(HAVE_LIBXC,1,[Defined if you have the LIBXC library.])
  $1
else
  if test $acx_libxc_ok != disable; then
    AC_MSG_ERROR([Could not find required libxc library.])
  fi
  $2
fi
])dnl ACX_LIBXC
