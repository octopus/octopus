AC_DEFUN([ACX_ETSF_IO], [
AC_REQUIRE([ACX_NETCDF])
acx_etsf_io_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(etsf-io-prefix, [AS_HELP_STRING([--with-etsf-io-prefix=DIR], [Directory where etsf_io was installed.])])
case $with_etsf_io_prefix in
  no ) acx_etsf_io_ok=disabled ;;
  "") if test x$FCFLAGS_ETSF_IO == x; then
    FCFLAGS_ETSF_IO="-I/usr/include"
  fi;;
  *) LIBS_ETSF_IO="-L$with_etsf_io_prefix/lib"; 
     FCFLAGS_ETSF_IO="$ax_cv_f90_modflag$with_etsf_io_prefix/include$fc_type" ;;
esac

AC_ARG_WITH(etsf-io-include, [AS_HELP_STRING([--with-etsf-io-include=DIR], [Directory where etsf_io Fortran headers were installed.])])
case $with_etsf_io_include in
  "") ;;
  *)  FCFLAGS_ETSF_IO="$ax_cv_f90_modflag$with_etsf_io_include" ;;
esac

dnl We cannot use etsf_io if netcdf is not found
if test "x$acx_netcdf_ok" != xyes; then
  acx_etsf_io_ok=disabled
fi

dnl Backup LIBS and FCFLAGS
acx_etsf_io_save_LIBS="$LIBS"
acx_etsf_io_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for etsf_io])
if test "$acx_etsf_io_ok" != disabled; then
  etsf_io_libs="$LIBS_ETSF_IO -letsf_io_utils -letsf_io"
  ABI_PROG_FC()
  for etsf_io_fcflags in "$FCFLAGS_ETSF_IO" "$FCFLAGS_ETSF_IO"/$fc_type ; do
    FCFLAGS="$etsf_io_fcflags $acx_etsf_io_save_FCFLAGS $FCFLAGS_NETCDF"
    LIBS="$etsf_io_libs $acx_etsf_io_save_LIBS $LIBS_NETCDF"
    AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
      use etsf_io
      type(etsf_vars) :: vars
      call etsf_io_vars_free(vars)
    ]), [acx_etsf_io_ok=yes; FCFLAGS_ETSF_IO="$etsf_io_fcflags"; LIBS_ETSF_IO="$etsf_io_libs"], [])
    if test $acx_netcdf_ok == yes; then 
      FCFLAGS_ETSF_IO=$etsf_io_fcflags
      break
    fi
  done
fi
AC_MSG_RESULT([$acx_etsf_io_ok ($FCFLAGS_ETSF_IO $LIBS_ETSF_IO)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_etsf_io_ok" = xyes; then
  AC_DEFINE(HAVE_ETSF_IO,1,[Defined if you have the ETSF_IO library.])
  $1
else
  AC_MSG_WARN([Could not find etsf_io library. 
           *** Will compile without ETSF I/O support])
  FCFLAGS_ETSF_IO=""
  LIBS_ETSF_IO=""
  $2
fi

AC_SUBST(FCFLAGS_ETSF_IO)
AC_SUBST(LIBS_ETSF_IO)
FCFLAGS="$acx_etsf_io_save_FCFLAGS"
LIBS="$acx_etsf_io_save_LIBS"
])dnl ACX_ETSF_IO
