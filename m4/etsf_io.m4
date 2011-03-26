AC_DEFUN([ACX_ETSF_IO], [
AC_REQUIRE([ACX_NETCDF])
acx_etsf_io_ok=no
FCFLAGS_ETSF_IO=""
LIBS_ETSF_IO=""

dnl We cannot use etsf_io if netcdf is not found
if test "x$acx_netcdf_ok" != xyes; then
  acx_etsf_io_ok=nonetcdf
  FCFLAGS_ETSF_IO=""
  LIBS_ETSF_IO=""
fi

dnl Check if the library was given in the command line
AC_ARG_WITH(etsf_io, [AS_HELP_STRING([--with-etsf_io=DIR], [http://www.etsf.eu/resources/software/libraries_and_tools])])
case $with_etsf_io in
  yes | "") ;;
  no ) acx_etsf_io_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_ETSF_IO="$with_etsf_io" ;;
  *) LIBS_ETSF_IO="-l$with_etsf_io" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_etsf_io_save_LIBS="$LIBS"
acx_etsf_io_save_FCFLAGS="$FCFLAGS"

dnl The tests
if test $acx_etsf_io_ok = no; then
  AC_MSG_CHECKING([for etsf_io])  
  # If LIBS_ETSF_IO has been passed with --with-etsf_io just test this
  if test "$LIBS_ETSF_IO"; then
    etsf_io_fcflags="$LIBS_ETSF_IO"; etsf_io_libs="$LIBS_ETSF_IO"
FCFLAGS="$etsf_io_fcflags $netcdf_fcflags $acx_etsf_io_save_FCFLAGS"
LIBS="$etsf_io_libs $netcdf_libs $acx_etsf_io_save_LIBS"
AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    use etsf_io
    type(etsf_vars) :: vars
    call etsf_io_vars_free(vars)
]), [acx_etsf_io_ok=yes; FCFLAGS_ETSF_IO="$etsf_io_fcflags"; LIBS_ETSF_IO="$etsf_io_libs"], [])
  else
    etsf_io_libs="-letsf_io_utils -letsf_io"
FCFLAGS="$etsf_io_fcflags $netcdf_fcflags $acx_etsf_io_save_FCFLAGS"
LIBS=" $acx_etsf_io_save_LIBS $etsf_io_libs $netcdf_libs"
AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    use etsf_io
    type(etsf_vars) :: vars
    call etsf_io_vars_free(vars)
]), [acx_etsf_io_ok=yes; FCFLAGS_ETSF_IO="$etsf_io_fcflags"; LIBS_ETSF_IO="$etsf_io_libs"], [])
  fi
  if test $acx_netcdf_ok = no -o -z "$FCFLAGS_ETSF_IO $LIBS_ETSF_IO"; then
    AC_MSG_RESULT([$acx_etsf_io_ok])
  else
    AC_MSG_RESULT([$acx_etsf_io_ok ($LIBS_ETSF_IO)])
  fi
fi

AC_SUBST(FCFLAGS_ETSF_IO)
AC_SUBST(LIBS_ETSF_IO)
FCFLAGS="$acx_etsf_io_save_FCFLAGS"
LIBS="$acx_etsf_io_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_etsf_io_ok" = xyes; then
  AC_DEFINE(HAVE_ETSF_IO,1,[Defined if you have the ETSF_IO library.])
  $1
else
  if test $acx_etsf_io_ok != disable; then
    AC_MSG_WARN([Could not find etsf_io library. 
                *** Will compile without etsf_io support])
  fi
  $2
fi
])dnl ACX_ETSF_IO
