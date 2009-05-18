AC_DEFUN([ACX_ETSF_IO], [
AC_REQUIRE([ACX_NETCDF])
acx_etsf_io_ok=no

dnl We cannot use etsf_io if netcdf is not found
if test "x$acx_netcdf_ok" != xyes; then
  acx_etsf_io_ok=nonetcdf
  FCFLAGS_ETSF_IO=""
  LIBS_ETSF_IO=""
fi

dnl Check if the library was given in the command line
if test $acx_etsf_io_ok = no; then
  AC_ARG_WITH(etsf_io, [AS_HELP_STRING([--with-etsf-io=DIR], [http://www.etsf.eu/resources/software/libraries_and_tools])])
  case $with_etsf_io in
    yes | "") ;;
    no ) acx_etsf_io_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_ETSF_IO="$with_etsf_io" ;;
    *) LIBS_ETSF_IO="-l$with_etsf_io" ;;
  esac
fi


dnl First, check LIBS_ETSF_IO environment variable
if test "x$LIBS_ETSF_IO" != x; then
  save_LIBS="$LIBS"; LIBS="$LIBS_ETSF_IO $LIBS_NETCDF"
  save_FCFLAGS="$FCFLAGS"; FCFLAGS="$FCFLAGS_ETSF_IO $FCFLAGS_NETCDF $FCFLAGS"
  AC_MSG_CHECKING([for etsf_io in $LIBS_ETSF_IO])
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    use etsf_io
    type(etsf_vars) :: vars
    call etsf_io_vars_free(vars)
    ]), [acx_etsf_io_ok=yes], [])
    AC_MSG_RESULT([$acx_etsf_io_ok])
    LIBS="$save_LIBS"
    FCFLAGS="$save_FCFLAGS"
    if test $acx_etsf_io_ok = no; then
      LIBS_ETSF_IO=""
      FCFLAGS_ETSF_IO=""
    fi
fi


dnl Second, just check in -letsf_io
if test $acx_etsf_io_ok = no; then
  save_LIBS="$LIBS"; LIBS="-letsf_io $LIBS_NETCDF"
  save_FCFLAGS="$FCFLAGS"; FCFLAGS="$FCFLAGS"
  AC_MSG_CHECKING([for etsf_io in -letsf_io])
    AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    use etsf_io
    type(etsf_vars) :: vars
    call etsf_io_vars_free(vars)
    ]), [acx_etsf_io_ok=yes; FCFLAGS_ETSF_IO=""; LIBS_ETSF_IO="-letsf_io"], [])
    AC_MSG_RESULT([$acx_etsf_io_ok])
  LIBS="$save_LIBS"
  FCFLAGS="$save_FCFLAGS"
  if test $acx_etsf_io_ok = no; then
    LIBS_ETSF_IO=""
    FCFLAGS_ETSF_IO=""
  fi
fi

AC_SUBST(FCFLAGS_ETSF_IO)
AC_SUBST(LIBS_ETSF_IO)


dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_etsf_io_ok" = xyes; then
  ifelse([$1],,AC_DEFINE(HAVE_ETSF_IO,1,[Define if you have ETSF_IO library.]
),[$1])
        :
else
  if test $acx_etsf_io_ok != disable; then
    AC_MSG_WARN([Could not find etsf_io library. 
                *** Will compile without etsf_io support])
  fi
  acx_etsf_io_ok=no
  $2
fi
])dnl ACX_ETSF_IO
