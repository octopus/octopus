AC_DEFUN([ACX_NETCDF], [
acx_netcdf_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(netcdf-prefix, [AS_HELP_STRING([--with-netcdf-prefix=DIR], [Directory where netcdf was installed.])])
case $with_netcdf_prefix in
  no ) acx_netcdf_ok=disabled ;;
  "") if test x$FCFLAGS_NETCDF == x; then
    FCFLAGS_NETCDF="-I/usr/include"
  fi;;
  *) LIBS_NETCDF="-L$with_netcdf_prefix/lib"; FCFLAGS_NETCDF="$ax_cv_f90_modflag$with_netcdf_prefix/include" ;;
esac

AC_ARG_WITH(netcdf-include, [AS_HELP_STRING([--with-netcdf-include=DIR], [Directory where netcdf Fortran headers were installed.])])
case $with_netcdf_include in
  "") ;;
  *)  FCFLAGS_NETCDF="$ax_cv_f90_modflag$with_netcdf_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_netcdf_save_LIBS="$LIBS"
acx_netcdf_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for netcdf])
if test "$acx_netcdf_ok" != disabled; then
  netcdf_fcflags="$FCFLAGS_NETCDF"
  FCFLAGS="$netcdf_fcflags $acx_netcdf_save_FCFLAGS"
  for netcdf_libsl in "" -lnetcdf "-lnetcdff -lnetcdf"; do
    netcdf_libs="$LIBS_NETCDF $netcdf_libsl"
    LIBS="$netcdf_libs $acx_netcdf_save_LIBS"
    AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
      use netcdf
      integer :: ncid
      integer :: status
      status = nf90_close(ncid)
    ]), [acx_netcdf_ok=yes; FCFLAGS_NETCDF="$netcdf_fcflags"; LIBS_NETCDF="$netcdf_libs"], [])
    if test $acx_netcdf_ok == yes; then 
      LIBS_NETCDF=$netcdf_libs
      break
    fi
  done
fi
AC_MSG_RESULT([$acx_netcdf_ok ($LIBS_NETCDF)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_netcdf_ok" = xyes; then
  AC_DEFINE(HAVE_NETCDF,1,[Defined if you have NETCDF library.])
  $1
else
  AC_MSG_WARN([Could not find NetCDF library. 
              *** Will compile without NetCDF and ETSF I/O support])
  FCFLAGS_NETCDF=""
  LIBS_NETCDF=""
  $2
fi

AC_SUBST(FCFLAGS_NETCDF)
AC_SUBST(LIBS_NETCDF)
FCFLAGS="$acx_netcdf_save_FCFLAGS"
LIBS="$acx_netcdf_save_LIBS"
])dnl ACX_NETCDF
