AC_DEFUN([ACX_NETCDF], [
acx_netcdf_ok=no
FCFLAGS_NETCDF=""
LIBS_NETCDF=""

dnl Check if the library was given in the command line
AC_ARG_WITH(netcdf-prefix, [AS_HELP_STRING([--with-netcdf-prefix=DIR], [Directory where netcdf was installed.])])
case $with_netcdf_prefix in
  no ) acx_netcdf_ok=disabled ;;
  "") LIBS_NETCDF=""; FCFLAGS_NETCDF="-I/usr/include" ;;
  *) LIBS_NETCDF="-L$with_netcdf_prefix/lib"; FCFLAGS_NETCDF="$ax_cv_f90_modflag$with_netcdf_prefix/include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_netcdf_save_LIBS="$LIBS"
acx_netcdf_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for netcdf])
if test "$acx_netcdf_ok" != disabled; then
  for netcdf_libsl in "-lnetcdff -lnetcdf" -lnetcdf; do
    netcdf_fcflags="$FCFLAGS_NETCDF"
    netcdf_libs="$LIBS_NETCDF $netcdf_libsl"
    FCFLAGS="$netcdf_fcflags $acx_netcdf_save_FCFLAGS"
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

AC_SUBST(FCFLAGS_NETCDF)
AC_SUBST(LIBS_NETCDF)
FCFLAGS="$acx_netcdf_save_FCFLAGS"
LIBS="$acx_netcdf_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_netcdf_ok" = xyes; then
  AC_DEFINE(HAVE_NETCDF,1,[Defined if you have NETCDF library.])
  $1
else
  AC_MSG_WARN([Could not find NetCDF library. 
              *** Will compile without NetCDF support])
  $2
fi
])dnl ACX_NETCDF
