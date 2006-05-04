AC_DEFUN([ACX_NETCDF], [
acx_netcdf_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(netcdf, [AS_HELP_STRING([--with-netcdf=DIR], [http://www.unidata.ucar.edu/packages/netcdf/])])
case $with_netcdf in
  yes) ;;
  no | "") acx_netcdf_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_NETCDF="$with_netcdf" ;;
  *) LIBS_NETCDF="-l$with_netcdf" ;;
esac

dnl Backup LIBS 
acx_netcdf_save_LIBS="$LIBS"

dnl First, check if it links
if test $acx_netcdf_ok = no; then
  LIBS="$LIBS_NETCDF $acx_netcdf_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for nf90_close])
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
use netcdf
integer :: ncid
integer :: status
status=nf90_close(ncid)]), acx_netcdf_ok=yes, [])
  AC_MSG_RESULT($acx_netcdf_ok)
fi

dnl NETCDF library?
if test $acx_netcdf_ok = no; then
  LIBS="$LIBS_NETCDF -lnetcdf $acx_netcdf_save_LIBS $FLIBS"

  AC_MSG_CHECKING([for nf90_close in -lnetcdf])
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
use netcdf
integer :: ncid
integer :: status
status=nf90_close(ncid)]), [acx_netcdf_ok=yes; LIBS_NETCDF="$LIBS_NETCDF -lnetcdf"], [])
  AC_MSG_RESULT($acx_netcdf_ok)
fi

AC_SUBST(LIBS_NETCDF)
LIBS="$acx_netcdf_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_netcdf_ok" = xyes; then
  AC_DEFINE(HAVE_NETCDF,1,[Defined if you have NETCDF library.])
  $1
else
  if test $acx_netcdf_ok != disable; then
    AC_MSG_WARN([Could not find netcdf library. 
                *** Will compile without netcdf support])
  fi
  $2
fi
])dnl ACX_NETCDF
