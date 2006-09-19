AC_DEFUN([ACX_TRLAN], [
acx_trlan_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(trlan, [AS_HELP_STRING([--with-trlan=DIR], [http://www.nersc.gov/research/SIMON/trlan.html])])
case $with_trlan in
  yes) ;;
  no | "") acx_trlan_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_TRLAN="$with_trlan" ;;
  *) LIBS_TRLAN="-l$with_trlan" ;;
esac

dnl Backup LIBS 
acx_trlan_save_LIBS="$LIBS"

dnl First, check if it links
if test $acx_trlan_ok = no; then
  LIBS="$LIBS_TRLAN $acx_trlan_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for trl_info module])
  AC_LINK_IFELSE(AC_LANG_PROGRAM([], [use trl_info]), acx_trlan_ok=yes, [])
  if test $acx_trlan_ok = no; then
    AC_MSG_RESULT([$acx_trlan_ok])
  else
    AC_MSG_RESULT([$acx_trlan_ok ($LIBS_TRLAN)])
  fi
fi

dnl Generic TRLAN library?
if test $acx_trlan_ok = no; then
  LIBS="$LIBS_TRLAN -ltrlan $acx_trlan_save_LIBS $FLIBS"

  AC_MSG_CHECKING([for trl_info module in -ltrlan])
  AC_LINK_IFELSE(AC_LANG_PROGRAM([], [use trl_info]), [acx_trlan_ok=yes; LIBS_TRLAN="$LIBS_TRLAN -ltrlan"], [])
  if test $acx_trlan_ok = no; then
    AC_MSG_RESULT([$acx_trlan_ok])
  else
    AC_MSG_RESULT([$acx_trlan_ok ($LIBS_TRLAN)])
  fi
fi

AC_SUBST(LIBS_TRLAN)
LIBS="$acx_trlan_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_trlan_ok" = xyes; then
  AC_DEFINE(HAVE_TRLAN,1,[Defined if you have TRLAN library.])
  $1
else
  if test $acx_trlan_ok != disable; then
    AC_MSG_WARN([Could not find trlan library. 
                *** Will compile without trlan support])
  fi
  $2
fi
])dnl ACX_TRLAN
