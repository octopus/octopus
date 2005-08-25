AC_DEFUN([ACX_METIS], [
acx_metis_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(metis, [AC_HELP_STRING([--with-metis=DIR], [http://www-users.cs.umn.edu/~karypis/metis/])])
case $with_metis in
  yes | "") ;;
  no) acx_metis_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_METIS="$with_metis" ;;
  *) LIBS_METIS="-l$with_metis" ;;
esac

dnl Backup LIBS 
acx_metis_save_LIBS="$LIBS"

dnl First, check if it links
if test $acx_metis_ok = no; then
  LIBS="$LIBS_METIS $acx_metis_save_LIBS"
  AC_CHECK_FUNC([METIS_PartGraphRecursive], [acx_trlan_ok=yes], [])
fi

dnl Generic METIS library?
if test $acx_metis_ok = no; then
  LIBS="$LIBS_METIS -lmetis -lm $acx_metis_save_LIBS"
	unset ac_cv_func_METIS_PartGraphRecursive
  AC_CHECK_FUNC([METIS_PartGraphRecursive], [acx_metis_ok=yes; LIBS_METIS="$LIBS_METIS -lmetis -lm"], [])
fi

AC_SUBST(LIBS_METIS)
LIBS="$acx_metis_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_metis_ok" = xyes; then
  AC_DEFINE(HAVE_METIS,1,[Defined if you have METIS library.])
  $1
else
  if test $acx_metis_ok != disable; then
    AC_MSG_WARN([Could not find metis library. 
                *** Will compile without metis support])
  fi
  $2
fi
])dnl ACX_METIS
