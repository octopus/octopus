dnl Copyright (c) 2011 David Strubbe

AC_DEFUN([ACX_BERKELEYGW], [
acx_berkeleygw_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(berkeleygw-prefix, [AS_HELP_STRING([--with-berkeleygw-prefix=DIR], [Directory where BerkeleyGW was installed.])])
if test "x$with_berkeleygw_prefix" = xno; then
  acx_berkeleygw_ok=disabled
fi

dnl Backup LIBS and FCFLAGS
acx_berkeleygw_save_LIBS="$LIBS"
acx_berkeleygw_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for BerkeleyGW])
if test "$acx_berkeleygw_ok" != disabled; then
  for location in "source" "install"; do
    if test "$location" = "source"; then
      LIBS_BERKELEYGW="-L$with_berkeleygw_prefix/library -lBGW_wfn"; 
      FCFLAGS_BERKELEYGW="$ax_cv_f90_modflag$with_berkeleygw_prefix/library";
    else
      LIBS_BERKELEYGW="-L$with_berkeleygw_prefix/lib -lBGW_wfn"; 
      FCFLAGS_BERKELEYGW="$ax_cv_f90_modflag$with_berkeleygw_prefix/include";
    fi
	
    FCFLAGS="$FCFLAGS_BERKELEYGW $acx_berkeleygw_save_FCFLAGS"
    LIBS="$LIBS_BERKELEYGW $acx_berkeleygw_save_LIBS"
    AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
      use wfn_rho_vxc_io_m
      call bgw_conf_test()
      ]), [acx_berkeleygw_ok=yes; break], [])
  done
fi
AC_MSG_RESULT([$acx_berkeleygw_ok ($FCFLAGS_BERKELEYGW $LIBS_BERKELEYGW)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_berkeleygw_ok" = xyes; then
  AC_DEFINE(HAVE_BERKELEYGW,1,[Defined if you have the BerkeleyGW library.])
else
  AC_MSG_WARN([Could not find BerkeleyGW library. 
           *** Will compile without BerkeleyGW support])
  FCFLAGS_BERKELEYGW=""
  LIBS_BERKELEYGW=""
fi

AC_SUBST(FCFLAGS_BERKELEYGW)
AC_SUBST(LIBS_BERKELEYGW)

FCFLAGS="$acx_berkeleygw_save_FCFLAGS"
LIBS="$acx_berkeleygw_save_LIBS"
])dnl ACX_BERKELEYGW
