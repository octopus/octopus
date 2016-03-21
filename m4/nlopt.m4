dnl Copyright (c) 2011 David Strubbe

AC_DEFUN([ACX_NLOPT], [
acx_nlopt_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(nlopt-prefix, [AS_HELP_STRING([--with-nlopt-prefix=DIR], [Directory where NLOPT was installed. The NLOPT library can be downloaded at http://ab-initio.mit.edu/wiki/index.php/Main_Page])])
if test "x$with_nlopt_prefix" = xno; then
  acx_nlopt_ok=disabled
fi

dnl Backup LIBS and FCFLAGS
acx_nlopt_save_LIBS="$LIBS"
acx_nlopt_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for NLOPT])
if test "$acx_nlopt_ok" != disabled; then
  for libdir in "lib" "lib64"; do
    if test "$libdir" = "lib"; then
      LIBS_NLOPT="-L$with_nlopt_prefix/lib -lnlopt"; 
      FCFLAGS_NLOPT="$ax_cv_f90_modflag$with_nlopt_prefix/include";
    else
      LIBS_NLOPT="-L$with_nlopt_prefix/lib64 -lnlopt"; 
      FCFLAGS_NLOPT="$ax_cv_f90_modflag$with_nlopt_prefix/include";
    fi
	
    FCFLAGS="$FCFLAGS_NLOPT $acx_nlopt_save_FCFLAGS"
    LIBS="$LIBS_NLOPT $acx_nlopt_save_LIBS"
    AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
      integer(8) :: nlopt
      integer :: nparams
      include 'nlopt.f'
      call nlo_create(opt, NLOPT_LD_LBFGS, nparams)
      ]), [acx_nlopt_ok=yes; break], [])
  done
fi
AC_MSG_RESULT([$acx_nlopt_ok ($FCFLAGS_NLOPT $LIBS_NLOPT)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_nlopt_ok" = xyes; then
  AC_DEFINE(HAVE_NLOPT,1,[Defined if you have the NLOPT library.])
else
  AC_MSG_WARN([Could not find NLOPT library. 
           *** Will compile without NLOPT support])
  FCFLAGS_NLOPT=""
  LIBS_NLOPT=""
fi

AC_SUBST(FCFLAGS_NLOPT)
AC_SUBST(LIBS_NLOPT)

FCFLAGS="$acx_nlopt_save_FCFLAGS"
LIBS="$acx_nlopt_save_LIBS"
])dnl ACX_NLOPT
