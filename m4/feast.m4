dnl Copyright (C) 2014, Ask Hjorth Larsen
dnl Shameless cargo-cult bastardization of adjacent .m4 files

AC_DEFUN([ACX_FEAST], [

dnl FEAST works without MPI, but let's not go there.
AC_REQUIRE([ACX_MPI])
AC_REQUIRE([ACX_BLAS])

acx_feast_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(feast, [AS_HELP_STRING([--with-feast=<lib>], [use FEAST eigensolver library <lib>. <lib> should be the path to the parallel library, typically libpfeast.a])])
case $with_feast in
  yes | "") ;;
  no) acx_feast_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_FEAST="$with_feast" ;;
  *) LIBS_FEAST="-l$with_feast" ;;
esac

dnl Backup LIBS
acx_feast_save_LIBS="$LIBS"

dnl The tests
AC_MSG_CHECKING([for FEAST])
if test "$acx_feast_ok" != disabled; then
  LIBS="$LIBS_FEAST $LIBS_BLAS $acx_feast_save_LIBS"
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    call zfeast_grcix()
    ]), [acx_feast_ok=yes], [])
fi
AC_MSG_RESULT([$acx_feast_ok ($LIBS_FEAST)])

LIBS="$acx_feast_save_LIBS"

if test x"$acx_feast_ok" = xyes; then
  AC_DEFINE(HAVE_FEAST,1,[Defined if you have the FEAST library.])
  AC_SUBST(LIBS_FEAST)
else
  AC_MSG_WARN([Could not find FEAST library.
      	      *** Will compile without FEAST support])
fi

])
