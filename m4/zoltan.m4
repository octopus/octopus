AC_DEFUN([ACX_ZOLTAN], [
compile_zoltan=no

AC_ARG_WITH(external_zoltan, [AS_HELP_STRING([--with-external_zoltan=LIB], [external installation of Zoltan mesh-partitioner])])
case $with_external_zoltan in
  yes | "") ;;
  no ) compile_zoltan=yes ;;
  -* | */* | *.a | *.so | *.so.* | *.o) LIBS_ZOLTAN="$with_external_zoltan";;
  *) LIBS_ZOLTAN="-l$with_external_zoltan" ;;
esac

dnl If in parallel, enable Zoltan compilation
if test x"$enable_mpi" != x"no" && test x"$compile_zoltan" == x"no"; then
  acx_zoltan_save_CPPFLAGS="$CPPFLAGS"
  acx_zoltan_save_LIBS="$LIBS"

  if test x"$CPPFLAGS_ZOLTAN" != x; then
    CPPFLAGS="$CPPFLAGS $CPPFLAGS_ZOLTAN"
  fi
  AC_CHECK_HEADERS(zoltan.h,[compile_zoltan=no],[compile_zoltan=yes])

  # if failed, try with a common location for trilinos installation
  if test x"$compile_zoltan" == x"yes"; then
    CPPFLAGS_ZOLTAN="-I/usr/include/trilinos"
    CPPFLAGS="$acx_zoltan_save_CPPFLAGS $CPPFLAGS_ZOLTAN"
    # remove cache so it actually checks again
    $as_unset AS_TR_SH([ac_cv_header_zoltan_h])
    AC_CHECK_HEADERS(zoltan.h,[compile_zoltan=no],[compile_zoltan=yes])
  fi

  # if couldn't find header, don't bother checking for function
  if test x"$compile_zoltan" == x"no"; then
    LIBS="$LIBS $LIBS_ZOLTAN"
    AC_CHECK_FUNCS(Zoltan_Initialize,[compile_zoltan=no],[compile_zoltan=yes])

    # if failed, try with a common location for trilinos installation
    if test x"$compile_zoltan" == x"yes"; then
      LIBS_ZOLTAN=-ltrilinos_zoltan
      LIBS="$acx_zoltan_save_LIBS $LIBS_ZOLTAN"
      # remove cache so it actually checks again
      $as_unset AS_TR_SH([ac_cv_func_Zoltan_Initialize])
      AC_CHECK_FUNCS(Zoltan_Initialize,[compile_zoltan=no],[compile_zoltan=yes])
    fi
  fi

  if test x"$compile_zoltan" == x"no"; then
    AC_DEFINE(HAVE_EXTERNAL_ZOLTAN, 1, [Zoltan is a external library])
    AC_MSG_NOTICE([Found a usable installation of Zoltan.])
    AC_SUBST(CPPFLAGS_ZOLTAN)
    AC_SUBST(LIBS_ZOLTAN)
  else
    AC_MSG_NOTICE([Could not find an installed version of Zoltan. It will be compiled.])
    CPPFLAGS_ZOLTAN=""	
    LIBS_ZOLTAN=""
  fi
  CPPFLAGS="$acx_zoltan_save_CPPFLAGS"
  LIBS="$acx_zoltan_save_LIBS"
fi

])dnl ACX_ZOLTAN
