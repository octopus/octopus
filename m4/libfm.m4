dnl NOT available from the GNU Autoconf Macro Archive at:
dnl https://trac.version.fz-juelich.de/libfm
dnl
AC_DEFUN([ACX_LIBFM], [
acx_libfm_ok=no

dnl We cannot use LIBFM if MPI is not found
if test "x$acx_mpi_ok" != xyes; then
  acx_libfm_ok=nompi
fi

dnl Get fortran linker name of LIBFM function to check for.
dnl if not compiling with fortran, convert the names
m4_if(_AC_LANG, Fortran, [fmm_init=fcs_run], [AC_FC_FUNC(fmm_init)])

dnl Check if the library was given in the command line
if test $acx_libfm_ok = no; then
  AC_ARG_WITH(libfm, [AS_HELP_STRING([--with-libfm=<lib>], [use LIBFM library <lib>])])
  case $with_libfm in
    yes | "") ;;
    no) acx_libfm_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_LIBFM="$with_libfm" ;;
    *) LIBS_LIBFM="-l$with_libfm" ;;
  esac
fi

dnl Backup LIBS 
acx_libfm_save_LIBS="$LIBS"
LIBS="$LIBS_LIBFM $LIBS $FLIBS"

dnl First, check LIBS_LIBFM environment variable
if test $acx_libfm_ok = no; then
  AC_MSG_CHECKING([for $fmm_init in $LIBS_LIBFM])
  AC_TRY_LINK_FUNC($fmm_init, [acx_libfm_ok=yes], [])
  if test $acx_libfm_ok = no; then
    AC_MSG_RESULT([$acx_libfm_ok ($LIBS_LIBFM)])
  else
    AC_MSG_RESULT([$acx_libfm_ok ($LIBS_LIBFM)])
  fi
fi

dnl Generic LIBFM library?
for libfm in fm_r64; do
  if test $acx_libfm_ok = no; then
    AC_CHECK_LIB($libfm , $fmm_init,
      [acx_libfm_ok=yes; LIBS_LIBFM="$LIBS_LIBFM -l$libfm"], [], [$FLIBS])
  fi
done

dnl Generic LIBFM library?
for libfm in fm_r64; do
  dnl if test x"$libfm" = xlibfm-openmpi; then       
  dnl   libfmCinit="libfmCinit-openmpi"
  dnl else
  dnl libfm="fm_r64"
  dnl fi
  if test $acx_libfm_ok = no; then
    AC_CHECK_LIB($libfm, $fmm_init,
      [acx_libfm_ok=yes; LIBS_LIBFM="$LIBS_LIBFM -l$libfm"], [], [$FLIBS])
  fi
done

AC_SUBST(LIBS_LIBFM)
LIBS="$acx_libfm_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libfm_ok" = xyes; then
  AC_DEFINE(HAVE_LIBFM,1,[Defined if you have LIBFM library.])
  $1
else
  AC_MSG_WARN([Could not find Libfm library. 
               *** Will compile without Libfm support])
  $2
fi
])dnl ACX_LIBFM
