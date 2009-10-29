AC_DEFUN([ACX_NEWUOA], [

dnl NEWUOA is included in the distribution

AC_ARG_ENABLE(newuoa, AS_HELP_STRING([--enable-newuoa], [Compile with internal NEWUOA optimization library.]),[acx_newuoa_ok=$enableval],[acx_newuoa_ok=no])
dnl NEWUOA is disabled by default

AC_MSG_CHECKING([whether NEWUOA is enabled])

if test x"$acx_newuoa_ok" = x"yes"; then
  HAVE_NEWUOA=1
  AC_DEFINE(HAVE_NEWUOA, 1, [This is defined when the NEWUOA code is to be compiled.])
fi

AC_MSG_RESULT([$acx_newuoa_ok])

])dnl ACX_NEWUOA
