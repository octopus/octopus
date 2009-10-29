AC_DEFUN([ACX_NEWUOA], [

dnl NEWUOA is included in the distribution

dnl LibUOA is disabled by default
acx_newuoa_ok=no

AC_ARG_ENABLE(newuoa, AS_HELP_STRING([--enable-newuoa], [Enable support for the newuoa optimization software.]))

AC_MSG_CHECKING([whether newuoa is enabled])

acx_newuoa_ok="$enable_newuoa"

if test x"$acx_newuoa_ok" = x"yes"; then
  HAVE_NEWUOA=1
  AC_DEFINE(HAVE_NEWUOA, 1, [This is defined when the newuoa code is to be compiled.])
fi

AC_MSG_RESULT([$acx_newuoa_ok])

])dnl ACX_LIBUOA
