AC_DEFUN([ACX_LIBNBC], [

dnl LibNBC is included in the distribution
acx_libnbc_ok=yes

dnl LibNBC is enabled by default
AC_ARG_ENABLE(libnbc, AS_HELP_STRING([--disable-libnbc], [Disable support for non-blocking collectives via LibNBC]), [acx_libnbc_ok=no])

AC_MSG_CHECKING([for libnbc])

if test x"$acx_libnbc_ok" = x"yes"; then
  HAVE_LIBNBC=1
  AC_DEFINE(HAVE_LIBNBC, 1, [This is defined when we should compile with LibNBC support.])
  AC_MSG_RESULT([enabled])
else
  AC_MSG_RESULT([disabled])
fi

])dnl ACX_LIBNBC
