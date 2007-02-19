AC_DEFUN([ACX_LIBNBC], [

dnl LibNBC is included in the distribution but still experimental,
dnl so disabled by default

acx_libnbc_ok=no

dnl LibNBC is only enabled, if explicitly requested by the user
AC_ARG_ENABLE(libnbc, AS_HELP_STRING([--enable-libnbc], [Include support for non-blocking collectives via LibNBC (experimental)]), [acx_libnbc_ok=yes])

if test x"$acx_libnbc_ok" = x"yes"; then
  HAVE_LIBNBC=1
  AC_DEFINE(HAVE_LIBNBC, 1, [This is defined when we should compile with LibNBC support (not default).])
fi

])dnl ACX_LIBNBC
