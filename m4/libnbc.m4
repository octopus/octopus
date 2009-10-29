AC_DEFUN([ACX_LIBNBC], [

dnl LibNBC is included in the distribution

AC_ARG_ENABLE(libnbc, AS_HELP_STRING([--disable-libnbc], [Do not compile with internal LibNBC library for non-blocking collectives.]),[acx_libnbc_ok=$enableval],[acx_libnbc_ok=yes])
dnl LibNBC is enabled by default

AC_MSG_CHECKING([whether LibNBC is enabled])

if test x"$acx_libnbc_ok" = x"yes"; then
   case $host_cpu in
   ia64) 
     acx_libnbc_ok=no
     AC_MSG_WARN([LibNBC does not work in ia64 systems: disabling it.])
     ;;
   sparc) 
     acx_libnbc_ok=no
     AC_MSG_WARN([LibNBC does not work in Sparc systems: disabling it.])
     ;;
   *)
     HAVE_LIBNBC=1
     AC_DEFINE(HAVE_LIBNBC, 1, [This is defined when we should compile with LibNBC support.])
     ;;
   esac
fi
AC_MSG_RESULT([$acx_libnbc_ok])

])dnl ACX_LIBNBC
