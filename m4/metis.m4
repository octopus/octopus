AC_DEFUN([ACX_METIS], [

dnl METIS is included in the distribution

dnl We disable METIS support only if the user is requesting this explicitly
AC_ARG_ENABLE(metis, AS_HELP_STRING([--disable-metis], [Do not compile with internal METIS domain-partitioning library.]),[acx_metis_ok=$enableval],[acx_metis_ok=yes])
dnl Since v2.0, METIS ships with octopus, so we enable it by default

AC_MSG_CHECKING([whether METIS is enabled])

AC_MSG_RESULT([$acx_metis_ok])

if test x"$acx_metis_ok" = xyes; then
  HAVE_METIS=1
  AC_DEFINE(HAVE_METIS, 1, [This is defined when we should compile with METIS support (default).])
fi

])dnl ACX_METIS
