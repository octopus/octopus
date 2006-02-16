AC_DEFUN([ACX_METIS], [

dnl Since v2.0, Metis ships with octopus, so we enable it by default
acx_metis_ok=yes

dnl We disable METIS support only if the user is requesting this explicitely
AC_ARG_ENABLE(metis, AC_HELP_STRING([--disable-metis], [Do not compile with internal Metis support]), [acx_metis_ok=no])

if test x"$acx_metis_ok" = xyes; then
  HAVE_METIS=1
  AC_DEFINE(HAVE_METIS, 1, [This is defined when we should compile with METIS support (default).])
else
  AC_MSG_WARN([METIS library has been disabled.
                *** octopus will not be able to use mesh partitioning.])
fi

])dnl ACX_METIS
