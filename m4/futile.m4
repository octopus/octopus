## Copyright (C) 2020 M. Oliveira
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
##

AC_DEFUN([ACX_FUTILE], [
AC_REQUIRE([ACX_YAML])
AC_REQUIRE([ACX_BLAS])
acx_futile_ok=no

dnl Check if the library was given in the command line
dnl if not, use environment variables or defaults
AC_ARG_WITH(futile-prefix, [AS_HELP_STRING([--with-futile-prefix=DIR], [Directory where Futile was installed.])])

# Set FCFLAGS_FUTILE only if not set from environment
if test x"$FCFLAGS_FUTILE" = x; then
  case $with_futile_prefix in
    "") FCFLAGS_FUTILE="$ax_cv_f90_modflag/usr/include" ;;
    *)  FCFLAGS_FUTILE="$ax_cv_f90_modflag$with_futile_prefix/include" ;;
  esac
fi

AC_ARG_WITH(futile-include, [AS_HELP_STRING([--with-futile-include=DIR], [Directory where Futile Fortran headers were installed.])])
case $with_futile_include in
  "") ;;
  *)  FCFLAGS_FUTILE="$ax_cv_f90_modflag$with_futile_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_futile_save_LIBS="$LIBS"
acx_futile_save_FCFLAGS="$FCFLAGS"


dnl The test
AC_MSG_CHECKING([for Futile])

testprogram="AC_LANG_PROGRAM([],[[
  use yaml_parse
  use yaml_output
  use f_utils
  use dynamic_memory
  use dictionaries

  call yaml_mapping_open('Test')
  call yaml_map('foo',1)
  call yaml_mapping_close()

]])"


FCFLAGS="$FCFLAGS_FUTILE $acx_futile_save_FCFLAGS"

if test -z "$LIBS_FUTILE"; then
  if test ! -z "$with_futile_prefix"; then
    LIBS_FUTILE="-L$with_futile_prefix/lib"
  else
    LIBS_FUTILE=""
  fi
  LIBS_FUTILE="$LIBS_FUTILE -lfutile-1"
fi

LIBS="$LIBS_FUTILE $LIBS_LIBYAML $LIBS_BLAS $acx_futile_save_LIBS"

AC_LINK_IFELSE($testprogram, [acx_futile_ok=yes], [])

AC_MSG_RESULT([$acx_futile_ok ($FCFLAGS_FUTILE $LIBS_FUTILE)])


dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_futile_ok" = xyes; then
  AC_DEFINE(HAVE_FUTILE, 1, [Defined if you have the Futile library.])
else
  AC_MSG_WARN([Could not find Futile library.])
  LIBS_FUTILE=""
  FCFLAGS_FUTILE=""
fi

AC_SUBST(LIBS_FUTILE)
AC_SUBST(FCFLAGS_FUTILE)
LIBS="$acx_futile_save_LIBS"
FCFLAGS="$acx_futile_save_FCFLAGS"

])dnl ACX_FUTILE
