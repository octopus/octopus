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

AC_DEFUN([ACX_YAML], [
acx_libyaml_ok=no

dnl Check if the library was given in the command line
dnl if not, use environment variables or defaults
AC_ARG_WITH(yaml-prefix, [AS_HELP_STRING([--with-yaml-prefix=DIR], [Directory where LibYAML was installed.])])

# Set CFLAGS_LIBYAML only if not set from environment
if test x"$CFLAGS_LIBYAML" = x; then
  case $with_yaml_prefix in
    "") CFLAGS_LIBYAML="-I/usr/include" ;;
    *)  CFLAGS_LIBYAML="-I$with_yaml_prefix/include" ;;
  esac
fi

dnl Backup LIBS and CFLAGS
acx_libyaml_save_LIBS="$LIBS"
acx_libyaml_save_CFLAGS="$CFLAGS"

if test -z "$LIBS_LIBYAML"; then
  case $with_yaml_prefix in
    "") LIBS_LIBYAML="" ;;
    *)  LIBS_LIBYAML="-L$with_yaml_prefix/lib" ;;
  esac
fi

LIBS="$LIBS_LIBYAML $acx_libyaml_save_LIBS"
AC_CHECK_LIB([yaml], [yaml_parser_parse], [acx_libyaml_ok=yes], [acx_libyaml_ok=no])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libyaml_ok" = xyes; then
  LIBS_LIBYAML="$LIBS_YAML -lyaml"
  AC_DEFINE([HAVE_YAML], 1, [Defined if you have the LibYAML library.])

else
  AC_MSG_WARN([Could not find LibYAML library.])
  LIBS_LIBYAML=""
  CFLAGS_LIBYAML=""
fi

AC_SUBST(LIBS_LIBYAML)
AC_SUBST(CFLAGS_LIBYAML)
LIBS="$acx_libyaml_save_LIBS"
CFLAGS="$acx_libyaml_save_CFLAGS"

])dnl ACX_LIBYAML
