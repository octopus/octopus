## Copyright (C) 2019 S. Ohlmann
## based on the ELPA m4 file in this directory authored by X. Andrade
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

AC_DEFUN([ACX_LIKWID],
[

  acx_likwid_ok=no

  dnl BACKUP LIBS AND FCFLAGS
  acx_likwid_save_LIBS="$LIBS"
  acx_likwid_save_FCFLAGS="$FCFLAGS"

  dnl Check if the library was given in the command line
  AC_ARG_WITH(likwid-prefix, [AS_HELP_STRING([--with-likwid-prefix=DIR], [Directory where likwid was installed.])])

  # Set FCFLAGS_LIKWID only if not set from environment
  if test x"$FCFLAGS_LIKWID" = x; then
    case $with_likwid_prefix in
      "") FCFLAGS_LIKWID="-I/usr/include" ;;
      *)  FCFLAGS_LIKWID="-I$with_likwid_prefix/include" ;;
    esac
  fi

  AC_MSG_CHECKING([for likwid])

  likwid_program="AC_LANG_PROGRAM([],[
    use :: likwid
    implicit none

    call likwid_markerInit()
  ])"

  FCFLAGS="$FCFLAGS_LIKWID $acx_likwid_save_FCFLAGS"

  if test ! -z "$with_likwid_prefix"; then
    LIBS_LIKWID="-L$with_likwid_prefix/lib -llikwid"
  else
    LIBS_LIKWID="-llikwid"
  fi

  LIBS="$LIBS_LIKWID $acx_likwid_save_LIBS $LIBS_LAPACK $LIBS_BLAS"
  AC_LINK_IFELSE($likwid_program, [acx_likwid_ok=yes], [acx_likwid_ok=no])

  AC_MSG_RESULT([$acx_likwid_ok ($FCFLAGS_LIKWID $LIBS_LIKWID)])

  if test x$acx_likwid_ok != xyes; then

    AC_MSG_WARN([Could not find the likwid library])

    FCFLAGS_LIKWID=""
    LIBS_LIKWID=""

  else

    AC_DEFINE(HAVE_LIKWID, 1, [Define if LIKWID is available])

  fi

  AC_SUBST(FCFLAGS_LIKWID)
  AC_SUBST(LIBS_LIKWID)

  FCFLAGS="$acx_likwid_save_FCFLAGS"
  LIBS="$acx_likwid_save_LIBS"

])
