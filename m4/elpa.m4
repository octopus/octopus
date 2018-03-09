## Copyright (C) 2016 X. Andrade
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

AC_DEFUN([ACX_ELPA],
[
  
  acx_elpa_ok=no
  acx_elpa_old_ok=no

  dnl BACKUP LIBS AND FCFLAGS
  acx_elpa_save_LIBS="$LIBS"
  acx_elpa_save_FCFLAGS="$FCFLAGS"

  dnl Check if the library was given in the command line
  AC_ARG_WITH(elpa-prefix, [AS_HELP_STRING([--with-elpa-prefix=DIR], [Directory where elpa was installed.])])

  # Set FCFLAGS_ELPA only if not set from environment
  if test x"$FCFLAGS_ELPA" = x; then
    case $with_elpa_prefix in
      "") FCFLAGS_ELPA="-I/usr/include" ;;
      *)  FCFLAGS_ELPA="-I$with_elpa_prefix/include" ;;
    esac
  fi

  AC_MSG_CHECKING([for elpa])

  elpa_program_old="AC_LANG_PROGRAM([],[
    use :: elpa1
    implicit none

    integer :: c1, c2, c3, err
    err = get_elpa_communicators(c1, 0, 0, c2, c3)

  ])"

  elpa_program="AC_LANG_PROGRAM([],[
    use :: elpa
    implicit none
    class(elpa_t), pointer :: e
    integer :: status

    status = elpa_init(20170403)
    e => elpa_allocate()

  ])"

  FCFLAGS="$FCFLAGS_ELPA $acx_elpa_save_FCFLAGS"

  if test ! -z "$with_elpa_prefix"; then
    LIBS_ELPA="-L$with_elpa_prefix/lib -lelpa"
  else
    LIBS_ELPA="-lelpa"
  fi

  LIBS="$LIBS_ELPA $acx_elpa_save_LIBS $LIBS_LAPACK $LIBS_BLAS"
  AC_LINK_IFELSE($elpa_program, [acx_elpa_ok=yes], [acx_elpa_ok=no])
  AC_LINK_IFELSE($elpa_program_old, [acx_elpa_old_ok=yes], [acx_elpa_old_ok=no])

  AC_MSG_RESULT([$acx_elpa_ok ($FCFLAGS_ELPA $LIBS_ELPA)])

  if test x$acx_elpa_ok != xyes; then

    if test x$acx_elpa_old_ok != xyes; then
      AC_MSG_WARN([Could not find the elpa library])
    else
      AC_MSG_WARN([Could only find an old version of the elpa library, compiling without elpa support])
    fi

    FCFLAGS_ELPA=""
    LIBS_ELPA=""

  else

    AC_DEFINE(HAVE_ELPA, 1, [Define if ELPA is available])

  fi

  AC_SUBST(FCFLAGS_ELPA)
  AC_SUBST(LIBS_ELPA)

  FCFLAGS="$acx_elpa_save_FCFLAGS"
  LIBS="$acx_elpa_save_LIBS"

])
