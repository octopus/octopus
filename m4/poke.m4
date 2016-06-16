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
## $Id: fft.m4 14470 2015-07-29 09:02:39Z dannert $
##

AC_DEFUN([ACX_POKE],
[
  
  acx_poke_ok=no

  dnl BACKUP LIBS AND FCFLAGS
  acx_poke_save_LIBS="$LIBS"
  acx_poke_save_FCFLAGS="$FCFLAGS"

  dnl Check if the library was given in the command line
  AC_ARG_WITH(poke-prefix, [AS_HELP_STRING([--with-poke-prefix=DIR], [Directory where poke was installed.])])

  # Set FCFLAGS_POKE only if not set from environment
  if test x"$FCFLAGS_POKE" = x; then
    case $with_poke_prefix in
      "") FCFLAGS_POKE="-I/usr/include" ;;
      *)  FCFLAGS_POKE="-I$with_poke_prefix/include" ;;
    esac
  fi

  AC_MSG_CHECKING([for poke])

  poke_program="AC_LANG_PROGRAM([],[
    use :: poke
    implicit none

    type(PokeGrid) :: grid

    grid = PokeGrid(1.0_8, (/10, 10, 10/))

  ])"

  FCFLAGS="$FCFLAGS_POKE $acx_poke_save_FCFLAGS"

  if test ! -z "$with_poke_prefix"; then
    LIBS_POKE="-L$with_poke_prefix/lib -lpoke"
  else
    LIBS_POKE="-lpoke"
  fi

  LIBS="$LIBS_POKE $acx_poke_save_LIBS $LIBS_FFTW"
  AC_LINK_IFELSE($poke_program, [acx_poke_ok=yes], [acx_poke_ok=no])

  AC_MSG_RESULT([$acx_poke_ok ($FCFLAGS_POKE $LIBS_POKE)])

  if test x$acx_poke_ok != xyes; then

    AC_MSG_WARN([Could not find the poke library])

    FCFLAGS_POKE=""
    LIBS_POKE=""

  else

    AC_DEFINE(HAVE_POKE, 1, [Define if POKE is available])

  fi

  AC_SUBST(FCFLAGS_POKE)
  AC_SUBST(LIBS_POKE)

  FCFLAGS="$acx_poke_save_FCFLAGS"
  LIBS="$acx_poke_save_LIBS"

])
