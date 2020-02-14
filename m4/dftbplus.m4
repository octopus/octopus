## Copyright (C) 2020 H. Appel, S. Ohlmann
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

AC_DEFUN([ACX_DFTBPLUS], [
acx_dftbplus_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(dftbplus-prefix, [AS_HELP_STRING([--with-dftbplus-prefix=<lib>], [https://www.dftbplus.org/])])

case $with_dftbplus_prefix in
  yes | "") ;;
  no) acx_dftbplus_ok=disable ;;
  *.a | *.so | *.so.* | *.o) LIBS_DFTBPLUS=$with_dftbplus_prefix ;
     xpath=${with_dftbplus_prefix%/lib/*}
     FCFLAGS_DFTBPLUS="$ax_cv_f90_modflag$xpath/include/dftb+/modfiles";;
  *) LIBS_DFTBPLUS="-L$with_dftbplus_prefix/lib -ldftbplus";
     FCFLAGS_DFTBPLUS="$ax_cv_f90_modflag$with_dftbplus_prefix/include/dftb+/modfiles" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_dftbplus_save_LIBS="$LIBS"
acx_dftbplus_save_FCFLAGS="$FCFLAGS"

FCFLAGS="$FCFLAGS_DFTBPLUS $acx_dftbplus_save_FCFLAGS"

testprogram="AC_LANG_PROGRAM([],[
    use dftbplus
    implicit none

    type(TDftbPlus) :: dftbp

    call TDftbPlus_init(dftbp)
  ])"

if test $acx_dftbplus_ok = no; then
  AC_MSG_CHECKING([for dftbplus library])
  LIBS="$LIBS_DFTBPLUS $LIBS_BLAS $LIBS_LAPACK $acx_dftbplus_save_LIB"
  AC_LINK_IFELSE($testprogram, [acx_dftbplus_ok=yes], [])

  if test $acx_dftbplus_ok = no; then
    AC_MSG_RESULT([$acx_dftbplus_ok])
  else
    AC_MSG_RESULT([$acx_dftbplus_ok ($LIBS_DFTBPLUS)])
  fi
fi

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_dftbplus_ok" = xyes; then
  AC_DEFINE(HAVE_DFTBPLUS,1,[Defined if you have DFTBPLUS library.])
else
  AC_MSG_WARN([Could not find DFTBPLUS library.
               *** Will compile without DFTBPLUS support])
  LIBS_DFTBPLUS=""
  FCFLAGS_DFTBPLUS=""
fi

AC_SUBST(LIBS_DFTBPLUS)
AC_SUBST(FCFLAGS_DFTBPLUS)
LIBS="$acx_dftbplus_save_LIBS"
FCFLAGS="$acx_dftbplus_save_FCFLAGS"

])dnl ACX_DFTBPLUS
