## Copyright (C) 2012 U. De Giovannini
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

dnl looks for libarpack.a
AC_DEFUN([ACX_ARPACK], [
AC_REQUIRE([ACX_BLAS])
acx_arpack_ok=no

dnl We cannot use ARPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
  acx_arpack_ok=noblas
fi

dnl Backup LIBS 
acx_arpack_save_LIBS="$LIBS"


dnl Check if the library was given in the command line
if test $acx_arpack_ok = no; then
  AC_ARG_WITH(arpack, [AS_HELP_STRING([--with-arpack=<lib>], [use ARPACK library <lib> http://forge.scilab.org/index.php/p/arpack-ng/])])
  case $with_arpack in
    yes | "") ;;
    no) acx_arpack_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_ARPACK="$with_arpack" ;;
    *) LIBS_ARPACK="-l$with_arpack" ;;
  esac
fi

testprog="AC_LANG_PROGRAM([],[call dsaupd])"

dnl First, check LIBS_ARPACK environment variable
if test $acx_arpack_ok = no; then
  LIBS="$LIBS_ARPACK $LIBS_LAPACK $LIBS_BLAS $acx_arpack_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for arpack library])
  AC_LINK_IFELSE($testprog, [acx_arpack_ok=yes], [])
  if test $acx_arpack_ok = no; then
    AC_MSG_RESULT([$acx_arpack_ok])
  else
    AC_MSG_RESULT([$acx_arpack_ok ($LIBS_ARPACK)])
  fi
fi

if test $acx_arpack_ok = no; then
  AC_MSG_CHECKING([for arpack library with -larpack])
  if test "$LIBS_ARPACK" = ""; then
    LIBS=" -larpack $LIBS_LAPACK $LIBS_BLAS $acx_arpack_save_LIBS $FLIBS"
    AC_LINK_IFELSE($testprog, [acx_arpack_ok=yes; LIBS_ARPACK=" -larpack"], [])
  else
    LIBS="-L$LIBS_ARPACK -larpack $LIBS_LAPACK $LIBS_BLAS $acx_arpack_save_LIBS $FLIBS"
    AC_LINK_IFELSE($testprog, [acx_arpack_ok=yes; LIBS_ARPACK="-L$LIBS_ARPACK -larpack"], [])  
  fi
  if test $acx_arpack_ok = no; then
    AC_MSG_RESULT([$acx_arpack_ok])
  else
    AC_MSG_RESULT([$acx_arpack_ok ($LIBS_ARPACK)])
  fi
fi

LIBS="$acx_arpack_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_arpack_ok" = xyes; then
  AC_SUBST(LIBS_ARPACK)
  AC_DEFINE(HAVE_ARPACK,1,[Defined if you have ARPACK library.])
  $1
else
    AC_MSG_WARN([Could not find ARPACK library. 
               *** Will compile without ARPACK support])
  $2
fi
])dnl ACX_ARPACK
