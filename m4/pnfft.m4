## Copyright (C) 2013 U. De Giovannini
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

AC_DEFUN([ACX_PNFFT], [
AC_REQUIRE([ACX_FFT])
AC_REQUIRE([ACX_PFFT])
acx_pnfft_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(pnfft-prefix, [AS_HELP_STRING([--with-pnfft-prefix=DIR], [http://www-user.tu-chemnitz.de/~mpip/software.php. It requires PFFT.])])
case $with_pnfft_prefix in
  yes | "") ;;
  no) acx_pnfft_ok=disable ;;
  *.a | *.so | *.so.* | *.o) LIBS_PFFT="-L$with_pfft_prefix" ;;
  *) LIBS_PNFFT="-L$with_pnfft_prefix/lib"; 
     FCFLAGS_PNFFT="$ax_cv_f90_modflag$with_pnfft_prefix/include" ;;  
esac

dnl The include dir must be specified when the library is given with a 
dnl specified file to be compiled static (i.e. *.a etc.)
AC_ARG_WITH(pnfft-include, [AS_HELP_STRING([--with-pnfft-include=DIR], [Directory where PNFFT Fortran include files were installed.])])
case $with_pnfft_include in
  "") if test "x$FCFLAGS_PNFFT" == x; then
  FCFLAGS_PNFFT="$ax_cv_f90_modflag /usr/include"
  fi;;
  *)  FCFLAGS_PNFFT="$ax_cv_f90_modflag$with_pnfft_include" ;;
esac


dnl We cannot use PNFFT if PFFT is not found
if test "x$acx_pfft_ok" != xyes; then
  acx_pnfft_ok=nopfft
fi

dnl Backup LIBS and FCFLAGS
acx_pnfft_save_LIBS="$LIBS"
acx_pnfft_save_FCFLAGS="$FCFLAGS"

FCFLAGS_PNFFT="$FCFLAGS_PNFFT $FCFLAGS_PFFT"
FCFLAGS="$FCFLAGS_PNFFT  $acx_pnfft_save_FCFLAGS"    


testprogram="AC_LANG_PROGRAM([],[ 
    use iso_c_binding
    include \"fftw3-mpi.f03\"
    include \"pfft.f03\"
    include \"pnfft.f03\"
    
    call pnfft_init()
  ])"


dnl First, check LIBS_PNFFT environment variable
if test x"$acx_pnfft_ok" = xno; then
  LIBS="$LIBS_PNFFT $LIBS_PFFT $LIBS"

  AC_MSG_CHECKING([for pnfft library])
  AC_LINK_IFELSE($testprogram, [acx_pnfft_ok=yes; LIBS_PNFFT="$LIBS_PNFFT $LIBS_PFFT"], [])

  if test $acx_pnfft_ok = no; then
    AC_MSG_RESULT([$acx_pnfft_ok])
  else
    AC_MSG_RESULT([$acx_pnfft_ok ($LIBS_PNFFT)])
  fi
fi


dnl Generic PNFFT library 
if test $acx_pnfft_ok = no; then
  AC_MSG_CHECKING([for pnfft library with -lpnfft])
  if test "$LIBS_PNFFT" = ""; then
    LIBS="-lpnfft $LIBS_PFFT $LIBS"
    AC_LINK_IFELSE($testprogram, [acx_pnfft_ok=yes; LIBS_PNFFT=" -lpnfft $LIBS_PFFT"], [])
  else
    LIBS="$LIBS_PNFFT -lpnfft $LIBS_PFFT $LIBS"
    AC_LINK_IFELSE($testprogram, [acx_pnfft_ok=yes; LIBS_PNFFT="$LIBS_PNFFT -lpnfft $LIBS_PFFT"], [])  
  fi
  if test $acx_pnfft_ok = no; then
    AC_MSG_RESULT([$acx_pnfft_ok])
  else
    AC_MSG_RESULT([$acx_pnfft_ok ($LIBS_PNFFT)])
  fi
fi


dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pnfft_ok" = xyes; then
  AC_DEFINE(HAVE_PNFFT,1,[Defined if you have PNFFT library.])
  $1
else
  AC_MSG_WARN([Could not find PNFFT library. 
               *** Will compile without PNFFT support])
  LIBS_PNFFT=""
  FCFLAGS_PNFFT=""
  $2
fi

AC_SUBST(LIBS_PNFFT)
AC_SUBST(FCFLAGS_PNFFT)
LIBS="$acx_pnfft_save_LIBS"
FCFLAGS="$acx_pnfft_save_FCFLAGS"

])dnl ACX_PNFFT
