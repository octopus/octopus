## Copyright (C) 2011 J. Alberdi
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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id: pfft.m4 6722 2010-06-13 12:44:43Z joseba $
##

AC_DEFUN([ACX_PFFT], [
AC_REQUIRE([ACX_FFT])
acx_pfft_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(pfft-prefix, [AS_HELP_STRING([--with-pfft-prefix=DIR], [http://www-user.tu-chemnitz.de/~mpip/software.php])])
case $with_pfft_prefix in
  no) acx_pfft_ok=disable ;;
  *) LIBS_PFFT="-L$with_pfft_prefix/lib -lpfft"; 
     FCFLAGS_PFFT="$ax_cv_f90_modflag$with_pfft_prefix/include" ;;
esac

dnl We only use the include file 'fftw3.f' with PFFT but not otherwise with FFTW, so fft.m4 did not check for it.
AC_ARG_WITH(fft-include, [AS_HELP_STRING([--with-fft-include=DIR], [Directory where FFTW Fortran include files were installed.])])
case $with_fft_include in
 "") FCFLAGS_FFT="$ax_cv_f90_modflag /usr/include";;
 *)  FCFLAGS_FFT="$ax_cv_f90_modflag$with_fft_include" ;;
esac

AC_ARG_WITH(pfft-include, [AS_HELP_STRING([--with-pfft-include=DIR], [Directory where PFFT Fortran include files were installed.])])
case $with_pfft_include in
  "") if test "x$FCFLAGS_PFFT" == x; then
  FCFLAGS_PFFT="$ax_cv_f90_modflag /usr/include"
  fi;;
  *)  FCFLAGS_PFFT="$ax_cv_f90_modflag$with_pfft_include" ;;
esac
FCFLAGS_PFFT="$FCFLAGS_PFFT $FCFLAGS_FFT"

dnl We cannot use PFFT if MPI is not found
if test "x$acx_mpi_ok" != xyes; then
  acx_pfft_ok=nompi
fi

dnl We cannot use PFFT if FFTW3 is not found
if test "x$acx_fft_ok" != xyes; then
  acx_pfft_ok=nofftw3
fi

dnl Backup LIBS and FCFLAGS
acx_pfft_save_LIBS="$LIBS"
acx_pfft_save_FCFLAGS="$FCFLAGS"

LIBS="$LIBS_PFFT $LIBS_FFT $LIBS $FLIBS"
FCFLAGS="$FCFLAGS_PFFT $FCFLAGS"

pfft_func="dpfft_plan_dft_3d"

dnl First, check LIBS_PFFT environment variable
if test x"$acx_pfft_ok" = xno; then
  AC_MSG_CHECKING([for $pfft_func in $FCFLAGS_PFFT $LIBS_PFFT])
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    include 'fftw3.f'
    include 'pfft.f'
    integer :: x = PFFT_REDFT00
    call dpfft_plan_dft_3d()
  ]), [acx_pfft_ok=yes], [])
fi

AC_MSG_RESULT([$acx_pfft_ok ($FCFLAGS_PFFT $LIBS_PFFT)])

dnl Generic PFFT library? is there such a thing?
for pfft in pfft; do
  if test x"$acx_pfft_ok" = xno; then
    AC_CHECK_LIB($pfft -lfftw3_mpi -lfftw3 -lfftw3_mpi -lm, $pfft_func,
      [acx_pfft_ok=yes; LIBS_PFFT="$LIBS_PFFT -lfftw3 -lfftw3_mpi -lm"], [], [$FLIBS])
  fi
done

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pfft_ok" = xyes; then
  AC_DEFINE(HAVE_PFFT,1,[Defined if you have PFFT library.])
  $1
else
  AC_MSG_WARN([Could not find PFFT library. 
               *** Will compile without PFFT support])
  LIBS_PFFT=""
  FCFLAGS_PFFT=""
  $2
fi

AC_SUBST(LIBS_PFFT)
AC_SUBST(FCFLAGS_PFFT)
LIBS="$acx_pfft_save_LIBS"
FCFLAGS="$acx_pfft_save_FCFLAGS"

])dnl ACX_PFFT
