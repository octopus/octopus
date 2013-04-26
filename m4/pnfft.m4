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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
##

AC_DEFUN([ACX_PNFFT], [
AC_REQUIRE([ACX_FFT])
AC_REQUIRE([ACX_PFFT])
acx_pnfft_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(pnfft-prefix, [AS_HELP_STRING([--with-pnfft-prefix=DIR], [http://www-user.tu-chemnitz.de/~mpip/software.php. It requires PFFT.])])
case $with_pnfft_prefix in
  no) acx_pnfft_ok=disable ;;
  *) LIBS_PNFFT="-L$with_pnfft_prefix/lib -lpnfft"; 
     FCFLAGS_PNFFT="$ax_cv_f90_modflag$with_pnfft_prefix/include" ;;
esac


dnl We cannot use PNFFT if MPI is not found
if test "x$acx_mpi_ok" != xyes; then
  acx_pnfft_ok=nompi
fi

dnl We cannot use PNFFT if FFTW3 is not found
if test "x$acx_fft_ok" != xyes; then
  acx_pnfft_ok=nofftw3
fi

dnl We cannot use PNFFT if PFFT is not found
if test "x$acx_pfft_ok" != xyes; then
  acx_pnfft_ok=nopfft
fi

dnl Backup LIBS and FCFLAGS
acx_pnfft_save_LIBS="$LIBS"
acx_pnfft_save_FCFLAGS="$FCFLAGS"


LIBS="$LIBS_PNFFT $LIBS_PFFT $LIBS"
FCFLAGS="$FCFLAGS_PNFFT $FCFLAGS_PFFT $FCFLAGS_FFT $FCFLAGS"


pnfft_func="pnfft_init"


dnl First, check LIBS_PNFFT environment variable
if test x"$acx_pnfft_ok" = xno; then
  AC_MSG_CHECKING([for $pnfft_func in $FCFLAGS_PNFFT $LIBS_PNFFT])
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[  
    use iso_c_binding
    include "fftw3-mpi.f03"
    include "pfft.f03"
    include "pnfft.f03"
    
    call pnfft_init()
  ]), [acx_pnfft_ok=yes], [])
fi


AC_MSG_RESULT([$acx_pnfft_ok ($FCFLAGS_PNFFT $LIBS_PNFFT)])

dnl Generic PNFFT library? is there such a thing?
for pnfft in pnfft; do
  if test x"$acx_pnfft_ok" = xno; then
    AC_CHECK_LIB($pnfft -lpfft -lfftw3_mpi -lfftw3 -lfftw3_mpi -lm, $pnfft_func,
      [acx_pnfft_ok=yes; LIBS_PNFFT="$LIBS_PNFFT -lpfft -lfftw3 -lfftw3_mpi -lm"], [], [$FLIBS])
  fi
done

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
