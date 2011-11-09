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
acx_pfft_ok=no

dnl We cannot use PFFT if MPI is not found
if test "x$acx_mpi_ok" != xyes; then
  acx_pfft_ok=nompi
fi

dnl We cannot use PFFT if FFTW3 is not found
if test "x$acx_fft_ok" != xyes; then
  acx_pfft_ok=nofftw3
fi

dnl Get fortran linker name of PFFT function to check for.
dnl if not compiling with fortran, convert the names
m4_if(_AC_LANG, Fortran, [blacs_pinfo=blacs_pinfo], [AC_FC_FUNC(blacs_pinfo)])

pfft_func="dpfft_plan_dft_3d"

dnl Check if the library was given in the command line
if test $acx_pfft_ok = no; then
  AC_ARG_WITH(pfft, [AS_HELP_STRING([--with-pfft=<lib>], [use PFFT library (http://www-user.tu-chemnitz.de/~mpip/software.php)])])
  case $with_pfft in
    yes | "") ;;
    no) acx_pfft_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_PFFT="$with_pfft" ;;
    *) LIBS_PFFT="-l$with_pfft" ;;
  esac
fi

dnl Backup LIBS 
acx_pfft_save_LIBS="$LIBS"
LIBS="$LIBS_PFFT $LIBS_FFT $LIBS $FLIBS"

dnl First, check LIBS_PFFT environment variable
if test $acx_pfft_ok = no; then
  AC_MSG_CHECKING([for $pfft_func in $LIBS_PFFT])
  AC_TRY_LINK_FUNC($pfft_func, [acx_pfft_ok=yes], [])
  if test $acx_pfft_ok = no; then
    AC_MSG_RESULT([$acx_pfft_ok ($LIBS_PFFT)])
  else
    AC_MSG_RESULT([$acx_pfft_ok ($LIBS_PFFT)])
  fi
fi

dnl Generic PFFT library?
for pfft in pfft; do
  if test $acx_pfft_ok = no; then
    AC_CHECK_LIB($pfft -lfftw3_mpi -lfftw3 -lfftw3_mpi -lm, $pfft_func,
      [acx_pfft_ok=yes; LIBS_PFFT="$LIBS_PFFT -lfftw3 -lfftw3_mpi -lm"], [], [$FLIBS])
  fi
done

AC_SUBST(LIBS_PFFT)
LIBS="$acx_pfft_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pfft_ok" = xyes; then
  AC_DEFINE(HAVE_PFFT,1,[Defined if you have PFFT library.])
  $1
else
  AC_MSG_WARN([Could not find PFFT library. 
               *** Will compile without PFFT support])
  $2
fi
])dnl ACX_PFFT
