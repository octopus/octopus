## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

AC_DEFUN([ACX_FFTW],
[
acx_fftw_ok=no

AC_ARG_WITH(fft, [  --with-fft=ARG    fft support
      ARG=no          disable FFT support
      ARG=fftw3       FFTW3 support through libfftw3 (default)
      ARG=fftw2       FFTW2 support through libfftw],
  [fft=$withval], [fft=fftw3])

case $fft in
  fftw3*)
    fft=3
    if test "x${SINGLE_PRECISION}" != x; then
      fftw_func="sfftw_plan_dft_1d"
      fftw_libs="fftw3f"
    else
      fftw_func="dfftw_plan_dft_1d"
      fftw_libs="fftw3"
    fi
    ;;
  fftw2*)
    fft=2
    if test "x${SINGLE_PRECISION}" != x; then
      AC_MSG_ERROR([Version 2 of FFTW can not be used with single precision])
    else
      fftw_func="fftw3d_f77_create_plan"
      fftw_libs="fftw dfftw"
    fi
    ;;
  *)
    fft=no
    ;;
esac

if test "${fft}" != "no"; then

  dnl Check if the library was given in the command line
  if test $acx_fftw_ok = no; then
    AC_ARG_WITH(fft-lib, [AC_HELP_STRING([--with-fft-lib=<lib>], [use FFT library <lib>])])
    case $with_fft_lib in
      yes | "") ;;
      no) acx_fftw_ok=disable ;;
      -* | */* | *.a | *.so | *.so.* | *.o) LIBS_FFT="$with_fft_lib" ;;
      *) LIBS_FFT="-l$with_fft_lib" ;;
    esac
  fi

  acx_fftw_save_LIBS="$LIBS"
  LIBS="$LIBS_FFT $LIBS $FLIBS"

  dnl First, check LIBS_FFTW environment variable
  if test $acx_fftw_ok = no; then
    if test "x$LIBS_FFTW" != x; then
      AC_MSG_CHECKING([for $fftw_func in $LIBS_FFTW])
      AC_TRY_LINK_FUNC($fftw_func, [acx_fftw_ok=yes], [LIBS_FFTW=""])
      AC_MSG_RESULT($acx_fftw_ok)
    fi
  fi

  dnl FFTW linked to by default?
  if test $acx_fftw_ok = no; then
    AC_CHECK_FUNC($fftw_func, [acx_fftw_ok=yes])
  fi

  dnl search libraries
  for fftw in $fftw_libs; do
    if test $acx_fftw_ok = no; then
      AC_CHECK_LIB($fftw, $fftw_func,
        [acx_fftw_ok=yes; LIBS_FFTW="$LIBS_FFTW -l$fftw"], [], [$FLIBS])
    fi
  done

  dnl if we have fftw2, search also for the real transforms
  if test "${fft}" = "2"; then
    if test $acx_fftw_ok = yes; then

      acx_fftw_ok=no;

      AC_CHECK_FUNC(rfftw3d_f77_create_plan, [acx_fftw_ok=yes])

      for fftw in rfftw drfftw; do
        if test $acx_fftw_ok = no; then
          AC_CHECK_LIB($fftw, rfftw3d_f77_create_plan,
            [acx_fftw_ok=yes; LIBS_FFTW="$LIBS_FFTW -l$fftw"], [], [$LIBS_FFTW $FLIBS])
        fi
      done
    fi
  fi

fi

AC_SUBST(LIBS_FFTW)
LIBS="$acx_fftw_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_fftw_ok" = xyes; then
  AC_DEFINE_UNQUOTED(HAVE_FFT, [$fft], [FFT library (fftw2 | fftw3)])
  $1
else
  if test $acx_fftw_ok != disable; then
    AC_MSG_WARN([Could not find fftw library. 
                *** Will compile without fftw support])
  fi
  $2
fi
])
