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

# Check size of a pointer
AC_DEFUN(ACX_POINTER_SIZE,
[AC_MSG_CHECKING([for the size of a pointer])
AC_REQUIRE([AC_PROG_CC])
if test -z "$POINTER_SIZE"; then
cat >pointertest.c <<EOF
#include <stdio.h>
void main()
{
  printf("%ld", sizeof(void *));
}
EOF
	ac_try='$CC $CFLAGS -o pointertest.x pointertest.c 1>&AC_FD_CC'
	if AC_TRY_EVAL(ac_try); then
  	ac_try=""
	else
  	echo "configure: failed program was:" >&AC_FD_CC
	  cat pointertest.c >&AC_FD_CC
  	rm -f pointertest*
	  AC_MSG_ERROR(failed to compile c program to find the size of a pointer)
	fi
	ac_pointersize=`./pointertest.x`;
	rm -f pointertest*
	AC_DEFINE_UNQUOTED(POINTER_SIZE, ${ac_pointersize}, [The size of a C pointer])
	AC_MSG_RESULT([${ac_pointersize} bytes])
fi
])

AC_DEFUN([ACX_CHECK_FUNC],
[
	if test -z "${$1}"; then
		# the space in " $3" is needed.
		AC_CHECK_FUNC([$2], [$1=" $3"], [$4], [$5])
  fi
])

AC_DEFUN([ACX_CHECK_LIB],
[
	if test -z "${$1}"; then
		AC_CHECK_LIB([$2], [$3], [$1="$4"], [$5], [$6])
	fi
])

AC_DEFUN([ACX_LIB_FFTW],
[
AC_ARG_ENABLE(fft, [  --enable-fft=ARG    enable fft support
      ARG=no          (the same as --disable-fft) disable FFT support
      ARG=fftw3       FFTW3 support through libfftw3 (default)
      ARG=fftw2       FFTW2 support through libfftw],
  [fft=$enableval], [fft=fftw3])

case "${fft}" in
  fftw3*)
    if test "${SINGLE_PRECISION}"; then
   echo HHHHHHHHHHHHHHHHHHHHHHHHHHHHh
      ACX_CHECK_FUNC([LIB_FFT], [sfftw_plan_dft_1d])
      ACX_CHECK_LIB([LIB_FFT], [fftw3f],  [sfftw_plan_dft_1d], [-lfftw3f])
    else
      ACX_CHECK_FUNC([LIB_FFT], [dfftw_plan_dft_1d])
      ACX_CHECK_LIB([LIB_FFT], [fftw3],  [dfftw_plan_dft_1d], [-lfftw3])
    fi

    AS_IF([test -z "${LIB_FFT}"], AC_MSG_ERROR([could not find fftw3 library]))
    fft=3
    ;;
  fftw2*)
      ACX_CHECK_FUNC([LIB_FFT], [fftw3d_f77_create_plan])
      ACX_CHECK_LIB([LIB_FFT], [fftw],  [fftw3d_f77_create_plan], [-lfftw])
      ACX_CHECK_LIB([LIB_FFT], [dfftw], [fftw3d_f77_create_plan], [-ldfftw])
      AS_IF([test -z "${LIB_FFT}"], AC_MSG_ERROR([could not find fftw library]))
	
      ACX_CHECK_FUNC([lib_rfft], [rfftw3d_f77_create_plan], [], [], [${LIB_FFT}])
      ACX_CHECK_LIB([lib_rfft], [rfftw], [rfftw3d_f77_create_plan], [-lrfftw], [], [${LIB_FFT}])
      ACX_CHECK_LIB([lib_rfft], [drfftw], [rfftw3d_f77_create_plan], [-ldrfftw], [], [${LIB_FFT}])
      AS_IF([test -z "${lib_rfft}"], AC_MSG_ERROR([could not find rfftw library]))

      LIB_FFT="${lib_rfft} ${LIB_FFT}"
      fft=2
      ;;
    *)
      fft=no
      ;;
  esac

if test -z "${LIB_FFT}"; then
  AC_MSG_WARN([Octopus will be compiled *without* fft support])
else
  AC_DEFINE_UNQUOTED(HAVE_FFT, $fft, [FFT library (fftw2 | fftw3)])
  AC_SUBST(LIB_FFT)
fi
])
