## Copyright (C) 2011 J. Alberdi-Rodriguez
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
## $Id: pfft.m4 6722 2010-06-13 12:44:43Z joseba $
##

AC_DEFUN([ACX_PFFT], [
AC_REQUIRE([ACX_FFT])
acx_pfft_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(pfft-prefix, [AS_HELP_STRING([--with-pfft-prefix=<lib>], [PFFT version 1.0.7 is required, linked with a patched version of FFTW3.3.3. Source and more information at http://www-user.tu-chemnitz.de/~mpip/software.php])])
case $with_pfft_prefix in
  yes | "") ;;
  no) acx_pfft_ok=disable ;;
  *.a | *.so | *.so.* | *.o) LIBS_PFFT=$with_pfft_prefix ;
     xpath=${with_pfft_prefix%/lib/*} 
     FCFLAGS_PFFT="$ax_cv_f90_modflag$xpath/include";;  
  *) LIBS_PFFT="-L$with_pfft_prefix/lib"; 
     FCFLAGS_PFFT="$ax_cv_f90_modflag$with_pfft_prefix/include" ;;
esac

dnl The include dir must be specified when the library is given with a 
dnl specified file to be compiled static (i.e. *.a etc.)
AC_ARG_WITH(pfft-include, [AS_HELP_STRING([--with-pfft-include=DIR], [PFFT Fortran include files directory])])
case $with_pfft_include in
  "") if test "x$FCFLAGS_PFFT" == x; then
        FCFLAGS_PFFT="$ax_cv_f90_modflag/usr/include"
      fi;;
  *)  FCFLAGS_PFFT="$ax_cv_f90_modflag$with_pfft_include" ;;
esac


dnl We need to link against the MPI FFTW3 used to compile PFFT 
AC_ARG_WITH(mpifftw-prefix, [AS_HELP_STRING([--with-mpifftw-prefix=DIR], [MPI FFTW3 libraries directory (the one used to build PFFT).])])
case $with_mpifftw_prefix in
 "")  case $LIBS_FFT in 
      *.a | *.so | *.so.* | *.o) xpath=${LIBS_FFT%/lib/*};
         LIBS_MPIFFT="-L$xpath/lib -lfftw3_mpi"; 
         FCFLAGS_MPIFFT="$ax_cv_f90_modflag$xpath/include";;
      *) LIBS_MPIFFT="-L/usr/lib -lfftw3_mpi"; 
         FCFLAGS_MPIFFT="$ax_cv_f90_modflag/usr/include";;
      esac
      ;;
 *.a | *.so | *.so.* | *.o) LIBS_MPIFFT=$with_mpifftw_prefix;
     xpath=${with_mpifftw_prefix%/lib/*} 
     FCFLAGS_MPIFFT="$ax_cv_f90_modflag$xpath/include";;
 *)  LIBS_MPIFFT="-L$with_mpifftw_prefix/lib -lfftw3_mpi -lfftw3";
     FCFLAGS_MPIFFT="$ax_cv_f90_modflag$with_mpifftw_prefix/include" ;;
esac


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

FCFLAGS_PFFT="$FCFLAGS_PFFT $FCFLAGS_MPIFFT"
FCFLAGS="$FCFLAGS_PFFT $acx_pfft_save_FCFLAGS"

# some symbols below will not be defined for version 1.0.4, making sure
# we have a version that is able to work in our code
testprogram="AC_LANG_PROGRAM([],[ 
    implicit none

    include 'fftw3.f'
    include 'pfft.f'

    integer                 :: comm_cart_2d
    integer(ptrdiff_t_kind) :: n(3)
    integer(8)              :: fftplan
    complex(8), allocatable :: data_in(:)
    complex(8), allocatable :: data_out(:)

    call dpfft_plan_dft_3d(fftplan, n, data_in, data_out, comm_cart_2d, &
     &     PFFT_FORWARD, PFFT_TRANSPOSED_OUT + PFFT_MEASURE + PFFT_DESTROY_INPUT)

  ])"


dnl First, check LIBS_PFFT environment variable
if test x"$acx_pfft_ok" = xno; then
  LIBS="$LIBS_PFFT $LIBS_MPIFFT $LIBS_FFT $acx_pfft_save_LIB"
  AC_MSG_CHECKING([for pfft library])
  AC_LINK_IFELSE($testprogram, [acx_pfft_ok=yes; LIBS_PFFT="$LIBS_PFFT $LIBS_MPIFFT $LIBS_FFT"], [])
  if test $acx_pfft_ok = no; then
    AC_MSG_RESULT([$acx_pfft_ok])
  else
    AC_MSG_RESULT([$acx_pfft_ok ($LIBS_PFFT)])
  fi
fi


dnl Generic PFFT library 
if test $acx_pfft_ok = no; then
  AC_MSG_CHECKING([for pfft library with -lpfft])
  if test "$LIBS_PFFT" = ""; then
    LIBS="-lpfft $LIBS_MPIFFT $LIBS_FFT $LIBS $acx_pfft_save_LIB"
    AC_LINK_IFELSE($testprogram, [acx_pfft_ok=yes; LIBS_PFFT="-lpfft $LIBS_MPIFFT $LIBS_FFT"], [])
  else
    LIBS="$LIBS_PFFT -lpfft $LIBS_MPIFFT $LIBS_FFT $acx_pfft_save_LIB"
    AC_LINK_IFELSE($testprogram, [acx_pfft_ok=yes; 
                                  LIBS_PFFT="$LIBS_PFFT -lpfft $LIBS_MPIFFT $LIBS_FFT"], [])  
  fi
  if test $acx_pfft_ok = no; then
    AC_MSG_RESULT([$acx_pfft_ok])
  else
    AC_MSG_RESULT([$acx_pfft_ok ($LIBS_PFFT)])
  fi
fi


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
