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
## $Id$
##

AC_DEFUN([ACX_FFT],
[
  acx_fft_ok=no

  if test "x${SINGLE_PRECISION}" != x; then
    fft_func="sfftw_plan_dft_1d"
    fft_lib="fftw3f"
  else
    fft_func="dfftw_plan_dft_1d"
    fft_lib="fftw3"
  fi

  dnl Check if the library was given in the command line
  AC_ARG_WITH(fft-lib, [AS_HELP_STRING([--with-fft-lib=DIR], [use FFT library prefix directory.])])
  case $with_fft_lib in
    yes | "") ;;
    *.a | *.so | *.so.* | *.o) LIBS_FFT="$with_fft_lib" ;;
    *) LIBS_FFT="-L$with_fft_lib/lib";
       FCFLAGS_FFT="$ax_cv_f90_modflag$with_fft_lib/include" ;;
  esac

  
  dnl The include dir must be specified when the library is given with a 
  dnl specified file to be compiled static (i.e. *.a etc.)
  AC_ARG_WITH(fft-include, [AS_HELP_STRING([--with-fft-include=DIR], [FFT library include directory.])])
  case $with_fft_include in
    "") if test "x$FCFLAGS_FFT" == x; then
          FCFLAGS_FFT="$ax_cv_f90_modflag/usr/include"
        fi;;
    *)  FCFLAGS_FFT="$ax_cv_f90_modflag$with_fft_include" ;;
  esac

  acx_fft_save_LIBS="$LIBS"
  acx_fft_save_FCFLAGS="$FCFLAGS"

  testprogram="AC_LANG_PROGRAM([],[ 
      call $fft_func
    ])"


  dnl First, check LIBS_FFT lib
  if test $acx_fft_ok = no; then
    LIBS="$LIBS_FFT $LIBS $FLIBS"
    AC_MSG_CHECKING([for $fft_func in $LIBS_FFT])
    AC_LINK_IFELSE($testprogram, [acx_fft_ok=yes], [])    

    if test $acx_fft_ok = no; then
      AC_MSG_RESULT([$acx_fft_ok])
    else
      AC_MSG_RESULT([$acx_fft_ok ($LIBS_FFT)])
    fi
    
  fi

  dnl FFTW linked to by default?
  if test $acx_fft_ok = no; then
    AC_CHECK_FUNC($fft_func, [acx_fft_ok=yes])
  fi

  dnl search libraries
  if test $acx_fft_ok = no; then
    AC_MSG_CHECKING([for fft library with -l$fft_lib])
    if test "$LIBS_FFT" = ""; then
      LIBS="-l$fft_lib $LIBS $acx_fft_save_LIB"
      AC_LINK_IFELSE($testprogram, [acx_fft_ok=yes; 
                                    LIBS_FFT="-l$fft_lib"], [])
    else
      LIBS="$LIBS_FFT -l$fft_lib $acx_fft_save_LIB"
      AC_LINK_IFELSE($testprogram, [acx_fft_ok=yes; 
                                    LIBS_FFT="$LIBS_FFT -l$fft_lib"], [])  

    fi
    if test $acx_fft_ok = no; then
      AC_MSG_RESULT([$acx_fft_ok])
    else
      AC_MSG_RESULT([$acx_fft_ok ($LIBS_FFT)])
    fi
  fi
  

  if test $acx_fft_ok = yes; then
    AC_CHECK_FUNC(dfftw_init_threads, 
                  AC_DEFINE(HAVE_FFTW3_THREADS, 1,[Define if the threaded version of FFTW3 is available.]))
  fi

  dnl check if the MPI version of FFTW3 is also available
  if test $acx_fft_ok = yes; then

    acx_fftmpi_ok=no
    
    testprogram="AC_LANG_PROGRAM([],[ 
        use iso_c_binding
        include 'fftw3-mpi.f03'

        call fftw_mpi_init()
      ])"
      
    FCFLAGS="$FCFLAGS_FFT $FCFLAGS"
    
    AC_MSG_CHECKING([for fftw_mpi_init])
    
    if test $acx_fftmpi_ok = no; then
      LIBS="$LIBS_FFT $LIBS $FLIBS"
      AC_LINK_IFELSE($testprogram, [acx_fftmpi_ok=yes; LIBS_FFTMPI=" "], [])
    fi
            
    if test $acx_fftmpi_ok = no; then
      LIBS="-lfftw3_mpi $LIBS_FFT $LIBS $FLIBS"
      AC_LINK_IFELSE($testprogram, [acx_fftmpi_ok=yes; LIBS_FFTMPI="-lfftw3_mpi"], [])
    fi
    
    if test $acx_fftmpi_ok = no; then
      AC_MSG_RESULT([$acx_fftmpi_ok])
    else
      AC_MSG_RESULT([$acx_fftmpi_ok ($LIBS_FFTMPI)])
    fi
  
  fi

  AC_SUBST(LIBS_FFT)
  AC_SUBST(FCFLAGS_FFT)
  LIBS="$acx_fft_save_LIBS"
  FCFLAGS="$acx_fft_save_FCFLAGS"

  # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
  if test x"$acx_fft_ok" != xyes; then
    if test $acx_fft_ok != disable; then
      AC_MSG_ERROR([Could not find required FFT library.])
    fi
    $2
  fi
])
