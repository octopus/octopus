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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$
##

AC_DEFUN([ACX_FFT],
[
  acx_fft_ok=no

  dnl Check if the library was given in the command line
  dnl Three possibilities to give the libraries:
  dnl 1. ".../../...libfftw3.a .../libfftw3-mpi.a ..." (static linking)
  dnl 2. "-Lpath/to/lib -lfftw3 -lfftw3-mpi ..." (dynamic linking)
  dnl 3. "/prefix/path/to/lib" (determine the subdir from configure)
  if test $acx_fft_ok = no; then
    AC_ARG_WITH(fft-lib, [AS_HELP_STRING([--with-fft-lib=<lib>], [use FFT library <lib>])])
    case $with_fft_lib in
      yes | "") LIBS_FFT=""
                FCFLAGS_FFT=""
		libinput=""
		;;
      -* |  *.so | *.so.*) shared_fft_libs=yes
                           LIBS_FFT="$with_fft_lib"
			   libinput="$with_fft_lib"
			   ;;
      *.a |*.o) shared_fft_libs=no
                LIBS_FFT="$with_fft_lib"
		libinput="$with_fft_lib"
		;;
      *) FFT_PREFIX="$with_fft_lib"
         FCFLAGS_FFT="$ax_cv_f90_modflag$FFT_PREFIX/include"
	 LIBS_FFT="-L$FFT_PREFIX/lib -lfftw3"
	 shared_fft_libs=yes
	 libinput=""
	 ;;
    esac
  fi


  acx_fft_save_FCFLAGS="$FCFLAGS"

  dnl The include dir must be specified when the library is given with a
  dnl specified file to be compiled static (i.e. *.a etc.)
  AC_ARG_WITH(fft-include, [AS_HELP_STRING([--with-fft-include=DIR], [FFTW3 Fortran include files directory])])
  if test "x$with_fft_include" != x; then
    FCFLAGS_FFT="$ax_cv_f90_modflag$with_fft_include"
  else
    if test "x$libinput" != x; then
      dnl split the libinput and use the -L... part
      dnl to determine the library path. Then
      dnl try libpath/../include as include path
      for x in $libinput; do
          case $x in
	    -L*) fft_libpath=$(sed -n "s/-L\(.*\)/\1/p" <<<"$x")
	      ;;
	    */*.a) fft_libpath=$(dirname "$x")
	      ;;
	  esac
      done
      if test "x$fft_libpath" != x; then
         FCFLAGS_FFT="$ax_cv_f90_modflag$fft_libpath/../include"
      fi
    fi
  fi


  testcompile="AC_LANG_PROGRAM([],[
    use, intrinsic :: iso_c_binding
    implicit none

    include 'fftw3.f03'
  ])"

  dnl First, check for the include path
  acx_fft_inc_ok=no
  if test x$acx_fft_inc_ok = xno; then
     AC_MSG_CHECKING([for fftw3.f03])
     AC_COMPILE_IFELSE($testcompile,[acx_fft_inc_ok=yes],[])
     AC_MSG_RESULT([$acx_fft_inc_ok])
     if test x$acx_fft_inc_ok != xyes; then
        FCFLAGS="$FCFLAGS_FFT $FCFLAGS"
	AC_MSG_CHECKING([for fftw3.f03 with $FCFLAGS_FFT])
        AC_COMPILE_IFELSE($testcompile,[acx_fft_inc_ok=yes],[])
	AC_MSG_RESULT([$acx_fft_inc_ok])
     fi
  fi

  if test $acx_fft_inc_ok != yes; then
      AC_MSG_ERROR([Could not find required fftw3.f03 include file.])
  fi

  dnl Now check for the library
  acx_fft_save_LIBS="$LIBS"
  LIBS="$LIBS_FFT $LIBS $FLIBS"

  testlink="AC_LANG_PROGRAM([],[
    use, intrinsic :: iso_c_binding
    implicit none

    include 'fftw3.f03'

    type(C_PTR) :: plan
    integer(C_INT) :: n0
    COMPLEX(C_DOUBLE_COMPLEX) :: in(10),out(10)
    integer(C_INT) :: sign=FFTW_FORWARD
    INTEGER(C_INT) :: flags

    plan = fftw_plan_dft_1d(n0, in, out, sign, flags)
  ])"

  dnl First, check LIBS_FFT environment variable
  if test $acx_fft_ok = no; then
    if test "x$LIBS_FFT" != x; then
      AC_MSG_CHECKING([for fftw_plan_dft_1d in $LIBS_FFT])
      AC_LINK_IFELSE($testlink, [acx_fft_ok=yes], [])
      if test x$acx_fft_ok = xno; then
        AC_MSG_RESULT([$acx_fft_ok])
      else
        AC_MSG_RESULT([$acx_fft_ok ($LIBS_FFT)])
      fi
    fi
  fi

  AC_DEFINE(HAVE_FFTW3, 1, [Define if FFTW3 is available])

  dnl FFTW linked to by default?
  dnl if test $acx_fft_ok = no; then
  dnl   AC_CHECK_FUNC($fft_func, [acx_fft_ok=yes])
  dnl fi

  dnl search libraries
  dnl for fftl in $fft_libs; do
  dnl  if test $acx_fft_ok = no; then
  dnl    AC_CHECK_LIB($fftl, $fft_func,
  dnl      [acx_fft_ok=yes; LIBS_FFT="$LIBS_FFT -l$fftl"], [], [$FLIBS])
  dnl  fi
  dnl done

  testprogram="AC_LANG_PROGRAM([],[
    use, intrinsic :: iso_c_binding
    implicit none

    include 'fftw3.f03'

    INTEGER(kind=C_INT) :: iret

    iret = fftw_init_threads()
  ])"

  AC_MSG_CHECKING([for fftw_init_threads in $LIBS_FFT])
  AC_LINK_IFELSE($testprogram, [acx_fft_threads_ok=yes], [acx_fft_threads_ok=no])
  AC_MSG_RESULT([$acx_fft_threads_ok])
  if test x"$acx_fft_threads_ok" != xyes; then
     if test "x$shared_fft_libs" == xyes; then
        dnl shared libraries
	testlibs="$LIBS_FFT -lfftw3_threads -lfftw3"
     else
        dnl static libraries
	testlibs="$fft_libpath/libfftw3_threads.a $LIBS_FFT"
     fi
     LIBS="$testlibs $LIBS"
     AC_MSG_CHECKING([for fftw_init_threads in $testlibs])
     AC_LINK_IFELSE($testprogram,[acx_fft_threads_ok=yes; LIBS_FFT="$testlibs"],
                                                            [acx_fft_threads_ok=no])
     AC_MSG_RESULT([$acx_fft_threads_ok])
  fi
  if test x"$acx_fft_threads_ok" = xyes; then  
    AC_DEFINE(HAVE_FFTW3_THREADS, 1,[Define if the threaded version of FFTW3 is available.])
  fi


  dnl Check for distributed FFTW3 library
  acx_fft_mpi_ok=no

testprogram="AC_LANG_PROGRAM([],[ 
    use, intrinsic :: iso_c_binding
    implicit none

    include 'fftw3-mpi.f03'

    call MPI_Init
    call fftw_mpi_init
    call MPI_Finalize
  ])"

  AC_MSG_CHECKING([for fftw_mpi_init in $LIBS_FFT])

  if test x"$acx_mpi_ok" != xyes; then
     AC_MSG_RESULT([requires MPI])
     acx_fft_mpi_ok=no
  else
	dnl We have MPI, check the link path
	AC_LINK_IFELSE($testprogram, [acx_fft_mpi_ok=yes], [acx_fft_mpi_ok=no])
  	AC_MSG_RESULT([$acx_fft_mpi_ok])
  	if test x"$acx_fft_mpi_ok" != xyes; then
	   if test "x$shared_fft_libs" == xyes; then
	      dnl shared libraries
	      testlibs="$LIBS_FFT -lfftw3_mpi -lfftw3"
	   else
	      dnl static libraries
	      testlibs="$fft_libpath/libfftw3_mpi.a $LIBS_FFT"
	   fi
	   LIBS="$testlibs $LIBS"
	   AC_MSG_CHECKING([for fftw_mpi_init in $testlibs])
           AC_LINK_IFELSE($testprogram,[acx_fft_mpi_ok=yes; LIBS_FFT="$testlibs"],
                                                            [acx_fft_mpi_ok=no])
    							    AC_MSG_RESULT([$acx_fft_mpi_ok])
        fi

  fi
  if test x"$acx_fft_mpi_ok" = xyes; then  
    AC_DEFINE(HAVE_FFTW3_MPI, 1,[Define if the distributed version of FFTW3 is available.])
  fi

  CFLAGS_FFT="$FCFLAGS_FFT"
  AC_SUBST(CFLAGS_FFT)
  AC_SUBST(FCFLAGS_FFT)
  FCFLAGS="$acx_fft_save_FCFLAGS"
  AC_SUBST(LIBS_FFT)
  LIBS="$acx_fft_save_LIBS"

  # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
  if test x"$acx_fft_ok" != xyes; then
    if test $acx_fft_ok != disable; then
      AC_MSG_ERROR([Could not find required FFT library.])
    fi
    $2
  fi
])
