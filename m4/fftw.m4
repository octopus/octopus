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
##

AC_DEFUN([ACX_FFTW],
[

acx_fftw_ok=no

dnl FIRST CHECK FOR THE BASE FFTW
AC_MSG_CHECKING([for fftw])

dnl Check if the library was given in the command line
AC_ARG_WITH(fftw-prefix, [AS_HELP_STRING([--with-fftw-prefix=DIR], [Directory where fftw was installed.])])


fftw_program="AC_LANG_PROGRAM([],[
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  type(C_PTR) :: plan
  integer(C_INT) :: n0 = 32
  COMPLEX(C_DOUBLE_COMPLEX) :: in(10),out(10)
  integer(C_INT) :: sign=FFTW_FORWARD
  INTEGER(C_INT) :: flags = 0

  plan = fftw_plan_dft_1d(n0, in, out, sign, flags)
])"

fftw_threads_program="AC_LANG_PROGRAM([],[
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  INTEGER(kind=C_INT) :: iret

  iret = fftw_init_threads()
])"

fftw_mpi_program="AC_LANG_PROGRAM([],[
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3-mpi.f03'

  call MPI_Init
  call fftw_mpi_init
  call MPI_Finalize
])"

dnl use LAPACK/BLAS libs
acx_fftw_save_LIBS="$LIBS"
acx_fftw_save_FCFLAGS="$FCFLAGS"
LIBS="$LIBS_LAPACK $LIBS_BLAS $LIBS $FLIBS"
if test x"$FCFLAGS_FFTW" = x; then
  dnl FCFLAGS_FFTW not specified
  dnl check if --with-fftw-prefix given
  if test x$with_fftw_prefix != x; then
    FCFLAGS="-I$with_fftw_prefix/include $acx_fftw_save_FCFLAGS"
    AC_LINK_IFELSE($fftw_program, [acx_fftw_ok=yes], [acx_fftw_ok=no])
    AC_LINK_IFELSE($fftw_threads_program, [acx_fftw_threads_ok=yes], [acx_fftw_threads_ok=no])
    if test x$acx_mpi_ok = xyes; then
      AC_LINK_IFELSE($fftw_mpi_program, [acx_fftw_mpi_ok=yes], [acx_fftw_mpi_ok=no])
    fi
    if test x"$acx_fftw_ok" = xyes; then
      AC_MSG_ERROR([--with-fftw-prefix given, but linking with BLAS/LAPACK successful! The specified FFTW would not be used. Please specify only the header location in FCFLAGS_FFTW.])
    fi
    FCFLAGS="$acx_fftw_save_FCFLAGS"
  fi
else
  dnl also use FCFLAGS_FFTW if given
  FCFLAGS="$FCFLAGS_FFTW $acx_fftw_save_FCFLAGS"
  AC_LINK_IFELSE($fftw_program, [acx_fftw_ok=yes], [acx_fftw_ok=no])
  AC_LINK_IFELSE($fftw_threads_program, [acx_fftw_threads_ok=yes], [acx_fftw_threads_ok=no])
  if test x$acx_mpi_ok = xyes; then
    AC_LINK_IFELSE($fftw_mpi_program, [acx_fftw_mpi_ok=yes], [acx_fftw_mpi_ok=no])
  fi
  dnl now complain if LIBS_FFTW given
  if test x"$acx_fftw_ok" = xyes; then
    if test x"$LIBS_FFTW" != x; then
      AC_MSG_ERROR([LIBS_FFTW given, but linking with BLAS/LAPACK successful! The specified FFTW would not be used. Please specify only the header location in FCFLAGS_FFTW.])
    fi
  else
    dnl now try with LIBS_FFTW instead of BLAS/LAPACK
    LIBS="$LIBS_FFTW"
    AC_LINK_IFELSE($fftw_program, [acx_fftw_ok=yes], [acx_fftw_ok=no])
    AC_LINK_IFELSE($fftw_threads_program, [acx_fftw_threads_ok=yes], [acx_fftw_threads_ok=no])
    if test x$acx_mpi_ok = xyes; then
      AC_LINK_IFELSE($fftw_mpi_program, [acx_fftw_mpi_ok=yes], [acx_fftw_mpi_ok=no])
    fi
  fi
fi

dnl reset LIBS and FCFLAGS
LIBS="$acx_fftw_save_LIBS"
FCFLAGS="$acx_fftw_save_FCFLAGS"


if test x"$acx_fftw_ok" = xyes; then
  dnl if fftw was already found, add the FCFLAGS/CFLAGS/LIBS if given
  AC_SUBST(FCFLAGS_FFTW)
  CFLAGS_FFTW="$FCFLAGS_FFTW"
  AC_SUBST(CFLAGS_FFTW)
  AC_SUBST(LIBS_FFTW)

  AC_MSG_RESULT([$acx_fftw_ok (no extra flags needed, provided by other library, e.g. MKL, or environment variables (FCFLAGS_FFTW and/or LIBS_FFTW))])
  AC_DEFINE(HAVE_FFTW3, 1, [Define if FFTW3 is available])

  AC_MSG_CHECKING([for fftw multithread support])
  AC_MSG_RESULT([$acx_fftw_threads_ok])
  if test x"$acx_fftw_threads_ok" != xno; then
    AC_DEFINE(HAVE_FFTW3_THREADS, 1,[Define if the threaded version of FFTW3 is available.])
  fi

  dnl NOW CHECK WHETHER THE MPI VERSION IS AVAILABLE
  if test x$acx_mpi_ok = xyes; then
    AC_MSG_CHECKING([whether fftw has MPI support])
    AC_MSG_RESULT([$acx_fftw_mpi_ok])
    if test x"$acx_fftw_mpi_ok" = xyes; then
      AC_DEFINE(HAVE_FFTW3_MPI, 1,[Define if the distributed version of FFTW3 is available.])
    fi
  fi

else
  dnl BACKUP LIBS AND FCFLAGS
  acx_fftw_save_LIBS="$LIBS"
  acx_fftw_save_FCFLAGS="$FCFLAGS"

  # Set FCFLAGS_FFTW only if not set from environment
  if test x"$FCFLAGS_FFTW" = x; then
    case $with_fftw_prefix in
      "") FCFLAGS_FFTW="-I/usr/include" ;;
      *)  FCFLAGS_FFTW="-I$with_fftw_prefix/include" ;;
    esac
  fi

  FCFLAGS="$FCFLAGS_FFTW $acx_fftw_save_FCFLAGS"

  if test ! -z "$with_fftw_prefix"; then
    LIBS_FFTW="-L$with_fftw_prefix/lib"
  else
    LIBS_FFTW=""
  fi

  dnl We do not append -lfftw3 at the end, as we might need to prefix other libraries
  LIBS="$LIBS_FFTW -lfftw3 $acx_fftw_save_LIBS"
  AC_LINK_IFELSE($fftw_program, [acx_fftw_ok=yes], [acx_fftw_ok=no])
  AC_MSG_RESULT([$acx_fftw_ok ($FCFLAGS_FFTW $LIBS_FFTW -lfftw3)])

  if test x$acx_fftw_ok != xyes; then
    AC_MSG_ERROR([Could not find required fftw library])
  fi

  AC_DEFINE(HAVE_FFTW3, 1, [Define if FFTW3 is available])

  dnl NOW CHECK WHETHER THE MULTITHREADED VERSION IS AVALIABLE
  AC_MSG_CHECKING([for fftw multithread support])
  dnl First try the OpenMP version
  LIBS="$LIBS_FFTW -lfftw3_omp -lfftw3 $acx_fftw_save_LIBS"
  AC_LINK_IFELSE($fftw_threads_program, [acx_fftw_threads_ok=fftw3_omp], [acx_fftw_threads_ok=no])

  if test x"$acx_fftw_threads_ok" == xfftw3_omp; then
    LIBS_FFTW="$LIBS_FFTW -lfftw3_omp"
  fi

  dnl Second try, the pthreads version
  if test x"$acx_fftw_threads_ok" == xno; then
    LIBS="$LIBS_FFTW -lfftw3_threads -lfftw3 $acx_fftw_save_LIBS"
    AC_LINK_IFELSE($fftw_threads_program, [acx_fftw_threads_ok=fftw3_threads], [acx_fftw_threads_ok=no])
  fi

  if test x"$acx_fftw_threads_ok" == xfftw3_threads; then
    LIBS_FFTW="$LIBS_FFTW -lfftw3_threads"
  fi

  AC_MSG_RESULT([$acx_fftw_threads_ok])
  if test x"$acx_fftw_threads_ok" != xno; then
    AC_DEFINE(HAVE_FFTW3_THREADS, 1,[Define if the threaded version of FFTW3 is available.])
  fi

  dnl NOW CHECK WHETHER THE MPI VERSION IS AVAILABLE
  if test x$acx_mpi_ok = xyes; then
    AC_MSG_CHECKING([whether fftw has MPI support])
    LIBS="$LIBS_FFTW -lfftw3_mpi -lfftw3 $acx_fftw_save_LIBS"
    AC_LINK_IFELSE($fftw_mpi_program, [acx_fftw_mpi_ok=yes], [acx_fftw_mpi_ok=no])
    AC_MSG_RESULT([$acx_fftw_mpi_ok])

    if test x"$acx_fftw_mpi_ok" = xyes; then
      LIBS_FFTW="$LIBS_FFTW -lfftw3_mpi"
      AC_DEFINE(HAVE_FFTW3_MPI, 1,[Define if the distributed version of FFTW3 is available.])
    fi
  fi

  dnl now we append -lfftw3
  LIBS_FFTW="$LIBS_FFTW -lfftw3"

  CFLAGS_FFTW="$FCFLAGS_FFTW"
  AC_SUBST(FCFLAGS_FFTW)
  AC_SUBST(CFLAGS_FFTW)
  AC_SUBST(LIBS_FFTW)

  FCFLAGS="$acx_fftw_save_FCFLAGS"
  LIBS="$acx_fftw_save_LIBS"
fi

])
