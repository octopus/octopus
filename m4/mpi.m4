AC_DEFUN([ACX_MPI], [
acx_mpi_ok=no

dnl Backup LIBS 
acx_mpi_save_LIBS="$LIBS"
LIBS="$LIBS_MPI $LIBS $FLIBS"

dnl First, check LIBS_MPI environment variable
if test $acx_mpi_ok = no; then
  AC_MSG_CHECKING([for MPI_init in $LIBS_MPI])
  AC_TRY_LINK_FUNC(MPI_Init, [acx_mpi_ok=yes], [])
  if test $acx_mpi_ok = no; then
    AC_MSG_RESULT([$acx_mpi_ok])
  else
    AC_MSG_RESULT([$acx_mpi_ok ($LIBS_MPI)])
  fi
fi

if test $acx_mpi_ok = no; then
  AC_CHECK_LIB(mpi, MPI_Init, [acx_mpi_ok=yes; LIBS_MPI="$LIBS_MPI -lmpi"])
fi

AC_SUBST(LIBS_MPI)
LIBS="$acx_mpi_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_mpi_ok" = xyes; then
  AC_DEFINE(HAVE_MPI,1,[Defined if you have MPI library.])
  $1
else
  $2
fi
])dnl ACX_MPI

AC_DEFUN([ACX_MPI_FC_MODULE], [
dnl let us see if we have a mpi module
AC_MSG_CHECKING([for MPI Fortran headers])
save_ldflags="$LDFLAGS"
AS_IF([test "$LIB_MPI"], [LDFLAGS="${LDFLAGS} -L${LIB_MPI}"])

AC_COMPILE_IFELSE(AC_LANG_PROGRAM([], [
include 'mpif.h'
integer :: ierr
call MPI_Init(ierr)
]), [HAVE_MPIF_H=1], [HAVE_MPIF_H=0])

if test "$HAVE_MPIF_H" = 1; then
  AC_DEFINE(MPI_H, 1, [have MPI Fortran header file])
  AC_MSG_RESULT([mpif.h])
else
  AC_COMPILE_IFELSE(AC_LANG_PROGRAM([], [
  use mpi
  integer	:: ierr
  call MPI_Init(ierr)
  ]), [HAVE_MPI_MOD=1], [HAVE_MPI_MOD=0])

  if test "$HAVE_MPI_MOD" = 1; then
    AC_DEFINE(MPI_MOD, 1, [have mpi module])
    AC_MSG_RESULT([mpi module])
  else
    AC_MSG_ERROR([Could neither find the mpi module nor mpif.h.])
  fi
fi
])dnl ACX_MPI_FC_MODULE

AC_DEFUN([ACX_MPI2], [
acx_mpi2_ok=no

AC_MSG_CHECKING([for MPI 2])

AC_ARG_ENABLE(mpi2, AS_HELP_STRING([--disable-mpi2], [Disable the usage of MPI 2 routines.]))

if test x"$enable_mpi2" != x"no"; then

if test "$HAVE_MPIF_H" = 1; then
AC_COMPILE_IFELSE(AC_LANG_PROGRAM([], [[
implicit none
include 'mpif.h'
integer :: aa, ierr
call MPI_Allreduce(MPI_IN_PLACE, aa, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
]]), [acx_mpi2_ok=yes], [])
else
AC_COMPILE_IFELSE(AC_LANG_PROGRAM([], [[
use mpi
implicit none
integer :: aa, ierr
call MPI_Allreduce(MPI_IN_PLACE, aa, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
]]), [acx_mpi2_ok=yes], [])

fi

if test $acx_mpi2_ok = yes; then
  AC_DEFINE(HAVE_MPI2, 1, [Defined if you have an MPI 2 implementation])
fi
AC_MSG_RESULT([$acx_mpi2_ok])

else

AC_MSG_RESULT([disabled])

fi

])

