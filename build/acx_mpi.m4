AC_DEFUN([ACX_MPI], [
acx_mpi_ok=no

dnl Backup LIBS 
acx_mpi_save_LIBS="$LIBS"
LIBS="$LIBS_MPI $LIBS $FLIBS"

dnl First, check LIBS_MPI environment variable
if test $acx_mpi_ok = no; then
  AC_MSG_CHECKING([for MPI_init in $LIBS_MPI])
  AC_TRY_LINK_FUNC(MPI_init, [acx_mpi_ok=yes], [])
  AC_MSG_RESULT($acx_mpi_ok)
fi

if test $acx_mpi_ok = no; then
  AC_CHECK_LIB(mpi, MPI_init, [acx_mpi_ok=yes; LIBS_MPI="$LIBS_MPI -l$mpi"])
fi

dnl let us see if we have a mpi module
save_ldflags="$LDFLAGS"
AS_IF([test "$LIB_MPI"], [LDFLAGS="${LDFLAGS} -L${LIB_MPI}"])
AC_COMPILE_IFELSE([
use mpi
integer :: ierr
call MPI_Init(ierr)
], [HAVE_MPI_MOD=1], [HAVE_MPI_MOD=0])

if test "$HAVE_MPI_MOD" = 1; then
  AC_DEFINE(MPI_MOD, 1, [have mpi module])
else
  AC_COMPILE_IFELSE([
include mpif.h
integer :: ierr
call MPI_Init(ierr)
  ], [HAVE_MPIF_H=1], [HAVE_MPIF_H=0])

  if test "$HAVE_MPIF_H" = 1; then
    AC_DEFINE(MPI_H, 1, [have mpi Fortran header file])
  else
    AC_MSG_WARN([Could not find neither the mpi module or mpif.h. 
                *** Continue at your own risk ;)])
  fi
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
