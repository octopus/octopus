dnl looks for PARPACK library 
AC_DEFUN([ACX_PARPACK], [
AC_REQUIRE([ACX_ARPACK])
acx_parpack_ok=no

dnl We cannot use PARPACK if MPI is not found
if test x$acx_mpi_ok != xyes; then
  acx_parpack_ok=nompi
fi 
dnl We cannot use PARPACK if ARPACKis not found
if test x$acx_arpack_ok != xyes; then
  acx_parpack_ok=noarpack
fi


dnl Backup LIBS 
acx_parpack_save_LIBS="$LIBS"


dnl Check if the library was given in the command line
if test $acx_parpack_ok = no; then
  AC_ARG_WITH(parpack, [AS_HELP_STRING([--with-parpack=<lib>], [use PARPACK library <lib>. Requires ARPACK.])])
  case $with_parpack in
    yes | "") ;;
    no) acx_parpack_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_PARPACK="$with_parpack" ;;
    *) LIBS_PARPACK="-l$with_parpack" ;;
  esac
fi

dnl PARPACK always wants ARPACK 
LIBS_PARPACK="$LIBS_PARPACK $LIBS_ARPACK"

dnl First, check LIBS_PARPACK environment variable
if test $acx_parpack_ok = no; then
  LIBS="$LIBS_PARPACK $LIBS_LAPACK $LIBS_BLAS $acx_parpack_save_LIBS $FLIBS"
  AC_MSG_CHECKING([for parpack library])
  AC_LINK_IFELSE([
	    program main
	    call pdsaupd
	    end program main], [acx_parpack_ok=yes], [])
  if test $acx_parpack_ok = no; then
    AC_MSG_RESULT([$acx_parpack_ok])
  else
    AC_MSG_RESULT([$acx_parpack_ok ($LIBS_PARPACK)])
  fi
fi

if test $acx_parpack_ok = no; then
  AC_MSG_CHECKING([for parpack library with -lparpack])
  if test "$LIBS_PARPACK" = ""; then
    LIBS="-lparpack  $LIBS_LAPACK $LIBS_BLAS $acx_parpack_save_LIBS $FLIBS" 
    AC_LINK_IFELSE([
      program main
      call pdsaupd
      end program main
    ], [acx_parpack_ok=yes; LIBS_PARPACK=" -lparpack"], [])
  else
    LIBS="-L$LIBS_PARPACK -lparpack  $LIBS_LAPACK $LIBS_BLAS $acx_parpack_save_LIBS $FLIBS" 
    AC_LINK_IFELSE([
      program main
      call pdsaupd
      end program main
    ], [acx_parpack_ok=yes; LIBS_PARPACK="-L$LIBS_PARPACK -lparpack"], [])  
  fi
  if test $acx_parpack_ok = no; then
    AC_MSG_RESULT([$acx_parpack_ok])
  else
    AC_MSG_RESULT([$acx_parpack_ok ($LIBS_PARPACK)])
  fi
fi


AC_SUBST(LIBS_PARPACK)
LIBS="$acx_parpack_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_parpack_ok" = xyes; then
  AC_DEFINE(HAVE_PARPACK,1,[Defined if you have PARPACK library.])
  $1
else
  if test $acx_parpack_ok = nompi; then
    AC_MSG_WARN([PARPACK requires MPI support.]) 
  fi
  if test $acx_parpack_ok = noarpack; then
    AC_MSG_WARN([PARPACK requires ARPACK.]) 
  fi   
  AC_MSG_WARN([Could not find PARPACK library. 
              *** Will compile without PARPACK support])
  $2
fi


])dnl ACX_PARPACK
