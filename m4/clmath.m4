AC_DEFUN([ACX_CLFFT], [
acx_clfft_ok=no

dnl Backup LIBS and CFLAGS
acx_clfft_save_LIBS="$LIBS"
acx_clfft_save_CFLAGS="$CFLAGS"

dnl Check if the library was given in the command line
AC_ARG_WITH(clfft-prefix, [AS_HELP_STRING([--with-clfft-prefix=DIR], [Directory where clFFT is installed.])])

# Set CFLAGS_CLFFT only if not set from environment
if test x"$CFLAGS_CLFFT" = x; then
  case $with_clfft_prefix in
    "") CFLAGS_CLFFT="-I/usr/include" ;;
    *)  CFLAGS_CLFFT="-I$with_clfft_prefix/include" ;;
  esac
fi

AC_ARG_WITH(clfft-include, [AS_HELP_STRING([--with-clfft-include=DIR], [Directory where clfft headers are installed.])])
case $with_clfft_include in
  "") ;;
  *)  CFLAGS_CLFFT="-I$with_clfft_include" ;;
esac

CFLAGS="$CFLAGS_CLFFT $acx_clfft_save_CFLAGS"

AC_MSG_CHECKING([for clfft])

CFLAGS="$CFLAGS_CLFFT $acx_clfft_save_CFLAGS"

if test ! -z "$LIBS_CLFFT"; then
  LIBS="$LIBS_CLFFT $acx_clfft_save_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <clFFT.h>]],[[
  cl_uint cl_major, cl_minor, cl_patch;
  clfftGetVersion(&cl_major, &cl_minor, &cl_patch);]])], [acx_clfft_ok=yes], [])
fi

if test ! -z "$with_clfft_prefix"; then
  if test x"$acx_clfft_ok" = xno; then
    LIBS_CLFFT="-L$with_clfft_prefix/lib64/ -lclFFT"
    LIBS="$LIBS_CLFFT $acx_clfft_save_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <clFFT.h>]],[[
    cl_uint cl_major, cl_minor, cl_patch;
    clfftGetVersion(&cl_major, &cl_minor, &cl_patch);]])], [acx_clfft_ok=yes], [])
  fi

  if test x"$acx_clfft_ok" = xno; then
    LIBS_CLFFT="-L$with_clfft_prefix/lib/ -lclFFT"
    LIBS="$LIBS_CLFFT $acx_clfft_save_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <clFFT.h>]],[[
    cl_uint cl_major, cl_minor, cl_patch;
    clfftGetVersion(&cl_major, &cl_minor, &cl_patch);]])], [acx_clfft_ok=yes], [])    
  fi

fi
  
AC_MSG_RESULT([$acx_clfft_ok ($CFLAGS_CLFFT $LIBS_CLFFT)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_clfft_ok" = xyes; then
  AC_DEFINE(HAVE_CLAMDFFT, 1, [Defined if you have the CLFFT library.])
else
  AC_MSG_ERROR([Could not find the required clfft library.])
fi

AC_SUBST([CFLAGS_CLFFT])
AC_SUBST([LIBS_CLFFT])
CFLAGS="$acx_clfft_save_CFLAGS"
LIBS="$acx_clfft_save_LIBS"
])dnl ACX_CLFFT


AC_DEFUN([ACX_CLBLAS], [
acx_clblas_ok=no

dnl Backup LIBS and CFLAGS
acx_clblas_save_LIBS="$LIBS"
acx_clblas_save_CFLAGS="$CFLAGS"

dnl Check if the library was given in the command line
AC_ARG_WITH(clblas-prefix, [AS_HELP_STRING([--with-clblas-prefix=DIR], [Directory where clBLAS is installed.])])

# Set CFLAGS_CLBLAS only if not set from environment
if test x"$CFLAGS_CLBLAS" = x; then
  case $with_clblas_prefix in
    "") CFLAGS_CLBLAS="-I/usr/include" ;;
    *)  CFLAGS_CLBLAS="-I$with_clblas_prefix/include" ;;
  esac
fi

AC_ARG_WITH(clblas-include, [AS_HELP_STRING([--with-clblas-include=DIR], [Directory where clblas headers are installed.])])
case $with_clblas_include in
  "") ;;
  *)  CFLAGS_CLBLAS="-I$with_clblas_include" ;;
esac

CFLAGS="$CFLAGS_CLBLAS $acx_clblas_save_CFLAGS"

AC_MSG_CHECKING([for clblas])

CFLAGS="$CFLAGS_CLBLAS $acx_clblas_save_CFLAGS"

if test ! -z "$LIBS_CLBLAS"; then
  LIBS="$LIBS_CLBLAS $acx_clblas_save_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <clBLAS.h>]],[[
  cl_uint cl_major, cl_minor, cl_patch;
  clblasGetVersion(&cl_major, &cl_minor, &cl_patch);]])], [acx_clblas_ok=yes], [])
fi

if test ! -z "$with_clblas_prefix"; then
  if test x"$acx_clblas_ok" = xno; then
    LIBS_CLBLAS="-L$with_clblas_prefix/lib64/ -lclBLAS"
    LIBS="$LIBS_CLBLAS $acx_clblas_save_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <clBLAS.h>]],[[
    cl_uint cl_major, cl_minor, cl_patch;
    clblasGetVersion(&cl_major, &cl_minor, &cl_patch);]])], [acx_clblas_ok=yes], [])
  fi

  if test x"$acx_clblas_ok" = xno; then
    LIBS_CLBLAS="-L$with_clblas_prefix/lib/ -lclBLAS"
    LIBS="$LIBS_CLBLAS $acx_clblas_save_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <clBLAS.h>]],[[
    cl_uint cl_major, cl_minor, cl_patch;
    clblasGetVersion(&cl_major, &cl_minor, &cl_patch);]])], [acx_clblas_ok=yes], [])    
  fi

fi
  
AC_MSG_RESULT([$acx_clblas_ok ($CFLAGS_CLBLAS $LIBS_CLBLAS)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_clblas_ok" = xyes; then
  AC_DEFINE(HAVE_CLAMDBLAS, 1, [Defined if you have the CLBLAS library.])
else
  AC_MSG_ERROR([Could not find the required clblas library.])
fi

AC_SUBST([CFLAGS_CLBLAS])
AC_SUBST([LIBS_CLBLAS])
CFLAGS="$acx_clblas_save_CFLAGS"
LIBS="$acx_clblas_save_LIBS"
])dnl ACX_CLBLAS
