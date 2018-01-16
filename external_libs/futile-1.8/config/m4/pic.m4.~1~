# Define a macro to test if the Fortran and C compiler support -fPIC
#
# Copyright (c) 2011-2013 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.

AC_DEFUN([AX_FLAG_PIC],
[
  AC_MSG_CHECKING([for position-independant code option flag])

  AC_LANG_PUSH(Fortran)
  FCFLAGS_SVG=$FCFLAGS
  FCFLAGS="$FCFLAGS -fPIC"
  AC_COMPILE_IFELSE([AC_LANG_SOURCE([
subroutine api_func(a)
  integer, intent(out) :: a

  a = 1
end subroutine api_func])],
    [ax_fc_pic="yes"], [ax_fc_pic="no"])
  FCFLAGS=$FCFLAGS_SVG
  AC_LANG_POP(Fortran)

  AC_LANG_PUSH(C)
  CFLAGS_SVG=$CFLAGS
  CFLAGS="$CFLAGS -fPIC"
  AC_COMPILE_IFELSE([AC_LANG_SOURCE([
int api_func()
{
 return 1;
}])],
    [ax_c_pic="yes"], [ax_c_pic="no"])
  CFLAGS=$CFLAGS_SVG
  AC_LANG_POP(C)

  if test $ax_fc_pic = "yes" -a $ax_c_pic = "yes" ; then
    ax_flag_pic="-fPIC"
  else
    ax_flag_pic=
  fi

  AC_MSG_RESULT([$ax_flag_pic])
])
