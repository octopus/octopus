dnl Check size of a pointer
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


m4_include([fortran90.m4])