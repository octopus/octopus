dnl This macros checks for the presence of GCC style vectorial expressions.
AC_DEFUN([ACX_GCC_VECTORS],[
    AC_MSG_CHECKING(for gcc style vectorial expressions)
    AC_LANG(C)
    AC_TRY_RUN(
      [
	typedef double v2df __attribute__ ((vector_size (16)));
	main() { 
	v2df a, b, c;
	a = a + b*c; 
 	return (sizeof(v2df) == 16) ? 0 : 1;
	}
      ],
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_GCC_VECTORS, 1,
                  [Define if the compiler supports gcc vectorial expressions])
      ],
      [
        AC_MSG_RESULT(no)
      ])
])

dnl This macro checks if the Fortran 90 system can link to another malloc(3) function.
dnl In order to do so, a fake malloc in C is compiled that sets a common block memory
dnl location when called. A Fortran 90 program is linked against this C object file,
dnl allocates an array, and returns the value of the common block location as exit code.
dnl This is 0 if the fake routine has been called, 1 otherwise.
AC_DEFUN([ACX_FAKE_MALLOC], [
    AC_MSG_CHECKING(if a fake malloc can be used)
    if test "x$cross_compiling" = "xyes"; then
       acx_fake_malloc_continue="no"
    else
       acx_fake_malloc_continue="yes"
    fi

    if test "x$acx_fake_malloc_continue" = "xyes"; then
       AC_LANG(C)
       AC_LANG_CONFTEST([
           AC_LANG_SOURCE([[
               #ifdef HAVE_STDLIB_H
     	       #include <stdlib.h>
     	       #endif
     	       extern struct {
     	         int status;
     	       } FC_FUNC_(fake_malloc, FAKE_MALLOC);

     	       void *malloc(size_t size) {
     	         FC_FUNC_(fake_malloc, FAKE_MALLOC).status = 0;
     	         return (void *)0;
     	       }

     	       void FC_FUNC(f_exit, F_EXIT)(int *status) {
     	         exit(*status);
               }
           ]])
       ])
       rm -f conftest.$ac_objext
       if _AC_EVAL_STDERR($ac_compile) &&
          AC_TRY_COMMAND([test -z "$ac_[]_AC_LANG_ABBREV[]_werror_flag"[]dnl
     	                  || test ! -s conftest.err]) &&
          AC_TRY_COMMAND([test -s conftest.$ac_objext]); then
          acx_fake_malloc_continue="yes"
     	  mv conftest.$ac_objext conftestc.$ac_objext
       else
     	  acx_fake_malloc_continue="no"
     	  rm -f conftest.err conftest.$ac_objext
       fi
       rm -f conftest.$ac_ext
    fi

    if test "x$acx_fake_malloc_continue" = "xyes"; then
        AC_LANG(Fortran)
        AC_LANG_CONFTEST([
            AC_LANG_PROGRAM([], [[
            	implicit none
            	common /fake_malloc/ status
            	integer              :: status
	    	integer              :: ierr
            	integer, allocatable :: a(:)
            	status = 1
            	allocate(a(10), stat=ierr)
            	call f_exit(status)
            ]])
       ])
       rm -f conftest.$ac_objext
       if _AC_EVAL_STDERR($ac_compile) &&
          AC_TRY_COMMAND([test -z "$ac_[]_AC_LANG_ABBREV[]_werror_flag"[]dnl
    		          || test ! -s conftest.err]) &&
          AC_TRY_COMMAND([test -s conftest.$ac_objext]); then
    	  acx_fake_malloc_continue="yes"
    	  mv conftest.$ac_objext conftestf.$ac_objext
       else
    	   acx_fake_malloc_continue="no"
    	   rm -f conftest.err conftest.$ac_objext conftestc.$ac_objext
       fi
       rm -f conftest.$ac_ext
    fi

    rm -f conftest$ac_exeext
    if test "x$acx_fake_malloc_continue" = "xyes"; then
       if _AC_EVAL_STDERR([$FC -o conftest$ac_exeext $FCFLAGS $LDFLAGS $FCFLAGS_SRCEXT \
                           conftestc.$ac_objext conftestf.$ac_objext $LIBS >&AS_MESSAGE_LOG_FD]) &&
          AC_TRY_COMMAND([test -z "$ac_[]_AC_LANG_ABBREV[]_werror_flag"[]dnl
                          || test ! -s conftest.err]) &&
          AC_TRY_COMMAND([test -s conftest$ac_exeext]); then
	  acx_fake_malloc_continue="yes"
       else
          acx_fake_malloc_continue="no"
	  rm -f conftest.err conftestc.$ac_objext conftestf.$ac_objext conftest$ac_exeext
       fi
    fi

    if test "x$acx_fake_malloc_continue" = "xyes"; then
        if AC_TRY_COMMAND(./conftest$ac_exeext); then
           acx_fake_malloc_continue="yes"
        else
	   acx_fake_malloc_continue="no"
        fi
        rm -f core *.core gmon.out bb.out conftest$ac_exeext \
            conftestc.$ac_objext conftestf.$ac_objext
    fi

    if test "x$acx_fake_malloc_continue" = "xyes"; then
        AC_MSG_RESULT([yes])
        AC_DEFINE(HAVE_FAKE_MALLOC, 1,
                  [Defined if it is possible to link an alternative version of malloc])
    else
        AC_MSG_RESULT([no])
    fi
]) #ACX_FAKE_MALLOC
