dnl A test, if malloc(3) returns addresses that are multiples of 16.
AC_DEFUN([ACX_MALLOC_ALIGNMENT], [
    AC_MSG_CHECKING(if malloc aligns memory to 16 byte boundaries)
    AC_LANG(C)
    AC_RUN_IFELSE([
        AC_LANG_SOURCE([[
            #ifdef HAVE_STDLIB_H
            #include <stdlib.h>
            #endif
            #ifdef HAVE_STDIO_H
            #include <stdio.h>
            #endif

            int main() {
              int i;
              void *ptr;
              int not_aligned;

              not_aligned = 0;

              for(i = 1; i < 500; i++) {
                ptr = malloc(i*500);
                if(((size_t)ptr)%16 != 0) {
                  not_aligned = 1;
                }
              }
              return not_aligned;
           }]])],
        [AC_MSG_RESULT(yes)
         AC_DEFINE(HAVE_16_BYTES_ALIGNED_MALLOC, 1,
                   [Defined if malloc aligns memory to 16 bytes])],
        [AC_MSG_RESULT(no)]
    )
]) # ACX_MALLOC_ALIGNMENT


dnl This macro checks if the Fortran 90 system links to the malloc(3) function.
dnl In order to do so, a fake malloc in C is compiled that sets a common block memory
dnl location when called. A Fortran 90 program is linked against this C object file,
dnl allocates an array, and returns the value of the common block location as exit code.
dnl This is 0 if the fake routine has been called, 1 otherwise.
AC_DEFUN([ACX_FC_USES_MALLOC], [
    AC_MSG_CHECKING(if the Fortran compiler uses malloc)
    if test "x$cross_compiling" = "xyes"; then
       acx_fc_uses_malloc_continue="no"
    else
       acx_fc_uses_malloc_continue="yes"
    fi

    if test "x$acx_fc_uses_malloc_continue" = "xyes"; then
       AC_LANG(C)
       AC_LANG_CONFTEST([
           AC_LANG_SOURCE([[
               #ifdef HAVE_STDLIB_H
	       #define _XOPEN_SOURCE 600
               #include <stdlib.h>
               #endif

               extern struct {
                 int status;
               } FC_FUNC_(fake_malloc, FAKE_MALLOC);

               void *malloc(size_t size) {
	         static char a[1048576];
                 static size_t current = 0;

                 FC_FUNC_(fake_malloc, FAKE_MALLOC).status = 0;
		 current += size;
                 return a + current - size;
               }

               void free(void * d) {
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
          acx_fc_uses_malloc_continue="yes"
          mv conftest.$ac_objext conftestc.$ac_objext
       else
          acx_fc_uses_malloc_continue="no"
          rm -f conftest.err conftest.$ac_objext
       fi
       rm -f conftest.$ac_ext
    fi

    if test "x$acx_fc_uses_malloc_continue" = "xyes"; then
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
          acx_fc_uses_malloc_continue="yes"
          mv conftest.$ac_objext conftestf.$ac_objext
       else
           acx_fc_uses_malloc_continue="no"
           rm -f conftest.err conftest.$ac_objext conftestc.$ac_objext
       fi
       rm -f conftest.$ac_ext
    fi

    rm -f conftest$ac_exeext
    if test "x$acx_fc_uses_malloc_continue" = "xyes"; then
       if _AC_EVAL_STDERR([$FC -o conftest$ac_exeext $FCFLAGS $LDFLAGS $FCFLAGS_SRCEXT \
                           conftestc.$ac_objext conftestf.$ac_objext $LIBS >&AS_MESSAGE_LOG_FD]) &&
          AC_TRY_COMMAND([test -z "$ac_[]_AC_LANG_ABBREV[]_werror_flag"[]dnl
                          || test ! -s conftest.err]) &&
          AC_TRY_COMMAND([test -s conftest$ac_exeext]); then
          acx_fc_uses_malloc_continue="yes"
       else
          acx_fc_uses_malloc_continue="no"
          rm -f conftest.err conftestc.$ac_objext conftestf.$ac_objext conftest$ac_exeext
       fi
    fi

    if test "x$acx_fc_uses_malloc_continue" = "xyes"; then
        if AC_TRY_COMMAND(./conftest$ac_exeext); then
           acx_fc_uses_malloc_continue="yes"
        else
           acx_fc_uses_malloc_continue="no"
        fi
        rm -f core *.core gmon.out bb.out conftest$ac_exeext \
            conftestc.$ac_objext conftestf.$ac_objext
    fi

    if test "x$acx_fc_uses_malloc_continue" = "xyes"; then
        AC_MSG_RESULT([yes])
        AC_DEFINE(FC_USES_MALLOC, 1,
                  [Defined if Fortran compiler uses malloc(3) and can link an alternative version])
    else
        AC_MSG_RESULT([no])
    fi
]) #ACX_FC_USES_MALLOC

dnl This macros checks for the presence of SSE2 support in the compiler
AC_DEFUN([ACX_C_SSE2],[
    AC_MSG_CHECKING(for sse2 support)
    AC_LANG(C)
    AC_TRY_COMPILE(
      [
        #if !defined(__SSE2__) && !defined(__amd64__) && !defined(__x86_64__)
	#error NO_SSE2
	#endif
      ],
      [
      ],     
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_C_SSE2, 1,
                  [Define if the compiler supports sse2 instructions])
      ],
      [
        AC_MSG_RESULT(no)
      ])
])