dnl This macros checks for the presence of GCC style vectorial expressions

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

AC_DEFUN([ACX_MALLOC_ALIGNMENT],[
    AC_MSG_CHECKING(if malloc memory alignes memory to 16 bytes)
    AC_LANG(C)
    AC_TRY_RUN(
      [
        #include <stdlib.h>
	#include <stdio.h>

	int main(){
	  int i;
	  void * ptr;
	  int not_aligned = 0;

	  for(i = 0; i < 100; i++){
	    ptr = malloc( i*500);
	    if ( ((size_t) ptr)%16 != 0 ) not_aligned = 1;
	  }

	  return not_aligned;
	}
      ],
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_16_BYTES_ALIGNED_MALLOC, 1,
                  [Define if malloc alignes memory to 16 bytes])
      ],
      [
        AC_MSG_RESULT(no)
      ])
])






AC_DEFUN([ACX_FAKE_MALLOC],[
    AC_MSG_CHECKING(if a fake malloc can be used)
    AC_LANG(C)
    AC_TRY_LINK(
      [
	void * malloc(int s){ return 0; }
      ],
      [
          void * p;
          p = malloc(1);
      ],
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_FAKE_MALLOC, 1,
                  [Define if it is possible to link an alternative version of malloc])
      ],
      [
        AC_MSG_RESULT(no)
      ])
])

