dnl These macros check for the presence of builtin compiler macros

AC_DEFUN([ACX_C_BUILTIN_EXPECT],[
    AC_MSG_CHECKING(for __builtin_expect)
    AC_LANG(C)
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [
          int f(){return 1;};
      ],
      [
          if (__builtin_expect(f(), 0));
      ])],
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_BUILTIN_EXPECT, 1,
                  [Define if the compiler provides __builtin_expect])
      ],
      [
        AC_MSG_RESULT(no)
      ])
])


AC_DEFUN([ACX_C_BUILTIN_PREFETCH],[
    AC_MSG_CHECKING(for __builtin_prefetch)
    AC_LANG(C)
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [
          int * f;
      ],
      [
          __builtin_prefetch(f, 0, 3);
      ])],
      [
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_BUILTIN_PREFETCH, 1,
                  [Define if the compiler provides __builtin_prefetch])
      ],
      [
        AC_MSG_RESULT(no)
      ])
])
