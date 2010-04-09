# Check whether zdotc routine in BLAS works properly.
# The test program here may either give nonzero results or crash.

AC_DEFUN([ACX_ZDOTC], [

  AC_MSG_CHECKING(whether zdotc works)
  AC_LANG_ASSERT(Fortran)

  acx_blas_save_LIBS="$LIBS"
  LIBS="$LIBS $LIBS_BLAS"

  acx_zdotc_ok=yes
  rm -f conf.zdotc

  AC_TRY_RUN([
program testzdotc

implicit none

complex*16, allocatable :: f1(:), f2(:)
complex*16 :: result1, result2
complex*16, external :: zdotc
integer :: nn, ii

nn=100

allocate(f1(nn), f2(nn))

result1 = cmplx(0.0d0,0.0d0)
do ii=1, nn
  f1(ii) = cmplx(exp(real(ii)/real(nn)), 0.0d0)
  f2(ii) = cmplx(cos(real(ii)/real(nn)), 0.0d0)
  result1 = result1 + exp(real(ii)/real(nn))*cos(real(ii)/real(nn))
end do

result2 = cmplx(0.0d0,0.0d0)
result2 = zdotc(nn,f1,1,f2,1)

print *, abs(result1-result2)

open(1, file='conf.zdotc')
if(abs(result1-result2) .lt. 1d-5) then
  write(1, '(a)') 'success'
endif

end program
], [], [acx_zdotc_ok=no], [echo $ac_n "cross compiling; assumed OK... $ac_c"])

  if test "x$acx_zdotc_ok" = "xyes"; then
    if test "x`cat conf.zdotc`" != "xsuccess"; then
      acx_zdotc_ok=no
      # program didn't crash, but gave wrong answer
    fi
  fi

  rm -f conf.zdotc
  AC_MSG_RESULT([$acx_zdotc_ok])
  LIBS="$acx_blas_save_LIBS"

  if test "x$acx_zdotc_ok" = "xno"; then
    AC_MSG_WARN([Substituting simple loop for zdotc routine.])
    AC_DEFINE_UNQUOTED(ZDOTC_BAD, 1, [Define if zdotc cannot be used.])
  fi
])
