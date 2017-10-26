#
# Copyright (C) 2006 Steven G. Johnson <stevenj@alum.mit.edu>
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or (at
#your option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
#USA.
#
#As a special exception, the respective Autoconf Macro's copyright
#owner gives unlimited permission to copy, distribute and modify the
#configure scripts that are the output of Autoconf when processing the
#Macro. You need not follow the terms of the GNU General Public License
#when using or distributing such scripts, even though portions of the
#text of the Macro appear in them. The GNU General Public License (GPL)
#does govern all other use of the material that constitutes the
#Autoconf Macro.
#
#This special exception to the GPL applies to versions of the Autoconf
#Macro released by the Autoconf Macro Archive. When you make and
#distribute a modified version of the Autoconf Macro, you may extend
#this special exception to the GPL to apply to your modified version as
#well.
#
# http://autoconf-archive.cryp.to/ax_openmp.html
#

AC_DEFUN([AX_OPENMP], [
AC_PREREQ(2.59) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for OpenMP flag of _AC_LANG compiler], ax_cv_[]_AC_LANG_ABBREV[]_openmp, [save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
ax_cv_[]_AC_LANG_ABBREV[]_openmp=unknown
# Flags to try:  -fopenmp (gcc), -openmp and -qopenmp (Intel), -mp (SGI & PGI),
#                -xopenmp (Sun), -omp (Tru64), -qsmp=omp (AIX), none
ax_openmp_flags="-qsmp=omp -qopenmp -openmp -mp=numa -mp=nonuma -mp -xopenmp -omp -fopenmp none"
if test "x$OPENMP_[]_AC_LANG_PREFIX[]FLAGS" != x; then
  ax_openmp_flags="$OPENMP_[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flags"
fi
for ax_openmp_flag in $ax_openmp_flags; do
  case $ax_openmp_flag in
    none) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[] $_AC_LANG_PREFIX[]FLAGS" ;;
    *) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flag";;
  esac
  AC_LINK_IFELSE([AC_LANG_CALL([], [omp_set_num_threads])],
        [ax_cv_[]_AC_LANG_ABBREV[]_openmp=$ax_openmp_flag; break])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
])
if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" = "xunknown"; then
  m4_default([$2],:)
else
  if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" != "xnone"; then
    OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ax_cv_[]_AC_LANG_ABBREV[]_openmp
  fi
  m4_default([$1], [AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
])dnl AX_OPENMP




AC_DEFUN([AX_OPENMP_SIMD], [
AC_PREREQ(2.59) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for OpenMP SIMP flag of _AC_LANG compiler], ax_cv_[]_AC_LANG_ABBREV[]_openmp_simd, [save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
ax_cv_[]_AC_LANG_ABBREV[]_openmp_simd=unknown
# Flags to try:  -fopenmp-simd (gcc), -qopenmp-simd (Intel),  (SGI & PGI),
#                 (Sun),  (Tru64), (AIX), none
ax_openmp_simd_flags=" -qopenmp-simd -fopenmp-simd none"
if test "x$OPENMP_SIMD_[]_AC_LANG_PREFIX[]FLAGS" != x; then
  ax_openmp_simd_flags="$OPENMP_SIMD_[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_simd_flags"
fi

testprog="AC_LANG_PROGRAM([],[
  use OMP_lib
  implicit none
  integer :: a(128)
  integer :: i
  
  !$OMP PARALLEL DO SIMD SAFELEN(4)
  do i = 1,128
    a(i) = i
  end do
  ])"

for ax_openmp_simd_flag in $ax_openmp_simd_flags; do
  case $ax_openmp_simd_flag in
    none) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[] $_AC_LANG_PREFIX[]FLAGS" ;;
    *) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_simd_flag";;
  esac
  AC_LINK_IFELSE($testprog,
        [ax_cv_[]_AC_LANG_ABBREV[]_openmp_simd=$ax_openmp_simd_flag; break])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
])
if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp_simd" = "xunknown"; then
  m4_default([$2],:)
else
  if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp_simd" != "xnone"; then
    OPENMP_SIMD_[]_AC_LANG_PREFIX[]FLAGS=$ax_cv_[]_AC_LANG_ABBREV[]_openmp_simd
  fi
  m4_default([$1], [AC_DEFINE(HAVE_OPENMP_SIMD,1,[Define if OpenMP SIMD is enabled])])
fi
])dnl AX_OPENMP_SIMD
