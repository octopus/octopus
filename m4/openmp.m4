# Taken from the autoconf source in 05/2020
# Copyright (C) 2001-2017, 2020 Free Software Foundation, Inc.

# This file is part of Autoconf.  This program is free
# software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Under Section 7 of GPL version 3, you are granted additional
# permissions described in the Autoconf Configure Script Exception,
# version 3.0, as published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# and a copy of the Autoconf Configure Script Exception along with
# this program; see the files COPYINGv3 and COPYING.EXCEPTION
# respectively.  If not, see <https://www.gnu.org/licenses/>.

# Written by David MacKenzie, with help from
# Akim Demaille, Paul Eggert,
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.


AC_DEFUN([AX_OPENMP],
[
  AC_PREREQ(2.62)
  OPENMP_[]_AC_LANG_PREFIX[]FLAGS=

  AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
    [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
    [AC_LINK_IFELSE([_AC_LANG_OPENMP],
       [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
       [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
        dnl Try these flags:
        dnl   GCC >= 4.2           -fopenmp
        dnl   Intel C              -qopenmp -openmp
        dnl   SGI C, PGI C         -mp
        dnl   Tru64 Compaq C       -omp
        dnl   IBM XL C (AIX, Linux) -qsmp=omp
        dnl   Cray CCE             -homp
        dnl   NEC SX               -Popenmp
        dnl   Lahey Fortran (Linux)  --openmp
        dnl   SunPRO C             -xopenmp
        dnl If in this loop a compiler is passed an option that it doesn't
        dnl understand or that it misinterprets, the AC_LINK_IFELSE test
        dnl will fail (since we know that it failed without the option),
        dnl therefore the loop will continue searching for an option, and
        dnl no output file called 'penmp' or 'mp' is created.
        for ac_option in -fopenmp -qopenmp -openmp -mp -omp -qsmp=omp -homp \
                         -Popenmp --openmp -xopenmp; do
          ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
          _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
          AC_LINK_IFELSE([_AC_LANG_OPENMP],
            [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
          _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
          if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
            break
          fi
        done])])
  case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
    "none needed" | unsupported)
      ;; #(
    *)
      OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp 
      AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled]) ;;
  esac

  AC_SUBST([OPENMP_]_AC_LANG_PREFIX[FLAGS])
])
