!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

#include "global.h"

module lib_oct
  implicit none

  ! Define the which routines can be seen from the outside
  private
  public :: loct_gamma, loct_incomplete_gamma, loct_bessel, loct_bessel_In
  public :: loct_sph_bessel, loct_asinh, loct_erf, loct_erfc, loct_ylm
  public :: loct_bessel_j0, loct_bessel_j1, loct_bessel_k0, loct_bessel_k1
  public :: loct_ran_init, loct_ran_end, loct_ran_gaussian
  public :: loct_clock, loct_getmem, loct_sysname, loct_getcwd
  public :: loct_mkdir, loct_rm, loct_number_of_lines
  public :: loct_fft_optimize, loct_wfs_list, loct_progress_bar, loct_printRecipe
  public :: write_iter_init, write_iter_clear, write_iter_flush, write_iter_end
  public :: write_iter_start, write_iter_string, write_iter_header_start, write_iter_header
  public :: write_iter_nl, write_iter_double, write_iter_int

  ! ------------------------------------------------------------
  ! write_iter functions
  interface
     subroutine write_iter_init(out,  iter, factor, file)
       integer(POINTER_SIZE) :: out
       integer               :: iter
       FLOAT                 :: factor
       character(len=*)      :: file
     end subroutine write_iter_init
     subroutine write_iter_clear(out)
       integer(POINTER_SIZE) :: out
     end subroutine write_iter_clear
     subroutine write_iter_flush(out)
       integer(POINTER_SIZE) :: out
     end subroutine write_iter_flush
     subroutine write_iter_end(out)
       integer(POINTER_SIZE) :: out
     end subroutine write_iter_end
     subroutine write_iter_start(out)
       integer(POINTER_SIZE) :: out
     end subroutine write_iter_start
     subroutine write_iter_string(out, string)
       integer(POINTER_SIZE) :: out
       character(len=*)      :: string
     end subroutine write_iter_string
     subroutine write_iter_header_start(out)
       integer(POINTER_SIZE) :: out
     end subroutine write_iter_header_start
     subroutine write_iter_header(out, string)
       integer(POINTER_SIZE) :: out
       character(len=*)      :: string
     end subroutine write_iter_header
     subroutine write_iter_nl(out)
       integer(POINTER_SIZE) :: out
     end subroutine write_iter_nl
  end interface

  interface write_iter_double
     subroutine write_iter_double_1(out, d, n)
       integer(POINTER_SIZE) :: out
       integer               :: n
       FLOAT                 :: d
     end subroutine write_iter_double_1
     subroutine write_iter_double_n(out, d, n)
       integer(POINTER_SIZE) :: out
       integer               :: n
       FLOAT                 :: d(n)
     end subroutine write_iter_double_n
  end interface
  interface write_iter_int
     subroutine write_iter_int_1(out, i, n)
       integer(POINTER_SIZE) :: out
       integer               :: n
       integer               :: i
     end subroutine write_iter_int_1
     subroutine write_iter_int_n(out, i, n)
       integer(POINTER_SIZE) :: out
       integer               :: n
       integer               :: i(n)
     end subroutine write_iter_int_n
  end interface

  ! ------------------------------------------------------------
  ! Special functions
  interface loct_gamma
    function oct_gamma(x)
      real(8) :: oct_gamma
      real(8), intent(in) :: x
    end function oct_gamma
    module procedure oct_gamma4
  end interface

  interface loct_incomplete_gamma
    function oct_incomplete_gamma(a, x)
      real(8) :: oct_incomplete_gamma
      real(8), intent(in) :: a, x
    end function oct_incomplete_gamma
    module procedure oct_incomplete_gamma4
  end interface

  interface loct_bessel
    function oct_bessel(n, x)
      real(8) :: oct_bessel
      integer, intent(in) :: n
      real(8), intent(in)  :: x
    end function oct_bessel
    module procedure oct_bessel4
  end interface

  interface loct_bessel_in
    function oct_bessel_in(n, x)
      real(8) :: oct_bessel_in
      integer, intent(in) :: n
      real(8), intent(in)  :: x
    end function oct_bessel_in
    module procedure oct_bessel_in4
  end interface

  interface loct_sph_bessel
    function oct_sph_bessel(l, x)
      real(8) :: oct_sph_bessel
      integer, intent(in) :: l
      real(8), intent(in)  :: x
    end function oct_sph_bessel
    module procedure oct_sph_bessel4
  end interface

  interface loct_bessel_j0
    function oct_bessel_j0(x)
      real(8) :: oct_bessel_j0
      real(8), intent(in)  :: x
    end function oct_bessel_j0
    module procedure oct_bessel_j04
  end interface
  
  interface loct_bessel_j1
    function oct_bessel_j1(x)
      real(8) :: oct_bessel_j1
      real(8), intent(in)  :: x
    end function oct_bessel_j1
    module procedure oct_bessel_j14
  end interface
  
  interface loct_bessel_k0
    function oct_bessel_k0(x)
      real(8) :: oct_bessel_k0
      real(8), intent(in)  :: x
    end function oct_bessel_k0
    module procedure oct_bessel_k04
  end interface
  
 interface loct_bessel_k1
    function oct_bessel_k1(x)
      real(8) :: oct_bessel_k1
      real(8), intent(in)  :: x
    end function oct_bessel_k1
    module procedure oct_bessel_k14
  end interface

  interface loct_asinh
    function oct_asinh(x)
      real(8) :: oct_asinh
      real(8), intent(in) :: x
    end function oct_asinh
    module procedure oct_asinh4
  end interface

  interface loct_erf
    function oct_erf(x)
      real(8) :: oct_erf
      real(8), intent(in) :: x
    end function oct_erf
    module procedure oct_erf4
  end interface

  interface loct_erfc
    function oct_erfc(x)
      real(8) oct_erfc
      real(8), intent(in) :: x
    end function oct_erfc
    module procedure oct_erfc4
  end interface

  interface loct_ylm
    function oct_ylm(x, y, z, l, m)
      real(8) :: oct_ylm
      real(8), intent(in) :: x, y, z
      integer, intent(in) :: l, m
    end function oct_ylm
    module procedure oct_ylm4
  end interface

  ! ------------------------------------------------------------
  ! Functions to generate random numbers
  interface loct_ran_init
    subroutine oct_ran_init(r)
      integer(POINTER_SIZE), intent(out) :: r
    end subroutine oct_ran_init
  end interface

  interface loct_ran_end
    subroutine oct_ran_end(r)
      integer(POINTER_SIZE), intent(out) :: r
    end subroutine oct_ran_end
  end interface

  interface loct_ran_gaussian
    function oct_ran_gaussian(r, sigma)
      real(8) :: oct_ran_gaussian
      integer(POINTER_SIZE), intent(in) :: r
      real(8), intent(in) :: sigma
    end function oct_ran_gaussian
    module procedure oct_ran_gaussian4
  end interface

  ! ------------------------------------------------------------
  ! System information (time, memory, sysname)
  interface loct_clock
    function oct_clock()
      real(8) :: oct_clock
    end function oct_clock
  end interface

  interface loct_getmem
    integer function oct_getmem()
    end function oct_getmem
  end interface

  interface loct_sysname
    subroutine oct_sysname(name)
      character(len=*), intent(out) :: name
    end subroutine oct_sysname
  end interface

  interface loct_getcwd
    subroutine oct_getcwd(name)
      character(len=*), intent(out) :: name
    end subroutine oct_getcwd
  end interface

  ! ------------------------------------------------------------
  ! File handling
  interface loct_mkdir
    subroutine oct_mkdir(name)
      character(len=*), intent(in) :: name
    end subroutine oct_mkdir
  end interface

  interface loct_rm
    subroutine oct_rm(name)
      character(len=*), intent(in) :: name
    end subroutine oct_rm
  end interface

  interface loct_number_of_lines
    integer function number_of_lines(filename)
      character(len=*), intent(in) :: filename
    end function number_of_lines
  end interface

  ! ------------------------------------------------------------
  ! Varia
  interface loct_fft_optimize
    subroutine oct_fft_optimize(n, p, par)
      integer, intent(inout) :: n
      integer, intent(in) :: p, par
    end subroutine oct_fft_optimize
  end interface

  interface loct_wfs_list
    subroutine oct_wfs_list(str, l)
      character(len=*), intent(in) :: str
      integer, intent(out) :: l(32)
    end subroutine oct_wfs_list
  end interface

  interface loct_progress_bar
    subroutine oct_progress_bar(a, max)
      integer, intent(in) :: a, max
    end subroutine oct_progress_bar
  end interface

  interface loct_printRecipe
    subroutine oct_printRecipe(dir, filename)
      character(len=*), intent(in)  :: dir
      character(len=*), intent(out) :: filename
    end subroutine oct_printRecipe
  end interface

contains
  ! single precision version of the functions
  real(4) function oct_gamma4(x)
    real(4), intent(in) :: x

    oct_gamma4 = real(oct_gamma(real(x, kind=8)), kind=4)
  end function oct_gamma4

  real(4) function oct_incomplete_gamma4(a, x)
    real(4), intent(in) :: a, x

    oct_incomplete_gamma4 = real(oct_incomplete_gamma(real(a, kind = 8), real(x, kind=8)), kind=4)
  end function oct_incomplete_gamma4
  
  real(4) function oct_bessel4(n, x)
    integer, intent(in) :: n
    real(4), intent(in)  :: x
    
    oct_bessel4 = real(oct_bessel(n, real(x, kind=8)), kind=4)
  end function oct_bessel4

  real(4) function oct_bessel_in4(n, x)
    integer, intent(in) :: n
    real(4), intent(in)  :: x
    
    oct_bessel_in4 = real(oct_bessel_in(n, real(x, kind=8)), kind=4)
  end function oct_bessel_In4

  real(4) function oct_sph_bessel4(l, x)
    integer, intent(in) :: l
    real(4), intent(in)  :: x
    
    oct_sph_bessel4 = real(oct_sph_bessel(l, real(x, kind=8)), kind=4)
  end function oct_sph_bessel4

  real(4) function oct_bessel_j04(x)
    real(4), intent(in)  :: x
    
    oct_bessel_j04 = real(oct_bessel_j0(real(x, kind=8)), kind=4)
  end function oct_bessel_j04

  real(4) function oct_bessel_j14(x)
    real(4), intent(in)  :: x
    
    oct_bessel_j14 = real(oct_bessel_j1(real(x, kind=8)), kind=4)
  end function oct_bessel_j14

  real(4) function oct_bessel_k04(x)
    real(4), intent(in)  :: x
    
    oct_bessel_k04 = real(oct_bessel_k0(real(x, kind=8)), kind=4)
  end function oct_bessel_k04

  real(4) function oct_bessel_k14(x)
    real(4), intent(in)  :: x
    
    oct_bessel_k14 = real(oct_bessel_k1(real(x, kind=8)), kind=4)
  end function oct_bessel_k14

  real(4) function oct_asinh4(x)
    real(4), intent(in) :: x

    oct_asinh4 = real(oct_asinh(real(x, kind=8)), kind=4)
  end function oct_asinh4
  
  real(4) function oct_erf4(x)
    real(4), intent(in)  :: x

    oct_erf4 = real(oct_erf(real(x, kind=8)), kind=4)
  end function oct_erf4

  real(4) function oct_erfc4(x)
    real(4), intent(in)  :: x

    oct_erfc4 = real(oct_erfc(real(x, kind=8)), kind=4)
  end function oct_erfc4

  real(4) function oct_ylm4(x, y, z, l, m)
    real(4), intent(in) :: x, y, z
    integer, intent(in) :: l, m
    
    oct_ylm4 = real(oct_ylm(real(x, kind=8), real(y, kind=8), real(z, kind=8), l, m), kind=4)
  end function oct_ylm4

  real(4) function oct_ran_gaussian4(r, sigma)
    integer(POINTER_SIZE), intent(in) :: r
    real(4), intent(in) :: sigma

    oct_ran_gaussian4 = real(oct_ran_gaussian(r, real(sigma, kind=8)), kind=4)
  end function oct_ran_gaussian4

end module lib_oct
