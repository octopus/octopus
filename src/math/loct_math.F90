!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module loct_math_oct_m
  implicit none

  !> Define which routines can be seen from the outside.
  private
  public ::                    &
    loct_gamma,                &
    loct_incomplete_gamma,     &
    loct_erf,                  &
    loct_erfc,                 &
    loct_ylm,                  &
    loct_bessel,               &
    loct_bessel_in,            &
    loct_bessel_j0,            &
    loct_bessel_j1,            &
    loct_bessel_k0,            &
    loct_bessel_k1,            &
    loct_sph_bessel,           &
    loct_legendre_sphplm,      &
    loct_sf_laguerre_n,        &
    loct_ran_init,             &
    loct_ran_end,              &
    loct_ran_gaussian,         &
    loct_ran_flat,             &
    loct_fft_optimize,         &
    loct_combination_init,     &
    loct_combination_end,      &
    loct_combination_next,     &
    loct_get_combination


  ! ---------------------------------------------------------
  !> Special functions
  interface loct_gamma
    function oct_gamma(x)
      implicit none
      real(8) :: oct_gamma
      real(8), intent(in) :: x
    end function oct_gamma
  end interface loct_gamma

  interface loct_incomplete_gamma
    function oct_incomplete_gamma(a, x)
      implicit none
      real(8) :: oct_incomplete_gamma
      real(8), intent(in) :: a, x
    end function oct_incomplete_gamma
  end interface loct_incomplete_gamma

  interface loct_bessel
    function oct_bessel(n, x)
      implicit none
      real(8) :: oct_bessel
      integer, intent(in) :: n
      real(8), intent(in) :: x
    end function oct_bessel
  end interface loct_bessel

  interface loct_bessel_in
    function oct_bessel_in(n, x)
      implicit none
      real(8) :: oct_bessel_in
      integer, intent(in) :: n
      real(8), intent(in) :: x
    end function oct_bessel_in
  end interface loct_bessel_in

  interface loct_sph_bessel
    function oct_sph_bessel(l, x)
      implicit none
      real(8) :: oct_sph_bessel
      integer, intent(in) :: l
      real(8), intent(in) :: x
    end function oct_sph_bessel
  end interface loct_sph_bessel

  interface loct_bessel_j0
    function oct_bessel_j0(x)
      implicit none
      real(8) :: oct_bessel_j0
      real(8), intent(in) :: x
    end function oct_bessel_j0
  end interface loct_bessel_j0

  interface loct_bessel_j1
    function oct_bessel_j1(x)
      implicit none
      real(8) :: oct_bessel_j1
      real(8), intent(in) :: x
    end function oct_bessel_j1
  end interface loct_bessel_j1

  interface loct_bessel_k0
    function oct_bessel_k0(x)
      implicit none
      real(8) :: oct_bessel_k0
      real(8), intent(in) :: x
    end function oct_bessel_k0
  end interface loct_bessel_k0

  interface loct_bessel_k1
    function oct_bessel_k1(x)
      implicit none
      real(8) :: oct_bessel_k1
      real(8), intent(in) :: x
    end function oct_bessel_k1
  end interface loct_bessel_k1

  interface loct_erf
    function oct_erf(x)
      implicit none
      real(8) :: oct_erf
      real(8), intent(in) :: x
    end function oct_erf
  end interface loct_erf

  interface loct_erfc
    function oct_erfc(x)
      implicit none
      real(8) oct_erfc
      real(8), intent(in) :: x
    end function oct_erfc
  end interface loct_erfc

  interface loct_legendre_sphplm
    function oct_legendre_sphplm(l, m, x)
      implicit none
      real(8) :: oct_legendre_sphplm
      integer, intent(in) :: l, m
      real(8), intent(in) :: x
    end function oct_legendre_sphplm
  end interface loct_legendre_sphplm

  interface loct_sf_laguerre_n
    function oct_sf_laguerre_n(n, a, x)
      implicit none
      real(8) :: oct_sf_laguerre_n
      integer, intent(in) :: n
      real(8), intent(in) :: a
      real(8), intent(in) :: x
    end function oct_sf_laguerre_n
  end interface loct_sf_laguerre_n

  interface loct_ylm
    subroutine oct_ylm(n, x, y, z, l, m, ylm)
      implicit none
      integer, intent(in)  :: n
      real(8), intent(in)  :: x, y, z
      integer, intent(in)  :: l, m
      real(8), intent(out) :: ylm
    end subroutine oct_ylm
  end interface loct_ylm
  
  ! ---------------------------------------------------------
  !> Functions to generate combinations
  interface loct_combination_init
    subroutine oct_combination_init(c, n, k)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: c
      integer,     intent(in)  :: n, k
    end subroutine oct_combination_init
  end interface loct_combination_init

  interface loct_combination_end
    subroutine oct_combination_end(c)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: c
    end subroutine oct_combination_end
  end interface loct_combination_end

  interface loct_combination_next
    subroutine oct_combination_next(c, next)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: c
      integer,     intent(out)  :: next
    end subroutine oct_combination_next
  end interface loct_combination_next

  interface
    subroutine oct_get_combination(c, comb)
      use iso_c_binding
      type(c_ptr), intent(in)  :: c
      integer,     intent(out) :: comb
    end subroutine oct_get_combination
  end interface

  ! ---------------------------------------------------------
  !> Functions to generate random numbers
  interface loct_ran_init
    subroutine oct_ran_init(r)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: r
    end subroutine oct_ran_init
  end interface loct_ran_init

  interface loct_ran_end
    subroutine oct_ran_end(r)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: r
    end subroutine oct_ran_end
  end interface loct_ran_end

  interface loct_ran_gaussian
    function oct_ran_gaussian(r, sigma)
      use iso_c_binding
      implicit none
      real(8) :: oct_ran_gaussian
      type(c_ptr), intent(in) :: r
      real(8),     intent(in) :: sigma
    end function oct_ran_gaussian
  end interface loct_ran_gaussian

  interface loct_ran_flat
    function oct_ran_flat(r, a, b)
      use iso_c_binding
      implicit none
      real(8) :: oct_ran_flat
      type(c_ptr), intent(in) :: r
      real(8),     intent(in) :: a
      real(8),     intent(in) :: b
    end function oct_ran_flat
  end interface loct_ran_flat

  interface loct_fft_optimize
    subroutine oct_fft_optimize(n, par)
      implicit none
      integer, intent(inout) :: n
      integer, intent(in)    :: par
    end subroutine oct_fft_optimize
  end interface loct_fft_optimize

contains
  
  subroutine loct_get_combination(c, comb)
    use iso_c_binding
    type(c_ptr),      intent(in)  :: c
    integer,          intent(out) :: comb(0:) !< Assume C-style array indices (i.e. start from 0) 

    call oct_get_combination(c, comb(0))
  end subroutine loct_get_combination

end module loct_math_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
