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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module loct_math_m
  implicit none

  ! Define which routines can be seen from the outside.
  private
  public ::                  &
    loct_gamma,              &
    loct_incomplete_gamma,   &
    loct_hypergeometric,     &
    loct_asinh,              &
    loct_erf,                &
    loct_erfc,               &
    loct_ylm,                &
    loct_bessel,             &
    loct_bessel_In,          &
    loct_bessel_j0,          &
    loct_bessel_j1,          &
    loct_bessel_k0,          &
    loct_bessel_k1,          &
    loct_sph_bessel,         &
    loct_legendre_sphplm,    &
    loct_sine_integral,      &
    loct_ran_init,           &
    loct_ran_end,            &
    loct_ran_gaussian,       &
    loct_ran_flat,           &
    loct_fft_optimize,       &
    loct_1dminimize,         &
    loct_minimize,           &
    loct_minimize_direct,    &
    loct_numerical_derivative


  integer, public, parameter ::     &
    MINMETHOD_STEEPEST_DESCENT = 1, &
    MINMETHOD_FR_CG            = 2, &
    MINMETHOD_PR_CG            = 3, &
    MINMETHOD_BFGS             = 4, &
    MINMETHOD_BFGS2            = 5, &
    MINMETHOD_NMSIMPLEX        = 6

  ! ---------------------------------------------------------
  !> Numerical derivative.
  interface loct_numerical_derivative
    subroutine oct_numerical_derivative(x, h, result, abserr, f)
      real(8) :: x, h, result, abserr
      interface
        subroutine f(x, fx)
          real(8), intent(in) :: x
          real(8), intent(inout) :: fx
        end subroutine f
      end interface
    end subroutine oct_numerical_derivative
  end interface loct_numerical_derivative


  ! ---------------------------------------------------------
  !> Special functions
  interface loct_gamma
    function oct_gamma(x)
      real(8) :: oct_gamma
      real(8), intent(in) :: x
    end function oct_gamma
    module procedure oct_gamma4
  end interface loct_gamma

  interface loct_incomplete_gamma
    function oct_incomplete_gamma(a, x)
      real(8) :: oct_incomplete_gamma
      real(8), intent(in) :: a, x
    end function oct_incomplete_gamma
    module procedure oct_incomplete_gamma4
  end interface loct_incomplete_gamma

  interface loct_hypergeometric
    function oct_hypergeometric(a, b, x)
      real(8) :: oct_hypergeometric
      real(8), intent(in) :: a, b, x
    end function oct_hypergeometric
    module procedure oct_hypergeometric4
  end interface loct_hypergeometric

  interface loct_bessel
    function oct_bessel(n, x)
      real(8) :: oct_bessel
      integer, intent(in) :: n
      real(8), intent(in)  :: x
    end function oct_bessel
    module procedure oct_bessel4
  end interface loct_bessel

  interface loct_bessel_in
    function oct_bessel_in(n, x)
      real(8) :: oct_bessel_in
      integer, intent(in) :: n
      real(8), intent(in)  :: x
    end function oct_bessel_in
    module procedure oct_bessel_in4
  end interface loct_bessel_in

  interface loct_sph_bessel
    function oct_sph_bessel(l, x)
      real(8) :: oct_sph_bessel
      integer, intent(in) :: l
      real(8), intent(in)  :: x
    end function oct_sph_bessel
    module procedure oct_sph_bessel4
  end interface loct_sph_bessel

  interface loct_bessel_j0
    function oct_bessel_j0(x)
      real(8) :: oct_bessel_j0
      real(8), intent(in)  :: x
    end function oct_bessel_j0
    module procedure oct_bessel_j04
  end interface loct_bessel_j0

  interface loct_bessel_j1
    function oct_bessel_j1(x)
      real(8) :: oct_bessel_j1
      real(8), intent(in)  :: x
    end function oct_bessel_j1
    module procedure oct_bessel_j14
  end interface loct_bessel_j1

  interface loct_bessel_k0
    function oct_bessel_k0(x)
      real(8) :: oct_bessel_k0
      real(8), intent(in)  :: x
    end function oct_bessel_k0
    module procedure oct_bessel_k04
  end interface loct_bessel_k0

  interface loct_bessel_k1
    function oct_bessel_k1(x)
      real(8) :: oct_bessel_k1
      real(8), intent(in)  :: x
    end function oct_bessel_k1
    module procedure oct_bessel_k14
  end interface loct_bessel_k1

  interface loct_asinh
    function oct_asinh(x)
      real(8) :: oct_asinh
      real(8), intent(in) :: x
    end function oct_asinh
    module procedure oct_asinh4
  end interface loct_asinh

  interface loct_erf
    function oct_erf(x)
      real(8) :: oct_erf
      real(8), intent(in) :: x
    end function oct_erf
    module procedure oct_erf4
  end interface loct_erf

  interface loct_erfc
    function oct_erfc(x)
      real(8) oct_erfc
      real(8), intent(in) :: x
    end function oct_erfc
    module procedure oct_erfc4
  end interface loct_erfc

  interface loct_legendre_sphplm
    function oct_legendre_sphplm(l, m, x)
      real(8) :: oct_legendre_sphplm
      integer, intent(in) :: l, m
      real(8), intent(in) :: x
    end function oct_legendre_sphplm
    module procedure oct_legendre_sphplm4
  end interface loct_legendre_sphplm

  interface loct_sine_integral
    function oct_sine_integral(x)
      real(8) :: oct_sine_integral
      real(8) :: x
    end function oct_sine_integral
    module procedure oct_sine_integral4
  end interface loct_sine_integral

  interface loct_ylm
    subroutine oct_ylm(n, x, y, z, l, m, ylm)
      integer, intent(in)  :: n
      real(8), intent(in)  :: x, y, z
      integer, intent(in)  :: l, m
      real(8), intent(out) :: ylm
    end subroutine oct_ylm
    module procedure oct_ylm4
  end interface loct_ylm

  ! ---------------------------------------------------------
  !> Functions to generate random numbers
  interface loct_ran_init
    subroutine oct_ran_init(r)
      use c_pointer_m
      type(c_ptr), intent(out) :: r
    end subroutine oct_ran_init
  end interface loct_ran_init

  interface loct_ran_end
    subroutine oct_ran_end(r)
      use c_pointer_m
      type(c_ptr), intent(inout) :: r
    end subroutine oct_ran_end
  end interface loct_ran_end

  interface loct_ran_gaussian
    function oct_ran_gaussian(r, sigma)
      use c_pointer_m
      real(8) :: oct_ran_gaussian
      type(c_ptr), intent(in) :: r
      real(8),   intent(in) :: sigma
    end function oct_ran_gaussian
    module procedure oct_ran_gaussian4
  end interface loct_ran_gaussian

  interface loct_ran_flat
    function oct_ran_flat(r, a, b)
      use c_pointer_m
      real(8) :: oct_ran_flat
      type(c_ptr), intent(in) :: r
      real(8),     intent(in) :: a
      real(8),     intent(in) :: b
    end function oct_ran_flat
    module procedure oct_ran_flat4
  end interface loct_ran_flat

  interface loct_1dminimize
    subroutine oct_1dminimize(a, b, m, f, status)
      real(8), intent(in) :: a, b, m
      interface
        subroutine f(x, fx)
          real(8), intent(in) :: x
          real(8), intent(out) :: fx
        end subroutine f
      end interface
      integer, intent(inout) :: status
    end subroutine oct_1dminimize
  end interface loct_1dminimize

  interface loct_minimize
    function oct_minimize(method, dim, x, step, tolgrad, toldr, maxiter, f, write_iter_info, minimum)
      integer :: oct_minimize
      integer, intent(in)    :: method
      integer, intent(in)    :: dim
      real(8), intent(inout) :: x
      real(8), intent(in)    :: step
      integer, intent(in)    :: maxiter
      real(8), intent(in)    :: tolgrad
      real(8), intent(in)    :: toldr
      real(8), intent(out)   :: minimum
      interface
        subroutine f(n, x, val, getgrad, grad)
          integer, intent(in) :: n
          real(8), intent(in) :: x(n)
          real(8), intent(inout) :: val
          integer, intent(in)  :: getgrad
          real(8), intent(inout) :: grad(n)
        end subroutine f
        subroutine write_iter_info(iter, n, val, maxdr, maxgrad, x)
          integer, intent(in) :: iter
          integer, intent(in) :: n
          real(8), intent(in) :: val
          real(8), intent(in) :: maxdr
          real(8), intent(in) :: maxgrad
          real(8), intent(in) :: x(n)
        end subroutine write_iter_info
      end interface
    end function oct_minimize
  end interface loct_minimize

  interface loct_minimize_direct
    function oct_minimize_direct(method, dim, x, step, toldr, maxiter, f, write_iter_info, minimum)
      integer :: oct_minimize_direct
      integer, intent(in)    :: method
      integer, intent(in)    :: dim
      real(8), intent(inout) :: x
      real(8), intent(in)    :: step
      integer, intent(in)    :: maxiter
      real(8), intent(in)    :: toldr
      real(8), intent(out)   :: minimum
      interface
        subroutine f(n, x, val)
          integer :: n
          real(8) :: x(n)
          real(8) :: val
        end subroutine f
        subroutine write_iter_info(iter, n, val, maxdr, x)
          integer, intent(in) :: iter
          integer, intent(in) :: n
          real(8), intent(in) :: val
          real(8), intent(in) :: maxdr
          real(8), intent(in) :: x(n)
        end subroutine write_iter_info
      end interface
    end function oct_minimize_direct
  end interface loct_minimize_direct

  interface loct_fft_optimize
    subroutine oct_fft_optimize(n, par)
      integer, intent(inout) :: n
      integer, intent(in) :: par
    end subroutine oct_fft_optimize
  end interface loct_fft_optimize

contains

  ! single-precision version of the functions
  real(4) function oct_gamma4(x)
    real(4), intent(in) :: x

    oct_gamma4 = real(oct_gamma(real(x, kind=8)), kind=4)
  end function oct_gamma4

  real(4) function oct_incomplete_gamma4(a, x)
    real(4), intent(in) :: a, x

    oct_incomplete_gamma4 = real(oct_incomplete_gamma(real(a, kind = 8), real(x, kind=8)), kind=4)
  end function oct_incomplete_gamma4

  real(4) function oct_hypergeometric4(a, b, x)
    real(4), intent(in) :: a, b, x

    oct_hypergeometric4 = real(oct_hypergeometric(real(a, kind = 8), real(b, kind = 8), real(x, kind = 8)), &
      kind = 4)
  end function oct_hypergeometric4

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

  real(4) function oct_legendre_sphplm4(l, m, x)
    integer, intent(in) :: l, m
    real(4), intent(in) :: x

    oct_legendre_sphplm4 = real(oct_legendre_sphplm(l, m, real(x, kind=8)), kind=4)
  end function oct_legendre_sphplm4

  real(4) function oct_sine_integral4(x)
    real(4), intent(in)  :: x

    oct_sine_integral4 = real(oct_sine_integral(real(x, kind=8)), kind=4)
  end function oct_sine_integral4

  subroutine oct_ylm4(n, x, y, z, l, m, ylm)
    integer, intent(in)  :: n
    real(4), intent(in)  :: x, y, z
    integer, intent(in)  :: l, m
    real(4), intent(out) :: ylm

    real(8) :: ylm8
    call oct_ylm(n, real(x, kind=8), real(y, kind=8), real(z, kind=8), l, m, ylm8)
    ylm = real(ylm, kind=4)
  end subroutine oct_ylm4

  real(4) function oct_ran_gaussian4(r, sigma)
    use c_pointer_m
    type(c_ptr), intent(in) :: r
    real(4),   intent(in) :: sigma

    oct_ran_gaussian4 = real(oct_ran_gaussian(r, real(sigma, kind=8)), kind=4)
  end function oct_ran_gaussian4

  real(4) function oct_ran_flat4(r, a, b)
    use c_pointer_m
    type(c_ptr), intent(in) :: r
    real(4),     intent(in) :: a
    real(4),     intent(in) :: b

    oct_ran_flat4 = real(oct_ran_flat(r, real(a, kind=8), real(b, kind=8)), kind=4)
  end function oct_ran_flat4

end module loct_math_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
