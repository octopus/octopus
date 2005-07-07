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
!!
!! $Id$

#include "global.h"

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains a data type (loct_spline_type) to contain 1D functions,
! along with a series of procedures to manage them (to define them, to obtain
! its values, to operate on them, etc). The internal representation of the
! functions is done through cubic splines, handled by the GSL library. For
! the user of the module, this internal representation is hidden, one just
! works with what are called hereafter "spline functions".
!
! WARNINGS:
!  * The problem of precision has been disregarded, so it probably does not
!    work if the code is compiled with single precision. To be fixed.
!
! To define a function, one must supply a set {x(i),y(i)} of pairs of values
! -- the abscissa and the value of the function.
!
! [*] DATA TYPE:
!     To define a spline function:
!     type(loct_spline_type) :: f
!
! [1] INITIALIZATION:
!     Before using any function, one should initialize it:
!
!     Interface:
!     subroutine loct_spline_init(spl)
!        type(loct_spline_type), intent(out) :: spl [or spl(:) or spl(:, :)]
!     end subroutine loct_spline_init
!
!     Usage:
!     call loct_spline_init(f)
!
! [2] FINALIZATION:
!     To empty any allocated space, one should finalize the function:
!
!     Interface:
!     subroutine loct_spline_end(spl)
!        type(loct_spline_type), intent(inout) :: spl [or spl(:) or spl(:, :)]
!     end subroutine loct_spline_end
!
!     Usage
!     type(loct_spline_type) :: f
!     call loct_spline_end(f)
!
! [3] TO DEFINE A FUNCTION:
!     To "fill" an initialized function f,use loct_spline_fit
!
!     Interface:
!     subroutine loct_spline_fit(n, x, y, spl)
!        integer, intent(in) :: nrc
!        real(X), intent(in) :: x(n), y(n)
!        type(loct_spline_type), intent(out) :: spl
!     end subroutine spline_fit8
!
!     (X may be 4 or eight, for single or double precision)
!
!     Usage:
!     call loct_spline_fit(n, x, y, f)
!     n is the number of values that are supplied, x the abscissas, and y
!     the value of the function to represent at each point.
!
! [4] FUNCTION VALUES:
!     To retrieve the value of a function at a given point:
!
!     Interface:
!
!     real(8) function loct_splint(spl, x)
!       type(loct_spline_type), intent(in) :: spl
!       real(8), intent(in) :: x
!     end function loct_splint
!
!     Usage:
!     To retrieve the value of function f at point x, and place it into
!     real value val:
!     val = loct_splint(f, x)
!
! [5] SUM:
!     If you have two defined functions, f and g, and an initialized function
!     h, you may sum f and g and place the result onto g. For this purpose,
!     the grid that defined the first operand, f, will be used -- the values
!     of g will be interpolated to get the sum and place it in h.
!
!     Interface:
!     subroutine loct_spline_sum(spl1, spl2, splsum)
!       type(loct_spline_type), intent(in)  :: spl1, spl2
!       type(loct_spline_type), intent(out) :: splsum
!     end subroutine loct_spline_sum
!
!     Usage:
!     call loct_spline_init(f)
!     call loct_spline_init(g)
!     call loct_spline_init(h)
!     call loct_spline_fit(npointsf, xf, yf, f)
!     call loct_spline_fit(npointsg, xg, yg, g)
!     call loct_spline_sum(f, g, h)
!
! [6] MULTIPLICATION BY A SCALAR
!     You may multiply a given spline-represented spline by a real number:
!
!     Interface:
!     subroutine loct_spline_times(a, spl)
!       type(loct_spline_type), intent(inout)  :: spl
!       real(8), intent(in) :: a
!     end subroutine loct_spline_times
!
!     Usage:
!     call loct_spline_init(f)
!     call loct_spline_fit(npoints, x, y, f) ! Fill f with y values at x points
!     call loct_spline_times(a, f) ! Now f contains a*y values at x points.
!
! [7] INTEGRAL:
!     Given a defined function, the function loct_spline_integral returns its
!     integral. The interval of integration may or may not be supplied.
!
!     Interface:
!     real(8) function loct_spline_integral(spl [,a,b])
!       type(loct_spline_integral), intent(in) :: spl
!       real(8), intent(in), optional :: a, b
!     end function loct_spline_integral
!
! [8] DOT PRODUCT:
!     Given two defined functions f and g, the function loct_spline_dotp returns
!     the value of their dot-product: int {dx f(x)*g(x)}. The mesh used to do
!     so the mesh of the fist-function (note that as a result the definition
!     is no longer conmutative).
!
!     Interface:
!     real(8) function loct_spline_dotp(spl1, spl2)
!       type(loct_spline_type), intent(in) :: spl1, spl2
!     end function loct_spline_dotp
!
! Note: The following routines, loct_spline_3dft, loct_spline_cut and loct_spline_filter,
! assume that the spline functions are the radial part of a 3 dimensional function with
! spherical symmetry. This is why the Fourier transform of F(\vec{r}) = f(r), is:
!       F(\vec{g}) = f(g) = \frac{4\pi}{g} \int_{0}^{\infty} { r*sin(g*r)*f(r) }
! which coincides with the inverse Fourier transform, except that the inverse Fourier
! transform should be multiplied by a (2*\pi)^{-3} factor.
!
! [9] FOURIER TRANSFORM:
!     If a spline function f is representing the radial part of a spherically
!     symmetric fuction F(\vec{r}), its Fourier transform is:
!       F(\vec{g}) = f(g) = \frac{4\pi}{g} \int_{0}^{\infty} { r*sin(g*r)*f(r) }
!     It is assumed that f is defined in some interval [0,rmax], and that it is
!     null at rmax and beyond. One may obtain f(g) by using loct_spline_3dft.
!     The result is placed on the spline data type splw. This has to be initialized,
!     and may or may not be filled. If it is filled, the abscissas that it has
!     are used to define the function. If it is not filled, an equally spaced
!     grid is constructed to define the faction, in the interval [0, gmax], where
!     gmax has to be supplied by the caller.
!
!     Interface:
!     subroutine loct_spline_3dft(spl, splw, gmax)
!       type(loct_spline_type), intent(in)    :: spl
!       type(loct_spline_type), intent(inout) :: splw
!       real(8), intent(in), optional :: gmax
!     end subroutine loct_spline_3dft
!
! [10] BESSEL TRANSFORM:
!
!
! [11] CUTTING A FUNCTION:
!     loct_spline_cut multiplies a given function by a cutoff-function, which
!     is defined to be one in [0, cutoff], and \exp\{-beta*(x/cutoff-1)^2\}
!
!     Interface:
!     subroutine loct_spline_cut(spl, cutoff, beta)
!       type(loct_spline_type), intent(in) :: spl
!       real(8), intent(in) :: cutoff, beta
!     end subroutine loct_spline_cut
!
! [12] FILTERING A FUNCTION, BOTH IN REAL AND FOURIER SPACE:
!     The function loct_spline_filter permits to filter out high-values
!     of a given spline function, either in real or in Fourier space.
!     If the optional argument fs is supplied, a filter in Fourier space
!     will be done by moving to the Fourier representation and then calling 
!     loct_spline_cut with cutoff = fs(1) and beta = fs(2). If the optional 
!     argument rs is supplied, a filter in real space will be done by calling 
!     loct_spline_cut with cutoff = rs(1) and beta = rs(2). If both arguments
!     are supplied, the Fourier filter will be applied *before*.
!
!     Interface:
!     subroutine loct_spline_filter(spl, fs, rs)
!       type(loct_spline_type), intent(inout) :: spl
!       real(8), intent(in), optional :: fs(2)
!       real(8), intent(in), optional :: rs(2)
!     end subroutine loct_spline_filter
!
! [13] PRINTING A FUNCTION:
!     It prints to a file the (x,y) values that were used to define a function.
!     The file is pointed by its Fortran unit given by argument iunit.
!
!     Interface:
!     subroutine loct_spline_print(spl, iunit)
!       type(loct_spline_type), intent(in) :: spl [ or spl(:) or spl(:, :)]
!       integer, intent(in) :: iunit
!     end subroutine loct_spline_print
!
! [14] DERIVATE A FUNCTION:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/!
module lib_oct_gsl_spline
  use lib_oct
  implicit none

  ! Define the which routines can be seen from the outside
  private
  public :: loct_spline_type,     & ! [*]
            loct_spline_init,     & ! [1]
            loct_spline_end,      & ! [2]
            loct_spline_fit,      & ! [3]
            loct_splint,          & ! [4]
            loct_spline_sum,      & ! [5]
            loct_spline_times,    & ! [6]
            loct_spline_integral, & ! [7]
            loct_spline_dotp,     & ! [8]
            loct_spline_3dft,     & ! [9]
            loct_spline_besselft, & ! [10]
            loct_spline_cut,      & ! [11]
            loct_spline_filter,   & ! [12]
            loct_spline_print,    & ! [13]
            loct_spline_der

  ! the basic spline datatype
  type loct_spline_type
    integer(POINTER_SIZE) :: spl, acc
  end type loct_spline_type

  ! Both the filling of the fuction, and the retrieval of the values
  ! may be done using single or double precision values.
  interface loct_spline_fit
    module procedure spline_fit4
    module procedure spline_fit8
  end interface
  interface loct_splint
    module procedure splint4
    module procedure splint8
  end interface

  ! The integral may be done with or without integration limits, but
  ! we want the interface to be common.
  interface loct_spline_integral
    module procedure loct_spline_integral_full
    module procedure loct_spline_integral_limits
  end interface

  ! Some operations may be done for one spline-function, or for an array of them
  interface loct_spline_init
    module procedure loct_spline_init_0
    module procedure loct_spline_init_1
    module procedure loct_spline_init_2
  end interface
  interface loct_spline_end
    module procedure loct_spline_end_0
    module procedure loct_spline_end_1
    module procedure loct_spline_end_2
  end interface
  interface loct_spline_print
    module procedure loct_spline_print_0
    module procedure loct_spline_print_1
    module procedure loct_spline_print_2
  end interface

  interface loct_spline_filter
    module procedure loct_spline_filter_ft
    module procedure loct_spline_filter_bessel
  end interface

  ! These are interfaces to functions defined in oct_gsl_f.c, which  in turn
  ! take care of calling the GSL library.
  interface
    subroutine oct_spline_end(spl, acc)
      integer(POINTER_SIZE), intent(inout) :: spl, acc
    end subroutine oct_spline_end

    subroutine oct_spline_fit(nrc, x, y, spl, acc)
      integer, intent(in) :: nrc
      real(8), intent(in) :: x, y
      integer(POINTER_SIZE), intent(inout) :: spl, acc
    end subroutine oct_spline_fit
    
    real(8) function oct_spline_eval(x, spl, acc)
      real(8), intent(in) :: x
      integer(POINTER_SIZE), intent(in) :: spl, acc
    end function oct_spline_eval

    real(8) function oct_spline_eval_der(x, spl, acc)
      real(8), intent(in) :: x
      integer(POINTER_SIZE), intent(in) :: spl, acc
    end function oct_spline_eval_der

    integer function oct_spline_npoints(spl)
      integer(POINTER_SIZE), intent(in) :: spl
    end function oct_spline_npoints

    subroutine oct_spline_x(spl, x)
      integer(POINTER_SIZE), intent(in) :: spl
      real(8), intent(out) :: x
    end subroutine oct_spline_x

    subroutine oct_spline_y(spl, y)
      integer(POINTER_SIZE), intent(in) :: spl
      real(8), intent(out) :: y
    end subroutine oct_spline_y

    real(8) function oct_spline_eval_integ(spl, a, b, acc)
      integer(POINTER_SIZE) :: spl, acc
      real(8) :: a, b
    end function oct_spline_eval_integ
  end interface

  real(8), parameter :: M_PI = CNST(3.141592653589793)

contains

  subroutine loct_spline_init_0(spl)
    type(loct_spline_type), intent(out) :: spl
    spl%spl = 0; spl%acc = 0
  end subroutine loct_spline_init_0

  subroutine loct_spline_init_1(spl)
    type(loct_spline_type), intent(out) :: spl(:)
    integer :: i
    do i = 1, size(spl)
       call loct_spline_init_0(spl(i))
    enddo
  end subroutine loct_spline_init_1

  subroutine loct_spline_init_2(spl)
    type(loct_spline_type), intent(out) :: spl(:, :)
    integer :: i, j
    do i = 1, size(spl, 1)
       do j = 1, size(spl, 2)
          call loct_spline_init_0(spl(i, j))
       enddo
    enddo
  end subroutine loct_spline_init_2

  subroutine loct_spline_end_0(spl)
    type(loct_spline_type), intent(inout) :: spl
    if(spl%spl.ne.0 .and. spl%acc.ne.0) then
      call oct_spline_end(spl%spl, spl%acc)
    end if
    spl%spl = 0; spl%acc = 0
  end subroutine loct_spline_end_0

  subroutine loct_spline_end_1(spl)
    type(loct_spline_type), intent(inout) :: spl(:)
    integer :: i
    do i = 1, size(spl)
       call loct_spline_end_0(spl(i))
    enddo
  end subroutine loct_spline_end_1

  subroutine loct_spline_end_2(spl)
    type(loct_spline_type), intent(inout) :: spl(:, :)
    integer :: i, j
    do i = 1, size(spl, 1)
       do j = 1, size(spl, 2)
          call loct_spline_end_0(spl(i, j))
       enddo
    enddo
  end subroutine loct_spline_end_2

  subroutine spline_fit8(nrc, rofi, ffit, spl)
    integer, intent(in) :: nrc
    real(8), intent(in) :: ffit(nrc), rofi(nrc)
    type(loct_spline_type), intent(out) :: spl
    call oct_spline_fit(nrc, rofi(1), ffit(1), spl%spl, spl%acc)
  end subroutine spline_fit8

  subroutine spline_fit4(nrc, rofi, ffit, spl)
    integer, intent(in) :: nrc
    real(4), intent(in) :: rofi(nrc), ffit(nrc)
    type(loct_spline_type), intent(out) :: spl
    
    real(8), allocatable :: rofi8(:), ffit8(:)

    allocate(rofi8(nrc), ffit8(nrc))
    rofi8 = real(rofi, kind=8)
    ffit8 = real(ffit, kind=8)
    call oct_spline_fit(nrc, rofi8(1), ffit8(1), spl%spl, spl%acc)
    deallocate(rofi8, ffit8)
  end subroutine spline_fit4

  real(8) function splint8(spl, x)
    type(loct_spline_type), intent(in) :: spl
    real(8), intent(in) :: x
  
    splint8 = oct_spline_eval(x, spl%spl, spl%acc)
  end function splint8

  real(4) function splint4(spl, x)
    type(loct_spline_type), intent(in) :: spl
    real(4), intent(in) :: x
  
    splint4 = real(oct_spline_eval(real(x, kind=8), spl%spl, spl%acc), kind=4)
  end function splint4

  subroutine loct_spline_sum(spl1, spl2, splsum)
    type(loct_spline_type), intent(in)  :: spl1, spl2
    type(loct_spline_type), intent(out) :: splsum
    integer :: npoints, i

    real(8), allocatable :: x(:), y(:), y2(:)

    npoints = oct_spline_npoints(spl1%spl)

    allocate(x(npoints), y(npoints), y2(npoints))

    call oct_spline_x(spl1%spl, x(1))
    call oct_spline_y(spl1%spl, y(1))

    do i = 1, npoints
       y2(i) = splint8(spl2, x(i))
    enddo

    y2 = y2 + y
    call oct_spline_fit(npoints, x(1), y2(1), splsum%spl, splsum%acc)

    deallocate(x, y, y2)
  end subroutine loct_spline_sum

  subroutine loct_spline_times(a, spl)
    real(8), intent(in) :: a
    type(loct_spline_type), intent(inout)  :: spl

    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    npoints = oct_spline_npoints(spl%spl)
    allocate(x(npoints), y(npoints))
    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))
    call oct_spline_end(spl%spl, spl%acc)  
    do i = 1, npoints
       y(i) = a*y(i)
    enddo
    call oct_spline_fit(npoints, x(1), y(1), spl%spl, spl%acc)

    deallocate(x, y)
  end subroutine loct_spline_times

  real(8) function loct_spline_integral_full(spl) result(res)
    type(loct_spline_type), intent(in) :: spl
    integer :: npoints
    real(8), allocatable :: x(:)
    npoints = oct_spline_npoints(spl%spl)
    allocate(x(npoints))
    call oct_spline_x(spl%spl, x(1))
    res = oct_spline_eval_integ(spl%spl, x(1), x(npoints), spl%acc)
    deallocate(x)
  end function loct_spline_integral_full

  real(8) function loct_spline_integral_limits(spl, a, b) result(res)
    type(loct_spline_type), intent(in) :: spl
    real(8), intent(in) :: a, b
    res = oct_spline_eval_integ(spl%spl, a, b, spl%acc)
  end function loct_spline_integral_limits

  real(8) function loct_spline_dotp(spl1, spl2) result (res)
    type(loct_spline_type), intent(in) :: spl1, spl2

    type(loct_spline_type) :: aux
    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    npoints = oct_spline_npoints(spl1%spl)    
    allocate(x(npoints), y(npoints))
    call oct_spline_x(spl1%spl, x(1))
    call oct_spline_y(spl1%spl, y(1))
    do i = 1, npoints
       y(i) = y(i)*oct_spline_eval(x(i), spl2%spl, spl2%acc)
    enddo
    call loct_spline_init(aux)
    call loct_spline_fit(npoints, x, y, aux)
    res = oct_spline_eval_integ(aux%spl, x(1), x(npoints), aux%acc)

  end function loct_spline_dotp

  subroutine loct_spline_3dft(spl, splw, gmax)
    type(loct_spline_type), intent(in)    :: spl
    type(loct_spline_type), intent(inout) :: splw
    real(8), intent(in), optional :: gmax
    
    type(loct_spline_type) :: aux
    real(8) :: g, dg
    integer :: np
    integer :: npoints, i, j
    real(8), allocatable :: x(:), y(:), y2(:), xw(:), yw(:)

    npoints = oct_spline_npoints(spl%spl)    
    allocate(x(npoints), y(npoints),y2(npoints))
    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))

    ! Check wether splw comes with a defined grid, or else build it.
    if(splw%spl.ne.0) then
      np = oct_spline_npoints(splw%spl)
      allocate(xw(np), yw(np))
      call oct_spline_x(splw%spl, xw(1))
      ! But now we need to kill the input:
      call loct_spline_end(splw)
    else
      np = 200 ! hard coded value
      dg = gmax/(np-1)
      allocate(xw(np), yw(np))     
      do i = 1, np
         g = (i-1)*dg
         xw(i) = g
      enddo
    endif

        ! The first point, xw(1) = 0.0 and it has to be treated separately.
        do j = 1, npoints
           y2(j) = CNST(4.0)*M_PI*y(j)*x(j)**2
        enddo
        call loct_spline_init(aux)
        call oct_spline_fit(npoints, x(1), y2(1), aux%spl, aux%acc)
        yw(1) = oct_spline_eval_integ(aux%spl, x(1), x(npoints), aux%acc)
        call loct_spline_end(aux)
    do i = 2, np
        do j = 1, npoints
           y2(j) = (CNST(4.0)*M_PI/xw(i))*y(j)*x(j)*sin(xw(i)*x(j))
        enddo
        call loct_spline_init(aux)
        call oct_spline_fit(npoints, x(1), y2(1), aux%spl, aux%acc)
        yw(i) = oct_spline_eval_integ(aux%spl, x(1), x(npoints), aux%acc)
        call loct_spline_end(aux)
    enddo

    call loct_spline_init(splw)
    call oct_spline_fit(np, xw(1), yw(1), splw%spl, splw%acc)

    deallocate(x, y, y2, xw, yw)
  end subroutine loct_spline_3dft

  subroutine loct_spline_besselft(spl, splw, l, gmax)
    type(loct_spline_type), intent(in)    :: spl
    type(loct_spline_type), intent(inout) :: splw
    integer, intent(in) :: l
    real(8), intent(in), optional :: gmax
    
    type(loct_spline_type) :: aux
    real(8) :: g, dg
    integer :: np
    integer :: npoints, i, j
    real(8), allocatable :: x(:), y(:), y2(:), xw(:), yw(:)

    npoints = oct_spline_npoints(spl%spl)    
    allocate(x(npoints), y(npoints),y2(npoints))
    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))

    ! Check wether splw comes with a defined grid, or else build it.
    if(splw%spl.ne.0) then
      np = oct_spline_npoints(splw%spl)
      allocate(xw(np), yw(np))
      call oct_spline_x(splw%spl, xw(1))
      ! But now we need to kill the input:
      call loct_spline_end(splw)
    else
      np = 1000 ! hard coded value
      dg = gmax/(np-1)
      allocate(xw(np), yw(np))     
      do i = 1, np
         g = (i-1)*dg
         xw(i) = g
      enddo
    endif

    do i = 1, np
        do j = 1, npoints
           y2(j) = y(j)*x(j)**2*loct_sph_bessel(l, x(j)*xw(i))
        enddo
        call loct_spline_init(aux)
        call oct_spline_fit(npoints, x(1), y2(1), aux%spl, aux%acc)
        yw(i) = sqrt(CNST(2.0)/M_PI)*oct_spline_eval_integ(aux%spl, x(1), x(npoints), aux%acc)
        call loct_spline_end(aux)
    enddo

    call loct_spline_init(splw)
    call oct_spline_fit(np, xw(1), yw(1), splw%spl, splw%acc)

    deallocate(x, y, y2, xw, yw)
  end subroutine loct_spline_besselft

  subroutine loct_spline_cut(spl, cutoff, beta)
    type(loct_spline_type), intent(inout) :: spl
    real(8), intent(in) :: cutoff, beta

    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    npoints = oct_spline_npoints(spl%spl)
    allocate(x(npoints), y(npoints))
    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))
    call oct_spline_end(spl%spl, spl%acc)
    do i = npoints, 1, -1
       if(x(i)<cutoff) exit
       y(i) = y(i) * exp(-beta*(x(i)/cutoff - CNST(1.0))**2)
    enddo
    call oct_spline_fit(npoints, x(1), y(1), spl%spl, spl%acc)

    deallocate(x, y)
  end subroutine loct_spline_cut

  subroutine loct_spline_filter_ft(spl, fs, rs)
    type(loct_spline_type), intent(inout) :: spl
    real(8), intent(in), optional :: fs(2)
    real(8), intent(in), optional :: rs(2)

    type(loct_spline_type) :: splw

    if(present(fs)) then
      call loct_spline_init(splw)
      call loct_spline_3dft(spl, splw, CNST(2.0)*fs(1))
      call loct_spline_cut(splw, fs(1), fs(2))
      call loct_spline_3dft(splw, spl)
      call loct_spline_times(CNST(1.0)/(CNST(2.0)*M_PI)**3, spl)
      call loct_spline_end(splw)
    endif
    if(present(rs)) then
      call loct_spline_cut(spl, rs(1), rs(2))  
    endif

  end subroutine loct_spline_filter_ft

  subroutine loct_spline_filter_bessel(spl, l, fs, rs)
    type(loct_spline_type), intent(inout) :: spl
    integer, intent(in) :: l
    real(8), intent(in), optional :: fs(2)
    real(8), intent(in), optional :: rs(2)

    type(loct_spline_type) :: splw

    if(present(fs)) then
      call loct_spline_init(splw)
      call loct_spline_besselft(spl, splw, l, CNST(2.0)*fs(1))
      call loct_spline_cut(splw, fs(1), fs(2))
      call loct_spline_besselft(splw, spl, l)
      call loct_spline_end(splw)
    endif
    if(present(rs)) then
      call loct_spline_cut(spl, rs(1), rs(2))  
    endif

  end subroutine loct_spline_filter_bessel

  subroutine loct_spline_der(spl, dspl)
    type(loct_spline_type), intent(in)    :: spl
    type(loct_spline_type), intent(inout) :: dspl


    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    ! Use the grid of dspl if it is present, otherwise use the same one of spl.
    if(dspl%spl == 0) then ! use the grid of spl
      npoints = oct_spline_npoints(spl%spl)
      allocate(x(npoints), y(npoints))
      call oct_spline_x(spl%spl, x(1))
    else ! use the grid of dspl
      npoints = oct_spline_npoints(dspl%spl)
      allocate(x(npoints), y(npoints))    
      call oct_spline_x(dspl%spl, x(1))
    endif
    do i = 1, npoints
       y(i) = oct_spline_eval_der(x(i), spl%spl, spl%acc)
    enddo
    call oct_spline_fit(npoints, x(1), y(1), dspl%spl, dspl%acc)

  end subroutine loct_spline_der

  subroutine loct_spline_print_0(spl, iunit)
    type(loct_spline_type), intent(in) :: spl
    integer, intent(in) :: iunit

    integer :: np, i
    real(8), allocatable :: x(:), y(:)

    np = oct_spline_npoints(spl%spl)
    allocate(x(np), y(np))
    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))
    do i = 1, np
       write(iunit, '(2f16.8)') x(i), y(i)
    enddo
    deallocate(x, y)
  end subroutine loct_spline_print_0

  subroutine loct_spline_print_1(spl, iunit)
    type(loct_spline_type), intent(in) :: spl(:)
    integer, intent(in) :: iunit

    character(len=4)  :: fm
    integer :: np, i, n, j
    real(8), allocatable :: x(:), y(:)

    n = size(spl)
    if(n<=0) return
    write(fm,'(i4)') n + 1; fm = adjustl(fm)
    np = oct_spline_npoints(spl(1)%spl)
    allocate(x(np), y(np))
    call oct_spline_x(spl(1)%spl, x(1))
    call oct_spline_y(spl(1)%spl, y(1))
    do i = 1, np
       write(iunit, '('//trim(fm)//'f16.8)') x(i), (loct_splint(spl(j), x(i)), j = 1, size(spl))
    enddo
    deallocate(x, y)

  end subroutine loct_spline_print_1

  subroutine loct_spline_print_2(spl, iunit)
    type(loct_spline_type), intent(in) :: spl(:, :)
    integer, intent(in) :: iunit

    character(len=4)  :: fm
    integer :: np, i, n1, n2, j, k
    real(8), allocatable :: x(:), y(:)

    n1 = size(spl, 1); n2 = size(spl, 2)
    if(n1*n2<=0) return
    write(fm,'(i4)') n1*n2 + 1; fm = adjustl(fm)
    np = oct_spline_npoints(spl(1, 1)%spl)
    allocate(x(np), y(np))
    call oct_spline_x(spl(1, 1)%spl, x(1))
    call oct_spline_y(spl(1, 1)%spl, y(1))
    do i = 1, np
       write(iunit, '('//trim(fm)//'f16.8)') x(i), &
            ((loct_splint(spl(j, k), x(i)), j = 1, n1), k = 1, n2)
    enddo
    deallocate(x, y)
  end subroutine loct_spline_print_2

end module lib_oct_gsl_spline
