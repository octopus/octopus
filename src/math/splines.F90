!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!/*----------------------------------------------------------------------------
! This module contains a data type (spline_t) to contain 1D functions,
! along with a series of procedures to manage them (to define them, to obtain
! its values, to operate on them, etc). The internal representation of the
! functions is done through cubic splines, handled by the GSL library. For
! the user of the module, this internal representation is hidden; one just
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
!     type(spline_t) :: f
!
! [1] INITIALIZATION:
!     Before using any function, one should initialize it:
!
!     Interface:
!     subroutine spline_init(spl)
!        type(spline_t), intent(out) :: spl [or spl(:) or spl(:, :)]
!     end subroutine spline_init
!
!     Usage:
!     call spline_init(f)
!
! [2] FINALIZATION:
!     To empty any allocated space, one should finalize the function:
!
!     Interface:
!     subroutine spline_end(spl)
!        type(spline_t), intent(inout) :: spl [or spl(:) or spl(:, :)]
!     end subroutine spline_end
!
!     Usage
!     type(spline_t) :: f
!     call spline_end(f)
!
! [3] TO DEFINE A FUNCTION:
!     To "fill" an initialized function f,use spline_fit
!
!     Interface:
!     subroutine spline_fit(n, x, y, spl)
!        integer, intent(in) :: nrc
!        real(X), intent(in) :: x(n), y(n)
!        type(spline_t), intent(out) :: spl
!     end subroutine spline_fit8
!
!     (X may be 4 or eight, for single or double precision)
!
!     Usage:
!     call spline_fit(n, x, y, f)
!     n is the number of values that are supplied, x the abscissas, and y
!     the value of the function to represent at each point.
!
! [4] FUNCTION VALUES:
!     To retrieve the value of a function at a given point:
!
!     Interface:
!
!     real(8) function spline_eval(spl, x)
!       type(spline_t), intent(in) :: spl
!       real(8), intent(in) :: x
!     end function spline_eval
!
!     Usage:
!     To retrieve the value of function f at point x, and place it into
!     real value val:
!     val = spline_eval(f, x)
!
! [5] SUM:
!     If you have two defined functions, f and g, and an initialized function
!     h, you may sum f and g and place the result onto g. For this purpose,
!     the grid that defined the first operand, f, will be used -- the values
!     of g will be interpolated to get the sum and place it in h.
!
!     Interface:
!     subroutine spline_sum(spl1, spl2, splsum)
!       type(spline_t), intent(in)  :: spl1, spl2
!       type(spline_t), intent(out) :: splsum
!     end subroutine spline_sum
!
!     Usage:
!     call spline_init(f)
!     call spline_init(g)
!     call spline_init(h)
!     call spline_fit(npointsf, xf, yf, f)
!     call spline_fit(npointsg, xg, yg, g)
!     call spline_sum(f, g, h)
!
! [6] MULTIPLICATION BY A SCALAR
!     You may multiply a given spline-represented spline by a real number:
!
!     Interface:
!     subroutine spline_times(a, spl)
!       type(spline_t), intent(inout)  :: spl
!       real(8), intent(in) :: a
!     end subroutine spline_times
!
!     Usage:
!     call spline_init(f)
!     call spline_fit(npoints, x, y, f) ! Fill f with y values at x points
!     call spline_times(a, f) ! Now f contains a*y values at x points.
!
! [7] INTEGRAL:
!     Given a defined function, the function spline_integral returns its
!     integral. The interval of integration may or may not be supplied.
!
!     Interface:
!     real(8) function spline_integral(spl [,a,b])
!       type(spline_integral), intent(in) :: spl
!       real(8), intent(in), optional :: a, b
!     end function spline_integral
!
! [8] DOT PRODUCT:
!     Given two defined functions f and g, the function spline_dotp returns
!     the value of their dot-product: int {dx f(x)*g(x)}. The mesh used to do
!     so the mesh of the fist-function (note that as a result the definition
!     is no longer conmutative).
!
!     Interface:
!     real(8) function spline_dotp(spl1, spl2)
!       type(spline_t), intent(in) :: spl1, spl2
!     end function spline_dotp
!
! Note: The following routines, spline_3dft, spline_cut and spline_filter,
! assume that the spline functions are the radial part of a 3 dimensional function with
! spherical symmetry. This is why the Fourier transform of F(\vec{r}) = f(r), is:
!       F(\vec{g}) = f(g) = \frac{4\pi}{g} \int_{0}^{\infty} { r*sin(g*r)*f(r) }
! which coincides with the inverse Fourier transform, except that the inverse Fourier
! transform should be multiplied by a (2*\pi)^{-3} factor.
!
! [9] FOURIER TRANSFORM:
!     If a spline function f is representing the radial part of a spherically
!     symmetric function F(\vec{r}), its Fourier transform is:
!       F(\vec{g}) = f(g) = \frac{4\pi}{g} \int_{0}^{\infty} { r*sin(g*r)*f(r) }
!     It is assumed that f is defined in some interval [0,rmax], and that it is
!     null at rmax and beyond. One may obtain f(g) by using spline_3dft.
!     The result is placed on the spline data type splw. This has to be initialized,
!     and may or may not be filled. If it is filled, the abscissas that it has
!     are used to define the function. If it is not filled, an equally spaced
!     grid is constructed to define the function, in the interval [0, gmax], where
!     gmax has to be supplied by the caller.
!
!     Interface:
!     subroutine spline_3dft(spl, splw, gmax)
!       type(spline_t), intent(in)    :: spl
!       type(spline_t), intent(inout) :: splw
!       real(8), intent(in), optional :: gmax
!     end subroutine spline_3dft
!
! [10] BESSEL TRANSFORM:
!
!
! [11] CUTTING A FUNCTION:
!     spline_cut multiplies a given function by a cutoff-function, which
!     is defined to be one in [0, cutoff], and \exp\{-beta*(x/cutoff-1)^2\}
!
!     Interface:
!     subroutine spline_cut(spl, cutoff, beta)
!       type(spline_t), intent(in) :: spl
!       real(8), intent(in) :: cutoff, beta
!     end subroutine spline_cut
!
! [13] PRINTING A FUNCTION:
!     It prints to a file the (x,y) values that were used to define a function.
!     The file is pointed to by its Fortran unit given by argument iunit.
!
!     Interface:
!     subroutine spline_print(spl, iunit)
!       type(spline_t), intent(in) :: spl [ or spl(:) or spl(:, :)]
!       integer, intent(in) :: iunit
!     end subroutine spline_print
!
! [14] DIFFERENTIATE A FUNCTION:
!
!----------------------------------------------------------------------------*/!
module splines_m
  use global_m
  use loct_math_m
  use messages_m
  use c_pointer_m
  use profiling_m

  implicit none

  ! Define which routines can be seen from the outside.
  private
  public ::               &
    spline_t,        & ! [*]
    spline_init,     & ! [1]
    spline_end,      & ! [2]
    spline_fit,      & ! [3]
    spline_eval,     & ! [4]
    spline_eval_vec, & ! [4]
    spline_sum,      & ! [5]
    spline_times,    & ! [6]
    spline_integral, & ! [7]
    spline_dotp,     & ! [8]
    spline_3dft,     & ! [9]
    spline_besselft, & ! [10]
    spline_cut,      & ! [11]
    spline_print,    & ! [13]
    spline_der,      &
    spline_der2,     &
    spline_div,      & 
    spline_mult,     & 
    spline_cutoff_radius

  ! the basic spline datatype
  type spline_t
    real(8)     :: x_limit(2)
    type(c_ptr) :: spl, acc
  end type spline_t

  ! Both the filling of the function, and the retrieval of the values
  ! may be done using single or double precision values.
  interface spline_fit
    module procedure spline_fit4
    module procedure spline_fit8
  end interface

  interface spline_eval
    module procedure spline_eval4
    module procedure spline_eval8
  end interface

  interface spline_eval_vec
    module procedure spline_eval4_array
    module procedure spline_eval8_array
  end interface

  ! The integral may be done with or without integration limits, but
  ! we want the interface to be common.
  interface spline_integral
    module procedure spline_integral_full
    module procedure spline_integral_limits
  end interface

  ! Some operations may be done for one spline-function, or for an array of them
  interface spline_init
    module procedure spline_init_0
    module procedure spline_init_1
    module procedure spline_init_2
  end interface
  interface spline_end
    module procedure spline_end_0
    module procedure spline_end_1
    module procedure spline_end_2
  end interface
  interface spline_print
    module procedure spline_print_0
    module procedure spline_print_1
    module procedure spline_print_2
  end interface

  ! These are interfaces to functions defined in oct_gsl_f.c, which in turn
  ! take care of calling the GSL library.
  interface
    subroutine oct_spline_end(spl, acc)
      use c_pointer_m

      type(c_ptr), intent(inout) :: spl, acc
    end subroutine oct_spline_end

    subroutine oct_spline_fit(nrc, x, y, spl, acc)
      use c_pointer_m

      integer,     intent(in) :: nrc
      real(8),     intent(in) :: x
      real(8),     intent(in) :: y
      type(c_ptr), intent(inout) :: spl
      type(c_ptr), intent(inout) :: acc
    end subroutine oct_spline_fit

    real(8) pure function oct_spline_eval(x, spl, acc)
      use c_pointer_m

      real(8),     intent(in) :: x
      type(c_ptr), intent(in) :: spl
      type(c_ptr), intent(in) :: acc      
    end function oct_spline_eval

    pure subroutine oct_spline_eval_array(nn, xf, spl, acc)
      use c_pointer_m

      integer,     intent(in)    :: nn
      real(8),     intent(inout) :: xf
      type(c_ptr), intent(in) :: spl
      type(c_ptr), intent(in) :: acc  
    end subroutine oct_spline_eval_array

    pure subroutine oct_spline_eval_array4(nn, xf, spl, acc)
      use c_pointer_m

      integer,     intent(in)    :: nn
      real(4),     intent(inout) :: xf
      type(c_ptr), intent(in)    :: spl
      type(c_ptr), intent(in)    :: acc  
    end subroutine oct_spline_eval_array4

    real(8) pure function oct_spline_eval_der(x, spl, acc)
      use c_pointer_m

      real(8),     intent(in) :: x
      type(c_ptr), intent(in) :: spl
      type(c_ptr), intent(in) :: acc
    end function oct_spline_eval_der

    real(8) pure function oct_spline_eval_der2(x, spl, acc)
      use c_pointer_m

      real(8),     intent(in) :: x
      type(c_ptr), intent(in) :: spl
      type(c_ptr), intent(in) :: acc
    end function oct_spline_eval_der2

    integer pure function oct_spline_npoints(spl)
      use c_pointer_m

      type(c_ptr), intent(in) :: spl
    end function oct_spline_npoints

    pure subroutine oct_spline_x(spl, x)
      use c_pointer_m

      type(c_ptr), intent(in)  :: spl
      real(8),     intent(out) :: x
    end subroutine oct_spline_x

    subroutine oct_spline_y(spl, y)
      use c_pointer_m

      type(c_ptr), intent(in) :: spl
      real(8),     intent(out) :: y
    end subroutine oct_spline_y

    real(8) pure function oct_spline_eval_integ(spl, a, b, acc)
      use c_pointer_m

      type(c_ptr), intent(in) :: spl
      real(8),     intent(in) :: a
      real(8),     intent(in) :: b
      type(c_ptr), intent(in) :: acc
    end function oct_spline_eval_integ
  end interface

contains

  !------------------------------------------------------------
  subroutine spline_init_0(spl)
    type(spline_t), intent(out) :: spl

    ! No PUSH SUB, called too often

    call set_null(spl%spl)
    call set_null(spl%acc)

  end subroutine spline_init_0


  !------------------------------------------------------------
  subroutine spline_init_1(spl)
    type(spline_t), intent(out) :: spl(:)

    integer :: i

    PUSH_SUB(spline_init_1)

    do i = 1, size(spl)
      call spline_init_0(spl(i))
    end do

    POP_SUB(spline_init_1)
  end subroutine spline_init_1


  !------------------------------------------------------------
  subroutine spline_init_2(spl)
    type(spline_t), intent(out) :: spl(:, :)

    integer :: i, j

    PUSH_SUB(spline_init_2)

    do i = 1, size(spl, 1)
      do j = 1, size(spl, 2)
        call spline_init_0(spl(i, j))
      end do
    end do

    POP_SUB(spline_init_2)
  end subroutine spline_init_2


  !------------------------------------------------------------
  subroutine spline_end_0(spl)
    type(spline_t), intent(inout) :: spl

    ! No PUSH SUB, called too often.

    if(c_associated(spl%spl) .and. c_associated(spl%acc)) then
      call oct_spline_end(spl%spl, spl%acc)
    end if
    call set_null(spl%spl)
    call set_null(spl%acc)

  end subroutine spline_end_0


  !------------------------------------------------------------
  subroutine spline_end_1(spl)
    type(spline_t), intent(inout) :: spl(:)

    integer :: i

    PUSH_SUB(spline_end_1)

    do i = 1, size(spl)
      call spline_end_0(spl(i))
    end do

    POP_SUB(spline_end_1)
  end subroutine spline_end_1


  !------------------------------------------------------------
  subroutine spline_end_2(spl)
    type(spline_t), intent(inout) :: spl(:, :)

    integer :: i, j

    PUSH_SUB(spline_end_2)

    do i = 1, size(spl, 1)
      do j = 1, size(spl, 2)
        call spline_end_0(spl(i, j))
      end do
    end do

    POP_SUB(spline_end_2)
  end subroutine spline_end_2


  !------------------------------------------------------------
  subroutine spline_copy(splout, splin)
    type(spline_t), intent(inout) :: splout
    type(spline_t), intent(in)    :: splin

    integer :: npoints
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_copy)

    npoints = oct_spline_npoints(splin%spl)

    SAFE_ALLOCATE( x(1:npoints))
    SAFE_ALLOCATE( y(1:npoints))

    call oct_spline_x(splin%spl, x(1))
    call oct_spline_y(splout%spl, y(1))

    call spline_fit(npoints, x, y, splout)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_copy)
  end subroutine spline_copy


  !------------------------------------------------------------
  subroutine spline_fit8(nrc, rofi, ffit, spl)
    integer,        intent(in)  :: nrc
    real(8),        intent(in)  :: rofi(:)
    real(8),        intent(in)  :: ffit(:)
    type(spline_t), intent(out) :: spl

    !No PUSH SUB, called too often

    spl%x_limit(1) = rofi(1)
    spl%x_limit(2) = rofi(nrc)
    call oct_spline_fit(nrc, rofi(1), ffit(1), spl%spl, spl%acc)

  end subroutine spline_fit8


  !------------------------------------------------------------
  subroutine spline_fit4(nrc, rofi, ffit, spl)
    integer,        intent(in)    :: nrc
    real(4),        intent(in)    :: rofi(:)
    real(4),        intent(in)    :: ffit(:)
    type(spline_t), intent(out)   :: spl

    real(8), allocatable :: rofi8(:), ffit8(:)

    PUSH_SUB(spline_fit4)

    SAFE_ALLOCATE(rofi8(1:nrc))
    SAFE_ALLOCATE(ffit8(1:nrc))

    rofi8 = real(rofi, kind=8)
    ffit8 = real(ffit, kind=8)

    call spline_fit8(nrc, rofi8, ffit8, spl)

    SAFE_DEALLOCATE_A(rofi8)
    SAFE_DEALLOCATE_A(ffit8)

    POP_SUB(spline_fit4)
  end subroutine spline_fit4


  !------------------------------------------------------------
  real(8) pure function spline_eval8(spl, x)
    type(spline_t), intent(in) :: spl
    real(8),        intent(in) :: x

    spline_eval8 = oct_spline_eval(x, spl%spl, spl%acc)
  end function spline_eval8


  !------------------------------------------------------------
  real(4) pure function spline_eval4(spl, x)
    type(spline_t), intent(in) :: spl
    real(4),        intent(in) :: x

    spline_eval4 = real(oct_spline_eval(real(x, kind=8), spl%spl, spl%acc), kind=4)
  end function spline_eval4


  !------------------------------------------------------------
  pure subroutine spline_eval8_array(spl, nn, xf)
    type(spline_t), intent(in)    :: spl
    integer,        intent(in)    :: nn
    real(8),        intent(inout) :: xf(:)

    call oct_spline_eval_array(nn, xf(1), spl%spl, spl%acc)
  end subroutine spline_eval8_array


  !------------------------------------------------------------
  pure subroutine spline_eval4_array(spl, nn, xf)
    type(spline_t), intent(in)    :: spl
    integer,        intent(in)    :: nn
    real(4),        intent(inout) :: xf(:)

    call oct_spline_eval_array4(nn, xf(1), spl%spl, spl%acc)
  end subroutine spline_eval4_array


  !------------------------------------------------------------
  subroutine spline_sum(spl1, spl2, splsum)
    type(spline_t), intent(in)  :: spl1
    type(spline_t), intent(in)  :: spl2
    type(spline_t), intent(out) :: splsum
    integer :: npoints, i

    real(8), allocatable :: x(:), y(:), y2(:)

    PUSH_SUB(spline_sum)

    npoints = oct_spline_npoints(spl1%spl)

    SAFE_ALLOCATE( x(1:npoints))
    SAFE_ALLOCATE( y(1:npoints))
    SAFE_ALLOCATE(y2(1:npoints))

    call oct_spline_x(spl1%spl, x(1))
    call oct_spline_y(spl1%spl, y(1))

    do i = 1, npoints
      y2(i) = spline_eval8(spl2, x(i))
    end do

    y2 = y2 + y
    call oct_spline_fit(npoints, x(1), y2(1), splsum%spl, splsum%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    SAFE_DEALLOCATE_A(y2)

    POP_SUB(spline_sum)
  end subroutine spline_sum


  !------------------------------------------------------------
  subroutine spline_times(a, spl)
    FLOAT,          intent(in)     :: a
    type(spline_t), intent(inout)  :: spl

    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_times)

    npoints = oct_spline_npoints(spl%spl)
    SAFE_ALLOCATE(x(1:npoints))
    SAFE_ALLOCATE(y(1:npoints))

    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))
    call oct_spline_end(spl%spl, spl%acc)
    do i = 1, npoints
      y(i) = a*y(i)
    end do
    call oct_spline_fit(npoints, x(1), y(1), spl%spl, spl%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_times)
  end subroutine spline_times


  !------------------------------------------------------------
  real(8) function spline_integral_full(spl) result(res)
    type(spline_t), intent(in) :: spl

    integer :: npoints
    real(8), allocatable :: x(:)

    PUSH_SUB(spline_integral_full)

    npoints = oct_spline_npoints(spl%spl)
    SAFE_ALLOCATE(x(1:npoints))
    call oct_spline_x(spl%spl, x(1))
    res = oct_spline_eval_integ(spl%spl, x(1), x(npoints), spl%acc)
    SAFE_DEALLOCATE_A(x)

    POP_SUB(spline_integral_full)
  end function spline_integral_full


  !------------------------------------------------------------
  real(8) pure function spline_integral_limits(spl, a, b) result(res)
    type(spline_t), intent(in) :: spl
    real(8),        intent(in) :: a
    real(8),        intent(in) :: b

    res = oct_spline_eval_integ(spl%spl, a, b, spl%acc)
  end function spline_integral_limits


  !------------------------------------------------------------
  real(8) function spline_dotp(spl1, spl2) result (res)
    type(spline_t), intent(in) :: spl1
    type(spline_t), intent(in) :: spl2
    
    type(spline_t) :: aux
    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_dotp)

    npoints = oct_spline_npoints(spl1%spl)
    SAFE_ALLOCATE(x(1:npoints))
    SAFE_ALLOCATE(y(1:npoints))

    call oct_spline_x(spl1%spl, x(1))
    call oct_spline_y(spl1%spl, y(1))
    do i = 1, npoints
      y(i) = y(i)*oct_spline_eval(x(i), spl2%spl, spl2%acc)
    end do
    call spline_init(aux)
    call spline_fit(npoints, x, y, aux)
    res = oct_spline_eval_integ(aux%spl, x(1), x(npoints), aux%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_dotp)
  end function spline_dotp


  !------------------------------------------------------------
  subroutine spline_3dft(spl, splw, gmax)
    type(spline_t),      intent(in)    :: spl
    type(spline_t),      intent(inout) :: splw
    FLOAT, optional,     intent(in)    :: gmax

    type(spline_t) :: aux
    real(8) :: g, dg
    integer :: np
    integer :: npoints, i, j
    real(8), allocatable :: x(:), y(:), y2(:), xw(:), yw(:)

    PUSH_SUB(spline_3dft)

    npoints = oct_spline_npoints(spl%spl)
    SAFE_ALLOCATE( x(1:npoints))
    SAFE_ALLOCATE( y(1:npoints))
    SAFE_ALLOCATE(y2(1:npoints))

    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))

    ! Check whether splw comes with a defined grid, or else build it.
    if(c_associated(splw%spl)) then
      np = oct_spline_npoints(splw%spl)
      SAFE_ALLOCATE(xw(1:np))
      SAFE_ALLOCATE(yw(1:np))
      call oct_spline_x(splw%spl, xw(1))
      ! But now we need to kill the input:
      call spline_end(splw)
    else
      np = 200 ! hard coded value
      dg = gmax/(np-1)
      SAFE_ALLOCATE(xw(1:np))
      SAFE_ALLOCATE(yw(1:np))
      do i = 1, np
        g = (i-1)*dg
        xw(i) = g
      end do
    end if

    ! The first point, xw(1) = 0.0 and it has to be treated separately.
    do j = 1, npoints
      y2(j) = CNST(4.0)*M_PI*y(j)*x(j)**2
    end do
    call spline_init(aux)
    call oct_spline_fit(npoints, x(1), y2(1), aux%spl, aux%acc)
    yw(1) = oct_spline_eval_integ(aux%spl, x(1), x(npoints), aux%acc)
    call spline_end(aux)

    do i = 2, np
      do j = 1, npoints
        y2(j) = (CNST(4.0)*M_PI/xw(i))*y(j)*x(j)*sin(xw(i)*x(j))
      end do
      call spline_init(aux)
      call oct_spline_fit(npoints, x(1), y2(1), aux%spl, aux%acc)
      yw(i) = oct_spline_eval_integ(aux%spl, x(1), x(npoints), aux%acc)
      call spline_end(aux)
    end do

    call spline_init(splw)
    call oct_spline_fit(np, xw(1), yw(1), splw%spl, splw%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    SAFE_DEALLOCATE_A(y2)
    SAFE_DEALLOCATE_A(xw)
    SAFE_DEALLOCATE_A(yw)

    POP_SUB(spline_3dft)
  end subroutine spline_3dft


  !------------------------------------------------------------
  subroutine spline_besselft(spl, splw, l, gmax)
    type(spline_t),    intent(in)    :: spl
    type(spline_t),    intent(inout) :: splw
    integer,           intent(in)    :: l
    FLOAT,   optional, intent(in)    :: gmax

    type(spline_t) :: aux
    real(8) :: g, dg
    integer :: np
    integer :: npoints, i, j
    real(8), allocatable :: x(:), y(:), y2(:), xw(:), yw(:)

    PUSH_SUB(spline_besselft)

    npoints = oct_spline_npoints(spl%spl)
    SAFE_ALLOCATE( x(1:npoints))
    SAFE_ALLOCATE( y(1:npoints))
    SAFE_ALLOCATE(y2(1:npoints))

    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))

    ! Check whether splw comes with a defined grid, or else build it.
    if(c_associated(splw%spl)) then
      np = oct_spline_npoints(splw%spl)
      SAFE_ALLOCATE(xw(1:np))
      SAFE_ALLOCATE(yw(1:np))
      call oct_spline_x(splw%spl, xw(1))
      ! But now we need to kill the input:
      call spline_end(splw)
    else
      ASSERT(present(gmax))
      np = 1000 ! hard coded value
      dg = gmax/(np-1)
      SAFE_ALLOCATE(xw(1:np))
      SAFE_ALLOCATE(yw(1:np))
      do i = 1, np
        g = real(i-1, 8)*dg
        xw(i) = g
      end do
    end if

    do i = 1, np
      do j = 1, npoints
        y2(j) = y(j) * x(j)**2 * loct_sph_bessel(l, x(j)*xw(i))
      end do
      call spline_init(aux)
      call oct_spline_fit(npoints, x(1), y2(1), aux%spl, aux%acc)
      yw(i) = sqrt(CNST(2.0)/M_PI)*oct_spline_eval_integ(aux%spl, x(1), x(npoints), aux%acc)
      call spline_end(aux)
    end do

    call spline_init(splw)
    call oct_spline_fit(np, xw(1), yw(1), splw%spl, splw%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    SAFE_DEALLOCATE_A(y2)
    SAFE_DEALLOCATE_A(xw)
    SAFE_DEALLOCATE_A(yw)

    POP_SUB(spline_besselft)
  end subroutine spline_besselft


  !------------------------------------------------------------
  subroutine spline_cut(spl, cutoff, beta)
    type(spline_t), intent(inout) :: spl
    FLOAT,          intent(in)    :: cutoff
    FLOAT,          intent(in)    :: beta

    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_cut)

    npoints = oct_spline_npoints(spl%spl)
    SAFE_ALLOCATE(x(1:npoints))
    SAFE_ALLOCATE(y(1:npoints))

    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))
    call oct_spline_end(spl%spl, spl%acc)
    do i = npoints, 1, -1
      if(x(i)<cutoff) then
        exit
      endif
      y(i) = y(i) * exp(-beta*(x(i)/cutoff - CNST(1.0))**2)
    end do
    call oct_spline_fit(npoints, x(1), y(1), spl%spl, spl%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_cut)
  end subroutine spline_cut


  !------------------------------------------------------------
  subroutine spline_div(spla, splb)
    type(spline_t),   intent(inout) :: spla
    type(spline_t),   intent(in)    :: splb

    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)
    real(8) :: aa

    PUSH_SUB(spline_div)

    npoints = oct_spline_npoints(spla%spl)

    SAFE_ALLOCATE(x(1:npoints))
    SAFE_ALLOCATE(y(1:npoints))

    call oct_spline_x(spla%spl, x(1))
    call oct_spline_y(spla%spl, y(1))
    call oct_spline_end(spla%spl, spla%acc)
  
    do i = npoints, 1, -1
      if(x(i) > splb%x_limit(2)) cycle
      aa = spline_eval(splb, x(i))
      y(i) = y(i)/aa
    end do
    
    call oct_spline_fit(npoints, x(1), y(1), spla%spl, spla%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_div)
  end subroutine spline_div


  !------------------------------------------------------------
  subroutine spline_mult(spla, splb)
    type(spline_t),  intent(inout) :: spla
    type(spline_t),  intent(in)    :: splb

    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)
    real(8) :: aa

    PUSH_SUB(spline_mult)

    npoints = oct_spline_npoints(spla%spl)

    SAFE_ALLOCATE(x(1:npoints))
    SAFE_ALLOCATE(y(1:npoints))

    call oct_spline_x(spla%spl, x(1))
    call oct_spline_y(spla%spl, y(1))
    call oct_spline_end(spla%spl, spla%acc)
  
    do i = npoints, 1, -1
      if(x(i) > splb%x_limit(2)) then
        aa = M_ZERO
      else
        aa = spline_eval(splb, x(i))
      end if
      y(i) = y(i)*aa
    end do
    
    call oct_spline_fit(npoints, x(1), y(1), spla%spl, spla%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_mult)
  end subroutine spline_mult


  !------------------------------------------------------------
  subroutine spline_der(spl, dspl)
    type(spline_t), intent(in)    :: spl
    type(spline_t), intent(inout) :: dspl

    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_der)

    ! Use the grid of dspl if it is present, otherwise use the same one of spl.
    if(.not. c_associated(dspl%spl)) then ! use the grid of spl
      npoints = oct_spline_npoints(spl%spl)
      SAFE_ALLOCATE(x(1:npoints))
      SAFE_ALLOCATE(y(1:npoints))
      call oct_spline_x(spl%spl, x(1))
    else ! use the grid of dspl
      npoints = oct_spline_npoints(dspl%spl)
      SAFE_ALLOCATE(x(1:npoints))
      SAFE_ALLOCATE(y(1:npoints))
      call oct_spline_x(dspl%spl, x(1))
    end if
    do i = 1, npoints
      y(i) = oct_spline_eval_der(x(i), spl%spl, spl%acc)
    end do
    call oct_spline_fit(npoints, x(1), y(1), dspl%spl, dspl%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_der)
  end subroutine spline_der


  !------------------------------------------------------------
  subroutine spline_der2(spl, dspl)
    type(spline_t), intent(in)    :: spl
    type(spline_t), intent(inout) :: dspl

    integer :: npoints, i
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_der2)

    ! Use the grid of dspl if it is present, otherwise use the same one of spl.
    if(.not. c_associated(dspl%spl)) then ! use the grid of spl
      npoints = oct_spline_npoints(spl%spl)
      SAFE_ALLOCATE(x(1:npoints))
      SAFE_ALLOCATE(y(1:npoints))
      call oct_spline_x(spl%spl, x(1))
    else ! use the grid of dspl
      npoints = oct_spline_npoints(dspl%spl)
      SAFE_ALLOCATE(x(1:npoints))
      SAFE_ALLOCATE(y(1:npoints))
      call oct_spline_x(dspl%spl, x(1))
    end if
    do i = 1, npoints
      y(i) = oct_spline_eval_der2(x(i), spl%spl, spl%acc)
    end do
    call oct_spline_fit(npoints, x(1), y(1), dspl%spl, dspl%acc)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_der2)
  end subroutine spline_der2


  !------------------------------------------------------------
  subroutine spline_print_0(spl, iunit)
    type(spline_t), intent(in) :: spl
    integer,        intent(in) :: iunit

    integer :: np, i
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_print_0)

    np = oct_spline_npoints(spl%spl)
    SAFE_ALLOCATE(x(1:np))
    SAFE_ALLOCATE(y(1:np))

    call oct_spline_x(spl%spl, x(1))
    call oct_spline_y(spl%spl, y(1))
    do i = 1, np
      write(iunit, '(2f16.8)') x(i), y(i)
    end do

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_print_0)
  end subroutine spline_print_0


  !------------------------------------------------------------
  subroutine spline_print_1(spl, iunit)
    type(spline_t), intent(in) :: spl(:)
    integer,        intent(in) :: iunit

    character(len=4)  :: fm
    integer :: np, i, n, j
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_print_1)

    n = size(spl)
    if(n<=0) then
      POP_SUB(spline_print_1)
      return
    endif

    write(fm,'(i4)') n + 1; fm = adjustl(fm)
    np = oct_spline_npoints(spl(1)%spl)
    SAFE_ALLOCATE(x(1:np))
    SAFE_ALLOCATE(y(1:np))

    call oct_spline_x(spl(1)%spl, x(1))
    call oct_spline_y(spl(1)%spl, y(1))
    do i = 1, np
      write(iunit, '('//trim(fm)//'f16.8)') x(i), (spline_eval(spl(j), x(i)), j = 1, size(spl))
    end do

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_print_1)
  end subroutine spline_print_1


  !------------------------------------------------------------
  subroutine spline_print_2(spl, iunit)
    type(spline_t), intent(in) :: spl(:, :)
    integer,        intent(in) :: iunit

    character(len=4)  :: fm
    integer :: np, i, n1, n2, j, k
    real(8), allocatable :: x(:), y(:)

    PUSH_SUB(spline_print_2)

    n1 = size(spl, 1); n2 = size(spl, 2)
    if(n1*n2<=0) then
      POP_SUB(spline_print_2)
      return
    endif

    write(fm,'(i4)') n1*n2 + 1; fm = adjustl(fm)
    np = oct_spline_npoints(spl(1, 1)%spl)

    SAFE_ALLOCATE(x(1:np))
    SAFE_ALLOCATE(y(1:np))

    call oct_spline_x(spl(1, 1)%spl, x(1))
    call oct_spline_y(spl(1, 1)%spl, y(1))
    do i = 1, np
      write(iunit, '('//trim(fm)//'f16.8)') x(i), &
        ((spline_eval(spl(j, k), x(i)), j = 1, n1), k = 1, n2)
    end do

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)

    POP_SUB(spline_print_2)
  end subroutine spline_print_2

 
  !------------------------------------------------------------
  FLOAT pure function spline_cutoff_radius(spl, threshold) result(r)
    type(spline_t), intent(in) :: spl
    FLOAT,          intent(in) :: threshold

    integer :: ii, jj
    FLOAT, parameter :: dx = CNST(0.01)
    
    jj = int(spl%x_limit(2)/dx) + 1

    do ii = jj, 1, -1
      r = dx*(ii-1)
      if(abs(spline_eval(spl, r)) > threshold) exit
    end do

  end function spline_cutoff_radius

end module splines_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
