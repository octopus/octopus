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
!! $Id$

#include "global.h"

module minimizer_m
  use global_m
  use profiling_m
  use messages_m
#if defined(HAVE_NEWUOA)
  use newuoa_m
#endif

  implicit none

  private
  public ::                    &
    loct_1dminimize,           &
    minimize_multidim,         &
    minimize_multidim_nograd

  integer, public, parameter ::      &
    MINMETHOD_STEEPEST_DESCENT =  1, &
    MINMETHOD_FR_CG            =  2, &
    MINMETHOD_PR_CG            =  3, &
    MINMETHOD_BFGS             =  4, &
    MINMETHOD_BFGS2            =  5, &
    MINMETHOD_NMSIMPLEX        =  6, &
    MINMETHOD_SD_NATIVE        = -1, &
    MINMETHOD_NEWUOA           =  7

  interface loct_1dminimize
    subroutine oct_1dminimize(a, b, m, f, status)
      real(8), intent(inout) :: a, b, m
      interface
        subroutine f(x, fx)
          real(8), intent(in)  :: x
          real(8), intent(out) :: fx
        end subroutine f
      end interface
      integer, intent(out) :: status
    end subroutine oct_1dminimize
  end interface loct_1dminimize

  interface loct_minimize
    integer function oct_minimize(method, dim, x, step, line_tol, &
      tolgrad, toldr, maxiter, f, write_iter_info, minimum)
      integer, intent(in)    :: method
      integer, intent(in)    :: dim
      real(8), intent(inout) :: x
      real(8), intent(in)    :: step
      real(8), intent(in)    :: line_tol
      real(8), intent(in)    :: tolgrad
      real(8), intent(in)    :: toldr
      integer, intent(in)    :: maxiter
      interface
        subroutine f(n, x, val, getgrad, grad)
          integer, intent(in)    :: n
          real(8), intent(in)    :: x(n)
          real(8), intent(inout) :: val
          integer, intent(in)    :: getgrad
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
      real(8), intent(out)   :: minimum
    end function oct_minimize
  end interface loct_minimize

  interface loct_minimize_direct
    function oct_minimize_direct(method, dim, x, step, toldr, maxiter, f, write_iter_info, minimum)
      integer :: oct_minimize_direct
      integer, intent(in)    :: method
      integer, intent(in)    :: dim
      real(8), intent(inout) :: x
      real(8), intent(in)    :: step
      real(8), intent(in)    :: toldr
      integer, intent(in)    :: maxiter
      !> No intents here is unfortunately required because the same dummy function will be passed
      !! also to newuoa routines in opt_control, and there the interface has no intents.
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
      real(8), intent(out)   :: minimum
    end function oct_minimize_direct
  end interface loct_minimize_direct

contains 

  subroutine minimize_multidim_nograd(method, dim, x, step, toldr, maxiter, f, write_iter_info, minimum, ierr)
    integer, intent(in)    :: method
    integer, intent(in)    :: dim
    real(8), intent(inout) :: x(:)
    real(8), intent(in)    :: step
    real(8), intent(in)    :: toldr
    integer, intent(in)    :: maxiter
    !> No intents here is unfortunately required because the same dummy function will be passed
    !! also to newuoa routines in opt_control, and there the interface has no intents.
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
    real(8), intent(out)   :: minimum
    integer, intent(out)   :: ierr
    
    integer :: npt, iprint, sizeofw
    REAL_DOUBLE, allocatable :: w(:)
    
    PUSH_SUB(minimize_multidim_nograd)

    ASSERT(ubound(x, dim = 1) >= dim)
    
    select case(method)
    case(MINMETHOD_NMSIMPLEX) 
      ierr = loct_minimize_direct(method, dim, x(1), step, toldr, maxiter, f, write_iter_info, minimum)
#if defined(HAVE_NEWUOA)
    case(MINMETHOD_NEWUOA)
      npt = 2*dim + 1
      iprint = 2
      sizeofw = (npt + 13)*(npt + dim) + 3 * dim*(dim + 3)/2 
      SAFE_ALLOCATE(w(1:sizeofw))
      w = M_ZERO
      call newuoa(dim, npt, x, step, toldr, iprint, maxiter, w, f)
      call f(dim, x, minimum)
      SAFE_DEALLOCATE_A(w)
      ierr = 0
#endif
    end select

    POP_SUB(minimize_multidim_nograd)

  end subroutine minimize_multidim_nograd

  !----------------------------------------------

  subroutine minimize_multidim(method, dim, x, step, line_tol, tolgrad, toldr, maxiter, f, write_iter_info, minimum, ierr)
    integer, intent(in)    :: method
    integer, intent(in)    :: dim
    real(8), intent(inout) :: x(:)
    real(8), intent(in)    :: step
    real(8), intent(in)    :: line_tol
    real(8), intent(in)    :: tolgrad
    real(8), intent(in)    :: toldr
    integer, intent(in)    :: maxiter
    interface
      subroutine f(n, x, val, getgrad, grad)
        integer, intent(in)    :: n
        real(8), intent(in)    :: x(n)
        real(8), intent(inout) :: val
        integer, intent(in)    :: getgrad
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
    real(8), intent(out)   :: minimum
    integer, intent(out)   :: ierr
    
    PUSH_SUB(minimize_multidim)

    ASSERT(ubound(x, dim = 1) >= dim)
    
    select case(method)
    case(MINMETHOD_SD_NATIVE)
      call minimize_sd(dim, x, step, tolgrad, toldr, maxiter, f, write_iter_info, minimum, ierr)
      
    case default
      ierr = loct_minimize(method, dim, x(1), step, line_tol, tolgrad, toldr, maxiter, f, write_iter_info, minimum)
      
    end select

    POP_SUB(minimize_multidim)
    
  end subroutine minimize_multidim

  !----------------------------------------------

  subroutine minimize_sd(dim, x, step, tolgrad, toldr, maxiter, f, write_iter_info, minimum, ierr)
    integer, intent(in)    :: dim
    real(8), intent(inout) :: x(:)
    real(8), intent(in)    :: step
    real(8), intent(in)    :: tolgrad
    real(8), intent(in)    :: toldr
    integer, intent(in)    :: maxiter
    interface
      subroutine f(n, x, val, getgrad, grad)
        integer, intent(in)    :: n
        real(8), intent(in)    :: x(n)
        real(8), intent(inout) :: val
        integer, intent(in)    :: getgrad
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
    real(8), intent(out)   :: minimum
    integer, intent(out)   :: ierr
    
    integer :: iter
    real(8), allocatable :: grad(:)
    real(8) :: step2, maxgrad
    
    PUSH_SUB(minimize_sd)

    SAFE_ALLOCATE(grad(1:dim))
    
    step2 = step*CNST(10.0)
    do iter = 1, maxiter
      call f(dim, x, minimum, 1, grad)
      
      maxgrad = maxval(abs(grad))
      
      call write_iter_info(iter, dim, minimum, maxgrad*step2, maxgrad, x)
      
      x(1:dim) = x(1:dim) - step2*grad(1:dim)
      
      step2 = step2*CNST(0.99)
    end do
    ierr = 0

    POP_SUB(minimize_sd)
  end subroutine minimize_sd
  
end module minimizer_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
