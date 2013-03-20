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

module minimizer_m
  implicit none

  private
  public ::                    &
    loct_1dminimize,           &
    minimize_multidim,         &
    minimize_multidim_nograd

  integer, public, parameter ::     &
    MINMETHOD_STEEPEST_DESCENT = 1, &
    MINMETHOD_FR_CG            = 2, &
    MINMETHOD_PR_CG            = 3, &
    MINMETHOD_BFGS             = 4, &
    MINMETHOD_BFGS2            = 5, &
    MINMETHOD_NMSIMPLEX        = 6

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
    integer function oct_minimize(method, dim, x, step, tolgrad, toldr, maxiter, f, write_iter_info, minimum)
      integer, intent(in)    :: method
      integer, intent(in)    :: dim
      real(8), intent(inout) :: x
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
      integer, intent(out)   :: ierr

      ierr = loct_minimize_direct(method, dim, x, step, toldr, maxiter, f, write_iter_info, minimum)

    end subroutine minimize_multidim_nograd

    !----------------------------------------------

    subroutine minimize_multidim(method, dim, x, step, tolgrad, toldr, maxiter, f, write_iter_info, minimum, ierr)
      integer, intent(in)    :: method
      integer, intent(in)    :: dim
      real(8), intent(inout) :: x
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

      ierr = loct_minimize(method, dim, x, step, tolgrad, toldr, maxiter, f, write_iter_info, minimum)

    end subroutine minimize_multidim

end module minimizer_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
