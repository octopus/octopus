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

module minimizer_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use profiling_oct_m
  use messages_oct_m
  use mpi_oct_m

  implicit none

  private
  public ::                    &
    loct_1dminimize,           &
    minimize_multidim,         &
    minimize_multidim_nograd,  &
    minimize_multidim_nlopt


  integer, public, parameter ::      &
    MINMETHOD_STEEPEST_DESCENT =  1, &
    MINMETHOD_FR_CG            =  2, &
    MINMETHOD_PR_CG            =  3, &
    MINMETHOD_BFGS             =  4, &
    MINMETHOD_BFGS2            =  5, &
    MINMETHOD_NMSIMPLEX        =  6, &
    MINMETHOD_SD_NATIVE        = -1, &
    MINMETHOD_NLOPT_BOBYQA     =  7, &
    MINMETHOD_FIRE             =  8, &
    MINMETHOD_NLOPT_LBFGS      =  9

  interface loct_1dminimize
    subroutine oct_1dminimize(a, b, m, f, status)
      implicit none
      real(8), intent(inout) :: a, b, m
      interface
        subroutine f(x, fx)
          implicit none
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
      implicit none
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
          implicit none
          integer, intent(in)    :: n
          real(8), intent(in)    :: x(n)
          real(8), intent(inout) :: val
          integer, intent(in)    :: getgrad
          real(8), intent(inout) :: grad(n)
        end subroutine f
        subroutine write_iter_info(iter, n, val, maxdr, maxgrad, x)
          implicit none
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
     implicit none
      integer :: oct_minimize_direct
      integer, intent(in)    :: method
      integer, intent(in)    :: dim
      real(8), intent(inout) :: x
      real(8), intent(in)    :: step
      real(8), intent(in)    :: toldr
      integer, intent(in)    :: maxiter
      !> No intents here is unfortunately required because the same dummy function will be passed
      !! also to newuoa routines in opt_control, and there the interface has no intents.
      !! UPDATE: The newuoa interfaces are gone, so probably this can be fixed.
      interface
        subroutine f(n, x, val)
          implicit none
          integer :: n
          real(8) :: x(n)
          real(8) :: val
        end subroutine f
        subroutine write_iter_info(iter, n, val, maxdr, x)
          implicit none
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
    !! UPDATE: The newoua interface is gone, and therefore probably this can be fixed.
    interface
      subroutine f(n, x, val)
        implicit none
        integer :: n
        real(8) :: x(n)
        real(8) :: val
      end subroutine f
      subroutine write_iter_info(iter, n, val, maxdr, x)
        implicit none
        integer, intent(in) :: iter
        integer, intent(in) :: n
        real(8), intent(in) :: val
        real(8), intent(in) :: maxdr
        real(8), intent(in) :: x(n)
      end subroutine write_iter_info
    end interface
    real(8), intent(out)   :: minimum
    integer, intent(out)   :: ierr
    
    PUSH_SUB(minimize_multidim_nograd)

    ASSERT(ubound(x, dim = 1) >= dim)
    
    select case(method)
    case(MINMETHOD_NMSIMPLEX) 
      ierr = loct_minimize_direct(method, dim, x(1), step, toldr, maxiter, f, write_iter_info, minimum)
    end select

    POP_SUB(minimize_multidim_nograd)

  end subroutine minimize_multidim_nograd


  subroutine minimize_multidim_nlopt(ierr, method, dim, x, step, toldr, maxiter, f, minimum, lb, ub)
    integer, intent(out)   :: ierr
    integer, intent(in)    :: method
    integer, intent(in)    :: dim
    real(8), intent(inout) :: x(:)
    real(8), intent(in)    :: step
    real(8), intent(in)    :: toldr
    integer, intent(in)    :: maxiter
    interface
      subroutine f(val, n, x, grad, need_gradient, f_data)
        FLOAT    :: val
        integer  :: n
        FLOAT    :: x(n), grad(n)
        integer  :: need_gradient
        FLOAT    :: f_data
      end subroutine f
    end interface
    real(8), intent(out)   :: minimum
    real(8), intent(in), optional :: lb(:), ub(:)
#if defined(HAVE_NLOPT)

    integer(8) :: opt
    integer :: ires
    include 'nlopt.f'

    select case(method)
    case(MINMETHOD_NLOPT_BOBYQA)
      call nlo_create(opt, NLOPT_LN_BOBYQA, dim)
    case(MINMETHOD_NLOPT_LBFGS)
      call nlo_create(opt, NLOPT_LD_LBFGS, dim)
    end select

    if(present(lb)) then
      call nlo_set_lower_bounds(ires, opt, lb)
    end if
    if(present(ub)) then
      call nlo_set_upper_bounds(ires, opt, ub)
    end if

    call nlo_set_min_objective(ires, opt, f, 0)
    ! This would set an inequality constraint (TODO)
    !call nlo_add_inequality_constraint(ires, opt, myconstraint, d1, CNST(1.0e-8))

    call nlo_set_xtol_abs(ires, opt, toldr)
    call nlo_set_initial_step1(ires, opt, step)
    call nlo_set_maxeval(ires, opt, maxiter)

    call nlo_optimize(ires, opt, x, minimum)
    ierr = ires
    call nlo_destroy(opt)
#else
    ierr = 0
#endif
  end subroutine minimize_multidim_nlopt


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
        implicit none
        integer, intent(in)    :: n
        real(8), intent(in)    :: x(n)
        real(8), intent(inout) :: val
        integer, intent(in)    :: getgrad
        real(8), intent(inout) :: grad(n)
      end subroutine f
      subroutine write_iter_info(iter, n, val, maxdr, maxgrad, x)
        implicit none
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
      call minimize_sd(dim, x, step, maxiter, f, write_iter_info, minimum, ierr)
      
    case default
      ierr = loct_minimize(method, dim, x(1), step, line_tol, tolgrad, toldr, maxiter, f, write_iter_info, minimum)
      
    end select

    POP_SUB(minimize_multidim)
    
  end subroutine minimize_multidim

  !----------------------------------------------

  subroutine minimize_sd(dim, x, step, maxiter, f, write_iter_info, minimum, ierr)
    integer, intent(in)    :: dim
    real(8), intent(inout) :: x(:)
    real(8), intent(in)    :: step
    integer, intent(in)    :: maxiter
    interface
      subroutine f(n, x, val, getgrad, grad)
        implicit none
        integer, intent(in)    :: n
        real(8), intent(in)    :: x(n)
        real(8), intent(inout) :: val
        integer, intent(in)    :: getgrad
        real(8), intent(inout) :: grad(n)
      end subroutine f
      subroutine write_iter_info(iter, n, val, maxdr, maxgrad, x)
        implicit none
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

end module minimizer_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
