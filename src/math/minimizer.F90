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
  use iso_c_binding
  use lalg_basic_oct_m
  use profiling_oct_m
  use messages_oct_m

  implicit none

  private
  public ::                    &
    loct_1dminimize,           &
    minimize_fire,             &
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
        use iso_c_binding
        real(c_double), intent(out) :: val
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(*)
        real(c_double), intent(out) :: grad(*)
        integer(c_int), intent(in)  :: need_gradient
        type(c_ptr),    intent(in)  :: f_data
      end subroutine f
    end interface
    real(8), intent(out)   :: minimum
    real(8), intent(in), optional :: lb(:), ub(:)
#if defined(HAVE_NLOPT)

    interface
      subroutine nlocreate(opt, alg, n)
        use iso_c_binding
        type(c_ptr),    intent(out) :: opt
        integer(c_int), intent(in)  :: alg
        integer(c_int), intent(in)  :: n
      end subroutine nlocreate

      subroutine nlo_set_lower_bounds(ret, opt, lower_bounds)
        use iso_c_binding
        integer(c_int), intent(out)   :: ret
        type(c_ptr),    intent(inout) :: opt
        real(c_double), intent(in)    :: lower_bounds(*)
      end subroutine nlo_set_lower_bounds

      subroutine nlo_set_upper_bounds(ret, opt, upper_bounds)
        use iso_c_binding
        integer(c_int), intent(out)   :: ret
        type(c_ptr),    intent(inout) :: opt
        real(c_double), intent(in)    :: upper_bounds(*)
      end subroutine nlo_set_upper_bounds

      subroutine nlo_set_min_objective(ret, opt, f, f_data)
        use iso_c_binding
        integer(c_int), intent(out)   :: ret
        type(c_ptr),    intent(inout) :: opt
        interface
          subroutine f(val, n, x, grad, need_gradient, f_data)
            use iso_c_binding
            real(c_double), intent(out) :: val
            integer(c_int), intent(in)  :: n
            real(c_double), intent(in)  :: x(*)
            real(c_double), intent(out) :: grad(*)
            integer(c_int), intent(in)  :: need_gradient
            type(c_ptr),    intent(in)  :: f_data
          end subroutine f
        end interface
        type(c_ptr),    intent(in)    :: f_data
      end subroutine nlo_set_min_objective

      subroutine nlo_set_xtol_abs1(ret, opt, xtol_abs)
        use iso_c_binding
        integer(c_int), intent(out)   :: ret
        type(c_ptr),    intent(inout) :: opt
        real(c_double), intent(in)    :: xtol_abs
      end subroutine nlo_set_xtol_abs1

      subroutine nlo_set_initial_step1(ret, opt, initial_step1)
        use iso_c_binding
        integer(c_int), intent(out)   :: ret
        type(c_ptr),    intent(inout) :: opt
        real(c_double), intent(in)    :: initial_step1
      end subroutine nlo_set_initial_step1

      subroutine nlo_set_maxeval(ret, opt, maxeval)
        use iso_c_binding
        integer(c_int), intent(out)   :: ret
        type(c_ptr),    intent(inout) :: opt
        integer(c_int), intent(in)    :: maxeval
      end subroutine nlo_set_maxeval

      subroutine nlo_optimize(ret, opt, x, optf)
        use iso_c_binding
        integer(c_int), intent(out)   :: ret
        type(c_ptr),    intent(inout) :: opt
        real(c_double), intent(inout) :: x(*)
        real(c_double), intent(out)   :: optf
      end subroutine nlo_optimize

      subroutine nlo_destroy(opt)
        use iso_c_binding
        type(c_ptr), intent(inout) :: opt
      end subroutine nlo_destroy
    end interface

    type(c_ptr) :: opt
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

    call nlo_set_min_objective(ires, opt, f, C_NULL_PTR)
    ! This would set an inequality constraint (TODO)
    !call nlo_add_inequality_constraint(ires, opt, myconstraint, d1, CNST(1.0e-8))

    call nlo_set_xtol_abs1(ires, opt, toldr)
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

  !----------------------------------------------

  subroutine minimize_fire(dim, x, step, tolgrad, maxiter, f, write_iter_info, en, ierr, mass, integrator)
    integer, intent(in)    :: dim
    real(8), intent(inout) :: x(:)
    real(8), intent(in)    :: step
    real(8), intent(in)    :: tolgrad
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
    real(8), intent(out)   :: en
    integer, intent(out)   :: ierr
    real(8), intent(in)    :: mass(:)
    integer, intent(in)    :: integrator

    integer :: n_iter
    real(8), allocatable :: grad(:)
    real(8) :: dt

    integer :: n_min
    real (8) :: alpha
    real (8) :: alpha_start
    real (8) :: f_alpha
    real (8) :: p_value
    real (8) :: f_inc
    real (8) :: f_dec
    real (8) :: dt_max

    integer :: p_times

    real (8), allocatable :: grad_atoms(:)
    real (8), allocatable :: vel(:)
    real (8), allocatable :: dr_i(:)
    real (8), allocatable :: x_new(:)
    real (8), allocatable :: dr_atoms(:)

    integer :: i_tmp

    real (8) :: mod_dr
    real (8) :: maxmove

    PUSH_SUB(minimize_fire)

    SAFE_ALLOCATE(grad_atoms(1:dim/3))
    SAFE_ALLOCATE(grad(1:dim))

    alpha_start = CNST(0.1)

    dt = step

    alpha = alpha_start

    p_times = 0

    f_alpha = CNST(0.99)
    n_min = 5
    f_inc = CNST(1.1)
    dt_max = CNST(10.0) * dt
    f_dec = CNST(0.5)

    maxmove = CNST(0.2) * P_Ang

    grad = M_ZERO

    SAFE_ALLOCATE(vel(1:dim))
    vel = M_ZERO

    SAFE_ALLOCATE(dr_atoms(1:dim/3))
    SAFE_ALLOCATE(x_new(1:dim))
    SAFE_ALLOCATE(dr_i(1:dim))

    x_new = x
    dr_i = M_ZERO

    n_iter = 1

    do while (n_iter <= maxiter)

      call f(dim, x_new, en, 1, grad)

      select case (integrator)
        ! Velocity verlet
      case (OPTION__GOFIREINTEGRATOR__VERLET)
        vel(1:dim) = vel(1:dim) - M_HALF*grad(1:dim)*dt/mass(1:dim)
      end select

      if (n_iter /= 1) then
        p_value = M_ZERO
        do i_tmp = 0, dim/3 - 1
          p_value = p_value - grad(3*i_tmp+1)*vel(3*i_tmp+1) - grad(3*i_tmp+2)*vel(3*i_tmp+2) - grad(3*i_tmp+3)*vel(3*i_tmp+3)
        end do

        if(p_value > M_ZERO) then
          vel(1:dim) = (M_ONE - alpha) * vel(1:dim) - alpha * grad(1:dim) * lalg_nrm2(dim,vel) / lalg_nrm2(dim,grad)
          if(p_times > n_min) then
            dt = min(dt * f_inc , dt_max)
            alpha = alpha * f_alpha
          end if

        else
          p_times = 0
          dt = dt * f_dec
          alpha = alpha_start
          vel = M_ZERO
        end if
      end if

      x(1:dim)=x_new(1:dim)


      select case (integrator)
      case (OPTION__GOFIREINTEGRATOR__VERLET)
        ! Velocity Verlet
        dr_i = vel(1:dim)*dt - M_HALF*grad(1:dim)*dt**2/mass(1:dim)
        vel(1:dim) = vel(1:dim) - M_HALF*grad(1:dim)*dt/mass(1:dim)
      case (OPTION__GOFIREINTEGRATOR__EULER)
        ! Euler method
        vel(1:dim) = vel(1:dim) - grad(1:dim)*dt/mass(1:dim)
        dr_i(1:dim) = vel(1:dim)*dt
      end select

      mod_dr = lalg_nrm2(dim, dr_i)
      if (mod_dr > maxmove) then
        dr_i = maxmove * dr_i / mod_dr
      end if

      x_new(1:dim) = x(1:dim) + dr_i(1:dim)

      do i_tmp = 0, dim/3 - 1
        grad_atoms(i_tmp+1) = sqrt(grad(3*i_tmp+1)**2 + grad(3*i_tmp+2)**2 + grad(3*i_tmp+3)**2)
        dr_atoms(i_tmp+1) =   sqrt(dr_i(3*i_tmp+1)**2+dr_i(3*i_tmp+2)**2+dr_i(3*i_tmp+3)**2)
      end do
      call write_iter_info(n_iter, dim, en, maxval(dr_atoms(1:)), maxval(abs(grad_atoms(1:))), x)

      if(maxval(abs(grad_atoms(1:))) < tolgrad) then
        ierr = 0
        n_iter = maxiter + 1
      else
        n_iter = n_iter + 1
      end if

    end do

    SAFE_DEALLOCATE_A(dr_atoms)
    SAFE_DEALLOCATE_A(x_new)
    SAFE_DEALLOCATE_A(dr_i)
    SAFE_DEALLOCATE_A(vel)
    SAFE_DEALLOCATE_A(grad)
    SAFE_DEALLOCATE_A(grad_atoms)

    POP_SUB(minimize_fire)

  end subroutine minimize_fire

end module minimizer_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
