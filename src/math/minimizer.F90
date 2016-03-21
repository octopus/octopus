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
    MINMETHOD_FIRE             =  8


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
    
    integer :: npt, iprint, sizeofw
    REAL_DOUBLE, allocatable :: w(:)
    
    PUSH_SUB(minimize_multidim_nograd)

    ASSERT(ubound(x, dim = 1) >= dim)
    
    select case(method)
    case(MINMETHOD_NMSIMPLEX) 
      ierr = loct_minimize_direct(method, dim, x(1), step, toldr, maxiter, f, write_iter_info, minimum)
    end select

    POP_SUB(minimize_multidim_nograd)

  end subroutine minimize_multidim_nograd


  subroutine minimize_multidim_nlopt(method, dim, x, step, toldr, maxiter, f, minimum, ierr)
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
    integer, intent(out)   :: ierr
#if defined(HAVE_NLOPT)

    integer(8) :: opt
    real(8) :: minf
    integer :: ires
    include 'nlopt.f'
 
    call nlo_create(opt, NLOPT_LN_BOBYQA, dim)

    ! this would set lower and upper bounds (TODO)
    !call nlo_set_lower_bounds(ires, opt, lb)
    !call nlo_set_upper_bounds(ires, opt, lb)
    call nlo_set_min_objective(ires, opt, f, 0)
    ! This would set an inequality constraint (TODO)
    !call nlo_add_inequality_constraint(ires, opt, myconstraint, d1, CNST(1.0e-8))

    call nlo_set_xtol_abs(ires, opt, toldr)
    call nlo_set_initial_step1(ires, opt, step)
    call nlo_set_maxeval(ires, opt, maxiter)

    call nlo_optimize(ires, opt, x, minimum)
    ierr = ires
    call nlo_destroy(opt)

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

  subroutine minimize_fire(dim, x, step, tolgrad, maxiter, f, write_iter_info, en, ierr, mass)
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

    integer :: n_iter
    real(8), allocatable :: grad(:)
    real(8) :: dt
    real(8) :: max_grad_atoms

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
    real (8), allocatable :: vec_delta_pos(:)
    real (8), allocatable :: dr_i(:)
    real (8), allocatable :: x_new(:)
    real (8), allocatable :: dr_atoms(:)

    integer :: dr_atom_iter
    integer :: i_tmp

    real (8) :: mod_vel
    real (8) :: mod_force

    PUSH_SUB(minimize_fire)

    if(mpi_grp_is_root(mpi_world)) then
      call messages_experimental('GOMethod = fire')
    end if

    SAFE_ALLOCATE(grad_atoms(1:dim/3))
    SAFE_ALLOCATE(grad(1:dim))

    alpha_start = CNST(0.1)

    dt = step

    alpha = alpha_start
    
    p_times = 0

    f_alpha = CNST(0.99)
    n_min = 5
    f_inc = CNST(1.1)
    dt_max = 10.0 * dt
    f_dec = CNST(0.5)

    SAFE_ALLOCATE(vec_delta_pos(1:dim))

    grad = 0.0

    SAFE_ALLOCATE(vel(1:dim))
    vel = 0.0

    SAFE_ALLOCATE(dr_atoms(1:dim/3))
    SAFE_ALLOCATE(x_new(1:dim))
    SAFE_ALLOCATE(dr_i(1:dim))

    x_new = x
    dr_i = 0.0

    n_iter = 1
    do while (n_iter <= maxiter)
      
      vel(1:dim) = vel(1:dim) - grad(1:dim)*dt/mass(1:dim)

      call f(dim, x_new, en, 1, grad)
      
      vel(1:dim) = vel(1:dim) - grad(1:dim)*dt/mass(1:dim)
      
      grad_atoms = 0.0

      do i_tmp = 0, dim/3 - 1
        grad_atoms(i_tmp+1) = sqrt(grad(3*i_tmp+1)**2 + grad(3*i_tmp+2)**2 + grad(3*i_tmp+3)**2)
      end do

      max_grad_atoms = maxval(abs(grad_atoms(1:)))

      p_value = 0.0
      do i_tmp = 0, dim/3 - 1
        p_value = p_value - grad(3*i_tmp+1)*vel(3*i_tmp+1) - grad(3*i_tmp+2)*vel(3*i_tmp+2) - grad(3*i_tmp+3)*vel(3*i_tmp+3)
      end do

      if(p_value > 0.0) then
        p_times = p_times + 1
        if(p_times > n_min) then
          dt = min(dt * f_inc , dt_max)
          alpha = alpha * f_alpha
        end if

      else
        p_times = 0
        dt = dt * f_dec
        alpha = alpha_start
        vel = 0.0
      end if

      x(1:dim)=x_new(1:dim)

      mod_force = lalg_nrm2(dim,grad)
      mod_vel = lalg_nrm2(dim,vel)
      
      vel(1:dim) = (1.0 - alpha) * vel(1:dim) - alpha * grad(1:dim) * mod_vel / mod_force
      
      ! Velocity Verlet displacement : vel*dt - 1/2 * grad * dt**2
      vec_delta_pos(1:dim) = vel(1:dim)*dt - grad(1:dim)*dt**2 / mass(1:dim)
      
      x_new(1:dim) = x(1:dim) + vec_delta_pos(1:dim)
      dr_i(1:dim) = sqrt((x_new(1:dim)-x(1:dim))**2)

      do dr_atom_iter = 0, dim/3 - 1
        dr_atoms(dr_atom_iter+1) = sqrt(dr_i(3*dr_atom_iter+1)**2+dr_i(3*dr_atom_iter+2)**2+dr_i(3*dr_atom_iter+3)**2)
      end do

      call write_iter_info(n_iter, dim, en, maxval(dr_atoms(1:)), max_grad_atoms, x)

      if(max_grad_atoms < tolgrad) then
        ierr = 0
        n_iter = maxiter + 1
      else
        n_iter = n_iter + 1
      end if
      
    end do

    SAFE_DEALLOCATE_A(dr_atoms)
    SAFE_DEALLOCATE_A(x_new)
    SAFE_DEALLOCATE_A(dr_i)
    SAFE_DEALLOCATE_A(vec_delta_pos)
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
