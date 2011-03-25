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
!! $Id: opt_control.F90 6950 2010-08-19 22:19:34Z dstrubbe $

  ! ---------------------------------------------------------
  ! The following routines are to be called by C routines, which in turn
  ! are called by the main procedure of this module, opt_control_run, which
  ! is below.
  ! ---------------------------------------------------------


  subroutine opt_control_function_forward(x, f)
    REAL_DOUBLE, intent(in)    :: x
    REAL_DOUBLE, intent(inout) :: f 

    FLOAT, allocatable :: theta(:), y(:)
    integer :: dof
    type(states_t) :: psi

    PUSH_SUB(opt_control_function_forward)

    dof = controlfunction_dof(par_)
    SAFE_ALLOCATE(theta(1:dof))
    SAFE_ALLOCATE(y(1:dof))

    theta = x_
    theta(index_) = x

    call controlfunction_set_theta(par_, theta)
    call controlfunction_theta_to_basis(par_)
    call states_copy(psi, initial_st)
    call propagate_forward(sys_, hm_, td_, par_, target, psi)
    f = - j1_functional(target, sys_%gr, psi, sys_%geo) - controlfunction_j2(par_)
    call states_end(psi)

    SAFE_DEALLOCATE_A(theta)
    SAFE_DEALLOCATE_A(y)
    POP_SUB(opt_control_function_forward)
  end subroutine opt_control_function_forward


  ! ---------------------------------------------------------
  subroutine opt_control_cg_calc(n, x, f, getgrad, df)
    integer,         intent(in)  :: n
    REAL_DOUBLE,     intent(in)  :: x(n)
    REAL_DOUBLE,     intent(inout) :: f
    integer,         intent(in)  :: getgrad
    REAL_DOUBLE,     intent(inout) :: df(n)

    integer :: j
    type(controlfunction_t) :: par_new
    FLOAT :: j1, dx, fmdf
    FLOAT, allocatable :: theta(:), abserr(:), dfn(:), dff(:)
    type(states_t) :: psi


    PUSH_SUB(opt_control_cg_calc)

    SAFE_ALLOCATE(theta(1:n))
    if(getgrad .eq. 1) then
      theta = x
      call controlfunction_set_theta(par_, theta)
      call controlfunction_theta_to_basis(par_)
      call controlfunction_copy(par_new, par_)
      call f_striter(sys_, hm_, td_, par_new, j1)
      f = - j1 - controlfunction_j2(par_)
      if(oct%dump_intermediate) call iterator_write(iterator, par_)
      call iteration_manager_direct(real(-f, REAL_PRECISION), par_, iterator, sys_)
      call controlfunction_set_rep(par_new)
      SAFE_ALLOCATE(dff(1:n))
      dff = df
      call controlfunction_gradient(real(x, REAL_PRECISION), par_, par_new, dff)
      df = dff

      ! Check if the gradient has been computed properly... This should be done only
      ! for debugging purposes.
      if(abs(oct%check_gradient) > M_ZERO) then
        dx = oct%check_gradient
        SAFE_ALLOCATE(dfn(1:n))
        SAFE_ALLOCATE(x_(1:n))
        SAFE_ALLOCATE(abserr(1:n))

        x_ = x
        do j = 1, n
          index_ = j
          call oct_numerical_derivative(x(j), dx, dfn(j), abserr(j), opt_control_function_forward)
        end do

        write(message(1), '(70(''#''))')
        write(message(2), *) &
          'GRADIENT (FORWARD-BACKWARD) |         GRADIENT (NUMERICAL)          |'
        call messages_info(2)
        do j = 1, n
          write(message(1), '(4x,es18.8,7x,a,3x,es18.8,a4,es8.1,6x,a)') &
            df(j), '|', dfn(j), ' +/-', abserr(j), '|'
          call messages_info(1)
        end do
        write(message(1), '(70(''-''))')
        write(message(2), '(a,es18.8,''                                        |'')') 'REL DIFF = ', &
          sqrt(dot_product(df-dfn,df-dfn))/sqrt(dot_product(dfn, dfn))
        write(message(3), '(70(''#''))')
        call messages_info(3)

        SAFE_DEALLOCATE_A(dfn)
        SAFE_DEALLOCATE_A(x_)
        SAFE_DEALLOCATE_A(abserr)
      end if

      call controlfunction_end(par_new)

    else
      theta = x
      call controlfunction_set_theta(par_, theta)
      call controlfunction_theta_to_basis(par_)
      call states_copy(psi, initial_st)
      call propagate_forward(sys_, hm_, td_, par_, target, psi)
      f = - j1_functional(target, sys_%gr, psi, sys_%geo) - controlfunction_j2(par_)
      call states_end(psi)
      if(oct%dump_intermediate) call iterator_write(iterator, par_)
      call iteration_manager_direct(real(-f, REAL_PRECISION), par_, iterator, sys_)
    end if

    SAFE_DEALLOCATE_A(theta)
    POP_SUB(opt_control_cg_calc)
  end subroutine opt_control_cg_calc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine opt_control_cg_messages_info(iter, n, val, maxdx, maxdf, x)
    integer,     intent(in) :: iter, n
    REAL_DOUBLE, intent(in) :: val, maxdx, maxdf
    REAL_DOUBLE, intent(in) :: x(n)

    FLOAT :: fluence, j1, j2, j

    PUSH_SUB(opt_control_cg_messages_info)

    j = - val
    fluence = controlfunction_fluence(par_)
    j2 = controlfunction_j2(par_)
    j1 = j - j2

    write(message(1), '(a,i5)') 'CG optimization iteration #', iter
    call messages_print_stress(stdout, trim(message(1)))

    write(message(1), '(6x,a,f12.5)')    " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)')    " => J        = ", j
    write(message(3), '(6x,a,f12.5)')    " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)')    " => Fluence  = ", fluence
    write(message(5), '(6x,a,f12.5)')    " => Delta    = ", maxdx
    call messages_info(5)
    call messages_print_stress(stdout)

    call iteration_manager_main(iterator, j, j1, j2, real(maxdx, REAL_PRECISION))

    POP_SUB(opt_control_cg_messages_info)
  end subroutine opt_control_cg_messages_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine opt_control_direct_calc(n, x, f)
    integer, intent(in)  :: n
    REAL_DOUBLE, intent(in)  :: x(n)
    REAL_DOUBLE, intent(out) :: f

    FLOAT :: j1, delta
    FLOAT, allocatable :: theta(:)
    type(states_t) :: psi
    type(controlfunction_t) :: par_new

    PUSH_SUB(opt_control_direct_calc)

    SAFE_ALLOCATE(theta(1:n))
    theta = x
    call controlfunction_set_theta(par_, theta)
    call controlfunction_theta_to_basis(par_)

    if(oct%delta == M_ZERO) then
      ! We only need the value of the target functional.
      call states_copy(psi, initial_st)
      call propagate_forward(sys_, hm_, td_, par_, target, psi)
      f = - j1_functional(target, sys_%gr, psi, sys_%geo) - controlfunction_j2(par_)
      if(oct%dump_intermediate) call iterator_write(iterator, par_)
      call iteration_manager_direct(real(-f, REAL_PRECISION), par_, iterator, sys_)
      call states_end(psi)
    else
      call controlfunction_copy(par_new, par_)
      call f_striter(sys_, hm_, td_, par_new, j1)
      delta = controlfunction_diff(par_, par_new)
      f = - oct%eta * j1 + oct%delta * delta
      if(oct%dump_intermediate) call iterator_write(iterator, par_)
      call iteration_manager_direct(real(-f, REAL_PRECISION), par_, iterator, sys_, delta)
      call controlfunction_end(par_new)
    end if

    SAFE_DEALLOCATE_A(theta)
    POP_SUB(opt_control_direct_calc)
  end subroutine opt_control_direct_calc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine opt_control_direct_messages_info(iter, n, val, maxdx, x)
    integer,     intent(in) :: iter, n
    REAL_DOUBLE, intent(in) :: val, maxdx
    REAL_DOUBLE, intent(in) :: x(n)

    FLOAT :: fluence, j1, j2, j
    FLOAT, allocatable :: theta(:)

    PUSH_SUB(opt_control_direct_messages_info)

    SAFE_ALLOCATE(theta(1:n))
    theta = x
    call controlfunction_set_theta(par_, theta)
    call controlfunction_theta_to_basis(par_)
    SAFE_DEALLOCATE_A(theta)

    j = - val
    fluence = controlfunction_fluence(par_)
    j2 = controlfunction_j2(par_)
    j1 = j - j2

    write(message(1), '(a,i5)') 'Direct optimization iteration #', iter
    call messages_print_stress(stdout, trim(message(1)))

    write(message(1), '(6x,a,f12.5)')    " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)')    " => J        = ", j
    write(message(3), '(6x,a,f12.5)')    " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)')    " => Fluence  = ", fluence
    write(message(5), '(6x,a,f12.5)')    " => Delta    = ", maxdx
    call messages_info(5)
    call messages_print_stress(stdout)

    call iteration_manager_main(iterator, j, j1, j2, real(maxdx, REAL_PRECISION))

    POP_SUB(opt_control_direct_messages_info)
  end subroutine opt_control_direct_messages_info
  ! ---------------------------------------------------------


