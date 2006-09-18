!! Copyright (C) 2005-2006 Heiko Appel
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ---------------------------------------------------------
subroutine DRIVER (sk, op, opt, sol, rhs, sk_work)
  type(sparskit_solver_t), intent(inout) :: sk
  R_TYPE, intent(in)    :: rhs(:)
  R_TYPE, intent(out)   :: sol(:)
  FLOAT,  intent(inout) :: sk_work(:)
#ifdef R_TREAL
  interface
    subroutine op(x, y)
      FLOAT, intent(in)  :: x(:)
      FLOAT, intent(out) :: y(:)
    end subroutine op
    subroutine opt(x, y)
      FLOAT, intent(in)  :: x(:)
      FLOAT, intent(out) :: y(:)
    end subroutine opt
  end interface
#endif
#ifdef R_TCOMPLEX
  interface
    subroutine op(xre, xim, yre, yim)
      FLOAT, intent(in)  :: xre(:), xim(:)
      FLOAT, intent(out) :: yre(:), yim(:)
    end subroutine op
    subroutine opt(xre, xim, yre, yim)
      FLOAT, intent(in)  :: xre(:), xim(:)
      FLOAT, intent(out) :: yre(:), yim(:)
    end subroutine opt
  end interface
#endif

  integer :: i

  call push_sub('sparskit_driver.sparskit_driver')
#ifdef R_TREAL
  sk_b = rhs
  ! initial guess
  sk_y = sol
#endif
#ifdef R_TCOMPLEX
  do i = 1, sk%size/2
    sk_b(i)           = real (rhs(i))
    sk_b(i+sk%size/2) = aimag(rhs(i))
    ! initial guess
    sk_y(i)           = real (sol(i))
    sk_y(i+sk%size/2) = aimag(sol(i))
  end do
#endif


  ! Start iterative solution of the linear system
  solver_iter: do i = 1, sk%maxiter
    ! Run actual solver
#ifdef HAVE_SPARSKIT
    if(in_debug_mode) call push_sub('sparskit_driver.sparskit_solver')

    call SOLVER(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)

    if(in_debug_mode) call pop_sub()
#else
    message(1) = 'Error: Sparskit library required for usage of SparskitSolver'
    message(2) = '       configure with --with-sparskit=DIR and remake.'
    call write_fatal(2)
#endif

    ! Evaluate reverse communication protocol
    select case(sk%ipar(1))
    case(1)
#ifdef R_TREAL
      call op(sk_work(sk%ipar(8):sk%ipar(8)+sk%size),sk_work(sk%ipar(9):sk%ipar(9)+sk%size))
#endif
#ifdef R_TCOMPLEX
      call op(sk_work(sk%ipar(8):sk%ipar(8)+sk%size/2),sk_work(sk%ipar(8)+sk%size/2:sk%ipar(8)+sk%size), &
        sk_work(sk%ipar(9):sk%ipar(9)+sk%size/2),sk_work(sk%ipar(9)+sk%size/2:sk%ipar(9)+sk%size))
#endif
    case(2)
      ! call atmux(n,w(sk%ipar(8)),w(sk%ipar(9)),a,ja,ia)
#ifdef R_TREAL
      call opt(sk_work(sk%ipar(8):sk%ipar(8)+sk%size),sk_work(sk%ipar(9):sk%ipar(9)+sk%size))
#endif
#ifdef R_TCOMPLEX
      call opt(sk_work(sk%ipar(8):sk%ipar(8)+sk%size/2),sk_work(sk%ipar(8)+sk%size/2:sk%ipar(8)+sk%size), &
        sk_work(sk%ipar(9):sk%ipar(9)+sk%size/2),sk_work(sk%ipar(9)+sk%size/2:sk%ipar(9)+sk%size))
#endif
    case(3)
      ! left preconditioner solver
      message(1) = 'Error: Preconditioning not implemented yet.'
      call write_fatal(1)
    case(4)
      ! left preconditioner transposed solve
      message(1) = 'Error: Preconditioning not implemented yet.'
      call write_fatal(1)
    case(5)
      ! right preconditioner solve
      message(1) = 'Error: Preconditioning not implemented yet.'
      call write_fatal(1)
    case(6)
      ! right preconditioner transposed solve
      message(1) = 'Error: Preconditioning not implemented yet.'
      call write_fatal(1)
    case(0)
      ! successful exit of solver
      exit solver_iter
    case(-1)
!      message(1) = 'Warning: Maximum iteration number "SparskitMaxIter" exceeded.'
!      call write_warning(1)
      exit solver_iter
    case(-2)
      message(1) = 'Error: Insufficient work space.'
      call write_fatal(1)
    case(-3)
      message(1) = 'Error: Anticipated break-down / divide by zero.'
      call write_fatal(1)
    case(-4)
      message(1) = 'Error: "SparskitRelTolerance" and "SparskitAbsTolerance" are'
      message(2) = '       both <= 0. Valid ranges are 0 <= SparskitRelTolerance < 1,'
      message(3) = '       0 <= SparskitAbsTolerance.'
      call write_fatal(3)
    case(-9)
      message(1) = 'Error: while trying to detect a break-down, an abnormal number is detected.'
      call write_fatal(1)
    case(-10)
      message(1) = 'Error: return due to some non-numerical reasons, e.g. invalid'
      message(2) = 'floating-point numbers etc.'
      call write_fatal(2)
    case default
      message(1) = 'Error: Unknown Sparskit return value. Exiting ...'
      call write_fatal(1)
    end select

    if(sk%iter_out > 0) then
      if(mod(i, sk%iter_out) == 0) then
        write(message(1), '(a,i7)') 'Sparskit Iter: ', i
        call write_info(1)
      end if
    end if
      
  end do solver_iter

  if(i.gt.sk%maxiter) then
!    message(1) = 'Warning: Maxiter reached'
!    call write_warning(1)
  end if

  ! set back to zero to initialize the solver for the next call
  sk%ipar(1) = 0
  ! store the number of iterations used
  sk%used_iter = i - 1
  ! reset 
  sk%ipar(7) = 0

  ! store current error norm
  sk%residual_norm = sk%fpar(6)

  ! output status info
  if(sk%verbose) then
    write(message(1), '(a,I5,a,E18.12)') 'Sparskit iter: ', sk%used_iter, ' residual norm: ', sk%residual_norm
    call write_info(1)
  end if

#ifdef R_TREAL
  sol = sk_y
#endif
#ifdef R_TCOMPLEX
  do i = 1, sk%size/2
    sol(i) = sk_y(i) + M_zI*sk_y(i+sk%size/2)
  end do
#endif

  call pop_sub()
end subroutine DRIVER
