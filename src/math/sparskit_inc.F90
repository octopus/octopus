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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

! ---------------------------------------------------------
subroutine X(sparskit_solver_run)(sk, op, opt, sol, rhs)
  type(sparskit_solver_t), intent(inout) :: sk
  R_TYPE, intent(in)  :: rhs(:)
  R_TYPE, intent(out) :: sol(:)

#ifdef R_TREAL
  interface
    subroutine op(x, y)
      implicit none
      FLOAT, intent(in)  :: x(:)
      FLOAT, intent(out) :: y(:)
    end subroutine op
    subroutine opt(x, y)
      implicit none
      FLOAT, intent(in)  :: x(:)
      FLOAT, intent(out) :: y(:)
    end subroutine opt
  end interface
#endif
#ifdef R_TCOMPLEX
  interface
    subroutine op(xre, xim, yre, yim)
      implicit none
      FLOAT, intent(in)  :: xre(:), xim(:)
      FLOAT, intent(out) :: yre(:), yim(:)
    end subroutine op
    subroutine opt(xre, xim, yre, yim)
      implicit none
      FLOAT, intent(in)  :: xre(:), xim(:)
      FLOAT, intent(out) :: yre(:), yim(:)
    end subroutine opt
  end interface
#endif

  integer :: iter

  PUSH_SUB(X(sparskit_solver_run))

  ! initialize counter
  sk%used_iter = 0

#ifdef R_TREAL
  ASSERT(.not. sk%is_complex)
  sk%sk_b = rhs
  ! initial guess
  sk%sk_y = sol
#endif
#ifdef R_TCOMPLEX
  ASSERT(sk%is_complex)
  do iter = 1, sk%size/2
    sk%sk_b(iter)           = real (rhs(iter))
    sk%sk_b(iter + sk%size/2) = aimag(rhs(iter))
    ! initial guess
    sk%sk_y(iter)           = real (sol(iter))
    sk%sk_y(iter + sk%size/2) = aimag(sol(iter))
  end do
#endif

  ! Start iterative solution of the linear system
  solver_iter: do iter = 1, sk%maxiter
    
    select case(sk%solver_type)
    case(SK_CG)
      call cg(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_CGNR)                          
      call cgnr(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_BCG)                           
      call bcg(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_DBCG)                          
      call dbcg(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_BCGSTAB)                       
      call bcgstab(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_TFQMR)                         
      call tfqmr(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_FOM)                           
      call fom(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_GMRES)                         
      call gmres(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_FGMRES)                        
      call fgmres(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case(SK_DQGMRES)                       
      call dqgmres(sk%size, sk%sk_b, sk%sk_y, sk%ipar, sk%fpar, sk%sk_work)
    case default
      write(message(1), '(a,i4,a)') "Input: '", sk%solver_type, &
           "' is not a valid SPARSKIT solver."
      message(2) = '( SPARSKITSolver =  cg | cgnr | bcg | dbcg | bcgstab | tfqmr | fom | gmres | fgmres | dqgmres )'
      call messages_fatal(2)
    end select
    
    ! Evaluate reverse communication protocol
    select case(sk%ipar(1))
    case(1)
#ifdef R_TREAL
      call op(sk%sk_work(sk%ipar(8):sk%ipar(8)+sk%size),sk%sk_work(sk%ipar(9):sk%ipar(9)+sk%size))
#endif
#ifdef R_TCOMPLEX
      call op(sk%sk_work(sk%ipar(8):sk%ipar(8)+sk%size/2),sk%sk_work(sk%ipar(8)+sk%size/2:sk%ipar(8)+sk%size), &
           sk%sk_work(sk%ipar(9):sk%ipar(9)+sk%size/2),sk%sk_work(sk%ipar(9)+sk%size/2:sk%ipar(9)+sk%size))
#endif
    case(2)
      ! call atmux(n,w(sk%ipar(8)),w(sk%ipar(9)),a,ja,ia)
#ifdef R_TREAL
      call opt(sk%sk_work(sk%ipar(8):sk%ipar(8)+sk%size),sk%sk_work(sk%ipar(9):sk%ipar(9)+sk%size))
#endif
#ifdef R_TCOMPLEX
      call opt(sk%sk_work(sk%ipar(8):sk%ipar(8)+sk%size/2),sk%sk_work(sk%ipar(8)+sk%size/2:sk%ipar(8)+sk%size), &
           sk%sk_work(sk%ipar(9):sk%ipar(9)+sk%size/2),sk%sk_work(sk%ipar(9)+sk%size/2:sk%ipar(9)+sk%size))
#endif
    case(3, 4, 5, 6)
      ! left preconditioner solver
      ! left preconditioner transposed solve
      ! right preconditioner solve
      ! right preconditioner transposed solve
      call messages_not_implemented('Sparskit preconditioning')
    case(0)
      ! successful exit of solver
      exit solver_iter
    case(-1)
!      message(1) = 'Maximum iteration number "SPARSKITMaxIter" exceeded.'
!      call messages_warning(1)
      exit solver_iter
    case(-2)
      message(1) = 'Insufficient work space.'
      call messages_fatal(1)
    case(-3)
      message(1) = 'Anticipated break-down / divide by zero.'
      call messages_fatal(1)
    case(-4)
      message(1) = '"SPARSKITRelTolerance" and "SPARSKITAbsTolerance" are'
      message(2) = 'both <= 0. Valid ranges are 0 <= SPARSKITRelTolerance < 1,'
      message(3) = '0 <= SPARSKITAbsTolerance.'
      call messages_fatal(3)
    case(-9)
      message(1) = 'While trying to detect a break-down, an abnormal number is detected.'
      call messages_fatal(1)
    case(-10)
      message(1) = 'Return due to some non-numerical reasons, e.g. invalid floating-point numbers etc.'
      call messages_fatal(1)
    case default
      message(1) = 'Unknown SPARSKIT return value. Exiting ...'
      call messages_fatal(1)
    end select

    if(sk%iter_out > 0) then
      if(mod(iter, sk%iter_out) == 0) then
        write(message(1), '(a,i7)') 'SPARSKIT Iter: ', iter
        call messages_info(1)
      end if
    end if
      
  end do solver_iter



  if(iter  > sk%maxiter) then
!    message(1) = 'Maxiter reached'
!    call messages_warning(1)
  end if

  ! set back to zero to initialize the solver for the next call
  sk%ipar(1) = 0
  ! store the number of iterations used
  sk%used_iter = iter - 1
  ! reset 
  sk%ipar(7) = 0

  ! store current error norm
  sk%residual_norm = sk%fpar(6)

  ! output status info
  if(sk%verbose) then
    write(message(1), '(a,I5,a,E19.12)') 'SPARSKIT iter: ', sk%used_iter, ' residual norm: ', sk%residual_norm
    call messages_info(1)
  end if

#ifdef R_TREAL
  sol = sk%sk_y
#endif
#ifdef R_TCOMPLEX
  do iter = 1, sk%size/2
    sol(iter) = sk%sk_y(iter) + M_zI*sk%sk_y(iter+sk%size/2)
  end do
#endif

  POP_SUB(X(sparskit_solver_run))

end subroutine X(sparskit_solver_run)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
