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


subroutine DRIVER (sk, rhs, sol, work, op) 
  type(sparskit_solver_type), intent(inout) :: sk
  FLOAT, intent(in)    :: rhs(:)
  FLOAT, intent(out)   :: sol(:)
  FLOAT, intent(inout) :: work(:)
  interface
     subroutine op(x, y)
       FLOAT, intent(in)  :: x(:)
       FLOAT, intent(out) :: y(:)
     end subroutine op
  end interface

  integer :: i

  call push_sub("SPARSKIT_DRIVER")

  solver_iter: do i = 1, sk%maxiter

     ! Run actual solver
#if defined(HAVE_SPARSKIT)
     call SOLVER(sk%size, rhs, sol, sk%ipar, sk%fpar, work)
#else
     message(1) = 'Error: Sparskit library required for usage of SparskitSolver'
     message(2) = '       configure with --with-sparskit=DIR and remake.'
     call write_fatal(2)  
#endif
     ! Evaluate reverse communication protocol
     select case(sk%ipar(1))
     case(1)
        call op(work(sk%ipar(8):sk%ipar(8)+sk%size),work(sk%ipar(9):sk%ipar(9)+sk%size))
     case(2)
        ! call atmux(n,w(sk%ipar(8)),w(sk%ipar(9)),a,ja,ia)
        message(1) = 'Error: Matrix vector multiplication with A^T not implemented yet.'
        call write_fatal(1)        
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
        message(1) = 'Error: Maximum iteration number "SparskitMaxIter" exceeded.'
        call write_fatal(1)
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

  enddo solver_iter

  if(i.gt.sk%maxiter) then
     message(1) = 'Warning: Maxiter reached'
     call write_warning(1)
  endif

  ! store the number of iterations used
  sk%used_iter = i

  call pop_sub()
end subroutine DRIVER
