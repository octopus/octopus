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
!! $Id$

! ---------------------------------------------------------
subroutine X(sparskit_solver_init)(n, sk)
  type(sparskit_solver_t), intent(out) :: sk
  integer, intent(in)  :: n

  integer :: workspace_size, m

  PUSH_SUB(X(sparskit_solver_init))

  !%Variable SPARSKITSolver
  !%Type integer
  !%Default sk_cg
  !%Section Math::General
  !%Description
  !% Specifies what kind of linear solver will be used.
  !%Option sk_cg 1
  !% Conjugate Gradient Method
  !%Option sk_cgnr 2
  !% Conjugate Gradient Method (Normal Residual equation)
  !%Option sk_bcg 3
  !% Bi-Conjugate Gradient Method
  !%Option sk_dbcg 4
  !% BCG with partial pivoting
  !%Option sk_bcgstab 5
  !% BCG stabilized
  !%Option sk_tfqmr 6
  !% Transpose-Free Quasi-Minimum Residual method
  !%Option sk_fom 7
  !% Full Orthogonalization Method
  !%Option sk_gmres 8
  !% Generalized Minimum Residual method
  !%Option sk_fgmres 9
  !% Flexible version of Generalized Minimum Residual method
  !%Option sk_dqgmres 10
  !% Direct versions of the Quasi-Generalized Minimum Residual method
  !%End
  call parse_integer(datasets_check('SPARSKITSolver'),          SK_CG, sk%solver_type)
  if ( sk%solver_type.lt.SK_MINVAL.or.sk%solver_type.gt.SK_MAXVAL ) then
    call input_error('SPARSKITSolver')
  end if

  !%Variable SPARSKITKrylovSubspaceSize
  !%Type integer
  !%Default 15
  !%Section Math::General
  !%Description
  !% Some of the SPARSKIT solvers are Krylov subspace methods.
  !% This variable determines what size the solver will use 
  !% for the subspace.
  !%End
  call parse_integer(datasets_check('SPARSKITKrylovSubspaceSize'), 15, sk%krylov_size)

  ! preconditioner not implemented
  sk%preconditioning = 0

  !%Variable SPARSKITMaxIter
  !%Type integer
  !%Default 0
  !%Section Math::General
  !%Description
  !% This variable controls the maximum number of iteration steps that
  !% will be performed by the (iterative) linear solver.
  !%End
  call parse_integer(datasets_check('SPARSKITMaxIter'),          5000, sk%maxiter)

  !%Variable SPARSKITIterOut
  !%Type integer
  !%Default 0
  !%Section Math::General
  !%Description
  !% Determines how often status info of the solver is printed.
  !%End
  call parse_integer(datasets_check('SPARSKITIterOut'),            -1, sk%iter_out)

  !%Variable SPARSKITRelTolerance
  !%Type float
  !%Default 1e-8
  !%Section Math::General
  !%Description
  !% Some SPARSKIT solvers use a relative tolerance as a stopping criterion 
  !% for the iterative solution process. This variable can be used to 
  !% specify the tolerance.
  !%End
  call parse_float(datasets_check('SPARSKITRelTolerance'), CNST(1e-5), sk%rel_tolerance)

  !%Variable SPARSKITAbsTolerance
  !%Type float
  !%Default 1e-8
  !%Section Math::General
  !%Description
  !% Some SPARSKIT solvers use an absolute tolerance as a stopping criterion 
  !% for the iterative solution process. This variable can be used to 
  !% specify the tolerance.
  !%End
  call parse_float(datasets_check('SPARSKITAbsTolerance'), CNST(1e-10), sk%abs_tolerance)

  !%Variable SPARSKITVerboseSolver
  !%Type logical
  !%Default no
  !%Section Math::General
  !%Description
  !% When set to yes, the SPARSKIT solver will write more detailed output.
  !%End
  call parse_logical(datasets_check('SPARSKITVerboseSolver'), .false., sk%verbose)

  ! size of the problem
  sk%size = n
#ifdef R_TCOMPLEX
  sk%size = 2*n
#endif

  ! initialize workspace size
  workspace_size = 0 

  ! Krylov subspace size
  m = sk%krylov_size
  if (mod(m, 2).ne.0) m = m + 1

  select case(sk%solver_type)
  case(SK_CG)
    message(1) = 'Info: SPARSKIT solver type: Conjugate Gradient Method'
    workspace_size = 5*sk%size
  case(SK_CGNR)
    message(1) = 'Info: SPARSKIT solver type: Conjugate Gradient Method (Normal Residual equation)'
    workspace_size = 5*sk%size
  case(SK_BCG)
    message(1) = 'Info: SPARSKIT solver type: Bi-Conjugate Gradient Method'
    workspace_size = 7*sk%size
  case(SK_DBCG)
    message(1) = 'Info: SPARSKIT solver type: BCG with partial pivoting'
    workspace_size = 11*sk%size
  case(SK_BCGSTAB)
    message(1) = 'Info: SPARSKIT solver type: BCG stabilized'
    workspace_size = 8*sk%size
  case(SK_TFQMR)
    message(1) = 'Info: SPARSKIT solver type: Transpose-Free Quasi-Minimum Residual method'
    workspace_size = 11*sk%size
  case(SK_FOM)
    message(1) = 'Info: SPARSKIT solver type: Full Orthogonalization Method'
    workspace_size = (sk%size+3)*(m+2) + (m+1)*m/2
  case(SK_GMRES)
    message(1) = 'Info: SPARSKIT solver type: Generalized Minimum Residual method'
    workspace_size = (sk%size+3)*(m+2) + (m+1)*m/2
  case(SK_FGMRES)
    message(1) = 'Info: SPARSKIT solver type: Flexible version of Generalized Minimum Residual method'
    workspace_size =  2*sk%size*(m+1) + (m+1)*m/2 + 3*m + 2
  case(SK_DQGMRES)
    message(1) = 'Info: SPARSKIT solver type: Direct versions of Quasi-Generalized Minimum Residual method'
    workspace_size = sk%size + (m+1) * (2*sk%size+4)
  case default
    write(message(1), '(a,i4,a)') "Input: '", sk%solver_type, &
      "' is not a valid SPARSKIT Solver"
    message(2) = '( SPARSKIT Solver =  cg | cgnr | bcg | dbcg | bcgstab | tfqmr | fom | gmres | fgmres | dqgmres )'
    call messages_fatal(2)
  end select
  call messages_info(1)

  ! Now we initialize the arrays for the reverse communication protocol
  sk%ipar = 0
  sk%fpar = 0

  ! A call to the solver with ipar(1) == 0 will initialize the iterative solver.
  sk%ipar(1) = 0

  ! Stopping criteria; use convergence test scheme 1
  sk%ipar(3) = 1

  sk%ipar(4) = workspace_size
  sk%ipar(5) = sk%krylov_size

  ! Maximum number of matrix-vector multiplies
  sk%ipar(6) = sk%maxiter

  ! Relative tolerance
  sk%fpar(1) = sk%rel_tolerance

  ! Absolute tolerance
  sk%fpar(2) = sk%abs_tolerance

  ! allocate and initialize work arrays
  SAFE_ALLOCATE(sk_b(1:sk%size))
  SAFE_ALLOCATE(sk_y(1:sk%size))
  SAFE_ALLOCATE(sk_work(1:workspace_size))
  sk_work = M_ZERO
  sk_y    = M_ZERO

  POP_SUB(X(sparskit_solver_init))
end subroutine X(sparskit_solver_init)


! ---------------------------------------------------------
subroutine X(sparskit_solver_run)(sk, op, opt, sol, rhs)
  type(sparskit_solver_t), intent(inout) :: sk
  R_TYPE, intent(in)  :: rhs(:)
  R_TYPE, intent(out) :: sol(:)

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

  integer :: iter

  PUSH_SUB(X(sparskit_solver_run))

  ! initialize counter
  sk%used_iter = 0

#ifdef R_TREAL
  sk_b = rhs
  ! initial guess
  sk_y = sol
#endif
#ifdef R_TCOMPLEX
  do iter = 1, sk%size/2
    sk_b(iter)           = real (rhs(iter))
    sk_b(iter + sk%size/2) = aimag(rhs(iter))
    ! initial guess
    sk_y(iter)           = real (sol(iter))
    sk_y(iter + sk%size/2) = aimag(sol(iter))
  end do
#endif

  ! Start iterative solution of the linear system
  solver_iter: do iter = 1, sk%maxiter
    
    select case(sk%solver_type)
    case(SK_CG)
      call cg(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_CGNR)                          
      call cgnr(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_BCG)                           
      call bcg(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_DBCG)                          
      call dbcg(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_BCGSTAB)                       
      call bcgstab(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_TFQMR)                         
      call tfqmr(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_FOM)                           
      call fom(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_GMRES)                         
      call gmres(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_FGMRES)                        
      call fgmres(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
    case(SK_DQGMRES)                       
      call dqgmres(sk%size, sk_b, sk_y, sk%ipar, sk%fpar, sk_work)
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
!      message(1) = 'Warning: Maximum iteration number "SPARSKITMaxIter" exceeded.'
!      call messages_warning(1)
      exit solver_iter
    case(-2)
      message(1) = 'Error: Insufficient work space.'
      call messages_fatal(1)
    case(-3)
      message(1) = 'Error: Anticipated break-down / divide by zero.'
      call messages_fatal(1)
    case(-4)
      message(1) = 'Error: "SPARSKITRelTolerance" and "SPARSKITAbsTolerance" are'
      message(2) = '       both <= 0. Valid ranges are 0 <= SPARSKITRelTolerance < 1,'
      message(3) = '       0 <= SPARSKITAbsTolerance.'
      call messages_fatal(3)
    case(-9)
      message(1) = 'Error: while trying to detect a break-down, an abnormal number is detected.'
      call messages_fatal(1)
    case(-10)
      message(1) = 'Error: return due to some non-numerical reasons, e.g. invalid'
      message(2) = 'floating-point numbers etc.'
      call messages_fatal(2)
    case default
      message(1) = 'Error: Unknown SPARSKIT return value. Exiting ...'
      call messages_fatal(1)
    end select

    if(sk%iter_out > 0) then
      if(mod(iter, sk%iter_out) == 0) then
        write(message(1), '(a,i7)') 'SPARSKIT Iter: ', iter
        call messages_info(1)
      end if
    end if
      
  end do solver_iter



  if(iter .gt.sk%maxiter) then
!    message(1) = 'Warning: Maxiter reached'
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
    write(message(1), '(a,I5,a,E18.12)') 'SPARSKIT iter: ', sk%used_iter, ' residual norm: ', sk%residual_norm
    call messages_info(1)
  end if

#ifdef R_TREAL
  sol = sk_y
#endif
#ifdef R_TCOMPLEX
  do iter = 1, sk%size/2
    sol(iter) = sk_y(iter) + M_zI*sk_y(iter+sk%size/2)
  end do
#endif

  POP_SUB(X(sparskit_solver_run))

end subroutine X(sparskit_solver_run)


! ---------------------------------------------------------
subroutine X(sparskit_solver_end)()
  PUSH_SUB(X(sparskit_solver_end))

  SAFE_DEALLOCATE_A(sk_b)
  SAFE_DEALLOCATE_A(sk_y)
  SAFE_DEALLOCATE_A(sk_work)

  POP_SUB(X(sparskit_solver_end))
end subroutine X(sparskit_solver_end)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
