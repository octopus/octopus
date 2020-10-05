!! Copyright (C) 2004-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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


!> The two following subroutines, sym_conjugate_gradients, and bi_conjugate_gradients,
!! must be called under a common interface: conjugate_gradients. It provides an
!! approximate solution to the linear system problem  Ax = b.
!! Solving a symmetric linear system, which is either real or complex symmetric or
!! Hermitian the best choice is sym_conjugate_gradients (does not need A^\dagger)
!! Solving a real unsymmetric or a complex non-Hermitian system bi_conjugate_gradients
!! has to be chosen (where one does need A^\dagger).
!!
!! Note: the complex-valued versions (both CG and BiCG) work only with
!!       symmetric operators but they may be non-Hermitian. This is a
!!       property of the CG algorithm. A comment on this can be found
!!       in chapter 4.1 of ftp://ftp.netlib.org/templates/templates.ps
!!
!! subroutine conjugate_gradients(np, x, b, op, [opt,] dotp, iter [, residue] [, threshold] )
!!    integer, intent(in)  :: np      => The dimension of the problem.
!!    FLOAT, intent(inout) :: x       => On input, an estimate to the solution.
!!                                    => On output, the approximate solution.
!!    FLOAT, intent(in)    :: b       => The inhomogeneous term of the equation.
!!    interface
!!      subroutine op(x, y)
!!         FLOAT, intent(in)  :: x(:)
!!         FLOAT, intent(out) :: y(:)
!!      end subroutine op
!!    end interface                   => This should be a procedure that
!!                                       computes Ax = y.
!!    interface
!!      subroutine opt(x, y)
!!         FLOAT, intent(in)  :: x(:)
!!         FLOAT, intent(out) :: y(:)
!!      end subroutine opt
!!    end interface                   => If present, this should be a procedure that
!!                                       computes A^\dagger x = y.
!!                                       Only useful for non-Hermitian operators.
!!    interface
!!      R_TYPE function dotp(x, y)
!!      R_TYPE, intent(inout) :: x(:)
!!      R_TYPE, intent(in)    :: y(:)
!!    end function dotp               => Calculates the <x | y>.
!!                                       Depending on the matrix A one should choose:
!!                                       complex symmetric: <x | y> = x^T * y
!!                                       hermitian:         <x | y> = x^\dagger * y
!!                                       general:           <x | y> = x^\dagger * y
!!    integer, intent(inout) :: iter  => On input, the maximum number of iterations that
!!                                       the procedure is allowed to take.
!!                                       On output, the iterations it actually did.
!!    FLOAT, intent(out) :: residue   => If present, it measures the final error:
!!                                       residue = <Ax - b | Ax - b>.
!!    FLOAT, intent(in)  :: threshold => If present, it sets the required accuracy
!!                                       threshold for the algorithm to stop. If not
!!                                       present, this is set to 1.0e-6. [The algorithm
!!                                       stops when <Ax - b | Ax - b> <= threshold, or
!!                                       iter iterations are reached.]
!! end subroutine conjugate_gradients
!!
!! (*) NOTE: The algorithm assumes that the vectors are given in an orthonormal basis.


! ---------------------------------------------------------
subroutine X(sym_conjugate_gradients)(np, x, b, op, dotp, iter, residue, threshold)
  integer, intent(in)    :: np
  R_TYPE,  intent(inout) :: x(:)
  R_TYPE,  intent(in)    :: b(:)
  interface
    subroutine op(x, y)
      implicit none
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine op
    R_TYPE function dotp(x, y) result(res)
      implicit none
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(in)    :: y(:)
    end function dotp
  end interface
  integer,          intent(inout) :: iter
  FLOAT,  optional, intent(in)    :: threshold
  R_TYPE, optional, intent(out)   :: residue

  R_TYPE, allocatable :: r(:), ax(:), p(:), ap(:)
  R_TYPE              :: alpha, beta, gamma
  FLOAT               :: threshold_
  integer             :: max_iter, ip

  PUSH_SUB(X(sym_conjugate_gradients))

  threshold_ = optional_default(threshold, CNST(1.0e-6))

  SAFE_ALLOCATE( r(1:np))
  SAFE_ALLOCATE(ax(1:np))
  SAFE_ALLOCATE( p(1:ubound(x, 1)))
  SAFE_ALLOCATE(ap(1:np))

  ! Initial residue.
  call op(x, ax)
  do ip = 1, np
    r(ip) = b(ip) - ax(ip)
  end do

  ! Initial search direction.
  call lalg_copy(np, r, p)

  max_iter = iter
  iter = 1
  do while(iter < max_iter)
    gamma = dotp(r, r)
    if(abs(gamma) < threshold_) exit
    call op(p, ap)
    alpha   = gamma/dotp(p, ap)
    call lalg_axpy(np, -alpha, ap, r)
    call lalg_axpy(np, alpha, p, x)
    beta    = dotp(r, r)/gamma
    do ip = 1, np
      p(ip) = r(ip) + beta*p(ip)
    end do
    iter    = iter + 1
  end do
  if(present(residue)) residue = gamma

  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(ax)
  SAFE_DEALLOCATE_A(p)
  SAFE_DEALLOCATE_A(ap)

  POP_SUB(X(sym_conjugate_gradients))
end subroutine X(sym_conjugate_gradients)


! ---------------------------------------------------------
subroutine X(bi_conjugate_gradients)(np, x, b, op, opt, dotp, iter, residue, threshold)
  integer, intent(in)    :: np
  R_TYPE,  intent(inout) :: x(:)
  R_TYPE,  intent(in)    :: b(:)
  interface
    subroutine op(x, y)
      implicit none
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine op
  end interface
  interface
    subroutine opt(x, y)
      implicit none
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine opt
    R_TYPE function dotp(x, y) result(res)
      implicit none
      R_TYPE, intent(in) :: x(:)
      R_TYPE, intent(in) :: y(:)
    end function dotp
  end interface
  integer,          intent(inout) :: iter
  R_TYPE, optional, intent(out)   :: residue
  FLOAT, optional,  intent(in)    :: threshold

  R_TYPE, allocatable :: r(:), rr(:), ax(:), p(:), pp(:), ap(:), atp(:)
  FLOAT               :: alpha, beta, gamma, err, threshold_
  integer             :: max_iter

  PUSH_SUB(X(bi_conjugate_gradients))

  threshold_ = optional_default(threshold, CNST(1.0e-6))

  SAFE_ALLOCATE(  r(1:np))
  SAFE_ALLOCATE( rr(1:np))
  SAFE_ALLOCATE( ax(1:np))
  SAFE_ALLOCATE(  p(1:ubound(x, 1)))
  SAFE_ALLOCATE( pp(1:ubound(x, 1)))
  SAFE_ALLOCATE( ap(1:np))
  SAFE_ALLOCATE(atp(1:np))

  ! Initial residue.
  call op(x, ax)
  ! r <- b - ax, rr <- r
  call lalg_copy(np, b, r)
  call lalg_axpy(np, -R_TOTYPE(M_ONE), ax, r)
  call lalg_copy(np, r, rr)

  ! Initial search direction.
  call lalg_copy(np, r, p)
  call lalg_copy(np, p, pp)

  max_iter = iter
  iter     = 1
  do while(iter < max_iter)
    gamma = R_REAL(dotp(rr, r))
    err   = dotp(r, r)
    if(abs(err) < threshold_) exit
    call op (p,  ap)
    call opt(pp, atp)
    alpha = gamma/R_REAL(dotp(pp, ap))
    call lalg_axpy(np, -alpha, ap, r)         ! r  <- r - alpha*ap
    call lalg_axpy(np, -alpha, atp, rr)       ! rr <- rr - alpha*atp
    call lalg_axpy(np, alpha, p, x)           ! x  <- x + alpha*p
    beta = R_REAL(dotp(rr, r))/gamma
    call lalg_scal(np, beta, p)               ! p  <- r + beta*p
    call lalg_axpy(np, R_TOTYPE(M_ONE), r, p)
    call lalg_scal(np, beta, pp)              ! pp <- rr + beta*pp
    call lalg_axpy(np, R_TOTYPE(M_ONE), rr, pp)
    iter = iter + 1
  end do
  if(present(residue)) residue = err

  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(rr)
  SAFE_DEALLOCATE_A(ax)
  SAFE_DEALLOCATE_A(p)
  SAFE_DEALLOCATE_A(pp)
  SAFE_DEALLOCATE_A(ap)
  SAFE_DEALLOCATE_A(atp)

  POP_SUB(X(bi_conjugate_gradients))
end subroutine X(bi_conjugate_gradients)

  ! ---------------------------------------------------------
  subroutine X(qmr_sym_spec_dotu)(np, x, b, op, prec, iter, residue, threshold, showprogress, converged)
    integer, target, intent(in)    :: np    !< number of points
    R_TYPE,          intent(inout) :: x(:)  !< initial guess and result
    R_TYPE,          intent(in)    :: b(:)  !< the right side
    interface
      subroutine op(x, y)                   !< the matrix A as operator, y <- Ax
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine prec(x, y)                 !< preconditioner
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prec
    end interface
    integer,           intent(inout) :: iter         !< (in) number of max iterations, (out) used iterations
    FLOAT, optional,   intent(out)   :: residue      !< residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold    !< convergence threshold
    logical, optional, intent(in)    :: showprogress !< show progress bar
    logical, optional, intent(out)   :: converged    !< is the algorithm converged

    PUSH_SUB(X(qmr_sym_spec_dotu))

    np_p => np
    call X(qmr_sym_gen_dotu)(np, x, b, op, X(dotu_qmr), X(nrm2_qmr), prec, iter, &
      residue, threshold, showprogress, converged)

    POP_SUB(X(qmr_sym_spec_dotu))
  end subroutine X(qmr_sym_spec_dotu)

  ! ---------------------------------------------------------
  subroutine X(qmr_spec_dotu)(np, x, b, op, opt, prec, prect, iter, residue, threshold, showprogress, converged)
    integer, target, intent(in)    :: np     !< number of points
    R_TYPE,          intent(inout) :: x(:)   !< initial guess and the result
    R_TYPE,          intent(in)    :: b(:)   !< the right side
    interface
      subroutine op(x, y)                    !< the matrix A as operator y <- A*x
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine opt(x, y)                   !< the transposed matrix A as operator y <- A^T*x
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine opt
    end interface
    interface
      subroutine prec(x, y)                  !< the preconditioner
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prec
    end interface
    interface
      subroutine prect(x, y)                 !< the transposed preconditioner
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prect
    end interface
    integer,           intent(inout) :: iter         !< (in) number of max iterations, (out) used iterations
    FLOAT, optional,   intent(out)   :: residue      !< residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold    !< convergence threshold
    logical, optional, intent(in)    :: showprogress !< show progress bar
    logical, optional, intent(out)   :: converged    !< is the algorithm converged

    PUSH_SUB(X(qmr_spec_dotu))

    np_p => np
    call X(qmr_gen_dotu)(np, x, b, op, opt, X(dotu_qmr), X(nrm2_qmr), &
      prec, prect, iter, residue, threshold, showprogress, converged)

    POP_SUB(X(qmr_spec_dotu))
  end subroutine X(qmr_spec_dotu)


  ! ---------------------------------------------------------
  FLOAT function X(nrm2_qmr)(x)
    R_TYPE, intent(in) :: x(:)

    PUSH_SUB(X(nrm2_qmr))

    X(nrm2_qmr) = lalg_nrm2(np_p, x(:))

    POP_SUB(X(nrm2_qmr))
  end function X(nrm2_qmr)


  ! ---------------------------------------------------------
  !> for complex symmetric matrices
  !! W Chen and B Poirier, J Comput Phys 219, 198-209 (2006)
  subroutine X(qmr_sym_gen_dotu)(np, x, b, op, dotu, nrm2, prec, iter, &
    residue, threshold, showprogress, converged)
    integer, intent(in)    :: np    !< number of points
    R_TYPE,  intent(inout) :: x(:)  !< the initial guess and the result
    R_TYPE,  intent(in)    :: b(:)  !< the right side
    interface
      subroutine op(x, y)           !< the matrix A as operator
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      R_TYPE function dotu(x, y)    !< the dot product (must be x^T*y, not daggered)
        implicit none
        R_TYPE, intent(in) :: x(:)
        R_TYPE, intent(in) :: y(:)
      end function dotu
    end interface
    interface
      FLOAT function nrm2(x)        !< the 2-norm of the vector x
        implicit none
        R_TYPE, intent(in) :: x(:)
      end function nrm2
    end interface
    interface
      subroutine prec(x, y)         !< the preconditioner
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prec
    end interface
    integer,           intent(inout) :: iter         !< (in) the maximum number of iterations, (out) used iterations
    FLOAT, optional,   intent(out)   :: residue      !< the residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold    !< convergence threshold
    logical, optional, intent(in)    :: showprogress !< should there be a progress bar
    logical, optional, intent(out)   :: converged    !< has the algorithm converged

    R_TYPE, allocatable :: r(:), v(:), z(:), q(:), p(:), deltax(:), deltar(:)
    R_TYPE              :: eta, delta, epsilon, beta, rtmp
    FLOAT               :: rho, xsi, gamma, alpha, theta, threshold_, res, oldtheta, oldgamma, oldrho, tmp, norm_b
    integer             :: max_iter, err, ip, ilog_res, ilog_thr
    logical             :: showprogress_

    PUSH_SUB(X(qmr_sym_gen_dotu))

    if(present(converged)) converged = .false.
    threshold_ = optional_default(threshold, CNST(1.0e-6))
    showprogress_ = optional_default(showprogress, .false.)

    SAFE_ALLOCATE(r(1:np))
    SAFE_ALLOCATE(v(1:np))
    SAFE_ALLOCATE(z(1:np))
    SAFE_ALLOCATE(q(1:ubound(x, 1)))
    SAFE_ALLOCATE(p(1:np))
    SAFE_ALLOCATE(deltax(1:np))
    SAFE_ALLOCATE(deltar(1:np))

    !The initial starting point is zero
    x = M_ZERO
    call lalg_copy(np, b, r)
    call lalg_copy(np, b, v)


    rho      = nrm2(v)
    norm_b   = nrm2(b)

    max_iter = iter
    iter     = 0
    err      = 0
    res      = rho

    ! If rho is basically zero we are already done.
    if(abs(rho) > M_EPSILON) then
      call prec(v, z)

      xsi = nrm2(z)

      gamma = M_ONE
      eta   = -M_ONE
      alpha = M_ONE
      theta = M_ZERO

      ! initialize progress bar
      if(showprogress_) then
        ilog_thr = max(M_ZERO, -CNST(100.0)*log(threshold_))
        call loct_progress_bar(-1, ilog_thr)
      end if

      do while(iter < max_iter)
        iter = iter + 1
        if((abs(rho) < M_EPSILON) .or. (abs(xsi) < M_EPSILON)) then
          err = 1
          exit
        end if

        alpha = alpha*xsi/rho

        call lalg_scal(np, M_ONE/rho, v)
  
        call lalg_scal(np, M_ONE/xsi, z)

        delta = dotu(v, z)

        if(abs(delta) < M_EPSILON) then
          err = 2
          exit
        end if

        if(iter == 1) then
          call lalg_copy(np, z, q)
        else
          rtmp = -rho*delta/epsilon
          do ip = 1, np
            q(ip) = rtmp*q(ip) + z(ip)
          end do
        end if

        call op(q, p)
        call lalg_scal(np, alpha, p)

        epsilon = dotu(q, p)

        if(abs(epsilon) < M_EPSILON) then
          err = 3
          exit
        end if

        beta = epsilon/delta
        do ip = 1, np
          v(ip) = -beta*v(ip) + p(ip)
        end do

        oldrho = rho
        rho = nrm2(v)

        call prec(v, z)
        call lalg_scal(np, M_ONE/alpha, z)

        xsi = nrm2(z)

        oldtheta = theta
        theta    = rho/(gamma*abs(beta))

        oldgamma = gamma
        gamma    = M_ONE/sqrt(M_ONE+theta**2)

        if(abs(gamma) < M_EPSILON) then
          err = 4
          exit
        end if

        eta = -eta*oldrho*gamma**2/(beta*oldgamma**2)

        rtmp = eta*alpha
        if(iter == 1) then

          do ip = 1, np
            deltax(ip) = rtmp*q(ip)
            x(ip) = x(ip) + deltax(ip)
          end do

          do ip = 1, np
            deltar(ip) = eta*p(ip)
            r(ip) = r(ip) - deltar(ip)
          end do

        else

          tmp  = (oldtheta*gamma)**2
          do ip = 1, np
            deltax(ip) = tmp*deltax(ip) + rtmp*q(ip)
            x(ip) = x(ip) + deltax(ip)
          end do

          do ip = 1, np
            deltar(ip) = tmp*deltar(ip) + eta*p(ip)
            r(ip) = r(ip) - deltar(ip)
          end do

        end if

        ! avoid divide by zero
        if(abs(norm_b) < M_EPSILON) then
          res = M_HUGE
        else
          res = nrm2(r)/norm_b
        end if

        if(showprogress_) then
          ilog_res = CNST(100.0)*max(M_ZERO, -log(res))
          call loct_progress_bar(ilog_res, ilog_thr)
        end if

        if(res < threshold_) exit
      end do
    end if

    select case(err)
    case(0)
      if(res < threshold_) then
        if (present(converged)) converged = .true.
      end if
    case(1)
      write(message(1), '(a)') "QMR breakdown, cannot continue: b or P*b is the zero vector!"
    case(2)
      write(message(1), '(a)') "QMR breakdown, cannot continue: v^T*z is zero!"
    case(3)
      write(message(1), '(a)') "QMR breakdown, cannot continue: q^T*p is zero!"
    case(4)
      write(message(1), '(a)') "QMR breakdown, cannot continue: gamma is zero!"
    case(5)
    end select

    if (err>0) then
      write(message(2), '(a)') "Try to change some system parameters (e.g. Spacing, TDTimeStep, ...)."
      call messages_fatal(2)
    end if

    if(showprogress_) write(*,*) ''

    if(present(residue)) residue = res

    SAFE_DEALLOCATE_A(r)
    SAFE_DEALLOCATE_A(v)
    SAFE_DEALLOCATE_A(z)
    SAFE_DEALLOCATE_A(q)
    SAFE_DEALLOCATE_A(p)
    SAFE_DEALLOCATE_A(deltax)
    SAFE_DEALLOCATE_A(deltar)

    POP_SUB(X(qmr_sym_gen_dotu))
  end subroutine X(qmr_sym_gen_dotu)

  ! ---------------------------------------------------------
  !> for general complex matrices
  !! taken from 'An Implementation of the QMR Method based on Coupled Two-Term Recurrences' by
  !! R. W. Freund and N. M. Nachtigal (page 25)
  !! http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950017192_1995117192.pdf
  subroutine X(qmr_gen_dotu)(np, x, b, op, opt, dotu, nrm2, prec, prect, iter, &
    residue, threshold, showprogress, converged)
    integer, intent(in)    :: np      !< number of points
    R_TYPE,  intent(inout) :: x(:)    !< initial guess and result
    R_TYPE,  intent(in)    :: b(:)    !< right side
    interface
      subroutine op(x, y)             !< the matrix A as operator: y <- A*x
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine opt(x, y)            !< the transposed matrix A as operator: y <- A^T*x
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine opt
    end interface
    interface
      R_TYPE function dotu(x, y)      !< the dot product
        implicit none
        R_TYPE, intent(in) :: x(:)
        R_TYPE, intent(in) :: y(:)
      end function dotu
    end interface
    interface
      FLOAT function nrm2(x)          !< the 2-norm of a vector
        implicit none
        R_TYPE, intent(in) :: x(:)
      end function nrm2
    end interface
    interface
      subroutine prec(x, y)           !< preconditioner
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prec
    end interface
    interface
      subroutine prect(x, y)          !< transposed preconditioner
        implicit none
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prect
    end interface
    integer,           intent(inout) :: iter         !< (in) the maximum number of iterations, (out) used iterations
    FLOAT, optional,   intent(out)   :: residue      !< the residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold    !< convergence threshold
    logical, optional, intent(in)    :: showprogress !< should there be a progress bar
    logical, optional, intent(out)   :: converged    !< has the algorithm converged

    R_TYPE, allocatable :: r(:), v(:), w(:), z(:), q(:), p(:), deltax(:), tmp(:)
    R_TYPE              :: eta, delta, epsilon, beta
    FLOAT               :: rho, xsi, gamma, theta, threshold_, res, oldtheta, oldgamma, oldrho, nw, norm_b
    integer             :: max_iter, err, ilog_res, ilog_thr
    logical             :: showprogress_

    PUSH_SUB(X(qmr_gen_dotu))

    if(present(converged)) converged = .false.
    threshold_ = optional_default(threshold, CNST(1.0e-6))
    showprogress_ = optional_default(showprogress, .false.)

    SAFE_ALLOCATE(r(1:np))
    SAFE_ALLOCATE(v(1:np))
    SAFE_ALLOCATE(w(1:np))
    SAFE_ALLOCATE(z(1:np))
    SAFE_ALLOCATE(q(1:np))
    SAFE_ALLOCATE(p(1:ubound(x, 1)))
    SAFE_ALLOCATE(tmp(1:np))
    SAFE_ALLOCATE(deltax(1:np))

    ! r = b-Ax
    call lalg_copy(np, b, r)
    call op(x, tmp)
    call lalg_axpy(np, -M_ONE, tmp, r)
    rho = nrm2(r)
    max_iter = iter
    iter     = 0
    err      = 0
    res      = rho

! If rho is basically zero we are already done.
    if(abs(rho) > M_EPSILON) then
      call lalg_copy(np, r, v)
      call lalg_scal(np, M_ONE/rho, v)
      call lalg_copy(np, r, w)
      call prect(w, z)
      xsi = nrm2(z)
      call lalg_scal(np, M_ONE/xsi, w)
      epsilon = M_ONE
      gamma = M_ONE
      xsi = M_ONE
      theta = M_ZERO
      eta   = -M_ONE

      ! initialize progress bar
      if(showprogress_) then
        ilog_thr = max(M_ZERO, -CNST(100.0)*log(threshold_))
        call loct_progress_bar(-1, ilog_thr)
      end if

      do while(iter < max_iter)
        iter = iter + 1
        call prec(v, z)
        delta = dotu(w, z)
        if(abs(delta) < M_EPSILON) then
          err = 2
          exit
        end if

        if(iter == 1) then
          call lalg_copy(np, z, p)
          call prect(w, q)
        else
          call lalg_scal(np, -xsi*delta/epsilon, p)
          call lalg_axpy(np, M_ONE, z, p)
          call prect(w, z)
          call lalg_scal(np, -rho*delta/epsilon, q)
          call lalg_axpy(np, M_ONE, z, q)
        end if

        call op(p,tmp)
        epsilon = dotu(q, tmp)
        if(abs(epsilon) < M_EPSILON) then
          err = 3
          exit
        end if
        beta = epsilon/delta
        call lalg_scal(np, -beta, v)
        call lalg_axpy(np, M_ONE, tmp, v)
        oldrho = rho
        rho = nrm2(v)
        nw = nrm2(w)
        call lalg_scal(np, -beta, w)
        call opt(q,tmp)
        call lalg_axpy(np, M_ONE, tmp, w)
        xsi = nrm2(z)

        oldtheta = theta
        theta    = nrm2(w)*rho/(gamma*abs(beta)*nw)
        theta    = rho/(gamma*abs(beta))
        oldgamma = gamma
        gamma    = M_ONE/sqrt(M_ONE+theta**2)
        if(abs(gamma) < M_EPSILON) then
          err = 4
          exit
        end if
        eta = -eta*oldrho*gamma**2/(beta*oldgamma**2)

        if(iter == 1) then
          call lalg_copy(np, p, deltax)
          call lalg_scal(np, eta, deltax)
          call lalg_copy(np, deltax, x)
        else
          call lalg_scal(np, (oldtheta*gamma)**2, deltax)
          call lalg_axpy(np, eta, p, deltax)
          call lalg_axpy(np, M_ONE, deltax, x)
        end if
        call lalg_scal(np, M_ONE/rho, v)
        call lalg_scal(np, M_ONE/xsi, w)

        call lalg_copy(np, b, r)
        call op(x, tmp)
        call lalg_axpy(np, -M_ONE, tmp, r)

        norm_b = nrm2(x)
        ! avoid divide by zero
        if(norm_b < M_EPSILON) then
          res = M_HUGE
        else
          res = nrm2(r)/norm_b
        end if

        if(showprogress_) then
          ilog_res = max(M_ZERO, -CNST(100.0)*log(res))
          call loct_progress_bar(ilog_res, ilog_thr)
        end if
        
        if(res < threshold_) exit
      end do
    end if

    select case(err)
    case(0)
      if(res < threshold_) then
        if (present(converged)) converged = .true.
      else
        write(message(1), '(a)') "QMR solver not converged!"
        write(message(2), '(a)') "Try increasing the maximum number of iterations or the tolerance."
        call messages_warning(2)
      end if
    case(1)
      write(message(1), '(a)') "QMR failure, can't continue: b or P*b is the zero vector!"
    case(2)
      write(message(1), '(a)') "QMR failure, can't continue: z^T*y is zero!"
    case(3)
      write(message(1), '(a)') "QMR failure, can't continue: q^T*p is zero!"
    case(4)
      write(message(1), '(a)') "QMR failure, can't continue: gamma is zero!"
    end select

    if (err>0) then
      write(message(2), '(a)') "Try to change some system parameters (e.g. Spacing, TDTimeStep, ...)."
      call messages_fatal(2)
    end if

    if(showprogress_) write(*,*) ''

    if(present(residue)) residue = res

    SAFE_DEALLOCATE_A(r)
    SAFE_DEALLOCATE_A(v)
    SAFE_DEALLOCATE_A(w)
    SAFE_DEALLOCATE_A(z)
    SAFE_DEALLOCATE_A(q)
    SAFE_DEALLOCATE_A(p)
    SAFE_DEALLOCATE_A(deltax)
    SAFE_DEALLOCATE_A(tmp)

    POP_SUB(X(qmr_gen_dotu))
  end subroutine X(qmr_gen_dotu)



!> This is the "Induced Dimension Reduction", IDR(s) (for s=4). IDR(s) is a robust and efficient short recurrence 
!> Krylov subspace method for solving large nonsymmetric systems of linear equations. It is described in 
!> [Peter Sonneveld and Martin B. van Gijzen, SIAM J. Sci. Comput. 31, 1035 (2008)]. 
!>
!> We have adapted the code released by M. B. van Gizjen [http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html].
!> That code is licenced under the MIT licence, which allows to re-release this modified version under the GPL.
!> 
!> The original copyright notice, and licence, follows:
!>
!> Copyright (c) July 2015, Martin van Gijzen
!>
!> Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
!> documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
!> the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
!> and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
!>
!> The above copyright notice and this permission notice shall be included in all copies or substantial 
!> portions of the Software.
!>
!> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
!> TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
!> THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
!> CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
!> DEALINGS IN THE SOFTWARE.
  function X(idrs)( b, s, &
                   preconditioner, matrixvector, &
                   ddotprod, zdotprod, &
                   tolerance, maximum_iterations, variant, &
                   flag, relres, iterations, &
                   x0, U0, omega, resvec, H)

    implicit none

    ! Required input parameters:
    R_TYPE, intent(in) :: b(:, :)  ! system rhs
    integer, intent(in) :: s       ! s parameter
    ! Solution:
    R_TYPE :: X(idrs) (size(b,1), size(b,2))
    ! Optional input parameters:
    FLOAT, optional, intent(in) :: tolerance 
    integer, optional, intent(in) :: maximum_iterations 
    character(len=8), optional, intent(in) :: variant     
   ! Optional output parameters:
    integer, optional, intent(out) :: flag        
    FLOAT, optional, intent(out) :: relres     
    integer, optional, intent(out) :: iterations
    ! Optional input arrays:
    R_TYPE, optional, intent(in) :: x0(:, :)  
    R_TYPE, optional, intent(in) :: U0(:, :, :) 
    R_TYPE, optional, intent(in) :: omega(:) 
    ! Optional output arrays:
    FLOAT, optional, intent(out) :: resvec(:) 
    R_TYPE, optional, intent(out) :: H(:, :)   

    interface
      function preconditioner( v )
        R_TYPE, dimension(:,:), intent(in)     :: v
        R_TYPE, dimension(size(v,1),size(v,2)) :: preconditioner
      end function preconditioner
      function matrixvector( v ) 
        R_TYPE, intent(in)       :: v(:,:)
        R_TYPE                   :: matrixvector(size(v,1),size(v,2))
      end function matrixvector
      FLOAT function ddotprod(a, b)
        FLOAT, intent(in) :: a(:), b(:)
        FLOAT, allocatable :: apsi(:, :), bpsi(:, :)
      end function ddotprod
      CMPLX function zdotprod(a, b)
        CMPLX, intent(in) :: a(:), b(:)
        CMPLX, allocatable :: apsi(:, :), bpsi(:, :)
      end function zdotprod
    end interface
    
    ! Local arrays:
    FLOAT, allocatable           :: P(:,:,:) 
    R_TYPE, allocatable          :: R0(:,:) 
    R_TYPE                       :: x(size(b,1),size(b,2))
    R_TYPE                       :: G(size(b,1),size(b,2),s)
    R_TYPE                       :: U(size(b,1),size(b,2),s)
    R_TYPE                       :: r(size(b,1),size(b,2)) 
    R_TYPE                       :: v(size(b,1),size(b,2))   
    R_TYPE                       :: t(size(b,1),size(b,2))  
    R_TYPE                       :: M(s,s), f(s), mu(s)
    R_TYPE                       :: alpha(s), beta(s), gamma(s)

    R_TYPE                       :: om, tr    
    FLOAT                        :: nr, nt, rho, kappa

    ! Declarations:
    integer               :: n                 ! dimension of the system
    integer               :: nrhs              ! Number of RHS-vectors
    integer               :: maxit             ! maximum number of iterations
    integer               :: method            ! which IDR(s) variant?
    FLOAT                 :: tol               ! actual tolerance
    integer               :: info              ! convergence indicator
    logical               :: out_flag          ! store flag
    logical               :: out_relres        ! store relres
    logical               :: out_iterations    ! store number of iterations
    logical               :: inispace          ! initial search space
    logical               :: user_omega        ! user defined omega present
    integer               :: n_omega           ! number of user defined omegas
    logical               :: out_resvec        ! store residual norms
    logical               :: out_H             ! store iteration parameters in H
    integer               :: nritz             ! Number of wanted ritz values

    integer               :: iter              ! number of iterations
    integer               :: ii                ! inner iterations index
    integer               :: jj                ! G-space index
    FLOAT                 :: normb, normr, tolb! for tolerance check
    integer               :: i,j,k,l           ! loop counters

    ! Problem size:
    n    = size(b,1)
    ! Number of right-hand-side vectors:
    nrhs = size(b,2)
      
    ! Check optional input parameters:
    if ( present(tolerance) ) then
      if ( tolerance < 0 ) then
        message(1) = 'The tolerance parameter in idrs routine must be non-negative'
        call messages_fatal(1)
      end if
      tol = tolerance 
    else
      tol = CNST(1e-6)
    endif

    maxit=min(2*n,1000)
    if ( present(maximum_iterations) ) maxit = maximum_iterations 
   
    method = 1 ! biortho   
    if ( present(variant) ) then
      if ( variant == 'minsync' ) then
        method = 2
      elseif ( variant == 'bicgstab' ) then 
        method = 3
      endif
    endif

    ! Initialize the output variables 
    out_flag       = present(flag)
    if ( out_flag )       flag = -1 
    out_relres     = present(relres)
    if ( out_relres)      relres = 1. 
    out_iterations = present(iterations)
    if ( out_iterations ) iterations = 0 

    ! Check optional input arrays:
    x = M_ZERO
    if ( present(x0) ) x = x0
      
    U = M_ZERO
    inispace =  present(U0)
    if ( inispace ) U = U0

    user_omega = present(omega)
    if ( user_omega ) then
      n_omega = size(omega)
    end if

    ! Check output arrays
    out_resvec     = present(resvec)
    if ( out_resvec ) then
      if ( maxit+1 > size(resvec) ) then
        message(1) = 'idrs: Length of vector with residual norms too small, should be maxit+1'
        call messages_fatal(1)
      end if
    end if

    out_H = present(H)
    if ( out_H ) then
      nritz = size(H,1)-1
      if ( size(H,2) /= nritz ) then
        message(1) = 'Second dimension of H incompatible, with first'
        call messages_fatal(1)
      end if
      H = M_ZERO
    end if

    ! compute initial residual, set absolute tolerance
    normb = X(frob_norm)(b)
    ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
    ! we remove this feature
    !tolb = tol * normb
    tolb = tol
    r = b - matrixvector(x)
    normr = X(frob_norm)(r)
    if ( out_resvec ) resvec(1)= normr

    ! check if the initial solution is not already a solution within the prescribed
    ! tolerance
    if (normr <= tolb) then      
      if ( out_iterations ) iterations = 0               
      if ( out_flag )       flag  = 0
      ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
      ! we remove this feature
      !if ( out_relres )     relres = normr/normb
      if ( out_relres )     relres = normr
      return
    end if

    ! Define P and kappa (depending on the method)
    if ( method == 1 ) then
      allocate( P(n,nrhs,s) )
      call RANDOM_SEED
      call RANDOM_NUMBER(P)
      do j = 1,s
         do k = 1,j-1
            alpha(k) = dtrace_dot( P(:,:,k),P(:,:,j) )
            P(:,:,j) = P(:,:,j) - alpha(k)*P(:,:,k)
         end do
         P(:,:,j) = P(:,:,j)/dfrob_norm(P(:,:,j))
      end do
      kappa = CNST(0.7)
    elseif ( method == 2 ) then
    ! P is piecewise constant, minimum residual for omega
      kappa = M_ZERO
    elseif ( method == 3 ) then
      !if ( s /= 1 ) stop "s=1 is required for variant bicgstab"
      ASSERT( s .ne. 1 )
      allocate( R0(n,nrhs) )
      R0 = r
      kappa = M_ZERO
    endif

    ! Initialize local variables:
    M = M_ZERO
    om = M_ONE
    iter = 0
    info = -1
    jj = 0
    ii = 0
    
    ! This concludes the initialisation phase

    ! Main iteration loop, build G-spaces:
    
    do while (  info < 0 )  ! start of iteration loop
     
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Generate s vectors in G_j
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ! New right-hand side for small system:
      f = X(p_dot)( P, R0, r, s )

      do k = 1, s

         ! Update inner iteration counter
         ii = ii + 1

         ! Compute new v
         v = r 
         if ( jj > 0 ) then

            ! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
            do i = k,s
               gamma(i) = f(i)
               do j = k, i-1
                  gamma(i) = gamma(i) - M(i,j)*gamma(j)
               end do
               gamma(i) = gamma(i)/M(i,i)
               v = v - gamma(i)*G(:,:,i)
            end do

            ! Compute new U(:,:,k)
            t = om * preconditioner(v)
            do i = k,s
               t = t + gamma(i)*U(:,:,i)
            end do
            U(:,:,k) = t

            ! Compute Hessenberg matrix?
            if ( out_H .and. ii <= nritz ) &
               H(ii-s:ii-k,ii)   = -gamma(k:s)/beta(k:s)

         else if ( .not. inispace ) then

            ! Updates for the first s iterations (in G_0):
            U(:, :, k) = preconditioner(v)

         end if

         ! Compute new G(:,:,k), G(:,:,k) is in space G_j
         G(:, :, k) = matrixvector(U(:, :, k))
        
         ! Bi-Orthogonalise the new basis vectors: 
         mu = X(p_dot)( P, R0, G(:,:,k), s )
         do i = 1,k-1
            alpha(i) = mu(i)
            do j = 1, i-1
               alpha(i) = alpha(i) - M(i,j)*alpha(j)
            end do
            alpha(i) = alpha(i)/M(i,i)
            G(:,:,k) = G(:,:,k) - G(:,:,i)*alpha(i)
            U(:,:,k) = U(:,:,k) - U(:,:,i)*alpha(i)
            mu(k:s)  = mu(k:s)  - M(k:s,i)*alpha(i)
         end do
         M(k:s,k) = mu(k:s)

         ! Compute Hessenberg matrix?
         if ( out_H .and. ii <= nritz .and. k  > 1 ) &
            H(ii-k+1:ii-1,ii) =  alpha(1:k-1)/beta(1:k-1)

         ! Break down?
         if ( abs(M(k,k)) <= tiny(tol) ) then
            info = 3
            exit
         end if

         ! Make r orthogonal to p_i, i = 1..k, update solution and residual 
         beta(k) = f(k)/M(k,k)
         r = r - beta(k)*G(:,:,k)
         x = x + beta(k)*U(:,:,k)

         ! New f = P(prime) *r (first k  components are zero)
         if ( k < s ) then
            f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k)
         end if

         ! Compute Hessenberg matrix?
         if ( out_H .and. ii <= nritz ) then     
            H(ii,ii) = 1./beta(k)
            l = max(1,ii-s)
            H(l+1:ii+1,ii) = (H(l+1:ii+1,ii) - H(l:ii,ii))
            H(l:ii+1,ii)   = H(l:ii+1,ii)/om
         end if

         ! Check for convergence
         normr = X(frob_norm)(r)
         iter = iter + 1
         if ( out_resvec ) resvec(iter + 1) = normr
         if ( normr < tolb ) then
            info = 0
            exit
         elseif ( iter == maxit ) then
            info = 1
            exit
         end if 

      end do ! Now we have computed s+1 vectors in G_j
      if ( info >= 0 )  then
         exit
      end if

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute first residual in G_j+1
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ! Update G-space counter
      jj = jj + 1

      ! Compute first residual in G_j+1
      ! Note: r is already perpendicular to P so v = r
 
      ! Preconditioning:
      v = preconditioner(r)
      t = matrixvector(v)


      ! Computation of a new omega
      if ( user_omega ) then
        i = mod(jj,n_omega)
        if ( i == 0 ) i = n_omega
        om = omega(i)
      elseif ( abs(kappa) <= M_EPSILON  ) then

        ! Minimal residual (same as in Bi-CGSTAB):
        om = X(trace_dot)(t,r)/X(trace_dot)(t,t)
      else

        ! 'Maintaining the convergence':
        nr = X(frob_norm)(r)
        nt = X(frob_norm)(t)
        tr = X(trace_dot)(t,r)
        rho = abs(tr/(nt*nr))
        om=tr/(nt*nt)
        if ( rho < kappa ) then
          om = om*kappa/rho
        end if
      end if
      if ( abs(om) <= epsilon(tol) ) then 
         info = 3
         exit
      end if 

      ! Update solution and residual
      r = r - om*t 
      x = x + om*v 

      ! Check for convergence
      normr = X(frob_norm)(r)
      iter = iter + 1
      if ( out_resvec ) resvec(iter + 1) = normr
      if ( normr < tolb ) then
        info = 0
      elseif ( iter == maxit ) then
        info = 1
      end if 

    end do ! end of while loop

    ! Set output parameters
    r = b - matrixvector(x)
    normr = X(frob_norm)(r)

    if ( info == 0 .and. normr > tolb ) info = 2
    if ( out_iterations ) iterations = iter
    ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
    ! we remove this feature
    !if ( out_relres )     relres=normr/normb
    if ( out_relres )     relres=normr
    if ( out_flag )       flag = info

    X(idrs) = x

    contains

    function dtrace_dot(v, w)
      ! Trace inner product of complex matrices
      FLOAT, intent(in)      :: v(:,:), w(:,:)
      FLOAT                  :: dtrace_dot
      integer k
      dtrace_dot = M_ZERO
      do k = 1, size(v, 2)
        dtrace_dot = dtrace_dot + ddotprod(v(:, k), w(:, k))
      end do
    end function dtrace_dot

    function ztrace_dot(v, w)
      ! Trace inner product of complex matrices
      CMPLX, intent(in)      :: v(:,:), w(:,:)
      CMPLX                  :: ztrace_dot
      integer :: k
      ztrace_dot = M_z0
      do k = 1, size(v, 2)
        ztrace_dot = ztrace_dot + zdotprod(v(:, k), w(:, k))
      end do
    end function ztrace_dot

    function X(p_dot)(P, R0, w, s)
     ! P inner product of complex matrices
     FLOAT,    allocatable, intent(in) :: P(:,:,:) 
     R_TYPE, allocatable, intent(in) :: R0(:,:)
     R_TYPE, intent(in)              :: w(:,:)
     integer, intent(in)             :: s
     R_TYPE                          :: X(p_dot)(s)

     R_TYPE                          :: v(s)
     integer                         :: j, k, N, low(s), up(s), step, nrhs

     if ( allocated(P) ) then
       ! Biortho: P has orthogonal random numbers
       do i = 1, s
         v(i) = R_TOTYPE(M_ZERO)
         do k = 1, size(w, 2)
#if defined(R_TREAL)
           v(i) = v(i) + ddotprod(P(:, k, i), w(:, k))
#else
           v(i) = v(i) + TOCMPLX(ddotprod(P(:, k, i), R_REAL(w(:, k))), ddotprod(P(:, k, i), R_AIMAG(w(:, k)))) 
#endif
         end do
       end do
     else if ( allocated(R0) ) then
       ! BiCGSTAB: shadow vector equal to initial residual
       v(1) = R_TOTYPE(M_ZERO)
       do k = 1, size(w, 2)
         v(1) = v(1) + X(dotprod)(R0(:, k), w(:, k))
       end do
     else
        ! Minsync: P is piecewise constant 
        ! WARNING: the integrals are done here in a peculiar way, not consistent with the
        ! definition of the dot product. So this is probably not working.
        N    = size(w,1)
        nrhs = size(w,2)
        step = N/s
        low(1) = 1
        do i = 1,s-1
            low(i+1) = i*step+1
            up(i) = i*step
        end do
        up(s) = N

        do i = 1,s
          v(i)  = M_ZERO
          do j = 1, nrhs
            v(i) = v(i) + sum( w(low(i):up(i),j) )
          end do
        end do
     end if
      
      X(p_dot) = v
    end function X(p_dot)

    function dfrob_norm(v)
     ! Frobenius norm of complex matrix
      FLOAT, intent(in)      :: v(:,:)
      FLOAT                     :: dfrob_norm
      integer :: k
      dfrob_norm = M_ZERO
      do k = 1, size(v, 2)
        dfrob_norm = dfrob_norm + ddotprod( v(:, k), v(:, k) )
      end do
      dfrob_norm = sqrt( dfrob_norm )
    end function dfrob_norm

    function zfrob_norm(v)
      ! Frobenius norm of complex matrix
      CMPLX, intent(in)      :: v(:,:)
      FLOAT                  :: zfrob_norm
      integer :: k
      zfrob_norm = M_ZERO
      do k = 1, size(v, 2)
        zfrob_norm = zfrob_norm + zdotprod( v(:, k), v(:, k) )
      end do
      zfrob_norm = sqrt( zfrob_norm )
    end function zfrob_norm
   
  end function X(idrs)



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
