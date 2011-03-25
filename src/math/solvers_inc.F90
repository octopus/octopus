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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$


! The two following subroutines, sym_conjugate_gradients, and bi_conjugate_gradients,
! must be called under a common interface: conjugate_gradients. It provides an
! approximate solution to the linear system problem  Ax = b.
! Solving a symmetric linear system, which is either real or complex symmetric or
! Hermitian the best choice is sym_conjugate_gradients (does not need A^\dagger)
! Solving a real unsymmetric or a complex non-Hermitian system bi_conjugate_gradients
! has to be chosen (where one does need A^\dagger).
!
! Note: the complex-valued versions (both CG and BiCG) work only with
!       symmetric operators but they may be non-Hermitian. This is a
!       property of the CG algorithm. A comment on this can be found
!       in chapter 4.1 of ftp://ftp.netlib.org/templates/templates.ps
!
! subroutine conjugate_gradients(np, x, b, op, [opt,] dotp, iter [, residue] [, threshold] )
!    integer, intent(in)  :: np      => The dimension of the problem.
!    FLOAT, intent(inout) :: x       => On input, an estimate to the solution.
!                                    => On output, the approximate solution.
!    FLOAT, intent(in)    :: b       => The inhomogeneous term of the equation.
!    interface
!      subroutine op(x, y)
!         FLOAT, intent(in)  :: x(:)
!         FLOAT, intent(out) :: y(:)
!      end subroutine op
!    end interface                   => This should be a procedure that
!                                       computes Ax = y.
!    interface
!      subroutine opt(x, y)
!         FLOAT, intent(in)  :: x(:)
!         FLOAT, intent(out) :: y(:)
!      end subroutine opt
!    end interface                   => If present, this should be a procedure that
!                                       computes A^\dagger x = y.
!                                       Only useful for non-Hermitian operators.
!    interface
!      R_TYPE function dotp(x, y)
!      R_TYPE, intent(inout) :: x(:)
!      R_TYPE, intent(in)    :: y(:)
!    end function dotp               => Calculates the <x | y>.
!                                       Depending on the matrix A one should choose:
!                                       complex symmetric: <x | y> = x^T * y
!                                       hermitian:         <x | y> = x^\dagger * y
!                                       general:           <x | y> = x^\dagger * y
!    integer, intent(inout) :: iter  => On input, the maximum number of iterations that
!                                       the procedure is allowed to take.
!                                       On output, the iterations it actually did.
!    FLOAT, intent(out) :: residue   => If present, it measures the final error:
!                                       residue = <Ax - b | Ax - b>.
!    FLOAT, intent(in)  :: threshold => If present, it sets the required accuracy
!                                       threshold for the algorithm to stop. If not
!                                       present, this is set to 1.0e-6. [The algorithm
!                                       stops when <Ax - b | Ax - b> <= threshold, or
!                                       iter iterations are reached.]
! end subroutine conjugate_gradients
!
! (*) NOTE: The algorithm assumes that the vectors are given in an orthonormal basis.


! ---------------------------------------------------------
subroutine X(sym_conjugate_gradients)(np_part, np, x, b, op, dotp, iter, residue, threshold)
  integer, intent(in)    :: np_part
  integer, intent(in)    :: np
  R_TYPE,  intent(inout) :: x(:)
  R_TYPE,  intent(in)    :: b(:)
  interface
    subroutine op(x, y)
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine op
    R_TYPE function dotp(x, y) result(res)
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

  if(present(threshold)) then
    threshold_ = threshold
  else
    threshold_ = CNST(1.0e-6)
  end if

  SAFE_ALLOCATE( r(1:np))
  SAFE_ALLOCATE(ax(1:np))
  SAFE_ALLOCATE( p(1:np_part))
  SAFE_ALLOCATE(ap(1:np))

  ! Initial residue.
  call op(x, ax)
  forall(ip = 1:np) r(ip) = b(ip) - ax(ip)

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
    forall(ip = 1:np) p(ip) = r(ip) + beta*p(ip)
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
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine op
  end interface
  interface
    subroutine opt(x, y)
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine opt
    R_TYPE function dotp(x, y) result(res)
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

  if(present(threshold)) then
    threshold_ = threshold
  else
    threshold_ = CNST(1.0e-6)
  end if

  SAFE_ALLOCATE(  r(1:np))
  SAFE_ALLOCATE( rr(1:np))
  SAFE_ALLOCATE( ax(1:np))
  SAFE_ALLOCATE(  p(1:np))
  SAFE_ALLOCATE( pp(1:np))
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
    integer, target,   intent(in)    :: np    ! number of points
    R_TYPE,             intent(inout) :: x(:)  ! initial guess and result
    R_TYPE,             intent(in)    :: b(:)  ! the right side
    interface
      subroutine op(x, y)                     ! the matrix A as operator, y <- Ax
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine prec(x, y)                   ! preconditioner
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prec
    end interface
    integer,           intent(inout) :: iter  ! [in] number of max iterations, [out] used iterations
    FLOAT, optional,   intent(out)   :: residue ! residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold ! convergence threshold
    logical, optional, intent(in)    :: showprogress ! show progress bar
    logical, optional, intent(out)   :: converged ! is the algorithm converged

    PUSH_SUB(X(qmr_sym_spec_dotu))

    np_p => np
    call X(qmr_sym_gen_dotu)(np, x, b, op, X(dotu_qmr), X(nrm2_qmr), prec, iter, &
      residue, threshold, showprogress, converged)

    POP_SUB(X(qmr_sym_spec_dotu))
  end subroutine X(qmr_sym_spec_dotu)

  ! ---------------------------------------------------------
  subroutine X(qmr_spec_dotu)(np, x, b, op, opt, prec, prect, iter, residue, threshold, showprogress, converged)
    integer, target,   intent(in)    :: np    ! number of points
    R_TYPE,             intent(inout) :: x(:)  ! initial guess and the result
    R_TYPE,             intent(in)    :: b(:)  ! the right side
    interface
      subroutine op(x, y)                     ! the matrix A as operator y <- A*x
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine opt(x, y)                    ! the transposed matrix A as operator y <- A^T*x
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine opt
    end interface
    interface
      subroutine prec(x, y)                   ! the preconditioner
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prec
    end interface
    interface
      subroutine prect(x, y)                  ! the transposed preconditioner
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prect
    end interface
    integer,           intent(inout) :: iter  ! [in] number of max iterations, [out] used iterations
    FLOAT, optional,   intent(out)   :: residue ! residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold ! convergence threshold
    logical, optional, intent(in)    :: showprogress ! show progress bar
    logical, optional, intent(out)   :: converged ! is the algorithm converged

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
  ! for complex symmetric matrices
  ! W Chen and B Poirier, J Comput Phys 219, 198-209 (2006)
  subroutine X(qmr_sym_gen_dotu)(np, x, b, op, dotu, nrm2, prec, iter, &
    residue, threshold, showprogress, converged)
    integer,           intent(in)    :: np    ! number of points
    R_TYPE,             intent(inout) :: x(:)  ! the initial guess and the result
    R_TYPE,             intent(in)    :: b(:)  ! the right side
    interface
      subroutine op(x, y)                     ! the matrix A as operator
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine op
    end interface
    interface 
      R_TYPE function dotu(x, y)               ! the dot product (must be x^T*y, not daggered)
        R_TYPE, intent(in) :: x(:)
        R_TYPE, intent(in) :: y(:)
      end function dotu
    end interface
    interface 
      FLOAT function nrm2(x)                  ! the 2-norm of the vector x
        R_TYPE, intent(in) :: x(:)
      end function nrm2
    end interface
    interface
      subroutine prec(x, y)                   ! the preconditioner
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prec
    end interface
    integer,           intent(inout) :: iter  ! [in] the maximum number of iterations, [out] used iterations
    FLOAT, optional,   intent(out)   :: residue ! the residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold ! convergence threshold
    logical, optional, intent(in)    :: showprogress ! should there be a progress bar
    logical, optional, intent(out)   :: converged ! has the algorithm converged

    R_TYPE, allocatable :: r(:), v(:), z(:), q(:), p(:), deltax(:), deltar(:)
    R_TYPE              :: eta, delta, epsilon, beta, rtmp
    FLOAT               :: rho, xsi, gamma, alpha, theta, threshold_, res, oldtheta, oldgamma, oldrho, tmp
    integer             :: max_iter, err, ip
    logical             :: showprogress_
    FLOAT               :: log_res, log_thr, norm_b
    integer             :: ilog_res, ilog_thr

    PUSH_SUB(X(qmr_sym_gen_dotu))

    if(present(converged)) then
      converged = .false.
    end if
    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = CNST(1.0e-6)
    end if
    if(present(showprogress)) then
      showprogress_ = showprogress
    else
      showprogress_ = .false.
    end if

    SAFE_ALLOCATE(r(1:np))
    SAFE_ALLOCATE(v(1:np))
    SAFE_ALLOCATE(z(1:np))
    SAFE_ALLOCATE(q(1:np))
    SAFE_ALLOCATE(p(1:np))
    SAFE_ALLOCATE(deltax(1:np))
    SAFE_ALLOCATE(deltar(1:np))

    ! use v as temp var
    call op(x, v)

    forall (ip = 1:np) 
      r(ip) = b(ip) - v(ip)
      v(ip) = r(ip)
    end forall

    rho      = nrm2(v)
    norm_b   = nrm2(b)

    max_iter = iter
    iter     = 0
    err      = 0
    res      = rho

    ! If rho is basically zero we are already done.
    if(abs(rho).gt.M_EPSILON) then
      call prec(v, z)

      xsi = nrm2(z)

      gamma = M_ONE
      eta   = -M_ONE
      alpha = M_ONE
      theta = M_ZERO

      ! initialize progress bar
      log_thr = -log(threshold_)
      ilog_thr = M_TEN**2*log_thr
      if(showprogress_) call loct_progress_bar(-1, ilog_thr)

      do while(iter < max_iter)
        iter = iter + 1
        if((abs(rho) < M_EPSILON) .or. (abs(xsi) < M_EPSILON)) then
          err = 1
          exit
        end if
        alpha = alpha*xsi/rho
        tmp = M_ONE/rho
        forall (ip = 1:np) v(ip) = tmp*v(ip)
        tmp = M_ONE/xsi
        forall (ip = 1:np) z(ip) = tmp*z(ip)

        delta = dotu(v, z)

        if(abs(delta) < M_EPSILON) then
          err = 2
          exit
        end if
        if(iter == 1) then
          forall (ip = 1:np) q(ip) = z(ip)
        else
          rtmp = -rho*delta/epsilon
          forall (ip = 1:np) q(ip) = rtmp*q(ip) + z(ip)
        end if
        call op(q, p)
        forall (ip = 1:np) p(ip) = alpha*p(ip)

        epsilon = dotu(q, p)

        if(abs(epsilon) < M_EPSILON) then
          err = 3
          exit
        end if
        beta = epsilon/delta
        forall (ip = 1:np) v(ip) = -beta*v(ip) + p(ip)        
        oldrho = rho

        rho = nrm2(v)

        call prec(v, z)
        tmp = M_ONE/alpha
        forall (ip = 1:np) z(ip) = tmp*z(ip)

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

          forall (ip = 1:np) 
            deltax(ip) = rtmp*q(ip)
            x(ip) = x(ip) + deltax(ip)
          end forall

          forall (ip = 1:np)
            deltar(ip) = eta*p(ip)
            r(ip) = r(ip) - deltar(ip)
          end forall

        else

          tmp  = (oldtheta*gamma)**2
          forall (ip = 1:np) 
            deltax(ip) = tmp*deltax(ip) + rtmp*q(ip)
            x(ip) = x(ip) + deltax(ip)
          end forall

          forall (ip = 1:np) 
            deltar(ip) = tmp*deltar(ip) + eta*p(ip)
            r(ip) = r(ip) - deltar(ip)
          end forall

        end if

        ! avoid divide by zero
        if(abs(norm_b) < M_EPSILON) then
          res = M_HUGE
        else
          res = nrm2(r)/norm_b
        endif

        log_res = -log(res)
        if(log_res < 0) log_res = 0
        ilog_res = M_TEN**2*log_res
        if(showprogress_)  call loct_progress_bar(ilog_res, ilog_thr)
        if(res < threshold_) exit
      end do
    end if

    if((err.eq.0).and.(iter.eq.max_iter)) err = 5

    select case(err)
    case(0)
      ! converged
      if(present(converged)) then
        converged = .true.
      end if
    case(1)
      write(message(1), '(a)') "QMR breakdown, cannot continue: b or P*b is the zero vector!"
      call messages_fatal(1)
    case(2)
      write(message(1), '(a)') "QMR breakdown, cannot continue: v^T*z is zero!"
      call messages_fatal(1)
    case(3)
      write(message(1), '(a)') "QMR breakdown, cannot continue: q^T*p is zero!"
      call messages_fatal(1)
    case(4)
      write(message(1), '(a)') "QMR breakdown, cannot continue: gamma is zero!"
      call messages_fatal(1)
    case(5)
      write(message(1), '(a)') "QMR solver not converged!"
      write(message(2), '(a)') "Try increasing the maximum number of iterations or the tolerance."
      call messages_warning(2)
      if(present(converged)) then
        converged = .false.
      end if
    end select
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
  ! for general complex matrices
  ! taken from 'An Implementation of the QMR Method based on Coupled Two-Term Recurrences' by
  ! R. W. Freund and N. M. Nachtigal (page 25)
  ! http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950017192_1995117192.pdf
  subroutine X(qmr_gen_dotu)(np, x, b, op, opt, dotu, nrm2, prec, prect, iter, &
    residue, threshold, showprogress, converged)
    integer,            intent(in)    :: np    ! number of points
    R_TYPE,             intent(inout) :: x(:)  ! initial guess and result
    R_TYPE,             intent(in)    :: b(:)  ! right side
    interface
      subroutine op(x, y)                     ! the matrix A as operator: y <- A*x
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine opt(x, y)                    ! the transposed matrix A as operator: y <- A^T*x
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine opt
    end interface
    interface 
      R_TYPE function dotu(x, y)               ! the dot product
        R_TYPE, intent(in) :: x(:)
        R_TYPE, intent(in) :: y(:)
      end function dotu
    end interface
    interface 
      FLOAT function nrm2(x)                  ! the 2-norm of a vector
        R_TYPE, intent(in) :: x(:)
      end function nrm2
    end interface
    interface
      subroutine prec(x, y)                   ! preconditioner
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prec
    end interface
    interface
      subroutine prect(x, y)                  ! transposed preconditioner
        R_TYPE, intent(in)  :: x(:)
        R_TYPE, intent(out) :: y(:)
      end subroutine prect
    end interface
    integer,           intent(inout) :: iter  ! [in] the maximum number of iterations, [out] used iterations
    FLOAT, optional,   intent(out)   :: residue ! the residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold ! convergence threshold
    logical, optional, intent(in)    :: showprogress ! should there be a progress bar
    logical, optional, intent(out)   :: converged ! has the algorithm converged

    R_TYPE, allocatable :: r(:), v(:), w(:), z(:), q(:), p(:), deltax(:), tmp(:)
    R_TYPE              :: eta, delta, epsilon, beta
    FLOAT               :: rho, xsi, gamma, theta, threshold_, res, oldtheta, oldgamma, oldrho, nw, norm_b
    integer             :: max_iter, err
    logical             :: showprogress_
    FLOAT               :: log_res, log_thr
    integer             :: ilog_res, ilog_thr

    PUSH_SUB(X(qmr_gen_dotu))

    if(present(converged)) then
      converged = .false.
    end if
    if(present(threshold)) then
      threshold_ = threshold
    else
      threshold_ = CNST(1.0e-6)
    end if
    if(present(showprogress)) then
      showprogress_ = showprogress
    else
      showprogress_ = .false.
    end if

    SAFE_ALLOCATE(r(1:np))
    SAFE_ALLOCATE(v(1:np))
    SAFE_ALLOCATE(w(1:np))
    SAFE_ALLOCATE(z(1:np))
    SAFE_ALLOCATE(q(1:np))
    SAFE_ALLOCATE(p(1:np))
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
    if(abs(rho).gt.M_EPSILON) then
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
      log_thr = -log(threshold_)
      ilog_thr = M_TEN**2*log_thr
      if(showprogress_) call loct_progress_bar(-1, ilog_thr)

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
        endif

        log_res = -log(res)
        if(log_res < 0) log_res = 0
        ilog_res = M_TEN**2*log_res
        if(showprogress_)  call loct_progress_bar(ilog_res, ilog_thr)
        if(res < threshold_) exit
      end do
    end if

    if((err.eq.0).and.(iter.eq.max_iter)) err = 5

    select case(err)
    case(0)
      ! converged
      if(present(converged)) then
        converged = .true.
      end if
    case(1)
      write(message(1), '(a)') "QMR failure, can't continue: b or P*b is the zero vector!"
      call messages_fatal(1)
    case(2)
      write(message(1), '(a)') "QMR failure, can't continue: z^T*y is zero!"
      call messages_fatal(1)
    case(3)
      write(message(1), '(a)') "QMR failure, can't continue: q^T*p is zero!"
      call messages_fatal(1)
    case(4)
      write(message(1), '(a)') "QMR failure, can't continue: gamma is zero!"
      call messages_fatal(1)
    case(5)
      write(message(1), '(a)') "QMR Solver not converged!"
      call messages_warning(1)
      if(present(converged)) then
        converged = .false.
      end if
    end select
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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
