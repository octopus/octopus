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
!! $Id$

#include "global.h"

! This module is intended to contain "only mathematical" functions
! and procedures.

module solvers_m
  use global_m
  use lalg_basic_m
  use loct_math_m
  use messages_m
  use loct_m
  use blas_m
  use lalg_adv_m

  implicit none

  private
  public ::                     &
    dconjugate_gradients,       &
    zconjugate_gradients,       &
    zqmr_sym,                   &
    zqmr

  ! ---------------------------------------------------------
  ! QMR (quasi-minimal residual) algorithm for complex symmetric matrices
  ! algorithm taken from:
  ! Parallel implementation of efficient preconditioned linear solver for
  ! grid-based applications in chemical physics. II: QMR linear solver
  ! Appendix A. Simplified QMR algorithm
  interface zqmr_sym
    module procedure zqmr_sym_spec_dotp, zqmr_sym_gen_dotp
  end interface
  integer, pointer :: np_p

  ! ---------------------------------------------------------
  ! QMR (quasi-minimal residual) algorithm for complex matrices
  ! algorithm taken from: An Implementation of the QMR Method based on
  ! Coupled Two-Term Recurrences by R. W. Freund and N. M. Nachtigal (page 25)
  interface zqmr
    module procedure zqmr_spec_dotp, zqmr_gen_dotp
  end interface

  interface dconjugate_gradients
    module procedure dsym_conjugate_gradients, dbi_conjugate_gradients
  end interface

  interface zconjugate_gradients
    module procedure zsym_conjugate_gradients, zbi_conjugate_gradients
  end interface

contains

  ! ---------------------------------------------------------
  subroutine zqmr_sym_spec_dotp(np, x, b, op, prec, iter, residue, threshold, showprogress, converged)
    integer, target,   intent(in)    :: np    ! number of points
    CMPLX,             intent(inout) :: x(:)  ! initial guess and result
    CMPLX,             intent(in)    :: b(:)  ! the right side
    interface
      subroutine op(x, y)                     ! the matrix A as operator, y <- Ax
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine prec(x, y)                   ! preconditioner
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine prec
    end interface
    integer,           intent(inout) :: iter  ! [in] number of max iterations, [out] used iterations
    FLOAT, optional,   intent(out)   :: residue ! residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold ! convergence threshold
    logical, optional, intent(in)    :: showprogress ! show progress bar
    logical, optional, intent(out)   :: converged ! is the algorithm converged

    call push_sub('math.zqmr_sym_spec_dotp')

    np_p => np
    call zqmr_sym_gen_dotp(np, x, b, op, dotu_qmr, nrm2_qmr, prec, iter, &
      residue, threshold, showprogress, converged)

    call pop_sub()
  end subroutine zqmr_sym_spec_dotp

  ! ---------------------------------------------------------
  subroutine zqmr_spec_dotp(np, x, b, op, opt, prec, prect, iter, residue, threshold, showprogress, converged)
    integer, target,   intent(in)    :: np    ! number of points
    CMPLX,             intent(inout) :: x(:)  ! initial guess and the result
    CMPLX,             intent(in)    :: b(:)  ! the right side
    interface
      subroutine op(x, y)                     ! the matrix A as operator y <- A*x
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine opt(x, y)                    ! the transposed matrix A as operator y <- A^T*x
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine opt
    end interface
    interface
      subroutine prec(x, y)                   ! the preconditioner
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine prec
    end interface
    interface
      subroutine prect(x, y)                  ! the transposed preconditioner
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine prect
    end interface
    integer,           intent(inout) :: iter  ! [in] number of max iterations, [out] used iterations
    FLOAT, optional,   intent(out)   :: residue ! residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold ! convergence threshold
    logical, optional, intent(in)    :: showprogress ! show progress bar
    logical, optional, intent(out)   :: converged ! is the algorithm converged

    call push_sub('math.zqmr_spec_dotp')

    np_p => np
    call zqmr_gen_dotp(np, x, b, op, opt, dotu_qmr, nrm2_qmr, &
      prec, prect, iter, residue, threshold, showprogress, converged)

    call pop_sub()
  end subroutine zqmr_spec_dotp


  ! ---------------------------------------------------------
  ! the complex dot product without conjugated vector
  CMPLX function dotu_qmr(x, y)
    CMPLX, intent(in) :: x(:)
    CMPLX, intent(in) :: y(:)

    call push_sub('math.dotu_qmr')

    dotu_qmr = lalg_dotu(np_p, x, y)

    call pop_sub()
  end function dotu_qmr


  ! ---------------------------------------------------------
  FLOAT function nrm2_qmr(x)
    CMPLX, intent(in) :: x(:)

    call push_sub('math.nrm2_qmr')

    nrm2_qmr = lalg_nrm2(np_p, x(:))

    call pop_sub()
  end function nrm2_qmr


  ! ---------------------------------------------------------
  ! for complex symmetric matrices
  subroutine zqmr_sym_gen_dotp(np, x, b, op, dotu, nrm2, prec, iter, &
    residue, threshold, showprogress, converged)
    integer,           intent(in)    :: np    ! number of points
    CMPLX,             intent(inout) :: x(:)  ! the initial guess and the result
    CMPLX,             intent(in)    :: b(:)  ! the right side
    interface
      subroutine op(x, y)                     ! the matrix A as operator
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine op
    end interface
    interface 
      CMPLX function dotu(x, y)               ! the dot product (must be x^T*y, not daggered)
        CMPLX, intent(in) :: x(:)
        CMPLX, intent(in) :: y(:)
      end function dotu
    end interface
    interface 
      FLOAT function nrm2(x)                  ! the 2-norm of the vector x
        CMPLX, intent(in) :: x(:)
      end function nrm2
    end interface
    interface
      subroutine prec(x, y)                   ! the preconditioner
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine prec
    end interface
    integer,           intent(inout) :: iter  ! [in] the maximum number of iterations, [out] used iterations
    FLOAT, optional,   intent(out)   :: residue ! the residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold ! convergence threshold
    logical, optional, intent(in)    :: showprogress ! should there be a progress bar
    logical, optional, intent(out)   :: converged ! has the algorithm converged

    CMPLX, allocatable  :: r(:), v(:), z(:), q(:), p(:), deltax(:), deltar(:)
    CMPLX               :: eta, delta, epsilon, beta
    FLOAT               :: rho, xsi, gamma, alpha, theta, threshold_, res, oldtheta, oldgamma, oldrho
    integer             :: max_iter, err, ip
    logical             :: showprogress_
    FLOAT               :: log_res, log_thr
    integer             :: ilog_res, ilog_thr

    call push_sub('math.zqmr_sym_gen_dotp')

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

    ALLOCATE(r(np), np)
    ALLOCATE(v(np), np)
    ALLOCATE(z(np), np)
    ALLOCATE(q(np), np)
    ALLOCATE(p(np), np)
    ALLOCATE(deltax(np), np)
    ALLOCATE(deltar(np), np)

    ! use v as temp var
    call op(x, v)

    forall (ip = 1:np) 
      r(ip) = b(ip) - v(ip)
      v(ip) = r(ip)
    end forall

    rho = nrm2(v)

    max_iter = iter
    iter     = 0
    err      = 0
    res      = rho

    ! If rho is basically zero we are already done.
    if(abs(rho).gt.M_EPSILON) then
      call prec(v, z)

      xsi = nrm2(z)

      gamma = M_ONE
      eta   = -M_z1
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
          forall (ip = 1:np) q(ip) = -rho*delta/epsilon*q(ip) + z(ip)
        end if
        call op(q, p)
        call lalg_scal(np, alpha, p)

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

        if(iter == 1) then

          forall (ip = 1:np) 
            deltax(ip) = eta*alpha*q(ip)
            x(ip) = x(ip) + deltax(ip)
          end forall

          forall (ip = 1:np)
            deltar(ip) = eta*p(ip)
            r(ip) = r(ip) - deltar(ip)
          end forall

        else

          forall (ip = 1:np) 
            deltax(ip) = (oldtheta*gamma)**2*deltax(ip) + eta*alpha*q(ip)
            x(ip) = x(ip) + deltax(ip)
          end forall

          forall (ip = 1:np) 
            deltar(ip) = (oldtheta*gamma)**2*deltar(ip) + eta*p(ip)
            r(ip) = r(ip) - deltar(ip)
          end forall

        end if

        res = nrm2(r)/nrm2(x)

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
      call write_fatal(1)
    case(2)
      write(message(1), '(a)') "QMR failure, can't continue: v^T*z is zero!"
      call write_fatal(1)
    case(3)
      write(message(1), '(a)') "QMR failure, can't continue: q^T*p is zero!"
      call write_fatal(1)
    case(4)
      write(message(1), '(a)') "QMR failure, can't continue: gamma is zero!"
      call write_fatal(1)
    case(5)
      write(message(1), '(a)') "QMR Solver not converged!"
      call write_warning(1)
      if(present(converged)) then
        converged = .false.
      end if
    end select
    if(showprogress_) write(*,*) ''

    if(present(residue)) residue = res

    deallocate(r, v, z, q, p, deltax, deltar)

    call pop_sub()
  end subroutine zqmr_sym_gen_dotp

  ! ---------------------------------------------------------
  ! for general complex matrices
  ! taken from 'An Implementation of the QMR Method based on Coupled Two-Term Recurrences' by
  ! R. W. Freund and N. M. Nachtigal (page 25)
  subroutine zqmr_gen_dotp(np, x, b, op, opt, dotu, nrm2, prec, prect, iter, &
    residue, threshold, showprogress, converged)
    integer,           intent(in)    :: np    ! number of points
    CMPLX,             intent(inout) :: x(:)  ! initial guess and result
    CMPLX,             intent(in)    :: b(:)  ! right side
    interface
      subroutine op(x, y)                     ! the matrix A as operator: y <- A*x
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine op
    end interface
    interface
      subroutine opt(x, y)                    ! the transposed matrix A as operator: y <- A^T*x
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine opt
    end interface
    interface 
      CMPLX function dotu(x, y)               ! the dot product
        CMPLX, intent(in) :: x(:)
        CMPLX, intent(in) :: y(:)
      end function dotu
    end interface
    interface 
      FLOAT function nrm2(x)                  ! the 2-norm of a vector
        CMPLX, intent(in) :: x(:)
      end function nrm2
    end interface
    interface
      subroutine prec(x, y)                   ! preconditioner
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine prec
    end interface
    interface
      subroutine prect(x, y)                  ! transposed preconditioner
        CMPLX, intent(in)  :: x(:)
        CMPLX, intent(out) :: y(:)
      end subroutine prect
    end interface
    integer,           intent(inout) :: iter  ! [in] the maximum number of iterations, [out] used iterations
    FLOAT, optional,   intent(out)   :: residue ! the residue = abs(Ax-b)
    FLOAT, optional,   intent(in)    :: threshold ! convergence threshold
    logical, optional, intent(in)    :: showprogress ! should there be a progress bar
    logical, optional, intent(out)   :: converged ! has the algorithm converged

    CMPLX, allocatable  :: r(:), v(:), w(:), z(:), q(:), p(:), deltax(:), tmp(:)
    CMPLX               :: eta, delta, epsilon, beta
    FLOAT               :: rho, xsi, gamma, theta, threshold_, res, oldtheta, oldgamma, oldrho, nw
    integer             :: max_iter, err
    logical             :: showprogress_
    FLOAT               :: log_res, log_thr
    integer             :: ilog_res, ilog_thr

    call push_sub('math.zqmr_gen_dotp')

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

    ALLOCATE(r(np), np)
    ALLOCATE(v(np), np)
    ALLOCATE(w(np), np)
    ALLOCATE(z(np), np)
    ALLOCATE(q(np), np)
    ALLOCATE(p(np), np)
    ALLOCATE(tmp(np), np)
    ALLOCATE(deltax(np), np)

    ! r = b-Ax
    call lalg_copy(np, b, r)
    call op(x, tmp)
    call lalg_axpy(np, -M_z1, tmp, r)
    rho = nrm2(r)
    max_iter = iter
    iter     = 0
    err      = 0
    res      = rho

! If rho is basically zero we are already done.
    if(abs(rho).gt.M_EPSILON) then
      call lalg_copy(np, r, v)
      call lalg_scal(np, M_z1/rho, v)
      call lalg_copy(np, r, w)
      call prect(w, z)
      xsi = nrm2(z)
      call lalg_scal(np, M_z1/xsi, w)
      epsilon = M_z1
      gamma = M_ONE
      xsi = M_ONE
      theta = M_ZERO
      eta   = -M_z1

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
          call lalg_axpy(np, M_z1, z, p)
          call prect(w, z)
          call lalg_scal(np, -rho*delta/epsilon, q)
          call lalg_axpy(np, M_z1, z, q)
        end if

        call op(p,tmp)
        epsilon = dotu(q, tmp)
        if(abs(epsilon) < M_EPSILON) then
          err = 3
          exit
        end if
        beta = epsilon/delta
        call lalg_scal(np, -beta, v)
        call lalg_axpy(np, M_z1, tmp, v)
        oldrho = rho
        rho = nrm2(v)
        nw = nrm2(w)
        call lalg_scal(np, -beta, w)
        call opt(q,tmp)
        call lalg_axpy(np, M_z1, tmp, w)
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
          call lalg_axpy(np, M_z1, deltax, x)
        end if
        call lalg_scal(np, M_z1/rho, v)
        call lalg_scal(np, M_z1/xsi, w)

        call lalg_copy(np, b, r)
        call op(x, tmp)
        call lalg_axpy(np, -M_z1, tmp, r)
        res = nrm2(r)/nrm2(x)

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
      call write_fatal(1)
    case(2)
      write(message(1), '(a)') "QMR failure, can't continue: z^T*y is zero!"
      call write_fatal(1)
    case(3)
      write(message(1), '(a)') "QMR failure, can't continue: q^T*p is zero!"
      call write_fatal(1)
    case(4)
      write(message(1), '(a)') "QMR failure, can't continue: gamma is zero!"
      call write_fatal(1)
    case(5)
      write(message(1), '(a)') "QMR Solver not converged!"
      call write_warning(1)
      if(present(converged)) then
        converged = .false.
      end if
    end select
    if(showprogress_) write(*,*) ''

    if(present(residue)) residue = res

    deallocate(r, v, w, z, q, p, deltax, tmp)

    call pop_sub()
  end subroutine zqmr_gen_dotp

#include "undef.F90"
#include "complex.F90"
#include "solvers_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "solvers_inc.F90"

end module solvers_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
