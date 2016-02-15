!! Copyright (C) 2004 Xavier Andrade, M. Marques
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
!> This subroutine calculates the solution of (H + shift) x = y
!! Typically shift = - eigenvalue + omega
! ---------------------------------------------------------
subroutine X(linear_solver_solve_HXeY) (this, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used, occ_response)
  type(linear_solver_t), target, intent(inout) :: this
  type(hamiltonian_t),   target, intent(in)    :: hm
  type(grid_t),          target, intent(in)    :: gr
  type(states_t),        target, intent(in)    :: st
  integer,                       intent(in)    :: ist
  integer,                       intent(in)    :: ik
  R_TYPE,                        intent(inout) :: x(:,:)   !< x(gr%mesh%np_part, d%dim)
  R_TYPE,                        intent(in)    :: y(:,:)   !< y(gr%mesh%np, d%dim)
  R_TYPE,                        intent(in)    :: shift
  FLOAT,                         intent(in)    :: tol
  FLOAT,                         intent(out)   :: residue
  integer,                       intent(out)   :: iter_used
  logical, optional,             intent(in)    :: occ_response

  logical :: occ_response_
  R_TYPE, allocatable :: z(:, :)

  PUSH_SUB(X(linear_solver_solve_HXeY))
  call profiling_in(prof, "LINEAR_SOLVER")

  occ_response_ = .true.
  if(present(occ_response)) occ_response_ = occ_response

  args%ls       => this
  args%hm       => hm
  args%gr       => gr 
  args%st       => st
  args%ist      = ist
  args%ik       = ik
  args%X(shift) = shift
  iter_used = this%max_iter

  select case(this%solver)

  case(OPTION__LINEARSOLVER__CG)
    call X(linear_solver_cg)       (this, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used)

  case(OPTION__LINEARSOLVER__BICGSTAB)
    call X(linear_solver_bicgstab) (this, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used, occ_response_)

  case(OPTION__LINEARSOLVER__BICGSTAB2)
#if defined(R_TCOMPLEX)
    message(1) = 'The BICGSTAB(2) method (LinearSolver = bicgstab2) is not available for complex matrices.'
    call messages_fatal(1)
#else
  
    call dlinear_solver_bicgstab2 (this, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used, l = 2)
#endif

  case(OPTION__LINEARSOLVER__MULTIGRID)
    call X(linear_solver_multigrid)(this, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used)

  case(OPTION__LINEARSOLVER__QMR_SYMMETRIC)   
    ! complex symmetric: for Sternheimer, only if wfns are real
    call X(qmr_sym_gen_dotu)(gr%mesh%np, x(:, 1), y(:, 1), &
      X(linear_solver_operator_na), X(mf_dotu_aux), X(mf_nrm2_aux), X(linear_solver_preconditioner), &
      iter_used, residue = residue, threshold = tol, showprogress = .false.)

  case(OPTION__LINEARSOLVER__QMR_SYMMETRIZED)
    ! symmetrized equation
    SAFE_ALLOCATE(z(1:gr%mesh%np, 1:1))
    call X(linear_solver_operator_t_na)(y(:, 1), z(:, 1))
    call X(qmr_sym_gen_dotu)(gr%mesh%np, x(:, 1), z(:, 1), &
      X(linear_solver_operator_sym_na), X(mf_dotu_aux), X(mf_nrm2_aux), X(linear_solver_preconditioner), &
      iter_used, residue = residue, threshold = tol, showprogress = .false.)

  case(OPTION__LINEARSOLVER__QMR_DOTP)
    ! using conjugated dot product
    call X(qmr_sym_gen_dotu)(gr%mesh%np, x(:, 1), y(:, 1), &
      X(linear_solver_operator_na), X(mf_dotp_aux), X(mf_nrm2_aux), X(linear_solver_preconditioner), &
      iter_used, residue = residue, threshold = tol, showprogress = .false.)

  case(OPTION__LINEARSOLVER__QMR_GENERAL)
    ! general algorithm
    call X(qmr_gen_dotu)(gr%mesh%np, x(:, 1), y(:, 1), X(linear_solver_operator_na), X(linear_solver_operator_t_na), &
      X(mf_dotu_aux), X(mf_nrm2_aux), X(linear_solver_preconditioner), X(linear_solver_preconditioner), &
      iter_used, residue = residue, threshold = tol, showprogress = .false.)

  case(OPTION__LINEARSOLVER__SOS)
    call X(linear_solver_sos)(hm, gr, st, ist, ik, x, y, shift, residue, iter_used)

  case default 
    write(message(1), '(a,i2)') "Unknown linear-response solver", this%solver
    call messages_fatal(1)

  end select

  call profiling_out(prof)
  POP_SUB(X(linear_solver_solve_HXeY))

end subroutine X(linear_solver_solve_HXeY)

! ---------------------------------------------------------

subroutine X(linear_solver_solve_HXeY_batch) (this, hm, gr, st, ik, xb, yb, shift, tol, residue, iter_used, occ_response)
  type(linear_solver_t), target, intent(inout) :: this
  type(hamiltonian_t),   target, intent(in)    :: hm
  type(grid_t),          target, intent(in)    :: gr
  type(states_t),        target, intent(in)    :: st
  integer,                       intent(in)    :: ik
  type(batch_t),                 intent(inout) :: xb
  type(batch_t),                 intent(in)    :: yb
  R_TYPE,                        intent(in)    :: shift(:)
  FLOAT,                         intent(in)    :: tol
  FLOAT,                         intent(out)   :: residue(:)
  integer,                       intent(out)   :: iter_used(:)
  logical, optional,             intent(in)    :: occ_response

  integer :: ii

  PUSH_SUB(X(linear_solver_solve_HXeY_batch))

  select case(this%solver)
  case(OPTION__LINEARSOLVER__QMR_DOTP)
    call profiling_in(prof_batch, "LINEAR_SOLVER_BATCH")
    call X(linear_solver_qmr_dotp)(this, hm, gr, st, ik, xb, yb, shift, iter_used, residue, tol)
    call profiling_out(prof_batch)

  case default
    do ii = 1, xb%nst
      call X(linear_solver_solve_HXeY) (this, hm, gr, st, xb%states(ii)%ist, ik, xb%states(ii)%X(psi), yb%states(ii)%X(psi), &
        shift(ii), tol, residue(ii), iter_used(ii), occ_response)
    end do

  end select

  POP_SUB(X(linear_solver_solve_HXeY_batch))

end subroutine X(linear_solver_solve_HXeY_batch)

! ---------------------------------------------------------
!> Conjugate gradients
subroutine X(linear_solver_cg) (ls, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used)
  type(linear_solver_t), intent(inout) :: ls
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(in)    :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:,:)   !< x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(in)    :: y(:,:)   !< y(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(in)    :: shift
  FLOAT,                 intent(in)    :: tol
  FLOAT,                 intent(out)   :: residue
  integer,               intent(out)   :: iter_used

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter, idim
  logical :: conv_last, conv

  PUSH_SUB(X(linear_solver_cg))

  SAFE_ALLOCATE( r(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE( p(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(Hp(1:gr%mesh%np, 1:st%d%dim))

  ! Initial residue
  call X(linear_solver_operator)(hm, gr, st, ist, ik, shift, x, Hp)
  r(1:gr%mesh%np, 1:st%d%dim) = y(1:gr%mesh%np, 1:st%d%dim) - Hp(1:gr%mesh%np, 1:st%d%dim)
  
  ! Initial search direction
  p(1:gr%mesh%np, 1:st%d%dim) = r(1:gr%mesh%np, 1:st%d%dim)
  p((gr%mesh%np+1):gr%mesh%np_part,1:st%d%dim) = M_ZERO
  
  conv_last = .false.
  gamma     = M_ONE
  do iter = 1, ls%max_iter
    gamma = X(mf_dotp)(gr%mesh, st%d%dim, r, r)

    conv = ( abs(gamma) < tol**2)
    if(conv.and.conv_last) exit
    conv_last = conv
    
    call X(linear_solver_operator)(hm, gr, st, ist, ik, shift, p, Hp)

    alpha = gamma/X(mf_dotp) (gr%mesh, st%d%dim, p, Hp)

    do idim = 1, st%d%dim
      !r = r - alpha*Hp
      call lalg_axpy(gr%mesh%np, -alpha, Hp(:, idim), r(:, idim))
      !x = x + alpha*p
      call lalg_axpy(gr%mesh%np,  alpha,  p(:, idim), x(:, idim))
    end do


    beta = X(mf_dotp)(gr%mesh, st%d%dim, r, r)/gamma

    p(1:gr%mesh%np, 1:st%d%dim) = r(1:gr%mesh%np, 1:st%d%dim) + beta*p(1:gr%mesh%np, 1:st%d%dim)

  end do
    
  iter_used = iter
  residue = sqrt(abs(gamma))

  if(.not. conv) then 
    write(message(1), '(a)') "CG solver not converged!"
    call messages_warning(1)
  end if

  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(p)
  SAFE_DEALLOCATE_A(Hp)

  POP_SUB(X(linear_solver_cg))
end subroutine X(linear_solver_cg)



#if defined(R_TREAL)
! ---------------------------------------------------------
!> BICONJUGATE GRADIENTS STABILIZED(2)
!>
!> This is the BICGSTAB(L) described in:
!> 
!> * G.L.G. Sleijpen and D.R. Fokkema, BiCGstab(l) for linear equations involving unsymmetric 
!>   matrices with complex spectrum, Electronic Transactions on Numer. Anal. (ETNA) 1, 11 (1993)
!>
!> Only for l=1, and l=2. For l=1 it is not fundamentally different from the bicgstab, although
!> it seems to break down less. This version contains two enhancements, described in:
!>
!> * G.Sleijpen and H.van der Vorst, Maintaining convergence properties of BiCGstab methods in finite 
!>   precision arithmetic, Numerical Algorithms 10, 203 (1995).
!> * G.Sleijpen and H.van der Vorst, Reliable updated residuals in hybrid BiCG methods, 
!>   Computing 56, 141 (1996).
!> 
!> This version is based on the bc2g.f routine of M. A. Botchev:
!> [http://www.staff.science.uu.nl/~vorst102/software.html]
!>
subroutine dlinear_solver_bicgstab2 (ls, hm, gr, st, ist, ik, x, rhsy, shift, tol, residue, iter_used, l)
  type(linear_solver_t), intent(inout) :: ls
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(in)    :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:, :)    !< x(gr%mesh%np_part, st%d%dim)
  R_TYPE,                intent(in)    :: rhsy(:, :) !< rhsy(gr%mesh%np, st%d%dim)
  FLOAT,                 intent(in)    :: shift
  FLOAT,                 intent(in)    :: tol
  FLOAT,                 intent(out)   :: residue
  integer,               intent(out)   :: iter_used
  integer,               intent(in)       :: l

  integer :: ii, i1, jj, kk, z, zz, y0, yl, y, rr, r, u, xp, bp, n, idim
  FLOAT, allocatable :: psi(:, :), zeta(:, :), xx(:), rhs(:)
  logical :: rcmp, xpdt
  FLOAT :: alpha, beta, hatgamma, kappa0, kappal, mxnrmr, mxnrmx, omega, rho0, rho1, residue0, &
           sigma, varrho, delta

  PUSH_SUB(dlinear_solver_bicgstab2)

  n = gr%mesh%np * st%d%dim
  SAFE_ALLOCATE(psi(1:n, 1:2*l+5))
  SAFE_ALLOCATE(zeta(1:l+1, 1:3+2*(l+1)))
  SAFE_ALLOCATE(xx(1:n))
  SAFE_ALLOCATE(rhs(1:n))
  delta = CNST(1.0e-2)

  do idim = 1, st%d%dim
    xx((idim-1)*gr%mesh%np+1 : st%d%dim*gr%mesh%np) = x(1:gr%mesh%np, idim)
    rhs((idim-1)*gr%mesh%np+1 : st%d%dim*gr%mesh%np) = rhsy(1:gr%mesh%np, idim)
  end do

  rr = 1
  r = rr + 1
  u = r + (l+1)
  xp = u + (l+1)
  bp = xp + 1

  z = 1
  zz = z + (l+1)
  y0 = zz + (l+1)
  yl = y0 + 1
  y = yl + 1

  ! Calculation of the initial residual.
  call op (xx, psi(:, r) )
  psi(:, r) = rhs(:) - psi(:, r)
  iter_used = 1

  psi(:, rr) = psi(:, r)
  psi(:, bp) = psi(:, r)
  psi(:, xp) = xx(:)
  xx(:) = M_ZERO

  residue0 = sqrt(dotp(psi(:, r), psi(:, r)))
  residue = residue0

  mxnrmx = residue0
  mxnrmr = residue0
  rcmp = .false.
  xpdt = .false.

  alpha = M_ZERO
  omega = M_ONE
  sigma = M_ONE
  rho0 =  M_ONE

  main_iteration : do while ( (residue > tol) .and.  (iter_used < ls%max_iter) )

    ! BICG
    rho0 = -omega * rho0
    do kk = 1, l
      rho1 = dotp(psi(:, rr), psi(:, r+kk-1))

      if (rho0 .eq. M_ZERO) then
        write(message(1), '(a)') "BiCGSTAB solver not converged!"
        call messages_warning(1)
        exit main_iteration
      end if
      beta = alpha*(rho1/rho0)
      rho0 = rho1
      do jj = 0, kk-1
        psi(:, u+jj) = psi(:, r+jj) - beta * psi(:, u+jj)
      end do

      call op(psi(:, u+kk-1), psi(:, u+kk))
      iter_used = iter_used + 1
      sigma = dotp(psi(:, rr), psi(:, u+kk))

      if (sigma .eq. M_ZERO) then
        write(message(1), '(a)') "BiCGSTAB solver not converged!"
        call messages_warning(1)
        exit main_iteration
      end if

      alpha = rho1/sigma
      xx(:) = alpha * psi(:, u) + xx(:)

      do jj = 0, kk-1
        psi(:, r+jj) = - alpha * psi(:, u+jj+1) + psi(:, r+jj)
      end do

      call op (psi(:, r+kk-1), psi(:, r+kk))
      iter_used = iter_used + 1

      residue = sqrt(dotp(psi(:, r), psi(:, r)))
      mxnrmx = max (mxnrmx, residue)
      mxnrmr = max (mxnrmr, residue)

    end do

    ! Calculation of the matrix Z = R^* R
    do i1 = 1, l+1
      do jj = i1-1, l
        zeta(jj+1, z+i1-1) = dotp(psi(:, r+jj), psi(:, r+i1-1))
        zeta(z+i1-1, jj+1) = zeta(jj+1, z+i1-1)
      end do
    end do

    do i1 = zz, zz+l
      do ii = 1, l+1
        zeta(ii, i1)   = zeta(ii, i1+(z-zz))
      end do
    end do

    zeta(1, y0) = - M_ONE
    zeta(2, y0) = zeta(2, z) / zeta(2, zz+1)
    zeta(l+1, y0) = M_ZERO

    zeta(1, yl) = M_ZERO
    zeta(2, yl) = zeta(2, z+l) / zeta(2, zz+1)
    zeta(l+1, yl) = - M_ONE

    ! Convex combination
    do ii = 1, l+1
      zeta(ii, y) = M_ZERO
    end do
    do jj = 1, l+1
      do ii = 1, l+1
        zeta(ii, y) = zeta(ii, y) + zeta(jj, yl) * zeta(ii, z+jj-1)
      end do
    end do
    kappal = sqrt( dot_product(zeta(1:l+1, yl), zeta(1:l+1, y)) ) 

    do ii = 1, l+1
      zeta(ii, y) = M_ZERO
    end do
    do jj = 1, l+1
      do ii = 1, l+1
        zeta(ii, y) = zeta(ii, y) + zeta(jj, y0) * zeta(ii, z+jj-1)
      end do
    end do
    kappa0 = sqrt(dot_product(zeta(1:l+1, y0), zeta(1:l+1, y)))

    varrho = dot_product(zeta(1:l+1, yl), zeta(1:l+1, y))
    varrho = varrho / (kappa0*kappal)

    hatgamma = sign( M_ONE, varrho) * max(abs(varrho),CNST(0.7)) * (kappa0/kappal)

    do ii=1,l+1
      zeta(ii, y0) = -hatgamma*zeta(ii, yl) + zeta(ii, y0)
    end do

    !   Update
    omega = zeta(l+1,y0)
    do jj=1,l
        psi(:, u) = psi(:, u) - zeta(1+jj, y0) * psi(:, u+jj)
        xx(:)     = xx(:)     + zeta(1+jj, y0) * psi(:, r+jj-1)
        psi(:,r) = psi(:, r) - zeta(1+jj, y0) * psi(:, r+jj)
    end do

    do ii = 1, l+1
      zeta(ii, y) = M_ZERO 
    end do
    do jj = 1, l+1
      do ii = 1, l+1
        zeta(ii, y) = zeta(ii, y) + zeta(jj, y0) * zeta(ii, z+jj-1)
      end do
    end do

    residue = sqrt(dot_product(zeta(1:l+1, y0), zeta(1:l+1, y)))

    ! Reliable update.
    mxnrmx = max (mxnrmx, residue)
    mxnrmr = max (mxnrmr, residue)
    xpdt = (residue < delta*residue0) .and. (residue0 < mxnrmx)
    rcmp = ( ( (residue < delta*mxnrmr) .and. (residue0 < mxnrmr) ) .or. xpdt)
    if (rcmp) then
      call op (xx, psi(:, r) )
      iter_used = iter_used + 1
      psi(:, r) = psi(:, bp) - psi(:, r)
      mxnrmr = residue
      if (xpdt) then
        psi(:, xp) = xx(:) + psi(:, xp)
        xx(:) = M_ZERO
        psi(:, bp) = psi(:, r)
        mxnrmx = residue
      end if
    end if

  end do main_iteration

  ! End of iterations
  xx(:) = psi(:, xp) + xx(:)

  ! Computation of the "true" residual
  ! call op (xx, psi(:, r) )
  ! psi(:, r) = rhs(:) - psi(:, r)
  ! residue = sqrt( dotp(psi(:, r), psi(:, r)) )

  if (residue .gt. tol) then 
    write(message(1), '(a)') "BiCGSTAB solver not converged!"
    call messages_warning(1)
  end if

  do idim = 1, st%d%dim
    x(1:gr%mesh%np, idim) = xx((idim-1)*gr%mesh%np+1 : st%d%dim*gr%mesh%np)
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(zeta)
  SAFE_DEALLOCATE_A(xx)
  SAFE_DEALLOCATE_A(rhs)
  POP_SUB(dlinear_solver_bicgstab2)
  return

  contains

    FLOAT function dotp(x, y)
      FLOAT, intent(in) :: x(:)
      FLOAT, intent(in) :: y(:)
      integer :: idim
      dotp = M_ZERO
      do idim = 1, st%d%dim
        dotp = dotp + &
          dmf_dotp(gr%mesh, x((idim-1)*gr%mesh%np+1:st%d%dim*gr%mesh%np), y((idim-1)*gr%mesh%np+1:st%d%dim*gr%mesh%np))
      end do
    end function dotp

    subroutine op(x, y)
      FLOAT, intent(in) :: x(:)
      FLOAT, intent(out) :: y(:)

      integer :: idim
      FLOAT, allocatable :: a(:, :), opa(:, :)

      SAFE_ALLOCATE(a(gr%mesh%np_part, st%d%dim))
      SAFE_ALLOCATE(opa(gr%mesh%np, st%d%dim))

      a = M_ZERO
      do idim = 1, st%d%dim
        a(1:gr%mesh%np, idim) = x((idim-1)*gr%mesh%np+1 : st%d%dim*gr%mesh%np)
      end do
      call dlinear_solver_operator (hm, gr, st, ist, ik, shift, a, opa)
      do idim = 1, st%d%dim
        y((idim-1)*gr%mesh%np+1 : st%d%dim*gr%mesh%np) = opa(1:gr%mesh%np, idim)
      end do

      SAFE_DEALLOCATE_A(a)
      SAFE_DEALLOCATE_A(opa)
    end subroutine op

end subroutine dlinear_solver_bicgstab2
#endif


! ---------------------------------------------------------
!> BICONJUGATE GRADIENTS STABILIZED
!! see http://math.nist.gov/iml++/bicgstab.h.txt
subroutine X(linear_solver_bicgstab) (ls, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used, occ_response)
  type(linear_solver_t), intent(inout) :: ls
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(in)    :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:,:)   !< x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(in)    :: y(:,:)   !< y(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(in)    :: shift
  FLOAT,                 intent(in)    :: tol
  FLOAT,                 intent(out)   :: residue
  integer,               intent(out)   :: iter_used
  logical,               intent(in)    :: occ_response

  R_TYPE, allocatable :: r(:,:), Hp(:,:), rs(:,:), Hs(:,:), p(:,:), s(:,:), psi(:, :)
  R_TYPE, pointer :: phat(:,:), shat(:,:)
  R_TYPE  :: alpha, beta, w, rho_1, rho_2
  logical :: conv_last, conv
  integer :: iter, idim, ip
  FLOAT :: gamma
  
  PUSH_SUB(X(linear_solver_bicgstab))

  SAFE_ALLOCATE( r(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE( p(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(rs(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE( s(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(Hp(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(Hs(1:gr%mesh%np, 1:st%d%dim))

  ! this will store the preconditioned functions
  SAFE_ALLOCATE(phat(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(shat(1:gr%mesh%np_part, 1:st%d%dim))

  ! Initial residue
  call X(linear_solver_operator) (hm, gr, st, ist, ik, shift, x, Hp)

  forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np) r(ip, idim) = y(ip, idim) - Hp(ip, idim)

  !re-orthogonalize r, this helps considerably with convergence
  if (occ_response) then
    SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))

    call states_get_state(st, gr%mesh, ist, ik, psi)
    
    alpha = X(mf_dotp)(gr%mesh, st%d%dim, psi, r)
    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, -alpha, psi(:, idim), r(:, idim))
    end do

    SAFE_DEALLOCATE_A(psi)
  else
    ! project RHS onto the unoccupied states
    call X(lr_orth_vector)(gr%mesh, st, r, ist, ik, shift + st%eigenval(ist, ik))
  end if
          
  do idim = 1, st%d%dim
    call lalg_copy(gr%mesh%np, r(:, idim), rs(:, idim))
  end do

  gamma = X(mf_nrm2)(gr%mesh, st%d%dim, r)

  conv_last = .false.
  do iter = 1, ls%max_iter

    rho_1 = X(mf_dotp) (gr%mesh, st%d%dim, rs, r)

    if( abs(rho_1) < M_EPSILON ) exit

    if( iter == 1 ) then
      do idim = 1, st%d%dim
        call lalg_copy(gr%mesh%np, r(:, idim), p(:, idim))
      end do
    else
      beta = rho_1/rho_2*alpha/w
      forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np) p(ip, idim) = r(ip, idim) + beta*(p(ip, idim) - w*Hp(ip, idim))
    end if

    ! preconditioning 
    call X(preconditioner_apply)(ls%pre, gr, hm, ik, p, phat, shift)
    call X(linear_solver_operator)(hm, gr, st, ist, ik, shift, phat, Hp)
    
    alpha = rho_1/X(mf_dotp)(gr%mesh, st%d%dim, rs, Hp)

    forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np) s(ip, idim) = r(ip, idim) - alpha*Hp(ip, idim)

    gamma = X(mf_nrm2) (gr%mesh, st%d%dim, s)

    conv = (gamma < tol)
    if( conv ) then
      do idim = 1, st%d%dim 
        call lalg_axpy(gr%mesh%np, alpha, phat(:, idim), x(:, idim))
      end do
      exit
    end if

    call X(preconditioner_apply)(ls%pre, gr, hm, ik, s, shat, shift)
    call X(linear_solver_operator)(hm, gr, st, ist, ik, shift, shat, Hs)

    w = X(mf_dotp)(gr%mesh, st%d%dim, Hs, s)/X(mf_dotp) (gr%mesh, st%d%dim, Hs, Hs)

    forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np)
      x(ip, idim) = x(ip, idim) + alpha*phat(ip, idim) + w*shat(ip, idim)
      r(ip, idim) = s(ip, idim) - w*Hs(ip, idim)
    end forall

    rho_2 = rho_1

    gamma = X(mf_nrm2)(gr%mesh, st%d%dim, r)
    conv = (gamma < tol)
    if( conv .and. conv_last ) then 
      exit
    end if
    conv_last = conv

    if( abs(w) < M_EPSILON ) exit

  end do

  iter_used = iter
  residue = gamma

  if(.not. conv) then 
    write(message(1), '(a)') "BiCGSTAB solver not converged!"
    call messages_warning(1)
  end if

  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(p)
  SAFE_DEALLOCATE_A(Hp)
  SAFE_DEALLOCATE_A(s)
  SAFE_DEALLOCATE_A(rs)
  SAFE_DEALLOCATE_A(Hs)
  SAFE_DEALLOCATE_P(phat)
  SAFE_DEALLOCATE_P(shat)

  POP_SUB(X(linear_solver_bicgstab))
end subroutine X(linear_solver_bicgstab)


! ---------------------------------------------------------
subroutine X(linear_solver_multigrid) (ls, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used)
  type(linear_solver_t), intent(inout) :: ls
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(in)    :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:,:)   ! x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(in)    :: y(:,:)   ! y(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(in)    :: shift
  FLOAT,                 intent(in)    :: tol
  FLOAT,                 intent(out)   :: residue
  integer,               intent(out)   :: iter_used

  R_TYPE, allocatable :: diag(:,:), hx(:,:), res(:,:), psi(:, :)
  integer :: iter

  PUSH_SUB(X(linear_solver_multigrid))

  SAFE_ALLOCATE(diag(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(  hx(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE( res(1:gr%mesh%np, 1:st%d%dim))

  call X(hamiltonian_diagonal)(hm, gr%der, diag, ik)
  diag(1:gr%mesh%np, 1:st%d%dim) = diag(1:gr%mesh%np, 1:st%d%dim) + shift

  do iter = 1, ls%max_iter

    call smoothing(3)

    call smoothing(3)

    !calculate the residue
    call X(linear_solver_operator)(hm, gr, st, ist, ik, shift, x, hx)
    res(1:gr%mesh%np, 1:st%d%dim) = hx(1:gr%mesh%np, 1:st%d%dim) - y(1:gr%mesh%np, 1:st%d%dim)
    residue = X(mf_nrm2)(gr%mesh, st%d%dim, res)

    if(residue < tol) exit

    if(debug%info) then

      SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))
      
      call states_get_state(st, gr%mesh, ist, ik, psi)
      write(message(1), *)  "Multigrid: iter ", iter,  residue, abs(X(mf_dotp)(gr%mesh, st%d%dim, psi, x))
      call messages_info(1)

      SAFE_DEALLOCATE_A(psi)
      
    end if

  end do

  iter_used = iter

  if(residue > tol) then 
    write(message(1), '(a)') "Multigrid solver not converged!"
    call messages_warning(1)
  end if

  POP_SUB(X(linear_solver_multigrid))

contains 

  subroutine smoothing(steps)
    integer, intent(in) :: steps

    integer :: ii, ip, idim
    R_TYPE  :: rr

    PUSH_SUB(X(linear_solver_multigrid).smoothing)

    do ii = 1, steps

      call X(linear_solver_operator)(hm, gr, st, ist, ik, shift, x, hx)

      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          rr = hx(ip, idim) - y(ip, idim)
          x(ip, idim) = x(ip, idim) - M_TWOTHIRD * rr / diag(ip, idim)
        end do
      end do

    end do

    call X(lr_orth_vector)(gr%mesh, st, x, ist, ik, shift + st%eigenval(ist, ik))

    POP_SUB(X(linear_solver_multigrid).smoothing)
  end subroutine smoothing

end subroutine X(linear_solver_multigrid)


! ---------------------------------------------------------
!> This routine applies the operator hx = [H (+ Q) + shift] x
subroutine X(linear_solver_operator) (hm, gr, st, ist, ik, shift, x, hx)
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(in)    :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:,:)   !<  x(gr%mesh%np_part, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:,:)  !< Hx(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(in)    :: shift

  integer :: idim, jst
  FLOAT   :: alpha_j
  R_TYPE  :: proj
  R_TYPE, allocatable :: psi(:, :)

  PUSH_SUB(X(linear_solver_operator))

  call X(hamiltonian_apply)(hm, gr%der, x, Hx, ist, ik)

  !Hx = Hx + shift*x
  do idim = 1, st%d%dim
    call lalg_axpy(gr%mesh%np, shift, x(:, idim), Hx(:, idim))
  end do

  if(st%smear%method == SMEAR_SEMICONDUCTOR .or. st%smear%integral_occs) then
    POP_SUB(X(linear_solver_operator))
    return
  end if

  ! This is the Q term in Eq. (11) of PRB 51, 6773 (1995)
  ASSERT(.not. st%parallel_in_states)
  do jst = 1, st%nst
    alpha_j = lr_alpha_j(st, jst, ik)
    if(alpha_j == M_ZERO) cycle

    SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))

    call states_get_state(st, gr%mesh, jst, ik, psi)
    
    proj = X(mf_dotp)(gr%mesh, st%d%dim, psi, x)
    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, alpha_j*proj, psi(:, idim), Hx(:, idim))
    end do

    SAFE_DEALLOCATE_A(psi)
    
  end do

  POP_SUB(X(linear_solver_operator))

end subroutine X(linear_solver_operator)

! ---------------------------------------------------------
subroutine X(linear_solver_operator_batch) (hm, gr, st, ik, shift, xb, hxb)
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(in)    :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ik
  R_TYPE,                intent(in)    :: shift(:)
  type(batch_t),         intent(inout) :: xb   
  type(batch_t),         intent(inout) :: hxb  

  integer :: ii
  R_TYPE, allocatable :: shift_ist_indexed(:)

  PUSH_SUB(X(linear_solver_operator_batch))

  if(st%smear%method == SMEAR_SEMICONDUCTOR .or. st%smear%integral_occs) then

    call X(hamiltonian_apply_batch)(hm, gr%der, xb, hxb, ik)
    
    SAFE_ALLOCATE(shift_ist_indexed(st%st_start:st%st_end))
    
    do ii = 1, xb%nst 
      shift_ist_indexed(xb%states(ii)%ist) = shift(ii)
    end do
    
    call batch_axpy(gr%mesh%np, shift_ist_indexed, xb, hxb)
    
    SAFE_DEALLOCATE_A(shift_ist_indexed)

  else

    do ii = 1, xb%nst
      call X(linear_solver_operator)(hm, gr, st, xb%states(ii)%ist, ik, shift(ii), &
        xb%states(ii)%X(psi), hxb%states(ii)%X(psi))
    end do

  end if
    
  POP_SUB(X(linear_solver_operator_batch))

end subroutine X(linear_solver_operator_batch)

! ---------------------------------------------------------
!> applies linear_solver_operator with other arguments implicit as global variables
subroutine X(linear_solver_operator_na) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !<  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:)  !< Hx(gr%mesh%np, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)

  SAFE_ALLOCATE(tmpx(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpy(1:args%gr%mesh%np, 1:1))

  call lalg_copy(args%gr%mesh%np, x, tmpx(:, 1))
  call X(linear_solver_operator)(args%hm, args%gr, args%st, args%ist, args%ik, args%X(shift), tmpx, tmpy)
  call lalg_copy(args%gr%mesh%np, tmpy(:, 1), hx)

  SAFE_DEALLOCATE_A(tmpx)
  SAFE_DEALLOCATE_A(tmpy)

end subroutine X(linear_solver_operator_na)


! ---------------------------------------------------------
!> applies transpose of linear_solver_operator with other arguments implicit as global variables
!! \f$ (H - shift)^T = H* - shift = (H - shift*)* \f$ 
subroutine X(linear_solver_operator_t_na) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:)  ! Hx(gr%mesh%np, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)

  SAFE_ALLOCATE(tmpx(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpy(1:args%gr%mesh%np, 1:1))

  call lalg_copy(args%gr%mesh%np, R_CONJ(x), tmpx(:, 1))
  call X(linear_solver_operator)(args%hm, args%gr, args%st, args%ist, args%ik, R_CONJ(args%X(shift)), tmpx, tmpy)
  call lalg_copy(args%gr%mesh%np, R_CONJ(tmpy(:, 1)), hx)

  SAFE_DEALLOCATE_A(tmpx)
  SAFE_DEALLOCATE_A(tmpy)

end subroutine X(linear_solver_operator_t_na)


! ---------------------------------------------------------
!> applies linear_solver_operator in symmetrized form: \f$  A^T A \f$ 
subroutine X(linear_solver_operator_sym_na) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !<  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:)  !< Hx(gr%mesh%np, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)
  R_TYPE, allocatable :: tmpz(:, :)

  SAFE_ALLOCATE(tmpx(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpy(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpz(1:args%gr%mesh%np_part, 1:1))

  call lalg_copy(args%gr%mesh%np, x, tmpx(:, 1))
  call X(linear_solver_operator)(args%hm, args%gr, args%st, args%ist, args%ik, args%X(shift), tmpx, tmpy)
  call X(linear_solver_operator_t_na)(tmpy(:, 1), tmpz(:, 1))
  call lalg_copy(args%gr%mesh%np, tmpz(:, 1), hx)

  SAFE_DEALLOCATE_A(tmpx)
  SAFE_DEALLOCATE_A(tmpy)
  SAFE_DEALLOCATE_A(tmpz)

end subroutine X(linear_solver_operator_sym_na)

! ---------------------------------------------------------
subroutine X(linear_solver_preconditioner) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !<  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: hx(:)  !< Hx(gr%mesh%np, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)

  PUSH_SUB(X(linear_solver_preconditioner))

  SAFE_ALLOCATE(tmpx(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpy(1:args%gr%mesh%np_part, 1:1))

  call lalg_copy(args%gr%mesh%np, x, tmpx(:, 1))
  call X(preconditioner_apply)(args%ls%pre, args%gr, args%hm, args%ik, tmpx, tmpy, args%X(shift))
  call lalg_copy(args%gr%mesh%np, tmpy(:, 1), hx)

  SAFE_DEALLOCATE_A(tmpx)
  SAFE_DEALLOCATE_A(tmpy)
  POP_SUB(X(linear_solver_preconditioner))

end subroutine X(linear_solver_preconditioner)

! ---------------------------------------------------------
subroutine X(linear_solver_sos) (hm, gr, st, ist, ik, x, y, shift, residue, iter_used)
  type(hamiltonian_t),            intent(in)    :: hm
  type(grid_t),                   intent(in)    :: gr
  type(states_t),                 intent(in)    :: st
  integer,                        intent(in)    :: ist
  integer,                        intent(in)    :: ik
  R_TYPE,                         intent(inout) :: x(:,:)   !< x(gr%mesh%np, st%d%dim)
  R_TYPE,                         intent(in)    :: y(:,:)   !< y(gr%mesh%np, st%d%dim)
  R_TYPE,                         intent(in)    :: shift
  FLOAT,                          intent(out)   :: residue
  integer,                        intent(out)   :: iter_used

  integer :: jst, idim
  R_TYPE  :: aa
  R_TYPE, allocatable  :: rr(:, :)
  R_TYPE, allocatable :: psi(:, :)
  
  PUSH_SUB(X(linear_solver_sos))

  x(1:gr%mesh%np, 1:st%d%dim) = M_ZERO

  SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))
  
  do jst = 1, st%nst
    if(ist == jst) cycle

    call states_get_state(st, gr%mesh, jst, ik, psi)
    
    aa = X(mf_dotp)(gr%mesh, st%d%dim, psi, y)
    aa = aa/(st%eigenval(jst, ik) + lr_alpha_j(st, jst, ik) + shift)
    ! Normally the expression in perturbation theory would have here
    ! denominator = st%eigenval(jst, ik) - st%eigenval(ist, ik)
    ! For solving this type of problem, -st%eigenval(ist, ik) is included in shift

    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, aa, psi(:, idim), x(:, idim))
    end do
  end do

  SAFE_DEALLOCATE_A(psi)

  ! calculate the residual
  SAFE_ALLOCATE(rr(1:gr%mesh%np, 1:st%d%dim))
  call X(linear_solver_operator)(hm, gr, st, ist, ik, shift, x, rr)

  do idim = 1, st%d%dim
    call lalg_axpy(gr%mesh%np, -M_ONE, y(:, idim), rr(:, idim))
  end do
  
  residue = X(mf_nrm2)(gr%mesh, st%d%dim, rr)
  iter_used = 1

  SAFE_DEALLOCATE_A(rr)
  POP_SUB(X(linear_solver_sos))

end subroutine X(linear_solver_sos)

! ---------------------------------------------------------
!> for complex symmetric matrices
!! W Chen and B Poirier, J Comput Phys 219, 198-209 (2006)
subroutine X(linear_solver_qmr_dotp)(this, hm, gr, st, ik, xb, bb, shift, iter_used, residue, threshold)
  type(linear_solver_t), intent(inout) :: this
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(in)    :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ik
  type(batch_t),         intent(inout) :: xb
  type(batch_t),         intent(in)    :: bb
  R_TYPE,                intent(in)    :: shift(:)
  integer,               intent(out)   :: iter_used(:) 
  FLOAT,                 intent(out)   :: residue(:)   !< the residue = abs(Ax-b)
  FLOAT,                 intent(in)    :: threshold    !< convergence threshold

  type(batch_t) :: vvb, res, zzb, qqb, ppb, deltax, deltar
  FLOAT               :: oldgamma
  integer             :: ii, iter
  FLOAT, allocatable  :: rho(:), oldrho(:), norm_b(:), xsi(:), gamma(:), alpha(:), theta(:), oldtheta(:), saved_res(:)
  R_TYPE, allocatable :: eta(:), beta(:), delta(:), eps(:), exception_saved(:, :, :)
  integer, allocatable :: status(:), saved_iter(:)

  integer, parameter ::        &
    QMR_NOT_CONVERGED    = 0,  &
    QMR_CONVERGED        = 1,  &
    QMR_RES_ZERO         = 2,  &
    QMR_B_ZERO           = 3,  &
    QMR_BREAKDOWN_PB     = 4,  &
    QMR_BREAKDOWN_VZ     = 5,  &
    QMR_BREAKDOWN_QP     = 6,  &
    QMR_BREAKDOWN_GAMMA  = 7

  PUSH_SUB(X(linear_solver_qmr_dotp))

  SAFE_ALLOCATE(rho(1:xb%nst))
  SAFE_ALLOCATE(oldrho(1:xb%nst))
  SAFE_ALLOCATE(norm_b(1:xb%nst))
  SAFE_ALLOCATE(xsi(1:xb%nst))
  SAFE_ALLOCATE(gamma(1:xb%nst))
  SAFE_ALLOCATE(alpha(1:xb%nst))
  SAFE_ALLOCATE(eta(1:xb%nst))
  SAFE_ALLOCATE(theta(1:xb%nst))
  SAFE_ALLOCATE(oldtheta(1:xb%nst))
  SAFE_ALLOCATE(beta(1:xb%nst))
  SAFE_ALLOCATE(delta(1:xb%nst))
  SAFE_ALLOCATE(eps(1:xb%nst))
  SAFE_ALLOCATE(saved_res(1:xb%nst))

  SAFE_ALLOCATE(status(1:xb%nst))
  SAFE_ALLOCATE(saved_iter(1:xb%nst))

  SAFE_ALLOCATE(exception_saved(1:gr%mesh%np, 1:st%d%dim, 1:xb%nst))

  call batch_copy(xb, vvb)
  call batch_copy(xb, res)
  call batch_copy(xb, zzb)
  call batch_copy(xb, qqb)
  call batch_copy(xb, ppb)
  call batch_copy(xb, deltax)
  call batch_copy(xb, deltar)

  call X(linear_solver_operator_batch)(hm, gr, st, ik, shift, xb, vvb)

  call batch_xpay(gr%mesh%np, bb, CNST(-1.0), vvb)
  call batch_copy_data(gr%mesh%np, vvb, res)

  call mesh_batch_nrm2(gr%mesh, vvb, rho)
  call mesh_batch_nrm2(gr%mesh, bb, norm_b)

  status = QMR_NOT_CONVERGED

  iter = 0

  do ii = 1, xb%nst

    residue(ii) = rho(ii)

    ! If rho(ii) is basically zero we are already done.
    if(abs(rho(ii)) <= M_EPSILON) then
      status(ii) = QMR_RES_ZERO
      residue(ii) = rho(ii)
      call batch_get_state(xb, ii, gr%mesh%np, exception_saved(:, :, ii))
      saved_iter(ii) = iter
      saved_res(ii) = residue(ii)
    end if

    ! if b is zero, the solution is trivial
    if(status(ii) == QMR_NOT_CONVERGED .and. abs(norm_b(ii)) <= M_EPSILON) then
      exception_saved = CNST(0.0)
      status(ii) = QMR_B_ZERO
      residue(ii) = norm_b(ii)
      saved_iter(ii) = iter
      saved_res(ii) = residue(ii)
    end if

  end do

  call X(preconditioner_apply_batch)(this%pre, gr, hm, ik, vvb, zzb, omega = shift)
  call mesh_batch_nrm2(gr%mesh, zzb, xsi)

  gamma = CNST(1.0)
  eta   = CNST(-1.0)
  alpha = CNST(1.0)
  theta = CNST(0.0)

  do while(iter < this%max_iter)
    iter = iter + 1

    if(all(status /= QMR_NOT_CONVERGED)) exit

    do ii = 1, xb%nst
      if(status(ii) == QMR_NOT_CONVERGED .and. (abs(rho(ii)) < M_EPSILON .or. abs(xsi(ii)) < M_EPSILON)) then
        call batch_get_state(xb, ii, gr%mesh%np, exception_saved(:, :, ii))
        status(ii) = QMR_BREAKDOWN_PB
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if

      alpha(ii) = alpha(ii)*xsi(ii)/rho(ii)
    end do

    call batch_scal(gr%mesh%np, CNST(1.0)/rho, vvb, a_full = .false.)
    call batch_scal(gr%mesh%np, CNST(1.0)/xsi, zzb, a_full = .false.)

    call X(mesh_batch_dotp_vector)(gr%mesh, vvb, zzb, delta)

    do ii = 1, xb%nst
      if(status(ii) == QMR_NOT_CONVERGED .and. abs(delta(ii)) < M_EPSILON) then
        call batch_get_state(xb, ii, gr%mesh%np, exception_saved(:, :, ii))
        status(ii) = QMR_BREAKDOWN_VZ
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if
    end do

    if(iter == 1) then
      call batch_copy_data(gr%mesh%np, zzb, qqb)
    else
      call batch_xpay(gr%mesh%np, zzb, -rho*delta/eps, qqb, a_full = .false.)
    end if

    call X(linear_solver_operator_batch)(hm, gr, st, ik, shift, qqb, ppb)

    call batch_scal(gr%mesh%np, alpha, ppb, a_full = .false.)

    call X(mesh_batch_dotp_vector)(gr%mesh, qqb, ppb, eps)

    do ii = 1, xb%nst
      if(status(ii) == QMR_NOT_CONVERGED .and. abs(eps(ii)) < M_EPSILON) then
        call batch_get_state(xb, ii, gr%mesh%np, exception_saved(:, :, ii))
        status(ii) = QMR_BREAKDOWN_QP
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if

      beta(ii) = eps(ii)/delta(ii)
    end do

    call batch_xpay(gr%mesh%np, ppb, -beta, vvb, a_full = .false.)

    forall (ii = 1:xb%nst) oldrho(ii) = rho(ii)

    call mesh_batch_nrm2(gr%mesh, vvb, rho)

    call X(preconditioner_apply_batch)(this%pre, gr, hm, ik, vvb, zzb, omega = shift)

    call batch_scal(gr%mesh%np, CNST(1.0)/alpha, zzb, a_full = .false.)

    call mesh_batch_nrm2(gr%mesh, zzb, xsi)

    do ii = 1, xb%nst
      oldtheta(ii) = theta(ii)
      theta(ii) = rho(ii)/(gamma(ii)*abs(beta(ii)))
      oldgamma = gamma(ii)
      gamma(ii) = CNST(1.0)/sqrt(CNST(1.0) + theta(ii)**2)

      if(status(ii) == QMR_NOT_CONVERGED .and. abs(gamma(ii)) < M_EPSILON) then
        call batch_get_state(xb, ii, gr%mesh%np, exception_saved(:, :, ii))
        status(ii) = QMR_BREAKDOWN_GAMMA
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if

      eta(ii) = -eta(ii)*oldrho(ii)*gamma(ii)**2/(beta(ii)*oldgamma**2)
    end do

    if(iter == 1) then

      call batch_copy_data(gr%mesh%np, qqb, deltax)
      call batch_scal(gr%mesh%np, eta*alpha, deltax, a_full = .false.)
      call batch_axpy(gr%mesh%np, CNST(1.0), deltax, xb)
      
      call batch_copy_data(gr%mesh%np, ppb, deltar)
      call batch_scal(gr%mesh%np, eta, deltar, a_full = .false.)
      call batch_axpy(gr%mesh%np, CNST(-1.0), deltar, res)

    else

      call batch_scal(gr%mesh%np, (oldtheta*gamma)**2, deltax, a_full = .false.)
      call batch_axpy(gr%mesh%np, eta*alpha, qqb, deltax, a_full = .false.)
      call batch_axpy(gr%mesh%np, CNST(1.0), deltax, xb)

      call batch_scal(gr%mesh%np, (oldtheta*gamma)**2, deltar, a_full = .false.)
      call batch_axpy(gr%mesh%np, eta, ppb, deltar, a_full = .false.)
      call batch_axpy(gr%mesh%np, CNST(-1.0), deltar, res)

    end if

    call mesh_batch_nrm2(gr%mesh, res, residue)
    forall(ii = 1:xb%nst) residue(ii) = residue(ii)/norm_b(ii)

    do ii = 1, xb%nst
      if(status(ii) == QMR_NOT_CONVERGED .and. residue(ii) < threshold) then
        status(ii) = QMR_CONVERGED
      end if
    end do

  end do

  do ii = 1, xb%nst
    if(status(ii) == QMR_NOT_CONVERGED .or. status(ii) == QMR_CONVERGED) then
      iter_used(ii) = iter
    else
      call batch_set_state(xb, ii, gr%mesh%np, exception_saved(:, :, ii))
      iter_used(ii) = saved_iter(ii)
      residue(ii) = saved_res(ii) 
    end if

    select case(status(ii))
    case(QMR_NOT_CONVERGED)
      write(message(1), '(a)') "QMR solver not converged!"
      write(message(2), '(a)') "Try increasing the maximum number of iterations or the tolerance."
      call messages_warning(2)
    case(QMR_BREAKDOWN_PB)
      write(message(1), '(a)') "QMR breakdown, cannot continue: b or P*b is the zero vector!"
      call messages_warning(1)
    case(QMR_BREAKDOWN_VZ)
      write(message(1), '(a)') "QMR breakdown, cannot continue: v^T*z is zero!"
      call messages_warning(1)
    case(QMR_BREAKDOWN_QP)
      write(message(1), '(a)') "QMR breakdown, cannot continue: q^T*p is zero!"
      call messages_warning(1)
    case(QMR_BREAKDOWN_GAMMA)
      write(message(1), '(a)') "QMR breakdown, cannot continue: gamma is zero!"
      call messages_warning(1)
    end select

  end do

  call batch_end(vvb)
  call batch_end(res)
  call batch_end(zzb)
  call batch_end(qqb)
  call batch_end(ppb)
  call batch_end(deltax)
  call batch_end(deltar)

  SAFE_DEALLOCATE_A(exception_saved)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(oldrho)
  SAFE_DEALLOCATE_A(norm_b)
  SAFE_DEALLOCATE_A(xsi)
  SAFE_DEALLOCATE_A(gamma)
  SAFE_DEALLOCATE_A(alpha)
  SAFE_DEALLOCATE_A(eta)
  SAFE_DEALLOCATE_A(theta)
  SAFE_DEALLOCATE_A(oldtheta)
  SAFE_DEALLOCATE_A(beta)
  SAFE_DEALLOCATE_A(delta)
  SAFE_DEALLOCATE_A(eps)

  SAFE_DEALLOCATE_A(status)
  SAFE_DEALLOCATE_A(saved_iter)

  POP_SUB(X(linear_solver_qmr_dotp))
end subroutine X(linear_solver_qmr_dotp)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
