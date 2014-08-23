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

  case(LS_CG)
    call X(linear_solver_cg)       (this, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used)

  case(LS_BICGSTAB)
    call X(linear_solver_bicgstab) (this, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used, occ_response_)

  case(LS_MULTIGRID)
    call X(linear_solver_multigrid)(this, hm, gr, st, ist, ik, x, y, shift, tol, residue, iter_used)

  case(LS_QMR_SYMMETRIC)   
    ! complex symmetric: for Sternheimer, only if wfns are real
    call X(qmr_sym_gen_dotu)(gr%mesh%np, x(:, 1), y(:, 1), &
      X(linear_solver_operator_na), X(mf_dotu_aux), X(mf_nrm2_aux), X(linear_solver_preconditioner), &
      iter_used, residue = residue, threshold = tol, showprogress = .false.)

  case(LS_QMR_SYMMETRIZED)
    ! symmetrized equation
    SAFE_ALLOCATE(z(1:gr%mesh%np, 1:1))
    call X(linear_solver_operator_t_na)(y(:, 1), z(:, 1))
    call X(qmr_sym_gen_dotu)(gr%mesh%np, x(:, 1), z(:, 1), &
      X(linear_solver_operator_sym_na), X(mf_dotu_aux), X(mf_nrm2_aux), X(linear_solver_preconditioner), &
      iter_used, residue = residue, threshold = tol, showprogress = .false.)

  case(LS_QMR_DOTP)
    ! using conjugated dot product
    call X(qmr_sym_gen_dotu)(gr%mesh%np, x(:, 1), y(:, 1), &
      X(linear_solver_operator_na), X(mf_dotp_aux), X(mf_nrm2_aux), X(linear_solver_preconditioner), &
      iter_used, residue = residue, threshold = tol, showprogress = .false.)

  case(LS_QMR_GENERAL)
    ! general algorithm
    call X(qmr_gen_dotu)(gr%mesh%np, x(:, 1), y(:, 1), X(linear_solver_operator_na), X(linear_solver_operator_t_na), &
      X(mf_dotu_aux), X(mf_nrm2_aux), X(linear_solver_preconditioner), X(linear_solver_preconditioner), &
      iter_used, residue = residue, threshold = tol, showprogress = .false.)

  case(LS_SOS)
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
  case(LS_QMR_DOTP)
    call X(linear_solver_qmr_dotp)(this, hm, gr, st, ik, xb, yb, shift, iter_used, residue, tol)

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
  endif

  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(p)
  SAFE_DEALLOCATE_A(Hp)

  POP_SUB(X(linear_solver_cg))
end subroutine X(linear_solver_cg)

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

  R_TYPE, allocatable :: r(:,:), Hp(:,:), rs(:,:), Hs(:,:), p(:,:), s(:,:)
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
    alpha = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), r)
    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, -alpha, st%X(psi)(:, idim, ist, ik), r(:, idim))
    end do
  else
    ! project RHS onto the unoccupied states
    call X(lr_orth_vector)(gr%mesh, st, r, ist, ik, shift + st%eigenval(ist, ik))
  endif
          
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
  endif

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

  R_TYPE, allocatable :: diag(:,:), hx(:,:), res(:,:)
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

    if(in_debug_mode) then 
      write(message(1), *)  "Multigrid: iter ", iter,  residue, abs(X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), x))
      call messages_info(1)
    end if

  end do

  iter_used = iter

  if(residue > tol) then 
    write(message(1), '(a)') "Multigrid solver not converged!"
    call messages_warning(1)
  endif

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
  FLOAT   :: alpha_j, proj

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
    
    proj = X(mf_dotp) (gr%mesh, st%d%dim, st%X(psi)(:, :, jst, ik), x)
    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, R_TOTYPE(alpha_j * proj), st%X(psi)(:, idim, jst, ik), Hx(:, idim))
    end do

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
  type(batch_t),         intent(out)   :: hxb  

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

  PUSH_SUB(X(linear_solver_sos))

  x(1:gr%mesh%np, 1:st%d%dim) = M_ZERO
  
  do jst = 1, st%nst
    if(ist == jst) cycle

    aa = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, jst, ik), y)
    aa = aa/(st%eigenval(jst, ik) + lr_alpha_j(st, jst, ik) + shift)
    ! Normally the expression in perturbation theory would have here
    ! denominator = st%eigenval(jst, ik) - st%eigenval(ist, ik)
    ! For solving this type of problem, -st%eigenval(ist, ik) is included in shift

    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, aa, st%X(psi)(:, idim, jst, ik), x(:, idim))
    end do
  end do

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

  type(batch_t) :: vvb, res, zzb, qqb, ppb, deltax, deltar, exception_saved
  R_TYPE              :: rtmp
  FLOAT               :: oldgamma, tmp
  integer             :: ip, ii, iter, idim
  FLOAT, allocatable  :: rho(:), oldrho(:), norm_b(:), xsi(:), gamma(:), alpha(:), theta(:), oldtheta(:), saved_res(:)
  R_TYPE, allocatable :: eta(:), beta(:), delta(:), eps(:)
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

  call batch_copy(xb, vvb, reference = .false.)
  call batch_copy(xb, res, reference = .false.)
  call batch_copy(xb, zzb, reference = .false.)
  call batch_copy(xb, qqb, reference = .false.)
  call batch_copy(xb, ppb, reference = .false.)
  call batch_copy(xb, deltax, reference = .false.)
  call batch_copy(xb, deltar, reference = .false.)
  call batch_copy(xb, exception_saved, reference = .false.)

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
      exception_saved%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim) = xb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)
      saved_iter(ii) = iter
      saved_res(ii) = residue(ii)
    end if

    ! if b is zero, the solution is trivial
    if(status(ii) == QMR_NOT_CONVERGED .and. abs(norm_b(ii)) <= M_EPSILON) then
      status(ii) = QMR_B_ZERO
      exception_saved%states(ii)%X(psi) = CNST(0.0)
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
        exception_saved%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim) = xb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)
        status(ii) = QMR_BREAKDOWN_PB
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if

      alpha(ii) = alpha(ii)*xsi(ii)/rho(ii)
    end do

    call batch_scal(gr%mesh%np, CNST(1.0)/rho, vvb, a_full = .false.)
    call batch_scal(gr%mesh%np, CNST(1.0)/xsi, zzb, a_full = .false.)
    
    do ii = 1, xb%nst
      do idim = 1, st%d%dim
!        call lalg_scal(gr%mesh%np, CNST(1.0)/rho(ii), vvb%states(ii)%X(psi)(:, idim))
!        call lalg_scal(gr%mesh%np, CNST(1.0)/xsi(ii), zzb%states(ii)%X(psi)(:, idim))
      end do
    end do

    call X(mesh_batch_dotp_vector)(gr%mesh, vvb, zzb, delta)

    do ii = 1, xb%nst
      if(status(ii) == QMR_NOT_CONVERGED .and. abs(delta(ii)) < M_EPSILON) then
        exception_saved%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim) = xb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)
        status(ii) = QMR_BREAKDOWN_VZ
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if
    end do

    do ii = 1, xb%nst
      if(iter == 1) then
        do idim = 1, st%d%dim
          call lalg_copy(gr%mesh%np, zzb%states(ii)%X(psi)(:, idim), qqb%states(ii)%X(psi)(:, idim))
        end do
      else
        rtmp = -rho(ii)*delta(ii)/eps(ii)
        forall (ip = 1:gr%mesh%np) qqb%states(ii)%X(psi)(ip, 1) = rtmp*qqb%states(ii)%X(psi)(ip, 1) + zzb%states(ii)%X(psi)(ip, 1)
      end if
    end do

    call X(linear_solver_operator_batch)(hm, gr, st, ik, shift, qqb, ppb)

    call batch_scal(gr%mesh%np, alpha, ppb, a_full = .false.)

    do ii = 1, xb%nst
      eps(ii) = X(mf_dotp)(gr%mesh, st%d%dim, qqb%states(ii)%X(psi), ppb%states(ii)%X(psi))

      if(status(ii) == QMR_NOT_CONVERGED .and. abs(eps(ii)) < M_EPSILON) then
        exception_saved%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim) = xb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)
        status(ii) = QMR_BREAKDOWN_QP
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if

      beta(ii) = eps(ii)/delta(ii)
    end do

    do ii = 1, xb%nst
      forall (ip = 1:gr%mesh%np) 
        vvb%states(ii)%X(psi)(ip, 1) = -beta(ii)*vvb%states(ii)%X(psi)(ip, 1) + ppb%states(ii)%X(psi)(ip, 1)
      end forall
      oldrho(ii) = rho(ii)
    end do

    do ii = 1, xb%nst
      rho(ii) = X(mf_nrm2)(gr%mesh, st%d%dim, vvb%states(ii)%X(psi))
    end do

    call X(preconditioner_apply_batch)(this%pre, gr, hm, ik, vvb, zzb, omega = shift)

    call batch_scal(gr%mesh%np, CNST(1.0)/alpha, zzb, a_full = .false.)

    do ii = 1, xb%nst
      xsi(ii) = X(mf_nrm2)(gr%mesh, st%d%dim, zzb%states(ii)%X(psi))
    end do

    do ii = 1, xb%nst
      oldtheta(ii) = theta(ii)
      theta(ii) = rho(ii)/(gamma(ii)*abs(beta(ii)))
      oldgamma = gamma(ii)
      gamma(ii) = CNST(1.0)/sqrt(CNST(1.0) + theta(ii)**2)

      if(status(ii) == QMR_NOT_CONVERGED .and. abs(gamma(ii)) < M_EPSILON) then
        exception_saved%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim) = xb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)
        status(ii) = QMR_BREAKDOWN_GAMMA
        saved_iter(ii) = iter
        saved_res(ii) = residue(ii)
      end if

      eta(ii) = -eta(ii)*oldrho(ii)*gamma(ii)**2/(beta(ii)*oldgamma**2)
    end do

    do ii = 1, xb%nst

      rtmp = eta(ii)*alpha(ii)

      if(iter == 1) then

        forall (ip = 1:gr%mesh%np)
          deltax%states(ii)%X(psi)(ip, 1) = rtmp*qqb%states(ii)%X(psi)(ip, 1)
          xb%states(ii)%X(psi)(ip, 1) = xb%states(ii)%X(psi)(ip, 1) + deltax%states(ii)%X(psi)(ip, 1)
        end forall

        forall (ip = 1:gr%mesh%np)
          deltar%states(ii)%X(psi)(ip, 1) = eta(ii)*ppb%states(ii)%X(psi)(ip, 1)
          res%states(ii)%X(psi)(ip, 1) = res%states(ii)%X(psi)(ip, 1) - deltar%states(ii)%X(psi)(ip, 1)
        end forall

      else

        tmp  = (oldtheta(ii)*gamma(ii))**2
        forall (ip = 1:gr%mesh%np)
          deltax%states(ii)%X(psi)(ip, 1) = tmp*deltax%states(ii)%X(psi)(ip, 1) + rtmp*qqb%states(ii)%X(psi)(ip, 1)
          xb%states(ii)%X(psi)(ip, 1) = xb%states(ii)%X(psi)(ip, 1) + deltax%states(ii)%X(psi)(ip, 1)
        end forall

        forall (ip = 1:gr%mesh%np)
          deltar%states(ii)%X(psi)(ip, 1) = tmp*deltar%states(ii)%X(psi)(ip, 1) + eta(ii)*ppb%states(ii)%X(psi)(ip, 1)
          res%states(ii)%X(psi)(ip, 1) = res%states(ii)%X(psi)(ip, 1) - deltar%states(ii)%X(psi)(ip, 1)
        end forall

      end if
    end do

    do ii = 1, xb%nst
      residue(ii) = X(mf_nrm2)(gr%mesh, st%d%dim, res%states(ii)%X(psi))/norm_b(ii)
    end do

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
      xb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim) = exception_saved%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)
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
  call batch_end(exception_saved)

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
