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
  type(grid_t),          target, intent(inout) :: gr
  type(states_t),        target, intent(in)    :: st
  integer,                       intent(in)    :: ist
  integer,                       intent(in)    :: ik
  R_TYPE,                        intent(inout) :: x(:,:)   !< x(gr%mesh%np, d%dim)
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
    call X(linear_solver_sos)(this, hm, gr, st, ist, ik, x, y, shift, residue, iter_used)

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
  type(grid_t),          target, intent(inout) :: gr
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
    call X(linear_solver_qmr_dotp)(this, hm, gr, st, ik, xb, yb, &
      shift, iter_used, residue, tol)

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
  type(grid_t),          intent(inout) :: gr
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
  type(grid_t),          intent(inout) :: gr
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
  type(grid_t),          intent(inout) :: gr
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
  type(grid_t),          intent(inout) :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:,:)   !<  x(gr%mesh%np, st%d%dim)
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
  type(grid_t),          intent(inout) :: gr
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
subroutine X(linear_solver_sos) (ls, hm, gr, st, ist, ik, x, y, shift, residue, iter_used)
  type(linear_solver_t),          intent(inout) :: ls
  type(hamiltonian_t),            intent(in)    :: hm
  type(grid_t),                   intent(inout) :: gr
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
subroutine X(linear_solver_qmr_dotp)(this, hm, gr, st, ik, xb, bb, shift, iter_used, &
  residue, threshold, showprogress, converged)
  type(linear_solver_t), intent(inout) :: this
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(inout) :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ik
  type(batch_t),         intent(inout) :: xb
  type(batch_t),         intent(in)    :: bb
  R_TYPE,                intent(in)    :: shift(:)
  integer,               intent(out)   :: iter_used(:) 
  FLOAT,                 intent(out)   :: residue(:)   !< the residue = abs(Ax-b)
  FLOAT,   optional,     intent(in)    :: threshold    !< convergence threshold
  logical, optional,     intent(in)    :: showprogress !< should there be a progress bar
  logical, optional,     intent(out)   :: converged    !< has the algorithm converged
  
  type(batch_t) :: vvb
  R_TYPE, allocatable :: x(:, :), b(:, :), r(:), v(:, :), z(:, :), q(:, :), p(:, :), deltax(:), deltar(:)
  R_TYPE              :: eta, delta, epsilon, beta, rtmp
  FLOAT               :: rho, xsi, gamma, alpha, theta, threshold_, res, oldtheta, oldgamma, oldrho, tmp, norm_b
  integer             :: err, ip, ilog_res, ilog_thr, ii, iter, idim
  logical             :: showprogress_

  integer, parameter ::        &
    QMR_NORMAL           = 0,  &
    QMR_BREAKDOWN_PB     = 1,  &
    QMR_BREAKDOWN_VZ     = 2,  &
    QMR_BREAKDOWN_QP     = 3,  &
    QMR_BREAKDOWN_GAMMA  = 4
    
  PUSH_SUB(X(linear_solver_qmr_dotp))

  if(present(converged)) converged = .false.
  threshold_ = optional_default(threshold, CNST(1.0e-6))
  showprogress_ = optional_default(showprogress, .false.)

  SAFE_ALLOCATE(x(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(b(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(r(1:gr%mesh%np))
  SAFE_ALLOCATE(v(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(z(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(q(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(p(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(deltax(1:gr%mesh%np))
  SAFE_ALLOCATE(deltar(1:gr%mesh%np))

  call batch_copy(xb, vvb, reference = .false.)

  call X(linear_solver_operator_batch)(hm, gr, st, ik, shift, xb, vvb)

  do ii = 1, xb%nst
    x(1:gr%mesh%np, 1:st%d%dim) = xb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)
    b(1:gr%mesh%np, 1:st%d%dim) = bb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)
    v(1:gr%mesh%np, 1:st%d%dim) = vvb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim)

    forall (ip = 1:gr%mesh%np)
      r(ip) = b(ip, 1) - v(ip, 1)
      v(ip, 1) = r(ip)
    end forall

    rho      = X(mf_nrm2)(gr%mesh, st%d%dim, v)
    norm_b   = X(mf_nrm2)(gr%mesh, st%d%dim, b)

    iter     = 0
    err      = QMR_NORMAL
    res      = rho

    ! If rho is basically zero we are already done.
    if(abs(rho) > M_EPSILON) then
      call X(preconditioner_apply)(this%pre, gr, hm, ik, v, z, omega = shift(ii))

      xsi = X(mf_nrm2)(gr%mesh, st%d%dim, z)

      gamma = M_ONE
      eta   = -M_ONE
      alpha = M_ONE
      theta = M_ZERO

      ! initialize progress bar
      if(showprogress_) then
        ilog_thr = max(M_ZERO, -CNST(100.0)*log(threshold_))
        call loct_progress_bar(-1, ilog_thr)
      end if

      do while(iter < this%max_iter)
        iter = iter + 1
        if((abs(rho) < M_EPSILON) .or. (abs(xsi) < M_EPSILON)) then
          err = QMR_BREAKDOWN_PB
          exit
        end if
        alpha = alpha*xsi/rho

        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, CNST(1.0)/rho, v(:, idim))
          call lalg_scal(gr%mesh%np, CNST(1.0)/xsi, z(:, idim))
        end do

        delta = X(mf_dotp)(gr%mesh, st%d%dim, v, z)

        if(abs(delta) < M_EPSILON) then
          err = QMR_BREAKDOWN_VZ
          exit
        end if

        if(iter == 1) then
          do idim = 1, st%d%dim
            call lalg_copy(gr%mesh%np, z(:, idim), q(:, idim))
          end do
        else
          rtmp = -rho*delta/epsilon
          forall (ip = 1:gr%mesh%np) q(ip, 1) = rtmp*q(ip, 1) + z(ip, 1)
        end if

        call X(linear_solver_operator)(hm, gr, st, xb%states(ii)%ist, ik, shift(ii), q, p)

        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, alpha, p(:, idim))
        end do

        epsilon = X(mf_dotp)(gr%mesh, st%d%dim, q, p)

        if(abs(epsilon) < M_EPSILON) then
          err = QMR_BREAKDOWN_QP
          exit
        end if

        beta = epsilon/delta
        forall (ip = 1:gr%mesh%np) v(ip, 1) = -beta*v(ip, 1) + p(ip, 1)
        oldrho = rho

        rho = X(mf_nrm2)(gr%mesh, st%d%dim, v)

        call X(preconditioner_apply)(this%pre, gr, hm, ik, v, z, omega = shift(ii))

        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, CNST(1.0)/alpha, z(:, idim))
        end do

        xsi = X(mf_nrm2)(gr%mesh, st%d%dim, z)

        oldtheta = theta
        theta    = rho/(gamma*abs(beta))
        oldgamma = gamma
        gamma    = M_ONE/sqrt(M_ONE+theta**2)

        if(abs(gamma) < M_EPSILON) then
          err = QMR_BREAKDOWN_GAMMA
          exit
        end if

        eta = -eta*oldrho*gamma**2/(beta*oldgamma**2)

        rtmp = eta*alpha

        if(iter == 1) then

          forall (ip = 1:gr%mesh%np)
            deltax(ip) = rtmp*q(ip, 1)
            x(ip, 1) = x(ip, 1) + deltax(ip)
          end forall

          forall (ip = 1:gr%mesh%np)
            deltar(ip) = eta*p(ip, 1)
            r(ip) = r(ip) - deltar(ip)
          end forall

        else

          tmp  = (oldtheta*gamma)**2
          forall (ip = 1:gr%mesh%np)
            deltax(ip) = tmp*deltax(ip) + rtmp*q(ip, 1)
            x(ip, 1) = x(ip, 1) + deltax(ip)
          end forall

          forall (ip = 1:gr%mesh%np)
            deltar(ip) = tmp*deltar(ip) + eta*p(ip, 1)
            r(ip) = r(ip) - deltar(ip)
          end forall

        end if

        ! avoid divide by zero
        if(abs(norm_b) < M_EPSILON) then
          res = M_HUGE
        else
          res = X(mf_nrm2)(gr%mesh, r)/norm_b
        endif

        if(showprogress_) then
          ilog_res = CNST(100.0)*max(M_ZERO, -log(res))
          call loct_progress_bar(ilog_res, ilog_thr)
        end if

        if(res < threshold_) exit
      end do
    end if
    
    iter_used(ii) = iter
    xb%states(ii)%X(psi)(1:gr%mesh%np, 1:st%d%dim) = x(1:gr%mesh%np, 1:st%d%dim)

    select case(err)
    case(QMR_NORMAL)
      if(res < threshold_) then
        if (present(converged)) converged = .true.
      else
        write(message(1), '(a)') "QMR solver not converged!"
        write(message(2), '(a)') "Try increasing the maximum number of iterations or the tolerance."
        call messages_warning(2)
      end if
    case(QMR_BREAKDOWN_PB)
      write(message(1), '(a)') "QMR breakdown, cannot continue: b or P*b is the zero vector!"
    case(QMR_BREAKDOWN_VZ)
      write(message(1), '(a)') "QMR breakdown, cannot continue: v^T*z is zero!"
    case(QMR_BREAKDOWN_QP)
      write(message(1), '(a)') "QMR breakdown, cannot continue: q^T*p is zero!"
    case(QMR_BREAKDOWN_GAMMA)
      write(message(1), '(a)') "QMR breakdown, cannot continue: gamma is zero!"
    end select

    if (err /= QMR_NORMAL) then
      write(message(2), '(a)') "Try to change some system parameters (e.g. Spacing, TDTimeStep, ...)."
      call messages_fatal(2)
    end if

    if(showprogress_) write(*,*) ''

    residue(ii) = res

  end do
  
  SAFE_DEALLOCATE_A(x)
  SAFE_DEALLOCATE_A(b)
  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(v)
  SAFE_DEALLOCATE_A(z)
  SAFE_DEALLOCATE_A(q)
  SAFE_DEALLOCATE_A(p)
  SAFE_DEALLOCATE_A(deltax)
  SAFE_DEALLOCATE_A(deltar)
  
  POP_SUB(X(linear_solver_qmr_dotp))
end subroutine X(linear_solver_qmr_dotp)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
