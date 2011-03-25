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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$


! ---------------------------------------------------------
! This subroutine calculates the solution of (H + shift) x = y
! Typically shift = - eigenvalue + omega
! ---------------------------------------------------------
subroutine X(solve_HXeY) (this, hm, gr, st, ist, ik, x, y, shift, tol, occ_response)
  type(linear_solver_t), target, intent(inout) :: this
  type(hamiltonian_t),   target, intent(in)    :: hm
  type(grid_t),          target, intent(inout) :: gr
  type(states_t),        target, intent(in)    :: st
  integer,                       intent(in)    :: ist
  integer,                       intent(in)    :: ik
  R_TYPE,                        intent(inout) :: x(:,:)   ! x(gr%mesh%np, d%dim)
  R_TYPE,                        intent(in)    :: y(:,:)   ! y(gr%mesh%np, d%dim)
  R_TYPE,                        intent(in)    :: shift
  FLOAT,                         intent(in)    :: tol
  logical, optional,             intent(in)    :: occ_response

  logical :: occ_response_
  R_TYPE, allocatable :: z(:, :)

  PUSH_SUB(X(solve_HXeY))
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
  this%iter = this%max_iter

  select case(this%solver)

  case(LS_CG)
    call X(ls_solver_cg)       (this, hm, gr, st, ist, ik, x, y, shift, tol)

  case(LS_BICGSTAB)
    call X(ls_solver_bicgstab) (this, hm, gr, st, ist, ik, x, y, shift, tol, occ_response_)

  case(LS_MULTIGRID)
    call X(ls_solver_multigrid)(this, hm, gr, st, ist, ik, x, y, shift, tol)

  case(LS_QMR_SYMMETRIC)   
    ! complex symmetric: for Sternheimer, only if wfns are real
    call X(qmr_sym)(gr%mesh%np, x(:, 1), y(:, 1), &
      X(ls_solver_operator_na), X(ls_dotu_qmr), X(ls_nrm2_qmr), X(ls_preconditioner), &
      this%iter, residue = this%abs_psi, threshold = tol, showprogress = .false.)

  case(LS_QMR_SYMMETRIZED)
    ! symmetrized equation
    SAFE_ALLOCATE(z(1:gr%mesh%np, 1:1))
    call X(ls_solver_operator_t_na)(y(:, 1), z(:, 1))
    call X(qmr_sym)(gr%mesh%np, x(:, 1), z(:, 1), &
      X(ls_solver_operator_sym_na), X(ls_dotu_qmr), X(ls_nrm2_qmr), X(ls_preconditioner), &
      this%iter, residue = this%abs_psi, threshold = tol, showprogress = .false.)

  case(LS_QMR_DOTP)
    ! using conjugated dot product
    call X(qmr_sym)(gr%mesh%np, x(:, 1), y(:, 1), &
      X(ls_solver_operator_na), X(ls_dotp_qmr), X(ls_nrm2_qmr), X(ls_preconditioner), &
      this%iter, residue = this%abs_psi, threshold = tol, showprogress = .false.)

  case(LS_QMR_GENERAL)
    ! general algorithm
    call X(qmr)(gr%mesh%np, x(:, 1), y(:, 1), X(ls_solver_operator_na), X(ls_solver_operator_t_na), &
      X(ls_dotu_qmr), X(ls_nrm2_qmr), X(ls_preconditioner), X(ls_preconditioner), &
      this%iter, residue = this%abs_psi, threshold = tol, showprogress = .false.)

  case(LS_SOS)
    call X(ls_solver_sos)(this, hm, gr, st, ist, ik, x, y, shift)

  case default 
    write(message(1), '(a,i2)') "Unknown linear-response solver", this%solver
    call messages_fatal(1)

  end select

  call profiling_out(prof)
  POP_SUB(X(solve_HXeY))

end subroutine X(solve_HXeY)


! ---------------------------------------------------------
FLOAT function X(ls_nrm2_qmr)(x)
  R_TYPE, intent(in) :: x(:)
  
  X(ls_nrm2_qmr) = X(mf_nrm2)(args%gr%mesh, x)
  
end function X(ls_nrm2_qmr)


! ---------------------------------------------------------
R_TYPE function X(ls_dotu_qmr)(x,y)
  R_TYPE, intent(in) :: x(:)
  R_TYPE, intent(in) :: y(:)
  
  X(ls_dotu_qmr) = X(mf_dotp)(args%gr%mesh, x, y, dotu = .true.)
  
end function X(ls_dotu_qmr)

! ---------------------------------------------------------
R_TYPE function X(ls_dotp_qmr)(x,y)
  R_TYPE, intent(in) :: x(:)
  R_TYPE, intent(in) :: y(:)
  
  X(ls_dotp_qmr) = X(mf_dotp)(args%gr%mesh, x, y)
  
end function X(ls_dotp_qmr)

! ---------------------------------------------------------
!Conjugate gradients
subroutine X(ls_solver_cg) (ls, hm, gr, st, ist, ik, x, y, shift, tol)
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

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter, idim
  logical :: conv_last, conv

  PUSH_SUB(X(ls_solver_cg))

  SAFE_ALLOCATE( r(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE( p(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(Hp(1:gr%mesh%np, 1:st%d%dim))

  ! Initial residue
  call X(ls_solver_operator)(hm, gr, st, ist, ik, shift, x, Hp)
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
    
    call X(ls_solver_operator)(hm, gr, st, ist, ik, shift, p, Hp)

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
    
  ls%iter = iter
  ls%abs_psi = sqrt(abs(gamma))

  if(.not. conv) then 
    write(message(1), '(a)') "CG solver not converged!"
    call messages_warning(1)
  endif

  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(p)
  SAFE_DEALLOCATE_A(Hp)

  POP_SUB(X(ls_solver_cg))
end subroutine X(ls_solver_cg)

! ---------------------------------------------------------
!BICONJUGATE GRADIENTS STABILIZED
!see http://math.nist.gov/iml++/bicgstab.h.txt
subroutine X(ls_solver_bicgstab) (ls, hm, gr, st, ist, ik, x, y, shift, tol, occ_response)
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
  logical,               intent(in)    :: occ_response

  R_TYPE, allocatable :: r(:,:), Hp(:,:), rs(:,:), Hs(:,:), p(:,:), s(:,:)
  R_TYPE, pointer :: phat(:,:), shat(:,:)
  R_TYPE  :: alpha, beta, w, rho_1, rho_2
  logical :: conv_last, conv
  integer :: iter, idim, ip
  FLOAT :: gamma
  
  PUSH_SUB(X(ls_solver_bicgstab))

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
  call X(ls_solver_operator) (hm, gr, st, ist, ik, shift, x, Hp)

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
    call X(ls_solver_operator)(hm, gr, st, ist, ik, shift, phat, Hp)
    
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
    call X(ls_solver_operator)(hm, gr, st, ist, ik, shift, shat, Hs)

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

  ls%iter = iter
  ls%abs_psi = gamma

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

  POP_SUB(X(ls_solver_bicgstab))
end subroutine X(ls_solver_bicgstab)


! ---------------------------------------------------------
subroutine X(ls_solver_multigrid) (ls, hm, gr, st, ist, ik, x, y, shift, tol)
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

  R_TYPE, allocatable :: diag(:,:), hx(:,:), res(:,:)
  integer :: iter

  PUSH_SUB(X(ls_solver_multigrid))

  SAFE_ALLOCATE(diag(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(  hx(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE( res(1:gr%mesh%np, 1:st%d%dim))

  call X(hamiltonian_diagonal)(hm, gr%der, diag, ik)
  diag(1:gr%mesh%np, 1:st%d%dim) = diag(1:gr%mesh%np, 1:st%d%dim) + shift

  do iter = 1, ls%max_iter

    call smoothing(3)

    call smoothing(3)

    !calculate the residue
    call X(ls_solver_operator)(hm, gr, st, ist, ik, shift, x, hx)
    res(1:gr%mesh%np, 1:st%d%dim) = hx(1:gr%mesh%np, 1:st%d%dim) - y(1:gr%mesh%np, 1:st%d%dim)
    ls%abs_psi = X(mf_nrm2)(gr%mesh, st%d%dim, res)

    if(ls%abs_psi < tol) exit

    if(in_debug_mode) then 
      write(message(1), *)  "Multigrid: iter ", iter,  ls%abs_psi, abs(X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), x))
      call messages_info(1)
    end if

  end do

  ls%iter = iter

  if(ls%abs_psi > tol) then 
    write(message(1), '(a)') "Multigrid solver not converged!"
    call messages_warning(1)
  endif

  POP_SUB(X(ls_solver_multigrid))

contains 

  subroutine smoothing(steps)
    integer, intent(in) :: steps

    integer :: ii, ip, idim
    R_TYPE  :: rr

    PUSH_SUB(X(ls_solver_multigrid).smoothing)

    do ii = 1, steps

      call X(ls_solver_operator)(hm, gr, st, ist, ik, shift, x, hx)

      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          rr = hx(ip, idim) - y(ip, idim)
          x(ip, idim) = x(ip, idim) - M_TWOTHIRD * rr / diag(ip, idim)
        end do
      end do

    end do

    call X(lr_orth_vector)(gr%mesh, st, x, ist, ik, shift + st%eigenval(ist, ik))

    POP_SUB(X(ls_solver_multigrid).smoothing)
  end subroutine smoothing

end subroutine X(ls_solver_multigrid)


! ---------------------------------------------------------
! This routine applies the operator hx = [H (+ Q) + shift] x
subroutine X(ls_solver_operator) (hm, gr, st, ist, ik, shift, x, hx)
  type(hamiltonian_t),   intent(in)    :: hm
  type(grid_t),          intent(inout) :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:,:)   !  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:,:)  ! Hx(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(in)    :: shift

  integer :: idim, jst
  FLOAT   :: alpha_j, proj

  PUSH_SUB(X(ls_solver_operator))

  call X(hamiltonian_apply)(hm, gr%der, x, Hx, ist, ik)

  !Hx = Hx + shift*x
  do idim = 1, st%d%dim
    call lalg_axpy(gr%mesh%np, shift, x(:, idim), Hx(:, idim))
  end do

  if(st%smear%method == SMEAR_SEMICONDUCTOR .or. st%smear%integral_occs) then
    POP_SUB(X(ls_solver_operator))
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

  POP_SUB(X(ls_solver_operator))

end subroutine X(ls_solver_operator)


! ---------------------------------------------------------
! applies ls_solver_operator with other arguments implicit as global variables
subroutine X(ls_solver_operator_na) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:)  ! Hx(gr%mesh%np, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)

  SAFE_ALLOCATE(tmpx(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpy(1:args%gr%mesh%np, 1:1))

  call lalg_copy(args%gr%mesh%np, x, tmpx(:, 1))
  call X(ls_solver_operator)(args%hm, args%gr, args%st, args%ist, args%ik, args%X(shift), tmpx, tmpy)
  call lalg_copy(args%gr%mesh%np, tmpy(:, 1), hx)

  SAFE_DEALLOCATE_A(tmpx)
  SAFE_DEALLOCATE_A(tmpy)

end subroutine X(ls_solver_operator_na)


! ---------------------------------------------------------
! applies transpose of ls_solver_operator with other arguments implicit as global variables
! (H - shift)^T = H* - shift = (H - shift*)*
subroutine X(ls_solver_operator_t_na) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:)  ! Hx(gr%mesh%np, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)

  SAFE_ALLOCATE(tmpx(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpy(1:args%gr%mesh%np, 1:1))

  call lalg_copy(args%gr%mesh%np, R_CONJ(x), tmpx(:, 1))
  call X(ls_solver_operator)(args%hm, args%gr, args%st, args%ist, args%ik, R_CONJ(args%X(shift)), tmpx, tmpy)
  call lalg_copy(args%gr%mesh%np, R_CONJ(tmpy(:, 1)), hx)

  SAFE_DEALLOCATE_A(tmpx)
  SAFE_DEALLOCATE_A(tmpy)

end subroutine X(ls_solver_operator_t_na)


! ---------------------------------------------------------
! applies ls_solver_operator in symmetrized form: A^T A
subroutine X(ls_solver_operator_sym_na) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:)  ! Hx(gr%mesh%np, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)
  R_TYPE, allocatable :: tmpz(:, :)

  SAFE_ALLOCATE(tmpx(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpy(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpz(1:args%gr%mesh%np_part, 1:1))

  call lalg_copy(args%gr%mesh%np, x, tmpx(:, 1))
  call X(ls_solver_operator)(args%hm, args%gr, args%st, args%ist, args%ik, args%X(shift), tmpx, tmpy)
  call X(ls_solver_operator_t_na)(tmpy(:, 1), tmpz(:, 1))
  call lalg_copy(args%gr%mesh%np, tmpz(:, 1), hx)

  SAFE_DEALLOCATE_A(tmpx)
  SAFE_DEALLOCATE_A(tmpy)
  SAFE_DEALLOCATE_A(tmpz)

end subroutine X(ls_solver_operator_sym_na)

! ---------------------------------------------------------
subroutine X(ls_preconditioner) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !  x(gr%mesh%np, st%d%dim)
  R_TYPE,                intent(out)   :: hx(:)  ! Hx(gr%mesh%np, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)

  PUSH_SUB(X(ls_preconditioner))

  SAFE_ALLOCATE(tmpx(1:args%gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmpy(1:args%gr%mesh%np_part, 1:1))

  call lalg_copy(args%gr%mesh%np, x, tmpx(:, 1))
  call X(preconditioner_apply)(args%ls%pre, args%gr, args%hm, args%ik, tmpx, tmpy, args%X(shift))
  call lalg_copy(args%gr%mesh%np, tmpy(:, 1), hx)

  SAFE_DEALLOCATE_A(tmpx)
  SAFE_DEALLOCATE_A(tmpy)
  POP_SUB(X(ls_preconditioner))

end subroutine X(ls_preconditioner)

! ---------------------------------------------------------
subroutine X(ls_solver_sos) (ls, hm, gr, st, ist, ik, x, y, shift)
  type(linear_solver_t),          intent(inout) :: ls
  type(hamiltonian_t),            intent(in)    :: hm
  type(grid_t),                   intent(inout) :: gr
  type(states_t),                 intent(in)    :: st
  integer,                        intent(in)    :: ist
  integer,                        intent(in)    :: ik
  R_TYPE,                         intent(inout) :: x(:,:)   ! x(gr%mesh%np, st%d%dim)
  R_TYPE,                         intent(in)    :: y(:,:)   ! y(gr%mesh%np, st%d%dim)
  R_TYPE,                         intent(in)    :: shift

  integer :: jst, idim
  R_TYPE  :: aa
  R_TYPE, allocatable  :: rr(:, :)

  PUSH_SUB(X(ls_solver_sos))

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
  call X(ls_solver_operator)(hm, gr, st, ist, ik, shift, x, rr)

  do idim = 1, st%d%dim
    call lalg_axpy(gr%mesh%np, -M_ONE, y(:, idim), rr(:, idim))
  end do
  
  ls%abs_psi = X(mf_nrm2)(gr%mesh, st%d%dim, rr)
  ls%iter = 1

  SAFE_DEALLOCATE_A(rr)
  POP_SUB(X(ls_solver_sos))

end subroutine X(ls_solver_sos)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
