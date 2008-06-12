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
! This subroutine calculates the solution using conjugate gradients
!    of (H + omega) x = y
! ---------------------------------------------------------
subroutine X(solve_HXeY) (this, h, gr, st, ist, ik, x, y, omega)
  type(linear_solver_t), target, intent(inout) :: this
  type(hamiltonian_t),   target, intent(inout) :: h
  type(grid_t),          target, intent(inout) :: gr
  type(states_t),        target, intent(in)    :: st
  integer,                       intent(in)    :: ist
  integer,                       intent(in)    :: ik
  R_TYPE,                        intent(inout) :: x(:,:)   ! x(NP, d%dim)
  R_TYPE,                        intent(in)    :: y(:,:)   ! y(NP, d%dim)
  R_TYPE,                        intent(in)    :: omega

  call profiling_in(prof, "LINEAR_SOLVER")

  select case(this%solver)

  case(LS_CG)
    call X(ls_solver_cg)       (this, h, gr, st, ist, ik, x, y, omega)

  case(LS_BICGSTAB)
    call X(ls_solver_bicgstab) (this, h, gr, st, ist, ik, x, y, omega)

  case(LS_MULTIGRID)
    call X(ls_solver_multigrid)(this, h, gr, st, ist, ik, x, y, omega)

#ifdef R_TCOMPLEX
  case(LS_QMR)
    args%ls       => this
    args%h        => h
    args%gr       => gr 
    args%st       => st
    args%ist      = ist
    args%ik       = ik
    args%X(omega) = omega

    this%iter = this%max_iter
    
    call zqmr_sym(NP, x(:, 1), y(:, 1), X(ls_solver_operator_na), X(ls_preconditioner), &
         this%iter, residue = this%abs_psi, threshold = this%tol, showprogress = .false.)
#endif

  case(LS_SOS)
    call X(ls_solver_sos)(this, h, gr, st, ist, ik, x, y, omega)

  case default 
    write(message(1), '(a,i2)') "Unknown linear response solver", this%solver
    call write_fatal(1)

  end select

  call profiling_out(prof)

end subroutine X(solve_HXeY)


! ---------------------------------------------------------
!Conjugate gradients
subroutine X(ls_solver_cg) (ls, h, gr, st, ist, ik, x, y, omega)
  type(linear_solver_t),          intent(inout) :: ls
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: ist
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,              intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,              intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter, idim
  logical :: conv_last, conv

  call push_sub('linear_response_solvers.Xls_solver_cg')

  ALLOCATE( r(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE( p(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(Hp(NP, st%d%dim),      NP     *st%d%dim)

  ! Initial residue
  call X(ls_solver_operator)(h, gr, st, ist, ik, omega, x, Hp)
  r(1:NP, 1:st%d%dim) = y(1:NP, 1:st%d%dim) - Hp(1:NP, 1:st%d%dim)
  
  ! Initial search direction
  p(1:NP, 1:st%d%dim) = r(1:NP, 1:st%d%dim)
  p((NP+1):NP_PART,1:st%d%dim) = M_ZERO
  
  conv_last = .false.
  do iter = 1, ls%max_iter
    gamma = X(states_dotp)(gr%m, st%d%dim, r, r)

    conv = ( abs(gamma) < ls%tol**2)
    if(conv.and.conv_last) exit
    conv_last = conv
    
    call X(ls_solver_operator)(h, gr, st, ist, ik, omega, p, Hp)

    alpha = gamma/X(states_dotp) (gr%m, st%d%dim, p, Hp)

    do idim = 1, st%d%dim
      !r = r - alpha*Hp
      call lalg_axpy(NP, -alpha, Hp(:, idim), r(:, idim))
      !x = x + alpha*p
      call lalg_axpy(NP,  alpha,  p(:, idim), x(:, idim))
    end do


    beta = X(states_dotp)(gr%m, st%d%dim, r, r)/gamma

    p(1:NP, 1:st%d%dim) = r(1:NP, 1:st%d%dim) + beta*p(1:NP, 1:st%d%dim)

  end do
    
  ls%iter = iter
  ls%abs_psi = sqrt(abs(gamma))

  deallocate(r, p, Hp)

  call pop_sub()
end subroutine X(ls_solver_cg)

! ---------------------------------------------------------
!BICONJUGATE GRADIENTS STABILIZED
!see http://math.nist.gov/iml++/bicgstab.h.txt
subroutine X(ls_solver_bicgstab) (ls, h, gr, st, ist, ik, x, y, omega)
  type(linear_solver_t),          intent(inout) :: ls
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: ist
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,              intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,              intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), Hp(:,:), rs(:,:), Hs(:,:), p(:,:), s(:,:)
  R_TYPE, pointer :: phat(:,:), shat(:,:)
  R_TYPE  :: alpha, beta, w, rho_1, rho_2
  logical :: conv_last, conv
  integer :: iter, idim, ip
  FLOAT :: gamma

  call push_sub('linear_response_solver.Xls_solver_bicgstab')

  ALLOCATE( r(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE( p(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(rs(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE( s(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(Hp(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE(Hs(NP, st%d%dim),      NP     *st%d%dim)

  !$omp parallel workshare
  r = R_TOTYPE(M_ZERO)
  p = R_TOTYPE(M_ZERO)
  s = R_TOTYPE(M_ZERO)
  rs = R_TOTYPE(M_ZERO)
  Hp = R_TOTYPE(M_ZERO)
  Hs = R_TOTYPE(M_ZERO)
  !$omp end parallel workshare

  ! this will store the preconditioned functions
  ALLOCATE(phat(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(shat(NP_PART, st%d%dim), NP_PART*st%d%dim)

  !$omp parallel workshare
  phat = R_TOTYPE(M_ZERO)
  shat = R_TOTYPE(M_ZERO)
  !$omp end parallel workshare

  ! Initial residue
  call X(ls_solver_operator) (h, gr, st, ist, ik, omega, x, Hp)

  !$omp parallel workshare
  r(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - Hp(1:NP,1:st%d%dim)
  rs(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)
  !$omp end parallel workshare

  gamma = X(states_nrm2)(gr%m, st%d%dim, r)

  conv_last = .false.

  do iter = 1, ls%max_iter

    rho_1 = X(states_dotp) (gr%m, st%d%dim, rs, r)

    if( abs(rho_1) < M_EPSILON ) exit

    if( iter == 1 ) then
      do idim = 1, st%d%dim
        call lalg_copy(NP, r(:, idim), p(:, idim))
      end do
    else
      beta = rho_1/rho_2*alpha/w
      do idim = 1, st%d%dim
        !$omp parallel do
        do ip = 1, NP
          p(ip, idim) = r(ip, idim) + beta*(p(ip, idim) - w*Hp(ip, idim))
        end do
        !$omp end parallel do
      end do
    end if

    ! preconditioning 
    call X(preconditioner_apply)(ls%pre, gr, h, p, phat, omega)
    call X(ls_solver_operator)(h, gr, st, ist, ik, omega, phat, Hp)
    
    alpha = rho_1/X(states_dotp)(gr%m, st%d%dim, rs, Hp)

    do idim = 1, st%d%dim
      !$omp parallel do
      do ip = 1, NP
        s(ip, idim) = r(ip, idim) - alpha*Hp(ip, idim)
      end do
      !$omp end parallel do
    end do

    gamma = X(states_nrm2) (gr%m, st%d%dim, s)

    if( gamma < ls%tol ) then
      do idim = 1, st%d%dim 
        call lalg_axpy(NP, alpha, phat(:, idim), x(:, idim))
      end do
      exit
    end if

    call X(preconditioner_apply)(ls%pre, gr, h, s, shat, omega)
    call X(ls_solver_operator)(h, gr, st, ist, ik, omega, shat, Hs)

    w = X(states_dotp)(gr%m, st%d%dim, Hs, s)/X(states_dotp) (gr%m, st%d%dim, Hs, Hs)

    do idim = 1, st%d%dim
      !$omp parallel do
      do ip = 1, NP
        x(ip, idim) = x(ip, idim) + alpha*phat(ip, idim) + w*shat(ip, idim)
        r(ip, idim) = s(ip, idim) - w*Hs(ip, idim)
      end do
      !$omp end parallel do
    end do

    rho_2 = rho_1

    gamma = X(states_nrm2)(gr%m, st%d%dim, r)
    conv = (gamma < ls%tol)

    if( conv .and. conv_last ) then 
      exit
    end if
    conv_last = conv

    if( abs(w) < M_EPSILON ) exit

  end do
    
  ls%iter = iter
  ls%abs_psi = gamma

  deallocate(r, p, Hp, s, rs, Hs)
  deallocate(phat, shat)

  call pop_sub()
end subroutine X(ls_solver_bicgstab)


! ---------------------------------------------------------
subroutine X(ls_solver_multigrid) (ls, h, gr, st, ist, ik, x, y, omega)
  type(linear_solver_t), intent(inout) :: ls
  type(hamiltonian_t),   intent(inout) :: h
  type(grid_t),          intent(inout) :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,                intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,                intent(in)    :: omega

  R_TYPE, allocatable :: diag(:,:), hx(:,:), res(:,:)
  integer :: iter

  ALLOCATE(diag(NP, st%d%dim), NP*st%d%dim)
  ALLOCATE(hx(NP, st%d%dim), NP*st%d%dim)
  ALLOCATE(res(NP, st%d%dim), NP*st%d%dim)

  call X(Hpsi_diag)(h, gr, diag, ik)
  diag(1:NP, 1:st%d%dim) = diag(1:NP, 1:st%d%dim) + omega

  do iter = 1, ls%max_iter

    call smoothing(3)

    call smoothing(3)

    !calculate the residue
    call X(ls_solver_operator)(h, gr, st, ist, ik, omega, x, hx)
    res(1:NP, 1:st%d%dim) = hx(1:NP, 1:st%d%dim) - y(1:NP, 1:st%d%dim)
    ls%abs_psi = X(states_nrm2)(gr%m, st%d%dim, res)

    if(ls%abs_psi < ls%tol) exit

    if(in_debug_mode) then 
      write(message(1), *)  "Multigrid: iter ", iter,  ls%abs_psi, abs(X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, ist, ik), x))
      call write_info(1)
    end if

  end do

  ls%iter = iter

contains 

  subroutine smoothing(steps)
    integer, intent(in) :: steps

    integer :: ii, ip, idim
    R_TYPE  :: rr

    do ii = 1, steps

      call X(ls_solver_operator)(h, gr, st, ist, ik, omega, x, hx)

      do idim = 1, st%d%dim
        do ip = 1, gr%m%np
          rr = hx(ip, idim) - y(ip, idim)
          x(ip, idim) = x(ip, idim) - CNST(0.666666) * rr / diag(ip, idim)
        end do
      end do

    end do

    call X(lr_orth_vector)(gr%m, st, x, ist, ik)

  end subroutine smoothing

end subroutine X(ls_solver_multigrid)


! ---------------------------------------------------------
! This routine applies the operator hx = [H (+ Q) + omega] x
subroutine X(ls_solver_operator) (h, gr, st, ist, ik, omega, x, hx)
  type(hamiltonian_t),   intent(inout) :: h
  type(grid_t),          intent(inout) :: gr
  type(states_t),        intent(in)    :: st
  integer,               intent(in)    :: ist
  integer,               intent(in)    :: ik
  R_TYPE,                intent(inout) :: x(:,:)   !  x(NP, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:,:)  ! Hx(NP, st%d%dim)
  R_TYPE,                intent(in)    :: omega

  integer :: idim, jst
  FLOAT   :: alpha_j, proj, dsmear

  call X(Hpsi)(h, gr, x, Hx, ist, ik)

  !Hx = Hx + omega*x
  do idim = 1, st%d%dim
    call lalg_axpy(NP, omega, x(:, idim), Hx(:, idim))
  end do

  ! This is the Q term in Eq. (11) of PRB 51, 6773 (1995)
  ASSERT(.not.st%parallel_in_states)

  dsmear = max(CNST(1e-14), st%smear%dsmear)
  do jst = 1, st%nst
    alpha_j = max(st%smear%e_fermi + M_THREE*dsmear - st%eigenval(jst, ik), M_ZERO)
    if(alpha_j == M_ZERO) cycle
      
    proj = X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:, :, jst, ik), x)
    do idim = 1, st%d%dim
      call lalg_axpy(NP, R_TOTYPE(alpha_j*proj), st%X(psi)(:, idim, jst, ik), Hx(:, idim))
    end do

  end do

end subroutine X(ls_solver_operator)


! ---------------------------------------------------------
subroutine X(ls_solver_operator_na) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !  x(NP, st%d%dim)
  R_TYPE,                intent(out)   :: Hx(:)  ! Hx(NP, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)

  ALLOCATE(tmpx(args%gr%m%np_part, 1), args%gr%m%np_part)
  ALLOCATE(tmpy(args%gr%m%np_part, 1), args%gr%m%np_part)

  call lalg_copy(args%gr%m%np, x, tmpx(:, 1))
  call X(ls_solver_operator)(args%h, args%gr, args%st, args%ist, args%ik, args%X(omega), tmpx, tmpy)
  call lalg_copy(args%gr%m%np, tmpy(:, 1), hx)

  deallocate(tmpx, tmpy)

end subroutine X(ls_solver_operator_na)


! ---------------------------------------------------------
subroutine X(ls_preconditioner) (x, hx)
  R_TYPE,                intent(in)    :: x(:)   !  x(NP, st%d%dim)
  R_TYPE,                intent(out)   :: hx(:)  ! Hx(NP, st%d%dim)

  R_TYPE, allocatable :: tmpx(:, :)
  R_TYPE, allocatable :: tmpy(:, :)

  ALLOCATE(tmpx(args%gr%m%np_part, 1), args%gr%m%np_part)
  ALLOCATE(tmpy(args%gr%m%np_part, 1), args%gr%m%np_part)

  call lalg_copy(args%gr%m%np, x, tmpx(:, 1))
  call X(preconditioner_apply)(args%ls%pre, args%gr, args%h, tmpx, tmpy, args%X(omega))
  call lalg_copy(args%gr%m%np, tmpy(:, 1), hx)

  deallocate(tmpx, tmpy)

end subroutine X(ls_preconditioner)

! ---------------------------------------------------------
subroutine X(ls_solver_sos) (this, h, gr, st, ist, ik, x, y, omega)
  type(linear_solver_t),          intent(inout) :: this
  type(hamiltonian_t),            intent(inout) :: h
  type(grid_t),                   intent(inout) :: gr
  type(states_t),                 intent(in)    :: st
  integer,                        intent(in)    :: ist
  integer,                        intent(in)    :: ik
  R_TYPE,                         intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,                         intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,                         intent(in)    :: omega

  integer :: ist2, idim
  R_TYPE  :: aa
  R_TYPE, allocatable  :: rr(:, :)

  x(1:NP, 1:st%d%dim) = M_ZERO
  
  do ist2 = 1, st%nst

    aa = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, ist2, ik), y)
    aa = aa/(st%eigenval(ist2, ik) - st%eigenval(ist, ik) + omega)

    do idim = 1, st%d%dim
      call lalg_axpy(NP, aa, st%X(psi)(:, idim, ist2, ik), x(:, idim))
    end do

  end do

  ! calculate the residual

  ALLOCATE(rr(1:NP, 1:st%d%dim), NP*st%d%dim)

  call X(ls_solver_operator)(h, gr, st, ist, ik, omega, x, rr)

  do idim = 1, st%d%dim
    call lalg_axpy(NP, -M_ONE, y(:, idim), rr(:, idim))
  end do
  
  this%abs_psi = X(states_nrm2)(gr%m, st%d%dim, rr)
  this%iter = 1

  deallocate(rr)

end subroutine X(ls_solver_sos)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
