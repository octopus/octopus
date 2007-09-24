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
! This subroutine calculates the solution of using conjugated gradients
!    (H + omega) x = y
! ---------------------------------------------------------

subroutine X(solve_HXeY) (this, h, gr, st, ik, x, y, omega)
  type(linear_solver_t), intent(inout) :: this
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: x(:,:)   ! x(NP, d%dim)
  R_TYPE,              intent(in)    :: y(:,:)   ! y(NP, d%dim)
  R_TYPE,              intent(in)    :: omega

  select case(this%solver)

  case(LS_CG)
    call X(ls_solver_cg)(this, h, gr, st, ik, x, y, omega)

  case(LS_BICGSTAB)
    call X(ls_solver_bicgstab)(this, h, gr, st, ik, x, y, omega)

  case default 
    message(1)="Unknown linear response solver"
    call write_fatal(1)

  end select

end subroutine X(solve_HXeY)


!Conjugated gradients
subroutine X(ls_solver_cg) (ls, h, gr, st, ik, x, y, omega)
  type(linear_solver_t),          intent(inout) :: ls
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,                 intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,                 intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter, idim
  logical :: conv_last, conv

  call push_sub('linear_response_solvers.Xls_solver_cg')

  ALLOCATE( r(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE( p(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(Hp(NP, st%d%dim),      NP     *st%d%dim)

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  r(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - ( Hp(1:NP,1:st%d%dim) + omega*x(1:NP,1:st%d%dim) )
  
  ! Initial search direction
  p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)
  p((NP+1):NP_PART,1:st%d%dim)=M_ZERO
  
  conv_last = .false.
  do iter = 1, ls%max_iter
    gamma = X(states_dotp) (gr%m, st%d%dim, r, r)

    conv = ( abs(gamma) < ls%tol**2)
    if(conv.and.conv_last) exit
    conv_last = conv
    
    call X(Hpsi)(h, gr, p, Hp, ik)
    !Hp = Hp + omega*p
    do idim = 1, st%d%dim
      call lalg_axpy(NP, omega, p(:, idim), Hp(:, idim))
    end do

    alpha = gamma/ X(states_dotp) (gr%m, st%d%dim, p, Hp)

    do idim = 1, st%d%dim
      !r = r - alpha*Hp
      call lalg_axpy(NP, -alpha, Hp(:, idim), r(:, idim))
      !x = x + alpha*p
      call lalg_axpy(NP,  alpha,  p(:, idim), x(:, idim))
    end do


    beta = X(states_dotp) (gr%m, st%d%dim, r, r)/gamma

    p(1:NP, 1:st%d%dim) = r(1:NP, 1:st%d%dim) + beta*p(1:NP, 1:st%d%dim)

  end do
    
  ls%iter = iter
  ls%abs_psi = sqrt(abs(gamma))

  deallocate(r, p, Hp)

  call pop_sub()
end subroutine X(ls_solver_cg)

!BICONJUGATED GRADIENTS STABILIZED
!see http://math.nist.gov/iml++/bicgstab.h.txt

subroutine X(ls_solver_bicgstab) (ls, h, gr, st, ik, x, y, omega)
  type(linear_solver_t),          intent(inout) :: ls
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,              intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,              intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), Hp(:,:), rs(:,:), Hs(:,:), p(:,:), s(:,:)
  R_TYPE, pointer :: phat(:,:), shat(:,:)
  R_TYPE  :: alpha, beta, w, rho_1, rho_2
  logical :: conv_last, conv
  integer :: iter, idim
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
  ALLOCATE( phat(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE( shat(NP_PART, st%d%dim), NP_PART*st%d%dim)
  !$omp parallel workshare
  phat = R_TOTYPE(M_ZERO)
  shat = R_TOTYPE(M_ZERO)
  !$omp end parallel workshare

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  !$omp parallel workshare
  r(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - ( Hp(1:NP,1:st%d%dim) + omega*x(1:NP,1:st%d%dim) )
  rs(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)
  !$omp end parallel workshare

  conv_last = .false.

  do iter = 1, ls%max_iter

    rho_1 = X(states_dotp) (gr%m, st%d%dim, rs, r)

    if( abs(rho_1) < M_EPSILON ) exit

    if( iter == 1 ) then
      !$omp parallel workshare
      p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)
      !$omp end parallel workshare
    else
      beta = rho_1/rho_2*alpha/w
      !$omp parallel workshare
      p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) + beta * (p(1:NP,1:st%d%dim) - w*Hp(1:NP,1:st%d%dim))
      !$omp end parallel workshare
    end if

    ! preconditioning 
    call X(preconditioner_apply)(ls%pre, gr, h, p, phat, omega)
    call X(Hpsi)(h, gr, phat, Hp, ik)

    !Hp = Hp + omega*phat
    do idim = 1, st%d%dim 
      call lalg_axpy(NP, omega, phat(:, idim), Hp(:, idim))
    end do
    
    alpha = rho_1/X(states_dotp) (gr%m, st%d%dim, rs, Hp)

    !$omp parallel workshare    
    s(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) - alpha * Hp(1:NP,1:st%d%dim)
    !$omp end parallel workshare

    gamma = X(states_nrm2) (gr%m, st%d%dim, s)

    if( gamma < ls%tol ) then
      !$omp parallel workshare
      x(1:NP,1:st%d%dim) = x(1:NP,1:st%d%dim) + alpha*phat(1:NP,1:st%d%dim)
      !$omp end parallel workshare
      exit
    end if

    call X(preconditioner_apply)(ls%pre, gr, h, s, shat, omega)
    call X(Hpsi)(h, gr, shat, Hs, ik)

    !Hs = Hs + omega*shat
    do idim = 1, st%d%dim 
      call lalg_axpy(NP, omega, shat(:, idim), Hs(:, idim))
    end do

    w = X(states_dotp) (gr%m, st%d%dim, Hs, s) / X(states_dotp) (gr%m, st%d%dim, Hs, Hs)

    !$omp parallel workshare
    x(1:NP,1:st%d%dim) = x(1:NP,1:st%d%dim) + alpha*phat(1:NP,1:st%d%dim) + w*shat(1:NP,1:st%d%dim)

    r(1:NP,1:st%d%dim) = s(1:NP,1:st%d%dim) - w * Hs(1:NP,1:st%d%dim)
    !$omp end parallel workshare

    rho_2=rho_1

    gamma = X(states_nrm2) (gr%m, st%d%dim, r)
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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
