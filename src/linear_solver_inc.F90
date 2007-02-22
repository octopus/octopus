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
!! -*- coding: utf-8 mode: f90 -*-
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

  case(LS_HX_FIXED)
    call X(ls_solver_hx)(this, h, gr, st, ik, x, y, omega, 1)

  case(LS_HX)
    call X(ls_solver_hx)(this, h, gr, st, ik, x, y, omega, 2)

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

  p=R_TOTYPE(M_ZERO)
  s=R_TOTYPE(M_ZERO)

  ! this will store the preconditioned functions
  ALLOCATE( phat(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE( shat(NP_PART, st%d%dim), NP_PART*st%d%dim)
  phat=R_TOTYPE(M_ZERO)
  shat=R_TOTYPE(M_ZERO)

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  r(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - ( Hp(1:NP,1:st%d%dim) + omega*x(1:NP,1:st%d%dim) )
  rs(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)

  conv_last = .false.

  do iter = 1, ls%max_iter


    rho_1 = X(states_dotp) (gr%m, st%d%dim, rs, r)
    
    if( rho_1 == M_ZERO ) then
      message(1)="rho_1 == MZERO"
      call write_fatal(1)
    end if

    if( iter == 1 ) then 
      p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)
    else
      beta = rho_1/rho_2*alpha/w
      p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) + &
           beta * (p(1:NP,1:st%d%dim) - w*Hp(1:NP,1:st%d%dim))
    end if

    ! preconditioning 
    call X(preconditioner_apply)(ls%pre, gr, h, p, phat, omega)
    call X(Hpsi)(h, gr, phat, Hp, ik)

    !Hp = Hp + omega*phat
    do idim = 1, st%d%dim 
      call lalg_axpy(NP, omega, phat(:, idim), Hp(:, idim))
    end do
    
    alpha = rho_1/X(states_dotp) (gr%m, st%d%dim, rs, Hp)
    
    s(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) - alpha * Hp(1:NP,1:st%d%dim)

    gamma = abs(X(states_dotp) (gr%m, st%d%dim, s, s))

    if( gamma < ls%tol**2 ) then 
      x(1:NP,1:st%d%dim) = x(1:NP,1:st%d%dim) + alpha*phat(1:NP,1:st%d%dim)
      exit
    end if

    call X(preconditioner_apply)(ls%pre, gr, h, s, shat, omega)
    call X(Hpsi)(h, gr, shat, Hs, ik)

    !Hs = Hs + omega*shat
    do idim = 1, st%d%dim 
      call lalg_axpy(NP, omega, shat(:, idim), Hs(:, idim))
    end do

    w = X(states_dotp) (gr%m, st%d%dim, Hs, s) / X(states_dotp) (gr%m, st%d%dim, Hs, Hs)

    x(1:NP,1:st%d%dim) = x(1:NP,1:st%d%dim) &
      + alpha*phat(1:NP,1:st%d%dim) + w*shat(1:NP,1:st%d%dim)

    r(1:NP,1:st%d%dim) = s(1:NP,1:st%d%dim) - w * Hs(1:NP,1:st%d%dim)

    rho_2=rho_1

    gamma = abs(X(states_dotp) (gr%m, st%d%dim, r, r))
    conv = (gamma < ls%tol**2)

    if( conv .and. conv_last ) then 
      exit
    end if
    conv_last = conv

    if( w == M_ZERO ) then
      message(1)="w == MZERO"
      call write_fatal(1)
    end if
        
  end do
    
  ls%iter = iter
  ls%abs_psi=sqrt(gamma)

  deallocate(r, p, Hp, s, rs, Hs)
  deallocate(phat, shat)

  call pop_sub()
end subroutine X(ls_solver_bicgstab)

!---------------------------------------------------------------------

subroutine X(ls_solver_hx) (ls, h, gr, st, ik, x, y, omega, mode)
  type(linear_solver_t), intent(inout) :: ls
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,                 intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,                 intent(in)    :: omega
  integer,                intent(in)    :: mode ! 1 => fixed 2 => var

  R_TYPE :: w, tau
  R_TYPE, allocatable :: x0(:,:), r(:,:), xi(:,:)
  
  logical :: x0_is_ok

  integer :: total_cg_iter, iter, max_iter

  !separate the frequency
  w = R_REAL(omega)
  tau = -M_zI*R_AIMAG(omega)

  ALLOCATE( x0(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE( xi(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE( r(NP, st%d%dim), NP*st%d%dim)

  x0 = M_ZERO
  xi = M_ZERO

  x0(1:NP_PART, 1:st%d%dim) = x(1:NP_PART, 1:st%d%dim)

  !calculate an approximation of the solution
  call X(Hpsi)(h, gr, x0, r, ik)
  r(1:NP,1:st%d%dim) =  ( r(1:NP,1:st%d%dim) + omega*x0(1:NP,1:st%d%dim) ) - y(1:NP,1:st%d%dim)

  ls%abs_psi = sqrt(abs(X(states_dotp) (gr%m, st%d%dim, r, r)))

  if ( ls%abs_psi < ls%tol ) then
    return
  end if
  
  if ( ls%abs_psi < 1e-02 ) then
    x0_is_ok = .true.
  else
    x0_is_ok = .false.
  end if
  !first term
  call X(ls_solver_cg)(ls, h, gr, st, ik, x, y, w)

  !othogonalize the result against the occupied states 
  !  call X(ls_orth_vector)(gr%m, st, x, ik)

#ifdef R_TREAL
  !we just solved the real case
  return
#endif

  r(1:NP, 1:st%d%dim) = x(1:NP, 1:st%d%dim)

  ASSERT( mode == 1 .or. mode == 2 )
  if( mode == 1 ) max_iter = 1
  if( mode == 2 ) max_iter = 10

  total_cg_iter = ls%iter
  do iter = 1, max_iter

    !estimate an initial value
    if(x0_is_ok) xi(1:NP_PART, 1:st%d%dim) = (x0(1:NP_PART, 1:st%d%dim) - x(1:NP_PART, 1:st%d%dim))/tau

    !multiply by A^-1
    call X(ls_solver_cg)(ls, h, gr, st, ik, xi, r, w)

    !accumulate the number of iterations
    total_cg_iter = total_cg_iter + ls%iter + 1
    
    !othogonalize the result against the occupied states
!    call X(ls_orth_vector)(gr%m, st, xi, ik)
    
    !accumulate
    x(1:NP, 1:st%d%dim) = x(1:NP, 1:st%d%dim) + tau*xi(1:NP, 1:st%d%dim)

    !calculate the residual
    call X(Hpsi)(h, gr, x, r, ik)
    r(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - ( r(1:NP,1:st%d%dim) + omega*x(1:NP,1:st%d%dim) )
    ls%abs_psi = sqrt(abs(X(states_dotp) (gr%m, st%d%dim, r, r)))

    !and check for convergency
    if ( ls%abs_psi < ls%tol ) then
      exit
    end if

    !copy xi to r and multiply by tau
    r(1:NP, 1:st%d%dim) = tau*xi(1:NP, 1:st%d%dim)
    
  end do !i

  ls%iter = total_cg_iter 

  deallocate(x0, xi, r)

end subroutine X(ls_solver_hx)

