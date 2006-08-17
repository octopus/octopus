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

subroutine X(lr_solve_HXeY) (lr, h, gr, st, ik, x, y, omega)
  type(lr_t),          intent(inout) :: lr
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: x(:,:)   ! x(NP, d%dim)
  R_TYPE,              intent(in)    :: y(:,:)   ! y(NP, d%dim)
  R_TYPE,              intent(in)    :: omega

  select case(lr%solver)

  case(LR_CG)
    call X(lr_solver_cg)(lr, h, gr, st, ik, x, y, omega)

  case(LR_HX_FIXED)
    call X(lr_solver_hx)(lr, h, gr, st, ik, x, y, omega, 1)

  case(LR_HX)
    call X(lr_solver_hx)(lr, h, gr, st, ik, x, y, omega, 2)

  case(LR_BCG)
    call X(lr_solver_bcg)(lr, h, gr, st, ik, x, y, omega)

  case(LR_BICGSTAB)
    call X(lr_solver_bicgstab)(lr, h, gr, st, ik, x, y, omega)

  case default 
    message(1)="Unknown linear response solver"
    call write_fatal(1)

  end select

end subroutine X(lr_solve_HXeY)


subroutine X(preconditioning)(lr, h, gr, st, ik, a, ahat, omega)
  type(lr_t),          intent(inout) :: lr
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: a(:,:)   ! x(NP, d%dim)
  R_TYPE,              intent(inout) :: ahat(:,:)   ! y(NP, d%dim)
  R_TYPE,              intent(in)    :: omega

  call push_sub('linear_response_solver.Xpreconditioning')

  select case(lr%preconditioner)

  case(LR_NONE) 
    ahat(1:NP, 1:st%d%dim)=a(1:NP, 1:st%d%dim)
    
  case(LR_DIAG)
    call X(Hpsi_diag)(h, gr, ahat, ik)
    ahat(1:NP, 1:st%d%dim)=a(1:NP, 1:st%d%dim)/(ahat+omega)
    
  case default
    message(1)="Invalid preconditioner"
    call write_fatal(1)
    
  end select
  
  call pop_sub()

end subroutine X(preconditioning)

!Conjugated gradients
subroutine X(lr_solver_cg) (lr, h, gr, st, ik, x, y, omega)
  type(lr_t),          intent(inout) :: lr
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,                 intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,                 intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter
  logical :: conv_last, conv

  call push_sub('linear_response_solvers.Xlr_solver_cg')

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
  do iter = 1, lr%max_iter
    gamma = X(states_dotp) (gr%m, st%d%dim, r, r)

    conv = ( abs(gamma) < lr%conv_abs_psi**2)
    if(conv.and.conv_last) exit
    conv_last = conv
    
    call X(Hpsi)(h, gr, p, Hp, ik)
    Hp(1:NP,1:st%d%dim) = Hp(1:NP,1:st%d%dim) + omega*p(1:NP,1:st%d%dim)
    
    alpha = gamma/ X(states_dotp) (gr%m, st%d%dim, p, Hp)
    
    r(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) - alpha*Hp(1:NP,1:st%d%dim)
    x(1:NP_PART,1:st%d%dim) = x(1:NP_PART,1:st%d%dim) + alpha* p(1:NP_PART,1:st%d%dim)
    
    beta = X(states_dotp) (gr%m, st%d%dim, r, r)/gamma
    p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) + beta*p(1:NP,1:st%d%dim)
    
  end do
    
  lr%iter = iter
  lr%abs_psi = sqrt(abs(gamma))

  deallocate(r, p, Hp)

  call pop_sub()
end subroutine X(lr_solver_cg)

!BICONJUGATED GRADIENTS
!Saad Page 210
subroutine X(lr_solver_bcg) (lr, h, gr, st, ik, x, y, omega)
  type(lr_t),          intent(inout) :: lr
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: x(:,:)   ! x(NP, st%d%dim)
  R_TYPE,              intent(in)    :: y(:,:)   ! y(NP, st%d%dim)
  R_TYPE,              intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:), rs(:,:), ps(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter
  logical :: conv_last, conv, orto

  call push_sub('linear_response_solvers.Xlr_solver_bcg')

  ALLOCATE( r(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE( p(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(rs(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE(ps(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(Hp(NP, st%d%dim),      NP     *st%d%dim)

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  r(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - ( Hp(1:NP,1:st%d%dim) + omega*x(1:NP,1:st%d%dim) )
  rs(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - ( Hp(1:NP,1:st%d%dim) + R_CONJ(omega)*x(1:NP,1:st%d%dim))

  ! Initial search direction
  p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)
  p((NP+1):NP_PART,1:st%d%dim)=M_ZERO

  ps(1:NP,1:st%d%dim) = rs(1:NP,1:st%d%dim)
  ps((NP+1):NP_PART,1:st%d%dim)=M_ZERO
  
  orto= .false.
  conv_last = .false.
  do iter = 1, lr%max_iter
    if( iter == lr%ort_min_step ) orto = .true.
    
    gamma = X(states_dotp) (gr%m, st%d%dim, r, rs)

    conv = (sqrt(abs(gamma)) < lr%conv_abs_psi)
    if(conv.and.conv_last) exit
    conv_last = conv
    
    if(orto) call X(lr_orth_vector)(gr%m, st, p, ik)

    call X(Hpsi)(h, gr, p, Hp, ik)
    Hp(1:NP,1:st%d%dim) = Hp(1:NP,1:st%d%dim) + omega*p(1:NP,1:st%d%dim)
    
    alpha = gamma / X(states_dotp) (gr%m, st%d%dim, Hp, ps)

    x(1:NP_PART,1:st%d%dim) = x(1:NP_PART,1:st%d%dim) + alpha*p(1:NP_PART,1:st%d%dim)
    r(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) - alpha*Hp(1:NP,1:st%d%dim)

    if(orto) call X(lr_orth_vector)(gr%m, st, ps, ik)   

    call X(Hpsi)(h, gr, ps, Hp, ik)
    Hp(1:NP,1:st%d%dim) = Hp(1:NP,1:st%d%dim) + R_CONJ(omega)*p(1:NP,1:st%d%dim)

    rs(1:NP,1:st%d%dim) = rs(1:NP,1:st%d%dim) - alpha*Hp(1:NP,1:st%d%dim)
    
    beta = X(states_dotp) (gr%m, st%d%dim, r, rs)/gamma

    p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) + beta*p(1:NP,1:st%d%dim)
    ps(1:NP,1:st%d%dim) = rs(1:NP,1:st%d%dim) + beta*ps(1:NP,1:st%d%dim)
    
  end do
    
  lr%iter     = iter
  lr%abs_psi = sqrt(abs(gamma))

  deallocate(r, p, rs, ps, Hp)

  call pop_sub()
end subroutine X(lr_solver_bcg)



!BICONJUGATED GRADIENTS STABILIZED
!see http://math.nist.gov/iml++/bicgstab.h.txt

subroutine X(lr_solver_bicgstab) (lr, h, gr, st, ik, x, y, omega)
  type(lr_t),          intent(inout) :: lr
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
  logical :: conv_last, conv, orto
  integer :: iter
  FLOAT :: gamma

  call push_sub('linear_response_solver.Xlr_solver_bicgstab')

  ALLOCATE( r(NP, st%d%dim), NP     *st%d%dim)
  ALLOCATE( p(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(rs(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE( s(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(Hp(NP, st%d%dim),      NP     *st%d%dim)
  ALLOCATE(Hs(NP, st%d%dim),      NP     *st%d%dim)

  p=M_ZERO
  s=M_ZERO

  if(precondition(lr)) then 
    ALLOCATE( phat(NP_PART, st%d%dim), NP_PART*st%d%dim)
    ALLOCATE( shat(NP_PART, st%d%dim), NP_PART*st%d%dim)
    phat=M_ZERO
    shat=M_ZERO
  end if

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  r(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - ( Hp(1:NP,1:st%d%dim) + omega*x(1:NP,1:st%d%dim) )
  rs(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)

  orto = .false.
  conv_last = .false.
  do iter = 1, lr%max_iter
    if( iter == lr%ort_min_step ) orto = .true.    

    rho_1 = X(states_dotp) (gr%m, st%d%dim, rs, r)
    
    if( rho_1 == M_ZERO ) then
      message(1)="rho_1 == MZERO"
      call write_fatal(1)
      exit 
    end if

    if( iter == 1 ) then 
      p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim)
    else
      beta = rho_1/rho_2*alpha/w
      p(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) + &
           beta * (p(1:NP,1:st%d%dim) - w*Hp(1:NP,1:st%d%dim))
    end if

    if(orto) call X(lr_orth_vector)(gr%m, st, p, ik)   

    if(precondition(lr)) then 
      call X(preconditioning)(lr, h, gr, st, ik, p, phat, omega)
      call X(Hpsi)(h, gr, phat, Hp, ik)
      Hp(1:NP,1:st%d%dim) = Hp(1:NP,1:st%d%dim) + omega*phat(1:NP,1:st%d%dim)
    else 
      call X(Hpsi)(h, gr, p, Hp, ik)
      Hp(1:NP,1:st%d%dim) = Hp(1:NP,1:st%d%dim) + omega*p(1:NP,1:st%d%dim)
    end if
      
    alpha = rho_1/X(states_dotp) (gr%m, st%d%dim, rs, Hp)
    
    s(1:NP,1:st%d%dim) = r(1:NP,1:st%d%dim) - alpha * Hp(1:NP,1:st%d%dim)

    gamma = abs(X(states_dotp) (gr%m, st%d%dim, s, s))

    if( gamma < lr%conv_abs_psi**2 ) then 
       if(precondition(lr)) then
         x(1:NP,1:st%d%dim) = x(1:NP,1:st%d%dim) + alpha*phat(1:NP,1:st%d%dim)
       else
         x(1:NP,1:st%d%dim) = x(1:NP,1:st%d%dim) + alpha*p(1:NP,1:st%d%dim)
       end if
       exit
    end if

    if(orto) call X(lr_orth_vector)(gr%m, st, s, ik, CNST(1e-4))

    if(precondition(lr)) then 
      call X(preconditioning)(lr, h, gr, st, ik, s, shat, omega)
      call X(Hpsi)(h, gr, shat, Hs, ik)
      Hs(1:NP,1:st%d%dim) = Hs(1:NP,1:st%d%dim) + omega*shat(1:NP,1:st%d%dim)
    else 
      call X(Hpsi)(h, gr, s, Hs, ik)
      Hs(1:NP,1:st%d%dim) = Hs(1:NP,1:st%d%dim) + omega*s(1:NP,1:st%d%dim)
    end if

    w = X(states_dotp) (gr%m, st%d%dim, Hs, s) / X(states_dotp) (gr%m, st%d%dim, Hs, Hs)

    if(precondition(lr)) then
      x(1:NP,1:st%d%dim) = x(1:NP,1:st%d%dim) &
           + alpha*phat(1:NP,1:st%d%dim) + w*shat(1:NP,1:st%d%dim)
    else
      x(1:NP,1:st%d%dim) = x(1:NP,1:st%d%dim) &
           + alpha*p(1:NP,1:st%d%dim) + w*s(1:NP,1:st%d%dim)
    end if

    r(1:NP,1:st%d%dim) = s(1:NP,1:st%d%dim) - w * Hs(1:NP,1:st%d%dim)

    rho_2=rho_1

    gamma = abs(X(states_dotp) (gr%m, st%d%dim, r, r))
    conv = (gamma < lr%conv_abs_psi**2)

    if( conv .and. conv_last ) then 
      exit
    end if
    conv_last = conv

    if( w == M_ZERO ) then
      message(1)="w == MZERO"
      call write_fatal(1)
      exit
    end if
        
  end do
    
  lr%iter = iter
  lr%abs_psi=sqrt(gamma)

  deallocate(r, p, Hp, s, rs, Hs)

  if(precondition(lr)) deallocate(phat, shat)

  call pop_sub()
end subroutine X(lr_solver_bicgstab)

!---------------------------------------------------------------------

subroutine X(lr_solver_hx) (lr, h, gr, st, ik, x, y, omega, mode)
  type(lr_t),          intent(inout) :: lr
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

  lr%abs_psi = sqrt(abs(X(states_dotp) (gr%m, st%d%dim, r, r)))

  if ( lr%abs_psi < lr%conv_abs_psi ) then
    return
  end if
  
  if ( lr%abs_psi < 1e-02 ) then
    x0_is_ok = .true.
  else
    x0_is_ok = .false.
  end if
  !first term
  call X(lr_solver_cg)(lr, h, gr, st, ik, x, y, w)

  !othogonalize the result against the occupied states 
  call X(lr_orth_vector)(gr%m, st, x, ik)

#ifdef R_TREAL
  !we just solved the real case
  return
#endif

  r(1:NP, 1:st%d%dim) = x(1:NP, 1:st%d%dim)

  ASSERT( mode == 1 .or. mode == 2 )
  if( mode == 1 ) max_iter = 1
  if( mode == 2 ) max_iter = 10

  total_cg_iter = lr%iter
  do iter = 1, max_iter

    !estimate an initial value
    if(x0_is_ok) xi(1:NP_PART, 1:st%d%dim) = (x0(1:NP_PART, 1:st%d%dim) - x(1:NP_PART, 1:st%d%dim))/tau

    !multiply by A^-1
    call X(lr_solver_cg)(lr, h, gr, st, ik, xi, r, w)

    !accumulate the number of iterations
    total_cg_iter = total_cg_iter + lr%iter + 1
    
    !othogonalize the result against the occupied states
    call X(lr_orth_vector)(gr%m, st, xi, ik)
    
    !accumulate
    x(1:NP, 1:st%d%dim) = x(1:NP, 1:st%d%dim) + tau*xi(1:NP, 1:st%d%dim)

    !calculate the residual
    call X(Hpsi)(h, gr, x, r, ik)
    r(1:NP,1:st%d%dim) = y(1:NP,1:st%d%dim) - ( r(1:NP,1:st%d%dim) + omega*x(1:NP,1:st%d%dim) )
    lr%abs_psi = sqrt(abs(X(states_dotp) (gr%m, st%d%dim, r, r)))

    !and check for convergency
    if ( lr%abs_psi < lr%conv_abs_psi ) then
      exit
    end if

    !copy xi to r and multiply by tau
    r(1:NP, 1:st%d%dim) = tau*xi(1:NP, 1:st%d%dim)
    
  end do !i

  lr%iter = total_cg_iter 

  deallocate(x0, xi, r)

end subroutine X(lr_solver_hx)

