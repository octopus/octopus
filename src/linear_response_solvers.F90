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
!! $Id: linear_response_inc.F90 2196 2006-06-12 15:59:20Z xavier $

!Conjugated gradients
subroutine X(lr_solver_cg) (lr, h, gr, d, ik, x, y, omega, st)
  type(lr_t),          intent(inout) :: lr
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_dim_t),  intent(in)    :: d
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: x(:,:)   ! x(NP, d%dim)
  R_TYPE,                 intent(in)    :: y(:,:)   ! y(NP, d%dim)
  R_TYPE,                 intent(in)    :: omega
  type(states_t), optional, intent(inout) :: st

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter, iunit
  logical :: conv_last, conv

  logical :: orto

  call push_sub('linear_response_inc.Xlr_solver_cg')

  !if we have the states, we reorthogonalize after each application of the hamiltonian
  orto = present(st) .and. lr%ort_each_step

  ALLOCATE( r(NP, d%dim),      NP     *d%dim)
  ALLOCATE( p(NP_PART, d%dim), NP_PART*d%dim)
  ALLOCATE(Hp(NP, d%dim),      NP     *d%dim)

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  r(1:NP,1:d%dim) = y(1:NP,1:d%dim) - ( Hp(1:NP,1:d%dim) + omega*x(1:NP,1:d%dim) )
  
  ! Initial search direction
  p(1:NP,1:d%dim) = r(1:NP,1:d%dim)
  p((NP+1):NP_PART,1:d%dim)=M_ZERO
  
  conv_last = .false.
  do iter = 1, lr%max_iter
    gamma = X(states_dotp) (gr%m, d%dim, r, r)

    ! we require more precision here than for the density
    conv = (sqrt(abs(gamma)) < lr%conv_abs_dens)
    if(conv.and.conv_last) exit
    conv_last = conv
    
    if(orto) call X(lr_orth_vector)(gr%m, st, p, ik)   

    call X(Hpsi)(h, gr, p, Hp, ik)
    Hp(1:NP,1:d%dim) = Hp(1:NP,1:d%dim) + omega*p(1:NP,1:d%dim)
    
    alpha = gamma/ X(states_dotp) (gr%m, d%dim, p, Hp)
    
    r(1:NP,1:d%dim) = r(1:NP,1:d%dim) - alpha*Hp(1:NP,1:d%dim)
    x(1:NP_PART,1:d%dim) = x(1:NP_PART,1:d%dim) + alpha* p(1:NP_PART,1:d%dim)
    
    beta = X(states_dotp) (gr%m, d%dim, r, r)/gamma
    p(1:NP,1:d%dim) = r(1:NP,1:d%dim) + beta*p(1:NP,1:d%dim)
    
  end do
    
  lr%iter     = iter
  lr%abs_dens = gamma

  deallocate(r, p, Hp)

  call pop_sub()
end subroutine X(lr_solver_cg)

!BICONJUGATED GRADIENTS
!Saad Page 210
subroutine X(lr_solver_bcg) (lr, h, gr, d, ik, x, y, omega, st)
  type(lr_t),          intent(inout) :: lr
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_dim_t),  intent(in)    :: d
  type(states_t), optional, intent(inout) :: st

  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: x(:,:)   ! x(NP, d%dim)
  R_TYPE,                 intent(in)    :: y(:,:)   ! y(NP, d%dim)
  R_TYPE,                 intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:), rs(:,:), ps(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter, iunit, id
  logical :: conv_last, conv

  logical :: orto

  call push_sub('linear_response_inc.Xlr_solver_bcg')

  !if we have the states, we reorthogonalize after each application of the hamiltonian
  orto = present(st) .and. lr%ort_each_step

  ALLOCATE( r(NP, d%dim),      NP     *d%dim)
  ALLOCATE( p(NP_PART, d%dim), NP_PART*d%dim)
  ALLOCATE(rs(NP, d%dim),      NP     *d%dim)
  ALLOCATE(ps(NP_PART, d%dim), NP_PART*d%dim)
  ALLOCATE(Hp(NP, d%dim),      NP     *d%dim)

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  r(1:NP,1:d%dim) = y(1:NP,1:d%dim) - ( Hp(1:NP,1:d%dim) + omega*x(1:NP,1:d%dim) )
  rs(1:NP,1:d%dim) = y(1:NP,1:d%dim) - ( Hp(1:NP,1:d%dim) + R_CONJ(omega)*x(1:NP,1:d%dim))

  ! Initial search direction
  p(1:NP,1:d%dim) = r(1:NP,1:d%dim)
  p((NP+1):NP_PART,1:d%dim)=M_ZERO

  ps(1:NP,1:d%dim) = rs(1:NP,1:d%dim)
  ps((NP+1):NP_PART,1:d%dim)=M_ZERO
  
  conv_last = .false.
  do iter = 1, lr%max_iter
    gamma = X(states_dotp) (gr%m, d%dim, r, rs)

    ! we require more precision here than for the density
    conv = (sqrt(abs(gamma)) < lr%conv_abs_dens)
    if(conv.and.conv_last) exit
    conv_last = conv
    
    if(orto) call X(lr_orth_vector)(gr%m, st, p, ik)

    call X(Hpsi)(h, gr, p, Hp, ik)
    Hp(1:NP,1:d%dim) = Hp(1:NP,1:d%dim) + omega*p(1:NP,1:d%dim)
    
    alpha = gamma / X(states_dotp) (gr%m, d%dim, Hp, ps)

    x(1:NP_PART,1:d%dim) = x(1:NP_PART,1:d%dim) + alpha*p(1:NP_PART,1:d%dim)
    r(1:NP,1:d%dim) = r(1:NP,1:d%dim) - alpha*Hp(1:NP,1:d%dim)

    if(orto) call X(lr_orth_vector)(gr%m, st, ps, ik)   

    call X(Hpsi)(h, gr, ps, Hp, ik)
    Hp(1:NP,1:d%dim) = Hp(1:NP,1:d%dim) + R_CONJ(omega)*p(1:NP,1:d%dim)

    rs(1:NP,1:d%dim) = rs(1:NP,1:d%dim) - alpha*Hp(1:NP,1:d%dim)
    
    beta = X(states_dotp) (gr%m, d%dim, r, rs)/gamma

    p(1:NP,1:d%dim) = r(1:NP,1:d%dim) + beta*p(1:NP,1:d%dim)
    ps(1:NP,1:d%dim) = rs(1:NP,1:d%dim) + beta*ps(1:NP,1:d%dim)
    
  end do
    
  lr%iter     = iter
  lr%abs_dens = gamma

  deallocate(r, p, Hp)

  call pop_sub()
end subroutine X(lr_solver_bcg)



!BICONJUGATED GRADIENTS STABILIZED
!see http://math.nist.gov/iml++/bicgstab.h.txt

subroutine X(lr_solver_bicgstab) (lr, h, gr, d, ik, x, y, omega, st)
  type(lr_t),          intent(inout) :: lr
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_dim_t),  intent(in)    :: d
  type(states_t), optional, intent(inout) :: st

  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: x(:,:)   ! x(NP, d%dim)
  R_TYPE,                 intent(in)    :: y(:,:)   ! y(NP, d%dim)
  R_TYPE,                 intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:), rs(:,:), s(:,:), Hs(:,:)
  R_TYPE  :: alpha, beta, gamma, w, normy, rho_1, rho_2
  integer :: iter, iunit, id
  logical :: conv_last, conv

  logical :: orto

  call push_sub('linear_response_inc.Xlr_solver_bcg')

  !if we have the states, we reorthogonalize after each application of the hamiltonian
  orto = present(st)  .and. lr%ort_each_step

  ALLOCATE( r(NP, d%dim), NP     *d%dim)
  ALLOCATE( p(NP_PART, d%dim), NP_PART*d%dim)
  ALLOCATE(rs(NP, d%dim),      NP     *d%dim)
  ALLOCATE( s(NP_PART, d%dim), NP_PART*d%dim)
  ALLOCATE(Hp(NP, d%dim),      NP     *d%dim)
  ALLOCATE(Hs(NP, d%dim),      NP     *d%dim)

  p=M_ZERO
  s=M_ZERO

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  r(1:NP,1:d%dim) = y(1:NP,1:d%dim) - ( Hp(1:NP,1:d%dim) + omega*x(1:NP,1:d%dim) )
  rs(1:NP,1:d%dim) = r(1:NP,1:d%dim)

  normy = X(states_dotp) (gr%m, d%dim, y, y)

  if(normy==M_ZERO) normy=M_ONE

  do iter = 1, lr%max_iter
    
    rho_1 = X(states_dotp) (gr%m, d%dim, rs, r)
    
    if( rho_1 == M_ZERO ) then
      message(1)="rho_1 == MZERO"
      call write_fatal(1)
      exit 
    end if

    if( iter == 1 ) then 
      p(1:NP,1:d%dim) = r(1:NP,1:d%dim)
    else
      beta = rho_1/rho_2*alpha/w
      p(1:NP,1:d%dim) = r(1:NP,1:d%dim) + &
           beta * (p(1:NP,1:d%dim) - w*Hp(1:NP,1:d%dim))
    end if

    if(orto) call X(lr_orth_vector)(gr%m, st, p, ik)   

    call X(Hpsi)(h, gr, p, Hp, ik)
    Hp(1:NP,1:d%dim) = Hp(1:NP,1:d%dim) + omega*p(1:NP,1:d%dim)
    
    alpha = rho_1/X(states_dotp) (gr%m, d%dim, rs, Hp)
    
    s(1:NP,1:d%dim) = r(1:NP,1:d%dim) - alpha * Hp(1:NP,1:d%dim)

    if( abs(X(states_dotp) (gr%m, d%dim, s, s)/normy) < lr%conv_abs_dens**2 ) then 
      lr%abs_dens = X(states_dotp) (gr%m, d%dim, s, s)/normy
      exit
    end if

    if(orto) call X(lr_orth_vector)(gr%m, st, s, ik)   
    call X(Hpsi)(h, gr, s, Hs, ik)
    Hs(1:NP,1:d%dim) = Hs(1:NP,1:d%dim) + omega*s(1:NP,1:d%dim)

    w = X(states_dotp) (gr%m, d%dim, Hs, s) / X(states_dotp) (gr%m, d%dim, Hs, Hs)

    x(1:NP,1:d%dim) = x(1:NP,1:d%dim) &
         + alpha*p(1:NP,1:d%dim) + w*s(1:NP,1:d%dim)
    
    r(1:NP,1:d%dim) = s(1:NP,1:d%dim) - w * Hs(1:NP,1:d%dim)
    
    rho_2=rho_1

    lr%abs_dens = abs(X(states_dotp) (gr%m, d%dim, r, r)/normy)

    if( lr%abs_dens < lr%conv_abs_dens**2 ) then 
      exit
    end if
    
    if( w == M_ZERO ) then
      message(1)="w == MZERO"
      call write_fatal(1)
      exit
    end if
        
  end do
    
  lr%iter = iter

  deallocate(r, p, Hp)

  call pop_sub()
end subroutine X(lr_solver_bicgstab)


