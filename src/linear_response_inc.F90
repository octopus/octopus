!! Copyright (C) 2004 E.S. Kadantsev, M. Marques
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
subroutine X(lr_alloc_fHxc) (st, m, lr)
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: m
  type(lr_t),     intent(out) :: lr

  lr%abs_dens = M_ZERO
  lr%iter     = 0

  ! allocate variables
  ALLOCATE(lr%X(dl_rho) (m%np, st%d%nspin), m%np*st%d%nspin)
  ALLOCATE(lr%X(dl_Vhar)(m%np), m%np)
  ALLOCATE(lr%dl_Vxc(m%np, st%d%nspin, st%d%nspin), m%np*st%d%nspin*st%d%nspin)

end subroutine X(lr_alloc_fHxc)


! ---------------------------------------------------------
subroutine X(lr_orth_vector) (m, st, v, ik)
  type(mesh_t),        intent(in)    :: m
  type(states_t),      intent(in)    :: st
  R_TYPE,                 intent(inout) :: v(:,:)
  integer,                intent(in)    :: ik

  R_TYPE  :: scalp
  integer :: ist

  call push_sub('linear_response_inc.Xlr_orth_vector')

  do ist = 1, st%nst
    if(st%occ(ist, ik) > M_ZERO) then
      scalp = X(states_dotp)(m, st%d%dim, st%X(psi)(:,:, ist, ik), v)
      v(:,:) = v(:,:) - scalp*st%X(psi)(:,:, ist, ik)
    end if
  end do

  call pop_sub()
end subroutine X(lr_orth_vector)


! ---------------------------------------------------------
! calculates
!    lr%dl_rho += sum_{i occ} psi_i^0 (r) * psi_i^1*(r)   <=   type=1
!    lr%dl_rho += sum_{i occ} psi_i^0*(r) * psi_i^1 (r)   <=   type=2
!    type 3 => type 1 + type 2
subroutine X(lr_build_dl_rho) (m, st, lr, type)
  type(mesh_t),   intent(in)    :: m
  type(states_t), intent(in)    :: st
  type(lr_t),     intent(inout) :: lr
  integer,           intent(in)    :: type

  integer :: i, p, ik, sp
  CMPLX   :: c
  R_TYPE  :: d(st%d%nspin)

  call push_sub('linear_response_inc.Xlr_build_dl_rho')

  sp = 1
  if(st%d%ispin == SPIN_POLARIZED) sp = 2

  do ik = 1, st%d%nik, sp
    do p  = st%st_start, st%st_end
      do i = 1, m%np
        d(1) = st%d%kweights(ik)*st%occ(p, ik) * &
          st%X(psi)(i, 1, p, ik)*R_CONJ(lr%X(dl_psi)(i, 1, p, ik))

        select case(st%d%ispin)
        case(SPIN_POLARIZED)
          d(2) = st%d%kweights(ik+1)*st%occ(p, ik+1) * &
            st%X(psi)(i, 1, p, ik+1)*R_CONJ(lr%X(dl_psi)(i, 1, p, ik+1))

        case(SPINORS)
          lr%X(dl_rho)(i, 2) = lr%X(dl_rho)(i, 2) + st%d%kweights(ik)*st%occ(p, ik) * &
            st%X(psi)(i, 2, p, ik)*R_CONJ(lr%X(dl_psi)(i, 2, p, ik))

          c = st%X(psi)(i, 1, p, ik) * R_CONJ(lr%X(dl_psi)(i, 2, p, ik))

          d(3) = st%d%kweights(ik)*st%occ(p, ik) * R_REAL(c)
          d(4) = st%d%kweights(ik)*st%occ(p, ik) * R_AIMAG(c)
        end select

        if(type == 2) d(:) = R_CONJ(d(:))
        if(type == 3) d(:) = R_CONJ(d(:))+d(:)
        lr%X(dl_rho)(i, :) = lr%X(dl_rho)(i, :) + d(:)
      end do
    end do
  end do

  call pop_sub()
end subroutine X(lr_build_dl_rho)


! ---------------------------------------------------------
! This subroutine calculates the solution of
!    (H - eps_{ist,ik} + omega) psi^1_{ist,ik} = y
! ---------------------------------------------------------
subroutine X(lr_solve_HXeY) (lr, h, gr, d, ik, x, y, omega)
  type(lr_t),          intent(inout) :: lr
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_dim_t),  intent(in)    :: d

  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: x(:,:)   ! x(NP, d%dim)
  R_TYPE,                 intent(in)    :: y(:,:)   ! y(NP, d%dim)
  R_TYPE, optional,       intent(in)    :: omega

  R_TYPE, allocatable :: r(:,:), p(:,:), Hp(:,:)
  R_TYPE  :: alpha, beta, gamma
  integer :: iter
  logical :: conv_last, conv

  call push_sub('linear_response_inc.Xlr_solve_HXeY')

  ALLOCATE( r(NP, d%dim),      NP     *d%dim)
  ALLOCATE( p(NP_PART, d%dim), NP_PART*d%dim)
  ALLOCATE(Hp(NP, d%dim),      NP     *d%dim)

  ! Initial residue
  call X(Hpsi)(h, gr, x, Hp, ik)
  r(:,:) = y(:,:) - Hp(:,:)
  if(present(omega)) r(:,:) = r(:,:) - omega*x(:,:)

  ! Initial search direction
  p = R_TOTYPE(M_ZERO)
  p(1:NP,:) = r(1:NP,:)

  conv_last = .false.
  do iter = 1, lr%max_iter
    gamma = X(states_dotp) (gr%m, d%dim, r, r)

    ! we require more precision here than for the density
    conv = (sqrt(abs(gamma)) < lr%conv_abs_dens)
    if(conv.and.conv_last) exit
    conv_last = conv

    call X(Hpsi)(h, gr, p, Hp, ik)
    if(present(omega)) Hp(:,:) = Hp(:,:) + omega*p(:,:)

    alpha = gamma/ X(states_dotp) (gr%m, d%dim, p, Hp)
    r(:,:) = r(:,:) - alpha*Hp(:,:)
    x(:,:) = x(:,:) + alpha* p(:,:)

    beta = X(states_dotp) (gr%m, d%dim, r, r)/gamma
    p(1:NP,:) = r(1:NP,:) + beta*p(1:NP,:)
  end do

  lr%iter     = iter
  lr%abs_dens = gamma
  deallocate(r, p, Hp)

  call pop_sub()
end subroutine X(lr_solve_HXeY)
