!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
subroutine X(vlaser_operator_quadratic) (laser, der, std, psi, hpsi)
  type(laser_t),       intent(in)    :: laser
  type(derivatives_t), intent(in)    :: der
  type(states_dim_t),  intent(in)    :: std
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(der%mesh%np_part, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(der%mesh%np_part, h%d%dim)

  integer :: ip
  logical :: vector_potential, magnetic_field

  FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), bb(1:MAX_DIM), b_prime(1:MAX_DIM)
  FLOAT, allocatable :: aa(:, :), a_prime(:, :)

  PUSH_SUB(X(vlaser_operator_quadratic))

  a_field = M_ZERO

  vector_potential = .false.
  magnetic_field = .false.

  select case(laser_kind(laser))
  case(E_FIELD_ELECTRIC) ! do nothing
  case(E_FIELD_MAGNETIC)
    if(.not. allocated(aa)) then 
      SAFE_ALLOCATE(aa(1:der%mesh%np_part, 1:der%mesh%sb%dim))
      aa = M_ZERO
      SAFE_ALLOCATE(a_prime(1:der%mesh%np_part, 1:der%mesh%sb%dim))
      a_prime = M_ZERO
    end if
    a_prime = M_ZERO
    call laser_vector_potential(laser, der%mesh, a_prime)
    aa = aa + a_prime
    b_prime = M_ZERO
    call laser_field(laser, b_prime(1:der%mesh%sb%dim))
    bb = bb + b_prime
    magnetic_field = .true.
  case(E_FIELD_VECTOR_POTENTIAL)
    a_field_prime = M_ZERO
    call laser_field(laser, a_field_prime(1:der%mesh%sb%dim))
    a_field = a_field + a_field_prime
    vector_potential = .true.
  end select

  if(magnetic_field) then
    do ip = 1, der%mesh%np
      hpsi(ip, :) = hpsi(ip, :) + M_HALF * &
        dot_product(aa(ip, 1:der%mesh%sb%dim), aa(ip, 1:der%mesh%sb%dim)) * psi(ip, :) / P_c**2
    end do
    SAFE_DEALLOCATE_A(aa)
    SAFE_DEALLOCATE_A(a_prime)
  end if
  if(vector_potential) then
    do ip = 1, der%mesh%np
      hpsi(ip, :) = hpsi(ip, :) + M_HALF * &
        dot_product(a_field(1:der%mesh%sb%dim), a_field(1:der%mesh%sb%dim))*psi(ip, :) / P_c**2
    end do
  end if

  POP_SUB(X(vlaser_operator_quadratic))
end subroutine X(vlaser_operator_quadratic)

! ---------------------------------------------------------
subroutine X(vlaser_operator_linear) (laser, der, std, psi, hpsi, ik, gyromagnetic_ratio, a_static)
  type(laser_t),       intent(in)    :: laser
  type(derivatives_t), intent(in)    :: der
  type(states_dim_t),  intent(in)    :: std
  R_TYPE,              intent(inout) :: psi(:,:) 
  R_TYPE,              intent(inout) :: hpsi(:,:)
  integer,             intent(in)    :: ik
  FLOAT,               intent(in)    :: gyromagnetic_ratio
  FLOAT, optional,     intent(in)    :: a_static(:,:)

  integer :: ip, idim
  logical :: electric_field, vector_potential, magnetic_field
  R_TYPE, allocatable :: grad(:, :, :), lhpsi(:, :)

  FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), bb(1:MAX_DIM), b_prime(1:MAX_DIM)

  FLOAT, allocatable :: vv(:), pot(:), aa(:, :), a_prime(:, :)

  PUSH_SUB(X(vlaser_operator_linear))

  a_field = M_ZERO

  electric_field = .false.
  vector_potential = .false.
  magnetic_field = .false.

  select case(laser_kind(laser))
  case(E_FIELD_SCALAR_POTENTIAL)
    if(.not. allocated(vv)) then 
      SAFE_ALLOCATE(vv(1:der%mesh%np))
    end if
    vv = M_ZERO
    call laser_potential(laser, der%mesh, vv)
    electric_field = .true.

  case(E_FIELD_ELECTRIC)
    if(.not. allocated(vv)) then 
      SAFE_ALLOCATE(vv(1:der%mesh%np))
      vv = M_ZERO
      SAFE_ALLOCATE(pot(1:der%mesh%np))
    end if
    pot = M_ZERO
    call laser_potential(laser, der%mesh, pot)
    vv = vv + pot
    electric_field = .true.
    SAFE_DEALLOCATE_A(pot)

  case(E_FIELD_MAGNETIC)
    if(.not. allocated(aa)) then 
      SAFE_ALLOCATE(aa(1:der%mesh%np_part, 1:der%mesh%sb%dim))
      aa = M_ZERO
      SAFE_ALLOCATE(a_prime(1:der%mesh%np_part, 1:der%mesh%sb%dim))
      a_prime = M_ZERO
    end if
    a_prime = M_ZERO
    call laser_vector_potential(laser, der%mesh, a_prime)
    aa = aa + a_prime
    b_prime = M_ZERO
    call laser_field(laser, b_prime(1:der%mesh%sb%dim))
    bb = bb + b_prime
    magnetic_field = .true.
  case(E_FIELD_VECTOR_POTENTIAL)
    a_field_prime = M_ZERO
    call laser_field(laser, a_field_prime(1:der%mesh%sb%dim))
    a_field = a_field + a_field_prime
    vector_potential = .true.
  end select

  if(electric_field) then
    do idim = 1, std%dim
      hpsi(1:der%mesh%np, idim)= hpsi(1:der%mesh%np, idim) + vv(1:der%mesh%np) * psi(1:der%mesh%np, idim)
    end do
    SAFE_DEALLOCATE_A(vv)
  end if

  if(magnetic_field) then

    SAFE_ALLOCATE(grad(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:std%dim))
 
    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim))
    end do

    ! If there is a static magnetic field, its associated vector potential is coupled with
    ! the time-dependent one defined as a "laser" (ideally one should just add them all and
    ! do the calculation only once...). Note that h%ep%a_static already has been divided
    ! by P_c, and therefore here we only divide by P_c, and not P_c**2.
    !
    ! We put a minus sign, since for the moment vector potential for
    ! lasers and for the static magnetic field use a different
    ! convetion.
    if(present(a_static)) then
      do ip = 1, der%mesh%np
        hpsi(ip, :) = hpsi(ip, :) - dot_product(aa(ip, 1:der%mesh%sb%dim), a_static(ip, 1:der%mesh%sb%dim)) * psi(ip, :) / P_c
      end do
    end if

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do ip = 1, der%mesh%np
        hpsi(ip, 1) = hpsi(ip, 1) - M_zI * dot_product(aa(ip, 1:der%mesh%sb%dim), grad(ip, 1:der%mesh%sb%dim, 1)) / P_c
      end do
    case (SPINORS)
      do ip = 1, der%mesh%np
        do idim = 1, std%dim
          hpsi(ip, idim) = hpsi(ip, idim) - M_zI * &
            dot_product(aa(ip, 1:der%mesh%sb%dim), grad(ip, 1:der%mesh%sb%dim, idim)) / P_c
        end do
      end do
    end select


    select case (std%ispin)
    case (SPIN_POLARIZED)
      SAFE_ALLOCATE(lhpsi(1:der%mesh%np, 1:std%dim))
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(1:der%mesh%np, 1) = - M_HALF / P_c * sqrt(dot_product(bb, bb)) * psi(1:der%mesh%np, 1)
      else
        lhpsi(1:der%mesh%np, 1) = + M_HALF / P_c * sqrt(dot_product(bb, bb)) * psi(1:der%mesh%np, 1)
      end if
      hpsi(1:der%mesh%np, :) = hpsi(1:der%mesh%np, :) + (gyromagnetic_ratio * M_HALF) * lhpsi(1:der%mesh%np, :)
      SAFE_DEALLOCATE_A(lhpsi)

    case (SPINORS)
      SAFE_ALLOCATE(lhpsi(1:der%mesh%np, 1:std%dim))
      lhpsi(1:der%mesh%np, 1) = M_HALF / P_c * (bb(3) * psi(1:der%mesh%np, 1) &
           + (bb(1) - M_zI * bb(2)) * psi(1:der%mesh%np, 2))
      lhpsi(1:der%mesh%np, 2) = M_HALF / P_c * (-bb(3) * psi(1:der%mesh%np, 2) &
           + (bb(1) + M_zI * bb(2)) * psi(1:der%mesh%np, 1))
      hpsi(1:der%mesh%np, :) = hpsi(1:der%mesh%np, :) + (gyromagnetic_ratio * M_HALF) * lhpsi(1:der%mesh%np, :)
      SAFE_DEALLOCATE_A(lhpsi)
    end select

    SAFE_DEALLOCATE_A(grad)
    SAFE_DEALLOCATE_A(aa)
    SAFE_DEALLOCATE_A(a_prime)
  end if

  if(vector_potential) then
    SAFE_ALLOCATE(grad(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:std%dim))

    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim))
    end do

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do ip = 1, der%mesh%np
        hpsi(ip, 1) = hpsi(ip, 1) - M_zI * dot_product(a_field(1:der%mesh%sb%dim), grad(ip, 1:der%mesh%sb%dim, 1)) / P_c
      end do
    case (SPINORS)
      do ip = 1, der%mesh%np
        do idim = 1, std%dim
          hpsi(ip, idim) = hpsi(ip, idim) - M_zI * &
            dot_product(a_field(1:der%mesh%sb%dim), grad(ip, 1:der%mesh%sb%dim, idim)) / P_c
        end do
      end do
    end select
    SAFE_DEALLOCATE_A(grad)
  end if

  POP_SUB(X(vlaser_operator_linear))
end subroutine X(vlaser_operator_linear)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
