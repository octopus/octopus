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
subroutine X(vlaser_operator_quadratic) (laser, gr, std, psi, hpsi)
  type(laser_t),       intent(in)    :: laser
  type(grid_t),        intent(inout) :: gr
  type(states_dim_t),  intent(in)    :: std
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)

  integer :: k
  logical :: vector_potential, magnetic_field

  FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), b(1:MAX_DIM), b_prime(1:MAX_DIM)
  FLOAT, allocatable :: a(:, :), a_prime(:, :)

  call push_sub('h_inc.Xvlaser_operator_quadratic')

  a_field = M_ZERO

  vector_potential = .false.
  magnetic_field = .false.

  select case(laser_kind(laser))
  case(E_FIELD_ELECTRIC) ! do nothing
  case(E_FIELD_MAGNETIC)
    if(.not. allocated(a)) then 
      ALLOCATE(a(NP_PART, gr%mesh%sb%dim), NP_PART*gr%mesh%sb%dim)
      a = M_ZERO
      ALLOCATE(a_prime(NP_PART, gr%mesh%sb%dim), NP_PART*gr%mesh%sb%dim)
      a_prime = M_ZERO
    end if
    call laser_vector_potential(laser, a_prime)
    a = a + a_prime
    call laser_field(gr%sb, laser, b_prime)
    b = b + b_prime
    magnetic_field = .true.
  case(E_FIELD_VECTOR_POTENTIAL)
    call laser_field(gr%sb, laser, a_field_prime)
    a_field = a_field + a_field_prime
    vector_potential = .true.
  end select

  if(magnetic_field) then
    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a(k, 1:gr%mesh%sb%dim), a(k, 1:gr%mesh%sb%dim))*psi(k, :) / P_c**2
    end do
    deallocate(a, a_prime)
  end if
  if(vector_potential) then
    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a_field(1:gr%mesh%sb%dim), a_field(1:gr%mesh%sb%dim))*psi(k, :) / P_c**2
    end do
  end if

  call pop_sub()
end subroutine X(vlaser_operator_quadratic)



! ---------------------------------------------------------
subroutine X(vlaser_operator_linear) (laser, gr, std, psi, hpsi, ik, gyromagnetic_ratio, a_static)
  type(laser_t),       intent(in)    :: laser
  type(grid_t),        intent(inout) :: gr
  type(states_dim_t),  intent(in)    :: std
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, std%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, std%dim)
  integer,             intent(in)    :: ik
  FLOAT,               intent(in)    :: gyromagnetic_ratio
  FLOAT,               pointer       :: a_static(:,:)

  integer :: k, idim
  logical :: electric_field, vector_potential, magnetic_field
  R_TYPE, allocatable :: grad(:, :, :), lhpsi(:, :)

  FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), b(1:MAX_DIM), b_prime(1:MAX_DIM)

  FLOAT, allocatable :: v(:), pot(:), a(:, :), a_prime(:, :)

  call push_sub('h_inc.Xvlaser_operator_linear')

  a_field = M_ZERO

  electric_field = .false.
  vector_potential = .false.
  magnetic_field = .false.

  select case(laser_kind(laser))
  case(E_FIELD_SCALAR_POTENTIAL)
    if(.not. allocated(v)) then 
      ALLOCATE(v(NP), NP)
      v = M_ZERO
    end if
    call laser_potential(gr%sb, laser, gr%mesh, v)
    electric_field = .true.
  case(E_FIELD_ELECTRIC)
    if(.not. allocated(v)) then 
      ALLOCATE(v(NP), NP)
      v = M_ZERO
      ALLOCATE(pot(NP), NP)
      pot = M_ZERO
    end if
    call laser_potential(gr%sb, laser, gr%mesh, pot)
    v = v + pot
    electric_field = .true.
    deallocate(pot)
  case(E_FIELD_MAGNETIC)
    if(.not. allocated(a)) then 
      ALLOCATE(a(NP_PART, gr%mesh%sb%dim), NP_PART*gr%mesh%sb%dim)
      a = M_ZERO
      ALLOCATE(a_prime(NP_PART, gr%mesh%sb%dim), NP_PART*gr%mesh%sb%dim)
      a_prime = M_ZERO
    end if
    call laser_vector_potential(laser, a_prime)
    a = a + a_prime
    call laser_field(gr%sb, laser, b_prime)
    b = b + b_prime
    magnetic_field = .true.
  case(E_FIELD_VECTOR_POTENTIAL)
    call laser_field(gr%sb, laser, a_field_prime)
    a_field = a_field + a_field_prime
    vector_potential = .true.
  end select

  if(electric_field) then
    do k = 1, std%dim
      hpsi(1:NP, k)= hpsi(1:NP, k) + v(1:NP) * psi(1:NP, k)
    end do
    deallocate(v)
  end if

  if(magnetic_field) then

    ALLOCATE(grad(NP_PART, gr%mesh%sb%dim, std%dim), NP_PART*std%dim*gr%mesh%sb%dim)
 
    do idim = 1, std%dim
      call X(derivatives_grad)(gr%der, psi(:, idim), grad(:, :, idim))
    end do

    ! If there is a static magnetic field, its associated vector potential is coupled with
    ! the time-dependent one defined as a "laser" (ideally one should just add them all and
    ! do the calculation only once...). Note that h%ep%a_static already has been divided
    ! by P_c, and therefore here we only divide by P_c, and not P_c**2.
    if(associated(a_static)) then
      do k = 1, NP
        hpsi(k, :) = hpsi(k, :) + dot_product(a(k, 1:gr%mesh%sb%dim), a_static(k, 1:gr%mesh%sb%dim))*psi(k, :) / P_c
      end do
    end if

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a(k, 1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, 1)) / P_c
      end do
    case (SPINORS)
      do k = 1, NP
        do idim = 1, std%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(a(k, 1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, idim)) / P_c
        end do
      end do
    end select


    select case (std%ispin)
    case (SPIN_POLARIZED)
      ALLOCATE(lhpsi(NP, std%dim), NP*std%dim)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(1:NP, 1) = - M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      else
        lhpsi(1:NP, 1) = + M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      end if
      hpsi(1:NP, :) = hpsi(1:NP, :) + (gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
      deallocate(lhpsi)
    case (SPINORS)
      ALLOCATE(lhpsi(NP, std%dim), NP*std%dim)
      lhpsi(1:NP, 1) = M_HALF/P_C*( b(3)*psi(1:NP, 1) &
           + (b(1) - M_zI*b(2))*psi(1:NP, 2))
      lhpsi(1:NP, 2) = M_HALF/P_C*(-b(3)*psi(1:NP, 2) &
           + (b(1) + M_zI*b(2))*psi(1:NP, 1))
      hpsi(1:NP, :) = hpsi(1:NP, :) + (gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
      deallocate(lhpsi)
    end select

    deallocate(grad)
    deallocate(a, a_prime)
  end if

  if(vector_potential) then
    ALLOCATE(grad(NP_PART, gr%mesh%sb%dim, std%dim), NP_PART*std%dim*gr%mesh%sb%dim)

    do idim = 1, std%dim
      call X(derivatives_grad)(gr%der, psi(:, idim), grad(:, :, idim))
    end do

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a_field(1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, 1)) / P_c
      end do
    case (SPINORS)
      do k = 1, NP
        do idim = 1, std%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(a_field(1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, idim)) / P_c
        end do
      end do
    end select
    deallocate(grad)
  end if

  call pop_sub()
end subroutine X(vlaser_operator_linear)


! ---------------------------------------------------------
subroutine X(vlasers) (lasers, nlasers, gr, std, psi, hpsi, grad, ik, t, gyromagnetic_ratio, a_static)
  type(laser_t),       intent(in)    :: lasers(:)
  integer,             intent(in)    :: nlasers
  type(grid_t),        intent(inout) :: gr
  type(states_dim_t),  intent(in)    :: std
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, std%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, std%dim)
  R_TYPE,              intent(in)    :: grad(: , :, :)
  integer,             intent(in)    :: ik
  FLOAT,    optional,  intent(in)    :: t
  FLOAT,               intent(in)    :: gyromagnetic_ratio
  FLOAT,               pointer       :: a_static(:,:)

  integer :: i, k, idim
  logical :: electric_field, vector_potential, magnetic_field
  R_TYPE, allocatable :: lhpsi(:, :)
  type(profile_t), save :: ext_fields_profile

  FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), b(1:MAX_DIM), b_prime(1:MAX_DIM)

  FLOAT, allocatable :: v(:), pot(:), a(:, :), a_prime(:, :)

  call push_sub('h_inc.Xvlasers')
  call profiling_in(ext_fields_profile, 'EXTERNAL_FIELDS')

  a_field = M_ZERO

  electric_field = .false.
  vector_potential = .false.
  magnetic_field = .false.

  do i = 1, nlasers

    select case(laser_kind(lasers(i)))
    case(E_FIELD_SCALAR_POTENTIAL, E_FIELD_ELECTRIC)
      if(.not. allocated(v)) then 
        ALLOCATE(v(NP), NP)
        v = M_ZERO
      end if
      ALLOCATE(pot(NP), NP)
      pot = M_ZERO
      call laser_potential(gr%sb, lasers(i), gr%mesh, pot, t)
      v = v + pot
      electric_field = .true.
      deallocate(pot)

    case(E_FIELD_MAGNETIC)
      if(.not. allocated(a)) then 
        ALLOCATE(a(NP_PART, gr%mesh%sb%dim), NP_PART*gr%mesh%sb%dim)
        a = M_ZERO
        ALLOCATE(a_prime(NP_PART, gr%mesh%sb%dim), NP_PART*gr%mesh%sb%dim)
        a_prime = M_ZERO
      end if

      call laser_vector_potential(lasers(i), a_prime, t)
      a = a + a_prime
      if(present(t)) then
        call laser_field(gr%sb, lasers(i), b_prime, t)
      else 
        call laser_field(gr%sb, lasers(i), b_prime)
      end if
      b = b + b_prime
      magnetic_field = .true.

    case(E_FIELD_VECTOR_POTENTIAL)
      call laser_field(gr%sb, lasers(i), a_field_prime, t)
      a_field = a_field + a_field_prime
      vector_potential = .true.

    end select

  end do

  if(electric_field) then
    do k = 1, std%dim
      hpsi(1:NP, k)= hpsi(1:NP, k) + v(1:NP) * psi(1:NP, k)
    end do
    deallocate(v)
  end if

  if(magnetic_field) then
    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a(k, 1:gr%mesh%sb%dim), a(k, 1:gr%mesh%sb%dim))*psi(k, :) / P_c**2
    end do

    ! If there is a static magnetic field, its associated vector potential is coupled with
    ! the time-dependent one defined as a "laser" (ideally one should just add them all and
    ! do the calculation only once...). Note that h%ep%a_static already has been divided
    ! by P_c, and therefore here we only divide by P_c, and not P_c**2.
    if(associated(a_static)) then
      do k = 1, NP
        hpsi(k, :) = hpsi(k, :) + dot_product(a(k, 1:gr%mesh%sb%dim), a_static(k, 1:gr%mesh%sb%dim))*psi(k, :) / P_c
      end do
    end if

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a(k, 1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, 1)) / P_c
      end do
    case (SPINORS)
      do k = 1, NP
        do idim = 1, std%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(a(k, 1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, idim)) / P_c
        end do
      end do
    end select


    select case (std%ispin)
    case (SPIN_POLARIZED)
      ALLOCATE(lhpsi(NP, std%dim), NP*std%dim)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(1:NP, 1) = - M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      else
        lhpsi(1:NP, 1) = + M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      end if
      hpsi(1:NP, :) = hpsi(1:NP, :) + (gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
      deallocate(lhpsi)
    case (SPINORS)
      ALLOCATE(lhpsi(NP, std%dim), NP*std%dim)
      lhpsi(1:NP, 1) = M_HALF/P_C*( b(3)*psi(1:NP, 1) &
           + (b(1) - M_zI*b(2))*psi(1:NP, 2))
      lhpsi(1:NP, 2) = M_HALF/P_C*(-b(3)*psi(1:NP, 2) &
           + (b(1) + M_zI*b(2))*psi(1:NP, 1))
      hpsi(1:NP, :) = hpsi(1:NP, :) + (gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
      deallocate(lhpsi)
    end select

    deallocate(a, a_prime)
  end if

  if(vector_potential) then
    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a_field(1:gr%mesh%sb%dim), a_field(1:gr%mesh%sb%dim))*psi(k, :) / P_c**2
    end do

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a_field(1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, 1)) / P_c
      end do
    case (SPINORS)
      do k = 1, NP
        do idim = 1, std%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(a_field(1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, idim)) / P_c
        end do
      end do
    end select
  end if

  call profiling_out(ext_fields_profile)
  call pop_sub()
end subroutine X(vlasers)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
