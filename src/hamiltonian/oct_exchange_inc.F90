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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! ---------------------------------------------------------
subroutine X(oct_exchange_operator)(this, namespace, mesh, hpsi, ist, ik)
  type(oct_exchange_t), intent(in)    :: this
  type(namespace_t),    intent(in)    :: namespace
  type(mesh_t),         intent(in)    :: mesh
  R_TYPE,               intent(inout) :: hpsi(:, :)
  integer,              intent(in)    :: ist
  integer,              intent(in)    :: ik

  integer :: ik2
  R_TYPE, allocatable :: psi(:, :), psi2(:, :)
  integer :: ip

  PUSH_SUB(X(oct_exchange_operator))

  SAFE_ALLOCATE(psi(1:mesh%np, 1:this%oct_st%d%dim))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:this%oct_st%d%dim))

  select case(this%oct_st%d%ispin)
  case(UNPOLARIZED)
    ASSERT(this%oct_st%d%nik  ==  1)
    call states_elec_get_state(this%oct_st, mesh, ist, 1, psi2)
    do ip = 1, mesh%np
      hpsi(ip, 1) = hpsi(ip, 1) + M_TWO*M_zI*psi2(ip, 1)*(this%oct_pot(ip, 1) + this%oct_fxc(ip, 1, 1)*this%oct_rho(ip, 1))
    end do

  case(SPIN_POLARIZED)
    ASSERT(this%oct_st%d%nik  ==  2)

    call states_elec_get_state(this%oct_st, mesh, ist, ik, psi2)

    do ik2 = 1, 2
      do ip = 1, mesh%np
        hpsi(ip, 1) = hpsi(ip, 1) + M_TWO * M_zI * this%oct_st%occ(ist, ik) * &
          psi2(ip, 1) * (this%oct_pot(ip, ik2) + this%oct_fxc(ip, ik, ik2)*this%oct_rho(ip, ik2))
       end do
     end do

  case(SPINORS)
    call messages_not_implemented("Function oct_exchange_operator for spinors", &
      namespace=namespace)
  end select

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(psi2)
  POP_SUB(X(oct_exchange_operator))
end subroutine X(oct_exchange_operator)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
