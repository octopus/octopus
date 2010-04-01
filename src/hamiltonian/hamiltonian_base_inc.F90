!! Copyright (C) 2009 X. Andrade
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
!! $Id: hamiltonian_base_inc.F90 3988 2008-03-31 15:06:50Z fnog $

subroutine X(hamiltonian_base_local_batch)(this, mesh, std, ispin, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(states_dim_t),          intent(in)    :: std
  integer,                     intent(in)    :: ispin
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)

  call profiling_in(prof_vlpsi, "VLPSI")
  call push_sub('hamiltonian_base_inc.Xhamiltonian_base_local_batch')

  if(associated(this%potential)) then

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      !$omp parallel do private(idim, ip)
      do ist = 1, psib%nst
        forall (idim = 1:psib%dim, ip = 1:mesh%np)
          vpsib%states(ist)%X(psi)(ip, idim) = this%potential(ip, ispin)*psib%states(ist)%X(psi)(ip, idim)
        end forall
      end do
      !$omp end parallel do

      call profiling_count_operations((R_MUL*psib%nst)*mesh%np)
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

    case(SPINORS)
      !the spinor case is more complicated since it mixes the two components.
      do ist = 1, psib%nst
        psi  => psib%states(ist)%X(psi)
        vpsi => vpsib%states(ist)%X(psi)

        forall(ip = 1:mesh%np)
          vpsi(ip, 1) = this%potential(ip, 1)*psi(ip, 1) + &
            (this%potential(ip, 3) + M_zI*this%potential(ip, 4))*psi(ip, 2)
          vpsi(ip, 2) = this%potential(ip, 2)*psi(ip, 2) + &
            (this%potential(ip, 3) - M_zI*this%potential(ip, 4))*psi(ip, 1)
        end forall

      end do
      call profiling_count_operations((6*R_ADD + 2*R_MUL)*mesh%np*psib%nst)

    end select

  end if

  call pop_sub('hamiltonian_base_inc.Xhamiltonian_base_local_batch')
  call profiling_out(prof_vlpsi)
end subroutine X(hamiltonian_base_local_batch)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
