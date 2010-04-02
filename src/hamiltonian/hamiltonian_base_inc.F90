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

subroutine X(hamiltonian_base_local)(this, mesh, std, ispin, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(states_dim_t),          intent(in)    :: std
  integer,                     intent(in)    :: ispin
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)

  call profiling_in(prof_vlpsi, "VLPSI")
  call push_sub('hamiltonian_base_inc.Xhamiltonian_base_local')

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

  call pop_sub('hamiltonian_base_inc.Xhamiltonian_base_local')
  call profiling_out(prof_vlpsi)
end subroutine X(hamiltonian_base_local)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_magnetic)(this, der, std, ep, ispin, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(derivatives_t),         intent(in)    :: der
  type(states_dim_t),          intent(in)    :: std
  type(epot_t),                intent(in)    :: ep
  integer,                     intent(in)    :: ispin
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)
  R_TYPE, allocatable :: grad(:, :, :)
  FLOAT :: a2, cc, b2, bb(1:MAX_DIM)
  CMPLX :: b12

  call profiling_in(prof_magnetic, "MAGNETIC")
  call push_sub('hamiltonian_base_inc.Xhamiltonian_base_magnetic')

  SAFE_ALLOCATE(grad(1:der%mesh%np, 1:MAX_DIM, 1:std%dim))

  do ist = 1, psib%nst
    psi  => psib%states(ist)%X(psi)
    vpsi => vpsib%states(ist)%X(psi)

    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
    grad(:, der%mesh%sb%dim + 1:MAX_DIM, :) = M_ZERO
 
    !! this is commmented for the moment, since this is applied by gauge_field_apply
    !    if(associated(this%uniform_vector_potential)) then
    !      
    !      a2 = sum(this%uniform_vector_potential(1:MAX_DIM)**2)
    !      
    !      forall (idim = 1:psib%dim, ip = 1:der%mesh%np)
    !        vpsi(ip, idim) = vpsi(ip, idim) + M_HALF*a2*psi(ip, idim) &
    !          + M_zI*dot_product(this%uniform_vector_potential(1:MAX_DIM), grad(ip, 1:MAX_DIM, idim))
    !      end forall
    !    end if
    
    if(associated(this%vector_potential)) then
      forall (idim = 1:std%dim, ip = 1:der%mesh%np)
        vpsi(ip, idim) = vpsi(ip, idim) + M_HALF*sum(this%vector_potential(1:MAX_DIM, ip)**2)*psi(ip, idim) &
          + M_zI*dot_product(this%vector_potential(1:MAX_DIM, ip), grad(ip, 1:MAX_DIM, idim))
      end forall
    end if

    if(associated(this%uniform_magnetic_field).and. std%ispin /= UNPOLARIZED) then
      ! Zeeman term
      cc = M_HALF/P_C*ep%gyromagnetic_ratio*M_HALF
      bb = this%uniform_magnetic_field
      b2 = sqrt(sum(this%uniform_magnetic_field**2))
      b12 = bb(1) - M_ZI*bb(2)

      select case (std%ispin)
      case (SPIN_POLARIZED)
        if(is_spin_down(ispin)) cc = -cc
        
        forall (ip = 1:der%mesh%np)
          vpsi(ip, 1) = vpsi(ip, 1) + cc*b2*psi(ip, 1)
        end forall
        
      case (SPINORS)
        forall (ip = 1:der%mesh%np)
          vpsi(ip, 1) = vpsi(ip, 1) + cc*(bb(3)*psi(ip, 1) + b12*psi(ip, 2))
          vpsi(ip, 2) = vpsi(ip, 2) + cc*(-bb(3)*psi(ip, 2) + conjg(b12)*psi(ip, 1))
        end forall
        
      end select
    end if
  end do
  
  SAFE_DEALLOCATE_A(grad)
  
  call pop_sub('hamiltonian_base_inc.Xhamiltonian_base_magnetic')
  call profiling_out(prof_magnetic)
end subroutine X(hamiltonian_base_magnetic)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
