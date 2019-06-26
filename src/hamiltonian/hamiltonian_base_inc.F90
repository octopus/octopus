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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

subroutine X(hamiltonian_base_local)(this, mesh, std, ispin, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(states_dim_t),          intent(in)    :: std
  integer,                     intent(in)    :: ispin
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  PUSH_SUB(X(hamiltonian_base_local))

  if(batch_status(psib) == BATCH_DEVICE_PACKED) then
    ASSERT(.not. allocated(this%Impotential))
    call X(hamiltonian_base_local_sub)(this%potential, mesh, std, ispin, &
      psib, vpsib, potential_opencl = this%potential_opencl)
  else
    if(allocated(this%Impotential)) then
      call X(hamiltonian_base_local_sub)(this%potential, mesh, std, ispin, &
        psib, vpsib, Impotential = this%Impotential)
    else
      call X(hamiltonian_base_local_sub)(this%potential, mesh, std, ispin, &
        psib, vpsib)
    end if
  end if

  POP_SUB(X(hamiltonian_base_local))
end subroutine X(hamiltonian_base_local)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_local_sub)(potential, mesh, std, ispin, psib, vpsib, Impotential, potential_opencl)
  FLOAT,                        intent(in)    :: potential(:,:)
  type(mesh_t),                 intent(in)    :: mesh
  type(states_dim_t),           intent(in)    :: std
  integer,                      intent(in)    :: ispin
  type(batch_t), target,        intent(in)    :: psib
  type(batch_t), target,        intent(inout) :: vpsib
  FLOAT, optional,              intent(in)    :: Impotential(:,:)
  type(accel_mem_t), optional,  intent(in)    :: potential_opencl

  integer :: ist, ip, dim2, dim3
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)
  R_TYPE  :: psi1, psi2
  FLOAT   :: vv, Imvv
  R_TYPE  :: pot(1:4) 
  logical :: pot_is_cmplx
  integer :: pnp, localsize

  call profiling_in(prof_vlpsi, "VLPSI")
  PUSH_SUB(X(hamiltonian_base_local_sub))

  pot_is_cmplx = .false.
  if(present(Impotential)) pot_is_cmplx = .true.

  if(batch_is_packed(psib) .or. batch_is_packed(vpsib)) then
    ASSERT(batch_is_packed(psib))
    ASSERT(batch_is_packed(vpsib))
  end if

  select case(batch_status(psib))
  case(BATCH_DEVICE_PACKED)
    ASSERT(.not. pot_is_cmplx) ! not implemented

    pnp = accel_padded_size(mesh%np)

    select case(std%ispin)

    case(UNPOLARIZED, SPIN_POLARIZED)
      call accel_set_kernel_arg(kernel_vpsi, 0, pnp*(ispin - 1))
      call accel_set_kernel_arg(kernel_vpsi, 1, mesh%np)
      call accel_set_kernel_arg(kernel_vpsi, 2, potential_opencl)
      call accel_set_kernel_arg(kernel_vpsi, 3, psib%pack%buffer)
      call accel_set_kernel_arg(kernel_vpsi, 4, log2(psib%pack%size_real(1)))
      call accel_set_kernel_arg(kernel_vpsi, 5, vpsib%pack%buffer)
      call accel_set_kernel_arg(kernel_vpsi, 6, log2(vpsib%pack%size_real(1)))

      localsize = accel_kernel_workgroup_size(kernel_vpsi)/psib%pack%size_real(1)

      dim3 = mesh%np/(accel_max_size_per_dim(2)*localsize) + 1
      dim2 = min(accel_max_size_per_dim(2)*localsize, pad(mesh%np, localsize))

      call accel_kernel_run(kernel_vpsi, (/psib%pack%size_real(1), dim2, dim3/), (/psib%pack%size_real(1), localsize, 1/))

    case(SPINORS)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 0, mesh%np)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 1, potential_opencl)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 2, pnp)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 3, psib%pack%buffer)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 4, psib%pack%size(1))
      call accel_set_kernel_arg(kernel_vpsi_spinors, 5, vpsib%pack%buffer)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 6, vpsib%pack%size(1))

      call accel_kernel_run(kernel_vpsi_spinors, (/psib%pack%size(1)/2, pnp/), &
        (/psib%pack%size(1)/2, 2*accel_max_workgroup_size()/psib%pack%size(1)/))

    end select

    call accel_finish()

    call profiling_count_operations((R_MUL*psib%nst)*mesh%np)
    call profiling_count_transfers(mesh%np, M_ONE)
    call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

  case(BATCH_PACKED)

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      if(pot_is_cmplx)then
        do ip = 1, mesh%np
          vv = potential(ip, ispin)
          Imvv = Impotential(ip, ispin)
          forall (ist = 1:psib%nst_linear)
            vpsib%pack%X(psi)(ist, ip) = vpsib%pack%X(psi)(ist, ip) + (vv+M_zI*Imvv)*psib%pack%X(psi)(ist, ip)
          end forall
        end do
      else
        !$omp parallel do private(vv, ist)
        do ip = 1, mesh%np
          vv = potential(ip, ispin)
          forall (ist = 1:psib%nst_linear)
            vpsib%pack%X(psi)(ist, ip) = vpsib%pack%X(psi)(ist, ip) + vv*psib%pack%X(psi)(ist, ip)
          end forall
        end do
        !$omp end parallel do
      end if
      call profiling_count_operations((2*R_ADD*psib%nst_linear)*mesh%np)
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst_linear, R_TOTYPE(M_ONE))

    case(SPINORS)
      ASSERT(mod(psib%nst_linear, 2) == 0)
      !the spinor case is more complicated since it mixes the two components.
      if(pot_is_cmplx)then
        !$omp parallel do private(psi1, psi2, ist, pot)
        do ip = 1, mesh%np
          do ist = 1, psib%nst_linear, 2
            psi1 = psib%pack%zpsi(ist    , ip)
            psi2 = psib%pack%zpsi(ist + 1, ip)
            pot(1:4) = potential(ip, 1:4) + M_zI * Impotential(ip, 1:4)
            vpsib%pack%zpsi(ist    , ip) = vpsib%pack%zpsi(ist    , ip) + &
                   pot(1)*psi1 + (pot(3) + M_zI*pot(4))*psi2
            vpsib%pack%zpsi(ist + 1, ip) = vpsib%pack%zpsi(ist + 1, ip) + &
                   pot(2)*psi2 + (pot(3) - M_zI*pot(4))*psi1            
          end do
        end do
        !$omp end parallel do
                
      else
        !$omp parallel do private(psi1, psi2, ist)
        do ip = 1, mesh%np
          do ist = 1, psib%nst_linear, 2
            psi1 = psib%pack%zpsi(ist    , ip)
            psi2 = psib%pack%zpsi(ist + 1, ip)
            vpsib%pack%zpsi(ist    , ip) = vpsib%pack%zpsi(ist    , ip) + &
              potential(ip, 1)*psi1 + (potential(ip, 3) + M_zI*potential(ip, 4))*psi2
            vpsib%pack%zpsi(ist + 1, ip) = vpsib%pack%zpsi(ist + 1, ip) + &
              potential(ip, 2)*psi2 + (potential(ip, 3) - M_zI*potential(ip, 4))*psi1            
          end do
        end do
        !$omp end parallel do
      end if
      
      call profiling_count_operations((6*R_ADD + 2*R_MUL)*mesh%np*psib%nst)

    end select

  case(BATCH_NOT_PACKED)

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      if(pot_is_cmplx)then
        do ist = 1, psib%nst
          forall (ip = 1:mesh%np)
            vpsib%states(ist)%X(psi)(ip, 1) = vpsib%states(ist)%X(psi)(ip, 1) + &
              (potential(ip, ispin)+ M_zI*Impotential(ip, ispin)) * psib%states(ist)%X(psi)(ip, 1) 
          end forall
        end do
      else
        !$omp parallel do private(ip)
        do ist = 1, psib%nst
          forall (ip = 1:mesh%np)
            vpsib%states(ist)%X(psi)(ip, 1) = vpsib%states(ist)%X(psi)(ip, 1) + &
              potential(ip, ispin) * psib%states(ist)%X(psi)(ip, 1)
          end forall
        end do
        !$omp end parallel do
      end if

      call profiling_count_operations((2*R_ADD*psib%nst)*mesh%np)
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

    case(SPINORS)
      !the spinor case is more complicated since it mixes the two components.
      if (pot_is_cmplx) then
        do ist = 1, psib%nst
          psi  => psib%states(ist)%X(psi)
          vpsi => vpsib%states(ist)%X(psi)
          
          do ip = 1, mesh%np
            pot(1:4) = potential(ip, 1:4) + M_zI * Impotential(ip, 1:4)
            vpsi(ip, 1) = vpsi(ip, 1) + pot(1)*psi(ip, 1) + &
                          (pot(3) + M_zI*pot(4))*psi(ip, 2)
            vpsi(ip, 2) = vpsi(ip, 2) + pot(2)*psi(ip, 2) + &
                          (pot(3) - M_zI*pot(4))*psi(ip, 1)
          end do
        end do
        
      else
        do ist = 1, psib%nst
          psi  => psib%states(ist)%X(psi)
          vpsi => vpsib%states(ist)%X(psi)

          forall(ip = 1:mesh%np)
            vpsi(ip, 1) = vpsi(ip, 1) + potential(ip, 1)*psi(ip, 1) + &
              (potential(ip, 3) + M_zI*potential(ip, 4))*psi(ip, 2)
            vpsi(ip, 2) = vpsi(ip, 2) + potential(ip, 2)*psi(ip, 2) + &
              (potential(ip, 3) - M_zI*potential(ip, 4))*psi(ip, 1)
          end forall
        end do
      end if
      call profiling_count_operations((6*R_ADD + 2*R_MUL)*mesh%np*psib%nst)

    end select

  end select

  call profiling_out(prof_vlpsi)
  POP_SUB(X(hamiltonian_base_local_sub))

end subroutine X(hamiltonian_base_local_sub)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_phase)(this, der, np, iqn, conjugate, psib, src)
  type(hamiltonian_base_t),              intent(in)    :: this
  type(derivatives_t),                   intent(in)    :: der
  integer,                               intent(in)    :: np
  integer,                               intent(in)    :: iqn
  logical,                               intent(in)    :: conjugate
  type(batch_t),                 target, intent(inout) :: psib
  type(batch_t),       optional, target, intent(in)    :: src

  integer :: ip, ii
  type(batch_t), pointer :: src_
  type(profile_t), save :: phase_prof
  CMPLX :: phase
  integer :: wgsize
  type(accel_kernel_t), save :: ker_phase

  PUSH_SUB(X(hamiltonian_base_phase))
  call profiling_in(phase_prof, "PBC_PHASE_APPLY")

  call profiling_count_operations(R_MUL*dble(np)*psib%nst_linear)

  ASSERT(np <= der%mesh%np_part)

  src_ => psib
  if(present(src)) src_ => src

  select case(batch_status(psib))
  case(BATCH_PACKED)

    if(conjugate) then

      !$omp parallel do private(ip, ii, phase)
      do ip = 1, np
        phase = conjg(this%phase(ip, iqn))
        do ii = 1, psib%nst_linear
          psib%pack%X(psi)(ii, ip) = phase*src_%pack%X(psi)(ii, ip)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(ip, ii, phase)
      do ip = 1, np
        phase = this%phase(ip, iqn)
        do ii = 1, psib%nst_linear
          psib%pack%X(psi)(ii, ip) = phase*src_%pack%X(psi)(ii, ip)
        end do
      end do
      !$omp end parallel do

    end if

  case(BATCH_NOT_PACKED)

    if(conjugate) then

      !$omp parallel private(ii, ip)
      do ii = 1, psib%nst_linear
        !$omp do
        do ip = 1, np
          psib%states_linear(ii)%X(psi)(ip) = conjg(this%phase(ip, iqn))*src_%states_linear(ii)%X(psi)(ip)
        end do
        !$omp end do nowait
      end do
      !$omp end parallel

    else
      !$omp parallel private(ii, ip)
      do ii = 1, psib%nst_linear
        !$omp do
        do ip = 1, np
          psib%states_linear(ii)%X(psi)(ip) = this%phase(ip, iqn)*src_%states_linear(ii)%X(psi)(ip)
        end do
        !$omp end do nowait
      end do
      !$omp end parallel

    end if

  case(BATCH_DEVICE_PACKED)
    call accel_kernel_start_call(ker_phase, 'phase.cl', 'phase_hamiltonian')

    if(conjugate) then
      call accel_set_kernel_arg(ker_phase, 0, 1_4)
    else
      call accel_set_kernel_arg(ker_phase, 0, 0_4)
    end if

    call accel_set_kernel_arg(ker_phase, 1, (iqn - this%buff_phase_qn_start)*der%mesh%np_part)
    call accel_set_kernel_arg(ker_phase, 2, np)
    call accel_set_kernel_arg(ker_phase, 3, this%buff_phase)
    call accel_set_kernel_arg(ker_phase, 4, src_%pack%buffer)
    call accel_set_kernel_arg(ker_phase, 5, log2(src_%pack%size(1)))
    call accel_set_kernel_arg(ker_phase, 6, psib%pack%buffer)
    call accel_set_kernel_arg(ker_phase, 7, log2(psib%pack%size(1)))

    wgsize = accel_kernel_workgroup_size(ker_phase)/psib%pack%size(1)

    call accel_kernel_run(ker_phase, (/psib%pack%size(1), pad(np, wgsize)/), (/psib%pack%size(1), wgsize/))

    call accel_finish()
  end select

  call profiling_out(phase_prof)
  POP_SUB(X(hamiltonian_base_phase))
end subroutine X(hamiltonian_base_phase)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_phase_spiral)(this, der, psib, ik)
  type(hamiltonian_base_t),              intent(in)    :: this
  type(derivatives_t),                   intent(in)    :: der
  type(batch_t),                         intent(inout) :: psib
  integer,                               intent(in)    :: ik

  integer :: ip, ii
  type(profile_t), save :: phase_prof

  PUSH_SUB(X(hamiltonian_base_phase_spiral))
  call profiling_in(phase_prof, "PBC_PHASE_SPIRAL")

  call profiling_count_operations(R_MUL*dble(der%mesh%np_part-der%mesh%np)*psib%nst_linear)


  ASSERT(.not. batch_status(psib) == BATCH_DEVICE_PACKED)

  ASSERT(der%boundaries%spiral)

  select case(batch_status(psib))
  case(BATCH_PACKED)

    !$omp parallel do private(ip, ii)
    do ip = der%mesh%np + 1, der%mesh%np_part
      do ii = 1, psib%nst_linear, 2
        if(this%spin(3,batch_linear_to_ist(psib, ii),ik)>0) then
          psib%pack%X(psi)(ii+1, ip) = psib%pack%X(psi)(ii+1, ip)*this%phase_spiral(ip, 1)
        else
          psib%pack%X(psi)(ii, ip) = psib%pack%X(psi)(ii, ip)*this%phase_spiral(ip, 2)
        end if
      end do
     end do
    !$omp end parallel do

  case(BATCH_NOT_PACKED)

    !$omp parallel private(ii, ip)
    do ii = 1, psib%nst_linear, 2
      if(this%spin(3,batch_linear_to_ist(psib, ii),ik)>0) then
        !$omp do
        do ip = der%mesh%np + 1, der%mesh%np_part
          psib%states_linear(ii+1)%X(psi)(ip) = psib%states_linear(ii+1)%X(psi)(ip)*this%phase_spiral(ip, 1)
        end do
        !$omp end do nowait
      else
        !$omp do
        do ip = der%mesh%np + 1, der%mesh%np_part
          psib%states_linear(ii)%X(psi)(ip) = psib%states_linear(ii)%X(psi)(ip)*this%phase_spiral(ip, 2)
        end do
        !$omp end do nowait
      end if
    end do
    !$omp end parallel

  end select

  call profiling_out(phase_prof)
  POP_SUB(X(hamiltonian_base_phase_spiral))
end subroutine X(hamiltonian_base_phase_spiral)


! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_rashba)(this, der, std, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(derivatives_t),         intent(in)    :: der
  type(states_dim_t),          intent(in)    :: std
  type(batch_t), target,       intent(in)    :: psib
  type(batch_t), target,       intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, allocatable :: psi(:, :), vpsi(:, :), grad(:, :, :)
  PUSH_SUB(X(hamiltonian_base_rashba))

  if(abs(this%rashba_coupling) < M_EPSILON) then
    POP_SUB(X(hamiltonian_base_rashba))
    return
  end if
  ASSERT(std%ispin == SPINORS)
  ASSERT(der%mesh%sb%dim == 2)

  SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:std%dim))
  SAFE_ALLOCATE(vpsi(1:der%mesh%np, 1:std%dim))
  SAFE_ALLOCATE(grad(1:der%mesh%np, 1:der%mesh%sb%dim, 1:std%dim))

  do ist = 1, psib%nst
    call batch_get_state(psib, ist, der%mesh%np_part, psi)
    call batch_get_state(vpsib, ist, der%mesh%np, vpsi)

    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
 
    if(allocated(this%vector_potential)) then
      forall(ip = 1:der%mesh%np)
        vpsi(ip, 1) = vpsi(ip, 1) + &
          (this%rashba_coupling) * (this%vector_potential(2, ip) + M_zI * this%vector_potential(1, ip)) * psi(ip, 2)
        vpsi(ip, 2) = vpsi(ip, 2) + &
          (this%rashba_coupling) * (this%vector_potential(2, ip) - M_zI * this%vector_potential(1, ip)) * psi(ip, 1)
      end forall
    end if

    forall(ip = 1:der%mesh%np)
      vpsi(ip, 1) = vpsi(ip, 1) - &
        this%rashba_coupling*( grad(ip, 1, 2) - M_zI*grad(ip, 2, 2) )
      vpsi(ip, 2) = vpsi(ip, 2) + &
        this%rashba_coupling*( grad(ip, 1, 1) + M_zI*grad(ip, 2, 1) )
    end forall

    call batch_set_state(vpsib, ist, der%mesh%np, vpsi)
  end do
  
  SAFE_DEALLOCATE_A(grad)
  SAFE_DEALLOCATE_A(vpsi)
  SAFE_DEALLOCATE_A(psi)
  
  POP_SUB(X(hamiltonian_base_rashba))
end subroutine X(hamiltonian_base_rashba)

! -----------------------------------------------------------------------------

subroutine X(hamiltonian_base_magnetic)(this, der, std, ep, ispin, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(derivatives_t),         intent(in)    :: der
  type(states_dim_t),          intent(in)    :: std
  type(epot_t),                intent(in)    :: ep
  integer,                     intent(in)    :: ispin
  type(batch_t), target,       intent(in)    :: psib
  type(batch_t), target,       intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, allocatable :: psi(:, :), vpsi(:, :), grad(:, :, :)
  FLOAT :: cc, b2, bb(1:MAX_DIM)
  CMPLX :: b12

  if(.not. hamiltonian_base_has_magnetic(this)) return

  call profiling_in(prof_magnetic, "MAGNETIC")
  PUSH_SUB(X(hamiltonian_base_magnetic))

  SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:std%dim))
  SAFE_ALLOCATE(vpsi(1:der%mesh%np, 1:std%dim))
  SAFE_ALLOCATE(grad(1:der%mesh%np, 1:der%mesh%sb%dim, 1:std%dim))

  do ist = 1, psib%nst
    call batch_get_state(psib, ist, der%mesh%np_part, psi)
    call batch_get_state(vpsib, ist, der%mesh%np, vpsi)

    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
 
    if(allocated(this%vector_potential)) then
      forall (idim = 1:std%dim, ip = 1:der%mesh%np)
        vpsi(ip, idim) = vpsi(ip, idim) + (M_HALF / this%mass) * &
          sum(this%vector_potential(1:der%mesh%sb%dim, ip)**2)*psi(ip, idim) &
          + (M_ONE / this%mass) * M_zI*dot_product(this%vector_potential(1:der%mesh%sb%dim, ip), grad(ip, 1:der%mesh%sb%dim, idim))
      end forall
    end if

    if(allocated(this%uniform_magnetic_field).and. std%ispin /= UNPOLARIZED) then
      ! Zeeman term
      cc = M_HALF/P_C*ep%gyromagnetic_ratio*M_HALF
      bb(1:max(der%mesh%sb%dim, 3)) = this%uniform_magnetic_field(1:max(der%mesh%sb%dim, 3))
      b2 = sqrt(sum(bb(1:max(der%mesh%sb%dim, 3))**2))
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

    call batch_set_state(vpsib, ist, der%mesh%np, vpsi)
  end do
  
  SAFE_DEALLOCATE_A(grad)
  SAFE_DEALLOCATE_A(vpsi)
  SAFE_DEALLOCATE_A(psi)
  
  POP_SUB(X(hamiltonian_base_magnetic))
  call profiling_out(prof_magnetic)
end subroutine X(hamiltonian_base_magnetic)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_start)(this, mesh, std, bnd, ik, psib, projection)
  type(hamiltonian_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_dim_t),               intent(in)    :: std
  type(boundaries_t),               intent(in)    :: bnd
  integer,                          intent(in)    :: ik
  type(batch_t),                    intent(in)    :: psib
  type(projection_t),               intent(out)   :: projection

  integer :: ist, ip, iproj, imat, nreal, iprojection
  integer :: npoints, nprojs, nst, maxnpoints
  integer, allocatable :: ind(:)
  type(projector_matrix_t), pointer :: pmat
  integer :: padnprojs, lnprojs, size, idim
  type(profile_t), save :: cl_prof
  type(accel_kernel_t), save, target :: ker_proj_bra, ker_proj_bra_phase
  type(accel_kernel_t), pointer :: kernel
  R_TYPE, allocatable :: lpsi(:, :)

  integer :: block_size
  
  if(.not. this%apply_projector_matrices) return

  call profiling_in(prof_vnlpsi_start, "VNLPSI_MAT_BRA")
  PUSH_SUB(X(hamiltonian_base_nlocal_start))

  nst = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

  if(batch_is_packed(psib) .and. accel_is_enabled()) then

    call accel_create_buffer(projection%buff_projection, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, &
      this%full_projection_size*psib%pack%size_real(1))

    call profiling_in(cl_prof, "CL_PROJ_BRA")
    ! only do this if we have some points of projector matrices
    if(this%max_npoints > 0) then

      if(allocated(this%projector_phases)) then
        call accel_kernel_start_call(ker_proj_bra_phase, 'projector.cl', 'projector_bra_phase')
        kernel => ker_proj_bra_phase
        size = psib%pack%size(1)
        ASSERT(R_TYPE_VAL == TYPE_CMPLX)
      else
        call accel_kernel_start_call(ker_proj_bra, 'projector.cl', 'projector_bra')
        kernel => ker_proj_bra
        size = psib%pack%size_real(1)
      end if

      call accel_set_kernel_arg(kernel, 0, this%nprojector_matrices)
      call accel_set_kernel_arg(kernel, 1, this%buff_offsets)
      call accel_set_kernel_arg(kernel, 2, this%buff_matrices)
      call accel_set_kernel_arg(kernel, 3, this%buff_maps)
      call accel_set_kernel_arg(kernel, 4, this%buff_scals)
      call accel_set_kernel_arg(kernel, 5, psib%pack%buffer)
      call accel_set_kernel_arg(kernel, 6, log2(size))
      call accel_set_kernel_arg(kernel, 7, projection%buff_projection)
      call accel_set_kernel_arg(kernel, 8, log2(size))

      if(allocated(this%projector_phases)) then
        call accel_set_kernel_arg(kernel, 9, this%buff_projector_phases)
        call accel_set_kernel_arg(kernel, 10, (ik - std%kpt%start)*this%total_points)
      end if

      padnprojs = pad_pow2(this%max_nprojs)
      lnprojs = min(accel_kernel_workgroup_size(kernel)/size, padnprojs)

      call accel_kernel_run(kernel, &
        (/size, padnprojs, this%nprojector_matrices/), (/size, lnprojs, 1/))

      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)

        npoints = pmat%npoints
        nprojs = pmat%nprojs

        call profiling_count_operations(nreal*nprojs*M_TWO*npoints + nst*nprojs)
      end do

      call accel_finish()
    end if

    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(projection%X(projection)(1:psib%pack%size_real(1), 1:this%full_projection_size))
      call accel_read_buffer(projection%buff_projection, &
        this%full_projection_size*psib%pack%size_real(1), projection%X(projection))
    end if

    call profiling_out(cl_prof)

    POP_SUB(X(hamiltonian_base_nlocal_start))
    call profiling_out(prof_vnlpsi_start)
    return
  end if

  ! This routine uses blocking to optimize cache usage. One block of
  ! |phi> is loaded in cache L1 and then then we calculate the dot
  ! product of it with the corresponding blocks of |psi_k>, next we
  ! load another block and do the same. This way we only have to load
  ! |psi> from the L2 or memory.
  block_size = hardware%X(block_size)


  SAFE_ALLOCATE(projection%X(projection)(1:nst, 1:this%full_projection_size))
  projection%X(projection) = M_ZERO

  SAFE_ALLOCATE(ind(1:this%nprojector_matrices))

  iprojection = 0
  maxnpoints = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    npoints = pmat%npoints
    maxnpoints = max(maxnpoints, npoints)
    nprojs = pmat%nprojs
    ind(imat) = iprojection
    iprojection = iprojection + nprojs
    call profiling_count_operations(nprojs*(R_ADD + R_MUL)*npoints + nst*nprojs)
  end do

  SAFE_ALLOCATE(lpsi(1:nst, 1:maxnpoints))

  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    iprojection = ind(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints == 0) cycle

    if(.not. allocated(this%projector_phases)) then
      if(batch_is_packed(psib)) then
        
        !$omp parallel do private(ist)
        do ip = 1, npoints
          forall(ist=1:nst)
            lpsi(ist, ip) = psib%pack%X(psi)(ist, pmat%map(ip))
          end forall
        end do
        
      else
        
        do ist = 1, nst
          !$omp parallel do
          do ip = 1, npoints
            lpsi(ist, ip) = psib%states_linear(ist)%X(psi)(pmat%map(ip))
          end do
        end do
        
      end if

    else
       if(.not. bnd%spiral) then 
        if(batch_is_packed(psib)) then
          !$omp parallel do private(ist)
          do ip = 1, npoints
            do ist = 1, nst
              lpsi(ist, ip) = psib%pack%X(psi)(ist, pmat%map(ip))*this%projector_phases(ip, imat, 1, ik)
            end do
          end do

        else

          do ist = 1, nst
            !$omp parallel do
            do ip = 1, npoints
              lpsi(ist, ip) = psib%states_linear(ist)%X(psi)(pmat%map(ip))*this%projector_phases(ip, imat, 1, ik)
            end do
          end do

        end if
      else
        if(batch_is_packed(psib)) then
         !$omp parallel do private(ist)
         do ip = 1, npoints
           do ist = 1, nst, 2
             if(this%spin(3,batch_linear_to_ist(psib, ist),ik)>0) then
               lpsi(ist, ip)   = psib%pack%X(psi)(ist, pmat%map(ip))*this%projector_phases(ip, imat, 1, ik)
               lpsi(ist+1, ip) = psib%pack%X(psi)(ist+1, pmat%map(ip))*this%projector_phases(ip, imat, 2, ik)
             else 
               lpsi(ist, ip)   = psib%pack%X(psi)(ist, pmat%map(ip))*this%projector_phases(ip, imat, 3, ik)
               lpsi(ist+1, ip) = psib%pack%X(psi)(ist+1, pmat%map(ip))*this%projector_phases(ip, imat, 1, ik)
             end if
           end do
         end do

       else

         do ist = 1, nst
           if(this%spin(3,batch_linear_to_ist(psib, ist),ik)>0 .and. batch_linear_to_idim(psib, ist)==2) then
             idim = 2
           else if(this%spin(3,batch_linear_to_ist(psib, ist),ik)<0 .and. batch_linear_to_idim(psib, ist)==1) then
             idim = 3
           else
             idim = 1
           end if
           !$omp parallel do
           do ip = 1, npoints
             lpsi(ist, ip) = psib%states_linear(ist)%X(psi)(pmat%map(ip))*this%projector_phases(ip, imat, idim, ik)
           end do
         end do

        end if

      end if
    end if

    call blas_gemm('N', 'N', nreal, nprojs, npoints, &
        M_ONE, lpsi(1, 1), nreal, pmat%projectors(1, 1), npoints, M_ZERO,  projection%X(projection)(1, iprojection + 1), nreal)
    call profiling_count_operations(nreal*nprojs*M_TWO*npoints)

    do iproj = 1, nprojs
      do ist = 1, nst
        projection%X(projection)(ist, iprojection + iproj) = projection%X(projection)(ist, iprojection + iproj)*pmat%scal(iproj)
       end do
    end do

  end do

  SAFE_DEALLOCATE_A(ind)
  SAFE_DEALLOCATE_A(lpsi)

  POP_SUB(X(hamiltonian_base_nlocal_start))
  call profiling_out(prof_vnlpsi_start)
end subroutine X(hamiltonian_base_nlocal_start)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_finish)(this, mesh, std, bnd, ik, projection, vpsib)
  type(hamiltonian_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_dim_t),               intent(in)    :: std
  type(boundaries_t),               intent(in)    :: bnd
  integer,                          intent(in)    :: ik
  type(projection_t),       target, intent(inout) :: projection
  type(batch_t),                    intent(inout) :: vpsib

  integer :: ist, ip, imat, nreal, iprojection
  integer :: npoints, nprojs, nst, idim
  CMPLX  :: phase, phase_pq, phase_mq
  R_TYPE, allocatable :: psi(:, :)
  type(projector_matrix_t), pointer :: pmat
  type(profile_t), save :: reduce_prof

  if(.not. this%apply_projector_matrices) return

  call profiling_in(prof_vnlpsi_finish, "VNLPSI_MAT_KET")
  PUSH_SUB(X(hamiltonian_base_nlocal_finish))

  nst = vpsib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

  ! reduce the projections
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, "VNLPSI_MAT_REDUCE")
    call comm_allreduce(mesh%vp%comm, projection%X(projection))
    call profiling_out(reduce_prof)
  end if

  if(batch_is_packed(vpsib) .and. accel_is_enabled()) then

    if(mesh%parallel_in_domains) then
      ! only do this if we have points of some projector matrices
      if(this%max_npoints > 0) then
        call accel_write_buffer(projection%buff_projection, &
          this%full_projection_size*vpsib%pack%size_real(1), projection%X(projection))
      end if
      SAFE_DEALLOCATE_A(projection%X(projection))
    end if

    call finish_opencl()
    call accel_release_buffer(projection%buff_projection)
    
    POP_SUB(X(hamiltonian_base_nlocal_finish))
    call profiling_out(prof_vnlpsi_finish)
    return
  end if

  ASSERT(allocated(projection%X(projection)))

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(allocated(pmat%mix)) then
      do ist = 1, nst
        projection%X(projection)(ist, iprojection + 1:iprojection + nprojs) = &
          matmul(pmat%mix(1:nprojs, 1:nprojs), projection%X(projection)(ist, iprojection + 1:iprojection + nprojs))
      end do
    end if
    
    if(npoints /=  0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints))

      ! Matrix-multiply again.
      ! the line below does: psi = matmul(projection, transpose(pmat%projectors))
      call blas_gemm('N', 'T', nreal, npoints, nprojs, &
        M_ONE, projection%X(projection)(1, iprojection + 1), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, psi(1, 1), nreal)
      
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints)

      call profiling_in(prof_scatter, "PROJ_MAT_SCATTER")

      if(.not. allocated(this%projector_phases)) then    
        ! and copy the points from the local buffer to its position
        if(batch_is_packed(vpsib)) then
          !$omp parallel do private(ip, ist) if(.not. this%projector_self_overlap)
          do ip = 1, npoints
            forall(ist = 1:nst)
              vpsib%pack%X(psi)(ist, pmat%map(ip)) = vpsib%pack%X(psi)(ist, pmat%map(ip)) + psi(ist, ip)
            end forall
          end do
          !$omp end parallel do
        else
          do ist = 1, nst
            !$omp parallel do if(.not. this%projector_self_overlap)
            do ip = 1, npoints
              vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) = vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) + psi(ist, ip)
            end do
            !$omp end parallel do
          end do
        end if
      else
        if(.not. bnd%spiral) then
          ! and copy the points from the local buffer to its position
          if(batch_is_packed(vpsib)) then
            !$omp parallel do private(ip, ist, phase)
            do ip = 1, npoints
              phase = conjg(this%projector_phases(ip, imat, 1, ik))
              forall(ist = 1:nst)
                vpsib%pack%X(psi)(ist, pmat%map(ip)) = vpsib%pack%X(psi)(ist, pmat%map(ip)) &
                            + psi(ist, ip)*phase
              end forall
            end do
            !$omp end parallel do
          else
            do ist = 1, nst
              !$omp parallel do
              do ip = 1, npoints
                vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) = vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) &
                    + psi(ist, ip)*conjg(this%projector_phases(ip, imat, 1, ik))
              end do
              !$omp end parallel do
            end do
          end if
        else
          ! and copy the points from the local buffer to its position
          if(batch_is_packed(vpsib)) then
            !$omp parallel do private(ip, ist) if(.not. this%projector_self_overlap)
            do ip = 1, npoints
              phase = conjg(this%projector_phases(ip, imat, 1, ik))
              phase_pq = conjg(this%projector_phases(ip, imat, 2, ik))
              phase_mq = conjg(this%projector_phases(ip, imat, 3, ik))
              do ist = 1, nst, 2
                if(this%spin(3,batch_linear_to_ist(vpsib, ist),ik)>0) then
                  vpsib%pack%X(psi)(ist, pmat%map(ip)) = vpsib%pack%X(psi)(ist, pmat%map(ip)) &
                              + psi(ist, ip)*phase
                  vpsib%pack%X(psi)(ist+1, pmat%map(ip)) = vpsib%pack%X(psi)(ist+1, pmat%map(ip)) &
                              + psi(ist+1, ip)*phase_pq
                else
                  vpsib%pack%X(psi)(ist, pmat%map(ip)) = vpsib%pack%X(psi)(ist, pmat%map(ip)) &
                              + psi(ist, ip)*phase_mq
                  vpsib%pack%X(psi)(ist+1, pmat%map(ip)) = vpsib%pack%X(psi)(ist+1, pmat%map(ip)) &
                              + psi(ist+1, ip)*phase
                end if
              end do
            end do
            !$omp end parallel do
          else
            do ist = 1, nst
              if(this%spin(3,batch_linear_to_ist(vpsib, ist),ik)>0 .and. batch_linear_to_idim(vpsib, ist)==2) then
                idim = 2
              else if(this%spin(3,batch_linear_to_ist(vpsib, ist),ik)<0 .and. batch_linear_to_idim(vpsib, ist)==1) then
                idim = 3
              else
                idim = 1
              end if
              !$omp parallel do if(.not. this%projector_self_overlap)
              do ip = 1, npoints
                vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) = vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) &
                    + psi(ist, ip)*conjg(this%projector_phases(ip, imat, idim, ik))
              end do
              !$omp end parallel do
            end do
          end if

        end if
      end if
      call profiling_count_operations(nst*npoints*R_ADD)
      call profiling_out(prof_scatter)
    end if
    
    SAFE_DEALLOCATE_A(psi)
    
    INCR(iprojection, nprojs)
  end do
  
  SAFE_DEALLOCATE_A(projection%X(projection))
  
  POP_SUB(X(hamiltonian_base_nlocal_finish))
  call profiling_out(prof_vnlpsi_finish)

contains

  subroutine finish_opencl()
    integer :: wgsize, imat, iregion, size, padnprojs, lnprojs
    type(profile_t), save :: cl_prof
    type(accel_kernel_t), save, target :: ker_proj_ket, ker_proj_ket_phase, ker_mix
    type(accel_kernel_t), pointer :: kernel
    type(accel_mem_t), pointer :: buff_proj
    
    PUSH_SUB(X(hamiltonian_base_nlocal_finish).finish_opencl)

    ! In this case we run one kernel per projector, since all write to
    ! the wave-function. Otherwise we would need to do atomic
    ! operations.

    call profiling_in(cl_prof, "CL_PROJ_KET")

    ! only do this if we have points of some projector matrices
    if(this%max_npoints > 0) then

      if(this%projector_mix) then

        SAFE_ALLOCATE(buff_proj)
        call accel_create_buffer(buff_proj, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, this%full_projection_size*vpsib%pack%size_real(1))

        call accel_kernel_start_call(ker_mix, 'projector.cl', 'projector_mix')
        
        call accel_set_kernel_arg(ker_mix, 0, this%nprojector_matrices)
        call accel_set_kernel_arg(ker_mix, 1, this%buff_offsets)
        call accel_set_kernel_arg(ker_mix, 2, this%buff_mix)
        call accel_set_kernel_arg(ker_mix, 3, projection%buff_projection)
        call accel_set_kernel_arg(ker_mix, 4, log2(vpsib%pack%size_real(1)))
        call accel_set_kernel_arg(ker_mix, 5, buff_proj)
        
        padnprojs = pad_pow2(this%max_nprojs)
        lnprojs = min(accel_kernel_workgroup_size(ker_mix)/vpsib%pack%size_real(1), padnprojs)
        
        call accel_kernel_run(ker_mix, &
          (/vpsib%pack%size_real(1), padnprojs, this%nprojector_matrices/), (/vpsib%pack%size_real(1), lnprojs, 1/))
        
        call accel_finish()

      else

        buff_proj => projection%buff_projection
        
      end if
      
      if(allocated(this%projector_phases)) then
        call accel_kernel_start_call(ker_proj_ket_phase, 'projector.cl', 'projector_ket_phase')
        kernel => ker_proj_ket_phase
        size = vpsib%pack%size(1)
        ASSERT(R_TYPE_VAL == TYPE_CMPLX)
      else
        call accel_kernel_start_call(ker_proj_ket, 'projector.cl', 'projector_ket')
        kernel => ker_proj_ket
        size = vpsib%pack%size_real(1)
      end if

      do iregion = 1, this%nregions
        
        call accel_set_kernel_arg(kernel, 0, this%nprojector_matrices)
        call accel_set_kernel_arg(kernel, 1, this%regions(iregion) - 1)
        call accel_set_kernel_arg(kernel, 2, this%buff_offsets)
        call accel_set_kernel_arg(kernel, 3, this%buff_matrices)
        call accel_set_kernel_arg(kernel, 4, this%buff_maps)
        call accel_set_kernel_arg(kernel, 5, buff_proj)
        call accel_set_kernel_arg(kernel, 6, log2(size))
        call accel_set_kernel_arg(kernel, 7, vpsib%pack%buffer)
        call accel_set_kernel_arg(kernel, 8, log2(size))

        if(allocated(this%projector_phases)) then
          call accel_set_kernel_arg(kernel, 9, this%buff_projector_phases)
          call accel_set_kernel_arg(kernel, 10, (ik - std%kpt%start)*this%total_points)
        end if

        wgsize = accel_kernel_workgroup_size(kernel)/size    

        call accel_kernel_run(kernel, &
          (/size, pad(this%max_npoints, wgsize), this%regions(iregion + 1) - this%regions(iregion)/), &
          (/size, wgsize, 1/))
        
        call accel_finish()
        
      end do
      
      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)
        npoints = pmat%npoints
        nprojs = pmat%nprojs
        call profiling_count_operations(nreal*nprojs*M_TWO*npoints)
        call profiling_count_operations(nst*npoints*R_ADD)
      end do

      call accel_finish()

      if(this%projector_mix) then
        call accel_release_buffer(buff_proj)
        SAFE_DEALLOCATE_P(buff_proj)
      end if
    end if
    
    call profiling_out(cl_prof)

    POP_SUB(X(hamiltonian_base_nlocal_finish).finish_opencl)
  end subroutine finish_opencl

end subroutine X(hamiltonian_base_nlocal_finish)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_force)(this, mesh, st, geo, iqn, ndim, psi1b, psi2b, force)
  type(hamiltonian_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_t),                   intent(in)    :: st
  type(geometry_t),                 intent(in)    :: geo
  integer,                          intent(in)    :: iqn
  integer,                          intent(in)    :: ndim
  type(batch_t),                    intent(in)    :: psi1b
  type(batch_t),                    intent(in)    :: psi2b(:)
  FLOAT,                            intent(inout) :: force(:, :)

  integer :: ii, ist, ip, iproj, imat, nreal, iprojection, iatom, idir
  integer :: npoints, nprojs, nst
  R_TYPE, allocatable :: psi(:, :, :), projs(:, :, :), ff(:)
  type(projector_matrix_t), pointer :: pmat

  if(.not. this%apply_projector_matrices) return
    
  call profiling_in(prof_matelement, "VNLPSI_MAT_ELEM")
  PUSH_SUB(X(hamiltonian_base_nlocal_force))

  ASSERT(psi1b%nst_linear == psi2b(1)%nst_linear)
  ASSERT(batch_status(psi1b) == batch_status(psi2b(1)))
  
  if(batch_is_packed(psi1b) .and. accel_is_enabled()) call messages_not_implemented('Accel non-local force')
  
  nst = psi1b%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif
  
  SAFE_ALLOCATE(projs(0:ndim, 1:nst, 1:this%full_projection_size))
  
  projs = CNST(0.0)

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints /= 0) then

      SAFE_ALLOCATE(psi(0:ndim, 1:nst, 1:npoints))
      
      call profiling_in(prof_matelement_gather, "PROJ_MAT_ELEM_GATHER")

      ! collect all the points we need in a continuous array
      if(batch_is_packed(psi1b)) then
        forall(ip = 1:npoints)
          forall(ist = 1:nst)
            psi(0, ist, ip) = psi1b%pack%X(psi)(ist, pmat%map(ip))
            forall(idir = 1:ndim) psi(idir, ist, ip) = psi2b(idir)%pack%X(psi)(ist, pmat%map(ip))
          end forall
        end forall
      else
        forall(ip = 1:npoints)
          forall(ist = 1:nst) 
            psi(0, ist, ip) = psi1b%states_linear(ist)%X(psi)(pmat%map(ip))
            forall(idir = 1:ndim) psi(idir, ist, ip) = psi2b(idir)%states_linear(ist)%X(psi)(pmat%map(ip))
          end forall
        end forall
      end if

      if(allocated(this%projector_phases)) then
        forall(ip = 1:npoints)
          forall(ist = 1:nst)
            forall(idir = 0:ndim)
              psi(idir, ist, ip) = this%projector_phases(ip, imat, batch_linear_to_idim(psi1b, ist), iqn)*psi(idir, ist, ip)
            end forall
          end forall
        end forall
      end if

      call profiling_out(prof_matelement_gather)
      
      ! Now matrix-multiply to calculate the projections. We can do all the matrix multiplications at once
      call blas_gemm('N', 'N', (ndim + 1)*nreal, nprojs, npoints, M_ONE, &
        psi(0, 1, 1), (ndim + 1)*nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, projs(0, 1, iprojection + 1), (ndim + 1)*nreal)
      
      call profiling_count_operations(nreal*(ndim + 1)*nprojs*M_TWO*npoints)

    else
      
      projs(0:ndim, 1:nst, iprojection + 1:iprojection + nprojs) = CNST(0.0)

    end if

    SAFE_DEALLOCATE_A(psi)

    INCR(iprojection, nprojs)

  end do

  if(mesh%parallel_in_domains) then
    call profiling_in(prof_matelement_reduce, "VNLPSI_MAT_ELEM_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, projs)
    call profiling_out(prof_matelement_reduce)
  end if
  
  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    
    npoints = pmat%npoints
    nprojs = pmat%nprojs
          
    iatom = this%projector_to_atom(imat)

    if(allocated(pmat%mix)) then
      do idir = 1, ndim
        do ist = 1, nst
          projs(idir, ist, iprojection + 1:iprojection + nprojs) = &
            matmul(pmat%mix(1:nprojs, 1:nprojs), projs(idir, ist, iprojection + 1:iprojection + nprojs))
        end do
      end do
    end if
    
    SAFE_ALLOCATE(ff(1:ndim))
    
    ff(1:ndim) = CNST(0.0)
    
    do ii = 1, psi1b%nst_linear
      ist = batch_linear_to_ist(psi1b, ii)
      if(st%d%kweights(iqn)*abs(st%occ(ist, iqn)) <= M_EPSILON) cycle
      do iproj = 1, nprojs
        do idir = 1, ndim
          ff(idir) = ff(idir) - CNST(2.0)*st%d%kweights(iqn)*st%occ(ist, iqn)*pmat%scal(iproj)*mesh%volume_element*&
            R_CONJ(projs(0, ii, iprojection + iproj))*projs(idir, ii, iprojection + iproj)
        end do
      end do
    end do
    
    force(1:ndim, iatom) = force(1:ndim, iatom) + ff(1:ndim)
    
    call profiling_count_operations((R_ADD + 2*R_MUL)*nst*ndim*nprojs)
    
    SAFE_DEALLOCATE_A(ff)

    INCR(iprojection, nprojs)

  end do

  SAFE_DEALLOCATE_A(projs)

  POP_SUB(X(hamiltonian_base_nlocal_force))
  call profiling_out(prof_matelement)
end subroutine X(hamiltonian_base_nlocal_force)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_position_commutator)(this, mesh, std, ik, psib, commpsib)
  type(hamiltonian_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_dim_t),               intent(in)    :: std
  integer,                          intent(in)    :: ik
  type(batch_t),                    intent(in)    :: psib
  type(batch_t),                    intent(inout) :: commpsib(:)

  integer :: ist, ip, iproj, imat, nreal, iprojection, idir
  integer :: npoints, nprojs, nst
  integer, allocatable :: ind(:)
  R_TYPE :: aa, bb, cc, dd
  R_TYPE, allocatable :: projections(:, :, :)
  R_TYPE, allocatable :: psi(:, :, :)
  CMPLX :: phase(2)
  type(projector_matrix_t), pointer :: pmat
  type(profile_t), save :: prof, reduce_prof
  integer :: wgsize, size

  if(.not. this%apply_projector_matrices) return

  PUSH_SUB(X(hamiltonian_base_nlocal_position_commutator))
  call profiling_in(prof, "COMMUTATOR")

  ASSERT(batch_is_packed(psib))

  nst = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

  if(batch_is_packed(psib) .and. accel_is_enabled()) then
    call X(commutator_opencl)()
    call profiling_out(prof)
    POP_SUB(X(hamiltonian_base_nlocal_position_commutator))
    return
  end if

  SAFE_ALLOCATE(projections(1:nst, 1:this%full_projection_size, 0:3))
  projections = M_ZERO

  SAFE_ALLOCATE(ind(1:this%nprojector_matrices))

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs
    ind(imat) = iprojection
    iprojection = iprojection + nprojs
    !    call profiling_count_operations(nprojs*(R_ADD + R_MUL)*npoints + nst*nprojs)
  end do

  !$omp parallel do private(imat, pmat, iprojection, npoints, nprojs, iproj, ist, aa, bb, cc, dd, ip, phase)
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    iprojection = ind(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints == 0) cycle

    if(.not. allocated(this%projector_phases)) then
      do iproj = 1, nprojs

        do ist = 1, nst
          aa = CNST(0.0)
          bb = CNST(0.0)
          cc = CNST(0.0)
          dd = CNST(0.0)
          do ip = 1, npoints
            aa = aa + pmat%projectors(ip, iproj)*psib%pack%X(psi)(ist, pmat%map(ip))
            bb = bb + pmat%projectors(ip, iproj)*pmat%position(1, ip)*psib%pack%X(psi)(ist, pmat%map(ip))
            cc = cc + pmat%projectors(ip, iproj)*pmat%position(2, ip)*psib%pack%X(psi)(ist, pmat%map(ip))
            dd = dd + pmat%projectors(ip, iproj)*pmat%position(3, ip)*psib%pack%X(psi)(ist, pmat%map(ip))
          end do
          projections(ist, iprojection + iproj, 0) = pmat%scal(iproj)*aa
          projections(ist, iprojection + iproj, 1) = pmat%scal(iproj)*bb
          projections(ist, iprojection + iproj, 2) = pmat%scal(iproj)*cc
          projections(ist, iprojection + iproj, 3) = pmat%scal(iproj)*dd
        end do

      end do

    else

      do iproj = 1, nprojs

        do ist = 1, nst
          aa = CNST(0.0)
          bb = CNST(0.0)
          cc = CNST(0.0)
          dd = CNST(0.0)
          do ip = 1, npoints
            phase(1) = this%projector_phases(ip, imat, batch_linear_to_idim(psib, ist), ik)
            aa = aa + pmat%projectors(ip, iproj)*psib%pack%X(psi)(ist, pmat%map(ip))*phase(1)
            bb = bb + pmat%projectors(ip, iproj)*pmat%position(1, ip)*psib%pack%X(psi)(ist, pmat%map(ip))*phase(1)
            cc = cc + pmat%projectors(ip, iproj)*pmat%position(2, ip)*psib%pack%X(psi)(ist, pmat%map(ip))*phase(1)
            dd = dd + pmat%projectors(ip, iproj)*pmat%position(3, ip)*psib%pack%X(psi)(ist, pmat%map(ip))*phase(1)
          end do
          projections(ist, iprojection + iproj, 0) = pmat%scal(iproj)*aa
          projections(ist, iprojection + iproj, 1) = pmat%scal(iproj)*bb
          projections(ist, iprojection + iproj, 2) = pmat%scal(iproj)*cc
          projections(ist, iprojection + iproj, 3) = pmat%scal(iproj)*dd
        end do

      end do
    end if

  end do

  ! reduce the projections
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, "COMMUTATOR_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, projections)
    call profiling_out(reduce_prof)
  end if

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(allocated(pmat%mix)) then
      do idir = 0, 3
        do ist = 1, nst
          projections(ist, iprojection + 1:iprojection + nprojs, idir) = &
            matmul(pmat%mix(1:nprojs, 1:nprojs), projections(ist, iprojection + 1:iprojection + nprojs, idir))
        end do
      end do
    end if
    
    if(npoints /=  0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints, 0:3))

      ! Matrix-multiply again.
      ! the line below does: psi = matmul(projection, transpose(pmat%projectors))
      call blas_gemm('N', 'T', nreal, npoints, nprojs, &
        M_ONE, projections(1, iprojection + 1, 0), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, psi(1, 1, 0), nreal)

      call blas_gemm('N', 'T', nreal, npoints, nprojs, &
        M_ONE, projections(1, iprojection + 1, 1), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, psi(1, 1, 1), nreal)

      call blas_gemm('N', 'T', nreal, npoints, nprojs, &
        M_ONE, projections(1, iprojection + 1, 2), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, psi(1, 1, 2), nreal)

      call blas_gemm('N', 'T', nreal, npoints, nprojs, &
        M_ONE, projections(1, iprojection + 1, 3), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, psi(1, 1, 3), nreal)
            
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints*4)

      if(allocated(this%projector_phases)) then
        do idir = 0, 3
          !$omp parallel do private(ip, ist, phase)
          do ip = 1, npoints
            phase(:) = conjg(this%projector_phases(ip, imat, :, ik))
            forall(ist = 1:nst)
              psi(ist, ip, idir) = phase(batch_linear_to_idim(psib, ist))*psi(ist, ip, idir)
            end forall
          end do
          !$omp end parallel do
        end do
      end if

      do idir = 1, 3
        do ip = 1, npoints
          forall(ist = 1:nst)
            commpsib(idir)%pack%X(psi)(ist, pmat%map(ip)) = commpsib(idir)%pack%X(psi)(ist, pmat%map(ip)) &
              - psi(ist, ip, idir) + pmat%position(idir, ip)*psi(ist, ip, 0)
          end forall
        end do
      end do
      
      call profiling_count_operations(nst*npoints*9*R_ADD)
    end if
    
    SAFE_DEALLOCATE_A(psi)
    
    INCR(iprojection, nprojs)
  end do

  SAFE_DEALLOCATE_A(ind)

  call profiling_out(prof)
  POP_SUB(X(hamiltonian_base_nlocal_position_commutator))

contains

  subroutine X(commutator_opencl)()
    type(accel_kernel_t), target, save :: ker_commutator_bra, ker_commutator_bra_phase, ker_mix
    type(accel_kernel_t), target, save :: ker_commutator_ket, ker_commutator_ket_phase
    type(accel_kernel_t), pointer :: kernel
    type(accel_mem_t), target :: buff_proj
    type(accel_mem_t), pointer :: buff_proj_copy
    integer :: padnprojs, lnprojs, iregion
    FLOAT, allocatable :: proj(:)

    padnprojs = pad_pow2(this%max_nprojs)

    call accel_create_buffer(buff_proj, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, 4*this%full_projection_size*psib%pack%size_real(1))

    if(allocated(this%projector_phases)) then
      call accel_kernel_start_call(ker_commutator_bra_phase, 'projector.cl', 'projector_commutator_bra_phase')
      kernel => ker_commutator_bra_phase
      size = psib%pack%size(1)
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
    else
      call accel_kernel_start_call(ker_commutator_bra, 'projector.cl', 'projector_commutator_bra')
      size = psib%pack%size_real(1)
      kernel => ker_commutator_bra
    end if
    
    call accel_set_kernel_arg(kernel,  0, this%nprojector_matrices)
    call accel_set_kernel_arg(kernel,  1, this%buff_offsets)
    call accel_set_kernel_arg(kernel,  2, this%buff_matrices)
    call accel_set_kernel_arg(kernel,  3, this%buff_maps)
    call accel_set_kernel_arg(kernel,  4, this%buff_scals)
    call accel_set_kernel_arg(kernel,  5, this%buff_position)
    call accel_set_kernel_arg(kernel,  6, psib%pack%buffer)
    call accel_set_kernel_arg(kernel,  7, log2(size))
    call accel_set_kernel_arg(kernel,  8, buff_proj)
    call accel_set_kernel_arg(kernel,  9, log2(size))

    if(allocated(this%projector_phases)) then
      call accel_set_kernel_arg(kernel, 10, this%buff_projector_phases)
      call accel_set_kernel_arg(kernel, 11, (ik - std%kpt%start)*this%total_points)
    end if
      
    lnprojs = min(accel_kernel_workgroup_size(kernel)/size, padnprojs)

    call accel_kernel_run(kernel, (/size, padnprojs, this%nprojector_matrices/), (/size, lnprojs, 1/))

    call accel_finish()

    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(proj(1:4*this%full_projection_size*psib%pack%size_real(1)))
      call accel_read_buffer(buff_proj, 4*this%full_projection_size*psib%pack%size_real(1), proj)
      call comm_allreduce(mesh%vp%comm, proj)
      call accel_write_buffer(buff_proj, 4*this%full_projection_size*psib%pack%size_real(1), proj)
    end if

    if(this%projector_mix) then

      SAFE_ALLOCATE(buff_proj_copy)
      
      call accel_create_buffer(buff_proj_copy, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, &
        4*this%full_projection_size*psib%pack%size_real(1))

      size = 4*psib%pack%size_real(1)
      
      call accel_kernel_start_call(ker_mix, 'projector.cl', 'projector_mix')
      
      call accel_set_kernel_arg(ker_mix, 0, this%nprojector_matrices)
      call accel_set_kernel_arg(ker_mix, 1, this%buff_offsets)
      call accel_set_kernel_arg(ker_mix, 2, this%buff_mix)
      call accel_set_kernel_arg(ker_mix, 3, buff_proj)
      call accel_set_kernel_arg(ker_mix, 4, log2(size))
      call accel_set_kernel_arg(ker_mix, 5, buff_proj_copy)
      
      padnprojs = pad_pow2(this%max_nprojs)
      lnprojs = min(accel_kernel_workgroup_size(ker_mix)/size, padnprojs)
      
      call accel_kernel_run(ker_mix, (/size, padnprojs, this%nprojector_matrices/), (/size, lnprojs, 1/))
      
      call accel_finish()

    else

      buff_proj_copy => buff_proj
      
    end if
    
    if(allocated(this%projector_phases)) then
      call accel_kernel_start_call(ker_commutator_ket_phase, 'projector.cl', 'projector_commutator_ket_phase')
      kernel => ker_commutator_ket_phase
      size = psib%pack%size(1)
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
    else
      call accel_kernel_start_call(ker_commutator_ket, 'projector.cl', 'projector_commutator_ket')
      kernel => ker_commutator_ket
      size = psib%pack%size_real(1)
    end if
    
    do iregion = 1, this%nregions
      
      call accel_set_kernel_arg(kernel,  0, this%nprojector_matrices)
      call accel_set_kernel_arg(kernel,  1, this%regions(iregion) - 1)
      call accel_set_kernel_arg(kernel,  2, this%buff_offsets)
      call accel_set_kernel_arg(kernel,  3, this%buff_matrices)
      call accel_set_kernel_arg(kernel,  4, this%buff_maps)
      call accel_set_kernel_arg(kernel,  5, this%buff_position)
      call accel_set_kernel_arg(kernel,  6, buff_proj_copy)
      call accel_set_kernel_arg(kernel,  7, log2(size))
      call accel_set_kernel_arg(kernel,  8, commpsib(1)%pack%buffer)
      call accel_set_kernel_arg(kernel,  9, commpsib(2)%pack%buffer)
      call accel_set_kernel_arg(kernel, 10, commpsib(3)%pack%buffer)
      call accel_set_kernel_arg(kernel, 11, log2(size))

      if(allocated(this%projector_phases)) then
        call accel_set_kernel_arg(kernel, 12, this%buff_projector_phases)
        call accel_set_kernel_arg(kernel, 13, (ik - std%kpt%start)*this%total_points)
      end if
      
      wgsize = accel_kernel_workgroup_size(kernel)/size    

      call accel_kernel_run(kernel, &
        (/size, pad(this%max_npoints, wgsize), this%regions(iregion + 1) - this%regions(iregion)/), &
        (/size, wgsize, 1/))
      
      call accel_finish()
      
    end do
    
    if(this%projector_mix) then
      call accel_release_buffer(buff_proj_copy)
      SAFE_ALLOCATE(buff_proj_copy)
    end if

    call accel_release_buffer(buff_proj)
    
  end subroutine X(commutator_opencl)
  
end subroutine X(hamiltonian_base_nlocal_position_commutator)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
