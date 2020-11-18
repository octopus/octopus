!! Copyright (C) 2009-2020 X. Andrade, N. Tancogne-Dejean, M. Lueders
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

subroutine X(hamiltonian_elec_base_local)(this, mesh, std, ispin, psib, vpsib)
  type(hamiltonian_elec_base_t),  intent(in)    :: this
  type(mesh_t),                   intent(in)    :: mesh
  type(states_elec_dim_t),        intent(in)    :: std
  integer,                        intent(in)    :: ispin
  type(wfs_elec_t),               intent(in)    :: psib
  type(wfs_elec_t),               intent(inout) :: vpsib

  PUSH_SUB(X(hamiltonian_elec_base_local))

  if(psib%status() == BATCH_DEVICE_PACKED) then
    ASSERT(.not. allocated(this%Impotential))
    call X(hamiltonian_elec_base_local_sub)(this%potential, mesh, std, ispin, &
      psib, vpsib, potential_opencl = this%potential_opencl)
  else
    if(allocated(this%Impotential)) then
      call X(hamiltonian_elec_base_local_sub)(this%potential, mesh, std, ispin, &
        psib, vpsib, Impotential = this%Impotential)
    else
      call X(hamiltonian_elec_base_local_sub)(this%potential, mesh, std, ispin, psib, vpsib)
    end if
  end if

  POP_SUB(X(hamiltonian_elec_base_local))
end subroutine X(hamiltonian_elec_base_local)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_local_sub)(potential, mesh, std, ispin, psib, vpsib, Impotential, potential_opencl)
  FLOAT,                        intent(in)    :: potential(:,:)
  type(mesh_t),                 intent(in)    :: mesh
  type(states_elec_dim_t),      intent(in)    :: std
  integer,                      intent(in)    :: ispin
  type(wfs_elec_t), target,     intent(in)    :: psib
  type(wfs_elec_t), target,     intent(inout) :: vpsib
  FLOAT, optional,              intent(in)    :: Impotential(:,:)
  type(accel_mem_t),  optional, intent(in)    :: potential_opencl

  integer :: ist, ip, dim2, dim3
#ifdef R_TCOMPLEX
  R_TYPE :: psi1, psi2
  FLOAT  :: Imvv
  CMPLX  :: pot(1:4)
  CMPLX, pointer :: psi(:, :), vpsi(:, :)
#endif
  FLOAT   :: vv
  logical :: pot_is_cmplx
  integer :: pnp, localsize

  call profiling_in(prof_vlpsi, TOSTRING(X(VLPSI)))
  PUSH_SUB(X(hamiltonian_elec_base_local_sub))

  pot_is_cmplx = .false.
  if(present(Impotential)) pot_is_cmplx = .true.

  call psib%check_compatibility_with(vpsib)

  select case(psib%status())
  case(BATCH_DEVICE_PACKED)
    ASSERT(.not. pot_is_cmplx) ! not implemented

    pnp = accel_padded_size(mesh%np)

    select case(std%ispin)

    case(UNPOLARIZED, SPIN_POLARIZED)
      call accel_set_kernel_arg(kernel_vpsi, 0, pnp*(ispin - 1))
      call accel_set_kernel_arg(kernel_vpsi, 1, mesh%np)
      call accel_set_kernel_arg(kernel_vpsi, 2, potential_opencl)
      call accel_set_kernel_arg(kernel_vpsi, 3, psib%ff_device)
      call accel_set_kernel_arg(kernel_vpsi, 4, log2(psib%pack_size_real(1)))
      call accel_set_kernel_arg(kernel_vpsi, 5, vpsib%ff_device)
      call accel_set_kernel_arg(kernel_vpsi, 6, log2(vpsib%pack_size_real(1)))

      localsize = accel_kernel_workgroup_size(kernel_vpsi)/psib%pack_size_real(1)

      dim3 = mesh%np/(accel_max_size_per_dim(2)*localsize) + 1
      dim2 = min(accel_max_size_per_dim(2)*localsize, pad(mesh%np, localsize))

      call accel_kernel_run(kernel_vpsi, (/psib%pack_size_real(1), dim2, dim3/), (/psib%pack_size_real(1), localsize, 1/))

    case(SPINORS)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 0, mesh%np)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 1, potential_opencl)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 2, pnp)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 3, psib%ff_device)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 4, psib%pack_size(1))
      call accel_set_kernel_arg(kernel_vpsi_spinors, 5, vpsib%ff_device)
      call accel_set_kernel_arg(kernel_vpsi_spinors, 6, vpsib%pack_size(1))

      call accel_kernel_run(kernel_vpsi_spinors, (/psib%pack_size(1)/2, pnp/), &
        (/psib%pack_size(1)/2, 2*accel_max_workgroup_size()/psib%pack_size(1)/))

    end select

    call accel_finish()

    call profiling_count_operations((R_MUL*psib%nst_linear)*mesh%np)
    call profiling_count_transfers(mesh%np, M_ONE)
    call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

  case(BATCH_PACKED)

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      if(pot_is_cmplx)then
#ifdef R_TCOMPLEX
        do ip = 1, mesh%np
          vv = potential(ip, ispin)
          Imvv = Impotential(ip, ispin)
          do ist = 1, psib%nst_linear
            vpsib%zff_pack(ist, ip) = vpsib%zff_pack(ist, ip) + (vv + M_zI*Imvv)*psib%zff_pack(ist, ip)
          end do
        end do
        call profiling_count_operations(2*((R_ADD+R_MUL)*psib%nst_linear)*mesh%np)
#else
        ! Complex potential can only be applied to complex batches
        ASSERT(.false.)
#endif
      else
        !$omp parallel do private(vv, ist)
        do ip = 1, mesh%np
          vv = potential(ip, ispin)
          do ist = 1, psib%nst_linear
            vpsib%X(ff_pack)(ist, ip) = vpsib%X(ff_pack)(ist, ip) + vv*psib%X(ff_pack)(ist, ip)
          end do
        end do
        !$omp end parallel do
        call profiling_count_operations(((R_ADD+R_MUL)*psib%nst_linear)*mesh%np)
      end if
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst_linear, R_TOTYPE(M_ONE))

    case(SPINORS)
#ifdef R_TCOMPLEX
      ASSERT(mod(psib%nst_linear, 2) == 0)
      !the spinor case is more complicated since it mixes the two components.
      if(pot_is_cmplx)then
        !$omp parallel do private(psi1, psi2, ist, pot)
        do ip = 1, mesh%np
          do ist = 1, psib%nst_linear, 2
            psi1 = psib%zff_pack(ist    , ip)
            psi2 = psib%zff_pack(ist + 1, ip)
            pot(1:4) = potential(ip, 1:4) + M_zI * Impotential(ip, 1:4)
            vpsib%zff_pack(ist    , ip) = vpsib%zff_pack(ist    , ip) + &
                   pot(1)*psi1 + (pot(3) + M_zI*pot(4))*psi2
            vpsib%zff_pack(ist + 1, ip) = vpsib%zff_pack(ist + 1, ip) + &
                   pot(2)*psi2 + (pot(3) - M_zI*pot(4))*psi1            
          end do
        end do
        !$omp end parallel do
        call profiling_count_operations((7*R_ADD + 7*R_MUL)*mesh%np*psib%nst)
      else
        !$omp parallel do private(psi1, psi2, ist)
        do ip = 1, mesh%np
          do ist = 1, psib%nst_linear, 2
            psi1 = psib%zff_pack(ist    , ip)
            psi2 = psib%zff_pack(ist + 1, ip)
            vpsib%zff_pack(ist    , ip) = vpsib%zff_pack(ist    , ip) + &
              potential(ip, 1)*psi1 + (potential(ip, 3) + M_zI*potential(ip, 4))*psi2
            vpsib%zff_pack(ist + 1, ip) = vpsib%zff_pack(ist + 1, ip) + &
              potential(ip, 2)*psi2 + (potential(ip, 3) - M_zI*potential(ip, 4))*psi1            
          end do
        end do
        !$omp end parallel do
        call profiling_count_operations((6*R_ADD + 6*R_MUL)*mesh%np*psib%nst)
      end if
#else
      ! Spinors always imply complex batches
      ASSERT(.false.)
#endif      
    end select

  case(BATCH_NOT_PACKED)

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      if(pot_is_cmplx)then
#ifdef R_TCOMPLEX
        do ist = 1, psib%nst
          do ip = 1, mesh%np
            vpsib%X(ff)(ip, 1, ist) = vpsib%X(ff)(ip, 1, ist) + &
              (potential(ip, ispin)+ M_zI*Impotential(ip, ispin)) * psib%X(ff)(ip, 1, ist)
          end do
        end do
        call profiling_count_operations(2*((R_ADD+R_MUL)*psib%nst)*mesh%np)
#else
        ! Complex potential can only be applied to complex batches
        ASSERT(.false.)
#endif
      else
        !$omp parallel do private(ip)
        do ist = 1, psib%nst
          do ip = 1, mesh%np
            vpsib%X(ff)(ip, 1, ist) = vpsib%X(ff)(ip, 1, ist) + &
              potential(ip, ispin) * psib%X(ff)(ip, 1, ist)
          end do
        end do
        !$omp end parallel do

        call profiling_count_operations(((R_ADD+R_MUL)*psib%nst)*mesh%np)
      end if

      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

    case(SPINORS)
#ifdef R_TCOMPLEX
      !the spinor case is more complicated since it mixes the two components.
      if (pot_is_cmplx) then
        do ist = 1, psib%nst
          psi  => psib%zff(:, :, ist)
          vpsi => vpsib%zff(:, :, ist)
          
          do ip = 1, mesh%np
            pot(1:4) = potential(ip, 1:4) + M_zI * Impotential(ip, 1:4)
            vpsi(ip, 1) = vpsi(ip, 1) + pot(1)*psi(ip, 1) + &
                          (pot(3) + M_zI*pot(4))*psi(ip, 2)
            vpsi(ip, 2) = vpsi(ip, 2) + pot(2)*psi(ip, 2) + &
                          (pot(3) - M_zI*pot(4))*psi(ip, 1)
          end do
        end do
        call profiling_count_operations((7*R_ADD + 7*R_MUL)*mesh%np*psib%nst)
        
      else
        do ist = 1, psib%nst
          psi  => psib%zff(:, :, ist)
          vpsi => vpsib%zff(:, :, ist)

          do ip = 1, mesh%np
            vpsi(ip, 1) = vpsi(ip, 1) + potential(ip, 1)*psi(ip, 1) + &
              (potential(ip, 3) + M_zI*potential(ip, 4))*psi(ip, 2)
            vpsi(ip, 2) = vpsi(ip, 2) + potential(ip, 2)*psi(ip, 2) + &
              (potential(ip, 3) - M_zI*potential(ip, 4))*psi(ip, 1)
          end do
        end do
        call profiling_count_operations((6*R_ADD + 6*R_MUL)*mesh%np*psib%nst)
      end if
#else
      ! Spinors always imply complex batches
      ASSERT(.false.)
#endif  
    end select

  end select

  call profiling_out(prof_vlpsi)
  POP_SUB(X(hamiltonian_elec_base_local_sub))

end subroutine X(hamiltonian_elec_base_local_sub)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_rashba)(this, mesh, der, std, psib, vpsib)
  type(hamiltonian_elec_base_t),  intent(in)    :: this
  type(mesh_t),                   intent(in)    :: mesh
  type(derivatives_t),            intent(in)    :: der
  type(states_elec_dim_t),        intent(in)    :: std
  type(wfs_elec_t), target,       intent(in)    :: psib
  type(wfs_elec_t), target,       intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, allocatable :: psi(:, :), vpsi(:, :), grad(:, :, :)
  PUSH_SUB(X(hamiltonian_elec_base_rashba))

  if(abs(this%rashba_coupling) < M_EPSILON) then
    POP_SUB(X(hamiltonian_elec_base_rashba))
    return
  end if
  ASSERT(std%ispin == SPINORS)
  ASSERT(mesh%sb%dim == 2)

  SAFE_ALLOCATE(psi(1:mesh%np_part, 1:std%dim))
  SAFE_ALLOCATE(vpsi(1:mesh%np, 1:std%dim))
  SAFE_ALLOCATE(grad(1:mesh%np, 1:mesh%sb%dim, 1:std%dim))

  do ist = 1, psib%nst
    call batch_get_state(psib, ist, mesh%np_part, psi)
    call batch_get_state(vpsib, ist, mesh%np, vpsi)

    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do

    if(allocated(this%vector_potential)) then
      do ip = 1, mesh%np
        vpsi(ip, 1) = vpsi(ip, 1) + &
          (this%rashba_coupling) * (this%vector_potential(2, ip) + M_zI * this%vector_potential(1, ip)) * psi(ip, 2)
        vpsi(ip, 2) = vpsi(ip, 2) + &
          (this%rashba_coupling) * (this%vector_potential(2, ip) - M_zI * this%vector_potential(1, ip)) * psi(ip, 1)
      end do
    end if

    do ip = 1, mesh%np
      vpsi(ip, 1) = vpsi(ip, 1) - &
        this%rashba_coupling*( grad(ip, 1, 2) - M_zI*grad(ip, 2, 2) )
      vpsi(ip, 2) = vpsi(ip, 2) + &
        this%rashba_coupling*( grad(ip, 1, 1) + M_zI*grad(ip, 2, 1) )
    end do

    call batch_set_state(vpsib, ist, mesh%np, vpsi)
  end do

  SAFE_DEALLOCATE_A(grad)
  SAFE_DEALLOCATE_A(vpsi)
  SAFE_DEALLOCATE_A(psi)
  
  POP_SUB(X(hamiltonian_elec_base_rashba))
end subroutine X(hamiltonian_elec_base_rashba)

! -----------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_magnetic)(this, mesh, der, std, ep, ispin, psib, vpsib)
  type(hamiltonian_elec_base_t),  intent(in)    :: this
  type(mesh_t),                   intent(in)    :: mesh
  type(derivatives_t),            intent(in)    :: der
  type(states_elec_dim_t),        intent(in)    :: std
  type(epot_t),                   intent(in)    :: ep
  integer,                        intent(in)    :: ispin
  type(wfs_elec_t), target,       intent(in)    :: psib
  type(wfs_elec_t), target,       intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, allocatable :: psi(:, :), vpsi(:, :), grad(:, :, :)
  FLOAT :: cc, b2, bb(1:MAX_DIM)
  CMPLX :: b12

  if(.not. hamiltonian_elec_base_has_magnetic(this)) return

  call profiling_in(prof_magnetic, TOSTRING(X(MAGNETIC)))
  PUSH_SUB(X(hamiltonian_elec_base_magnetic))

  SAFE_ALLOCATE(psi(1:mesh%np_part, 1:std%dim))
  SAFE_ALLOCATE(vpsi(1:mesh%np, 1:std%dim))
  SAFE_ALLOCATE(grad(1:mesh%np, 1:mesh%sb%dim, 1:std%dim))

  do ist = 1, psib%nst
    call batch_get_state(psib, ist, mesh%np_part, psi)
    call batch_get_state(vpsib, ist, mesh%np, vpsi)

    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do

    if(allocated(this%vector_potential)) then
      do idim = 1, std%dim
        do ip = 1, mesh%np
          vpsi(ip, idim) = vpsi(ip, idim) + (M_HALF / this%mass) * &
            sum(this%vector_potential(1:mesh%sb%dim, ip)**2)*psi(ip, idim) &
            + (M_ONE / this%mass) * M_zI*dot_product(this%vector_potential(1:mesh%sb%dim, ip), grad(ip, 1:mesh%sb%dim, idim))
        end do
      end do
    end if

    if(allocated(this%uniform_magnetic_field).and. std%ispin /= UNPOLARIZED) then
      ! Zeeman term
      cc = M_HALF/P_C*ep%gyromagnetic_ratio*M_HALF
      bb(1:max(mesh%sb%dim, 3)) = this%uniform_magnetic_field(1:max(mesh%sb%dim, 3))
      b2 = sqrt(sum(bb(1:max(mesh%sb%dim, 3))**2))
      b12 = bb(1) - M_ZI*bb(2)

      select case (std%ispin)
      case (SPIN_POLARIZED)
        if(is_spin_down(ispin)) cc = -cc

        do ip = 1, mesh%np
          vpsi(ip, 1) = vpsi(ip, 1) + cc*b2*psi(ip, 1)
        end do

      case (SPINORS)
        do ip = 1, mesh%np
          vpsi(ip, 1) = vpsi(ip, 1) + cc*(bb(3)*psi(ip, 1) + b12*psi(ip, 2))
          vpsi(ip, 2) = vpsi(ip, 2) + cc*(-bb(3)*psi(ip, 2) + conjg(b12)*psi(ip, 1))
        end do

      end select
    end if

    call batch_set_state(vpsib, ist, mesh%np, vpsi)
  end do

  SAFE_DEALLOCATE_A(grad)
  SAFE_DEALLOCATE_A(vpsi)
  SAFE_DEALLOCATE_A(psi)
  
  POP_SUB(X(hamiltonian_elec_base_magnetic))
  call profiling_out(prof_magnetic)
end subroutine X(hamiltonian_elec_base_magnetic)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_nlocal_start)(this, mesh, std, bnd, psib, projection)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_elec_dim_t),          intent(in)    :: std
  type(boundaries_t),               intent(in)    :: bnd
  type(wfs_elec_t),                 intent(in)    :: psib
  type(projection_t),               intent(out)   :: projection

  integer :: ist, ip, iproj, imat, nreal, iprojection
  integer :: npoints, nprojs, nst_linear, maxnpoints
  integer, allocatable :: ind(:)
  type(projector_matrix_t), pointer :: pmat
  integer :: padnprojs, lnprojs, size, nphase
  type(profile_t), save :: cl_prof
  type(accel_kernel_t), save, target :: ker_proj_bra, ker_proj_bra_phase, ker_proj_bra_phase_spiral
  type(accel_kernel_t), pointer :: kernel
  integer, allocatable :: spin_to_phase(:)
  R_TYPE, allocatable :: lpsi(:, :)
#ifdef R_TCOMPLEX
  integer :: iphase
  CMPLX, allocatable :: tmp_proj(:, :)
#endif

  integer :: block_size
  integer :: size_unfolded

  if(.not. this%has_non_local_potential) return
  
  ASSERT(this%apply_projector_matrices)

  call profiling_in(prof_vnlpsi_start, TOSTRING(X(VNLPSI_MAT_BRA)))
  PUSH_SUB(X(hamiltonian_elec_base_nlocal_start))

  nst_linear = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst_linear
#else
  nreal = nst_linear
#endif
  nphase = 1

  if(psib%has_phase) then
    ASSERT(allocated(this%projector_phases))
  end if

  if(psib%status() == BATCH_DEVICE_PACKED) then

    call accel_create_buffer(projection%buff_projection, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, &
      this%full_projection_size*psib%pack_size_real(1))

    call profiling_in(cl_prof, TOSTRING(X(CL_PROJ_BRA)))
    ! only do this if we have some points of projector matrices
    if(this%max_npoints > 0) then

      if(allocated(this%projector_phases)) then

        if(bnd%spiral) then

          nphase = 3

          SAFE_ALLOCATE(spin_to_phase(psib%pack_size(1)))
          call accel_create_buffer(projection%buff_spin_to_phase, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, psib%pack_size(1))
    
          do ist = 1, nst_linear
            if(this%spin(3, psib%linear_to_ist(ist), psib%ik) > 0 .and. psib%linear_to_idim(ist)==2) then
              spin_to_phase(ist) = 1
            else if(this%spin(3, psib%linear_to_ist(ist), psib%ik) < 0 .and. psib%linear_to_idim(ist)==1) then
              spin_to_phase(ist) = 2
            else
              spin_to_phase(ist) = 0
            end if 
          end do
          ! This might not be necessary:
          do ist = nst_linear+1, psib%pack_size(1)
            spin_to_phase(ist) = 0
          end do

          call accel_write_buffer(projection%buff_spin_to_phase, psib%pack_size(1), spin_to_phase)

          call accel_kernel_start_call(ker_proj_bra_phase_spiral, 'projector.cl', 'projector_bra_phase_spiral')
          kernel => ker_proj_bra_phase_spiral
          SAFE_DEALLOCATE_A(spin_to_phase)
        else 
          call accel_kernel_start_call(ker_proj_bra_phase, 'projector.cl', 'projector_bra_phase')
          kernel => ker_proj_bra_phase
        end if
        size = psib%pack_size(1)
        ASSERT(R_TYPE_VAL == TYPE_CMPLX)  

      else
        call accel_kernel_start_call(ker_proj_bra, 'projector.cl', 'projector_bra')
        kernel => ker_proj_bra
        size = psib%pack_size_real(1)
      end if

      call accel_set_kernel_arg(kernel, 0, this%nprojector_matrices)
      call accel_set_kernel_arg(kernel, 1, this%buff_offsets)
      call accel_set_kernel_arg(kernel, 2, this%buff_matrices)
      call accel_set_kernel_arg(kernel, 3, this%buff_maps)
      call accel_set_kernel_arg(kernel, 4, this%buff_scals)
      call accel_set_kernel_arg(kernel, 5, psib%ff_device)
      call accel_set_kernel_arg(kernel, 6, log2(size))
      call accel_set_kernel_arg(kernel, 7, projection%buff_projection)
      call accel_set_kernel_arg(kernel, 8, log2(size))

      if(allocated(this%projector_phases)) then
        call accel_set_kernel_arg(kernel, 9, this%buff_projector_phases)
        ! Note: we need to use this%nphase, as the kernel might be called with spiral=false, but 
        !       the phases been built with spiralBC=true
        call accel_set_kernel_arg(kernel, 10, (psib%ik - std%kpt%start)*this%total_points*this%nphase)
        if(bnd%spiral) then
          call accel_set_kernel_arg(kernel, 11, projection%buff_spin_to_phase)
        end if
      end if

      ! In case of CUDA we use an optimized kernel, in which the loop over npoints is broken
      ! further into chunks, in order to parallelize over the threads within a warp.
      ! Therefore we need to launch warp_size * size kernels. The size of each block needs to 
      ! have multiples of warp_size as x-dimension.

      size_unfolded = size * accel%warp_size
      padnprojs = pad_pow2(this%max_nprojs)
      lnprojs = min(accel_kernel_workgroup_size(kernel)/accel%warp_size, padnprojs)

      call accel_kernel_run(kernel, &
        (/size_unfolded, padnprojs, this%nprojector_matrices/), (/accel%warp_size, lnprojs, 1/))

      do imat = 1, this%nprojector_matrices
        pmat => this%projector_matrices(imat)

        npoints = pmat%npoints
        nprojs = pmat%nprojs

        !! update number of operations for nphase !! 
        call profiling_count_operations(nreal*nprojs*M_TWO*npoints + nst_linear*nprojs)
      end do

      call accel_finish()
    end if

    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(projection%X(projection)(1:psib%pack_size_real(1), 1:this%full_projection_size))
      call accel_read_buffer(projection%buff_projection, &
        this%full_projection_size*psib%pack_size_real(1), projection%X(projection))
    end if


    call profiling_out(cl_prof)

    POP_SUB(X(hamiltonian_elec_base_nlocal_start))
    call profiling_out(prof_vnlpsi_start)
    return
  end if

  ! This routine uses blocking to optimize cache usage. One block of
  ! |phi> is loaded in cache L1 and then then we calculate the dot
  ! product of it with the corresponding blocks of |psi_k>, next we
  ! load another block and do the same. This way we only have to load
  ! |psi> from the L2 or memory.
  block_size = hardware%X(block_size)


  SAFE_ALLOCATE(projection%X(projection)(1:nst_linear, 1:this%full_projection_size))
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
    call profiling_count_operations(nprojs*(R_ADD + R_MUL)*npoints + nst_linear*nprojs)
    if(allocated(this%projector_phases)) then
      call profiling_count_operations(R_MUL*npoints*nst_linear)
    end if
  end do

  SAFE_ALLOCATE(lpsi(1:nst_linear, 1:maxnpoints))

  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    iprojection = ind(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints == 0) cycle

    if(.not. allocated(this%projector_phases)) then
      if(psib%status() == BATCH_PACKED) then
        
        !$omp parallel do private(ist, ip)
        do ip = 1, npoints
          do ist= 1, nst_linear
            lpsi(ist, ip) = psib%X(ff_pack)(ist, pmat%map(ip))
          end do
        end do
        
      else
        
        do ist = 1, nst_linear
          !$omp parallel do
          do ip = 1, npoints
            lpsi(ist, ip) = psib%X(ff_linear)(pmat%map(ip), ist)
          end do
        end do
        
      end if

    else
      if(.not. bnd%spiral) then 
        if(psib%status() == BATCH_PACKED) then
#ifdef R_TCOMPLEX
          !$omp parallel do private(ist)
          do ip = 1, npoints
            do ist = 1, nst_linear
              lpsi(ist, ip) = psib%zff_pack(ist, pmat%map(ip))*this%projector_phases(ip, 1, imat, psib%ik)
            end do
          end do
#else
        ! Phases not allowed for real batches
        ASSERT(.false.)
#endif
        else
#ifdef R_TCOMPLEX
          do ist = 1, nst_linear
            !$omp parallel do
            do ip = 1, npoints
              lpsi(ist, ip) = psib%zff_linear(pmat%map(ip), ist)*this%projector_phases(ip, 1, imat, psib%ik)
            end do
          end do
#else
          ! Phases not allowed for real batches
          ASSERT(.false.)
#endif
        end if
      else
        if(psib%status() == BATCH_PACKED) then
#ifdef R_TCOMPLEX
         !$omp parallel do private(ist)
         do ip = 1, npoints
           do ist = 1, nst_linear, 2
             if(this%spin(3,psib%linear_to_ist(ist), psib%ik)>0) then
               lpsi(ist, ip)   = psib%zff_pack(ist,   pmat%map(ip))*this%projector_phases(ip, 1, imat, psib%ik)
               lpsi(ist+1, ip) = psib%zff_pack(ist+1, pmat%map(ip))*this%projector_phases(ip, 2, imat, psib%ik)
             else 
               lpsi(ist, ip)   = psib%zff_pack(ist,   pmat%map(ip))*this%projector_phases(ip, 3, imat, psib%ik)
               lpsi(ist+1, ip) = psib%zff_pack(ist+1, pmat%map(ip))*this%projector_phases(ip, 1, imat, psib%ik)
             end if
           end do
         end do
#else
         ! Phases not allowed for real batches
         ASSERT(.false.)
#endif
       else
#ifdef R_TCOMPLEX
         do ist = 1, nst_linear
           if(this%spin(3, psib%linear_to_ist(ist), psib%ik) > 0 .and. psib%linear_to_idim(ist)==2) then
             iphase = 2
           else if(this%spin(3, psib%linear_to_ist(ist), psib%ik) < 0 .and. psib%linear_to_idim(ist)==1) then
             iphase = 3
           else
             iphase = 1
           end if
           !$omp parallel do
           do ip = 1, npoints
             lpsi(ist, ip) = psib%zff_linear(pmat%map(ip), ist)*this%projector_phases(ip, iphase, imat, psib%ik)
           end do
         end do
#else
         ! Phases not allowed for real batches
         ASSERT(.false.)
#endif
        end if

      end if
    end if

    if(pmat%is_cmplx) then
#ifdef R_TCOMPLEX
      SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:nst_linear))
      call blas_gemm('C', 'T', nprojs, nst_linear, npoints, &
          M_z1, pmat%zprojectors(1, 1), npoints, lpsi(1, 1), nst_linear, M_z0, tmp_proj(1,1), nprojs)
      !$omp parallel do private(iproj, ist)
      do iproj = 1, nprojs
        do ist = 1, nst_linear
          projection%X(projection)(ist, iprojection + iproj) = tmp_proj(iproj, ist)
        end do
      end do
      SAFE_DEALLOCATE_A(tmp_proj)
      call profiling_count_operations(nst_linear*nprojs*M_TWO*npoints)
#else
      ! Complex projection matrix not allowed for real batches
      ASSERT(.false.)
#endif
    else
      call blas_gemm('N', 'N', nreal, nprojs, npoints, &
          M_ONE, lpsi(1, 1), nreal, pmat%dprojectors(1, 1), npoints, M_ZERO,  projection%X(projection)(1, iprojection + 1), nreal)
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints)
    end if

    !$omp parallel do private(iproj, ist)
    do iproj = 1, nprojs
      do ist = 1, nst_linear
        projection%X(projection)(ist, iprojection + iproj) = projection%X(projection)(ist, iprojection + iproj)*pmat%scal(iproj)
      end do
    end do

  end do

  SAFE_DEALLOCATE_A(ind)
  SAFE_DEALLOCATE_A(lpsi)

  POP_SUB(X(hamiltonian_elec_base_nlocal_start))
  call profiling_out(prof_vnlpsi_start)
end subroutine X(hamiltonian_elec_base_nlocal_start)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_nlocal_finish)(this, mesh, bnd, std, projection, vpsib)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_elec_dim_t),          intent(in)    :: std
  type(boundaries_t),               intent(in)    :: bnd
  type(projection_t),       target, intent(inout) :: projection
  class(wfs_elec_t),                intent(inout) :: vpsib

  integer :: ist, ip, imat, nreal, iprojection
  integer :: npoints, nprojs, nst_linear, nphase
  R_TYPE, allocatable :: psi(:, :)
  type(projector_matrix_t), pointer :: pmat
  type(profile_t), save :: reduce_prof
#ifdef R_TCOMPLEX
  integer :: iproj, idim, iphase
  CMPLX  :: phase, phase_pq, phase_mq
  CMPLX, allocatable :: tmp_proj(:, :, :)
#endif

  if(.not. this%has_non_local_potential) return

  ASSERT(this%apply_projector_matrices)

  call profiling_in(prof_vnlpsi_finish, TOSTRING(X(VNLPSI_MAT_KET)))
  PUSH_SUB(X(hamiltonian_elec_base_nlocal_finish))

  nst_linear = vpsib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst_linear
#else
  nreal = nst_linear
#endif
  nphase = 1
  if(bnd%spiral) nphase = 3

  if(vpsib%has_phase) then
    ASSERT(allocated(this%projector_phases))
  end if

  ! reduce the projections
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, TOSTRING(X(VNLPSI_MAT_REDUCE)))
    call comm_allreduce(mesh%vp%comm, projection%X(projection))
    call profiling_out(reduce_prof)
  end if

  if(vpsib%status() == BATCH_DEVICE_PACKED) then

    if(mesh%parallel_in_domains) then
      ! only do this if we have points of some projector matrices
      if(this%max_npoints > 0) then
        call accel_write_buffer(projection%buff_projection, &
          this%full_projection_size*vpsib%pack_size_real(1), projection%X(projection))
      end if
      SAFE_DEALLOCATE_A(projection%X(projection))
    end if

    call finish_opencl()
    call accel_release_buffer(projection%buff_projection)
    if(bnd%spiral) then
      call accel_release_buffer(projection%buff_spin_to_phase)
    end if

    
    POP_SUB(X(hamiltonian_elec_base_nlocal_finish))
    call profiling_out(prof_vnlpsi_finish)
    return
  end if

  ASSERT(allocated(projection%X(projection)))

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(allocated(pmat%zmix)) then
#ifdef R_TCOMPLEX
      SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:vpsib%nst, 1:std%dim))

      do ist = 1, vpsib%nst
        tmp_proj(1:nprojs, ist, 1) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 1), &
                   projection%X(projection)((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs)) &
                                   + matmul(pmat%zmix(1:nprojs, 1:nprojs, 3), &
                   projection%X(projection)((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs))
        tmp_proj(1:nprojs, ist, 2) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 2), &
                   projection%X(projection)((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs)) &
                                   + matmul(pmat%zmix(1:nprojs, 1:nprojs, 4), &
                   projection%X(projection)((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs))
      end do

      do ist = 1, vpsib%nst
        do idim = 1, std%dim
          do iproj = 1, nprojs
            projection%X(projection)((ist-1)*std%dim+idim, iprojection + iproj) = tmp_proj(iproj, ist, idim) 
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(tmp_proj)
#else
      ! Complex projection matrix not allowed for real batches
      ASSERT(.false.)
#endif
   else if(allocated(pmat%dmix)) then

     if(allocated(pmat%dmix)) then
        do ist = 1, nst_linear
          projection%X(projection)(ist, iprojection + 1:iprojection + nprojs) = &
            matmul(pmat%dmix(1:nprojs, 1:nprojs), projection%X(projection)(ist, iprojection + 1:iprojection + nprojs))
        end do
      end if

    end if
    
    if(npoints /=  0) then

      SAFE_ALLOCATE(psi(1:nst_linear, 1:npoints))

      ! Matrix-multiply again.
      ! the line below does: psi = matmul(projection, transpose(pmat%projectors))
      if(.not.pmat%is_cmplx) then
        call blas_gemm('N', 'T', nreal, npoints, nprojs, &
          M_ONE, projection%X(projection)(1, iprojection + 1), nreal, pmat%dprojectors(1, 1), npoints, &
          M_ZERO, psi(1, 1), nreal)
          call profiling_count_operations(nreal*nprojs*M_TWO*npoints)
      else
#ifdef R_TCOMPLEX
        call blas_gemm('N', 'T', nst_linear, npoints, nprojs, &
          M_z1, projection%X(projection)(1, iprojection + 1), nst_linear, pmat%zprojectors(1, 1), npoints, &
          M_z0, psi(1, 1), nst_linear)
        call profiling_count_operations(nst_linear*nprojs*(R_ADD+R_MUL)*npoints)
#else
      ! Complex projection matrix not allowed for real batches
      ASSERT(.false.)
#endif
      end if
      
      call profiling_in(prof_scatter, TOSTRING(X(PROJ_MAT_SCATTER)))

      if(.not. allocated(this%projector_phases)) then    
        ! and copy the points from the local buffer to its position
        if(vpsib%status() == BATCH_PACKED) then
          !$omp parallel do private(ip, ist) if(.not. this%projector_self_overlap)
          do ip = 1, npoints
            do ist = 1, nst_linear\
              vpsib%X(ff_pack)(ist, pmat%map(ip)) = vpsib%X(ff_pack)(ist, pmat%map(ip)) + psi(ist, ip)
            end do
          end do
          !$omp end parallel do
        else
          do ist = 1, nst_linear
            !$omp parallel do if(.not. this%projector_self_overlap)
            do ip = 1, npoints
              vpsib%X(ff_linear)(pmat%map(ip), ist) = vpsib%X(ff_linear)(pmat%map(ip), ist) + psi(ist, ip)
            end do
            !$omp end parallel do
          end do
        end if
        call profiling_count_operations(nst_linear*npoints*R_ADD)
      else
#ifdef R_TCOMPLEX
        if(.not. bnd%spiral) then
          ! and copy the points from the local buffer to its position
          if(vpsib%status() == BATCH_PACKED) then
            !$omp parallel do private(ip, ist, phase) if(.not. this%projector_self_overlap)
            do ip = 1, npoints
              phase = conjg(this%projector_phases(ip, 1, imat, vpsib%ik))
              do ist = 1, nst_linear
                vpsib%zff_pack(ist, pmat%map(ip)) = vpsib%zff_pack(ist, pmat%map(ip)) &
                            + psi(ist, ip)*phase
              end do
            end do
            !$omp end parallel do
          else
            do ist = 1, nst_linear
              !$omp parallel do if(.not. this%projector_self_overlap)
              do ip = 1, npoints
                vpsib%zff_linear(pmat%map(ip), ist) = vpsib%zff_linear(pmat%map(ip), ist) &
                    + psi(ist, ip)*conjg(this%projector_phases(ip, 1, imat, vpsib%ik))
              end do
              !$omp end parallel do
            end do
          end if
          call profiling_count_operations(nst_linear*npoints*(R_ADD+R_MUL))
        else
          ! and copy the points from the local buffer to its position
          if(vpsib%status() == BATCH_PACKED) then
            !$omp parallel do private(ip, ist, phase, phase_pq, phase_mq) if(.not. this%projector_self_overlap)
            do ip = 1, npoints
              phase = conjg(this%projector_phases(ip, 1, imat, vpsib%ik))
              phase_pq = conjg(this%projector_phases(ip, 2, imat, vpsib%ik))
              phase_mq = conjg(this%projector_phases(ip, 3, imat, vpsib%ik))
              do ist = 1, nst_linear, 2
                if(this%spin(3, vpsib%linear_to_ist(ist), vpsib%ik) > 0) then
                  vpsib%zff_pack(ist, pmat%map(ip)) = vpsib%zff_pack(ist, pmat%map(ip)) &
                              + psi(ist, ip)*phase
                  vpsib%zff_pack(ist+1, pmat%map(ip)) = vpsib%zff_pack(ist+1, pmat%map(ip)) &
                              + psi(ist+1, ip)*phase_pq
                else
                  vpsib%zff_pack(ist, pmat%map(ip)) = vpsib%zff_pack(ist, pmat%map(ip)) &
                              + psi(ist, ip)*phase_mq
                  vpsib%zff_pack(ist+1, pmat%map(ip)) = vpsib%zff_pack(ist+1, pmat%map(ip)) &
                              + psi(ist+1, ip)*phase
                end if
              end do
            end do
            !$omp end parallel do
          else
            do ist = 1, nst_linear
              if(this%spin(3, vpsib%linear_to_ist(ist), vpsib%ik) > 0 .and. vpsib%linear_to_idim(ist) == 2) then
                iphase = 2
              else if(this%spin(3, vpsib%linear_to_ist(ist), vpsib%ik) < 0 .and. vpsib%linear_to_idim(ist) == 1) then
                iphase = 3
              else
                iphase = 1
              end if
              !$omp parallel do if(.not. this%projector_self_overlap)
              do ip = 1, npoints
                vpsib%zff_linear(pmat%map(ip), ist) = vpsib%zff_linear(pmat%map(ip), ist) &
                    + psi(ist, ip)*conjg(this%projector_phases(ip, iphase, imat, vpsib%ik))
              end do
              !$omp end parallel do
            end do
          end if
          call profiling_count_operations(nst_linear*npoints*(R_ADD+R_MUL))
        end if
#else
        ! Phases not allowed for real batches
        ASSERT(.false.)
#endif
      end if
      call profiling_out(prof_scatter)
    end if
    
    SAFE_DEALLOCATE_A(psi)
    
    INCR(iprojection, nprojs)
  end do
  
  SAFE_DEALLOCATE_A(projection%X(projection))
  
  POP_SUB(X(hamiltonian_elec_base_nlocal_finish))
  call profiling_out(prof_vnlpsi_finish)

contains

  subroutine finish_opencl()
    integer :: wgsize, imat, iregion, size, padnprojs, lnprojs
    type(profile_t), save :: cl_prof
    type(accel_kernel_t), save, target :: ker_proj_ket, ker_proj_ket_phase, ker_proj_ket_phase_spiral, ker_mix
    type(accel_kernel_t), pointer :: kernel
    type(accel_mem_t), pointer :: buff_proj
    
    PUSH_SUB(X(hamiltonian_elec_base_nlocal_finish).finish_opencl)

    ! In this case we run one kernel per projector, since all write to
    ! the wave-function. Otherwise we would need to do atomic
    ! operations.

    call profiling_in(cl_prof, TOSTRING(X(CL_PROJ_KET)))

    ! only do this if we have points of some projector matrices
    if(this%max_npoints > 0) then

      if(this%projector_mix) then

        SAFE_ALLOCATE(buff_proj)
        call accel_create_buffer(buff_proj, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, this%full_projection_size*vpsib%pack_size_real(1))

        call accel_kernel_start_call(ker_mix, 'projector.cl', 'projector_mix')
        
        call accel_set_kernel_arg(ker_mix, 0, this%nprojector_matrices)
        call accel_set_kernel_arg(ker_mix, 1, this%buff_offsets)
        call accel_set_kernel_arg(ker_mix, 2, this%buff_mix)
        call accel_set_kernel_arg(ker_mix, 3, projection%buff_projection)
        call accel_set_kernel_arg(ker_mix, 4, log2(vpsib%pack_size_real(1)))
        call accel_set_kernel_arg(ker_mix, 5, buff_proj)
        
        padnprojs = pad_pow2(this%max_nprojs)
        lnprojs = min(accel_kernel_workgroup_size(ker_mix)/vpsib%pack_size_real(1), padnprojs)
        
        call accel_kernel_run(ker_mix, &
          (/vpsib%pack_size_real(1), padnprojs, this%nprojector_matrices/), (/vpsib%pack_size_real(1), lnprojs, 1/))
        
        call accel_finish()

      else

        buff_proj => projection%buff_projection
        
      end if
      
      if(allocated(this%projector_phases)) then
        if(bnd%spiral) then
          call accel_kernel_start_call(ker_proj_ket_phase_spiral, 'projector.cl', 'projector_ket_phase_spiral')
          kernel => ker_proj_ket_phase_spiral
        else
          call accel_kernel_start_call(ker_proj_ket_phase, 'projector.cl', 'projector_ket_phase')
          kernel => ker_proj_ket_phase
        endif
        size = vpsib%pack_size(1)
        ASSERT(R_TYPE_VAL == TYPE_CMPLX)
      else
        call accel_kernel_start_call(ker_proj_ket, 'projector.cl', 'projector_ket')
        kernel => ker_proj_ket
        size = vpsib%pack_size_real(1)
      end if

      do iregion = 1, this%nregions
        
        call accel_set_kernel_arg(kernel, 0, this%nprojector_matrices)
        call accel_set_kernel_arg(kernel, 1, this%regions(iregion) - 1)
        call accel_set_kernel_arg(kernel, 2, this%buff_offsets)
        call accel_set_kernel_arg(kernel, 3, this%buff_matrices)
        call accel_set_kernel_arg(kernel, 4, this%buff_maps)
        call accel_set_kernel_arg(kernel, 5, buff_proj)
        call accel_set_kernel_arg(kernel, 6, log2(size))
        call accel_set_kernel_arg(kernel, 7, vpsib%ff_device)
        call accel_set_kernel_arg(kernel, 8, log2(size))

        if(allocated(this%projector_phases)) then
          call accel_set_kernel_arg(kernel, 9, this%buff_projector_phases)
          ! Note: we need to use this%nphase, as the kernel might be called with spiral=false, but 
          !       the phases been built with spiralBC=true
          call accel_set_kernel_arg(kernel, 10, (vpsib%ik - std%kpt%start)*this%total_points*this%nphase)
          if(bnd%spiral) then
            call accel_set_kernel_arg(kernel, 11, projection%buff_spin_to_phase)
          end if  
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
        call profiling_count_operations(nst_linear*npoints*R_ADD)
      end do

      call accel_finish()

      if(this%projector_mix) then
        call accel_release_buffer(buff_proj)
        SAFE_DEALLOCATE_P(buff_proj)
      end if
    end if
    
    call profiling_out(cl_prof)

    POP_SUB(X(hamiltonian_elec_base_nlocal_finish).finish_opencl)
  end subroutine finish_opencl

end subroutine X(hamiltonian_elec_base_nlocal_finish)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_nlocal_force)(this, mesh, st, bnd, iqn, ndim, psi1b, psi2b, force)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  type(mesh_t),                          intent(in)    :: mesh
  type(states_elec_t),                   intent(in)    :: st
  type(boundaries_t),                    intent(in)    :: bnd
  integer,                               intent(in)    :: iqn
  integer,                               intent(in)    :: ndim
  type(wfs_elec_t),                      intent(in)    :: psi1b
  type(wfs_elec_t),                      intent(in)    :: psi2b(:)
  FLOAT,                                 intent(inout) :: force(:, :)

  integer :: ii, ist, ip, iproj, imat, nreal, iprojection, iatom, idir
  integer :: npoints, nprojs, nst
  R_TYPE, allocatable :: psi(:, :, :), projs(:, :, :), ff(:)
  type(projector_matrix_t), pointer :: pmat
#ifdef R_TCOMPLEX
  integer :: idim
  CMPLX, allocatable :: tmp_proj(:, :, :)
#endif

  if(.not. this%has_non_local_potential) return

  ASSERT(this%apply_projector_matrices)
    
  call profiling_in(prof_matelement, TOSTRING(X(VNLPSI_MAT_ELEM)))
  PUSH_SUB(X(hamiltonian_elec_base_nlocal_force))

  ASSERT(psi1b%nst_linear == psi2b(1)%nst_linear)
  ASSERT(psi1b%status() == psi2b(1)%status())
  ASSERT(.not. (psi1b%status() == BATCH_DEVICE_PACKED))
  
  ASSERT(.not.bnd%spiral)

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

      call profiling_in(prof_matelement_gather, TOSTRING(X(PROJ_MAT_ELEM_GATHER)))

      ! collect all the points we need in a continuous array
      if(psi1b%status() == BATCH_PACKED) then
        do ip = 1, npoints
          do ist = 1, nst
            psi(0, ist, ip) = psi1b%X(ff_pack)(ist, pmat%map(ip))
            do idir = 1, ndim
              psi(idir, ist, ip) = psi2b(idir)%X(ff_pack)(ist, pmat%map(ip))
            end do
          end do
        end do
      else
        do ip = 1, npoints
          do ist = 1, nst
            psi(0, ist, ip) = psi1b%X(ff_linear)(pmat%map(ip), ist)
            do idir = 1, ndim
              psi(idir, ist, ip) = psi2b(idir)%X(ff_linear)(pmat%map(ip), ist)
            end do
          end do
        end do
      end if

      if(allocated(this%projector_phases)) then
#ifdef R_TCOMPLEX
        do ip = 1, npoints
          do ist = 1, nst
            do idir = 0, ndim
              psi(idir, ist, ip) = this%projector_phases(ip, 1, imat, psi1b%ik)*psi(idir, ist, ip)
            end do
          end do
        end do
#else
        ! Phases not allowed for real batches
        ASSERT(.false.)
#endif
      end if

      call profiling_out(prof_matelement_gather)

      ! Now matrix-multiply to calculate the projections. We can do all the matrix multiplications at once
      if(.not. pmat%is_cmplx) then
        call blas_gemm('N', 'N', (ndim + 1)*nreal, nprojs, npoints, M_ONE, &
          psi(0, 1, 1), (ndim + 1)*nreal, pmat%dprojectors(1, 1), npoints, &
          M_ZERO, projs(0, 1, iprojection + 1), (ndim + 1)*nreal)

        call profiling_count_operations(nreal*(ndim + 1)*nprojs*M_TWO*npoints)
      else
#ifdef R_TCOMPLEX
        SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:nst*(ndim + 1), 1:1))
        call blas_gemm('C', 'T', nprojs, (ndim + 1)*nst, npoints, &
          M_z1, pmat%zprojectors(1, 1), npoints, psi(0, 1, 1), (ndim + 1)*nst, &
          M_z0, tmp_proj(1,1,1), nprojs)
        do iproj = 1, nprojs
          do ist = 1, nst
            do idir = 0, ndim 
              projs(idir , ist, iprojection + iproj) = tmp_proj(iproj, (ist-1)*(ndim+1)+idir+1, 1)
            end do
          end do
        end do
        SAFE_DEALLOCATE_A(tmp_proj)

        call profiling_count_operations(nst*(ndim + 1)*nprojs*(R_ADD+R_MUL)*npoints)
#endif
      end if

    else
      
      projs(0:ndim, 1:nst, iprojection + 1:iprojection + nprojs) = CNST(0.0)

    end if

    SAFE_DEALLOCATE_A(psi)

    INCR(iprojection, nprojs)

  end do

  if(mesh%parallel_in_domains) then
    call profiling_in(prof_matelement_reduce, TOSTRING(X(VNLPSI_MAT_ELEM_REDUCE)))
    call comm_allreduce(mesh%mpi_grp%comm, projs)
    call profiling_out(prof_matelement_reduce)
  end if
  
  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    
    npoints = pmat%npoints
    nprojs = pmat%nprojs
          
    iatom = this%projector_to_atom(imat)

    if(allocated(pmat%zmix)) then
#ifdef R_TCOMPLEX
      SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:psi1b%nst, 1:st%d%dim))

      do idir = 1, ndim
        do ist = 1, psi1b%nst
          tmp_proj(1:nprojs, ist, 1) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 1), &
                   projs(idir, (ist-1)*st%d%dim+1, iprojection + 1:iprojection + nprojs)) &
                                   + matmul(pmat%zmix(1:nprojs, 1:nprojs, 3), &
                   projs(idir, (ist-1)*st%d%dim+2, iprojection + 1:iprojection + nprojs))
          tmp_proj(1:nprojs, ist, 2) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 2), &
                   projs(idir, (ist-1)*st%d%dim+2, iprojection + 1:iprojection + nprojs)) &
                                   + matmul(pmat%zmix(1:nprojs, 1:nprojs, 4), &
                   projs(idir, (ist-1)*st%d%dim+1, iprojection + 1:iprojection + nprojs))
        end do

        do ist = 1, psi1b%nst
          do idim = 1, st%d%dim
            do iproj = 1, nprojs
              projs(idir, (ist-1)*st%d%dim+idim, iprojection + iproj) = tmp_proj(iproj, ist, idim)  
            end do
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(tmp_proj)
#else
      ! Complex projector matrix not allowed for real batches
      ASSERT(.false.)
#endif
    else if(allocated(pmat%dmix)) then

      do idir = 1, ndim
        do ist = 1, nst
          projs(idir, ist, iprojection + 1:iprojection + nprojs) = &
            matmul(pmat%dmix(1:nprojs, 1:nprojs), projs(idir, ist, iprojection + 1:iprojection + nprojs))
        end do
      end do
    end if
    
    SAFE_ALLOCATE(ff(1:ndim))
    
    ff(1:ndim) = CNST(0.0)
    
    do ii = 1, psi1b%nst_linear
      ist = psi1b%linear_to_ist(ii)
      if(st%d%kweights(iqn)*abs(st%occ(ist, iqn)) <= M_EPSILON) cycle
      do iproj = 1, nprojs
        do idir = 1, ndim
          ff(idir) = ff(idir) - M_TWO*st%d%kweights(iqn)*st%occ(ist, iqn)*pmat%scal(iproj)*mesh%volume_element*&
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

  POP_SUB(X(hamiltonian_elec_base_nlocal_force))
  call profiling_out(prof_matelement)
end subroutine X(hamiltonian_elec_base_nlocal_force)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_elec_base_nlocal_position_commutator)(this, mesh, std, bnd, psib, commpsib)
  type(hamiltonian_elec_base_t), target, intent(in)    :: this
  type(mesh_t),                          intent(in)    :: mesh
  type(states_elec_dim_t),               intent(in)    :: std
  type(boundaries_t),                    intent(in)    :: bnd
  type(wfs_elec_t),                      intent(in)    :: psib
  type(wfs_elec_t),                      intent(inout) :: commpsib(:)

  integer :: ist, ip, iproj, imat, nreal, iprojection, idir
  integer :: npoints, nprojs, nst
  integer, allocatable :: ind(:)
  R_TYPE :: aa, bb, cc, dd
  R_TYPE, allocatable :: projections(:, :, :)
  R_TYPE, allocatable :: psi(:, :, :), lpsi(:,:)
  type(projector_matrix_t), pointer :: pmat
  type(profile_t), save :: prof, reduce_prof
  integer :: wgsize, size
#ifdef R_TCOMPLEX
  integer :: idim
  CMPLX :: phase
  CMPLX, allocatable :: tmp_proj(:, :, :)
#endif

  if(.not. this%has_non_local_potential) return

  ASSERT(this%apply_projector_matrices)

  PUSH_SUB(X(hamiltonian_elec_base_nlocal_position_commutator))
  call profiling_in(prof, TOSTRING(X(COMMUTATOR)))

  ASSERT(psib%is_packed())
  ASSERT(.not.bnd%spiral)

  nst = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

  if(psib%status() == BATCH_DEVICE_PACKED) then
    call X(commutator_opencl)()
    call profiling_out(prof)
    POP_SUB(X(hamiltonian_elec_base_nlocal_position_commutator))
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

#ifdef R_TCOMPLEX
  !$omp parallel do private(imat, pmat, iprojection, npoints, nprojs, iproj, ist, aa, bb, cc, dd, ip, phase, lpsi)
#else
  !$omp parallel do private(imat, pmat, iprojection, npoints, nprojs, iproj, ist, aa, bb, cc, dd, ip, lpsi)
#endif
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)
    iprojection = ind(imat)
    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints == 0) cycle

    SAFE_ALLOCATE(lpsi(1:npoints, 1:nst))
    if(.not. allocated(this%projector_phases)) then
      do ist = 1, nst
        do ip = 1, npoints
          lpsi(ip, ist) = psib%X(ff_pack)(ist, pmat%map(ip))
        end do
      end do
    else
#ifdef R_TCOMPLEX
      do ist = 1, nst
        do ip = 1, npoints
          lpsi(ip, ist) = psib%zff_pack(ist, pmat%map(ip)) &
                            *this%projector_phases(ip, 1, imat, psib%ik)
        end do
      end do
#else
      ! Phases not allowed for real batches
      ASSERT(.false.)
#endif
    end if

    do iproj = 1, nprojs

      if(pmat%is_cmplx) then
#ifdef R_TCOMPLEX
        do ist = 1, nst
          aa = CNST(0.0)
          bb = CNST(0.0)
          cc = CNST(0.0)
          dd = CNST(0.0)
          do ip = 1, npoints
            aa = aa + R_CONJ(pmat%zprojectors(ip, iproj))*lpsi(ip, ist)
            bb = bb + R_CONJ(pmat%zprojectors(ip, iproj))*pmat%position(1, ip)*lpsi(ip, ist)
            cc = cc + R_CONJ(pmat%zprojectors(ip, iproj))*pmat%position(2, ip)*lpsi(ip, ist)
            dd = dd + R_CONJ(pmat%zprojectors(ip, iproj))*pmat%position(3, ip)*lpsi(ip, ist)
          end do
          projections(ist, iprojection + iproj, 0) = pmat%scal(iproj)*aa
          projections(ist, iprojection + iproj, 1) = pmat%scal(iproj)*bb
          projections(ist, iprojection + iproj, 2) = pmat%scal(iproj)*cc
          projections(ist, iprojection + iproj, 3) = pmat%scal(iproj)*dd
        end do
#else
      ! Complex projection matrix not allowed for real batches
      ASSERT(.false.)
#endif
      else
        do ist = 1, nst
          aa = CNST(0.0)
          bb = CNST(0.0)
          cc = CNST(0.0)
          dd = CNST(0.0)
          do ip = 1, npoints
            aa = aa + pmat%dprojectors(ip, iproj)*lpsi(ip, ist)
            bb = bb + pmat%dprojectors(ip, iproj)*pmat%position(1, ip)*lpsi(ip, ist)
            cc = cc + pmat%dprojectors(ip, iproj)*pmat%position(2, ip)*lpsi(ip, ist)
            dd = dd + pmat%dprojectors(ip, iproj)*pmat%position(3, ip)*lpsi(ip, ist)
          end do
          projections(ist, iprojection + iproj, 0) = pmat%scal(iproj)*aa
          projections(ist, iprojection + iproj, 1) = pmat%scal(iproj)*bb
          projections(ist, iprojection + iproj, 2) = pmat%scal(iproj)*cc
          projections(ist, iprojection + iproj, 3) = pmat%scal(iproj)*dd
        end do

      end if
    end do

    SAFE_DEALLOCATE_A(lpsi)
  end do

  ! reduce the projections
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, TOSTRING(X(COMMUTATOR_REDUCE)))
    call comm_allreduce(mesh%mpi_grp%comm, projections)
    call profiling_out(reduce_prof)
  end if

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(allocated(pmat%zmix)) then
#ifdef R_TCOMPLEX
      SAFE_ALLOCATE(tmp_proj(1:nprojs, 1:psib%nst, 1:std%dim))

      do idir = 0, 3
        do ist = 1, psib%nst
          tmp_proj(1:nprojs, ist, 1) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 1), &
                     projections((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs, idir)) &
                                     + matmul(pmat%zmix(1:nprojs, 1:nprojs, 3), &
                     projections((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs, idir))
          tmp_proj(1:nprojs, ist, 2) = matmul(pmat%zmix(1:nprojs, 1:nprojs, 2), &
                     projections((ist-1)*std%dim+2, iprojection + 1:iprojection + nprojs, idir)) &
                                     + matmul(pmat%zmix(1:nprojs, 1:nprojs, 4), &
                     projections((ist-1)*std%dim+1, iprojection + 1:iprojection + nprojs, idir))
        end do

        do ist = 1, psib%nst
          do idim = 1, std%dim
            do iproj = 1, nprojs
              projections((ist-1)*std%dim+idim, iprojection + iproj, idir) = tmp_proj(iproj, ist, idim)
            end do
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(tmp_proj)
#else
      ! Complex projection matrix not allowed for real batches
      ASSERT(.false.)
#endif
   else if(allocated(pmat%dmix)) then
      do idir = 0, 3
        do ist = 1, nst
          projections(ist, iprojection + 1:iprojection + nprojs, idir) = &
            matmul(pmat%dmix(1:nprojs, 1:nprojs), projections(ist, iprojection + 1:iprojection + nprojs, idir))
        end do
      end do
    end if
    
    if(npoints /=  0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints, 0:3))

      ! Matrix-multiply again.
      ! the line below does: psi = matmul(projection, transpose(pmat%projectors))

      if(.not. pmat%is_cmplx) then
        do idir = 0, 3
          call blas_gemm('N', 'T', nreal, npoints, nprojs, &
            M_ONE, projections(1, iprojection + 1, idir), nreal, pmat%dprojectors(1, 1), npoints, &
            M_ZERO, psi(1, 1, idir), nreal)
        end do
        call profiling_count_operations(nreal*nprojs*M_TWO*npoints*4)

      else
 #ifdef R_TCOMPLEX
        do idir = 0, 3
          call blas_gemm('N', 'T', nst, npoints, nprojs, &
            M_z1, projections(1, iprojection + 1, idir), nst, pmat%zprojectors(1, 1), npoints, &
            M_z0, psi(1, 1, idir), nst)
        end do 
#endif
        call profiling_count_operations(nst*nprojs*(R_ADD+R_MUL)*npoints*4)
      end if

      if(allocated(this%projector_phases)) then
#ifdef R_TCOMPLEX
        do idir = 0, 3
          !$omp parallel do private(ip, ist, phase)
          do ip = 1, npoints
            phase = conjg(this%projector_phases(ip, 1, imat, psib%ik))
            do ist = 1, nst
              psi(ist, ip, idir) = phase*psi(ist, ip, idir)
            end do
          end do
          !$omp end parallel do
        end do
        call profiling_count_operations(nst*npoints*3*R_MUL)
#else
        ! Phases not allowed for real batches
        ASSERT(.false.)
#endif
      end if

      do idir = 1, 3
        do ip = 1, npoints
          do ist = 1, nst
            commpsib(idir)%X(ff_pack)(ist, pmat%map(ip)) = commpsib(idir)%X(ff_pack)(ist, pmat%map(ip)) &
              - psi(ist, ip, idir) + pmat%position(idir, ip)*psi(ist, ip, 0)
          end do
        end do
      end do
      
      call profiling_count_operations(nst*npoints*3*(2*R_ADD+R_MUL))
    end if

    SAFE_DEALLOCATE_A(psi)

    INCR(iprojection, nprojs)
  end do

  SAFE_DEALLOCATE_A(ind)

  call profiling_out(prof)
  POP_SUB(X(hamiltonian_elec_base_nlocal_position_commutator))

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

    call accel_create_buffer(buff_proj, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, 4*this%full_projection_size*psib%pack_size_real(1))

    if(allocated(this%projector_phases)) then
      call accel_kernel_start_call(ker_commutator_bra_phase, 'projector.cl', 'projector_commutator_bra_phase')
      kernel => ker_commutator_bra_phase
      size = psib%pack_size(1)
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
    else
      call accel_kernel_start_call(ker_commutator_bra, 'projector.cl', 'projector_commutator_bra')
      size = psib%pack_size_real(1)
      kernel => ker_commutator_bra
    end if
    
    call accel_set_kernel_arg(kernel,  0, this%nprojector_matrices)
    call accel_set_kernel_arg(kernel,  1, this%buff_offsets)
    call accel_set_kernel_arg(kernel,  2, this%buff_matrices)
    call accel_set_kernel_arg(kernel,  3, this%buff_maps)
    call accel_set_kernel_arg(kernel,  4, this%buff_scals)
    call accel_set_kernel_arg(kernel,  5, this%buff_position)
    call accel_set_kernel_arg(kernel,  6, psib%ff_device)
    call accel_set_kernel_arg(kernel,  7, log2(size))
    call accel_set_kernel_arg(kernel,  8, buff_proj)
    call accel_set_kernel_arg(kernel,  9, log2(size))

    if(allocated(this%projector_phases)) then
      call accel_set_kernel_arg(kernel, 10, this%buff_projector_phases)
      call accel_set_kernel_arg(kernel, 11, (psib%ik - std%kpt%start)*this%total_points)
    end if
      
    lnprojs = min(accel_kernel_workgroup_size(kernel)/size, padnprojs)

    call accel_kernel_run(kernel, (/size, padnprojs, this%nprojector_matrices/), (/size, lnprojs, 1/))

    call accel_finish()

    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(proj(1:4*this%full_projection_size*psib%pack_size_real(1)))
      call accel_read_buffer(buff_proj, 4*this%full_projection_size*psib%pack_size_real(1), proj)
      call comm_allreduce(mesh%vp%comm, proj)
      call accel_write_buffer(buff_proj, 4*this%full_projection_size*psib%pack_size_real(1), proj)
      SAFE_DEALLOCATE_A(proj)
    end if

    if(this%projector_mix) then

      SAFE_ALLOCATE(buff_proj_copy)
      
      call accel_create_buffer(buff_proj_copy, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, &
        4*this%full_projection_size*psib%pack_size_real(1))

      size = 4*psib%pack_size_real(1)
      
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
      size = psib%pack_size(1)
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
    else
      call accel_kernel_start_call(ker_commutator_ket, 'projector.cl', 'projector_commutator_ket')
      kernel => ker_commutator_ket
      size = psib%pack_size_real(1)
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
      call accel_set_kernel_arg(kernel,  8, commpsib(1)%ff_device)
      call accel_set_kernel_arg(kernel,  9, commpsib(2)%ff_device)
      call accel_set_kernel_arg(kernel, 10, commpsib(3)%ff_device)
      call accel_set_kernel_arg(kernel, 11, log2(size))

      if(allocated(this%projector_phases)) then
        call accel_set_kernel_arg(kernel, 12, this%buff_projector_phases)
        call accel_set_kernel_arg(kernel, 13, (psib%ik - std%kpt%start)*this%total_points)
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
  
end subroutine X(hamiltonian_elec_base_nlocal_position_commutator)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
