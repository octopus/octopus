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

  integer :: ist, ip
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)
  R_TYPE  :: psi1, psi2
  FLOAT   :: vv
#ifdef HAVE_OPENCL
  integer :: pnp, iprange
#endif

  if(.not. associated(this%potential)) then
    return
  end if

  call profiling_in(prof_vlpsi, "VLPSI")
  PUSH_SUB(X(hamiltonian_base_local))

  if(batch_is_packed(psib) .or. batch_is_packed(vpsib)) then
    ASSERT(batch_is_packed(psib))
    ASSERT(batch_is_packed(vpsib))
    call batch_pack_was_modified(vpsib)
  end if

  select case(batch_status(psib))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    pnp = opencl_padded_size(mesh%np)

    select case(std%ispin)

    case(UNPOLARIZED, SPIN_POLARIZED)
      call opencl_set_kernel_arg(kernel_vpsi, 0, pnp*(ispin - 1))
      call opencl_set_kernel_arg(kernel_vpsi, 1, this%potential_opencl)
      call opencl_set_kernel_arg(kernel_vpsi, 2, psib%pack%buffer)
      call opencl_set_kernel_arg(kernel_vpsi, 3, log2(psib%pack%size_real(1)))
      call opencl_set_kernel_arg(kernel_vpsi, 4, vpsib%pack%buffer)
      call opencl_set_kernel_arg(kernel_vpsi, 5, log2(vpsib%pack%size_real(1)))

      iprange = opencl_max_workgroup_size()/psib%pack%size_real(1)

      call opencl_kernel_run(kernel_vpsi, (/psib%pack%size_real(1), pnp/), (/psib%pack%size_real(1), iprange/))

    case(SPINORS)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 0, this%potential_opencl)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 1, pnp)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 2, psib%pack%buffer)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 3, psib%pack%size(1))
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 4, vpsib%pack%buffer)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 5, vpsib%pack%size(1))

      call opencl_kernel_run(kernel_vpsi_spinors, (/psib%pack%size(1)/2, pnp/), &
        (/psib%pack%size(1)/2, 2*opencl_max_workgroup_size()/psib%pack%size(1)/))

    end select

    call opencl_finish()

    call profiling_count_operations((R_MUL*psib%nst)*mesh%np)
    call profiling_count_transfers(mesh%np, M_ONE)
    call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

#endif
  case(BATCH_PACKED)

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      !$omp parallel do private(vv, ist)
      do ip = 1, mesh%np
        vv = this%potential(ip, ispin)
        forall (ist = 1:psib%nst_linear)
          vpsib%pack%X(psi)(ist, ip) = vpsib%pack%X(psi)(ist, ip) + vv*psib%pack%X(psi)(ist, ip)
        end forall
      end do

      call profiling_count_operations((2*R_ADD*psib%nst_linear)*mesh%np)
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst_linear, R_TOTYPE(M_ONE))

    case(SPINORS)
      ASSERT(mod(psib%nst_linear, 2) == 0)
      !the spinor case is more complicated since it mixes the two components.

      !$omp parallel do private(psi1, psi2, ist)
      do ip = 1, mesh%np
        do ist = 1, psib%nst_linear, 2
          psi1 = psib%pack%zpsi(ist    , ip)
          psi2 = psib%pack%zpsi(ist + 1, ip)
          vpsib%pack%zpsi(ist    , ip) = vpsib%pack%zpsi(ist    , ip) + &
            this%potential(ip, 1)*psi1 + (this%potential(ip, 3) + M_zI*this%potential(ip, 4))*psi2
          vpsib%pack%zpsi(ist + 1, ip) = vpsib%pack%zpsi(ist + 1, ip) + &
            this%potential(ip, 2)*psi2 + (this%potential(ip, 3) - M_zI*this%potential(ip, 4))*psi1
        end do
      end do

      call profiling_count_operations((6*R_ADD + 2*R_MUL)*mesh%np*psib%nst)

    end select

  case(BATCH_NOT_PACKED)

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      !$omp parallel do private(ip)
      do ist = 1, psib%nst
        forall (ip = 1:mesh%np)
          vpsib%states(ist)%X(psi)(ip, 1) = vpsib%states(ist)%X(psi)(ip, 1) + this%potential(ip, ispin)*psib%states(ist)%X(psi)(ip, 1)
        end forall
      end do
      !$omp end parallel do

      call profiling_count_operations((2*R_ADD*psib%nst)*mesh%np)
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

    case(SPINORS)
      !the spinor case is more complicated since it mixes the two components.
      do ist = 1, psib%nst
        psi  => psib%states(ist)%X(psi)
        vpsi => vpsib%states(ist)%X(psi)

        forall(ip = 1:mesh%np)
          vpsi(ip, 1) = vpsi(ip, 1) + this%potential(ip, 1)*psi(ip, 1) + &
            (this%potential(ip, 3) + M_zI*this%potential(ip, 4))*psi(ip, 2)
          vpsi(ip, 2) = vpsi(ip, 2) + this%potential(ip, 2)*psi(ip, 2) + &
            (this%potential(ip, 3) - M_zI*this%potential(ip, 4))*psi(ip, 1)
        end forall

      end do
      call profiling_count_operations((6*R_ADD + 2*R_MUL)*mesh%np*psib%nst)

    end select

  end select

  call profiling_out(prof_vlpsi)
  POP_SUB(X(hamiltonian_base_local))

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
  FLOAT :: cc, b2, bb(1:MAX_DIM)
  CMPLX :: b12

  if(.not. hamiltonian_base_has_magnetic(this)) return

  call profiling_in(prof_magnetic, "MAGNETIC")
  PUSH_SUB(X(hamiltonian_base_magnetic))

  SAFE_ALLOCATE(grad(1:der%mesh%np, 1:MAX_DIM, 1:std%dim))

  do ist = 1, psib%nst
    psi  => psib%states(ist)%X(psi)
    vpsi => vpsib%states(ist)%X(psi)

    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
    grad(:, der%mesh%sb%dim + 1:MAX_DIM, :) = M_ZERO
 
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
  
  POP_SUB(X(hamiltonian_base_magnetic))
  call profiling_out(prof_magnetic)
end subroutine X(hamiltonian_base_magnetic)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_start)(this, mesh, std, ik, psib, projection)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(states_dim_t),          intent(in)    :: std
  integer,                     intent(in)    :: ik
  type(batch_t),               intent(in)    :: psib
  type(projection_t),          intent(out)   :: projection

  integer :: ist, ip, iproj, imat, nreal, iprojection
  integer :: npoints, nprojs, nst
  R_TYPE, allocatable :: psi(:, :)
  type(projector_matrix_t), pointer :: pmat
#ifdef HAVE_OPENCL
  integer :: padnprojs, wgsize
  type(profile_t), save :: cl_prof
#endif
  if(.not. this%apply_projector_matrices) return

  call profiling_in(prof_vnlpsi, "VNLPSI_MAT")
  PUSH_SUB(X(hamiltonian_base_nlocal_start))

  nst = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

#ifdef HAVE_OPENCL
  if(batch_is_packed(psib) .and. opencl_is_enabled()) then
   
    call opencl_create_buffer(projection%buff_projection, CL_MEM_READ_WRITE, R_TYPE_VAL, &
      this%full_projection_size*psib%pack%size_real(1))
    
    call profiling_in(cl_prof, "CL_PROJ_BRA")

    call opencl_set_kernel_arg(kernel_projector_bra, 0, this%nprojector_matrices)
    call opencl_set_kernel_arg(kernel_projector_bra, 1, this%buff_sizes)
    call opencl_set_kernel_arg(kernel_projector_bra, 2, this%buff_offsets)
    call opencl_set_kernel_arg(kernel_projector_bra, 3, this%buff_matrices)
    call opencl_set_kernel_arg(kernel_projector_bra, 4, this%buff_maps)
    call opencl_set_kernel_arg(kernel_projector_bra, 5, this%buff_scals)
    call opencl_set_kernel_arg(kernel_projector_bra, 6, psib%pack%buffer)
    call opencl_set_kernel_arg(kernel_projector_bra, 7, log2(psib%pack%size_real(1)))
    call opencl_set_kernel_arg(kernel_projector_bra, 8, projection%buff_projection)
    call opencl_set_kernel_arg(kernel_projector_bra, 9, log2(psib%pack%size_real(1)))

    padnprojs = pad_pow2(this%max_nprojs)
    wgsize = min(32, opencl_kernel_workgroup_size(kernel_projector_bra)/(psib%pack%size_real(1)*padnprojs))

    call opencl_kernel_run(kernel_projector_bra, &
      (/psib%pack%size_real(1), padnprojs, pad(this%nprojector_matrices, wgsize)/), (/psib%pack%size_real(1), padnprojs, wgsize/))

    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)
      
      npoints = pmat%npoints
      nprojs = pmat%nprojs
      
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints + nst*nprojs)
    end do

    call opencl_finish()

    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(projection%X(projection)(1:psib%pack%size_real(1), 1:this%full_projection_size))
      call opencl_read_buffer(projection%buff_projection, &
        this%full_projection_size*psib%pack%size_real(1), projection%X(projection))
    end if

    call profiling_out(cl_prof)

    POP_SUB(X(hamiltonian_base_nlocal_start))
    call profiling_out(prof_vnlpsi)
    return
  end if
#endif

  SAFE_ALLOCATE(projection%X(projection)(1:nst, 1:this%full_projection_size))
  projection%X(projection) = M_ZERO
  iprojection = 0

  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints /= 0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints))

      call profiling_in(prof_gather, "PROJ_MAT_GATHER")

      ! collect all the points we need in a continous array
      ! MISSING: phases
      if(batch_is_packed(psib)) then
        forall(ip = 1:npoints)
          forall(ist = 1:nst)
            psi(ist, ip) = psib%pack%X(psi)(ist, pmat%map(ip))
          end forall
        end forall
      else
        forall(ist = 1:nst, ip = 1:npoints)
          psi(ist, ip) = psib%states_linear(ist)%X(psi)(pmat%map(ip))
        end forall
      end if

      if(associated(this%projector_phases)) then
        forall(ip = 1:npoints)
          forall(ist = 1:nst)
            psi(ist, ip) = this%projector_phases(ip, imat, ik)*psi(ist, ip)
          end forall
        end forall
      end if

      call profiling_out(prof_gather)
      
      ! Now matrix-multiply to calculate the projections.
      ! the line below does: projection = matmul(psi, pmat%projectors)
      call dgemm('N', 'N', nreal, nprojs, npoints, M_ONE, psi(1, 1), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO,  projection%X(projection)(1, iprojection + 1), nreal)

      ! apply the scale
      forall(ist = 1:nst, iproj = 1:nprojs)
        projection%X(projection)(ist, iprojection + iproj) = projection%X(projection)(ist, iprojection + iproj)*pmat%scal(iproj)
      end forall

      call profiling_count_operations(nreal*nprojs*M_TWO*npoints + nst*nprojs)

    end if

    SAFE_DEALLOCATE_A(psi)

    INCR(iprojection, nprojs)
  end do

  POP_SUB(X(hamiltonian_base_nlocal_start))
  call profiling_out(prof_vnlpsi)
end subroutine X(hamiltonian_base_nlocal_start)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_finish)(this, mesh, std, ik, projection, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(states_dim_t),          intent(in)    :: std
  integer,                     intent(in)    :: ik
  type(projection_t),          intent(inout) :: projection
  type(batch_t),               intent(inout) :: vpsib

  integer :: ist, ip, iproj, imat, nreal, iprojection
  integer :: npoints, nprojs, nst, d1
  R_TYPE, allocatable :: psi(:, :)
  type(projector_matrix_t), pointer :: pmat
#ifdef HAVE_MPI
  R_TYPE, allocatable :: projection_red(:, :)
  type(profile_t), save :: reduce_prof
#endif

  if(.not. this%apply_projector_matrices) return

  call profiling_in(prof_vnlpsi, "VNLPSI_MAT")
  PUSH_SUB(X(hamiltonian_base_nlocal_finish))

  nst = vpsib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

  ! reduce the projections
#ifdef HAVE_MPI
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, "VNLPSI_MAT_REDUCE")
    d1 = ubound(projection%X(projection), dim = 1)
    SAFE_ALLOCATE(projection_red(1:d1, 1:this%full_projection_size))
    projection_red = projection%X(projection)
    call MPI_Allreduce(projection_red(1, 1), projection%X(projection)(1, 1), d1*this%full_projection_size, &
      R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    SAFE_DEALLOCATE_A(projection_red)
    call profiling_out(reduce_prof)
  end if
#endif

#ifdef HAVE_OPENCL
  if(batch_is_packed(vpsib) .and. opencl_is_enabled()) then

    if(mesh%parallel_in_domains) then
      call opencl_write_buffer(projection%buff_projection, &
        this%full_projection_size*vpsib%pack%size_real(1), projection%X(projection))
      SAFE_DEALLOCATE_P(projection%X(projection))
    end if

    call finish_opencl()
    call opencl_release_buffer(projection%buff_projection)

    POP_SUB(X(hamiltonian_base_nlocal_finish))
    call profiling_out(prof_vnlpsi)
    return
  end if
#endif

  ASSERT(associated(projection%X(projection)))

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints /=  0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints))

      ! Matrix-multiply again.
      ! the line below does: psi = matmul(projection, transpose(pmat%projectors))
      call dgemm('N', 'T', nreal, npoints, nprojs, &
        M_ONE, projection%X(projection)(1, iprojection + 1), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO,  psi(1, 1), nreal)
      
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints)

      call profiling_in(prof_scatter, "PROJ_MAT_SCATTER")

      if(associated(this%projector_phases)) then
        forall(ip = 1:npoints)
          forall(ist = 1:nst)
            psi(ist, ip) = conjg(this%projector_phases(ip, imat, ik))*psi(ist, ip)
          end forall
        end forall
      end if
      
      ! and copy the points from the local buffer to its position
      ! MISSING: phases
      if(batch_is_packed(vpsib)) then
        forall(ip = 1:npoints)
          forall(ist = 1:nst)
            vpsib%pack%X(psi)(ist, pmat%map(ip)) = vpsib%pack%X(psi)(ist, pmat%map(ip)) + psi(ist, ip)
          end forall
        end forall

        call batch_pack_was_modified(vpsib)
      else
        do ist = 1, nst
          forall(ip = 1:npoints)
            vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) = vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) + psi(ist, ip)
          end forall
        end do
      end if
      call profiling_count_operations(nst*npoints*R_ADD)
      call profiling_out(prof_scatter)
    end if
    
    SAFE_DEALLOCATE_A(psi)
    
    INCR(iprojection, nprojs)
  end do
  
  SAFE_DEALLOCATE_P(projection%X(projection))
  
  POP_SUB(X(hamiltonian_base_nlocal_finish))
  call profiling_out(prof_vnlpsi)

contains

  subroutine finish_opencl()
#ifdef HAVE_OPENCL
    integer :: wgsize, imat
    type(profile_t), save :: cl_prof
    type(opencl_mem_t) :: lpsi

    ! In this case we run one kernel per projector, since all write to
    ! the wave-function. Otherwise we would need to do atomic
    ! operations.

    call profiling_in(cl_prof, "CL_PROJ_KET")

    call opencl_create_buffer(lpsi, CL_MEM_READ_WRITE, TYPE_FLOAT, this%total_points*vpsib%pack%size_real(1))

    call opencl_set_kernel_arg(kernel_projector_ket, 0, this%nprojector_matrices)
    call opencl_set_kernel_arg(kernel_projector_ket, 1, this%buff_sizes)
    call opencl_set_kernel_arg(kernel_projector_ket, 2, this%buff_offsets)
    call opencl_set_kernel_arg(kernel_projector_ket, 3, this%buff_matrices)
    call opencl_set_kernel_arg(kernel_projector_ket, 4, this%buff_maps)
    call opencl_set_kernel_arg(kernel_projector_ket, 5, projection%buff_projection)
    call opencl_set_kernel_arg(kernel_projector_ket, 6, log2(vpsib%pack%size_real(1)))
    call opencl_set_kernel_arg(kernel_projector_ket, 7, lpsi)
    call opencl_set_kernel_arg(kernel_projector_ket, 8, log2(vpsib%pack%size_real(1)))
    
    wgsize = opencl_max_workgroup_size()/vpsib%pack%size_real(1)    

    call opencl_kernel_run(kernel_projector_ket, &
      (/vpsib%pack%size_real(1), pad(this%max_npoints, wgsize), this%nprojector_matrices/), &
      (/vpsib%pack%size_real(1), wgsize, 1/))

    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)
      npoints = pmat%npoints
      nprojs = pmat%nprojs
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints)
      call profiling_count_operations(nst*npoints*R_ADD)
    end do

    call opencl_finish()

    call opencl_set_kernel_arg(kernel_projector_ket_copy, 0, mesh%np)
    call opencl_set_kernel_arg(kernel_projector_ket_copy, 1, this%buff_pos)
    call opencl_set_kernel_arg(kernel_projector_ket_copy, 2, this%buff_invmap)
    call opencl_set_kernel_arg(kernel_projector_ket_copy, 3, lpsi)
    call opencl_set_kernel_arg(kernel_projector_ket_copy, 4, log2(vpsib%pack%size_real(1)))
    call opencl_set_kernel_arg(kernel_projector_ket_copy, 5, vpsib%pack%buffer)
    call opencl_set_kernel_arg(kernel_projector_ket_copy, 6, log2(vpsib%pack%size_real(1)))
    
    wgsize = opencl_max_workgroup_size()/vpsib%pack%size_real(1)
    
    call opencl_kernel_run(kernel_projector_ket_copy, &
      (/vpsib%pack%size_real(1), pad(mesh%np, wgsize)/), (/vpsib%pack%size_real(1), wgsize/))
    
    call batch_pack_was_modified(vpsib)
    call opencl_finish()

    call opencl_release_buffer(lpsi)
    call profiling_out(cl_prof)

#endif
  end subroutine finish_opencl

end subroutine X(hamiltonian_base_nlocal_finish)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
