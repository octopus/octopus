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

  integer :: ist, idim, ip, iprange
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)
#ifdef HAVE_OPENCL
  integer :: pnp
#endif

  call profiling_in(prof_vlpsi, "VLPSI")
  call push_sub('hamiltonian_base_inc.Xhamiltonian_base_local')

  if(associated(this%potential)) then
#ifdef HAVE_OPENCL
    if(opencl_is_enabled() .and. batch_is_in_buffer(psib)) then
      ASSERT(batch_is_in_buffer(vpsib))
      
      pnp = opencl_padded_size(mesh%np)

      select case(std%ispin)

      case(UNPOLARIZED, SPIN_POLARIZED)
        call opencl_set_kernel_arg(kernel_vpsi, 0, pnp*(ispin - 1))
        call opencl_set_kernel_arg(kernel_vpsi, 1, this%potential_opencl)
        call opencl_set_kernel_arg(kernel_vpsi, 2, psib%buffer)
        call opencl_set_kernel_arg(kernel_vpsi, 3, log2(psib%ubound_real(1)))
        call opencl_set_kernel_arg(kernel_vpsi, 4, vpsib%buffer)
        call opencl_set_kernel_arg(kernel_vpsi, 5, log2(vpsib%ubound_real(1)))

        iprange = opencl_max_workgroup_size()/psib%ubound_real(1)

        call opencl_kernel_run(kernel_vpsi, (/psib%ubound_real(1), pnp/), (/psib%ubound_real(1), iprange/))

      case(SPINORS)
        call opencl_set_kernel_arg(kernel_vpsi_spinors, 0, this%potential_opencl)
        call opencl_set_kernel_arg(kernel_vpsi_spinors, 1, pnp)
        call opencl_set_kernel_arg(kernel_vpsi_spinors, 2, psib%buffer)
        call opencl_set_kernel_arg(kernel_vpsi_spinors, 3, batch_buffer_ubound(psib))
        call opencl_set_kernel_arg(kernel_vpsi_spinors, 4, vpsib%buffer)
        call opencl_set_kernel_arg(kernel_vpsi_spinors, 5, batch_buffer_ubound(vpsib))

        call opencl_kernel_run(kernel_vpsi_spinors, (/psib%ubound(1)/2, pnp/), &
          (/psib%ubound(1)/2, 2*opencl_max_workgroup_size()/psib%ubound(1)/))

      end select

      call batch_buffer_was_modified(vpsib)

      call profiling_count_operations((R_MUL*psib%nst)*mesh%np)
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))
      call profiling_out(prof_vlpsi)
      call pop_sub('hamiltonian_base_inc.Xhamiltonian_base_local')
      return
    end if
#endif
    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      !$omp parallel do private(idim, ip)
      do ist = 1, psib%nst
        forall (idim = 1:psib%dim, ip = 1:mesh%np)
          vpsib%states(ist)%X(psi)(ip, idim) = vpsib%states(ist)%X(psi)(ip, idim) + &
            this%potential(ip, ispin)*psib%states(ist)%X(psi)(ip, idim)
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
  FLOAT :: cc, b2, bb(1:MAX_DIM)
  CMPLX :: b12

  if(.not. hamiltonian_base_has_magnetic(this)) return

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

subroutine X(hamiltonian_base_non_local)(this, mesh, std, ik, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(states_dim_t),          intent(in)    :: std
  integer,                     intent(in)    :: ik
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  integer :: ist, ip, iproj, imat
  integer :: npoints, nprojs, nst
  R_TYPE, allocatable :: psi(:, :), projection(:,  :)
  type(projector_matrix_t), pointer :: pmat
#ifdef HAVE_MPI
  R_TYPE, allocatable :: projection_red(:, :)
#endif
#ifdef HAVE_OPENCL
  type(opencl_mem_t) :: buff_map
#endif

  if(.not. this%apply_projector_matrices) return

  call profiling_in(prof_vnlpsi, "VNLPSI_MAT")
  call push_sub('hamiltonian_base_inc.Xhamiltonian_base_non_local')

  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs
    nst = psib%nst_linear

#ifdef HAVE_OPENCL
    if(batch_is_in_buffer(psib) .or. batch_is_in_buffer(vpsib)) then
      call opencl_create_buffer(buff_map, CL_MEM_READ_ONLY, TYPE_INTEGER, npoints)
      call opencl_write_buffer(buff_map, npoints, pmat%map)
    end if
#endif

    SAFE_ALLOCATE(projection(1:nprojs, 1:nst))

    if(npoints == 0) then
      ! With domain parallelization it might happen that this
      ! processor does not have any point. In this case we only need
      ! to call MPI_Allreduce with zero.

      projection = M_ZERO

    else

      SAFE_ALLOCATE(psi(1:npoints, 1:nst))

      call profiling_in(prof_gather, "PROJ_MAT_GATHER")
      if(batch_is_in_buffer(psib)) then
        call gather_opencl()
      else
        ! collect all the points we need in a continous array
        forall(ist = 1:nst, ip = 1:npoints)
          psi(ip, ist) = psib%states_linear(ist)%X(psi)(pmat%map(ip))
        !MISSING: phases
        end forall
      endif
      call profiling_out(prof_gather)

      ! Now matrix-multiply to calculate the projections. Since the
      ! projector matrix is always real we can only use blas in the
      ! real case.
#ifdef R_TREAL
      call blas_gemm('T', 'N', nprojs, nst, npoints, &
        M_ONE, pmat%projectors(1, 1), npoints, psi(1, 1), npoints, &
        M_ZERO, projection(1, 1), nprojs)
#else
      projection(1:nprojs, 1:nst) = matmul(transpose(pmat%projectors(1:npoints, 1:nprojs)), psi(1:npoints, 1:nst))
#endif

      ! apply the scale
      forall(ist = 1:nst, iproj = 1:nprojs) projection(iproj, ist) = projection(iproj, ist)*pmat%scal(iproj)

    end if

    ! now reduce the projections
#ifdef HAVE_MPI
    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(projection_red(1:nprojs, 1:nst))
      forall(ist = 1:nst, iproj = 1:nprojs) projection_red(iproj, ist) = projection(iproj, ist)
      call MPI_Allreduce(projection_red(1, 1), projection(1, 1), nst*nprojs, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
      SAFE_DEALLOCATE_A(projection_red)
    end if
#endif

    if(npoints ==  0) then
      SAFE_DEALLOCATE_A(projection)
      cycle
    end if

    ! Matrix-multiply again.
#ifdef R_TREAL
    call blas_gemm('N', 'N', npoints, nst, nprojs, &
      M_ONE, pmat%projectors(1, 1), npoints, projection(1, 1), nprojs, &
      M_ZERO, psi(1, 1), npoints)
#else
    psi(1:npoints, 1:nst) = matmul(pmat%projectors(1:npoints, 1:nprojs), projection(1:nprojs, 1:nst))
#endif

    call profiling_count_operations(M_TWO*nst*nprojs*M_TWO*npoints*R_ADD + nst*nprojs*R_ADD)

    call profiling_in(prof_scatter, "PROJ_MAT_SCATTER")

    if(batch_is_in_buffer(vpsib)) then
      call scatter_opencl()
    else
      ! and copy the points from the local buffer to its position
      do ist = 1, nst
        forall(ip = 1:npoints)
          vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) = vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) + psi(ip, ist)
          !MISSING: phases
        end forall
      end do
    end if

    call profiling_count_operations(nst*npoints*R_ADD)
    call profiling_out(prof_scatter)

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(projection)

#ifdef HAVE_OPENCL
    if(batch_is_in_buffer(psib) .or. batch_is_in_buffer(vpsib)) then
      call opencl_release_buffer(buff_map)
    end if
#endif
  end do

  call pop_sub('hamiltonian_base_inc.Xhamiltonian_base_non_local')
  call profiling_out(prof_vnlpsi)
contains

  subroutine gather_opencl()
#ifdef HAVE_OPENCL
    type(opencl_mem_t) :: buff_psi
    integer :: wgsize

    call opencl_create_buffer(buff_psi, CL_MEM_WRITE_ONLY, R_TYPE_VAL, npoints*psib%nst_linear)
    
    call opencl_set_kernel_arg(X(projector_gather), 0, npoints)
    call opencl_set_kernel_arg(X(projector_gather), 1, buff_map)
    call opencl_set_kernel_arg(X(projector_gather), 2, psib%buffer)
    call opencl_set_kernel_arg(X(projector_gather), 3, batch_buffer_ubound(psib))
    call opencl_set_kernel_arg(X(projector_gather), 4, buff_psi)
    call opencl_set_kernel_arg(X(projector_gather), 5, npoints)

    wgsize = opencl_max_workgroup_size()

    call opencl_kernel_run(X(projector_gather), (/pad(npoints, wgsize), psib%nst_linear/), (/wgsize, 1/))
    
    call opencl_read_buffer(buff_psi, npoints*psib%nst_linear, psi)
    call opencl_release_buffer(buff_psi)
#endif
  end subroutine gather_opencl

  subroutine scatter_opencl()
#ifdef HAVE_OPENCL
    type(opencl_mem_t) :: buff_psi
    integer :: wgsize

    call opencl_create_buffer(buff_psi, CL_MEM_READ_ONLY, R_TYPE_VAL, npoints*vpsib%nst_linear)
    call opencl_write_buffer(buff_psi, npoints*vpsib%nst_linear, psi)
    
    call opencl_set_kernel_arg(X(projector_scatter), 0, npoints)
    call opencl_set_kernel_arg(X(projector_scatter), 1, buff_map)
    call opencl_set_kernel_arg(X(projector_scatter), 2, buff_psi)
    call opencl_set_kernel_arg(X(projector_scatter), 3, npoints)
    call opencl_set_kernel_arg(X(projector_scatter), 4, vpsib%buffer)
    call opencl_set_kernel_arg(X(projector_scatter), 5, batch_buffer_ubound(vpsib))

    wgsize = opencl_max_workgroup_size()

    call opencl_kernel_run(X(projector_scatter), (/pad(npoints, wgsize), vpsib%nst_linear/), (/wgsize, 1/))
    
    call batch_buffer_was_modified(vpsib)

    call opencl_release_buffer(buff_psi)
#endif
  end subroutine scatter_opencl

end subroutine X(hamiltonian_base_non_local)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
