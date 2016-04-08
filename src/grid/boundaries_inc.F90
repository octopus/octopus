!! Copyright (C) 2005-2006 Florian Lorenzen, Heiko Appel
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
!! $Id$

! ---------------------------------------------------------
!> Updates ghost points of every node. A vector suitable
!! for non-local operations contains local values and
!! ghost point values.
!! Length of v_local must be
!! vp%np_local+vp%np_ghost
subroutine X(vec_ghost_update)(vp, v_local)
  type(pv_t), intent(in)    :: vp
  R_TYPE,     intent(inout) :: v_local(:)

  R_TYPE,  allocatable :: ghost_send(:)
  integer              :: nsend
  
  call profiling_in(prof_update, "GHOST_UPDATE")

  PUSH_SUB(X(vec_ghost_update))

  nsend = subarray_size(vp%ghost_spoints)
  SAFE_ALLOCATE(ghost_send(1:nsend))
  call X(subarray_gather)(vp%ghost_spoints, v_local, ghost_send)

#ifdef HAVE_MPI
  call mpi_debug_in(vp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(ghost_send(1), vp%ghost_scounts(1), vp%ghost_sdispls(1), R_MPITYPE, &
       v_local(vp%np_local+1), vp%ghost_rcounts(1), vp%ghost_rdispls(1), R_MPITYPE, &
       vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_ALLTOALLV)
#endif

  SAFE_DEALLOCATE_A(ghost_send)

  POP_SUB(X(vec_ghost_update))

  call profiling_out(prof_update)
end subroutine X(vec_ghost_update)

! ---------------------------------------------------------

subroutine X(ghost_update_batch_start)(vp, v_local, handle)
  type(pv_t),    target,    intent(in)    :: vp
  type(batch_t), target,    intent(inout) :: v_local
  type(pv_handle_batch_t),  intent(out)   :: handle

  integer :: ipart, pos, ii, tag, nn

  call profiling_in(prof_start, "GHOST_UPDATE_START")
  PUSH_SUB(X(ghost_update_batch_start))

  ASSERT(v_local%nst_linear > 0)

  handle%nnb = 0
  handle%v_local => v_local
  handle%vp => vp

  SAFE_ALLOCATE(handle%requests(1:2*vp%npart*v_local%nst_linear))

  ! first post the receptions
  select case(batch_status(v_local))

  case(BATCH_CL_PACKED)
    SAFE_ALLOCATE(handle%X(recv_buffer)(1:v_local%pack%size(1)*vp%np_ghost))

    do ipart = 1, vp%npart
      if(vp%ghost_rcounts(ipart) == 0) cycle
      
      handle%nnb = handle%nnb + 1
      tag = 0
      pos = 1 + vp%ghost_rdispls(ipart)*v_local%pack%size(1)
#ifdef HAVE_MPI
      call MPI_Irecv(handle%X(recv_buffer)(pos), vp%ghost_rcounts(ipart)*v_local%pack%size(1), R_MPITYPE, &
           ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
#endif
    end do

  case(BATCH_PACKED)
    !In this case, data from different vectors is contiguous. So we can use one message per partition.
    do ipart = 1, vp%npart
      if(vp%ghost_rcounts(ipart) == 0) cycle
      
      handle%nnb = handle%nnb + 1
      tag = 0
      pos = vp%np_local + 1 + vp%ghost_rdispls(ipart)
#ifdef HAVE_MPI
      call MPI_Irecv(v_local%pack%X(psi)(1, pos), vp%ghost_rcounts(ipart)*v_local%pack%size(1), R_MPITYPE, &
           ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
#endif
    end do

  case(BATCH_NOT_PACKED)
    do ii = 1, v_local%nst_linear
      do ipart = 1, vp%npart
        if(vp%ghost_rcounts(ipart) == 0) cycle
        
        handle%nnb = handle%nnb + 1
        tag = ii
        pos = vp%np_local + 1 + vp%ghost_rdispls(ipart)
#ifdef HAVE_MPI
        call MPI_Irecv(v_local%states_linear(ii)%X(psi)(pos), vp%ghost_rcounts(ipart), R_MPITYPE, ipart - 1, tag, &
        vp%comm, handle%requests(handle%nnb), mpi_err)
#endif
      end do
    end do

  end select


  call batch_init(handle%ghost_send, 1, v_local%nst_linear)
  call X(batch_allocate)(handle%ghost_send, 1, v_local%nst_linear, subarray_size(vp%ghost_spoints))

  if(batch_is_packed(v_local)) call batch_pack(handle%ghost_send, copy = .false.)

  !now collect the data for sending
  call X(subarray_gather_batch)(vp%ghost_spoints, v_local, handle%ghost_send)

  if(batch_status(v_local) == BATCH_CL_PACKED) then
    nn = product(handle%ghost_send%pack%size(1:2))
    SAFE_ALLOCATE(handle%X(send_buffer)(1:nn))
#ifdef HAVE_OPENCL
    call opencl_read_buffer(handle%ghost_send%pack%buffer, nn, handle%X(send_buffer))
#endif
  end if

  select case(batch_status(v_local))

  case(BATCH_CL_PACKED)
    do ipart = 1, vp%npart
      if(vp%ghost_scounts(ipart) == 0) cycle
      handle%nnb = handle%nnb + 1
      tag = 0
#ifdef HAVE_MPI
      call MPI_Isend(handle%X(send_buffer)(1 + (vp%ghost_sendpos(ipart) - 1)*v_local%pack%size(1)), &
        vp%ghost_scounts(ipart)*v_local%pack%size(1), &
        R_MPITYPE, ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
#endif
    end do

  case(BATCH_PACKED)
    do ipart = 1, vp%npart
      if(vp%ghost_scounts(ipart) == 0) cycle
      handle%nnb = handle%nnb + 1
      tag = 0
#ifdef HAVE_MPI
      call MPI_Isend(handle%ghost_send%pack%X(psi)(1, vp%ghost_sendpos(ipart)), &
        vp%ghost_scounts(ipart)*v_local%pack%size(1), &
        R_MPITYPE, ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
#endif
    end do

  case(BATCH_NOT_PACKED)
    do ii = 1, v_local%nst_linear
      do ipart = 1, vp%npart
        if(vp%ghost_scounts(ipart) == 0) cycle
        handle%nnb = handle%nnb + 1
        tag = ii
#ifdef HAVE_MPI
        call MPI_Isend(handle%ghost_send%states_linear(ii)%X(psi)(vp%ghost_sendpos(ipart)), &
             vp%ghost_scounts(ipart), R_MPITYPE, ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
#endif
      end do
    end do
  end select

  POP_SUB(X(ghost_update_batch_start))
  call profiling_out(prof_start)

end subroutine X(ghost_update_batch_start)

! ---------------------------------------------------------

subroutine X(ghost_update_batch_finish)(handle)
  type(pv_handle_batch_t),  intent(inout)   :: handle

  integer, allocatable :: status(:, :)

  call profiling_in(prof_wait, "GHOST_UPDATE_WAIT")
  PUSH_SUB(X(ghost_update_batch_finish))
  
  ASSERT(handle%nnb > 0)

#ifdef HAVE_MPI
  SAFE_ALLOCATE(status(1:MPI_STATUS_SIZE, 1:handle%nnb))
  call MPI_Waitall(handle%nnb, handle%requests(1), status(1, 1), mpi_err)
#endif

  SAFE_DEALLOCATE_A(status)
  SAFE_DEALLOCATE_P(handle%requests)

  if(batch_status(handle%v_local) == BATCH_CL_PACKED) then
#ifdef HAVE_OPENCL
    call opencl_write_buffer(handle%v_local%pack%buffer, handle%v_local%pack%size(1)*handle%vp%np_ghost, &
      handle%X(recv_buffer), offset = handle%v_local%pack%size(1)*handle%vp%np_local)
#endif
    SAFE_DEALLOCATE_P(handle%X(send_buffer))
    SAFE_DEALLOCATE_P(handle%X(recv_buffer))
  end if

  call batch_end(handle%ghost_send, copy = .false.)

  call profiling_out(prof_wait)
  POP_SUB(X(ghost_update_batch_finish))
end subroutine X(ghost_update_batch_finish)

! ---------------------------------------------------------
!> Set all boundary points in ffb to zero to implement zero
!! boundary conditions for the derivatives, in finite system;
!! or set according to periodic boundary conditions.
subroutine X(boundaries_set_batch)(boundaries, ffb)
  type(boundaries_t),    intent(in)    :: boundaries
  type(batch_t), target, intent(inout) :: ffb

  integer :: bndry_start, bndry_end

  PUSH_SUB(X(boundaries_set_batch))
  call profiling_in(set_bc_prof, 'SET_BC')
  
  ASSERT(batch_type(ffb) == R_TYPE_VAL)

  ! The boundary points are at different locations depending on the presence
  ! of ghost points due to domain parallelization.
  if(boundaries%mesh%parallel_in_domains) then
    bndry_start = boundaries%mesh%vp%np_local + boundaries%mesh%vp%np_ghost + 1
    bndry_end   = boundaries%mesh%vp%np_local + boundaries%mesh%vp%np_ghost + boundaries%mesh%vp%np_bndry
  else
    bndry_start = boundaries%mesh%np + 1
    bndry_end   = boundaries%mesh%np_part
  end if
    
  if(boundaries%mesh%sb%periodic_dim < boundaries%mesh%sb%dim) call zero_boundaries()
  if(boundaries%mesh%sb%mr_flag) call multiresolution()
  if(boundaries%mesh%sb%periodic_dim > 0)     call periodic()

  call profiling_out(set_bc_prof)
  POP_SUB(X(boundaries_set_batch))

contains

  ! ---------------------------------------------------------
  subroutine zero_boundaries()
    integer :: ist, ip
#ifdef HAVE_OPENCL
    integer :: np
#endif

    PUSH_SUB(X(boundaries_set_batch).zero_boundaries)

    select case(batch_status(ffb))
    case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
      np = ffb%pack%size(1)*(bndry_end - bndry_start + 1)
      call opencl_set_buffer_to_zero(ffb%pack%buffer, batch_type(ffb), np, offset = ffb%pack%size(1)*(bndry_start - 1))
      call opencl_finish()
#endif
    case(BATCH_PACKED)
      forall(ip = bndry_start:bndry_end) 
        forall(ist = 1:ffb%nst_linear)
          ffb%pack%X(psi)(ist, ip) = R_TOTYPE(M_ZERO)
        end forall
      end forall

    case(BATCH_NOT_PACKED)
      do ist = 1, ffb%nst_linear
        forall (ip = bndry_start:bndry_end) ffb%states_linear(ist)%X(psi)(ip) = R_TOTYPE(M_ZERO)
      end do

    end select

    call batch_pack_was_modified(ffb)

    POP_SUB(X(boundaries_set_batch).zero_boundaries)
  end subroutine zero_boundaries


  ! ---------------------------------------------------------
  subroutine multiresolution()
    integer :: ist, ip
    integer :: ii, jj, kk, ix, iy, iz, dx, dy, dz, i_lev
    FLOAT :: weight
    R_TYPE, pointer :: ff(:)

    PUSH_SUB(X(boundaries_set_batch).multiresolution)

    do ist = 1, ffb%nst_linear
      ff => ffb%states_linear(ist)%X(psi)

      do ip = bndry_start, bndry_end
        ix = boundaries%mesh%idx%lxyz(ip, 1)
        iy = boundaries%mesh%idx%lxyz(ip, 2)
        iz = boundaries%mesh%idx%lxyz(ip, 3)

        i_lev = boundaries%mesh%resolution(ix,iy,iz)

        ! resolution is 2**num_radii for outer boundary points, but now we want inner boundary points
        if(i_lev /= 2**boundaries%mesh%sb%hr_area%num_radii) then
          dx = abs(mod(ix, 2**(i_lev)))
          dy = abs(mod(iy, 2**(i_lev)))
          dz = abs(mod(iz, 2**(i_lev)))

          do ii = 1, boundaries%mesh%sb%hr_area%interp%nn
            do jj = 1, boundaries%mesh%sb%hr_area%interp%nn
              do kk = 1, boundaries%mesh%sb%hr_area%interp%nn
                weight = boundaries%mesh%sb%hr_area%interp%ww(ii) * &
                  boundaries%mesh%sb%hr_area%interp%ww(jj) *        &
                  boundaries%mesh%sb%hr_area%interp%ww(kk)

                ff(ip) = ff(ip) + weight * ff(boundaries%mesh%idx%lxyz_inv(   &
                  ix + boundaries%mesh%sb%hr_area%interp%posi(ii) * dx,       &
                  iy + boundaries%mesh%sb%hr_area%interp%posi(jj) * dy,       &
                  iz + boundaries%mesh%sb%hr_area%interp%posi(kk) * dz))
              end do
            end do
          end do
        end if

      end do ! ip
    end do ! ist

    POP_SUB(X(boundaries_set_batch).multiresolution)
  end subroutine multiresolution


  ! ---------------------------------------------------------
  subroutine periodic()
    integer :: ip, ist, ip_bnd, ip_inn
    R_TYPE, pointer :: ff(:)

#ifdef HAVE_MPI
    R_TYPE, allocatable :: sendbuffer(:, :, :)
    R_TYPE, allocatable :: recvbuffer(:, :, :)
    integer, allocatable :: send_disp(:), send_count(:)
    integer, allocatable :: recv_disp(:), recv_count(:)
    integer :: ipart, npart, maxsend, maxrecv, ldbuffer, ip2
#endif
#ifdef HAVE_OPENCL
    type(octcl_kernel_t), save :: kernel_send, kernel_recv, kernel
    type(cl_kernel) :: kernel_ref
    integer :: wgsize
    type(opencl_mem_t) :: buff_send
    type(opencl_mem_t) :: buff_recv
#endif

    PUSH_SUB(X(boundaries_set_batch).periodic)

#ifdef HAVE_MPI
    if(boundaries%mesh%parallel_in_domains) then

      call profiling_in(set_bc_precomm_prof, 'SET_BC_PRECOMM')

      npart = boundaries%mesh%vp%npart
      maxsend = maxval(boundaries%nsend(1:npart))
      maxrecv = maxval(boundaries%nrecv(1:npart))

      ldbuffer = ffb%nst_linear
      if(batch_status(ffb) == BATCH_CL_PACKED) ldbuffer = ffb%pack%size(1)
      SAFE_ALLOCATE(sendbuffer(1:ldbuffer, 1:maxsend, 1:npart))

      select case(batch_status(ffb))

      case(BATCH_NOT_PACKED)

        do ipart = 1, npart
          do ip = 1, boundaries%nsend(ipart)
            ip2 = boundaries%per_send(ip, ipart)
            do ist = 1, ffb%nst_linear
              sendbuffer(ist, ip, ipart) = ffb%states_linear(ist)%X(psi)(ip2)
            end do
          end do
        end do

      case(BATCH_PACKED)

        do ipart = 1, npart
          !$omp parallel do private(ip, ip2, ist)
          do ip = 1, boundaries%nsend(ipart)
            ip2 = boundaries%per_send(ip, ipart)
            do ist = 1, ffb%nst_linear
              sendbuffer(ist, ip, ipart) = ffb%pack%X(psi)(ist, ip2)
            end do
          end do
        end do

      case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
        call opencl_create_buffer(buff_send, CL_MEM_WRITE_ONLY, R_TYPE_VAL, ffb%pack%size(1)*maxsend*npart)

        call octcl_kernel_start_call(kernel_send, 'boundaries.cl', 'boundaries_periodic_send')
        kernel_ref = octcl_kernel_get_ref(kernel_send)

        call opencl_set_kernel_arg(kernel_ref, 0, maxsend)
        call opencl_set_kernel_arg(kernel_ref, 1, boundaries%buff_nsend)
        call opencl_set_kernel_arg(kernel_ref, 2, boundaries%buff_per_send)
        call opencl_set_kernel_arg(kernel_ref, 3, ffb%pack%buffer)
        call opencl_set_kernel_arg(kernel_ref, 4, log2(ffb%pack%size_real(1)))
        call opencl_set_kernel_arg(kernel_ref, 5, buff_send)

        wgsize = opencl_kernel_workgroup_size(kernel_ref)/ffb%pack%size_real(1)

        call opencl_kernel_run(kernel_ref, (/ffb%pack%size_real(1), pad(maxsend, wgsize), npart/), &
          (/ffb%pack%size_real(1), wgsize, 1/))

        call opencl_finish()

        call opencl_read_buffer(buff_send, ffb%pack%size(1)*maxsend*npart, sendbuffer)
        call opencl_release_buffer(buff_send)
#endif
      end select


      SAFE_ALLOCATE(send_count(1:npart))
      SAFE_ALLOCATE(send_disp(1:npart))
      SAFE_ALLOCATE(recv_count(1:npart))
      SAFE_ALLOCATE(recv_disp(1:npart))

      do ipart = 1, npart
        send_count(ipart) = ldbuffer*boundaries%nsend(ipart)
        send_disp(ipart)  = ldbuffer*maxsend*(ipart - 1)
        recv_count(ipart) = ldbuffer*boundaries%nrecv(ipart)
        recv_disp(ipart)  = ldbuffer*maxrecv*(ipart - 1)
      end do

      ASSERT(send_count(boundaries%mesh%vp%partno) == 0)
      ASSERT(recv_count(boundaries%mesh%vp%partno) == 0)

      SAFE_ALLOCATE(recvbuffer(1:ldbuffer, 1:maxrecv, 1:npart))

      call profiling_out(set_bc_precomm_prof)

      call profiling_in(set_bc_comm_prof, 'SET_BC_COMM')

      call mpi_debug_in(boundaries%mesh%vp%comm, C_MPI_ALLTOALLV)
      call MPI_Alltoallv(sendbuffer, send_count, send_disp, R_MPITYPE, &
        recvbuffer, recv_count, recv_disp, R_MPITYPE, boundaries%mesh%vp%comm, mpi_err)
      call mpi_debug_out(boundaries%mesh%vp%comm, C_MPI_ALLTOALLV)

      call profiling_count_transfers(sum(boundaries%nsend(1:npart) + boundaries%nrecv(1:npart))*ffb%nst_linear, &
        R_TOTYPE(M_ONE))

      call profiling_out(set_bc_comm_prof)

      call profiling_in(set_bc_postcomm_prof, 'SET_BC_POSTCOMM')

      SAFE_DEALLOCATE_A(send_count)
      SAFE_DEALLOCATE_A(send_disp)
      SAFE_DEALLOCATE_A(recv_count)
      SAFE_DEALLOCATE_A(recv_disp)
      SAFE_DEALLOCATE_A(sendbuffer)

      select case(batch_status(ffb))

      case(BATCH_NOT_PACKED)

        do ipart = 1, npart
          do ip = 1, boundaries%nrecv(ipart)
            ip2 = boundaries%per_recv(ip, ipart)
            do ist = 1, ffb%nst_linear
              ffb%states_linear(ist)%X(psi)(ip2) = recvbuffer(ist, ip, ipart)
            end do
          end do
        end do

      case(BATCH_PACKED)

        do ipart = 1, npart
          !$omp parallel do private(ip, ip2, ist)
          do ip = 1, boundaries%nrecv(ipart)
            ip2 = boundaries%per_recv(ip, ipart)
            do ist = 1, ffb%nst_linear
              ffb%pack%X(psi)(ist, ip2) = recvbuffer(ist, ip, ipart)
            end do
          end do
        end do

      case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
        call opencl_create_buffer(buff_recv, CL_MEM_READ_ONLY, R_TYPE_VAL, ffb%pack%size(1)*maxrecv*npart)
        call opencl_write_buffer(buff_recv, ffb%pack%size(1)*maxrecv*npart, recvbuffer)

        call octcl_kernel_start_call(kernel_recv, 'boundaries.cl', 'boundaries_periodic_recv')
        kernel_ref = octcl_kernel_get_ref(kernel_recv)

        call opencl_set_kernel_arg(kernel_ref, 0, maxrecv)
        call opencl_set_kernel_arg(kernel_ref, 1, boundaries%buff_nrecv)
        call opencl_set_kernel_arg(kernel_ref, 2, boundaries%buff_per_recv)
        call opencl_set_kernel_arg(kernel_ref, 3, ubound(boundaries%per_recv, dim = 1))
        call opencl_set_kernel_arg(kernel_ref, 4, buff_recv)
        call opencl_set_kernel_arg(kernel_ref, 5, ffb%pack%buffer)
        call opencl_set_kernel_arg(kernel_ref, 6, log2(ffb%pack%size_real(1)))

        wgsize = opencl_kernel_workgroup_size(kernel_ref)/ffb%pack%size_real(1)

        call opencl_kernel_run(kernel_ref, (/ffb%pack%size_real(1), pad(maxrecv, wgsize), npart/), &
          (/ffb%pack%size_real(1), wgsize, 1/))

        call opencl_finish()

        call opencl_release_buffer(buff_recv)
#endif
      end select

      SAFE_DEALLOCATE_A(recvbuffer)        

      call profiling_out(set_bc_postcomm_prof)

    end if
#endif

    select case(batch_status(ffb))

    case(BATCH_NOT_PACKED)

      do ist = 1, ffb%nst_linear
        ff => ffb%states_linear(ist)%X(psi)
        forall (ip = 1:boundaries%nper)
          ff(boundaries%per_points(POINT_BOUNDARY, ip)) = ff(boundaries%per_points(POINT_INNER, ip))
        end forall
      end do

    case(BATCH_PACKED)

      !$omp parallel do private(ip, ip_bnd, ip_inn, ist)
      do ip = 1, boundaries%nper
        ip_bnd = boundaries%per_points(POINT_BOUNDARY, ip)
        ip_inn = boundaries%per_points(POINT_INNER, ip)
        forall(ist = 1:ffb%nst_linear) ffb%pack%X(psi)(ist, ip_bnd) = ffb%pack%X(psi)(ist, ip_inn)
      end do

    case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL

      call octcl_kernel_start_call(kernel, 'boundaries.cl', 'boundaries_periodic')
      kernel_ref = octcl_kernel_get_ref(kernel)

      call opencl_set_kernel_arg(kernel_ref, 0, boundaries%nper)
      call opencl_set_kernel_arg(kernel_ref, 1, boundaries%buff_per_points)
      call opencl_set_kernel_arg(kernel_ref, 2, ffb%pack%buffer)
      call opencl_set_kernel_arg(kernel_ref, 3, log2(ffb%pack%size_real(1)))

      wgsize = opencl_kernel_workgroup_size(kernel_ref)/ffb%pack%size_real(1)

      call opencl_kernel_run(kernel_ref, (/ffb%pack%size_real(1), pad(boundaries%nper, wgsize)/), &
        (/ffb%pack%size_real(1), wgsize/))

      call opencl_finish()

#endif

    end select

    POP_SUB(X(boundaries_set_batch).periodic)
  end subroutine periodic

end subroutine X(boundaries_set_batch)

! ---------------------------------------------------------

subroutine X(boundaries_set_single)(boundaries, ff)
  type(boundaries_t),  intent(in)    :: boundaries
  R_TYPE, target,      intent(inout) :: ff(:) !< target for batch_add_state

  type(batch_t) :: batch_ff

  PUSH_SUB(X(boundaries_set_single))

  call batch_init     (batch_ff, 1)
  call batch_add_state(batch_ff, ff)

  ASSERT(batch_is_ok(batch_ff))

  call X(boundaries_set_batch)(boundaries, batch_ff)

  call batch_end(batch_ff)
  POP_SUB(X(boundaries_set_single))

end subroutine X(boundaries_set_single)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
