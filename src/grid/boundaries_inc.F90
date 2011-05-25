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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

! ---------------------------------------------------------
!> Updates ghost points of every node. A vector suitable
!! for non local operations contains local values and
!! ghost point values.
!! Length of v_local must be
!! vp%np_local(vp%partno)+vp%np_ghost(vp%partno)
subroutine X(vec_ghost_update)(vp, v_local)
  type(pv_t), intent(in)    :: vp
  R_TYPE,     intent(inout) :: v_local(:)

  R_TYPE,  allocatable :: ghost_send(:)
  integer              :: nsend
  
  call profiling_in(prof_update, "GHOST_UPDATE")

  PUSH_SUB(X(vec_ghost_update))

  nsend = subarray_size(vp%sendpoints)
  SAFE_ALLOCATE(ghost_send(1:nsend))
  call X(subarray_gather)(vp%sendpoints, v_local, ghost_send)

  call mpi_debug_in(vp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(ghost_send(1), vp%np_ghost_neigh(1, vp%partno), vp%sdispls(1),           &
    R_MPITYPE, v_local(vp%np_local(vp%partno)+1), vp%rcounts(1), vp%rdispls(1), R_MPITYPE,    &
    vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_ALLTOALLV)

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
#ifdef HAVE_OPENCL
  case(BATCH_CL_PACKED)
    SAFE_ALLOCATE(handle%X(recv_buffer)(1:v_local%pack%size(1)*vp%np_ghost(vp%partno)))

    do ipart = 1, vp%npart
      if(vp%np_ghost_neigh(vp%partno, ipart) == 0) cycle
      
      handle%nnb = handle%nnb + 1
      tag = 0
      pos = 1 + vp%rdispls(ipart)*v_local%pack%size(1)
      call MPI_Irecv(handle%X(recv_buffer)(pos), vp%rcounts(ipart)*v_local%pack%size(1), R_MPITYPE, ipart - 1, tag, &
        vp%comm, handle%requests(handle%nnb), mpi_err)
    end do
#endif

  case(BATCH_PACKED)
    !In this case, data from different vectors is contiguous. So we can use one message per partition.
    do ipart = 1, vp%npart
      if(vp%np_ghost_neigh(vp%partno, ipart) == 0) cycle
      
      handle%nnb = handle%nnb + 1
      tag = 0
      pos = vp%np_local(vp%partno) + 1 + vp%rdispls(ipart)
      call MPI_Irecv(v_local%pack%X(psi)(1, pos), vp%rcounts(ipart)*v_local%pack%size(1), R_MPITYPE, ipart - 1, tag, &
        vp%comm, handle%requests(handle%nnb), mpi_err)
    end do

  case(BATCH_NOT_PACKED)
    do ii = 1, v_local%nst_linear
      do ipart = 1, vp%npart
        if(vp%np_ghost_neigh(vp%partno, ipart) == 0) cycle
        
        handle%nnb = handle%nnb + 1
        tag = ii
        pos = vp%np_local(vp%partno) + 1 + vp%rdispls(ipart)
        call MPI_Irecv(v_local%states_linear(ii)%X(psi)(pos), vp%rcounts(ipart), R_MPITYPE, ipart - 1, tag, &
        vp%comm, handle%requests(handle%nnb), mpi_err)
      end do
    end do

  end select


  call batch_init(handle%ghost_send, 1, v_local%nst_linear)
  call X(batch_new)(handle%ghost_send, 1, v_local%nst_linear, subarray_size(vp%sendpoints))

  if(batch_is_packed(v_local)) call batch_pack(handle%ghost_send, copy = .false.)

  !now collect the data for sending
  call X(subarray_gather_batch)(vp%sendpoints, v_local, handle%ghost_send)

#ifdef HAVE_OPENCL
  if(batch_status(v_local) == BATCH_CL_PACKED) then
    nn = product(handle%ghost_send%pack%size(1:2))
    SAFE_ALLOCATE(handle%X(send_buffer)(1:nn))
    call opencl_read_buffer(handle%ghost_send%pack%buffer, nn, handle%X(send_buffer))
  end if
#endif

  select case(batch_status(v_local))
#ifdef HAVE_OPENCL
  case(BATCH_CL_PACKED)
    do ipart = 1, vp%npart
      if(vp%np_ghost_neigh(ipart, vp%partno) == 0) cycle
      handle%nnb = handle%nnb + 1
      tag = 0
      call MPI_Isend(handle%X(send_buffer)(1 + (vp%sendpos(ipart) - 1)*v_local%pack%size(1)), &
        vp%np_ghost_neigh(ipart, vp%partno)*v_local%pack%size(1), &
        R_MPITYPE, ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
    end do
#endif

  case(BATCH_PACKED)
    do ipart = 1, vp%npart
      if(vp%np_ghost_neigh(ipart, vp%partno) == 0) cycle
      handle%nnb = handle%nnb + 1
      tag = 0
      call MPI_Isend(handle%ghost_send%pack%X(psi)(1, vp%sendpos(ipart)), &
        vp%np_ghost_neigh(ipart, vp%partno)*v_local%pack%size(1), &
        R_MPITYPE, ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
    end do

  case(BATCH_NOT_PACKED)
    do ii = 1, v_local%nst_linear
      do ipart = 1, vp%npart
        if(vp%np_ghost_neigh(ipart, vp%partno) == 0) cycle
        handle%nnb = handle%nnb + 1
        tag = ii
        call MPI_Isend(handle%ghost_send%states_linear(ii)%X(psi)(vp%sendpos(ipart)), vp%np_ghost_neigh(ipart, vp%partno), &
          R_MPITYPE, ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
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

  SAFE_ALLOCATE(status(1:MPI_STATUS_SIZE, 1:handle%nnb))
  call MPI_Waitall(handle%nnb, handle%requests(1), status(1, 1), mpi_err)
  SAFE_DEALLOCATE_A(status)
  SAFE_DEALLOCATE_P(handle%requests)

#ifdef HAVE_OPENCL
  if(batch_status(handle%v_local) == BATCH_CL_PACKED) then
    call opencl_write_buffer(handle%v_local%pack%buffer, handle%v_local%pack%size(1)*handle%vp%np_ghost(handle%vp%partno), &
      handle%X(recv_buffer), offset = handle%v_local%pack%size(1)*handle%vp%np_local(handle%vp%partno))

    SAFE_DEALLOCATE_P(handle%X(send_buffer))
    SAFE_DEALLOCATE_P(handle%X(recv_buffer))
  end if
#endif

  if(batch_is_packed(handle%ghost_send)) call batch_unpack(handle%ghost_send, copy = .false.)
  call batch_end(handle%ghost_send)

  call profiling_out(prof_wait)
  POP_SUB(X(ghost_update_batch_finish))
end subroutine X(ghost_update_batch_finish)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
