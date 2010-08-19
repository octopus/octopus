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
!! $Id: mesh_alltoall_inc.F90 6794 2010-07-08 08:54:56Z xavier $

! ---------------------------------------------------------

subroutine X(mesh_alltoall_batch_start)(this, send, recv_offset, recv)
  type(mesh_alltoall_t),    intent(inout) :: this
  type(batch_t),            intent(in)    :: send
  type(batch_t),            intent(inout) :: recv
  integer, optional,        intent(in)    :: recv_offset

  integer :: offset
  integer :: ipart, pos, ii, tag

  call push_sub('mesh_alltoall_inc.Xmesh_alltoall_batch_start')

  ASSERT(send%nst_linear > 0)
  ASSERT(batch_status(send) == batch_status(recv))
  ASSERT(.not. this%in_progress)
  this%in_progress = .true.

  offset = 0
  if(present(recv_offset)) offset = recv_offset

  this%nrequests = 0
  SAFE_ALLOCATE(this%requests(1:2*this%mesh%vp%npart*send%nst_linear))

  ! first post the receptions
  select case(batch_status(send))
  case(BATCH_CL_PACKED)
    ASSERT(.false.)

  case(BATCH_PACKED)
    !In this case, data from different vectors is contiguous. So we can use one message per partition.
    pos = 1
    do ipart = 1, this%mesh%vp%npart
      if(this%nrecv(ipart) == 0) cycle
      INCR(this%nrequests, 1)
      tag = 0
      call MPI_Irecv(recv%pack%X(psi)(1, pos + offset), this%nrecv(ipart)*recv%pack%size(1), R_MPITYPE, ipart - 1, tag, &
        this%mesh%mpi_grp%comm, this%requests(this%nrequests), mpi_err)
      INCR(pos, this%nrecv(ipart))
    end do

  case(BATCH_NOT_PACKED)
    do ii = 1, send%nst_linear
      pos = 1
      do ipart = 1, this%mesh%vp%npart
        if(this%nrecv(ipart) == 0) cycle
        INCR(this%nrequests, 1)
        tag = ii
        call MPI_Irecv(send%states_linear(ii)%X(psi)(pos + offset), this%nrecv(ipart), R_MPITYPE, ipart - 1, tag, &
          this%mesh%mpi_grp%comm, this%requests(this%nrequests), mpi_err)
        INCR(pos, this%nrecv(ipart))
      end do
    end do

  end select

  call batch_init(this%send_buffer, 1, send%nst_linear)
  call X(batch_new)(this%send_buffer, 1, send%nst_linear, subarray_size(this%sendpoints))
  if(batch_is_packed(send)) call batch_pack(this%send_buffer, copy = .false.)

  !now collect the data for sending
  call X(subarray_gather_batch)(this%sendpoints, send, this%send_buffer)

  select case(batch_status(send))
  case(BATCH_CL_PACKED)
    ASSERT(.false.)

  case(BATCH_PACKED)
    pos = 1
    do ipart = 1, this%mesh%vp%npart
      if(this%nsend(ipart) == 0) cycle
      INCR(this%nrequests, 1)
      tag = 0
      call MPI_Isend(this%send_buffer%pack%X(psi)(1, pos), this%nsend(ipart)*send%pack%size(1), &
        R_MPITYPE, ipart - 1, tag, this%mesh%mpi_grp%comm, this%requests(this%nrequests), mpi_err)
      INCR(pos, this%nsend(ipart))
    end do
    
  case(BATCH_NOT_PACKED)
    do ii = 1, send%nst_linear
      pos = 1
      do ipart = 1, this%mesh%vp%npart
        if(this%nsend(ipart) == 0) cycle
        INCR(this%nrequests, 1)
        tag = ii
        call MPI_Isend(this%send_buffer%states_linear(ii)%X(psi)(pos), this%nsend(ipart), &
          R_MPITYPE, ipart - 1, tag, this%mesh%mpi_grp%comm, this%requests(this%nrequests), mpi_err)
        INCR(pos, this%nsend(ipart))
      end do
    end do
  end select

  call pop_sub('mesh_alltoall_inc.Xmesh_alltoall_batch_start')
end subroutine X(mesh_alltoall_batch_start)

! -------------------------------------------------------------------

subroutine X(mesh_alltoall_batch_finish)(this)
  type(mesh_alltoall_t),    intent(inout) :: this

  integer, allocatable :: status(:, :)

  call push_sub('mesh_alltoall_inc.Xmesh_alltoall_batch_finish')
  
  ASSERT(this%in_progress)
  ASSERT(this%nrequests > 0)
  this%in_progress = .false.

  SAFE_ALLOCATE(status(1:MPI_STATUS_SIZE, 1:this%nrequests))
  call MPI_Waitall(this%nrequests, this%requests(1), status(1, 1), mpi_err)
  SAFE_DEALLOCATE_A(status)
  SAFE_DEALLOCATE_P(this%requests)

  if(batch_is_packed(this%send_buffer)) call batch_unpack(this%send_buffer, copy = .false.)
  call X(batch_delete)(this%send_buffer)
  call batch_end(this%send_buffer)

  call pop_sub('mesh_alltoall_inc.Xmesh_alltoall_batch_finish')
end subroutine X(mesh_alltoall_batch_finish)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
