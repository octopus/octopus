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

  call push_sub('par_vec_inc.Xvec_ghost_update')

  nsend = subarray_size(vp%sendpoints)
  SAFE_ALLOCATE(ghost_send(1:nsend))
  call X(subarray_gather)(vp%sendpoints, v_local, ghost_send)

  call mpi_debug_in(vp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(ghost_send(1), vp%np_ghost_neigh(1, vp%partno), vp%sdispls(1),           &
    R_MPITYPE, v_local(vp%np_local(vp%partno)+1), vp%rcounts(1), vp%rdispls(1), R_MPITYPE,    &
    vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_ALLTOALLV)

  SAFE_DEALLOCATE_A(ghost_send)

  call pop_sub('par_vec_inc.Xvec_ghost_update')

  call profiling_out(prof_update)
end subroutine X(vec_ghost_update)

! ---------------------------------------------------------

subroutine X(ghost_update_batch_start)(vp, v_local, handle)
  type(pv_t),               intent(in)    :: vp
  type(batch_t),            intent(inout) :: v_local
  type(pv_handle_batch_t),  intent(out)   :: handle

  integer :: ipart, pos, nsend, ii, tag

  call profiling_in(prof_start, "GHOST_UPDATE_START")
  call push_sub('par_vec_inc.Xghost_update_batch_start')

  ASSERT(v_local%nst_linear > 0)
  ASSERT(.not. batch_is_packed(v_local))

  call batch_init(handle%ghost_send, 1, v_local%nst_linear)

  nsend = subarray_size(vp%sendpoints)

  call X(batch_new)(handle%ghost_send, 1, v_local%nst_linear, nsend)

  handle%nnb = 0

  ! first post the receptions
  
  SAFE_ALLOCATE(handle%requests(1:2*vp%npart*v_local%nst_linear))
  
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
  
  !now pack the data for sending
  call X(subarray_gather_batch)(vp%sendpoints, v_local, handle%ghost_send)
  
  do ii = 1, v_local%nst_linear
    do ipart = 1, vp%npart
      if(vp%np_ghost_neigh(ipart, vp%partno) == 0) cycle
      handle%nnb = handle%nnb + 1
      tag = ii
      call MPI_Isend(handle%ghost_send%states_linear(ii)%X(psi)(vp%sendpos(ipart)), vp%np_ghost_neigh(ipart, vp%partno), &
        R_MPITYPE, ipart - 1, tag, vp%comm, handle%requests(handle%nnb), mpi_err)
    end do
  end do

  call pop_sub('par_vec_inc.Xghost_update_batch_start')
  call profiling_out(prof_start)

end subroutine X(ghost_update_batch_start)

! ---------------------------------------------------------

subroutine X(ghost_update_batch_finish)(handle)
  type(pv_handle_batch_t),  intent(inout)   :: handle

  integer, allocatable :: status(:, :)

  call profiling_in(prof_wait, "GHOST_UPDATE_WAIT")
  call push_sub('par_vec_inc.Xghost_update_batch_finish')
  
  ASSERT(handle%nnb > 0)

  SAFE_ALLOCATE(status(1:MPI_STATUS_SIZE, 1:handle%nnb))
  call MPI_Waitall(handle%nnb, handle%requests(1), status(1, 1), mpi_err)
  SAFE_DEALLOCATE_A(status)
  SAFE_DEALLOCATE_P(handle%requests)

  call X(batch_delete)(handle%ghost_send)
  call batch_end(handle%ghost_send)

  call profiling_out(prof_wait)
  call pop_sub('par_vec_inc.Xghost_update_batch_finish')
end subroutine X(ghost_update_batch_finish)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
