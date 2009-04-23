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

! Generally:
! Xvec_gather and Xvec_scatter only consider inner points.
! Xvec_scatter_bndry takes care of boundary points (there is
! no Xvec_gather_bndry as they are only written and not read).
! Xvec_scatter_all is Xvec_scatter followd by Xvec_scatter_bndry.

! ---------------------------------------------------------
! Scatters a vector v to all nodes in vp with respect to
! to point -> node mapping in vp.
! v_local has at least to be of size vp%np_local(vp%partno).
subroutine X(vec_scatter)(vp, v, v_local)
  type(pv_t), intent(in)  :: vp
  R_TYPE,     intent(in)  :: v(:)
  R_TYPE,     intent(out) :: v_local(:)

  integer              :: i         ! Counter.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Send buffer.

  call push_sub('par_vec_inc.Xvec_scatter')

  ! Skip the MPI call if domain parallelization is not used.
  if(vp%npart.lt.2) then
    v_local(1:vp%np) = v(1:vp%np)
    call pop_sub(); return
  end if

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal-1

  SAFE_ALLOCATE(v_tmp(1:1))
  if(vp%rank.eq.vp%root) then
  ! Fill send buffer.
    SAFE_DEALLOCATE_A(v_tmp)
    SAFE_ALLOCATE(v_tmp(1:vp%np))

    ! Rearrange copy of v. All points of node r are in
    ! v_tmp(xlocal(r):xlocal(r)+np_local(r)-1).
    do i = 1, vp%np
      v_tmp(i) = v(vp%local(i))
    end do
  end if

  ! Careful: MPI rank numbers range from 0 to mpiv%numprocs-1
  ! But partition numbers from 1 to vp%npart with usually
  ! vp%npart = mpiv%numprocs.
  call mpi_debug_in(vp%comm, C_MPI_SCATTERV)
  call MPI_Scatterv(v_tmp, vp%np_local, displs, R_MPITYPE, v_local, &
                    vp%np_local(vp%partno), R_MPITYPE,              &
                    vp%root, vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_SCATTERV)

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call pop_sub()

end subroutine X(vec_scatter)


! ---------------------------------------------------------
! v_local has to be of length np_local+np_ghost+np_bndry
! for this to work.
! And v has to be of length np_part.
subroutine X(vec_scatter_bndry)(vp, v, v_local)
  type(pv_t), intent(in)  :: vp
  R_TYPE,     intent(in)  :: v(:)
  R_TYPE,     intent(out) :: v_local(:)

  integer              :: i         ! Counter.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Send buffer.

  call push_sub('par_vec_inc.Xvec_scatter_bndry')

  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xbndry-1

  ! Fill send buffer.
  SAFE_ALLOCATE(v_tmp(1:1))
  if(vp%rank.eq.vp%root) then
    SAFE_DEALLOCATE_A(v_tmp)
    SAFE_ALLOCATE(v_tmp(1:vp%np_enl))

    ! Rearrange copy of v. All points of node r are in
    ! v_tmp(xlocal(r):xlocal(r)+np_local(r)-1).
    do i = 1, vp%np_enl
      v_tmp(i) = v(vp%bndry(i)+vp%np)
    end do
  end if

  ! Careful: MPI rank numbers range from 0 to mpiv%numprocs-1
  ! But partition numbers from 1 to vp%npart with usually
  ! vp%npart = mpiv%numprocs.
  call mpi_debug_in(vp%comm, C_MPI_SCATTERV)
  call MPI_Scatterv(v_tmp, vp%np_bndry, displs, R_MPITYPE,                     &
                    v_local(vp%np_local(vp%partno)+vp%np_ghost(vp%partno)+1:), &
                    vp%np_bndry(vp%partno), R_MPITYPE, vp%root, vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_SCATTERV)

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call pop_sub()

end subroutine X(vec_scatter_bndry)


! ---------------------------------------------------------
! Xvec_scatter followed by Xvec_scatter_bndry.
subroutine X(vec_scatter_all)(vp, v, v_local)
  type(pv_t), intent(in)  :: vp
  R_TYPE,     intent(in)  :: v(:)
  R_TYPE,     intent(out) :: v_local(:)

  call push_sub('par_vec_inc.Xvec_scatter_all')

  call X(vec_scatter)(vp, v, v_local)
  call X(vec_scatter_bndry)(vp, v, v_local)

  call pop_sub()

end subroutine X(vec_scatter_all)


! ---------------------------------------------------------
! Reverse operation of Xvec_scatter.
! All v_locals from the nodes are packed together
! into v on node vp%root in correct order.
subroutine X(vec_gather)(vp, v, v_local)
  type(pv_t), intent(in)  :: vp
  R_TYPE,     intent(out) :: v(:)
  R_TYPE,     intent(in)  :: v_local(:)

  integer              :: i         ! Counter.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Receive buffer.

  call push_sub('par_vec_inc.Xvec_gather')

  ! Skip the MPI call if domain parallelization is not used.
  if(vp%npart.lt.2) then
    v(1:vp%np) = v_local(1:vp%np)
    call pop_sub(); return
  end if

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal-1

  SAFE_ALLOCATE(v_tmp(1:vp%np))

  call mpi_debug_in(vp%comm, C_MPI_GATHERV)
  call MPI_Gatherv(v_local, vp%np_local(vp%partno), R_MPITYPE, v_tmp, &
                   vp%np_local, displs, R_MPITYPE,                    &
                   vp%root, vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_GATHERV)

  ! Copy values from v_tmp to their original position in v.
  if(vp%rank.eq.vp%root) then
    do i = 1, vp%np
      v(vp%local(i)) = v_tmp(i)
    end do

  end if

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call pop_sub()

end subroutine X(vec_gather)


! ---------------------------------------------------------
! Like Xvec_gather but the result is gathered
! on all nodes, i. e. v has to be a properly
! allocated array on all nodes.
subroutine X(vec_allgather)(vp, v, v_local)
  type(pv_t), intent(in)  :: vp
  R_TYPE,     intent(out) :: v(:)
  R_TYPE,     intent(in)  :: v_local(:)

  integer              :: i         ! Counter.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Receive buffer.

  call push_sub('par_vec_inc.Xvec_allgather')

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal-1

  SAFE_ALLOCATE(v_tmp(1:vp%np))

  call mpi_debug_in(vp%comm, C_MPI_ALLGATHERV)
  call MPI_Allgatherv(v_local, vp%np_local(vp%partno), R_MPITYPE, v_tmp, &
                      vp%np_local, displs, R_MPITYPE,                    &
                      vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_ALLGATHERV)

  ! Copy values from v_tmp to their original position in v.
  do i = 1, vp%np
    v(vp%local(i)) = v_tmp(i)
  end do

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call pop_sub()

end subroutine X(vec_allgather)


! ---------------------------------------------------------
! Updates ghost points of every node. A vector suitable
! for non local operations contains local values and
! ghost point values.
! Length of v_local must be
! vp%np_local(vp%partno)+vp%np_ghost(vp%partno)
subroutine X(vec_ghost_update)(vp, v_local)
  type(pv_t), intent(in)    :: vp
  R_TYPE,     intent(inout) :: v_local(:)

  R_TYPE,  pointer :: ghost_send(:)          ! Send buffer.

  call profiling_in(C_PROFILING_GHOST_UPDATE, "GHOST_UPDATE")

  call push_sub('par_vec_inc.Xvec_ghost_update')

  call X(vec_ghost_update_prepare)(vp, v_local, ghost_send)

  ! Bring it on the way.
  ! It has to examined whether it is better to use several point to
  ! point send operations (thus, non-neighbour sends could explicitly be
  ! skipped) or to use Alltoallv. In my opinion, Alltoallv just skips
  ! any I/O if some sendcount is 0. This means, there is actually
  ! no need to code this again.
  ! Note: It is actually possible to do pair-wise communication. Roughly
  ! p/2 pairs can communicate at the same time (as long as there is
  ! a switched network). It then takes as much rounds of
  ! communication as the maximum neighbour number of a particular node
  ! is. If this speeds up communication depends on the MPI
  ! implementation: A performant implementation might already work this
  ! way. So only touch this code if profiling with many processors
  ! shows this spot is a serious bottleneck.
  call mpi_debug_in(vp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(ghost_send, vp%np_ghost_neigh(1, vp%partno), vp%sdispls(1),           &
    R_MPITYPE, v_local(vp%np_local(vp%partno)+1), vp%rcounts(1), vp%rdispls(1), R_MPITYPE, &
    vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_ALLTOALLV)

  call X(vec_ghost_update_finish)(ghost_send)

  call pop_sub()

  call profiling_out(C_PROFILING_GHOST_UPDATE)
end subroutine X(vec_ghost_update)


! ---------------------------------------------------------
! The same as Xvec_ghost_update but in a non-blocking fashion.
! The handle is an NBC_Handle to be used in an NBC_Wait call.
subroutine X(vec_ighost_update)(vp, v_local, handle)
  type(pv_t),         intent(in)    :: vp
  R_TYPE,             intent(inout) :: v_local(:)
  type(pv_handle_t),  intent(inout) :: handle

  integer :: ipart, pos

  call profiling_in(C_PROFILING_GHOST_UPDATE, "GHOST_UPDATE")

  call push_sub('par_vec_inc.Xvec_ighost_update')

  select case(handle%comm_method)
#ifdef HAVE_LIBNBC
  case(NON_BLOCKING_COLLECTIVE)
    ! use a collective non-blocking call
    
    nullify(handle%ighost_send, handle%dghost_send, handle%zghost_send)
    call X(vec_ghost_update_prepare)(vp, v_local, handle%X(ghost_send))
    
    call NBCF_Ialltoallv(handle%X(ghost_send), vp%np_ghost_neigh(1, vp%partno), vp%sdispls(1),  &
         R_MPITYPE, v_local(vp%np_local(vp%partno)+1), vp%rcounts(1), vp%rdispls(1), R_MPITYPE, &
         vp%comm, handle%nbc_h, mpi_err)
    
#endif
  case(NON_BLOCKING)
    ! use a series of p2p non-blocking calls
    
    handle%nnb = 0
    do ipart = 1, vp%npart
      if(vp%np_ghost_neigh(ipart, vp%partno) == 0) cycle
      
      handle%nnb = handle%nnb + 1
      call MPI_Isend(v_local(1), 1, vp%X(send_type)(ipart), ipart - 1, 0, &
           vp%comm, handle%requests(handle%nnb), mpi_err)
      
    end do
    
    do ipart = 1, vp%npart
      if(vp%np_ghost_neigh(vp%partno, ipart) == 0) cycle
      
      handle%nnb = handle%nnb + 1
      pos = vp%np_local(vp%partno) + 1 + vp%rdispls(ipart)
      call MPI_Irecv(v_local(pos), vp%rcounts(ipart), R_MPITYPE, ipart - 1, 0, &
           vp%comm, handle%requests(handle%nnb), mpi_err)
    end do
    
  end select

  call pop_sub()

  call profiling_out(C_PROFILING_GHOST_UPDATE)

end subroutine X(vec_ighost_update)

subroutine X(vec_ghost_update_prepare) (vp, v_local, ghost_send)
  type(pv_t),       intent(in)    :: vp
  R_TYPE,           intent(in)    :: v_local(:)
  R_TYPE,           pointer       :: ghost_send(:)          ! Send buffer.


  integer :: i, j, k, r ! Counters.
  integer :: total      ! Total number of ghost points to send away.

  call push_sub('par_vec_inc.Xvec_ghost_update_prepare')

  
  ! Calculate number of ghost points current node
  ! has to send to neighbours and allocate send buffer.
  total = sum(vp%np_ghost_neigh(:, vp%partno))
  SAFE_ALLOCATE(ghost_send(1:total))
  
  ! Collect all local points that have to be sent to neighbours.
  j = 1
  ! Iterate over all possible receivers.
  do r = 1, vp%npart
    ! Iterate over all ghost points that r wants.
    do i = 0, vp%np_ghost_neigh(r, vp%partno)-1
      ! Get global number k of i-th ghost point.
      k = vp%ghost(vp%xghost_neigh(r, vp%partno)+i)
      ! Lookup up local number of point k and put
      ! value from v_local to send buffer.
      ghost_send(j) = v_local(vec_global2local(vp, k, vp%partno))
      j             = j + 1
    end do
  end do
  
  call pop_sub()
end subroutine X(vec_ghost_update_prepare)

subroutine X(vec_ghost_update_finish)(ghost_send)
  R_TYPE,  pointer :: ghost_send(:)

  SAFE_DEALLOCATE_P(ghost_send)

end subroutine X(vec_ghost_update_finish)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
