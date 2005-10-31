!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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


! Scatters a vector v to all nodes in vp with respect to
! to point -> node mapping in vp.
! v_local has at least to be of size vp%np_local(vp%partno).
subroutine X(vec_scatter)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(in)  :: v(:)
  R_TYPE,        intent(out) :: v_local(:)

  integer              :: i         ! Counter.
  integer              :: ierr      ! MPI errorcode.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Send buffer.

  call push_sub('par_vec.Xvec_scatter')

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  allocate(displs(vp%p))
  displs = vp%xlocal-1

  ! Fill send buffer.
  if(vp%rank.eq.vp%root) then
    allocate(v_tmp(vp%np))

    ! Rearrange copy of v. All points of node r are in
    ! v_tmp(xlocal(r):xlocal(r)+np_local(r)-1).
    do i = 1, vp%np
      v_tmp(i) = v(vp%local(i))
    end do
  end if

  ! Careful: MPI rank numbers range from 0 to mpiv%numprocs-1
  ! But partition numbers from 1 to vp%p with usually
  ! vp%p = mpiv%numprocs.
  call MPI_Debug_IN(vp%comm, C_MPI_SCATTERV)
  call MPI_Scatterv(v_tmp, vp%np_local, displs, R_MPITYPE, v_local, &
                    vp%np_local(vp%partno), R_MPITYPE,              &
                    vp%root, vp%comm, ierr)
  call MPI_Debug_OUT(vp%comm, C_MPI_SCATTERV)

  if(vp%rank.eq.vp%root) deallocate(v_tmp)

  deallocate(displs)

  call pop_sub()

end subroutine X(vec_scatter)


! v_local has to be of length np_local+np_ghost+np_bndry
! for this to work.
! And v has to be of length np_part.
subroutine X(vec_scatter_bndry)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(in)  :: v(:)
  R_TYPE,        intent(out) :: v_local(:)

  integer              :: i         ! Counter.
  integer              :: ierr      ! MPI errorcode.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Send buffer.

  call push_sub('par_vec.Xvec_scatter_bndry')

  allocate(displs(vp%p))
  displs = vp%xbndry-1

  ! Fill send buffer.
  if(vp%rank.eq.vp%root) then
    allocate(v_tmp(vp%np_enl))

    ! Rearrange copy of v. All points of node r are in
    ! v_tmp(xlocal(r):xlocal(r)+np_local(r)-1).
    do i = 1, vp%np_enl
      v_tmp(i) = v(vp%bndry(i)+vp%np)
    end do
  end if

  ! Careful: MPI rank numbers range from 0 to mpiv%numprocs-1
  ! But partition numbers from 1 to vp%p with usually
  ! vp%p = mpiv%numprocs.
  call MPI_Debug_IN(vp%comm, C_MPI_SCATTERV)
  call MPI_Scatterv(v_tmp, vp%np_bndry, displs, R_MPITYPE,                     &
                    v_local(vp%np_local(vp%partno)+vp%np_ghost(vp%partno)+1:), &
                    vp%np_bndry(vp%partno), R_MPITYPE, vp%root, vp%comm, ierr)
  call MPI_Debug_OUT(vp%comm, C_MPI_SCATTERV)

  if(vp%rank.eq.vp%root) deallocate(v_tmp)

  deallocate(displs)

  call pop_sub()

end subroutine X(vec_scatter_bndry)


! Xvec_scatter followed by Xvec_scatter_bndry.
subroutine X(vec_scatter_all)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(in)  :: v(:)
  R_TYPE,        intent(out) :: v_local(:)

  call push_sub('par_vec.Xvec_scatter_all')

  call X(vec_scatter)(vp, v, v_local)
  call X(vec_scatter_bndry)(vp, v, v_local)

  call pop_sub()

end subroutine X(vec_scatter_all)


! Reverse operation of Xvec_scatter.
! All v_locals from the nodes are packed together
! into v on node vp%root in correct order.
subroutine X(vec_gather)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(out) :: v(:)
  R_TYPE,        intent(in)  :: v_local(:)

  integer              :: i         ! Counter.
  integer              :: ierr      ! MPI errorcode.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Receive buffer.

  call push_sub('par_vec.Xvec_gather')

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  allocate(displs(vp%p))
  displs = vp%xlocal-1

  if(vp%rank.eq.vp%root) allocate(v_tmp(vp%np))

  call MPI_Debug_IN(vp%comm, C_MPI_GATHERV)
  call MPI_Gatherv(v_local, vp%np_local(vp%partno), R_MPITYPE, v_tmp, &
                   vp%np_local, displs, R_MPITYPE,                    &
                   vp%root, vp%comm, ierr)
  call MPI_Debug_OUT(vp%comm, C_MPI_GATHERV)

  ! Copy values from v_tmp to their original position in v.
  if(vp%rank.eq.vp%root) then
    do i = 1, vp%np
      v(vp%local(i)) = v_tmp(i)
    end do

    deallocate(v_tmp)
  end if

  deallocate(displs)

  call pop_sub()

end subroutine X(vec_gather)


! Like Xvec_gather but the result is gathered
! on all nodes, i. e. v has to be a properly
! allocated aray on all nodes.
subroutine X(vec_allgather)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(out) :: v(:)
  R_TYPE,        intent(in)  :: v_local(:)

  integer              :: i         ! Counter.
  integer              :: ierr      ! MPI errorcode.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Receive buffer.

  call push_sub('par_vec.Xvec_allgather')

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  allocate(displs(vp%p))
  displs = vp%xlocal-1

  allocate(v_tmp(vp%np))

  call MPI_Debug_IN(vp%comm, C_MPI_ALLGATHERV)
  call MPI_Allgatherv(v_local, vp%np_local(vp%partno), R_MPITYPE, v_tmp, &
                      vp%np_local, displs, R_MPITYPE,                    &
                      vp%comm, ierr)
  call MPI_Debug_OUT(vp%comm, C_MPI_ALLGATHERV)

  ! Copy values from v_tmp to their original position in v.
  do i = 1, vp%np
    v(vp%local(i)) = v_tmp(i)
  end do

  deallocate(v_tmp)
  deallocate(displs)

  call pop_sub()

end subroutine X(vec_allgather)


! Updates ghost points of every node. A vector suitable
! for non local operations contains local values and
! ghost point values.
! Length of v_local must be
! vp%np_local(vp%partno)+vp%np_ghost(vp%partno)
subroutine X(vec_ghost_update)(vp, v_local)
  type(pv_type), intent(in)    :: vp
  R_TYPE,        intent(inout) :: v_local(:)

  integer              :: i, j, k, r             ! Counters.
  integer              :: ierr                   ! MPI errorcode.
  integer              :: total                  ! Total number of ghost
                                                 ! points to send away.
  integer, allocatable :: sdispls(:), rdispls(:) ! Displacements for
                                                 ! MPI_Alltoallv.
  R_TYPE,  allocatable :: ghost_send(:)          ! Send buffer.

  call profiling_in(C_PROFILING_GHOST_UPDATE)

  call push_sub('par_vec.Xvec_ghost_update')

  ! Calculate number of ghost points current node
  ! has to send to neighbours and allocate send buffer.
  total = sum(vp%np_ghost_neigh(:, vp%partno))
  allocate(ghost_send(total))

  ! Send and receive displacements.
  ! Send displacement cannot directly be calculated
  ! from vp%xghost_neigh because those are indices for
  ! vp%np_ghost_neigh(vp%partno, :) and not
  ! vp%np_ghost_neigh(:, vp%partno) (rank being fixed).
  ! So what gets done is to pick out the number of ghost points
  ! each partition r wants to have from the current partiton
  ! vp%partno.
  allocate(sdispls(vp%p))
  sdispls(1) = 0
  do r = 2, vp%p
    sdispls(r) = sdispls(r-1)+vp%np_ghost_neigh(r-1, vp%partno)
  end do

  ! This is like in vec_scatter/gather.
  allocate(rdispls(vp%p))
  rdispls = vp%xghost_neigh(vp%partno, :)-vp%xghost(vp%partno)

  ! Collect all local points that have to be sent to neighbours.
  j = 1
  ! Iterate over all possible receivers.
  do r = 1, vp%p
    ! Iterate over all ghost points that r wants.
    do i = 0, vp%np_ghost_neigh(r, vp%partno)-1
      ! Get global number k of i-th ghost point.
      k = vp%ghost(vp%xghost_neigh(r, vp%partno)+i)
      ! Lookup up local number of point k and put
      ! value from v_local to send buffer.
      ghost_send(j) = v_local(vp%global(k, vp%partno))
      j             = j + 1
    end do
  end do

  ! Bring it on the way.
  ! It has to examined wheather it is better to use several point to
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
  call MPI_Debug_IN(vp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(ghost_send, vp%np_ghost_neigh(:, vp%partno), sdispls, &
                     R_MPITYPE, v_local(vp%np_local(vp%partno)+1:),        &
                     vp%np_ghost_neigh(vp%partno, :), rdispls, R_MPITYPE,  &
                     vp%comm, ierr)
  call MPI_Debug_OUT(vp%comm, C_MPI_ALLTOALLV)

  deallocate(sdispls, rdispls)
  deallocate(ghost_send)

  call pop_sub()

  call profiling_out(C_PROFILING_GHOST_UPDATE)

end subroutine X(vec_ghost_update)


! Sums over all elements of a vector v which was
! scattered to several v_local with vec_scatter.
! It is of no relevance, if v_local contains
! ghost points or not, because only the first
! vp%np_local(vp%partno) elements are considered.
R_TYPE function X(vec_integrate)(vp, v_local) result(s)
  type(pv_type), intent(in) :: vp
  R_TYPE,        intent(in) :: v_local(:)

  integer :: ierr
  R_TYPE  :: s_local ! Sum of v_local(i).

  call profiling_in(C_PROFILING_VEC_INTEGRATE)
  call push_sub('par_vec.Xvec_integrate')

  s_local = sum(v_local(:vp%np_local(vp%partno)))

  call MPI_Debug_IN(vp%comm, C_MPI_ALLREDUCE)
  call MPI_Allreduce(s_local, s, 1, R_MPITYPE, MPI_SUM, vp%comm, ierr)
  call MPI_Debug_OUT(vp%comm, C_MPI_ALLREDUCE)

  call pop_sub()
  call profiling_out(C_PROFILING_VEC_INTEGRATE)

end function X(vec_integrate)
