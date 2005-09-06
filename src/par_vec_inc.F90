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
! vec_gather and vec_scatter only consider inner points.
! vec_scatter_bndry takes care of boundary points (there is
! no vec_gather_bndry as they are only written and not read).
! vec_scatter_all is vec_scatter followd by vec_scatter_bndry.


! Scatters a vector v to all nodes in vp with respect to 
! to point -> node mapping in vp.
! v_local has at least to be of size vp%np_local(rank+1).
subroutine X(vec_scatter)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(in)  :: v(:)
  R_TYPE,        intent(out) :: v_local(:)

  integer              :: i         ! Counter.
  integer              :: ierr      ! MPI errorcode.
  integer              :: rank      ! Rank of node.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Send buffer.
  
  call push_sub('par_vec.Xvec_scatter')

  call MPI_Comm_rank(vp%comm, rank, ierr)
 
  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  allocate(displs(vp%p))
  displs = vp%xlocal-1

  ! Fill send buffer.
  if(rank.eq.vp%root) then
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
  call MPI_Scatterv(v_tmp, vp%np_local, displs, R_MPITYPE, v_local, &
                    vp%np_local(rank+1), R_MPITYPE,                 &
                    vp%root, vp%comm, ierr)

  if(rank.eq.vp%root) deallocate(v_tmp)

  deallocate(displs)

  call pop_sub()

end subroutine X(vec_scatter)


! v_local has to be of length np_local+np_ghost+np_bndry
! for this to work.
! And v has to be of length np_tot.
subroutine X(vec_scatter_bndry)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(in)  :: v(:)
  R_TYPE,        intent(out) :: v_local(:)

  integer              :: i         ! Counter.
  integer              :: ierr      ! MPI errorcode.
  integer              :: rank      ! Rank of node.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Send buffer.

  call push_sub('par_vec.Xvec_scatter_bndry')
  
  call MPI_Comm_rank(vp%comm, rank, ierr)

  allocate(displs(vp%p))
  displs = vp%xbndry-1

  ! Fill send buffer.
  if(rank.eq.vp%root) then
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
  call MPI_Scatterv(v_tmp, vp%np_bndry, displs, R_MPITYPE,  &
    v_local(vp%np_local(rank+1)+vp%np_ghost(rank+1)+1:),    &
    vp%np_bndry(rank+1), R_MPITYPE, vp%root, vp%comm, ierr)

  if(rank.eq.vp%root) deallocate(v_tmp)

  deallocate(displs)

  call pop_sub()

end subroutine X(vec_scatter_bndry)


! vec_scatter followed by vec_scatter_bndry.
subroutine X(vec_scatter_all)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(in)  :: v(:)
  R_TYPE,        intent(out) :: v_local(:)
  
  call push_sub('par_vec.Xvec_scatter_all')
  
  call X(vec_scatter)(vp, v, v_local)
  call X(vec_scatter_bndry)(vp, v, v_local)

  call pop_sub()

end subroutine X(vec_scatter_all)


! Reverse operation of vec_scatter.
! All v_locals from the nodes are packed together
! into v on node vp%root in correct order.
subroutine X(vec_gather)(vp, v, v_local)
  type(pv_type), intent(in)  :: vp
  R_TYPE,        intent(out) :: v(:)
  R_TYPE,        intent(in)  :: v_local(:)

  integer              :: i         ! Counter.
  integer              :: ierr      ! MPI errorcode.
  integer              :: rank      ! Rank of node.
  integer, allocatable :: displs(:) ! Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  ! Receive buffer.
  
  call push_sub('par_vec.Xvec_gather')
  
  call MPI_Comm_rank(vp%comm, rank, ierr)

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  allocate(displs(vp%p))
  displs = vp%xlocal-1

  if(rank.eq.vp%root) allocate(v_tmp(vp%np))

  call MPI_Gatherv(v_local, vp%np_local(rank+1), R_MPITYPE, v_tmp, &
                   vp%np_local, displs, R_MPITYPE,                 &
                   vp%root, vp%comm, ierr)

  ! Copy values from v_tmp to their original position in v.
  if(rank.eq.vp%root) then
    do i = 1, vp%np
      v(vp%local(i)) = v_tmp(i)
    end do
    
    deallocate(v_tmp)
  end if

  deallocate(displs)

  call pop_sub()

end subroutine X(vec_gather)


! Updates ghost points of every node. A vector suitable
! for non local operations contains local values and
! ghost point values.
! Length of v_local must be
! vp%np_local(rank+1)+vp%np_ghost(rank+1)
subroutine X(vec_ghost_update)(vp, v_local)
  type(pv_type), intent(in)    :: vp
  R_TYPE,        intent(inout) :: v_local(:)
  
  integer              :: i, j, k, r             ! Counters.
  integer              :: ierr                   ! MPI errorcode.
  integer              :: rank                   ! Rank of current node.
  integer              :: total                  ! Total number of ghost
                                                 ! points to send away.
  integer, allocatable :: sdispls(:), rdispls(:) ! Displacements for
                                                 ! MPI_Alltoallv.
  R_TYPE,  allocatable :: ghost_send(:)          ! Send buffer.

  call push_sub('par_vec.Xvec_ghost_update')

  call MPI_Comm_rank(vp%comm, rank, ierr)

  ! Calculate number of ghost points current node
  ! has to send to neighbours and allocate send buffer.
  total = sum(vp%np_ghost_neigh(:, rank+1))
  allocate(ghost_send(total))
  
  ! Send and receive displacements.
  ! Send displacement cannot directly be calculated
  ! from vp%xghost_neigh because those are indices for
  ! vp%np_ghost_neigh(rank+1, :) and not
  ! vp%np_ghost_neigh(:, rank+1) (rank being fixed).
  ! So what gets done is to pick out the number of ghost points
  ! each partition r wants to have from the current partiton
  ! rank+1.
  allocate(sdispls(vp%p))
  sdispls(1) = 0
  do r = 2, vp%p
    sdispls(r) = sdispls(r-1)+vp%np_ghost_neigh(r-1, rank+1)
  end do
  
  ! This is like in vec_scatter/gather.
  allocate(rdispls(vp%p))
  rdispls = vp%xghost_neigh(rank+1, :)-vp%xghost(rank+1)

  ! Collect all local points that have to be sent to neighbours.
  j = 1
  ! Iterate over all possible receivers.
  do r = 1, vp%p
    ! Iterate over all ghost points that r wants.
    do i = 0, vp%np_ghost_neigh(r, rank+1)-1
      ! Get global number k of i-th ghost point.
      k = vp%ghost(vp%xghost_neigh(r, rank+1)+i)
      ! Lookup up local number of point k and put
      ! value from v_local to send buffer.
      ghost_send(j) = v_local(vp%global(k, rank+1))
      j             = j + 1
    end do
  end do
  
  ! Bring it on the way.
  ! It has to examined wheather it is better to use several point to
  ! point send operations (thus, non-neighbour sends could explicitly be
  ! skipped) or to use Alltoallv. In my opinion, Alltoallv just skips
  ! any I/O if some sendcount is 0. This means, there is actually
  ! no need to code this again.
  ! FIXME: It is actually possible to do pair-wise communication. Roughly
  ! p/2 pairs can communicate at the same time (as long as there is
  ! a switched network). It then takes as much rounds of
  ! communication as the maximum neighbour number of a particular node
  ! is. This will speed up exchange and will probably be implemented later.
  call MPI_Alltoallv(ghost_send, vp%np_ghost_neigh(:, rank+1), sdispls, &
                     R_MPITYPE, v_local(vp%np_local(rank+1)+1:),        &
                     vp%np_ghost_neigh(rank+1, :), rdispls, R_MPITYPE,  &
                     vp%comm, ierr)
                     
  deallocate(sdispls, rdispls)
  deallocate(ghost_send)
  
  call pop_sub()
  
end subroutine X(vec_ghost_update)


! Sums over all elements of a vector v which was
! scattered to several v_local with vec_scatter.
! It is of no relevance, if v_local contains
! ghost points or not, because only the first
! vp%np_local(rank+1) elements are considered.
R_TYPE function X(vec_integrate)(vp, v_local) result(s)
  type(pv_type), intent(in) :: vp
  R_TYPE,        intent(in) :: v_local(:)

  integer :: rank
  integer :: ierr
  R_TYPE  :: s_local ! Sum of v_local(i).

  call push_sub('par_vec.Xvec_integrate')
  
  call MPI_Comm_rank(vp%comm, rank, ierr)

  s_local = sum(v_local(:vp%np_local(rank+1)))

  call MPI_Allreduce(s_local, s, 1, R_MPITYPE, MPI_SUM, vp%comm, ierr)

  call pop_sub()

end function X(vec_integrate)
