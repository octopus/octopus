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

!> Generally:
!! Xvec_gather and Xvec_scatter only consider inner points.
!! Xvec_scatter_bndry takes care of boundary points (there is
!! no Xvec_gather_bndry as they are only written and not read).
!! Xvec_scatter_all is Xvec_scatter followd by Xvec_scatter_bndry.

!! ---------------------------------------------------------
!! Scatters a vector v to all nodes in vp with respect to
!! to point -> node mapping in vp.
!! v_local has at least to be of size vp%np_local(vp%partno).
subroutine X(vec_scatter)(vp, root, v, v_local)
  type(pv_t), intent(in)  :: vp
  integer,    intent(in)  :: root
  R_TYPE,     intent(in)  :: v(:)
  R_TYPE,     intent(out) :: v_local(:)

  integer              :: i         !< Counter.
  integer, allocatable :: displs(:) !< Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  !< Send buffer.

  PUSH_SUB(X(vec_scatter))
  call profiling_in(prof_scatter, "VEC_SCATTER")

  ! Skip the MPI call if domain parallelization is not used.
  if(vp%npart.lt.2) then
    v_local(1:vp%np) = v(1:vp%np)
    POP_SUB(X(vec_scatter))
    return
  end if

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal - 1

  SAFE_ALLOCATE(v_tmp(1:1))
  if(vp%rank.eq.root) then
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
  call MPI_Scatterv(v_tmp(1), vp%np_local, displs(1), R_MPITYPE, v_local(1), &
                    vp%np_local(vp%partno), R_MPITYPE,              &
                    root, vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_SCATTERV)

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call profiling_out(prof_scatter)
  POP_SUB(X(vec_scatter))

end subroutine X(vec_scatter)


! ---------------------------------------------------------
!> v_local has to be of length np_local+np_ghost+np_bndry
!! for this to work.
!! And v has to be of length np_part.
subroutine X(vec_scatter_bndry)(vp, root, v, v_local)
  type(pv_t), intent(in)  :: vp
  integer,    intent(in)  :: root
  R_TYPE,     intent(in)  :: v(:)
  R_TYPE,     intent(out) :: v_local(:)

  integer              :: i         !< Counter.
  integer, allocatable :: displs(:) !< Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  !< Send buffer.

  PUSH_SUB(X(vec_scatter_bndry))

  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xbndry-1

  ! Fill send buffer.
  SAFE_ALLOCATE(v_tmp(1:1))
  if(vp%rank.eq.root) then
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
                    vp%np_bndry(vp%partno), R_MPITYPE, root, vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_SCATTERV)

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  POP_SUB(X(vec_scatter_bndry))

end subroutine X(vec_scatter_bndry)


! ---------------------------------------------------------
! Xvec_scatter followed by Xvec_scatter_bndry.
subroutine X(vec_scatter_all)(vp, root, v, v_local)
  type(pv_t), intent(in)  :: vp
  integer,    intent(in)  :: root
  R_TYPE,     intent(in)  :: v(:)
  R_TYPE,     intent(out) :: v_local(:)

  PUSH_SUB(X(vec_scatter_all))

  call X(vec_scatter)(vp, root, v, v_local)
  call X(vec_scatter_bndry)(vp, root, v, v_local)

  POP_SUB(X(vec_scatter_all))

end subroutine X(vec_scatter_all)


! ---------------------------------------------------------
!> Reverse operation of Xvec_scatter.
!! All v_locals from the nodes are packed together
!! into v on node root in correct order.
subroutine X(vec_gather)(vp, root, v, v_local)
  type(pv_t), intent(in)  :: vp
  integer,    intent(in)  :: root
  R_TYPE,     intent(out) :: v(:)
  R_TYPE,     intent(in)  :: v_local(:)

  integer              :: i         !< Counter.
  integer, allocatable :: displs(:) !< Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  !< Receive buffer.

  PUSH_SUB(X(vec_gather))

  ! Skip the MPI call if domain parallelization is not used.
  if(vp%npart.lt.2) then
    v(1:vp%np) = v_local(1:vp%np)
    POP_SUB(X(vec_gather))
    return
  end if

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal-1

  SAFE_ALLOCATE(v_tmp(1:vp%np))

  call mpi_debug_in(vp%comm, C_MPI_GATHERV)
  call MPI_Gatherv(v_local(1), vp%np_local(vp%partno), R_MPITYPE, v_tmp(1), &
                   vp%np_local, displs(1), R_MPITYPE,                        &
                   root, vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_GATHERV)

  ! Copy values from v_tmp to their original position in v.
  if(vp%rank.eq.root) then
    do i = 1, vp%np
      v(vp%local(i)) = v_tmp(i)
    end do

  end if

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  POP_SUB(X(vec_gather))

end subroutine X(vec_gather)

! ---------------------------------------------------------
!> Like Xvec_gather but the result is gathered
!! on all nodes, i. e. v has to be a properly
!! allocated array on all nodes.
subroutine X(vec_allgather)(vp, v, v_local)
  type(pv_t), intent(in)  :: vp
  R_TYPE,     intent(out) :: v(:)
  R_TYPE,     intent(in)  :: v_local(:)

  integer              :: i         !< Counter.
  integer, allocatable :: displs(:) !< Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  !< Receive buffer.

  PUSH_SUB(X(vec_allgather))
  call profiling_in(prof_allgather, "VEC_ALLGATHER")
  
  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal-1

  SAFE_ALLOCATE(v_tmp(1:vp%np))

  call mpi_debug_in(vp%comm, C_MPI_ALLGATHERV)
  call MPI_Allgatherv(v_local(1), vp%np_local(vp%partno), R_MPITYPE, v_tmp(1), &
                      vp%np_local, displs(1), R_MPITYPE,                       &
                      vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_ALLGATHERV)

  ! Copy values from v_tmp to their original position in v.
  do i = 1, vp%np
    v(vp%local(i)) = v_tmp(i)
  end do

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call profiling_out(prof_allgather)
  POP_SUB(X(vec_allgather))

end subroutine X(vec_allgather)

!--------------------------------------------------------

!> This function collects points from the array src in nodes and puts
!! them in the arrat dst in the node with MPI rank root, the points to
!! collect are given by a list of global indices (of size nn).
subroutine X(vec_selective_gather)(this, nn, list, root, src, dst)
  type(pv_t),       intent(in)    :: this    !< The par_vec_t object.
  integer,          intent(in)    :: nn      !< The number of points in the list, must be the same on all nodes.
  integer,          intent(in)    :: list(:) !< The list of global points that will be collected, must be the same on all nodes.
  integer,          intent(in)    :: root    !< The node that will receive the data, in MPI convention.
  R_TYPE,           intent(in)    :: src(:)  !< The points, defined in all nodes.
  R_TYPE,           intent(out)   :: dst(:)  !< The destination array, only defined on the root node.
  
  integer, allocatable :: recv_count(:), copy_count(:)
  R_TYPE, allocatable  :: send_buffer(:), recv_buffer(:, :)
  integer :: send_count, ip, owner, ii, ipart
  logical :: i_am_root
  integer, save :: tag = 0

  i_am_root = (this%rank == root)

  if(i_am_root) then
    SAFE_ALLOCATE(recv_count(1:this%npart))
    SAFE_ALLOCATE(copy_count(1:this%npart))
    recv_count(1:this%npart) = 0
  else
    SAFE_ALLOCATE(send_buffer(1:nn))
    send_count = 0
  end if

  do ii = 1, nn
    ip = list(ii)
    owner = this%part(ip)

    if(owner == this%partno .and. i_am_root) then
      ! it is here, copy directly
      dst(ii) = src(vec_global2local(this, ip, owner))
    else if(owner == this%partno) then
      ! prepare to send
      send_count = send_count + 1
      send_buffer(send_count) = src(vec_global2local(this, ip, owner))
    else if(i_am_root) then
      ! prepare to recv the point
      recv_count(owner) = recv_count(owner) + 1
    end if
  end do

  if(i_am_root) then
    SAFE_ALLOCATE(recv_buffer(1:maxval(recv_count), 1:this%npart))

    !recv from all nodes
    do ipart = 1, this%npart
      if(ipart == this%partno) cycle
      if(recv_count(ipart) == 0) cycle
      call MPI_Recv(recv_buffer(1, ipart), recv_count(ipart), R_MPITYPE, &
        ipart - 1, tag, this%comm, MPI_STATUS_IGNORE, mpi_err)
    end do

    !now put the values in the right place
    copy_count = 0
    do ii = 1, nn
      owner = this%part(list(ii))
      if(owner /= this%partno) then
        copy_count(owner) = copy_count(owner) + 1
        dst(ii) = recv_buffer(copy_count(owner), owner)
      end if
    end do

    ASSERT(all(copy_count == recv_count))

    SAFE_DEALLOCATE_A(recv_buffer)
  else
    if(send_count > 0) then
      ! just send the values
      call MPI_Send(send_buffer(1), send_count, R_MPITYPE, root, tag, this%comm, mpi_err)
    end if
    SAFE_DEALLOCATE_A(send_buffer)
  end if
  
  tag = tag + 1

end subroutine X(vec_selective_gather)

! --------------------------------------------------------------------

subroutine X(vec_selective_scatter)(this, nn, list, root, src, dst)
  type(pv_t),       intent(in)    :: this
  integer,          intent(in)    :: nn
  integer,          intent(in)    :: list(:)
  integer,          intent(in)    :: root
  R_TYPE,           intent(in)    :: src(:)
  R_TYPE,           intent(out)   :: dst(:)
 
  integer, allocatable :: send_count(:)
  R_TYPE, allocatable  :: send_buffer(:), recv_buffer(:)
  integer :: recv_count, ip, owner, ii, ipart, jj
  logical :: i_am_root
  integer, save :: tag = 0

  i_am_root = (this%rank == root)

  if(i_am_root) then
    SAFE_ALLOCATE(send_count(1:this%npart))
    send_count(1:this%npart) = 0
  else
    recv_count = 0
  end if

  do ii = 1, nn
    ip = list(ii)
    owner = this%part(ip)

    if(owner == this%partno .and. i_am_root) then
      ! it is here, copy directly
      dst(vec_global2local(this, ip, owner)) = src(ii)
    else if(owner == this%partno) then
      recv_count = recv_count + 1
    else if(i_am_root) then
      send_count(owner) = send_count(owner) + 1
    end if
  end do

  if(i_am_root) then
    SAFE_ALLOCATE(send_buffer(1:maxval(send_count)))

    !copy to a buffer and send to all nodes
    do ipart = 1, this%npart
      if(ipart == this%partno) cycle
      if(send_count(ipart) == 0) cycle
      jj = 0
      do ii = 1, nn
        if(ipart /= this%part(list(ii))) cycle
        jj = jj + 1
        send_buffer(jj) = src(ii)
      end do

      ASSERT(jj == send_count(ipart))
      call MPI_Send(send_buffer(1), send_count(ipart), R_MPITYPE, ipart - 1, tag, this%comm, mpi_err)
    end do

    SAFE_DEALLOCATE_A(send_buffer)
    SAFE_DEALLOCATE_A(send_count)
  else
    SAFE_ALLOCATE(recv_buffer(1:recv_count))
    if(recv_count > 0) then

      call MPI_Recv(recv_buffer(1), recv_count, R_MPITYPE, root, tag, this%comm, MPI_STATUS_IGNORE, mpi_err)

      jj = 0
      do ii = 1, nn
        ip = list(ii)
        owner = this%part(ip)
        if(owner /= this%partno) cycle
        jj = jj + 1
        dst(vec_global2local(this, ip, owner)) = recv_buffer(jj)
      end do

      ASSERT(jj == recv_count)

    end if
    SAFE_DEALLOCATE_A(recv_buffer)
  end if
  
  tag = tag + 1

end subroutine X(vec_selective_scatter)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
