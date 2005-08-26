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

#include "global.h"

module par_vec

  ! Usage example for par_vec routines.
  !
  ! integer            :: np_local, np_ghost
  ! FLOAT              :: s
  ! FLOAT              :: u(np), v(np)
  ! FLOAT, allocatable :: ul(:), vl(:), w(:)
  ! type(pv_type)      :: vp
  ! type(mesh_type)    :: m
  !
  ! ! Fill u, v with sensible values.
  ! ! ...
  !
  ! ! Initialize parallelization with e. g.
  ! ! m          = sys%gr%m
  ! ! stencil    = op%stencil
  ! ! np_stencil = op%n
  ! call vec_init_default(m, stencil, np_stencil, vp, np_local, np_ghost)
  !
  ! ! Allocate space for local vectors.
  ! allocate(ul(np_local+np_ghost))
  ! allocate(vl(np_local+np_ghost))
  ! allocate(wl(np_local+np_ghost))
  !
  ! ! Distribute vectors.
  ! call vec_scatter(vp, u, ul)
  ! call vec_scatter(vp, v, vl)
  !
  ! ! Calculate scalar product.
  ! wl(:np_local) = ul*vl
  ! ! vec_integrate ignores ghost points (i. e. v(np_local+1:)).
  ! s = vec_integrate(vp, wl)
  !
  ! ! Compute some operator op: vl = op ul
  ! call vec_ghost_update(vp, ul)
  ! call vec_nl_operator(vp, op, ul, vl)
  ! ! Gather result of op in one vector v.
  ! call vec_gather(vp, v, vl)
  !
  ! ! Clean up.
  ! call vec_end(vp)

  use global
  use mesh
  use messages
  use nl_operator
#ifdef DEBUG
  use io
#endif

  implicit none

  private
  public :: pv_type,          &
            vec_init,         &
            vec_init_default, &
            vec_end,          &
            vec_scatter,      &
            vec_gather,       &
            vec_ghost_update, &
            vec_integrate,    &
            vec_dnl_operator


  ! Describes mesh distribution to nodes.
  ! There is actually some redundancy: p and np are
  ! also stored in type(mesh_type) but it is needed in
  ! vec-routines too and I see no use in additionally passing
  ! a type(mesh_type) argument in all vec-calls.
  type pv_type
    integer          :: comm                 ! MPI communicator to use.
    integer          :: root                 ! The master node.
    integer          :: p                    ! Number of partitions.
    integer          :: np                   ! Total number of points.
    integer, pointer :: np_local(:)          ! How many points has partition r?
    integer, pointer :: xlocal(:)            ! Points of partition r start at
                                             ! xlocal(r) in local.
    integer, pointer :: local(:)             ! Partition r has points
                                             ! local(xlocal(r):
                                             ! xlocal(r)+np_local(r)-1).
    integer, pointer :: global(:, :)         ! global(i, r) is local number
                                             ! of point i in partition r
                                             ! (if this is 0, i is neither
                                             ! a ghost point nor local to r).
    integer          :: total                ! Total number of ghost points.
    integer, pointer :: np_ghost(:)          ! How many ghost points has
                                             ! partition r?
    integer, pointer :: np_ghost_neigh(:, :) ! Number of ghost points per
                                             ! neighbour per partition.
    integer, pointer :: xghost(:)            ! Like xpoint.
    integer, pointer :: xghost_neigh(:, :)   ! Like xghost for neighbours.
    integer, pointer :: ghost(:)             ! Global indices of all local
                                             ! ghost points.
  end type pv_type


contains

  ! Wrapper: calls vec_init with MPI_COMM_WORLD and master node 0.
  subroutine vec_init_default(m, stencil, np_stencil, vp, &
                      np_local, np_ghost)
    type(mesh_type), intent(in)  :: m            ! The current mesh.
    integer,         intent(in)  :: stencil(:,:) ! The stencil for which to
                                                 ! calculate ghost points.
    integer,         intent(in)  :: np_stencil   ! Num. of points in stencil.
    type(pv_type),   intent(out) :: vp           ! Description of partition.
    ! Those two are shortcuts.
    integer,         intent(out) :: np_local     ! vp%np_local(rank).
    integer,         intent(out) :: np_ghost     ! vp%np_ghost(rank).

    call push_sub('par_vec.vec_init_default')
    
    call vec_init(MPI_COMM_WORLD, 0, m, stencil, np_stencil, vp, &
                  np_local, np_ghost)

    call pop_sub()

  end subroutine vec_init_default
 

  ! Initializes a pv_type object (parallel vector).
  ! It computes the local to global and global to local index tables
  ! and the ghost point exchange.
  ! The format for the stencil is: stencil(i, 3) for i=1, ..., np_stencil.
  ! For example a stencil like (in x-y-plane)
  !          .
  !        .....
  !          .
  ! is coded as
  !   stencil(1, :) = (/ 0,  1,  0/)
  !   stencil(2, :) = (/-2,  0,  0/)
  !   stencil(3, :) = (/-1,  0,  0/)
  !   stencil(4, :) = (/ 0,  0,  0/)
  !   stencil(5, :) = (/ 1,  0,  0/)
  !   stencil(6, :) = (/ 2,  0,  0/)
  !   stencil(7, :) = (/ 0, -1,  0/)
  ! The points are relative to the "application point" of the stencil.
  subroutine vec_init(comm, root, m, stencil, np_stencil, vp, &
                      np_local, np_ghost)
    integer,         intent(in)  :: comm         ! Communicator to use.
    integer,         intent(in)  :: root         ! The master node.
    type(mesh_type), intent(in)  :: m            ! The current mesh.
    integer,         intent(in)  :: stencil(:,:) ! The stencil for which to
                                                 ! calculate ghost points.
    integer,         intent(in)  :: np_stencil   ! Num. of points in stencil.
    type(pv_type),   intent(out) :: vp           ! Description of partition.
    ! Those two are shortcuts.
    integer,         intent(out) :: np_local     ! vp%np_local(rank).
    integer,         intent(out) :: np_ghost     ! vp%np_ghost(rank).
    
    ! Careful: MPI counts node ranks from 0 to numproc-1.
    ! Partition numbers from METIS range from 1 to numproc.
    ! For this reason, all ranks are incremented by one.
    integer              :: p                ! Number of partitions.
    integer              :: np               ! Number of points.
    integer              :: i, j, k, r       ! Counters.
    integer, allocatable :: ir(:), irr(:, :) ! Counters.
    integer              :: rank             ! Rank of current node.
    integer              :: ierr             ! MPI errorcode.
    integer              :: p1(3), p2(3)     ! Points.
    integer, allocatable :: ghost_flag(:, :) ! To remember ghost pnts.
#ifdef DEBUG
    integer              :: iunit            ! For debug output to files.
    character(len=3)     :: filenum
#endif

    call push_sub('par_vec.vec_init')

    p  = m%npart
    np = m%np
    call MPI_Comm_rank(comm, rank, ierr)
    call mpierr(ierr)
    
    allocate(ghost_flag(np, p))
    allocate(ir(p), irr(p, p))

    allocate(vp%np_local(p))
    allocate(vp%xlocal(p))
    allocate(vp%local(np))
    allocate(vp%global(np, p))
    allocate(vp%np_ghost(p))
    allocate(vp%np_ghost_neigh(p, p))
    allocate(vp%xghost(p))
    allocate(vp%xghost_neigh(p, p))
    
    ! Count number of points for each node.
    do i = 1, np
      vp%np_local(m%part(i)) = vp%np_local(m%part(i))+1
    end do

    ! Set up local to global index table (np_local, xlocal, local)
    ! and global to local index table (global).
    vp%xlocal(1) = 1
    do r = 2, p
      vp%xlocal(r) = vp%xlocal(r-1)+vp%np_local(r-1)
    end do
    ir = 0
    do i = 1, np
      vp%local(vp%xlocal(m%part(i))+ir(m%part(i))) = i
      ir(m%part(i))                                = ir(m%part(i))+1
    end do

    ! Format of ghost:
    !
    ! np_ghost_neigh, np_ghost, xghost_neigh, xghost are components of vp, the vp% is
    ! ommited due to space constraints.
    !
    ! The following figure shows, how ghost points of node r are put into ghost:
    !
    !  |<--------------------------------np_ghost(r)---------------------------------->|
    !  |                                                                               |
    !  |<-np_ghost_neigh(r,1)->|     |<-np_ghost_neigh(r,p-1)->|<-np_ghost_neigh(r,p)->|
    !  |                       |     |                         |                       |
    ! -----------------------------------------------------------------------------------
    !  |                       | ... |                         |                       |
    ! -----------------------------------------------------------------------------------
    !  ^                             ^                         ^
    !  |                             |                         |
    !  xghost_neigh(r,1)             xghost_neigh(r,p-1)       xghost_neigh(r,p)
    !  |
    !  xghost(r)
   
    ! Mark and count ghost points and neighbours
    ! (set vp%np_ghost_neigh, vp%np_ghost, ghost_flag).
    vp%total          = 0
    ghost_flag        = 0
    vp%np_ghost_neigh = 0
    vp%np_ghost       = 0
    ! Check all nodes.
    do r = 1, p
      ! Check all points of this node.
      do i = vp%xlocal(r), vp%xlocal(r)+vp%np_local(r)-1
        ! Get coordinates of current point.
        p1 = m%Lxyz(vp%local(i), :)
        ! For all points in stencil.
        do j = 1, np_stencil
          ! Get index k of possible ghost point.
          p2 = p1+stencil(j, :)
          ! Check, whether p2 is in the box.
          ! FIXME: How about boundary conditions?
          ! If p2 is out of the box, ignore it for now.
          if(m%nr(1, 1).gt.p2(1).or.p2(1).gt.m%nr(2, 1).or. &
             m%nr(1, 2).gt.p2(2).or.p2(2).gt.m%nr(2, 2).or. &
             m%nr(1, 3).gt.p2(3).or.p2(3).gt.m%nr(2, 3)) cycle
          ! If it is in the box, get its point number.
          k = m%Lxyz_inv(p2(1), p2(2), p2(3))
          ! If p2 is in the box but does not belong to the geometry,
          ! its point number k should be greater than np (correct me
          ! if I am wrong). If this is the case, ignore it for now.
          if(k.gt.m%np) cycle
          ! At this point, it is sure that point number k is a
          ! relevant point.
          ! If this index k does not belong to partition of node r,
          ! then k is a ghost point for r with m%part(k) now being
          ! a neighbour of r.
          if(m%part(k).ne.r) then
            ! Only mark and count this ghost point, if it is not
            ! done yet. Otherwise, points would possibly be registered
            ! more than once.
            if(ghost_flag(k, r).eq.0) then
              ! Mark point i as ghost point for r from m%part(k).
              ghost_flag(k, r)                = m%part(k)
              ! Increase number of ghost points of r from m%part(k).
              vp%np_ghost_neigh(r, m%part(k)) = vp%np_ghost_neigh(r, m%part(k))+1
              ! Increase total number of ghostpoints of r.
              vp%np_ghost(r)                  = vp%np_ghost(r)+1
              ! One more ghost point.
              vp%total                        = vp%total+1
            end if
          end if
        end do
      end do
    end do

    ! Set index tables xghost and xghost_neigh.
    vp%xghost(1) = 1
    do r = 2, p
      vp%xghost(r) = vp%xghost(r-1)+vp%np_ghost(r-1)
    end do
    do r = 1, p
      vp%xghost_neigh(r, 1) = vp%xghost(r)
      do j = 2, p
        vp%xghost_neigh(r, j) = vp%xghost_neigh(r, j-1)    &
                                +vp%np_ghost_neigh(r, j-1)
      end do
    end do
    
    ! Get space for ghost point vector.
    allocate(vp%ghost(vp%total))

    ! Fill ghost as described above.
    irr = 0
    do i = 1, np
      do r = 1, p
        j = ghost_flag(i, r)
        ! If point i is a ghost point for r from j, save this
        ! information.
        if(j.ne.0) then
          vp%ghost(vp%xghost_neigh(r, j)+irr(r, j)) = i
          irr(r, j)                                 = irr(r, j)+1
        end if
      end do
    end do

    ! Write information about amount of ghost points.
    message(1) = 'Info: Total number of ghostpoints of each node:'
    write(message(2), '(a,100i7)') 'Info:', vp%np_ghost
    call write_info(2)

#ifdef DEBUG
    ! Write numbers and coordinates of each nodes ghost points
    ! to a single file (like in mesh_partition_init) called
    ! debug/mesh_partition/ghost_points.###.
    if(mpiv%node.eq.0) then
      call io_mkdir('debug/mesh_partition')
      do r = 1, p 
        write(filenum, '(i3.3)') r
        iunit = io_open('debug/mesh_partition/ghost_points.'//filenum, &
                        action='write')
        do i = 1, vp%np_ghost(r)
          j = vp%ghost(vp%xghost(r)+i-1)
          write(iunit, '(i8,3f12.8)') j, m%x(j, :)
        end do
      call io_close(iunit)
      end do
    end if
#endif

    ! Create reverse (global to local) lookup.
    ! Given a global point number i and a vector v_local of
    ! length vp%np_local(r)+vp%np_ghost(r) global(i, r) gives
    ! the index of point i in v_local as long, as this point is
    ! local to r or a ghost point for r (if global(i, r) > vp%np_local(r)
    ! it is a ghost point). If global(i, r) is 0 then i is neither local
    ! to r nor a ghost point of r. This indicates a serious error.
    vp%global = 0
    do r = 1, p
      ! Local points.
      do i = 1, vp%np_local(r)
        vp%global(vp%local(vp%xlocal(r)+i-1), r) = i
      end do
      ! Ghost points.
      do i = 1, vp%np_ghost(r)
        vp%global(vp%ghost(vp%xghost(r)+i-1), r) = vp%np_local(r)+i
      end do
    end do

    ! Complete entries in vp.
    vp%comm = comm
    vp%root = root
    vp%np   = m%np
    vp%p    = p

    ! Return shortcuts.
    np_local = vp%np_local(rank+1)
    np_ghost = vp%np_ghost(rank+1)
    call pop_sub()
    
  end subroutine vec_init


  ! Deallocate memory used by vp.
  subroutine vec_end(vp)
    type(pv_type), intent(inout) :: vp

    call push_sub('par_vec.vec_end')
    
    if(associated(vp%np_local)) then
      deallocate(vp%np_local)
      nullify(vp%np_local)
    endif
    if(associated(vp%xlocal)) then
      deallocate(vp%xlocal)
      nullify(vp%xlocal)
    endif
    if(associated(vp%local)) then
      deallocate(vp%local)
      nullify(vp%local)
    endif
    if(associated(vp%global)) then
      deallocate(vp%global)
      nullify(vp%global)
    endif
    if(associated(vp%np_ghost)) then
      deallocate(vp%np_ghost)
      nullify(vp%np_ghost)
    endif
    if(associated(vp%np_ghost_neigh)) then
      deallocate(vp%np_ghost_neigh)
      nullify(vp%np_ghost_neigh)
    endif
    if(associated(vp%xghost)) then
      deallocate(vp%xghost)
      nullify(vp%xghost)
    endif
    if(associated(vp%xghost_neigh)) then
      deallocate(vp%xghost_neigh)
      nullify(vp%xghost_neigh)
    endif
    if(associated(vp%ghost)) then
      deallocate(vp%ghost)
      nullify(vp%ghost)
    endif

    call pop_sub()

  end subroutine vec_end


  ! Scatters a vector v to all nodes in vp with respect to 
  ! to point -> node mapping in vp.
  ! v_local has at least to be of size vp%np_local(mpiv%node).
  subroutine vec_scatter(vp, v, v_local)
    type(pv_type), intent(in)  :: vp
    FLOAT,         intent(in)  :: v(:)
    FLOAT,         intent(out) :: v_local(:)

    integer              :: i         ! Counter.
    integer              :: ierr      ! MPI errorcode.
    integer              :: rank      ! Rank of node.
    integer, allocatable :: displs(:) ! Displacements for scatter.
    FLOAT,   allocatable :: v_tmp(:)  ! Send buffer.
    
    call push_sub('par_vec.vec_scatter')

    call MPI_Comm_rank(vp%comm, rank, ierr)
    call mpierr(ierr)
   
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
    call MPI_Scatterv(v_tmp, vp%np_local, displs, MPI_FLOAT, v_local, &
                      vp%np_local(rank+1), MPI_FLOAT,                 &
                      vp%root, vp%comm, ierr)
    call mpierr(ierr)

    if(mpiv%node.eq.vp%root) deallocate(v_tmp)

    deallocate(displs)

    call pop_sub()

  end subroutine vec_scatter


  ! Reverse operation of vec_scatter.
  ! All v_locals from the nodes are packed together
  ! into v on node vp%root in correct order.
  subroutine vec_gather(vp, v, v_local)
    type(pv_type), intent(in)  :: vp
    FLOAT,         intent(out) :: v(:)
    FLOAT,         intent(in)  :: v_local(:)

    integer              :: i         ! Counter.
    integer              :: ierr      ! MPI errorcode.
    integer              :: rank      ! Rank of node.
    integer, allocatable :: displs(:) ! Displacements for scatter.
    FLOAT,   allocatable :: v_tmp(:)  ! Receive buffer.
    
    call push_sub('par_vec.vec_gather')
    
    call MPI_Comm_rank(vp%comm, rank, ierr)
    call mpierr(ierr)

    ! Unfortunately, vp%xlocal ist not quite the required
    ! displacement vector.
    allocate(displs(vp%p))
    displs = vp%xlocal-1

    if(rank.eq.vp%root) allocate(v_tmp(vp%np))

    call MPI_Gatherv(v_local, vp%np_local(rank+1), MPI_FLOAT, v_tmp, &
                     vp%np_local, vp%xlocal, MPI_FLOAT,              &
                     vp%root, vp%comm, ierr)
    call mpierr(ierr)

    ! Copy values from v_tmp to their original position in v.
    if(rank.eq.vp%root) then
      do i = 1, vp%np
        v(vp%local(i)) = v_tmp(i)
      end do
      
      deallocate(v_tmp)
    end if

    call pop_sub()

  end subroutine vec_gather


  ! Updates ghost points of every node. A vector suitable
  ! for non local operations contains local values and
  ! ghost point values.
  ! Length of v_local must be
  ! vp%np_local(mpiv%node+1)+vp%np_ghost(mpiv%node+1)
  subroutine vec_ghost_update(vp, v_local)
    type(pv_type), intent(in)    :: vp
    FLOAT,         intent(inout) :: v_local(:)
    
    integer              :: i, j, k, r             ! Counters.
    integer              :: ierr                   ! MPI errorcode.
    integer              :: rank                   ! Rank of current node.
    integer              :: total                  ! Total number of ghost
                                                   ! points to send away.
    integer, allocatable :: sdispls(:), rdispls(:) ! Displacements for
                                                   ! MPI_Alltoallv.
    FLOAT,   allocatable :: ghost_send(:)          ! Send buffer.

    call push_sub('par_vec.vec_ghost_update')

    call MPI_Comm_rank(vp%comm, rank, ierr)
    call mpierr(ierr)

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
    call MPI_Alltoallv(ghost_send, vp%np_ghost_neigh(:, rank+1), sdispls, &
                       MPI_FLOAT, v_local(vp%np_local(rank+1)+1:),        &
                       vp%np_ghost_neigh(rank+1, :), rdispls, MPI_FLOAT,  &
                       vp%comm, ierr)
    call mpierr(ierr)
                       
    deallocate(sdispls, rdispls)
    deallocate(ghost_send)
    
    call pop_sub()
    
  end subroutine vec_ghost_update


  ! Sums over all elements of a vector v which was
  ! scattered to several v_local with vec_scatter.
  ! It is of no relevance, if v_local contains
  ! ghost points or not, because only the first
  ! vp%np_local(rank+1) elements are considered.
  FLOAT function vec_integrate(vp, v_local) result(s)
    type(pv_type), intent(in) :: vp
    FLOAT,         intent(in) :: v_local(:)

    integer :: rank
    integer :: ierr
    FLOAT   :: s_local ! Sum of v_local(i).

    call push_sub('par_vec.vec_integrate')
    
    call MPI_Comm_rank(vp%comm, rank, ierr)
    call mpierr(ierr)

    s_local = sum(v_local(:vp%np_local(rank+1)))

    call MPI_Allreduce(s_local, s, 1, MPI_FLOAT, MPI_SUM, vp%comm, ierr)
    call mpierr(ierr)

    call pop_sub()

  end function vec_integrate

  
  ! Apply parallel non local operator (fo = op fi).
  ! fi, fo have to be of length vp%np_local(rank+1)+vp%np_ghost(rank+1)
  ! and vec_ghost_update should have been called on the input vector
  ! in advance to have sane values at the ghost points.
  subroutine vec_dnl_operator(vp, op, fi, fo)
    type(pv_type),          intent(in)  :: vp
    type(nl_operator_type), intent(in)  :: op
    FLOAT,                  intent(in)  :: fi(:)
    FLOAT,                  intent(out) :: fo(:)

    integer              :: i, k, n
    integer              :: rank
    integer              :: ierr
    FLOAT, allocatable   :: w_re(:)
    !integer, allocatable :: i_local(:, :)

    call push_sub('par_vec.vec_dnl_operator')

    call MPI_Comm_rank(vp%comm, rank, ierr)
    call mpierr(ierr)
    
    n = op%n

    ! Maybe it is faster to create a local index first.
    !allocate(i_local(n, vp%np))
    !do i = 1, vp%np
    !  i_local(1:n, i) = vp%global(op%i(1:n, i), rank+1)
    !end do
    
    if(op%const_w) then
       allocate(w_re(n))
       w_re(1:n) = op%w_re(1:n, 1)
       ! Only operate on local points.
       do i = 1, vp%np_local(rank+1)
          ! Get local number of global point to lookup up
          ! stencil points in op%i. Those indices are then
          ! translated into local indices (i locally corresponds to
          ! k globally).
          k = vp%local(vp%xlocal(rank+1)+i-1)
          fo(i) = sum(w_re(1:n)*fi(vp%global(op%i(1:n, k), rank+1)))
          ! With local index:
          !fo(i) = sum(w_re(1:n)*fi(i_local(1:n, k))
       end do
       deallocate(w_re)
    else
       do i = 1, vp%np_local(rank+1)
          k = vp%local(vp%xlocal(rank+1)+i-1)
          fo(i) = sum(op%w_re(1:n,k)*fi(vp%global(op%i(1:n, k), rank+1)))
          !fo(i) = sum(op%w_re(1:n,k)*fi(i_local(1:n, k)))
       end do
    end if

    !deallocate(i_local)

    call pop_sub()

  end subroutine vec_dnl_operator

  
  ! (Very) basic error handling for MPI calls.
  ! Checks, if ierr indicates an MPI error and terminates
  ! cleanly if so.
  subroutine mpierr(ierr)
    integer :: ierr

    if(ierr.ne.MPI_SUCCESS) then
      write(message(1), '(a,i5)') 'MPI error, errorcode', ierr
      call write_fatal(1)
    end if
  end subroutine mpierr

end module par_vec
