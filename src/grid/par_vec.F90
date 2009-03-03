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

#include "global.h"

module par_vec_m

  ! Some general things and nomenclature:
  !
  ! - Points that are stored only on one node are
  !   called local points.
  ! - Local points, that are stored redundantly on
  !   another node because of the partitioning are
  !   called ghost points.
  ! - Points from the enlargement are only stored
  !   once on the corresponding node and are called
  !   boundary points.
  ! - np ist the total number of inner points.
  !
  ! A globally defined vector v has two parts:
  ! - v(1:np) are the inner points
  ! - v(np+1:np_part) are the boundary points
  ! In the typical case of zero boundary conditions
  ! v(np+1:np_part) is 0.
  ! The two parts are split according to the partitions.
  ! The result of this split are local vectors vl on each node
  ! which consist of three parts:
  ! - vl(1:np_local)                                     local points.
  ! - vl(np_local+1:np_local+np_ghost)                   ghost points.
  ! - vl(np_local+np_ghost+1:np_local+np_ghost+np_bndry) boundary points.
  !
  !
  ! Usage example for par_vec routines.
  !
  ! ! Initialize parallelization with mesh m and operator op
  ! ! initialized and given.
  ! ! m          = sys%gr%mesh
  ! ! stencil    = op%stencil
  !
  ! FLOAT              :: s
  ! FLOAT              :: u(np_global), v(np_global)
  ! FLOAT, allocatable :: ul(:), vl(:), wl(:)
  ! type(mesh_t)    :: m
  !
  ! ! Fill u, v with sensible values.
  ! ! ...
  !
  ! ! Allocate space for local vectors.
  ! allocate(ul(np_part))
  ! allocate(vl(np_part))
  ! allocate(wl(np_part))
  !
  ! ! Distribute vectors.
  ! call X(vec_scatter)(vp, u, ul)
  ! call X(vec_scatter)(vp, v, vl)
  !
  ! ! Compute some operator op: vl = op ul
  ! call X(vec_ghost_update)(vp, ul)
  ! call X(nl_operator_operate)(op, ul, vl)
  ! ! Gather result of op in one vector v.
  ! call X(vec_gather)(vp, v, vl)
  !
  ! ! Clean up.
  ! deallocate(ul, vl, wl)

  use c_pointer_m
  use global_m
  use iihash_m
  use io_m
  use math_m
  use index_m
  use messages_m
  use mpi_debug_m
  use mpi_m
  use profiling_m
  use stencil_m

  implicit none

  private
  public :: pv_t

  type pv_t
    ! The content of these members is node dependent.
    integer          :: rank                 ! Our rank in the communicator.
    integer          :: partno               ! Partition number of the
                                             ! current node
    integer, pointer :: isend_type(:)        ! The datatypes to send
    integer, pointer :: dsend_type(:)        ! ghost points
    integer, pointer :: zsend_type(:)
    integer, pointer :: rdispls(:)
    integer, pointer :: sdispls(:)
    integer, pointer :: rcounts(:)

    ! The following members are set independent of the nodes.
    integer                 :: npart                    ! Number of partitions.
    integer                 :: root                 ! The master node.
    integer                 :: comm                 ! MPI communicator to use.
    integer                 :: np                   ! Number of points in mesh.
    integer                 :: np_enl               ! Number of points in enlargement.
    integer, pointer        :: part(:)              ! Point -> partition.
    integer, pointer        :: np_local(:)          ! How many points has partition r?
    integer, pointer        :: xlocal(:)            ! Points of partition r start at
                                                    ! xlocal(r) in local.
    integer, pointer        :: local(:)             ! Partition r has points
                                                    ! local(xlocal(r):
                                                    ! xlocal(r)+np_local(r)-1).
    integer, pointer        :: np_bndry(:)          ! Number of boundary points.
    integer, pointer        :: xbndry(:)            ! Index of bndry(:).
    integer, pointer        :: bndry(:)             ! Global numbers of boundary
                                                    ! points.
    type(iihash_t), pointer :: global(:)            ! global(r) contains the global ->
                                                    ! local mapping for partition r.
    integer                 :: total                ! Total number of ghost points.
    integer, pointer        :: np_ghost(:)          ! How many ghost points has
                                                    ! partition r?
    integer, pointer        :: np_ghost_neigh(:, :) ! Number of ghost points per
                                                    ! neighbour per partition.
    integer, pointer        :: xghost(:)            ! Like xlocal.
    integer, pointer        :: xghost_neigh(:, :)   ! Like xghost for neighbours.
    integer, pointer        :: ghost(:)             ! Global indices of all local
  end type pv_t

#if defined(HAVE_MPI)

  integer :: SEND = 1, RECV = 2

  integer, public, parameter :: BLOCKING = 1, NON_BLOCKING = 2, NON_BLOCKING_COLLECTIVE = 3

  type pv_handle_t
    private
    integer          :: comm_method
    type(c_ptr)      :: nbc_h
    integer          :: nnb
    integer, pointer :: requests(:)
    integer, pointer :: status(:, :)
    integer, pointer :: ighost_send(:)
    FLOAT,   pointer :: dghost_send(:)
    CMPLX,   pointer :: zghost_send(:)
  end type pv_handle_t

  type(profile_t), save :: C_PROFILING_GHOST_UPDATE

  public ::              &
    pv_handle_t,         &
    pv_handle_init,      &
    pv_handle_test,      &
    pv_handle_wait,      &
    pv_handle_end

  public ::              &
    vec_init,            &
    vec_end,             &
    vec_global2local,    &
    dvec_scatter,        &
    zvec_scatter,        &
    ivec_scatter,        &
    dvec_scatter_bndry,  &
    zvec_scatter_bndry,  &
    ivec_scatter_bndry,  &
    dvec_scatter_all,    &
    zvec_scatter_all,    &
    ivec_scatter_all,    &
    dvec_gather,         &
    zvec_gather,         &
    ivec_gather,         &
    dvec_allgather,      &
    zvec_allgather,      &
    ivec_allgather,      &
    dvec_ghost_update,   &
    zvec_ghost_update,   &
    ivec_ghost_update,   &
    dvec_ighost_update,  &
    zvec_ighost_update,  &
    ivec_ighost_update

contains

  ! Initializes a pv_type object (parallel vector).
  ! It computes the local to global and global to local index tables
  ! and the ghost point exchange.
  !
  ! Note: we can not pass in the i(:, :) array from the stencil
  ! because it is not yet computed (it is local to a node and
  ! must be initialized some time after vec_init is run).
  ! Warning: The naming scheme for the np_ variables if different
  ! from how it is in the rest of the code (for historical reasons
  ! and also because the vec_init has more a global than a local point
  ! of view on the mesh): See the comments in the parameter list.
  subroutine vec_init(comm, root, part, np, np_part, idx, stencil, dim, vp)
    integer,         intent(in)  :: comm         ! Communicator to use.
    integer,         intent(in)  :: root         ! The master node.

    ! The next seven entries come from the mesh.
    integer,         intent(in)  :: part(:)      ! Point -> partition.
    integer,         intent(in)  :: np           ! m%np_global
    integer,         intent(in)  :: np_part      ! m%np_part_global
    type(index_t),   intent(in)  :: idx
    type(stencil_t), intent(in)  :: stencil      ! The stencil for which to calculate ghost points.
    integer,         intent(in)  :: dim          ! Number of dimensions.
    type(pv_t),      intent(out) :: vp           ! Description of partition.

    ! Careful: MPI counts node ranks from 0 to numproc-1.
    ! Partition numbers from METIS range from 1 to numproc.
    ! For this reason, all ranks are incremented by one.
    integer                     :: p                ! Number of partitions.
    integer                     :: np_enl           ! Number of points in enlargement.
    integer                     :: i, j, k, r       ! Counters.
    integer, allocatable        :: ir(:), irr(:, :) ! Counters.
    integer                     :: rank             ! Rank of current node.
    integer                     :: p1(MAX_DIM)      ! Points.
    type(iihash_t), allocatable :: ghost_flag(:)    ! To remember ghost pnts.
    integer                     :: iunit            ! For debug output to files.
    character(len=3)            :: filenum
    integer                     :: tmp
    logical                     :: found

    call push_sub('par_vec.vec_init')

    ! Shortcuts.
    call MPI_Comm_Size(comm, p, mpi_err)
    np_enl = np_part-np

    ! Store partition number and rank for later reference.
    ! Having both variables is a bit redundant but makes the code readable.
    call MPI_Comm_Rank(comm, rank, mpi_err)
    vp%rank   = rank
    vp%partno = rank + 1

    ALLOCATE(ghost_flag(p),            p)
    ALLOCATE(ir(p),                    p)
    ALLOCATE(irr(p, p),                p*p)
    ALLOCATE(vp%part(np+np_enl),       np+np_enl)
    ALLOCATE(vp%np_local(p),           p)
    ALLOCATE(vp%xlocal(p),             p)
    ALLOCATE(vp%local(np),             np)
    ALLOCATE(vp%np_bndry(p),           p)
    ALLOCATE(vp%xbndry(p),             p)
    ALLOCATE(vp%bndry(np_enl),         np_enl)
    ALLOCATE(vp%global(p),             p)
    ALLOCATE(vp%np_ghost(p),           p)
    ALLOCATE(vp%np_ghost_neigh(p, p),  p*p)
    ALLOCATE(vp%xghost(p),             p)
    ALLOCATE(vp%xghost_neigh(p, p),    p*p)

    ! Count number of points for each node.
    ! Local points.
    vp%np_local = 0
    do i = 1, np
      vp%np_local(part(i)) = vp%np_local(part(i))+1
    end do
    ! Boundary points.
    vp%np_bndry = 0
    do i = 1, np_enl
      vp%np_bndry(part(i+np)) = vp%np_bndry(part(i+np))+1
    end do

    ! Set up local to global index table for local points
    ! (xlocal, local) and for boundary points (xbndry, bndry).
    vp%xlocal(1) = 1
    vp%xbndry(1) = 1
    do r = 2, p
      vp%xlocal(r) = vp%xlocal(r-1)+vp%np_local(r-1)
      vp%xbndry(r) = vp%xbndry(r-1)+vp%np_bndry(r-1)
    end do
    ir = 0
    do i = 1, np
      vp%local(vp%xlocal(part(i))+ir(part(i))) = i
      ir(part(i))                              = ir(part(i))+1
    end do
    ir = 0
    do i = np+1, np+np_enl
      vp%bndry(vp%xbndry(part(i))+ir(part(i))) = i
      ir(part(i))                              = ir(part(i))+1
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
    do r = 1, p
      call iihash_init(ghost_flag(r), vp%np_local(r))
    end do
    vp%total          = 0
    vp%np_ghost_neigh = 0
    vp%np_ghost       = 0
    ! Check all nodes.
    do r = 1, p
      ! Check all points of this node.
      do i = vp%xlocal(r), vp%xlocal(r)+vp%np_local(r)-1
        ! Get coordinates of current point.
        call index_to_coords(idx, dim, vp%local(i), p1)

        ! For all points in stencil.
        do j = 1, stencil%size
          ! Get point number of possible ghost point.
          k = index_from_coords(idx, dim, p1(:) + stencil%points(:, j))
          ASSERT(k.ne.0)
          ! If this index k does not belong to partition of node r,
          ! then k is a ghost point for r with part(k) now being
          ! a neighbour of r.
          if(part(k).ne.r) then
            ! Only mark and count this ghost point, if it is not
            ! done yet. Otherwise, points would possibly be registered
            ! more than once.
            tmp = iihash_lookup(ghost_flag(r), k, found)
            if(.not.found) then
              ! Mark point i as ghost point for r from part(k).
              call iihash_insert(ghost_flag(r), k, part(k))
              ! Increase number of ghost points of r from part(k).
              vp%np_ghost_neigh(r, part(k)) = vp%np_ghost_neigh(r, part(k))+1
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
    ALLOCATE(vp%ghost(vp%total), vp%total)

    ! Fill ghost as described above.
    irr = 0
    do i = 1, np+np_enl
      do r = 1, p
        j = iihash_lookup(ghost_flag(r), i, found)
        ! If point i is a ghost point for r from j, save this
        ! information.
        if(found) then
          vp%ghost(vp%xghost_neigh(r, j)+irr(r, j)) = i
          irr(r, j)                                 = irr(r, j)+1
        end if
      end do
    end do

    if(in_debug_mode) then
      ! Write numbers and coordinates of each nodes ghost points
      ! to a single file (like in mesh_partition_init) called
      ! debug/mesh_partition/ghost_points.###.
      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir('debug/mesh_partition')
        do r = 1, p
          write(filenum, '(i3.3)') r
          iunit = io_open('debug/mesh_partition/ghost_points.'//filenum, &
            action='write')
          do i = 1, vp%np_ghost(r)
            j = vp%ghost(vp%xghost(r)+i-1)
            write(iunit, '(4i8)') j, idx%Lxyz(j, :)
          end do
          call io_close(iunit)
        end do
      end if
    end if

    ! Set up the global to local point number mapping.
    do r = 1, p
      ! Create hash table.
      call iihash_init(vp%global(r), vp%np_local(r)+vp%np_ghost(r)+vp%np_bndry(r))
      ! Insert local points.
      do i = 1, vp%np_local(r)
        call iihash_insert(vp%global(r), vp%local(vp%xlocal(r)+i-1), i)
      end do
      ! Insert ghost points.
      do i = 1, vp%np_ghost(r)
        call iihash_insert(vp%global(r), vp%ghost(vp%xghost(r)+i-1), i+vp%np_local(r))
      end do
      ! Insert boundary points.
      do i = 1, vp%np_bndry(r)
        call iihash_insert(vp%global(r), vp%bndry(vp%xbndry(r)+i-1), i+vp%np_local(r)+vp%np_ghost(r))
      end do
    end do
    
    ! Complete entries in vp.
    vp%comm   = comm
    vp%root   = root
    vp%np     = np
    vp%np_enl = np_enl
    vp%npart      = p
    vp%part   = part

    call init_mpi_datatypes

    do r = 1, p
      call iihash_end(ghost_flag(r))
    end do

    call pop_sub()

  contains
    
    subroutine init_mpi_datatypes
      integer, allocatable :: blocklengths(:), displacements(:), offsets(:)
      integer :: ii, kk, ipart, total, ierr, nblocks

      ALLOCATE(vp%isend_type(1:vp%npart), vp%npart)
      ALLOCATE(vp%dsend_type(1:vp%npart), vp%npart)
      ALLOCATE(vp%zsend_type(1:vp%npart), vp%npart)

      ! Iterate over all possible receivers.
      do ipart = 1, vp%npart
        total = vp%np_ghost_neigh(ipart, vp%partno)

        if(total == 0) cycle

        ALLOCATE(blocklengths(total), total)
        ALLOCATE(offsets(total), total)
        ALLOCATE(displacements(total), total)
        
        ! Collect all local points that have to be sent to neighbours.
        
        ! Iterate over all ghost points that ipart wants.
        do ii = 0, vp%np_ghost_neigh(ipart, vp%partno) - 1
          ! Get global number kk of i-th ghost point.
          kk = vp%ghost(vp%xghost_neigh(ipart, vp%partno) + ii)
          ! Lookup up local number of point kk
          displacements(ii + 1) = vec_global2local(vp, kk, vp%partno) - 1
        end do

        call get_blocks(total, displacements, nblocks, blocklengths, offsets)

        call MPI_Type_indexed(nblocks, blocklengths, offsets, MPI_INTEGER, vp%isend_type(ipart), ierr)
        call MPI_Type_indexed(nblocks, blocklengths, offsets, MPI_FLOAT,   vp%dsend_type(ipart), ierr)
        call MPI_Type_indexed(nblocks, blocklengths, offsets, MPI_CMPLX,   vp%zsend_type(ipart), ierr)
        
        call MPI_Type_commit(vp%isend_type(ipart), ierr)
        call MPI_Type_commit(vp%dsend_type(ipart), ierr)
        call MPI_Type_commit(vp%zsend_type(ipart), ierr)

        deallocate(blocklengths, displacements, offsets)
        
      end do

      ! Send and receive displacements.
      ! Send displacement cannot directly be calculated
      ! from vp%xghost_neigh because those are indices for
      ! vp%np_ghost_neigh(vp%partno, :) and not
      ! vp%np_ghost_neigh(:, vp%partno) (rank being fixed).
      ! So what gets done is to pick out the number of ghost points
      ! each partition r wants to have from the current partiton
      ! vp%partno.
      
      ALLOCATE(vp%sdispls(1:vp%npart), vp%npart)
      ALLOCATE(vp%rdispls(1:vp%npart), vp%npart)
      ALLOCATE(vp%rcounts(1:vp%npart), vp%npart)

      vp%sdispls(1) = 0
      do ipart = 2, vp%npart
        vp%sdispls(ipart) = vp%sdispls(ipart - 1) + vp%np_ghost_neigh(ipart - 1, vp%partno)
      end do

      ! This is like in vec_scatter/gather.
      vp%rdispls(1:vp%npart) = vp%xghost_neigh(vp%partno, 1:vp%npart) - vp%xghost(vp%partno)
      
      vp%rcounts(1:vp%npart) = vp%np_ghost_neigh(vp%partno, 1:vp%npart)

    end subroutine init_mpi_datatypes

  end subroutine vec_init


  ! Deallocate memory used by vp.
  subroutine vec_end(vp)
    type(pv_t), intent(inout) :: vp

    integer :: ipart, r

    call push_sub('par_vec.vec_end')

    deallocate(vp%rdispls)
    deallocate(vp%sdispls)
    deallocate(vp%rcounts)

    if(associated(vp%isend_type)) then

      do ipart = 1, vp%npart
        if(vp%np_ghost_neigh(ipart, vp%partno) == 0) cycle
        call MPI_Type_free(vp%isend_type(ipart), mpi_err)
        call MPI_Type_free(vp%dsend_type(ipart), mpi_err)
        call MPI_Type_free(vp%zsend_type(ipart), mpi_err)
      end do
      deallocate(vp%isend_type)
      deallocate(vp%dsend_type)
      deallocate(vp%zsend_type)
      nullify(vp%isend_type)
      nullify(vp%dsend_type)
      nullify(vp%zsend_type)
    end if

    if(associated(vp%part)) then
      deallocate(vp%part)
      nullify(vp%part)
    end if
    if(associated(vp%np_local)) then
      deallocate(vp%np_local)
      nullify(vp%np_local)
    end if
    if(associated(vp%xlocal)) then
      deallocate(vp%xlocal)
      nullify(vp%xlocal)
    end if
    if(associated(vp%local)) then
      deallocate(vp%local)
      nullify(vp%local)
    end if
    if(associated(vp%np_bndry)) then
      deallocate(vp%np_bndry)
      nullify(vp%np_bndry)
    end if
    if(associated(vp%xbndry)) then
      deallocate(vp%xbndry)
      nullify(vp%xbndry)
    end if
    if(associated(vp%bndry)) then
      deallocate(vp%bndry)
      nullify(vp%bndry)
    end if
    if(associated(vp%np_ghost)) then
      deallocate(vp%np_ghost)
      nullify(vp%np_ghost)
    end if
    if(associated(vp%np_ghost_neigh)) then
      deallocate(vp%np_ghost_neigh)
      nullify(vp%np_ghost_neigh)
    end if
    if(associated(vp%xghost)) then
      deallocate(vp%xghost)
      nullify(vp%xghost)
    end if
    if(associated(vp%xghost_neigh)) then
      deallocate(vp%xghost_neigh)
      nullify(vp%xghost_neigh)
    end if
    if(associated(vp%ghost)) then
      deallocate(vp%ghost)
      nullify(vp%ghost)
    end if
    if(associated(vp%global)) then
      do r = 1, vp%npart
        call iihash_end(vp%global(r))
      end do
      deallocate(vp%global)
      nullify(vp%global)
    end if

    call pop_sub()

  end subroutine vec_end

  subroutine pv_handle_init(this, vp, comm_method)
    type(pv_handle_t), intent(out) :: this
    type(pv_t),        intent(in)  :: vp
    integer,           intent(in)  :: comm_method

    this%comm_method = comm_method

    select case(this%comm_method)
#ifdef HAVE_LIBNBC
    case(NON_BLOCKING_COLLECTIVE)
      call NBCF_Newhandle(this%nbc_h)
#endif
    case(NON_BLOCKING)
      ALLOCATE(this%requests(1:vp%npart*2), vp%npart*2)
      ALLOCATE(this%status(MPI_STATUS_SIZE, 1:vp%npart*2), vp%npart*2)
    end select
    nullify(this%ighost_send, this%dghost_send, this%zghost_send)
  end subroutine pv_handle_init

  subroutine pv_handle_end(this)
    type(pv_handle_t), intent(inout) :: this

    select case(this%comm_method)
#ifdef HAVE_LIBNBC
    case(NON_BLOCKING_COLLECTIVE)
      call NBCF_Freehandle(this%nbc_h)
#endif
    case(NON_BLOCKING)
      deallocate(this%requests, this%status)
    end select

  end subroutine pv_handle_end

  subroutine pv_handle_test(this)
    type(pv_handle_t), intent(inout) :: this
    
    select case(this%comm_method)
#ifdef HAVE_LIBNBC
    case(NON_BLOCKING_COLLECTIVE)
      call NBCF_Test(this%nbc_h, mpi_err)
#endif
    case(NON_BLOCKING)
    end select

  end subroutine pv_handle_test

  subroutine pv_handle_wait(this)
    type(pv_handle_t), intent(inout) :: this

    type(profile_t), save :: prof
    
    call profiling_in(prof, "GHOST_UPDATE_WAIT")

    select case(this%comm_method)
#ifdef HAVE_LIBNBC
    case(NON_BLOCKING_COLLECTIVE)
      call NBCF_Wait(this%nbc_h, mpi_err)
#endif
    case(NON_BLOCKING)
      call MPI_Waitall(this%nnb, this%requests, this%status, mpi_err)
    end select

    if(associated(this%ighost_send)) then
      deallocate(this%ighost_send)
      nullify(this%ighost_send)
    end if
    if(associated(this%dghost_send)) then
      deallocate(this%dghost_send)
      nullify(this%dghost_send)
    end if
    if(associated(this%zghost_send)) then
      deallocate(this%zghost_send)
      nullify(this%zghost_send)
    end if
    call profiling_out(prof)
  end subroutine pv_handle_wait


  ! ---------------------------------------------------------
  ! Returns local local number of global point i on partition r.
  ! If the result is zero, the point is neither a local nor a ghost
  ! point on r.
  integer function vec_global2local(vp, i, r)
    type(pv_t), intent(in) :: vp
    integer,    intent(in) :: i
    integer,    intent(in) :: r

    integer :: n
    logical :: found

    vec_global2local = 0
    if(associated(vp%global)) then
      n = iihash_lookup(vp%global(r), i, found)
      if(found) vec_global2local = n
    end if

  end function vec_global2local

#include "undef.F90"
#include "complex.F90"
#include "par_vec_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "par_vec_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "par_vec_inc.F90"

#endif
end module par_vec_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
