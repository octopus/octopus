!! Copyright (C) 2005-2006 Florian Lorenzen, Heiko Appel, J. Alberdi-Rodriguez
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

#include "global.h"
 
  !> Some general things and nomenclature:
  !!
  !! - Points that are stored only on one process are
  !!   called local points.
  !! - Local points that are stored redundantly on
  !!   another process because of the partitioning are
  !!   called ghost points.
  !! - Points from the enlargement are only stored
  !!   once on the corresponding process and are called
  !!   boundary points.
  !! - np is the total number of inner points.
  !!
  !! A globally defined vector v has two parts:
  !! - v(1:np) are the inner points
  !! - v(np+1:np_part) are the boundary points
  !! In the typical case of zero boundary conditions
  !! v(np+1:np_part) is 0.
  !! The two parts are split according to the partitions.
  !! The result of this split are local vectors vl on each process
  !! which consist of three parts:
  !! - vl(1:np_local_vec)                                     local points.
  !! - vl(np_local_vec+1:np_local_vec+np_ghost)                   ghost points.
  !! - vl(np_local_vec+np_ghost+1:np_local_vec+np_ghost+np_bndry) boundary points.
  !!
  !!
  !! Usage example for par_vec routines.
  !! \verbatim
  !!
  !! ! Initialize parallelization with mesh and operator op
  !! ! initialized and given.
  !! ! mesh       = sys%gr%mesh
  !! ! stencil    = op%stencil
  !!
  !! FLOAT              :: uu(np_global), vv(np_global)
  !! FLOAT, allocatable :: ul(:), vl(:), wl(:)
  !! type(mesh_t)       :: mesh
  !!
  !! ! Fill uu, vv with sensible values.
  !! ! ...
  !!
  !! ! Allocate space for local vectors.
  !! allocate(ul(np_part))
  !! allocate(vl(np_part))
  !! allocate(wl(np_part))
  !!
  !! ! Distribute vectors.
  !! call X(vec_scatter)(vp, uu, ul)
  !! call X(vec_scatter)(vp, vv, vl)
  !!
  !! ! Compute some operator op: vl = op ul
  !! call X(vec_ghost_update)(vp, ul)
  !! call X(nl_operator_operate)(op, ul, vl)
  !! !! Gather result of op in one vector vv.
  !! call X(vec_gather)(vp, vv, vl)
  !!
  !! ! Clean up.
  !! deallocate(ul, vl, wl)
  !! \endverbatim
module par_vec_m
  use global_m
  use iihash_m
  use index_m
  use io_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use partition_m
  use profiling_m
  use stencil_m
  use subarray_m

  implicit none

  private

  public ::                        &
    pv_t !< parallel information

#if defined(HAVE_MPI)
  public ::                        &
    vec_init,                      &
    vec_end,                       &
    vec_global2local,              &
    vec_index2local,               &
    dvec_scatter,                  &
    zvec_scatter,                  &
    ivec_scatter,                  &
    dvec_gather,                   &
    zvec_gather,                   &
    ivec_gather,                   &
    dvec_allgather,                &
    zvec_allgather,                &
    ivec_allgather
#endif
  !> Parallel information
  type pv_t
    ! The content of these members is process-dependent.
    integer          :: rank                 !< Our rank in the communicator. 
    !> Partition number of the
    !! current process
    integer          :: partno              
    type(subarray_t) :: sendpoints
    integer, pointer :: sendpos(:)

    integer, pointer :: rdispls(:)
    integer, pointer :: sdispls(:)
    integer, pointer :: rcounts(:)

    ! The following members are set independent of the processs.
    integer                 :: npart                !< Number of partitions.
    integer                 :: root                 !< The master process.
    integer                 :: comm                 !< MPI communicator to use.
    integer                 :: np_global            !< Number of points in mesh.
    integer, pointer        :: part_vec(:)          !< Global point        -> partition.
    integer, pointer        :: part_local(:)        !< Local point         -> partition
    integer, pointer        :: part_local_rev(:)    !< Local point`s value -> partition

    integer, pointer        :: np_local_vec(:)      !< How many points has partition r?
                                                    !! Global vector; npart elements.
    integer                 :: np_local             !< How many points has running partition? 
                                                    !! Local value.
    integer, pointer        :: xlocal_vec(:)        !< Points of partition r start at
                                                    !! xlocal_vec(r) in local. Global start point
                                                    !! of the local index.  
                                                    !! Global vector; npart elements.
    integer                 :: xlocal               !< Starting index of running process in local(:) vector.
                                                    !! Local value.
          
    integer, pointer        :: local_vec(:)         !< Partition r has points
                                                    !! local_vec(xlocal_vec(r):
                                                    !! xlocal_vec(r)+np_local_vec(r)-1). 
                                                    !! Global vector; np_global elements    
    integer, pointer        :: local(:)             !< Local points of running process
                                                    !! Local vector; np_local elements
    integer, pointer        :: recv_count(:)        !< Number of points to receive from all the other processes
    integer, pointer        :: send_count(:)        !< Number of points to send to all the other processes
                                                    !! in a MPI_Alltoallv.
    integer, pointer        :: recv_disp(:)         !< Displacement of points to receive from all the other processes
    integer, pointer        :: send_disp(:)         !< Displacement of points to send to all the other processes
                                                    !! in a MPI_Alltoallv.
    integer                 :: np_bndry             !< Number of boundary points.
                                                    !! Local value
    integer                 :: xbndry               !< Starting index of running process in bndry(:) 
                                                    !! Local value
 
    integer, pointer        :: bndry(:)             !< Global numbers of boundary points.
                                                    !! Global vector; np_enl elements
      
    type(iihash_t), pointer :: global(:)            !< global(r) contains the global ->
                                                    !! local mapping for partition r.   
    integer                 :: total                !< Total number of ghost points. 

    integer                 :: np_ghost                 !< How many ghost points has partition r?
                                                        !! Local value
    integer, pointer        :: np_ghost_neigh_partno(:) !< Number of the neighbours ghost points of the actual process
    integer                 :: xghost                   !< Starting index of running procces in ghost(:) vector.
    integer, pointer        :: ghost(:)                 !< Global indices of all local points.
                                                        !! Global vector; vp%total elements
  end type pv_t

#if defined(HAVE_MPI)

  type(profile_t), save :: prof_scatter
  type(profile_t), save :: prof_allgather

contains

  !> Initializes a pv_type object (parallel vector).
  !! It computes the local-to-global and global-to-local index tables
  !! and the ghost point exchange.
  !!
  !! Note: we cannot pass in the i(:, :) array from the stencil
  !! because it is not yet computed (it is local to a process and
  !! must be initialized some time after vec_init is run).
  !! \warning The naming scheme for the np_ variables is different
  !! from how it is in the rest of the code (for historical reasons
  !! and also because the vec_init has more a global than local point
  !! of view on the mesh): See the comments in the parameter list.
  subroutine vec_init(comm, root, np_global, np_part_global, idx, stencil, dim, periodic_dim, &
       inner_partition, bndry_partition, vp)
    integer,         intent(in)  :: comm         !< Communicator to use.
    integer,         intent(in)  :: root         !< The master process.

    !> The next seven entries come from the mesh.
    integer,          intent(in)    :: np_global      !< mesh%np_global
    integer,          intent(in)    :: np_part_global !< mesh%np_part_global
    type(index_t),    intent(in)    :: idx
    type(stencil_t),  intent(in)    :: stencil        !< The stencil for which to calculate ghost points.
    integer,          intent(in)    :: dim            !< Number of dimensions.
    integer,          intent(in)    :: periodic_dim   !< Number of periodic dimensions
    type(partition_t),intent(in)    :: inner_partition
    type(partition_t),intent(in)    :: bndry_partition
    type(pv_t),       intent(inout) :: vp             !< Description of partition.

    ! Careful: MPI counts process ranks from 0 to numproc-1.
    ! Partition numbers from METIS range from 1 to numproc.
    ! For this reason, all ranks are incremented by one.
    integer                     :: npart            !< Number of partitions.
    integer                     :: np_enl           !< Number of points in enlargement.
    integer                     :: gip, ip, jp, kp, jj, index, inode, jnode !< Counters.
    integer, allocatable        :: irr(:)           !< Counter.
    integer                     :: p1(MAX_DIM)      !< Points.
    type(iihash_t), allocatable :: ghost_flag(:)    !< To remember ghost pnts.
    integer                     :: iunit            !< For debug output to files.
    character(len=3)            :: filenum
    integer                     :: tmp, init, size, ii
    logical                     :: found
   
    integer                     :: idir, ipart, np_inner, np_bndry
    integer, pointer            :: np_ghost_tmp(:), np_bndry_tmp(:)
    !> Number of ghost points per
    !! neighbour per partition.      
    integer, pointer            :: np_ghost_neigh(:, :) 
    integer, pointer            :: xbndry_tmp(:)    !< Starting index of process i in bndry(:). 
    integer, pointer            :: xghost_tmp(:)  
    integer, pointer            :: xghost_neigh_partno(:)   !< Like xghost for neighbours.
    integer, pointer            :: xghost_neigh_back(:)     !< Same as previous, but outward
    integer, pointer            :: points(:), points_bndry(:), part_bndry(:), part_inner(:)

    PUSH_SUB(vec_init)

    ! Shortcuts.
    call MPI_Comm_Size(comm, npart, mpi_err)
    np_enl = np_part_global - np_global

    ! Store partition number and rank for later reference.
    ! Having both variables is a bit redundant but makes the code readable.
    call MPI_Comm_Rank(comm, vp%rank, mpi_err)
    vp%partno = vp%rank + 1

    SAFE_ALLOCATE(ghost_flag(1:npart))
    SAFE_ALLOCATE(irr(1:npart))
    SAFE_ALLOCATE(vp%np_local_vec(1:npart))
    SAFE_ALLOCATE(vp%xlocal_vec(1:npart))
    SAFE_ALLOCATE(np_bndry_tmp(1:npart))
    SAFE_ALLOCATE(xbndry_tmp(1:npart))
    SAFE_ALLOCATE(vp%global(1:npart))
    SAFE_ALLOCATE(vp%np_ghost_neigh_partno(1:npart))

    ! Count number of points for each process.
    ! Local points.
    call partition_get_np_local(inner_partition, vp%np_local_vec)
    vp%np_local = vp%np_local_vec(vp%partno)

    ! Boundary points.
    call partition_get_np_local(bndry_partition, np_bndry_tmp)
    vp%np_bndry = np_bndry_tmp(vp%partno)

    ! Set up local-to-global index table for local points
    ! (xlocal_vec, local) and for boundary points (xbndry, bndry).
    vp%xlocal_vec(1) = 1
    xbndry_tmp(1) = 1
    ! Set the starting point of local and boundary points
    do inode = 2, npart
      vp%xlocal_vec(inode) = vp%xlocal_vec(inode - 1) + vp%np_local_vec(inode - 1)
      xbndry_tmp(inode) = xbndry_tmp(inode - 1) + np_bndry_tmp(inode - 1)
    end do
    vp%xlocal = vp%xlocal_vec(vp%partno)
    ! Set the local and boundary points
    call init_local

    ! Format of ghost:
    !
    ! np_ghost_neigh, np_ghost, xghost_neigh, xghost are components of vp.
    ! The vp% is omitted due to space constraints.
    !
    ! The following figure shows how ghost points of process "inode" are put into ghost:
    !
    !  |<-------------------------------------------np_ghost(inode)--------------------------------------->|
    !  |                                                                                                   |
    !  |<-np_ghost_neigh(inode,1)->|     |<-np_ghost_neigh(inode,npart-1)->|<-np_ghost_neigh(inode,npart)->|
    !  |                           |     |                                 |                               |
    ! -------------------------------------------------------------------------------------------------------
    !  |                           | ... |                                 |                               |
    ! -------------------------------------------------------------------------------------------------------
    !  ^                                 ^                                 ^
    !  |                                 |                                 |
    !  xghost_neigh(inode,1)             xghost_neigh(inode,npart-1)       xghost_neigh(inode,npart)
    !  |
    !  xghost(inode)

    ! Mark and count ghost points and neighbours
    ! (set vp%np_ghost_neigh, vp%np_ghost, ghost_flag).
    do inode = 1, npart
      call iihash_init(ghost_flag(inode), vp%np_local_vec(inode))
    end do

    do jj = 1, stencil%size
      ASSERT(all(stencil%points(1:dim, jj) <= idx%enlarge(1:dim)))
    end do
    
    SAFE_ALLOCATE(vp%send_count(1:npart))
    vp%send_count = 0

    vp%total                 = 0
    vp%np_ghost_neigh_partno = 0
    vp%np_ghost              = 0
    ip                       = 0
    inode = vp%partno

    SAFE_ALLOCATE(points(1:vp%np_local))
    SAFE_ALLOCATE(vp%part_local(1:vp%np_local))
    ! Check all points of this node and create the local partition matrix
    do gip = vp%xlocal, vp%xlocal + vp%np_local - 1
      ip = ip + 1
      points(ip) = gip
    end do
    call partition_get_partition_number(inner_partition,vp%np_local, &
         points, vp%part_local)
    SAFE_DEALLOCATE_P(points)

    ip       = 0
    np_inner = 0
    np_bndry = 0
    SAFE_ALLOCATE(points(1:vp%np_local*stencil%size))
    SAFE_ALLOCATE(points_bndry(1:vp%np_bndry*stencil%size))
    points_bndry = 0
    do gip = vp%xlocal, vp%xlocal + vp%np_local - 1
      ip = ip + 1
      ! Update the receiving point
      ipart = vp%part_local(ip)
   
      vp%send_count(ipart) = vp%send_count(ipart) + 1
      ! Get coordinates of current point.
      call index_to_coords(idx, dim, vp%local(gip), p1)
      
      ! For all points in stencil.
      do jj = 1, stencil%size
        ! Get point number of possible ghost point.
        index = index_from_coords(idx, dim, p1(:) + stencil%points(:, jj))
        ASSERT(index /= 0)
        ! Global index can be either in the mesh or in the boundary.
        ! Different treatment is needed for each case.
        if (index > np_global) then
          np_bndry = np_bndry + 1
          points_bndry(np_bndry) = index - np_global
        else
          np_inner = np_inner + 1
          points(np_inner) = index
        end if
      end do
    end do
    
    SAFE_ALLOCATE(part_inner(1:np_inner))
    SAFE_ALLOCATE(part_bndry(1:np_bndry))

    call partition_get_partition_number(inner_partition, np_inner, &
         points, part_inner)
    call partition_get_partition_number(bndry_partition, np_bndry, &
         points_bndry, part_bndry)

    vp%total = 0
    do ip = 1, np_inner
      ! If this index does not belong to partition of working node "inode",
      ! then index is a ghost point for "inode" with part(index) now being
      ! a neighbour of "inode".
      if ( part_inner(ip) /= inode ) then
        ! Only mark and count this ghost point, if it is not
        ! done yet. Otherwise, points would possibly be registered
        ! more than once.
        tmp = iihash_lookup(ghost_flag(inode), points(ip), found)
        if(.not.found) then
          ! Mark point ip as ghost point for inode from part(index).
          call iihash_insert(ghost_flag(inode), points(ip), part_inner(ip))
          ! Increase number of ghost points of inode from part(index).
          vp%np_ghost_neigh_partno(part_inner(ip)) = vp%np_ghost_neigh_partno(part_inner(ip))+1
          ! Increase total number of ghostpoints of inode.
          vp%np_ghost = vp%np_ghost + 1
          ! One more ghost point.
          vp%total = vp%total + 1
        end if
      end if
    end do
    
    ! The same for boundary points
    do ip = 1, np_bndry
      if ( part_bndry(ip) /= inode ) then
        tmp = iihash_lookup(ghost_flag(inode), &
             points_bndry(ip)+np_global, found)
        if(.not.found) then
          call iihash_insert(ghost_flag(inode), &
               points_bndry(ip)+np_global, part_bndry(ip))
          vp%np_ghost_neigh_partno(part_bndry(ip)) = vp%np_ghost_neigh_partno(part_bndry(ip))+1
          vp%np_ghost = vp%np_ghost + 1
          vp%total = vp%total + 1
        end if
      end if
    end do

    SAFE_DEALLOCATE_P(points)
    SAFE_DEALLOCATE_P(points_bndry)
    SAFE_DEALLOCATE_P(part_inner)
    SAFE_DEALLOCATE_P(part_bndry)

    call init_MPI_Alltoall
    tmp=0
    call MPI_Allreduce(vp%total, tmp, 1, MPI_INTEGER, MPI_SUM, comm, mpi_err)
    vp%total = tmp    
    SAFE_ALLOCATE(np_ghost_neigh(1:npart, 1:npart))
    ! Distribute local data to all processes
    call MPI_Allgather(vp%np_ghost_neigh_partno(1),npart,MPI_INTEGER, &
         np_ghost_neigh(1,1),npart,MPI_INTEGER, &
         comm, mpi_err)

    ! Transpose data
    tmp = 0
    do inode = 1, npart-1
      do jnode = inode + 1, npart
        tmp = np_ghost_neigh(jnode, inode)
        np_ghost_neigh(jnode, inode) = np_ghost_neigh(inode, jnode)
        np_ghost_neigh(inode, jnode) = tmp
      end do
    end do
    vp%np_ghost_neigh_partno(1:npart) = np_ghost_neigh(1:npart, vp%partno)
    SAFE_ALLOCATE(vp%rcounts(1:npart))
    vp%rcounts(1:npart) = np_ghost_neigh(vp%partno, 1:npart)
    ! Set index tables xghost and xghost_neigh. 
    SAFE_ALLOCATE(np_ghost_tmp(1:npart))
    call MPI_Allgather(vp%np_ghost, 1, MPI_INTEGER, &
         np_ghost_tmp(1), 1, MPI_INTEGER, &
         comm, mpi_err)
   
    SAFE_ALLOCATE(xghost_tmp(1:npart))
    xghost_tmp(1) = 1
    do inode = 2, npart
      xghost_tmp(inode) = xghost_tmp(inode - 1) + np_ghost_tmp(inode - 1)
    end do
    vp%xghost = xghost_tmp(vp%partno)

    SAFE_ALLOCATE(xghost_neigh_partno(1:npart))
    SAFE_ALLOCATE(xghost_neigh_back(1:npart))
    tmp = 0
    do inode = 1, npart
      tmp = xghost_tmp(inode)
      xghost_neigh_partno(inode) = xghost_tmp(inode) 
      if (inode == vp%partno) then
        xghost_neigh_back(1)   = xghost_tmp(inode)
      end if
      do jnode = 2, npart
        tmp = tmp + np_ghost_neigh(inode, jnode - 1)
        if (jnode == vp%partno) then
          xghost_neigh_partno(inode) = tmp
        end if
        if (inode == vp%partno) then
          xghost_neigh_back(jnode) = tmp
        end if
      end do
    end do
    SAFE_DEALLOCATE_P(np_ghost_neigh)
    
    ! Get space for ghost point vector.
    SAFE_ALLOCATE(vp%ghost(1:vp%total))

    ! Fill ghost as described above.
    irr = 0
    do ip = 1, np_global+np_enl
      jnode = iihash_lookup(ghost_flag(vp%partno), ip, found)
      ! If point ip is a ghost point for vp%partno from jnode, save this
      ! information.
      if(found) then
        vp%ghost(xghost_neigh_back(jnode) + irr(jnode)) = ip
        irr(jnode) = irr(jnode) + 1
      end if
    end do

    do inode =  1, npart
      do jnode = 1, npart
        if(inode /= jnode) then
          init = xghost_neigh_partno(inode)
          size = vp%np_ghost_neigh_partno(inode)
          call MPI_Bcast(init, 1, MPI_INTEGER, jnode-1, comm, mpi_err)
          call MPI_Bcast(size, 1, MPI_INTEGER, jnode-1, comm, mpi_err)
          call MPI_Bcast(vp%ghost(init), size, MPI_INTEGER, &
               inode-1, comm, mpi_err)
        end if
      end do
    end do

    if(in_debug_mode) then
      ! Write numbers and coordinates of each process` ghost points
      ! to a single file (like in mesh_partition_init) called
      ! debug/mesh_partition/ghost_points.###.
      call io_mkdir('debug/mesh_partition')
      
      write(filenum, '(i3.3)') vp%partno
      iunit = io_open('debug/mesh_partition/ghost_points.'//filenum, action='write')
      do ip = 1, vp%np_ghost
        jp = vp%ghost(xghost_tmp(vp%partno) + ip - 1)
        write(iunit, '(4i8)') jp, (idx%lxyz(jp, idir), idir = 1, MAX_DIM)
      end do

      call io_close(iunit)
    end if

    
    SAFE_ALLOCATE(points(1:np_enl))
    SAFE_ALLOCATE(part_bndry(1:np_enl))
    do ii = 1, np_enl
      points(ii) = ii
    end do
    call partition_get_partition_number(bndry_partition, np_enl, &
         points, part_bndry)
    
    ! Set up the global-to-local point number mapping
    if (periodic_dim /= 0) then
      ip = 1
      jp = npart
      SAFE_ALLOCATE(vp%bndry(1:np_enl))
      irr = 0
      do ii = 1, np_enl
        vp%bndry(xbndry_tmp(part_bndry(ii)) + irr(part_bndry(ii))) = ii + np_global
        irr(part_bndry(ii)) = irr(part_bndry(ii)) + 1 ! increment the counter
      end do
    else
      ip = vp%partno
      jp = vp%partno
      ! initialize to zero all input
      do inode = 1, npart
        if (inode /= vp%partno) then
          call iihash_init(vp%global(inode),1)
        end if
      end do
      ii = xbndry_tmp(vp%partno) + np_bndry_tmp(vp%partno)
      SAFE_ALLOCATE(vp%bndry(xbndry_tmp(vp%partno):ii))
      tmp = 0
      do ii = 1, np_enl
        if(part_bndry(ii) == vp%partno) then
          vp%bndry(xbndry_tmp(part_bndry(ii)) + tmp) = ii + np_global
          tmp = tmp + 1 ! increment the counter
        end if
      end do
    end if  
    SAFE_DEALLOCATE_P(part_bndry)
    SAFE_DEALLOCATE_P(points)
    SAFE_DEALLOCATE_A(irr)
    
    do inode = ip, jp
      ! Create hash table.
      call iihash_init(vp%global(inode), vp%np_local_vec(inode) + &
           np_ghost_tmp(inode) + np_bndry_tmp(inode))
      ! Insert local points.
      do kp = 1, vp%np_local_vec(inode)
        call iihash_insert(vp%global(inode), vp%local_vec(vp%xlocal_vec(inode) + kp - 1), kp)
      end do
      ! Insert ghost points.
      do kp = 1, np_ghost_tmp(inode)
        call iihash_insert(vp%global(inode), vp%ghost(xghost_tmp(inode) + kp - 1), kp + vp%np_local_vec(inode))
      end do
      ! Insert boundary points.
      do kp = 1, np_bndry_tmp(inode)
        call iihash_insert(vp%global(inode), vp%bndry(xbndry_tmp(inode) + kp - 1), &
             kp + vp%np_local_vec(inode) + np_ghost_tmp(inode))
      end do
    end do
    vp%xbndry = xbndry_tmp(vp%partno)
    SAFE_DEALLOCATE_P(np_ghost_tmp)
    SAFE_DEALLOCATE_P(np_bndry_tmp)
    SAFE_DEALLOCATE_P(xbndry_tmp)
    SAFE_DEALLOCATE_P(xghost_tmp)    
    
    ! Complete entries in vp.
    vp%comm      = comm
    vp%root      = root
    vp%np_global = np_global
    vp%npart     = npart

    call init_send_points

    do inode = 1, npart
      call iihash_end(ghost_flag(inode))
    end do
    SAFE_DEALLOCATE_A(ghost_flag)

    POP_SUB(vec_init)

  contains
    subroutine init_local
      integer :: sp, ep, np_tmp
      integer, allocatable :: local_tmp(:), xlocal_tmp(:)
      PUSH_SUB(vec_init.init_local)
      
      sp = vp%xlocal
      ep = vp%xlocal + vp%np_local + 1
      SAFE_ALLOCATE(vp%local(sp:ep))

      sp = 1
      ep = np_global
      SAFE_ALLOCATE(vp%local_vec(sp:ep))

      ! Calculate the local vector in parallel
      call partition_get_local(inner_partition, local_tmp, np_tmp)

      ! Add padding to the calculated local vector
      do ip = 1, np_tmp
        vp%local(vp%xlocal + ip - 1) = local_tmp(ip)
      end do
      
      SAFE_ALLOCATE(xlocal_tmp(1:npart))
      xlocal_tmp = vp%xlocal_vec - 1
      ! Gather all the local vectors in a unique big one
      call mpi_debug_in(comm, C_MPI_ALLGATHERV)
      call MPI_Allgatherv(vp%local(vp%xlocal), vp%np_local, MPI_INTEGER, &
                          vp%local_vec, vp%np_local_vec, xlocal_tmp,  MPI_INTEGER, &
                          comm, mpi_err)
      call mpi_debug_out(comm, C_MPI_GATHERV)
      SAFE_DEALLOCATE_A(xlocal_tmp)

      POP_SUB(vec_init.init_local)
    end subroutine init_local

    subroutine init_MPI_Alltoall
      integer :: ipg
      PUSH_SUB(vec_init.init_MPI_Alltoall)
      
      SAFE_ALLOCATE(vp%recv_count(1:npart))
      SAFE_ALLOCATE(points(1:vp%np_local))
      vp%recv_count = 0
      do ip = 1, vp%np_local
        ! Get the temporally global point
        ipg = vp%local(vp%xlocal + ip - 1)
        ! Get the destination global point
        points(ip) = vp%local_vec(ipg)
      end do

      SAFE_ALLOCATE(vp%part_local_rev(1:vp%np_local))
      ! Get the destination partitions
      call partition_get_partition_number(inner_partition, vp%np_local, points, vp%part_local_rev)
      SAFE_DEALLOCATE_P(points)

      do ip = 1, vp%np_local
        ipart = vp%part_local_rev(ip)
        vp%recv_count(ipart) = vp%recv_count(ipart) + 1
      end do
      
      SAFE_ALLOCATE(vp%send_disp(1:npart))
      SAFE_ALLOCATE(vp%recv_disp(1:npart))

      vp%send_disp(1) = 0
      vp%recv_disp(1) = 0
      do ipart = 2, npart
        vp%send_disp(ipart) = vp%send_disp(ipart - 1) + vp%send_count(ipart - 1)
        vp%recv_disp(ipart) = vp%recv_disp(ipart - 1) + vp%recv_count(ipart - 1)
      end do
      
      POP_SUB(vec_init.init_MPI_Alltoall)
    end subroutine init_MPI_Alltoall
    
    subroutine init_send_points
      integer, allocatable :: displacements(:)
      integer :: ii, jj, kk, ipart, total

      PUSH_SUB(vec_init.init_send_points)

      SAFE_ALLOCATE(vp%sendpos(1:vp%npart))

      total = sum(vp%np_ghost_neigh_partno(1:vp%npart))

      SAFE_ALLOCATE(displacements(1:total))
        
      jj = 0
      ! Iterate over all possible receivers.
      do ipart = 1, vp%npart
        vp%sendpos(ipart) = jj + 1
        ! Iterate over all ghost points that ipart wants.
        do ii = 0, vp%np_ghost_neigh_partno(ipart) - 1
          ! Get global number kk of i-th ghost point.
          kk = vp%ghost(xghost_neigh_partno(ipart) + ii)
          ! Lookup up local number of point kk
          jj = jj + 1
          displacements(jj) = vec_global2local(vp, kk, vp%partno)
        end do
      end do

      call subarray_init(vp%sendpoints, total, displacements)

      SAFE_DEALLOCATE_A(displacements)
      SAFE_DEALLOCATE_P(xghost_neigh_partno)
        
      ! Send and receive displacements.
      ! Send displacement cannot directly be calculated
      ! from vp%xghost_neigh because those are indices for
      ! vp%np_ghost_neigh(vp%partno, :) and not
      ! vp%np_ghost_neigh(:, vp%partno) (rank being fixed).
      ! So what gets done is to pick out the number of ghost points
      ! each partition r wants to have from the current partiton
      ! vp%partno.
      
      SAFE_ALLOCATE(vp%sdispls(1:vp%npart))
      SAFE_ALLOCATE(vp%rdispls(1:vp%npart))

      vp%sdispls(1) = 0
      do ipart = 2, vp%npart
        vp%sdispls(ipart) = vp%sdispls(ipart - 1) + vp%np_ghost_neigh_partno(ipart - 1)
      end do

      ! This is like in vec_scatter/gather.
      vp%rdispls(1:vp%npart) = xghost_neigh_back(1:vp%npart) - vp%xghost
      SAFE_DEALLOCATE_P(xghost_neigh_back)
      
      POP_SUB(vec_init.init_send_points)
    end subroutine init_send_points

  end subroutine vec_init


  ! ---------------------------------------------------------
  !> Deallocate memory used by vp.
  subroutine vec_end(vp)
    type(pv_t), intent(inout) :: vp

    integer :: ipart

    PUSH_SUB(vec_end)

    call subarray_end(vp%sendpoints)

    SAFE_DEALLOCATE_P(vp%rdispls)
    SAFE_DEALLOCATE_P(vp%sdispls)
    SAFE_DEALLOCATE_P(vp%rcounts)
    SAFE_DEALLOCATE_P(vp%sendpos)    
    SAFE_DEALLOCATE_P(vp%send_disp)
    SAFE_DEALLOCATE_P(vp%recv_disp)
    SAFE_DEALLOCATE_P(vp%part_vec)
    SAFE_DEALLOCATE_P(vp%part_local)
    SAFE_DEALLOCATE_P(vp%part_local_rev)
    SAFE_DEALLOCATE_P(vp%np_local_vec)
    SAFE_DEALLOCATE_P(vp%xlocal_vec)
    SAFE_DEALLOCATE_P(vp%local)
    SAFE_DEALLOCATE_P(vp%local_vec)
    SAFE_DEALLOCATE_P(vp%send_count)
    SAFE_DEALLOCATE_P(vp%recv_count)
    SAFE_DEALLOCATE_P(vp%bndry)
    SAFE_DEALLOCATE_P(vp%np_ghost_neigh_partno)
    SAFE_DEALLOCATE_P(vp%ghost)

    if(associated(vp%global)) then 
      do ipart = 1, vp%npart
        call iihash_end(vp%global(ipart))
      end do
      SAFE_DEALLOCATE_P(vp%global)
    end if

    POP_SUB(vec_end)

  end subroutine vec_end

  ! ---------------------------------------------------------
  !> Returns local number of global point ip on partition inode.
  !! If the result is zero, the point is neither a local nor a ghost
  !! point on inode.
  integer function vec_global2local(vp, ip, inode)
    type(pv_t), intent(in) :: vp
    integer,    intent(in) :: ip
    integer,    intent(in) :: inode

    integer :: nn
    logical :: found

! no push_sub because called too frequently

    vec_global2local = 0
    if(associated(vp%global)) then
      nn = iihash_lookup(vp%global(inode), ip, found)
      if(found) vec_global2local = nn
    end if

  end function vec_global2local

  !> Change the value of one dimension (1=x, 2=y, 3=z) 
  !! according to the given value and return the local point
  !! \todo does not work if only one process is used
  integer function vec_index2local(vp, idx, ix, dim_pad, pad)
    type(pv_t),    intent(in) :: vp      !< All the required information
    type(index_t), intent(in) :: idx     !< Index information
    integer,       intent(in) :: ix(1:3) !< Global x,y,z indices
    integer,       intent(in) :: dim_pad !< The dimension that has to be changed
    integer,       intent(in) :: pad     !< How much we want to change that dimension

    integer :: global_point, local_point
    integer :: jx(1:3)

    ! no PUSH_SUB, called too often
    jx = ix
    jx(dim_pad) = jx(dim_pad) + pad
    global_point = index_from_coords(idx, 3, jx)
    
    if (mpi_world%size == 1) then 
      local_point = global_point
    else
      local_point = vec_global2local(vp, global_point, vp%partno)
      if (local_point == 0) then
        write(message(1), '(a)') "You are trying to access a neighbour that does not exist."
        write(message(2), '(a, i5)') "Global point = ", global_point
        write(message(3), '(a, 3i5)') "x,y,z point  = ", jx
        call messages_warning(3)
      end if
    end if
    
    vec_index2local = local_point
    
  end function vec_index2local
  

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
