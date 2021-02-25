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
  !! call vec_scatter(vp, uu, ul)
  !! call vec_scatter(vp, vv, vl)
  !!
  !! ! Compute some operator op: vl = op ul
  !! call X(vec_ghost_update)(vp, ul)
  !! call X(nl_operator_operate)(op, ul, vl)
  !! !! Gather result of op in one vector vv.
  !! call vec_gather(vp, vv, vl)
  !!
  !! ! Clean up.
  !! deallocate(ul, vl, wl)
  !! \endverbatim
module par_vec_oct_m
  use global_oct_m
  use iihash_oct_m
  use index_oct_m
  use io_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use namespace_oct_m
  use partition_oct_m
  use profiling_oct_m
  use space_oct_m
  use stencil_oct_m
  use subarray_oct_m

  implicit none

  private

  public ::            &
    pv_t,              & !< parallel information
    vec_init,          &
    vec_end,           &
    vec_scatter,       &
    vec_gather,        &
    vec_allgather,     &
    vec_global2local    
  
  !> Parallel information
  type pv_t
    ! Components are public by default

    ! The content of these members is process-dependent.
    integer              :: rank                 !< Our rank in the communicator. 
    !> Partition number of the
    !! current process
    integer              :: partno              
    type(subarray_t)     :: ghost_spoints     !< Ghost points that actual process must send to the neighbours
    integer, allocatable :: ghost_sendpos(:)  !< The positions of the points for the ghost communication

    integer, allocatable :: ghost_rdispls(:)  !< Ghost points receive displacements 
    integer, allocatable :: ghost_sdispls(:)  !< Ghost points send displacements 
    integer, allocatable :: ghost_rcounts(:)  !< Number of ghost points to receive
    integer, allocatable :: ghost_scounts(:)  !< Number of ghost points to send
    integer              :: ghost_scount
    integer, allocatable :: ghost_sendmap(:)  !< map for packing ghost points
    integer, allocatable :: ghost_recvmap(:)  !< map for unpacking ghost points

    ! The following members are set independent of the processs.
    integer                 :: npart                !< Number of partitions.
    integer                 :: comm                 !< MPI communicator to use.
    integer                 :: np_global            !< Number of points in mesh.
    integer, allocatable    :: part_local(:)        !< Local point         -> partition
    integer, allocatable    :: part_local_rev(:)    !< Local point`s value -> partition

    integer, allocatable    :: np_local_vec(:)      !< How many points has partition r?
                                                    !! Global vector; npart elements.
    integer                 :: np_local             !< How many points has running partition? 
                                                    !! Local value.
    integer, allocatable    :: xlocal_vec(:)        !< Points of partition r start at
                                                    !! xlocal_vec(r) in local. Global start point
                                                    !! of the local index.  
                                                    !! Global vector; npart elements.
    integer                 :: xlocal               !< Starting index of running process in local(:) vector.
                                                    !! Local value.
          
    integer, allocatable    :: local_vec(:)         !< Partition r has points
                                                    !! local_vec(xlocal_vec(r):
                                                    !! xlocal_vec(r)+np_local_vec(r)-1). 
                                                    !! Global vector; np_global elements    
    integer, allocatable    :: local(:)             !< Local points of running process
                                                    !! Local vector; np_local elements
    integer, allocatable    :: recv_count(:)        !< Number of points to receive from all the other processes
    integer, allocatable    :: send_count(:)        !< Number of points to send to all the other processes
                                                    !! in a MPI_Alltoallv.
    integer, allocatable    :: recv_disp(:)         !< Displacement of points to receive from all the other processes
    integer, allocatable    :: send_disp(:)         !< Displacement of points to send to all the other processes
                                                    !! in a MPI_Alltoallv.
    integer                 :: np_bndry             !< Number of boundary points.
                                                    !! Local value
    integer                 :: xbndry               !< Starting index of running process in bndry(:) 
                                                    !! Local value
 
    integer, allocatable    :: bndry(:)             !< Global numbers of boundary points.
                                                    !! Global vector; np_enl elements
      
    type(iihash_t), private :: global               !< global(r) contains the global ->
                                                        !! local mapping
    integer, private        :: total                    !< Total number of ghost points. 

    integer                 :: np_ghost                 !< How many ghost points has partition r?
                                                        !! Local value
    integer                 :: xghost                   !< Starting index of running procces in ghost(:) vector.
    integer, allocatable    :: ghost(:)                 !< Global indices of all local points.
                                                        !! Global vector; vp%total elements
  end type pv_t

  interface vec_scatter
    module procedure dvec_scatter
    module procedure zvec_scatter
    module procedure ivec_scatter
  end interface vec_scatter

  interface vec_gather
    module procedure dvec_gather
    module procedure zvec_gather
    module procedure ivec_gather
  end interface vec_gather

  interface vec_allgather
    module procedure dvec_allgather
    module procedure zvec_allgather
    module procedure ivec_allgather
  end interface vec_allgather
  
contains

  !> Initializes a pv_type object (parallel vector).
  !! It computes the local-to-global and global-to-local index tables
  !! and the ghost point exchange.
  !! \warning The naming scheme for the np_ variables is different
  !! from how it is in the rest of the code (for historical reasons
  !! and also because the vec_init has more a global than local point
  !! of view on the mesh): See the comments in the parameter list.
  subroutine vec_init(comm, np_global, np_part_global, idx, stencil, space, inner_partition, bndry_partition, vp, namespace)
    integer,         intent(in)  :: comm         !< Communicator to use.

    !> The next seven entries come from the mesh.
    integer,           intent(in)    :: np_global      !< mesh%np_global
    integer,           intent(in)    :: np_part_global !< mesh%np_part_global
    type(index_t),     intent(in)    :: idx
    type(stencil_t),   intent(in)    :: stencil        !< The stencil for which to calculate ghost points.
    type(space_t),     intent(in)    :: space
    type(partition_t), intent(in)    :: inner_partition
    type(partition_t), intent(in)    :: bndry_partition
    type(pv_t),        intent(inout) :: vp             !< Description of partition.
    type(namespace_t), intent(in)    :: namespace

    ! Careful: MPI counts process ranks from 0 to numproc-1.
    ! Partition numbers from METIS range from 1 to numproc.
    ! For this reason, all ranks are incremented by one.
    integer                     :: npart            !< Number of partitions.
    integer                     :: np_enl           !< Number of points in enlargement.
    integer                     :: gip, ip, jp, kp, jj, index, inode, jnode !< Counters.
    integer, allocatable        :: irr(:)           !< Counter.
    integer                     :: p1(MAX_DIM)      !< Points.
    type(iihash_t), allocatable :: ghost_flag(:)    !< To remember ghost pnts.
    type(iihash_t)              :: boundary, boundary_inv
    type(iihash_t)              :: ghost, ghost_inv
    integer                     :: iunit            !< For debug output to files.
    character(len=6)            :: filenum
    integer                     :: tmp, init, size, ii
    integer, allocatable        :: init_v(:), size_v(:), init_recv(:), sbuffer(:)
    logical                     :: found
   
    integer                     :: idir, ipart, np_inner, np_bndry
    integer, allocatable        :: np_ghost_tmp(:), np_bndry_tmp(:)
    !> Number of ghost points per
    !! neighbour per partition.      
    integer, allocatable        :: xbndry_tmp(:)    !< Starting index of process i in bndry(:). 
    integer, allocatable        :: xghost_tmp(:)
    integer, allocatable        :: xghost_neigh_partno(:)   !< Like xghost for neighbours.
    integer, allocatable        :: xghost_neigh_back(:)     !< Same as previous, but outward
    integer, allocatable        :: points(:), points_bndry(:), part_bndry(:), part_inner(:), part_ghost(:)

    PUSH_SUB(vec_init)

    ! Shortcuts.
#ifdef HAVE_MPI
    call MPI_Comm_Size(comm, npart, mpi_err)
#endif
    np_enl = np_part_global - np_global

    ! Store partition number and rank for later reference.
    ! Having both variables is a bit redundant but makes the code readable.
#ifdef HAVE_MPI
    call MPI_Comm_Rank(comm, vp%rank, mpi_err)
#endif
    vp%partno = vp%rank + 1

    vp%comm      = comm
    vp%np_global = np_global
    vp%npart     = npart


    SAFE_ALLOCATE(ghost_flag(1:npart))
    SAFE_ALLOCATE(irr(1:npart))
    SAFE_ALLOCATE(vp%np_local_vec(1:npart))
    SAFE_ALLOCATE(vp%xlocal_vec(1:npart))
    SAFE_ALLOCATE(vp%ghost_rcounts(1:npart))
    SAFE_ALLOCATE(vp%ghost_scounts(1:npart))

    ! Count number of points for each process.
    ! Local points.
    call partition_get_np_local(inner_partition, vp%np_local_vec)
    vp%np_local = vp%np_local_vec(vp%partno)


    ! Set up local-to-global index table for local points
    ! (xlocal_vec, local) and for boundary points (xbndry, bndry).
    vp%xlocal_vec(1) = 1
    ! Set the starting point of local and boundary points
    do inode = 2, npart
      vp%xlocal_vec(inode) = vp%xlocal_vec(inode - 1) + vp%np_local_vec(inode - 1)
    end do
    vp%xlocal = vp%xlocal_vec(vp%partno)

    call init_local()

    ! Create hash table.
    call iihash_init(vp%global)
    ! Insert local points.
    do ip = 1, vp%np_local
      call iihash_insert(vp%global, vp%local(vp%xlocal + ip - 1), ip)
    end do

    SAFE_ALLOCATE(points(1:vp%np_local))
    SAFE_ALLOCATE(vp%part_local(1:vp%np_local))
    do ip = 1, vp%np_local
      points(ip) = vp%xlocal + ip - 1
    end do
    call partition_get_partition_number(inner_partition, vp%np_local, &
         points, vp%part_local)
    SAFE_DEALLOCATE_A(points)

    call iihash_init(boundary)
    call iihash_init(ghost)
    call iihash_init(boundary_inv)
    call iihash_init(ghost_inv)
    vp%np_ghost = 0
    np_bndry = 0
    do gip = vp%xlocal, vp%xlocal + vp%np_local - 1
      ! Get coordinates of current point.
      call index_to_coords(idx, vp%local(gip), p1)
      ! For all points in stencil.
      do jj = 1, stencil%size
        ! Get point number of possible ghost point.
        index = index_from_coords(idx, p1(:) + stencil%points(:, jj))
        ASSERT(index /= 0)
        ! check if this point is a local point
        tmp = iihash_lookup(vp%global, index, found)
        if (found) cycle
        ! now check if the point is a potential ghost or boundary point
        if (index > np_global) then
          tmp = iihash_lookup(boundary, index, found)
          if (found) cycle
          vp%np_bndry = vp%np_bndry + 1
          call iihash_insert(boundary, index, vp%np_bndry)
          call iihash_insert(boundary_inv, vp%np_bndry, index)
        else
          tmp = iihash_lookup(ghost, index, found)
          if (found) cycle
          vp%np_ghost = vp%np_ghost + 1
          call iihash_insert(ghost, index, vp%np_ghost)
          call iihash_insert(ghost_inv, vp%np_ghost, index)
        end if
      end do
    end do

    SAFE_ALLOCATE(vp%bndry(1:vp%np_bndry))
    do ip = 1, vp%np_bndry
      vp%bndry(ip) = iihash_lookup(boundary_inv, ip, found)
      ASSERT(found)
    end do

    SAFE_ALLOCATE(vp%ghost(1:vp%np_ghost))
    do ip = 1, vp%np_ghost
      vp%ghost(ip) = iihash_lookup(ghost_inv, ip, found)
      ASSERT(found)
    end do
    call iihash_end(ghost_inv)
    call iihash_end(boundary_inv)
    call iihash_end(ghost)
    call iihash_end(boundary)

    SAFE_ALLOCATE(part_ghost(1:vp%np_ghost))
    call partition_get_partition_number(inner_partition, vp%np_ghost, &
         vp%ghost, part_ghost)

    vp%ghost_rcounts(:) = 0
    do ip = 1, vp%np_ghost
      ipart = part_ghost(ip)
      vp%ghost_rcounts(ipart) = vp%ghost_rcounts(ipart)+1
    end do
    ASSERT(sum(vp%ghost_rcounts) == vp%np_ghost)

    call init_MPI_Alltoall()
    
#ifdef HAVE_MPI
    call MPI_Alltoall(vp%ghost_rcounts(1), 1, MPI_INTEGER, &
         vp%ghost_scounts(1), 1, MPI_INTEGER, &
         comm, mpi_err)
#endif

    SAFE_ALLOCATE(vp%ghost_sdispls(1:vp%npart))
    SAFE_ALLOCATE(vp%ghost_rdispls(1:vp%npart))

    vp%ghost_sdispls(1) = 0
    vp%ghost_rdispls(1) = 0
    do ipart = 2, vp%npart
      vp%ghost_sdispls(ipart) = vp%ghost_sdispls(ipart - 1) + vp%ghost_scounts(ipart - 1)
      vp%ghost_rdispls(ipart) = vp%ghost_rdispls(ipart - 1) + vp%ghost_rcounts(ipart - 1)
    end do
    vp%ghost_scount = sum(vp%ghost_scounts)

    SAFE_ALLOCATE(vp%ghost_recvmap(1:vp%np_ghost))
    SAFE_ALLOCATE(points(1:vp%npart))
    points = 0
    do ip = 1, vp%np_ghost
      ipart = part_ghost(ip)
      points(ipart) =  points(ipart) + 1
      vp%ghost_recvmap(ip) = vp%ghost_rdispls(ipart) + points(ipart)
    end do
    SAFE_DEALLOCATE_A(points)

    SAFE_ALLOCATE(points(1:vp%np_ghost))
    do ip = 1, vp%np_ghost
      points(vp%ghost_recvmap(ip)) = vp%ghost(ip)
    end do
    SAFE_ALLOCATE(vp%ghost_sendmap(1:vp%ghost_scount))
#ifdef HAVE_MPI
    call MPI_Alltoallv(points(1), vp%ghost_rcounts(1), vp%ghost_rdispls(1), MPI_INTEGER, &
         vp%ghost_sendmap(1), vp%ghost_scounts(1), vp%ghost_sdispls(1), MPI_INTEGER, &
         vp%comm, mpi_err)
#endif
    print*, vp%partno, vp%ghost_scount, vp%np_ghost
    do ip = 1, vp%ghost_scount
      ! get local index
      index = iihash_lookup(vp%global, vp%ghost_sendmap(ip), found)
      ASSERT(found)
      vp%ghost_sendmap(ip) = index
    end do

    if(debug%info) then
      ! Write numbers and coordinates of each process` ghost points
      ! to a single file (like in mesh_partition_init) called
      ! debug/mesh_partition/ghost_points.###.
      call io_mkdir('debug/mesh_partition', namespace)
      
      write(filenum, '(i6.6)') vp%partno
      iunit = io_open('debug/mesh_partition/ghost_points.'//filenum, namespace, action='write')
      do ip = 1, vp%np_ghost
        jp = vp%ghost(ip)
        call index_to_coords(idx, jp, p1)
        write(iunit, '(99i8)') jp, (p1(idir), idir = 1, space%dim)
      end do

      call io_close(iunit)
    end if

    ! Insert ghost points.
    do ip = 1, vp%np_ghost
      call iihash_insert(vp%global, vp%ghost(ip), ip + vp%np_local)
    end do
    ! Insert boundary points.
    do ip = 1, vp%np_bndry
      call iihash_insert(vp%global, vp%bndry(ip), ip + vp%np_local + vp%np_ghost)
    end do
    
    !call init_send_points()

    SAFE_DEALLOCATE_A(part_ghost)

    POP_SUB(vec_init)

  contains
    subroutine init_local()

      integer :: sp, ep, np_tmp
      integer, allocatable :: xlocal_tmp(:)

      PUSH_SUB(vec_init.init_local)
      
      sp = vp%xlocal
      ep = vp%xlocal + vp%np_local + 1
      SAFE_ALLOCATE(vp%local(sp:ep))

      sp = 1
      ep = np_global
      SAFE_ALLOCATE(vp%local_vec(sp:ep))

      ! Calculate the local vector in parallel
      call partition_get_local(inner_partition, vp%local(vp%xlocal:), np_tmp)

      SAFE_ALLOCATE(xlocal_tmp(1:npart))
      xlocal_tmp = vp%xlocal_vec - 1
      ! Gather all the local vectors in a unique big one
      call mpi_debug_in(comm, C_MPI_ALLGATHERV)
#ifdef HAVE_MPI
      call MPI_Allgatherv(vp%local(vp%xlocal), vp%np_local, MPI_INTEGER, &
                          vp%local_vec, vp%np_local_vec, xlocal_tmp,  MPI_INTEGER, &
                          comm, mpi_err)
#endif
      call mpi_debug_out(comm, C_MPI_GATHERV)
      SAFE_DEALLOCATE_A(xlocal_tmp)

      POP_SUB(vec_init.init_local)
    end subroutine init_local

    subroutine init_MPI_Alltoall()

      integer :: ipg

      PUSH_SUB(vec_init.init_MPI_Alltoall)
      
      SAFE_ALLOCATE(vp%recv_count(1:npart))
      SAFE_ALLOCATE(vp%send_count(1:npart))
      SAFE_ALLOCATE(points(1:vp%np_local))
      vp%send_count = 0
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
      SAFE_DEALLOCATE_A(points)

      do ip = 1, vp%np_local
        ipart = vp%part_local_rev(ip)
        vp%recv_count(ipart) = vp%recv_count(ipart) + 1
        ipart = vp%part_local(ip)
        vp%send_count(ipart) = vp%send_count(ipart) + 1
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
    
    subroutine init_send_points()

      integer, allocatable :: displacements(:)
      integer :: ii, jj, kk, ipart, total

      PUSH_SUB(vec_init.init_send_points)

      SAFE_ALLOCATE(vp%ghost_sendpos(1:vp%npart))

      total = sum(vp%ghost_scounts(1:vp%npart))

      SAFE_ALLOCATE(displacements(1:total))
        
      jj = 0
      !! Iterate over all possible receivers.
      !do ipart = 1, vp%npart
      !  vp%ghost_sendpos(ipart) = jj + 1
      !  ! Iterate over all ghost points that ipart wants.
      !  do ii = 0, vp%ghost_scounts(ipart) - 1
      !    ! Get global number kk of i-th ghost point.
      !    !kk = vp%ghost(xghost_neigh_partno(ipart) + ii)
      !    kk = 1
      !    ! Lookup up local number of point kk
      !    jj = jj + 1
      !    displacements(jj) = vec_global2local(vp, kk)
      !  end do
      !end do

      !call subarray_init(vp%ghost_spoints, total, displacements)

      SAFE_DEALLOCATE_A(displacements)
      SAFE_DEALLOCATE_A(xghost_neigh_partno)

      POP_SUB(vec_init.init_send_points)
    end subroutine init_send_points

  end subroutine vec_init


  ! ---------------------------------------------------------
  !> Deallocate memory used by vp.
  subroutine vec_end(vp)
    type(pv_t), intent(inout) :: vp

    integer :: ipart

    PUSH_SUB(vec_end)

    call subarray_end(vp%ghost_spoints)

    SAFE_DEALLOCATE_A(vp%ghost_rdispls)
    SAFE_DEALLOCATE_A(vp%ghost_sdispls)
    SAFE_DEALLOCATE_A(vp%ghost_rcounts)
    SAFE_DEALLOCATE_A(vp%ghost_scounts)
    SAFE_DEALLOCATE_A(vp%ghost_sendpos)    
    SAFE_DEALLOCATE_A(vp%send_disp)
    SAFE_DEALLOCATE_A(vp%recv_disp)
    SAFE_DEALLOCATE_A(vp%part_local)
    SAFE_DEALLOCATE_A(vp%part_local_rev)
    SAFE_DEALLOCATE_A(vp%np_local_vec)
    SAFE_DEALLOCATE_A(vp%xlocal_vec)
    SAFE_DEALLOCATE_A(vp%local)
    SAFE_DEALLOCATE_A(vp%local_vec)
    SAFE_DEALLOCATE_A(vp%send_count)
    SAFE_DEALLOCATE_A(vp%recv_count)
    SAFE_DEALLOCATE_A(vp%bndry)
    SAFE_DEALLOCATE_A(vp%ghost)

    call iihash_end(vp%global)

    POP_SUB(vec_end)

  end subroutine vec_end


  ! ---------------------------------------------------------
  !> Returns local number of global point ip on the local node
  !! If the result is zero, the point is not available on the local node
  integer function vec_global2local(vp, ip)
    type(pv_t), intent(in) :: vp
    integer,    intent(in) :: ip

#ifdef HAVE_MPI
    integer :: nn
    logical :: found
#endif
    
! no push_sub because called too frequently

#ifdef HAVE_MPI
    
    vec_global2local = 0
    nn = iihash_lookup(vp%global, ip, found)
    if(found) vec_global2local = nn

#else

    vec_global2local = ip

#endif

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

end module par_vec_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
