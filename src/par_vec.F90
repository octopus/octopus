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
  ! When working with non-periodic boundary conditions
  ! a globally defined vector v has two parts:
  ! - v(:np) are the inner points
  ! - v(np+1:np_tot) are the boundary points
  ! In the typical case of zero boundary conditions
  ! v(np+1:np_tot) is 0.
  ! The two parts are split according to the partitions.
  ! The result of this slit are local vectors vl on each node
  ! which consist of three parts:
  ! - v(:np_local)                                      local points.
  ! - v(np_local+1:np_local+np_ghost)                   ghost points.
  ! - v(np_local+np_ghost+1:np_local+np_ghost+np_bndry) boundary points.
  !
  ! 
  ! Usage example for par_vec routines.
  !
  ! ! Initialize parallelization with mesh m and operator op
  ! ! initialized and given.
  ! ! m          = sys%gr%m
  ! ! stencil    = op%stencil
  ! ! np_stencil = op%n
  !
  ! FLOAT              :: s
  ! FLOAT              :: u(np_glob), v(np_glob)
  ! FLOAT, allocatable :: ul(:), vl(:), wl(:)
  ! type(mesh_type)    :: m
  !
  ! ! Fill u, v with sensible values.
  ! ! ...
  !
  ! ! Allocate space for local vectors.
  ! allocate(ul(np_tot))
  ! allocate(vl(np_tot))
  ! allocate(wl(np_tot))
  !
  ! ! Distribute vectors.
  ! call X(vec_scatter)(vp, u, ul)
  ! call X(vec_scatter)(vp, v, vl)
  !
  ! ! Calculate scalar product s=<u|v>.
  ! wl = ul*vl
  ! ! vec_integrate ignores ghost points (i. e. v(np_local+1:)).
  ! s = X(vec_integrate)(vp, wl)
  !
  ! ! Compute some operator op: vl = op ul
  ! call X(vec_ghost_update)(vp, ul)
  ! call X(nl_operator_operate)(op, ul, vl)
  ! ! Gather result of op in one vector v.
  ! call X(vec_gather)(vp, v, vl)
  !
  ! ! Clean up.
  ! deallocate(ul, vl, wl)

  use global
  use messages
#ifdef DEBUG
  use io
#endif

  implicit none

#if defined(HAVE_MPI) && defined(HAVE_METIS)
  private
  public :: pv_type,            &
            vec_init,           &
            vec_init_default,   &
            vec_end,            &
            dvec_scatter,       &
            zvec_scatter,       &
            dvec_scatter_bndry, &
            zvec_scatter_bndry, &
            dvec_scatter_all,   &
            zvec_scatter_all,   &
            dvec_gather,        &
            zvec_gather,        &
            dvec_ghost_update,  &
            zvec_ghost_update,  &
            dvec_integrate,     &
            zvec_integrate

  type pv_type
    integer          :: comm                 ! MPI communicator to use.
    integer          :: root                 ! The master node.
    integer          :: p                    ! Number of partitions.
    integer          :: np                   ! Number of points in mesh.
    integer          :: np_enl               ! Number of points in enlargement.
    integer, pointer :: np_local(:)          ! How many points has partition r?
    integer, pointer :: xlocal(:)            ! Points of partition r start at
                                             ! xlocal(r) in local.
    integer, pointer :: local(:)             ! Partition r has points
                                             ! local(xlocal(r):
                                             ! xlocal(r)+np_local(r)-1).
    integer, pointer :: np_bndry(:)          ! Number of boundary points.
    integer, pointer :: xbndry(:)            ! Index of bndry(:).
    integer, pointer :: bndry(:)             ! Global numbers of boundary
                                             ! points.
    integer, pointer :: global(:, :)         ! global(i, r) is local number
                                             ! of point i in partition r
                                             ! (if this is 0, i is neither
                                             ! a ghost point nor local to r).
    integer          :: total                ! Total number of ghost points.
    integer, pointer :: np_ghost(:)          ! How many ghost points has
                                             ! partition r?
    integer, pointer :: np_ghost_neigh(:, :) ! Number of ghost points per
                                             ! neighbour per partition.
    integer, pointer :: xghost(:)            ! Like xlocal.
    integer, pointer :: xghost_neigh(:, :)   ! Like xghost for neighbours.
    integer, pointer :: ghost(:)             ! Global indices of all local
                                             ! ghost points.
  end type pv_type


contains

  ! Wrapper: calls vec_init with MPI_COMM_WORLD and master node 0.
  subroutine vec_init_default(p, part, np, np_tot, nr,                 &
                              Lxyz_inv, Lxyz, stencil, np_stencil, vp)
    ! The next seven entries come from the mesh.
    integer,       intent(in)  :: p            ! m%npart
    integer,       intent(in)  :: part(:)      ! m%part
    integer,       intent(in)  :: np           ! m%np_glob
    integer,       intent(in)  :: np_tot       ! m%np_tot_glob
    integer,       intent(in)  :: nr(2, 3)     ! m%nr
    integer,       intent(in)  :: Lxyz_inv(nr(1,1):nr(2,1), &
                                           nr(1,2):nr(2,2), &
                                           nr(1,3):nr(2,3))
                                               ! m%Lxyz_inv
    integer,       intent(in)  :: Lxyz(:, :)   ! m%Lxyz
    integer,       intent(in)  :: stencil(:,:) ! The stencil for which to
                                               ! calculate ghost points.
    integer,       intent(in)  :: np_stencil   ! Num. of points in stencil.
    type(pv_type), intent(out) :: vp           ! Description of partition.

    ! Throw them away.
    integer :: np_local ! vp%np_local(rank+1).
    integer :: np_ghost ! vp%np_ghost(rank+1).
    integer :: np_bndry ! vp%np_bndry(rank+1).

    call push_sub('par_vec.vec_init_default')
    
    call vec_init(MPI_COMM_WORLD, 0, p, part, np, np_tot, nr,  &
                  Lxyz_inv, Lxyz, stencil, np_stencil, vp,     &
                  np_local, np_ghost, np_bndry)

    call pop_sub()

  end subroutine vec_init_default
 

  ! Initializes a pv_type object (parallel vector).
  ! It computes the local to global and global to local index tables
  ! and the ghost point exchange.
  ! The format for the stencil is: stencil(3, i) for i=1, ..., np_stencil
  ! (as for type(nl_operator_type)).
  ! For example a stencil like (in x-y-plane)
  !          .
  !        .....
  !          .
  ! is coded as
  !   stencil(:, 1) = (/ 0,  1,  0/)
  !   stencil(:, 2) = (/-2,  0,  0/)
  !   stencil(:, 3) = (/-1,  0,  0/)
  !   stencil(:, 3) = (/ 0,  0,  0/)
  !   stencil(:, 5) = (/ 1,  0,  0/)
  !   stencil(:, 6) = (/ 2,  0,  0/)
  !   stencil(:, 7) = (/ 0, -1,  0/)
  ! The points are relative to the "application point" of the stencil.
  ! (This is the same format as used in type(nl_operator_type), so
  ! just passing op%stencil is possible.)
  subroutine vec_init(comm, root, p, part, np, np_tot, nr,     &
                      Lxyz_inv, Lxyz, stencil, np_stencil, vp, &
                      np_local, np_ghost, np_bndry)
    integer,       intent(in)  :: comm         ! Communicator to use.
    integer,       intent(in)  :: root         ! The master node.
    ! The next seven entries come from the mesh.
    integer,       intent(in)  :: p            ! m%npart
    integer,       intent(in)  :: part(:)      ! m%part
    integer,       intent(in)  :: np           ! m%np_glob
    integer,       intent(in)  :: np_tot       ! m%np_tot_glob
    integer,       intent(in)  :: nr(2, 3)     ! m%nr
    integer,       intent(in)  :: Lxyz_inv(nr(1,1):nr(2,1), &
                                           nr(1,2):nr(2,2), &
                                           nr(1,3):nr(2,3))
                                               ! m%Lxyz_inv
    integer,       intent(in)  :: Lxyz(:, :)   ! m%Lxyz
    integer,       intent(in)  :: stencil(:,:) ! The stencil for which to
                                               ! calculate ghost points.
    integer,       intent(in)  :: np_stencil   ! Num. of points in stencil.
    type(pv_type), intent(out) :: vp           ! Description of partition.
    ! Those three are shortcuts.
    integer,       intent(out) :: np_local     ! vp%np_local(rank+1).
    integer,       intent(out) :: np_ghost     ! vp%np_ghost(rank+1).
    integer,       intent(out) :: np_bndry     ! vp%np_bndry(rank+1).
    
    ! Careful: MPI counts node ranks from 0 to numproc-1.
    ! Partition numbers from METIS range from 1 to numproc.
    ! For this reason, all ranks are incremented by one.
    integer              :: np_enl           ! Number of points in enlargement.
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

    ! Shortcuts.
    np_enl = np_tot-np
    call MPI_Comm_rank(comm, rank, ierr)
    
    allocate(ghost_flag(np, p))
    allocate(ir(p), irr(p, p))

    allocate(vp%np_local(p))
    allocate(vp%xlocal(p))
    allocate(vp%local(np))
    allocate(vp%np_bndry(p))
    allocate(vp%xbndry(p))
    allocate(vp%bndry(np_enl))
    allocate(vp%global(np+np_enl, p))
    allocate(vp%np_ghost(p))
    allocate(vp%np_ghost_neigh(p, p))
    allocate(vp%xghost(p))
    allocate(vp%xghost_neigh(p, p))
    
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
      ir(part(i))                                = ir(part(i))+1
    end do
    ir = 0
    do i = np+1, np+np_enl
      vp%bndry(vp%xbndry(part(i))+ir(part(i))) = i
      ir(part(i))                                = ir(part(i))+1
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
        p1 = Lxyz(vp%local(i), :)
        ! For all points in stencil.
        do j = 1, np_stencil
          ! Get coordinates of possible ghost point.
          p2 = p1+stencil(:, j)
          ! Check, whether p2 is in the box.
          ! If p2 is out of the box, ignore it.
          ! Actually, the box should be big enough, that
          ! this case does not arise.
          if(nr(1, 1).gt.p2(1).or.p2(1).gt.nr(2, 1).or. &
             nr(1, 2).gt.p2(2).or.p2(2).gt.nr(2, 2).or. &
             nr(1, 3).gt.p2(3).or.p2(3).gt.nr(2, 3)) cycle
          ! If it is in the box, get its (global) point number.
          k = Lxyz_inv(p2(1), p2(2), p2(3))
          ! If p2 is in the box but does not belong to the inner mesh,
          ! its point number k is greater than np.
          ! If this is the case, it is a boundary point and is handled
          ! elsewhere.
          if(k.gt.np) cycle
          ! At this point, it is sure that point number k is a
          ! relevant point.
          ! If this index k does not belong to partition of node r,
          ! then k is a ghost point for r with part(k) now being
          ! a neighbour of r.
          if(part(k).ne.r) then
            ! Only mark and count this ghost point, if it is not
            ! done yet. Otherwise, points would possibly be registered
            ! more than once.
            if(ghost_flag(k, r).eq.0) then
              ! Mark point i as ghost point for r from part(k).
              ghost_flag(k, r)                = part(k)
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
          write(iunit, '(4i8)') j, Lxyz(j, :)
        end do
      call io_close(iunit)
      end do
    end if
#endif

    ! Create reverse (global to local) lookup.
    ! Given a global point number i and a vector v_local of
    ! length vp%np_local(r)+vp%np_ghost(r)+vp%np_bndry(r) global(i, r) gives
    ! the index of point i in v_local as long as this point is
    ! local to r or a ghost point for r (if vp%np_bndry(r) >
    ! global(i, r) > vp%np_local(r) it is a ghost point).
    ! If global(i, r) is 0 then i is neither local
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
      ! Boundary points.
      do i = 1, vp%np_bndry(r)
        vp%global(vp%bndry(vp%xbndry(r)+i-1), r) =  vp%np_local(r)   &
                                                   +vp%np_ghost(r)+i
      end do
    end do

    ! Complete entries in vp.
    vp%comm   = comm
    vp%root   = root
    vp%np     = np
    vp%np_enl = np_enl
    vp%p      = p

    ! Return shortcuts.
    np_local = vp%np_local(rank+1)
    np_ghost = vp%np_ghost(rank+1)
    np_bndry = vp%np_bndry(rank+1)
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
    if(associated(vp%np_bndry)) then
      deallocate(vp%np_bndry)
      nullify(vp%np_bndry)
    endif
    if(associated(vp%xbndry)) then
      deallocate(vp%xbndry)
      nullify(vp%xbndry)
    endif
    if(associated(vp%bndry)) then
      deallocate(vp%bndry)
      nullify(vp%bndry)
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

#include "undef.F90"
#include "complex.F90"
#include "par_vec_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "par_vec_inc.F90"

#endif
end module par_vec
