!! Copyright (C) 2010 X. Andrade
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
!! $Id: mesh_alltoall.F90 6792 2010-07-08 08:06:12Z xavier $

#include "global.h"

module mesh_alltoall_m
  use batch_m
  use global_m
  use messages_m
  use mesh_m
  use mpi_m
  use mpi_debug_m
  use par_vec_m
  use profiling_m
  use subarray_m

  implicit none
  
  private

  ! WARNING: This module has not been tested. So there is little
  ! probability that it actually works.

#if defined(HAVE_MPI)
  public ::                        &
    mesh_alltoall_init,            &
    mesh_alltoall_end,             &
    dmesh_alltoall_batch_start,    &
    zmesh_alltoall_batch_start,    &
    dmesh_alltoall_batch_finish,   &
    zmesh_alltoall_batch_finish

  integer :: SEND = 1, RECV = 2

  type mesh_alltoall_t
    private
    type(batch_t)            :: send_buffer
    type(subarray_t)         :: sendpoints
    type(mesh_t),    pointer :: mesh
    integer,         pointer :: requests(:)
    integer                  :: nrequests
    logical                  :: in_progress
    integer,         pointer :: nsend(:)
    integer,         pointer :: nrecv(:)
  end type mesh_alltoall_t
      
contains
  
  ! ---------------------------------------------
  subroutine mesh_alltoall_init(this, mesh, npoints, points, reordered)
    type(mesh_alltoall_t), intent(out)   :: this
    type(mesh_t),  target, intent(in)    :: mesh
    integer,               intent(in)    :: npoints
    integer,               intent(inout) :: points(:)
    logical,               intent(out)   :: reordered

    integer :: npoints_max, ipart, npart, ii, ipoint, ntot
    integer, allocatable :: allpoints(:, :), points_copy(:), sendlist(:), sender(:), pos(:)
    
    this%mesh => mesh
    this%in_progress = .false.

    !prepare the send
    npart = this%mesh%mpi_grp%size

    call MPI_Allreduce(npoints, npoints_max, 1, MPI_INTEGER, MPI_MAX, this%mesh%mpi_grp%comm, mpi_err)

    SAFE_ALLOCATE(points_copy(1:npoints_max))
    points_copy(1:npoints) = points(1:npoints)
    points_copy(npoints + 1:npoints_max) = 0

    SAFE_ALLOCATE(allpoints(1:npoints_max, 1:npart))

    call MPI_Allgather(points_copy(1), npoints_max, MPI_INTEGER, allpoints(1, 1), npoints_max, MPI_INTEGER, &
      this%mesh%mpi_grp%comm, mpi_err)
    
    SAFE_ALLOCATE(this%nsend(1:npart))
    SAFE_ALLOCATE(sendlist(1:npoints_max*npart))

    ntot = 0
    do ipart = 1, npart
      this%nsend(ipart) = 0
      do ii = 1, npoints_max
        ipoint = allpoints(ii, ipart)

        if(ipoint == 0) exit
        
        if(mesh%vp%part(ipoint) == mesh%vp%partno) then
          INCR(this%nsend(ipart), 1)
          INCR(ntot, 1)
          sendlist(ntot) = vec_global2local(mesh%vp, ipoint, mesh%vp%partno)
          
          ASSERT(sendlist(ntot) /= 0)
        end if

      end do
    end do
    
    SAFE_DEALLOCATE_A(sendlist)
    SAFE_DEALLOCATE_A(allpoints)

    call subarray_init(this%sendpoints, ntot, sendlist)
    
    !prepare the reception
    SAFE_ALLOCATE(sender(1:npoints))
    SAFE_ALLOCATE(this%nrecv(1:npart))

    reordered = .false.
    this%nrecv = 0
    do ii = 1, npoints
      sender(ii) = mesh%vp%part(points(ii))
      ! count the points we receive from each node
      INCR(this%nrecv(sender(ii)), 1)
      ! and check if the input list is ordered by partition
      if(ii > 1) then
        if(sender(ii) < sender(ii - 1)) reordered = .true.
      end if
    end do

    if(reordered) then
      ! we need to return a list of points in the order that we will get them
      SAFE_ALLOCATE(pos(1:npart))

      pos(1) = 1
      do ipart = 2, npart
        pos(ipart) = pos(ipart - 1) + this%nrecv(ipart - 1)
      end do

      ! we insert each value in the correct position (this is in fact
      ! an O(N) sort algorithm)
      do ii = 1, npoints
        ipart = sender(ii)
        points(pos(ipart)) = points_copy(ii)
        INCR(pos(ipart), 1)
      end do

      SAFE_DEALLOCATE_A(pos)
    end if

    SAFE_DEALLOCATE_A(points_copy)
    SAFE_DEALLOCATE_A(sender)
    
  end subroutine mesh_alltoall_init
  
  ! ---------------------------------------------

  subroutine mesh_alltoall_end(this)
    type(mesh_alltoall_t), intent(inout) :: this

    nullify(this%mesh)
    ASSERT(.not. this%in_progress)
    
    SAFE_DEALLOCATE_P(this%nsend)
    SAFE_DEALLOCATE_P(this%nrecv)

    call subarray_end(this%sendpoints)
    
  end subroutine mesh_alltoall_end
  
  ! ---------------------------------------------

#include "undef.F90"
#include "complex.F90"
#include "mesh_alltoall_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "mesh_alltoall_inc.F90"

#endif
end module mesh_alltoall_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
