!! Copyright (C) 2005-2006 Heiko Appel, Florian Lorenzen
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

module blacs_proc_grid_m
  use global_m
  use blacs_m
  use mpi_m
  use messages_m
  use profiling_m

  implicit none

  private

#ifdef HAVE_SCALAPACK
  public ::                      &
    blacs_proc_grid_t,           &
    blacs_proc_grid_nullify,     &
    blacs_proc_grid_init,        &
    blacs_proc_grid_end,         &
    blacs_proc_grid_copy,        &
    blacs_proc_grid_null

  type blacs_proc_grid_t
    integer          :: context       !< The blacs context, -1 is object is null.
    integer          :: nprocs        !< Number of processors.
    integer          :: nprow         !< Number of processors per row.
    integer          :: npcol         !< Number of processors per column.
    integer          :: iam           !< Process indentifier.
    integer          :: myrow         !< The row of the processor in the processor grid.
    integer          :: mycol         !< The column of the processor in the processor grid.
    integer, pointer :: usermap(:, :) !< The index of each processor in the grid.
  end type blacs_proc_grid_t

contains

  ! ----------------------------------------------------

  subroutine blacs_proc_grid_nullify(this)
    type(blacs_proc_grid_t), intent(inout) :: this

    PUSH_SUB(blacs_proc_grid_nullify)

    this%context = -1

    POP_SUB(blacs_proc_grid_nullify)
  end subroutine blacs_proc_grid_nullify

  ! -----------------------------------------------------------------------
  !> Initializes a blacs context from an MPI communicator with
  !! topological information.
  !!
  !! \Warning: For the moment this function only works if mpi_grp holds
  !! all the nodes of mpi_world.
  subroutine blacs_proc_grid_init(this, mpi_grp)
    type(blacs_proc_grid_t), intent(out) :: this
    type(mpi_grp_t),         intent(in)  :: mpi_grp

    integer, parameter :: maxdims = 2
    integer :: dims(1:2), topo, coords(1:2), ix, iy, id, xy(2)
    logical :: periods(1:2)
    integer :: mpi_err
    integer :: comm
    logical :: reorder
    integer, allocatable :: procmap(:)

    PUSH_SUB(blacs_proc_grid_init)

    call MPI_Topo_test(mpi_grp%comm, topo, mpi_err)

    if(topo /= MPI_CART) then
      ! We create a new communicator with Cartesian topology
      dims(1) = mpi_grp%size
      dims(2) = 1
      periods = .false.
      reorder = .false.
      call MPI_Cart_create(mpi_grp%comm, 2, dims(1), periods(1), reorder, comm, mpi_err)
    else
      comm = mpi_grp%comm
    end if

    call blacs_pinfo(this%iam, this%nprocs)

    ! The process ID from ScaLAPACK is not always the
    ! same as MPI, so we need to construct a map.
    SAFE_ALLOCATE(procmap(0:mpi_grp%size - 1))
    call MPI_Allgather(this%iam, 1, MPI_INTEGER, procmap(0), 1, MPI_INTEGER, comm, mpi_err)

    ASSERT(this%nprocs == mpi_grp%size)
    ASSERT(this%iam == procmap(mpi_grp%rank))

    dims = 1
    coords = 0
    
    call MPI_Cart_get(comm, maxdims, dims(1), periods(1), coords(1), mpi_err)

    SAFE_ALLOCATE(this%usermap(1:dims(1), 1:dims(2)))
    
    do ix = 1, dims(1)
      xy(1) = ix - 1
      do iy = 1, dims(2)
        xy(2) = iy - 1
        call MPI_Cart_rank(comm, xy, id, mpi_err)
        this%usermap(ix, iy) = procmap(id)
      end do
    end do

    ! get the default system context
    call blacs_get(-1, what = 0, val = this%context)

    ! now get the context associated with the map
    call blacs_gridmap(this%context, this%usermap(1, 1), dims(1), dims(1), dims(2))

    ! and fill the rest of the structure
    call blacs_gridinfo(this%context, this%nprow, this%npcol, this%myrow, this%mycol)

    !check that Blacs and MPI are consistent
    ASSERT(this%nprow == dims(1))
    ASSERT(this%npcol == dims(2))
    ASSERT(this%myrow == coords(1))
    ASSERT(this%mycol == coords(2))

    if(topo /= MPI_CART) then
      call MPI_Comm_free(comm, mpi_err)
    end if

    SAFE_DEALLOCATE_A(procmap)

    POP_SUB(blacs_proc_grid_init)
  end subroutine blacs_proc_grid_init
  
  ! ----------------------------------------------------

  subroutine blacs_proc_grid_end(this)
    type(blacs_proc_grid_t), intent(inout) :: this

    PUSH_SUB(blacs_proc_grid_end)

    if(this%context /= -1) then
      call blacs_gridexit(this%context)
      SAFE_DEALLOCATE_P(this%usermap)
    end if

    POP_SUB(blacs_proc_grid_end)
  end subroutine blacs_proc_grid_end

  ! ----------------------------------------------------

  subroutine blacs_proc_grid_copy(cin, cout)
    type(blacs_proc_grid_t), intent(in)  :: cin
    type(blacs_proc_grid_t), intent(out) :: cout

    PUSH_SUB(blacs_proc_grid_copy)

    cout%context = cin%context 
    cout%nprocs  = cin%nprocs
    cout%nprow   = cin%nprow
    cout%npcol   = cin%npcol
    cout%iam     = cin%iam 
    cout%myrow   = cin%myrow
    cout%mycol   = cin%mycol
    
    if(cout%context /= -1) then
      ! we have to create a new context
      call blacs_get(-1, what = 0, val = cout%context)
      SAFE_ALLOCATE(cout%usermap(1:cout%nprow, 1:cout%npcol))
      cout%usermap = cin%usermap
      call blacs_gridmap(cout%context, cout%usermap(1, 1), cout%nprow, cout%nprow, cout%npcol)
    end if
    
    POP_SUB(blacs_proc_grid_copy)
  end subroutine blacs_proc_grid_copy

  ! ----------------------------------------------------

  logical pure function blacs_proc_grid_null(this)
    type(blacs_proc_grid_t), intent(in) :: this

    blacs_proc_grid_null = this%context == -1
  end function blacs_proc_grid_null
  

#endif

end module blacs_proc_grid_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
