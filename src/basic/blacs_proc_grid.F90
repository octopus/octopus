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
    blacs_proc_grid_init,        &
    blacs_proc_grid_end,         &
    blacs_proc_grid_copy

  type blacs_proc_grid_t
    integer :: context !< blacs context
    integer :: nprocs !< number of processors
    integer :: nprow !< number of processors per row
    integer :: npcol !< number of processors per column
    integer :: iam !< process indentifier
    integer :: myrow !< the row of the processor in the processor grid 
    integer :: mycol !< the column of the processor in the processor grid
  end type blacs_proc_grid_t

contains

  ! -----------------------------------------------------------------------
  !> Initializes a blacs context from an MPI communicator with
  !! topological information.
  !!
  !! Warning: For the moment this function only works if mpi_grp holds
  !! all the nodes of mpi_world.
  subroutine blacs_proc_grid_init(this, mpi_grp)
    type(blacs_proc_grid_t), intent(out) :: this
    type(mpi_grp_t),         intent(in)  :: mpi_grp

    integer, parameter :: maxdims = 2
    integer :: dims(1:2), topo, coords(1:2), ix, iy
    logical :: periods(1:2)
    integer, allocatable :: usermap(:, :)
    integer :: mpi_err
    integer :: comm
    logical :: reorder

    PUSH_SUB(blacs_proc_grid_init)

    call MPI_Topo_test(mpi_grp%comm, topo, mpi_err)

    if(topo /= MPI_CART) then
      ! We create a new communicator with cartesian topology
      dims(1) = mpi_grp%size
      dims(2) = 1
      periods = .false.
      reorder = .false.
      call MPI_Cart_create(mpi_grp%comm, 2, dims(1), periods(1), reorder, comm, mpi_err)
    else
      comm = mpi_grp%comm
    end if

    dims = 1
    coords = 0
    
    call MPI_Cart_get(comm, maxdims, dims(1), periods(1), coords(1), mpi_err)

    SAFE_ALLOCATE(usermap(1:dims(1), 1:dims(2)))
    
    do ix = 1, dims(1)
      do iy = 1, dims(2)
        call MPI_Cart_rank(comm, (/ix - 1, iy - 1/), usermap(ix, iy), mpi_err)
      end do
    end do

    ! get the default system context
    call blacs_get(-1, what = 0, val = this%context)

    ! now get the context associated with the map
    call blacs_gridmap(this%context, usermap(1, 1), dims(1), dims(1), dims(2))

    ! and fill the rest of the structure
    this%nprocs = mpi_grp%size
    this%nprow = dims(1)
    this%npcol = dims(2)
    this%iam = mpi_grp%rank
    this%myrow = coords(1) + 1
    this%mycol = coords(2) + 1

    SAFE_DEALLOCATE_A(usermap)

    if(topo /= MPI_CART) then
      call MPI_Comm_free(comm, mpi_err)
    end if

    POP_SUB(blacs_proc_grid_init)
  end subroutine blacs_proc_grid_init
  
  ! ----------------------------------------------------

  subroutine blacs_proc_grid_end(this)
    type(blacs_proc_grid_t), intent(inout) :: this

    PUSH_SUB(blacs_proc_grid_end)

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

    POP_SUB(blacs_proc_grid_copy)
  end subroutine blacs_proc_grid_copy

  ! ----------------------------------------------------

#endif

end module blacs_proc_grid_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
