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

module mpi_m
#if defined(MPI_MOD)
  use mpi
#endif
  use blacs_m

  implicit none

  ! I do not make this module private on purpose, so that the symbols defined either in
  ! module mpi, or in mpif.h are exported

  ! some machines do not have a mpi module, but a mpi input file
#if defined(MPI_H)
include "mpif.h"
#endif

  ! This is defined even when running serial
  type mpi_grp_t
    integer :: comm !< copy of the mpi communicator
    integer :: size !< size of comm (defined also in serial mode)
    integer :: rank !< rank of comm (defined also in serial mode)
  end type mpi_grp_t
 
  integer, parameter :: BLACS_DLEN = 9

  type blacs_proc_grid_t
    integer :: context !< blacs context
    integer :: nprocs !< number of processors
    integer :: nprow !< number of processors per row
    integer :: npcol !< number of processors per column
    integer :: iam !< process indentifier
    integer :: myrow !< the row of the processor in the processor grid 
    integer :: mycol !< the column of the processor in the processor grid
  end type blacs_proc_grid_t

  type(mpi_grp_t), public :: mpi_world
  type (blacs_proc_grid_t), public :: blacs
contains
  ! ---------------------------------------------------------
  subroutine mpi_mod_init()
#if defined(HAVE_MPI)
    integer :: mpi_err
#ifdef HAVE_OPENMP
    integer :: provided
#endif

    ! initialize MPI
#if defined(HAVE_OPENMP) && defined(HAVE_MPI2)
    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, mpi_err)
#else
    call MPI_Init(mpi_err)
#endif
    call mpi_grp_init(mpi_world, MPI_COMM_WORLD)
    call MPI_Barrier(mpi_world%comm, mpi_err)
#else
    call mpi_grp_init(mpi_world, -1)
#endif

#ifdef HAVE_SCALAPACK

    ! Initialize Blacs to be able to use ScaLAPACK 
    ! Determine my process number and the number of processes in machine
    call BLACS_PINFO(blacs%iam, blacs%nprocs)
    ! If machine needs additional set up, do it now
    if( blacs%nprocs < 1 ) then
    !  if( blacs%iam == 0 )&
     !   blacs%nprocs = blacs%nprow*blacs%npcol 
      call BLACS_SETUP(blacs%iam,mpi_world%size)
    end if  
    
#endif
  end subroutine mpi_mod_init


  ! ---------------------------------------------------------
  subroutine mpi_mod_end()
#if defined(HAVE_MPI)
    integer :: mpi_err
#ifdef HAVE_SCALAPACK 
    ! I think BLACS context is also released whem MPI_Finalize is called
    integer,parameter :: continue = 1
#endif
    ! end MPI
    call MPI_Finalize(mpi_err)
#endif 
 
  end subroutine mpi_mod_end


  ! ---------------------------------------------------------
  subroutine mpi_grp_init(grp, comm)
    type(mpi_grp_t), intent(out)  :: grp   !< information about this MPI group
    integer,         intent(in)   :: comm  !< the communicator that defined the group

#if defined(HAVE_MPI)
    integer :: mpi_err

    if(comm .ne. -1 .and. comm .ne. MPI_COMM_NULL) then
      grp%comm = comm
      call MPI_Comm_rank(grp%comm, grp%rank, mpi_err)
      call MPI_Comm_size(grp%comm, grp%size, mpi_err)
    else
#endif
      grp%comm = -1
      grp%rank = 0
      grp%size = 1
#if defined(HAVE_MPI)
    end if
#endif
  end subroutine mpi_grp_init

  subroutine mpi_grp_copy_equal()
    stop "mpi_grp_copy_equal"
  end subroutine mpi_grp_copy_equal

  ! ---------------------------------------------------------
  subroutine mpi_grp_copy(mpi_grp_out, mpi_grp_in)
    type(mpi_grp_t), intent(out) :: mpi_grp_out
    type(mpi_grp_t), intent(in)  :: mpi_grp_in

    mpi_grp_out%comm = mpi_grp_in%comm
    mpi_grp_out%size = mpi_grp_in%size
    mpi_grp_out%rank = mpi_grp_in%rank
  end subroutine mpi_grp_copy

  ! ---------------------------------------------------------
  logical function mpi_grp_is_root(grp)
    type(mpi_grp_t), intent(in) :: grp
    
    mpi_grp_is_root = (grp%rank == 0)
  end function mpi_grp_is_root
  
  ! ---------------------------------------------------------

  !> Initializes a blacs context from an MPI communicator with
  !> topological information.
  subroutine blacs_proc_grid_from_mpi(this, mpi_grp)
    type(blacs_proc_grid_t), intent(out) :: this
    type(mpi_grp_t),         intent(in)  :: mpi_grp

    integer, parameter :: maxdims = 2
    integer :: dims(1:2), topo, coords(1:2), ix, iy
    logical :: periods(1:2)
    integer, allocatable :: usermap(:, :)
    integer :: mpi_err
    
    call MPI_Topo_test(mpi_grp%comm, topo, mpi_err)

    if(topo /= MPI_CART) then
      ! We cannot use anything more elegant here.
      stop "Not implemented."
    end if

    dims = 1
    coords = 0
    
    call MPI_Cart_get(mpi_grp%comm, maxdims, dims(1), periods(1), coords(1), mpi_err)
    
    !cannot use SAFE ALLOCATE here
    allocate(usermap(1:dims(1), 1:dims(2)))
    
    do ix = 1, dims(1)
      do iy = 1, dims(2)
        call MPI_Cart_rank(mpi_grp%comm, (/ix - 1, iy - 1/), usermap(ix, iy), mpi_err)
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

    deallocate(usermap)

  end subroutine blacs_proc_grid_from_mpi

end module mpi_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
