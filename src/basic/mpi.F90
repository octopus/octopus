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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module mpi_oct_m
#if defined(MPI_MOD)
  use mpi
#endif
  use blacs_oct_m

  implicit none

  ! I do not make this module private on purpose, so that the symbols defined either in
  ! module mpi, or in mpif.h are exported

  ! some machines do not have a mpi module, but a mpi input file
#if defined(MPI_H)
include "mpif.h"
#endif

  !> This is defined even when running serial
  type mpi_grp_t
    ! Components are public by default
    integer :: comm !< copy of the mpi communicator
    integer :: size !< size of comm (defined also in serial mode)
    integer :: rank !< rank of comm (defined also in serial mode)
  end type mpi_grp_t
 
  type(mpi_grp_t), public :: mpi_world

  !> used to store return values of mpi calls
  integer, public :: mpi_err

contains

  ! ---------------------------------------------------------
  subroutine mpi_mod_init(is_serial)
    logical, intent(in) :: is_serial

#if defined(HAVE_MPI)
#ifdef HAVE_OPENMP
    integer :: provided
#endif
#ifdef HAVE_SCALAPACK
    integer :: iam, nprocs
    integer :: blacs_default_system_context !< for blacs/openmpi bug workaround
#endif
#endif

    if(is_serial) then
      call mpi_grp_init(mpi_world, -1)
      return
    end if

#if defined(HAVE_MPI)
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
    call blacs_pinfo(iam, nprocs)

    ! If machine needs additional set up, do it now
    if( nprocs < 1 ) then
      call blacs_setup(iam, mpi_world%size)
    end if

    !> ensure there was at least one call to blacs_gridinit() or blacs_gridmap()
    !> without it, blacs_exit() triggers an error
    !>
    !> *** An error occurred in MPI_Type_free
    !> *** MPI_ERR_TYPE: invalid datatype
    !>
    !> in openmpi
    call blacs_get(-1,0, blacs_default_system_context)
    call blacs_gridinit(blacs_default_system_context, 'R', 1, 1)

#endif
  end subroutine mpi_mod_init


  ! ---------------------------------------------------------
  subroutine mpi_mod_end()

#ifdef HAVE_SCALAPACK
    if(mpi_world%comm /= -1) call blacs_exit(1)
#endif

#if defined(HAVE_MPI)
    ! end MPI, if we started it
    if(mpi_world%comm /= -1) call MPI_Finalize(mpi_err)
#endif 

  end subroutine mpi_mod_end


  ! ---------------------------------------------------------
  subroutine mpi_grp_init(grp, comm)
    type(mpi_grp_t), intent(out)  :: grp   !< information about this MPI group
    integer,         intent(in)   :: comm  !< the communicator that defined the group

    grp%comm = comm
#if defined(HAVE_MPI)
    if (grp%comm == MPI_COMM_NULL) grp%comm = -1
#endif

    if (grp%comm == -1) then
      grp%rank = 0
      grp%size = 1
#if defined(HAVE_MPI)
    else
      call MPI_Comm_rank(grp%comm, grp%rank, mpi_err)
      call MPI_Comm_size(grp%comm, grp%size, mpi_err)
#endif
    end if

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

end module mpi_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
