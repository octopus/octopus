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

  !> This is defined even when running serial
  type mpi_grp_t
    integer :: comm !< copy of the mpi communicator
    integer :: size !< size of comm (defined also in serial mode)
    integer :: rank !< rank of comm (defined also in serial mode)
  end type mpi_grp_t
 
  type(mpi_grp_t), public :: mpi_world

contains
  ! ---------------------------------------------------------
  subroutine mpi_mod_init()
#if defined(HAVE_MPI)
    integer :: mpi_err
#ifdef HAVE_OPENMP
    integer :: provided
#endif
#ifdef HAVE_SCALAPACK
    integer :: iam, nprocs
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
    call blacs_pinfo(iam, nprocs)

    ! If machine needs additional set up, do it now
    if( nprocs < 1 ) then
      call blacs_setup(iam, mpi_world%size)
    end if

#endif
  end subroutine mpi_mod_init


  ! ---------------------------------------------------------
  subroutine mpi_mod_end()
#if defined(HAVE_MPI)
    integer :: mpi_err

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

end module mpi_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
