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

module mpi_m
#if defined(MPI_MOD)
  use mpi
#endif

  implicit none

  ! I do not make this module private on purpose, so that the symbols defined either in
  ! module mpi, or in mpif.h are exported

  ! some machines do not have a mpi module, but a mpi input file
#if defined(MPI_H)
# include "mpif.h"
#endif

  ! This is defined even when running serial
  type mpi_grp_t
    integer :: comm ! copy of the mpi communicator
    integer :: size ! size of comm (defined also in serial mode)
    integer :: rank ! rank of comm (defined also in serial mode)
  end type mpi_grp_t

  type(mpi_grp_t), public :: mpi_world

contains
  ! ---------------------------------------------------------
  subroutine mpi_mod_init()
#if defined(HAVE_MPI)
    integer :: mpi_err
    
    ! initialize MPI
    call MPI_INIT(mpi_err)
    call mpi_grp_init(mpi_world, MPI_COMM_WORLD)
    call MPI_BARRIER(mpi_world%comm, mpi_err)
#else
    call mpi_grp_init(mpi_world, -1)
#endif
  end subroutine mpi_mod_init


  ! ---------------------------------------------------------
  subroutine mpi_mod_end()
#if defined(HAVE_MPI)
    integer :: mpi_err

    ! end MPI
    call MPI_FINALIZE(mpi_err)
#endif  
  end subroutine mpi_mod_end


  ! ---------------------------------------------------------
  subroutine mpi_grp_init(grp, comm)
    type(mpi_grp_t), intent(out) :: grp   ! information about this MPI group
    integer                         :: comm  ! the communicator that defined the group

#if defined(HAVE_MPI)
    integer :: mpi_err

    if(comm .ne. -1) then
      grp%comm = comm
      call MPI_COMM_RANK(grp%comm, grp%rank, mpi_err)
      call MPI_COMM_SIZE(grp%comm, grp%size, mpi_err)
    else
#endif
      grp%comm = -1
      grp%rank = 0
      grp%size = 1
#if defined(HAVE_MPI)
    end if
#endif
  end subroutine mpi_grp_init


  ! ---------------------------------------------------------
  logical function mpi_grp_is_root(grp)
    type(mpi_grp_t), intent(in) :: grp
    
    mpi_grp_is_root = (grp%rank == 0)
  end function mpi_grp_is_root

end module mpi_m
