!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: eigen.F90 3030 2007-06-25 16:45:05Z marques $

! This module contains some common usage patterns of MPI routines.

#include "global.h"

module mpi_lib_m
  use global_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use profiling_m

  implicit none

  private

#if !defined(HAVE_MPI)
  integer, public :: mpi_lib_dummy !< this avoids compilers complaining about empty module
#else
  public ::              &
    lmpi_gen_allgatherv, &
    lmpi_translate_rank

  interface lmpi_gen_allgatherv
    module procedure dlmpi_gen_allgatherv, zlmpi_gen_allgatherv, ilmpi_gen_allgatherv
  end interface
#endif

contains

#if defined(HAVE_MPI)
  ! ---------------------------------------------------------
  !> Returns the rank number of the node rank in from_comm for the
  !! to_comm communicator.
  integer function lmpi_translate_rank(from_comm, to_comm, rank)
    integer, intent(in) :: from_comm
    integer, intent(in) :: to_comm
    integer, intent(in) :: rank

    integer :: from_group, to_group, from_rank(1), to_rank(1)

    PUSH_SUB(lmpi_translate_rank)

    call MPI_Comm_group(from_comm, from_group, mpi_err)
    call MPI_Comm_group(to_comm, to_group, mpi_err)

    from_rank(1) = rank
    call MPI_Group_translate_ranks(from_group, 1, from_rank, to_group, to_rank, mpi_err)

    lmpi_translate_rank = to_rank(1)

    POP_SUB(lmpi_translate_rank)
  end function lmpi_translate_rank


#include "undef.F90"
#include "real.F90"
#include "mpi_lib_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mpi_lib_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "mpi_lib_inc.F90"
#else
  subroutine this_module_is_not_empty()
    integer :: neither_is_this_subroutine
    neither_is_this_subroutine = 0
  end subroutine this_module_is_not_empty
#endif
end module mpi_lib_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
