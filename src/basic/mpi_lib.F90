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
  use messages_m
  use mpi_m
  use mpi_debug_m

  implicit none

  private

#if defined(HAVE_MPI)
  public :: &
    lmpi_gen_alltoallv

  interface lmpi_gen_alltoallv
    module procedure dlmpi_gen_alltoallv, zlmpi_gen_alltoallv, ilmpi_gen_alltoallv
  end interface
#endif

contains

#if defined(HAVE_MPI)
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
