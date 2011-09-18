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
!! $Id: comm.F90 3587 2007-11-22 16:43:00Z xavier $

#include "global.h"

module comm_m
  use global_m
  use messages_m
  use mpi_m
  use types_m
  use profiling_m

  implicit none 

  private

  public ::                       &
    comm_allreduce

  interface comm_allreduce
    module procedure dcomm_allreduce_0, zcomm_allreduce_0, icomm_allreduce_0
    module procedure dcomm_allreduce_1, zcomm_allreduce_1, icomm_allreduce_1
    module procedure dcomm_allreduce_2, zcomm_allreduce_2, icomm_allreduce_2
    module procedure dcomm_allreduce_4, zcomm_allreduce_4, icomm_allreduce_4
  end interface

contains
  
#include "undef.F90"
#include "real.F90"
#include "comm_inc.F90"
  
#include "undef.F90"
#include "complex.F90"
#include "comm_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "comm_inc.F90"

end module comm_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
