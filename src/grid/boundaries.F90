!! Copyright (C) 2005-2010 Florian Lorenzen, Heiko Appel, X. Andrade
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

module boundaries_m
  use batch_m
  use global_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use par_vec_m
  use profiling_m
  use subarray_m

  implicit none
  
  private

#if defined(HAVE_MPI)
  public ::                        &
    pv_handle_batch_t,             &
    dvec_ghost_update,             &
    zvec_ghost_update,             &
    dghost_update_batch_start,     &
    zghost_update_batch_start,     &
    dghost_update_batch_finish,    &
    zghost_update_batch_finish

  integer :: SEND = 1, RECV = 2

  type pv_handle_batch_t
    private
    type(batch_t)        :: ghost_send
    integer,     pointer :: requests(:)
    integer              :: nnb
  end type pv_handle_batch_t

  type(profile_t), save :: prof_start
  type(profile_t), save :: prof_wait
  type(profile_t), save :: prof_update
    
contains

#include "undef.F90"
#include "complex.F90"
#include "boundaries_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "boundaries_inc.F90"

#endif
end module boundaries_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
