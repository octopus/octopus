!! Copyright (C) 2009 N. Helbig and M. Verstraete
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
!! $Id$

#include "global.h"

module modelmb_exchange_syms_m

  use batch_m
  use geometry_m
  use global_m
  use grid_m
  use hypercube_m
  use index_m
  use io_m
  use lalg_adv_m
  use loct_m
  use math_m
  use mesh_batch_m
  use mesh_function_m
  use messages_m
  use modelmb_particles_m
  use modelmb_density_matrix_m
  use mpi_m
  use mpi_lib_m
  use parser_m
  use permutations_m
  use profiling_m
  use states_m
  use young_m

  implicit none

  private

  public :: modelmb_sym_all_states
  public :: dmodelmb_sym_state, zmodelmb_sym_state
  public :: dmodelmb_sym_all_states, zmodelmb_sym_all_states

contains

#include "real.F90"
#include "modelmb_exchange_syms_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "modelmb_exchange_syms_inc.F90"
#include "undef.F90"

subroutine modelmb_sym_all_states (gr, st, geo)
  type(states_t),         intent(inout) :: st
  type(grid_t),           intent(inout) :: gr
  type(geometry_t),       intent(in)    :: geo

  if (states_are_complex(st)) then
    call zmodelmb_sym_all_states (gr, st, geo)
  else
    call dmodelmb_sym_all_states (gr, st, geo)
  end if
  
end subroutine modelmb_sym_all_states

end module modelmb_exchange_syms_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
