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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: states.F90 5022 2009-03-03 17:47:58Z nitsche $

#include "global.h"

module modelmb_exchange_syms_m

  use batch_m
  use datasets_m
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
  use modelmb_1part_m
  use mpi_m
  use mpi_lib_m
  use parser_m
  use permutations_m
  use profiling_m
  use states_m
  use young_m

  implicit none

  private

  public :: dmodelmb_sym_state, zmodelmb_sym_state

contains

#include "real.F90"
#include "modelmb_exchange_syms_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "modelmb_exchange_syms_inc.F90"
#include "undef.F90"

end module modelmb_exchange_syms_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
