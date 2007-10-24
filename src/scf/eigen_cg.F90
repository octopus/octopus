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
!! $Id$

#include "global.h"

module eigen_cg_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use lalg_adv_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use mpi_lib_m
  use preconditioners_m
  use profiling_m
  use states_m
  use states_block_m
  use units_m

  implicit none

  private
  public ::                  &
    deigen_solver_cg2,       &
    zeigen_solver_cg2,       &
    deigen_solver_cg2_new,   &
    zeigen_solver_cg2_new,   &
    deigen_solver_lobpcg,    &
    zeigen_solver_lobpcg
  
  contains

#include "real.F90"
#include "eigen_cg_inc.F90"
#include "eigen_lobpcg_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "eigen_cg_inc.F90"
#include "eigen_lobpcg_inc.F90"
end module eigen_cg_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
