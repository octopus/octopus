!! Copyright (C) 2018 N. Tancogne-Dejean 
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

module loewdin_oct_m
  use distributed_oct_m
  use global_oct_m
  use io_oct_m
  use kpoints_oct_m
  use lalg_adv_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use orbitalbasis_oct_m
  use orbitalset_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use types_oct_m  
 
  implicit none

  private

  public ::                            &
        dloewdin_orthogonalize,        &
        zloewdin_orthogonalize,        &
        dloewdin_overlap,              &
        zloewdin_overlap,              &
        dloewdin_info,                 &
        zloewdin_info

  contains

#include "undef.F90"
#include "real.F90"
#include "loewdin_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "loewdin_inc.F90"

end module loewdin_oct_m
