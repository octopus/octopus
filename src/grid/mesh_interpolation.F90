!! Copyright (C) 2014 X. Andrade
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

#include "global.h"

module mesh_interpolation_oct_m
  use comm_oct_m
  use global_oct_m
  use index_oct_m
  use iso_c_binding
  use loct_math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use par_vec_oct_m
  use simul_box_oct_m
  
  implicit none
  private

  public ::                           &
    mesh_interpolation_t,             &
    mesh_interpolation_init,          &
    mesh_interpolation_end,           &
    mesh_interpolation_evaluate,      &
    dmesh_interpolation_test,         &
    zmesh_interpolation_test

  type mesh_interpolation_t
    private
    
    type(mesh_t), pointer :: mesh
  end type mesh_interpolation_t

  interface mesh_interpolation_evaluate
    module procedure dmesh_interpolation_evaluate
    module procedure zmesh_interpolation_evaluate
    module procedure dmesh_interpolation_evaluate_vec
    module procedure zmesh_interpolation_evaluate_vec
  end interface mesh_interpolation_evaluate

contains

  subroutine mesh_interpolation_init(this, mesh)
    type(mesh_interpolation_t), intent(out)   :: this
    type(mesh_t), target,       intent(in)    :: mesh
    
    PUSH_SUB(mesh_interpolation_init)

    ASSERT(.not. mesh%use_curvilinear)
    
    this%mesh => mesh

    POP_SUB(mesh_interpolation_init)
  end subroutine mesh_interpolation_init
  
  ! ---------------------------------------------------
  subroutine mesh_interpolation_end(this)
    type(mesh_interpolation_t), intent(inout)   :: this

    PUSH_SUB(mesh_interpolation_end)

    nullify(this%mesh)

    POP_SUB(mesh_interpolation_end)  
  end subroutine mesh_interpolation_end


#include "undef.F90"
#include "real.F90"
#include "mesh_interpolation_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_interpolation_inc.F90"

end module mesh_interpolation_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
