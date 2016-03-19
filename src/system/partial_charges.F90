!! Copyright (C) 2015 X. Andrade
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
!! $Id:$

#include "global.h"

module partial_charges_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use derivatives_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use hamiltonian_base_oct_m
  use hirshfeld_oct_m
  use io_oct_m
  use io_function_oct_m
  use lalg_basic_oct_m
  use logrid_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use ps_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use species_oct_m
  use splines_oct_m
  use states_oct_m
  use states_dim_oct_m
  use submesh_oct_m
  use symmetries_oct_m
  use symmetrizer_oct_m
  use symm_op_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private

  type partial_charges_t
    integer :: dummy
  end type partial_charges_t

  public ::                             &
    partial_charges_t,                            &
    partial_charges_init,                         &
    partial_charges_end,                          &
    partial_charges_calculate

contains

  subroutine partial_charges_init(this)
    type(partial_charges_t), intent(out)   :: this

    PUSH_SUB(partial_charges_init)

    this%dummy = 0
    
    POP_SUB(partial_charges_init)
  end subroutine partial_charges_init

  !----------------------------------------------

  subroutine partial_charges_calculate(this, mesh, st, geo, hirshfeld_charges)
    type(partial_charges_t), intent(in)    :: this
    type(mesh_t),            intent(in)    :: mesh
    type(states_t),          intent(in)    :: st
    type(geometry_t),        intent(in)    :: geo
    FLOAT, optional,         intent(out)   :: hirshfeld_charges(:)

    integer :: iatom
    type(profile_t), save :: prof
    type(hirshfeld_t) :: hirshfeld
    
    PUSH_SUB(partial_charges_calculate)
    call profiling_in(prof, 'PARTIAL_CHARGES')

    if(present(hirshfeld_charges)) then

      call hirshfeld_init(hirshfeld, mesh, geo, st)
      
      do iatom = 1, geo%natoms
        call hirshfeld_charge(hirshfeld, iatom, st%rho, hirshfeld_charges(iatom))
      end do
      
      call hirshfeld_end(hirshfeld)
    end if
    
    call profiling_out(prof)
    POP_SUB(partial_charges_calculate)

  end subroutine partial_charges_calculate

  ! ---------------------------------------------------------

  subroutine partial_charges_end(this)
    type(partial_charges_t), intent(inout) :: this

    PUSH_SUB(partial_charges_end)

    POP_SUB(partial_charges_end)
  end subroutine partial_charges_end

end module partial_charges_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
