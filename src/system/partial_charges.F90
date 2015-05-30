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

module partial_charges_m
  use batch_m
  use batch_ops_m
  use comm_m
  use cube_m
  use cube_function_m
  use derivatives_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use hamiltonian_base_m
  use hirshfeld_m
  use io_m
  use io_function_m
  use lalg_basic_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
  use poisson_m
  use profiling_m
  use projector_m
  use ps_m
  use restart_m
  use simul_box_m
  use species_m
  use splines_m
  use states_m
  use states_dim_m
  use submesh_m
  use symmetries_m
  use symmetrizer_m
  use symm_op_m
  use unit_m
  use unit_system_m
  use varinfo_m

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

end module partial_charges_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
