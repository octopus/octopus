!! Copyright (C) 2008 X. Andrade
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

module gauge_field_force_m
  use batch_m
  use batch_ops_m
  use datasets_m
  use derivatives_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use lalg_basic_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
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

  public ::                               &
    gauge_field_get_force

contains

  ! ---------------------------------------------------------
  subroutine gauge_field_get_force(gr, hm, geo, pj, phases, st, force)
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(in)    :: hm
    type(geometry_t),     intent(in)    :: geo
    type(projector_t),    intent(in)    :: pj(:)
    CMPLX,                intent(in)    :: phases(:, :)
    type(states_t),       intent(inout) :: st
    type(gauge_force_t),  intent(out)   :: force

    integer :: idir

    PUSH_SUB(gauge_field_get_force)

    ASSERT(st%d%nspin == 1)

    do idir = 1, gr%sb%dim
      force%vecpot(idir) = CNST(4.0)*M_PI*P_c/gr%sb%rcell_volume*dmf_integrate(gr%mesh, st%current(:, idir, 1))
    end do

    POP_SUB(gauge_field_get_force)
  end subroutine gauge_field_get_force

end module gauge_field_force_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
