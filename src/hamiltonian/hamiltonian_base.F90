!! Copyright (C) 2009 X. Andrade
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
!! $Id: hamiltonian_base.F90 3988 2008-03-31 15:06:50Z fnog $

#include "global.h"

module hamiltonian_base_m
  use batch_m
  use datasets_m
  use derivatives_m
  use geometry_m
  use global_m
  use grid_m
  use hardware_m
  use io_m
  use lalg_basic_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use nl_operator_m
  use parser_m
  use profiling_m
  use projector_m
  use simul_box_m
  use solids_m
  use species_m
  use splines_m
  use states_m
  use states_dim_m
  use submesh_m
  use varinfo_m

  implicit none

  private

  public ::                                    &
    hamiltonian_base_t,                        &
    dhamiltonian_base_local_batch,             &
    zhamiltonian_base_local_batch,             &
    hamiltonian_base_init,                     &
    hamiltonian_base_end

  ! This object stores and applies an electromagnetic potential that
  ! can be represented by different types of potentials.

  type hamiltonian_base_t
    integer                      :: nspin
    type(nl_operator_t), pointer :: kinetic
    FLOAT, pointer               :: potential(:, :)
    type(projector_t),   pointer :: nlproj(:)
  end type hamiltonian_base_t

  integer, public ::                     &
    TERM_ALL                 = HUGE(1),  &
    TERM_KINETIC             =   1,      &
    TERM_LOCAL_POTENTIAL     =   2,      & 
    TERM_NON_LOCAL_POTENTIAL =   4,      &
    TERM_OTHERS              =   8,      &
    TERM_LOCAL_EXTERNAL      =  16

  type(profile_t), save :: prof_vlpsi

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_base_init(this, mesh, nspin)
    type(hamiltonian_base_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(in)    :: nspin

    call push_sub('hamiltonian_base.hamiltonian_base_init')

    this%nspin = nspin
    
    call pop_sub('hamiltonian_base.hamiltonian_base_init')
  end subroutine hamiltonian_base_init

  ! ---------------------------------------------------------
  subroutine hamiltonian_base_end(this)
    type(hamiltonian_base_t), intent(inout) :: this

    call push_sub('hamiltonian_base.hamiltonian_base_end')

    SAFE_DEALLOCATE_P(this%potential)

    call pop_sub('hamiltonian_base.hamiltonian_base_end')
  end subroutine hamiltonian_base_end

#include "undef.F90"
#include "real.F90"
#include "hamiltonian_base_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hamiltonian_base_inc.F90"

end module hamiltonian_base_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
