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

#include "global.h"

module modelmb_exchange_syms_oct_m
  use batch_oct_m
  use global_oct_m
  use grid_oct_m
  use index_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use modelmb_particles_oct_m
  use permutations_oct_m
  use profiling_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use young_oct_m

  implicit none

  private

  public :: &
    modelmb_sym_all_states, &
    dmodelmb_sym_state, &
    zmodelmb_sym_state, &
    dmodelmb_sym_all_states, &
    zmodelmb_sym_all_states

contains

#include "real.F90"
#include "modelmb_exchange_syms_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "modelmb_exchange_syms_inc.F90"
#include "undef.F90"

subroutine modelmb_sym_all_states (gr, st)
  type(states_elec_t),    intent(inout) :: st
  type(grid_t),           intent(in)    :: gr

  PUSH_SUB(modelmb_sym_all_states)

  if (st%parallel_in_states) then
    call messages_not_implemented("Model MB parallel in states")
  end if

  if (states_are_complex(st)) then
    call zmodelmb_sym_all_states (gr, st)
  else
    call dmodelmb_sym_all_states (gr, st)
  end if

  POP_SUB(modelmb_sym_all_states)
end subroutine modelmb_sym_all_states

end module modelmb_exchange_syms_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
