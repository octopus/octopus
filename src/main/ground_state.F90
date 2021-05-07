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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module ground_state_oct_m
  use calc_mode_par_oct_m
  use electrons_oct_m
  use electrons_ground_state_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use output_oct_m
  use states_elec_oct_m
  use v_ks_oct_m

  implicit none

  private
  public ::                       &
    ground_state_run_init,        &
    ground_state_run

contains

  subroutine ground_state_run_init()

    PUSH_SUB(ground_state_run_init)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
#ifdef HAVE_SCALAPACK
    call calc_mode_par_set_scalapack_compat()
#endif

    POP_SUB(ground_state_run_init)
  end subroutine ground_state_run_init

  ! ---------------------------------------------------------
  subroutine ground_state_run(system, from_scratch)
    class(*),        intent(inout) :: system
    logical,         intent(inout) :: from_scratch

    PUSH_SUB(ground_state_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = gs not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call ground_state_run_legacy(system, from_scratch)
    end select

    POP_SUB(ground_state_run)
  end subroutine ground_state_run

  subroutine ground_state_run_legacy(electrons, from_scratch)
    class(electrons_t), intent(inout) :: electrons
    logical,            intent(inout) :: from_scratch

    PUSH_SUB(ground_state_run_legacy)

    call electrons_ground_state_run(electrons%namespace, electrons%mc, electrons%gr, electrons%ions, electrons%st, electrons%ks, &
      electrons%hm, electrons%outp, electrons%space, from_scratch)

    POP_SUB(ground_state_run_legacy)
  end subroutine ground_state_run_legacy

end module ground_state_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
