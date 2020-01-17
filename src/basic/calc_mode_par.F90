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

!> This module uses a module-scope global object to allow calculation modes
!! to set the available parallelization strategies and whether the layout
!! must be compatible with ScaLAPACK, and to allow this information to be
!! accessed elsewhere. It does not, and should not, contain the definitions
!! of the calculation modes themselves, to avoid writing code explicitly
!! dependent on the calculation mode elsewhere.
module calc_mode_par_oct_m
  use global_oct_m
  use messages_oct_m
  use multicomm_oct_m

  implicit none

  private
  public ::                             &
    calc_mode_par_t,                     &
    calc_mode_par_init,                  &
    calc_mode_par_end,                   &
    calc_mode_par_set_parallelization,   &
    calc_mode_par_unset_parallelization, &
    calc_mode_par_parallel_mask,         &
    calc_mode_par_default_parallel_mask, &
    calc_mode_par_set_scalapack_compat,  &
    calc_mode_par_scalapack_compat

  type calc_mode_par_t
    private
    integer :: par_mask
    integer :: def_par_mask
    logical :: scalapack_compat
  end type calc_mode_par_t

  type(calc_mode_par_t) :: this

contains

  ! ----------------------------------------------------------
  !> Set domains and kpoints as possible parallelization strategies.
  !> Set domains as default parallelization strategy.
  subroutine calc_mode_par_init()
    ! no push_sub because this routine is called before everything
    ! is fully initialized for the debugging stack

    this%par_mask = 0
    this%par_mask = ibset(this%par_mask, P_STRATEGY_DOMAINS - 1)
    this%par_mask = ibset(this%par_mask, P_STRATEGY_KPOINTS - 1)

    this%def_par_mask = 0
    this%def_par_mask = ibset(this%def_par_mask, P_STRATEGY_DOMAINS - 1)
    this%def_par_mask = ibset(this%def_par_mask, P_STRATEGY_KPOINTS - 1)

    this%scalapack_compat = .false.
  end subroutine calc_mode_par_init

  ! -----------------------------------------------------

  subroutine calc_mode_par_end()

  end subroutine calc_mode_par_end

  ! -----------------------------------------------------
  !> Add a parallelization strategy to the list of possible ones.
  !> Make it default also if default = .true.
  subroutine calc_mode_par_set_parallelization(par, default)
    integer, intent(in) :: par
    logical, intent(in) :: default

    this%par_mask = ibset(this%par_mask, par - 1)
    if(default) this%def_par_mask = ibset(this%def_par_mask, par - 1)

  end subroutine calc_mode_par_set_parallelization

  ! -----------------------------------------------------
  !> Remove a parallelization strategy from the list of possible ones.
  !> It will also be removed from the default.
  subroutine calc_mode_par_unset_parallelization(par)
    integer, intent(in) :: par

    this%par_mask = ibclr(this%par_mask, par - 1)
    this%def_par_mask = ibclr(this%def_par_mask, par - 1)

  end subroutine calc_mode_par_unset_parallelization

  ! -----------------------------------------------------

  !> Defines that the current run mode requires division of states
  !! and domains to be compatible with scalapack.
  subroutine calc_mode_par_set_scalapack_compat()
    this%scalapack_compat = .true.
  end subroutine calc_mode_par_set_scalapack_compat

  ! -----------------------------------------------------

  !> Whether the current run mode requires divisions compatible with
  !! scalapack.
  logical pure function calc_mode_par_scalapack_compat() result(compat)
    compat = this%scalapack_compat
  end function calc_mode_par_scalapack_compat

  ! -----------------------------------------------------

  integer function calc_mode_par_parallel_mask() result(par_mask)
    PUSH_SUB(calc_mode_par_parallel_mask)

    par_mask = this%par_mask

    POP_SUB(calc_mode_par_parallel_mask)
  end function calc_mode_par_parallel_mask

  ! -----------------------------------------------------
  !> This function returns the default modes used for a calculation,
  !! that might be different from the modes available.
  integer function calc_mode_par_default_parallel_mask() result(par_mask)
    PUSH_SUB(calc_mode_par_default_parallel_mask)

    par_mask = this%def_par_mask

    POP_SUB(calc_mode_par_default_parallel_mask)
  end function calc_mode_par_default_parallel_mask

end module calc_mode_par_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
