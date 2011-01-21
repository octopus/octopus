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

module calc_mode_m
  use global_m
  use messages_m
  use multicomm_m

  public ::                            &
       calc_mode_t,                    &
       calc_mode_init,                 &
       calc_mode_set_parallelization,  &
       calc_mode_is,                   &
       calc_mode_parallel_mask,        &
       calc_mode_default_parallel_mask

  type calc_mode_t
    integer :: mode
    integer :: par_mask
    integer :: def_par_mask
  end type calc_mode_t

  type(calc_mode_t) :: this

  integer, public, parameter ::   &
    CM_NONE               =   0,  &
    CM_GS                 =   1,  &
    CM_UNOCC              =   2,  &
    CM_TD                 =   3,  &
    CM_GEOM_OPT           =   5,  &
    CM_OPT_CONTROL        =   7,  &
    CM_LR_POL             =   8,  &
    CM_CASIDA             =   9,  &
    CM_VDW                =  11,  &
    CM_PHONONS_LR         =  12,  &
    CM_RAMAN              =  13,  &
    CM_ONE_SHOT           =  14,  &
    CM_KDOTP              =  15,  &
    CM_GCM                =  16,  &
    CM_MEMORY             =  17,  &
    CM_INVERTKDS          =  18,  &
    CM_PULPO_A_FEIRA      =  99

  contains

    ! ----------------------------------------------------------

    subroutine calc_mode_init(mode)
      integer, optional, intent(in) :: mode

      ! no push_sub because this routine is called before everything
      ! is fully initialized for the debugging stack

      if(present(mode)) then
        this%mode = mode
      else
        this%mode = CM_NONE
      end if

      this%par_mask = 0
      this%par_mask = ibset(this%par_mask, P_STRATEGY_DOMAINS - 1)
      this%par_mask = ibset(this%par_mask, P_STRATEGY_KPOINTS - 1)

      this%def_par_mask = 0
      this%def_par_mask = ibset(this%def_par_mask, P_STRATEGY_DOMAINS - 1)

    end subroutine calc_mode_init

    ! -----------------------------------------------------

    subroutine calc_mode_end()

    end subroutine calc_mode_end

    ! -----------------------------------------------------

    subroutine calc_mode_set_parallelization(par, default)
      integer, intent(in) :: par
      logical, intent(in) :: default

      this%par_mask = ibset(this%par_mask, par - 1)
      if(default) this%def_par_mask = ibset(this%def_par_mask, par - 1)

    end subroutine calc_mode_set_parallelization

    ! ----------------------------------------------------- 
    !
    ! This function is deprecated. Do not use it. The behaviour of
    ! the code should not depend on the calculation mode.
    logical function calc_mode_is(mode)
      integer, intent(in) :: mode
      
      PUSH_SUB(calc_mode_is)

      calc_mode_is = (this%mode == mode)

      POP_SUB(calc_mode_is)
    end function calc_mode_is

    ! -----------------------------------------------------
    integer function calc_mode_parallel_mask() result(par_mask)
      PUSH_SUB(calc_mode_parallel_mask)

      par_mask = this%par_mask

      POP_SUB(calc_mode_parallel_mask)
    end function calc_mode_parallel_mask

    ! -----------------------------------------------------
    ! This function returns the default modes used for a calculation,
    ! that might be different from the modes available.
    integer function calc_mode_default_parallel_mask() result(par_mask)
      PUSH_SUB(calc_mode_default_parallel_mask)

      par_mask = this%def_par_mask

      POP_SUB(calc_mode_default_parallel_mask)
    end function calc_mode_default_parallel_mask

end module calc_mode_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
