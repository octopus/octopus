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
  use messages_m
  use multicomm_m

  public ::                            &
       calc_mode,                      &
       calc_mode_is,                   &
       calc_mode_set,                  &
       calc_mode_parallel_mask,        &
       calc_mode_default_parallel_mask

  integer :: calc_mode_id

  integer, public, parameter ::  &
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

    subroutine calc_mode_set(mode)
      integer, intent(in) :: mode

      call push_sub('calc_mode.calc_mode_set')
      calc_mode_id = mode

      call pop_sub()
    end subroutine calc_mode_set

    integer function calc_mode()

      call push_sub('calc_mode.calc_mode')
      calc_mode = calc_mode_id

      call pop_sub()
    end function calc_mode

    logical function calc_mode_is(mode)
      integer, intent(in) :: mode
      
      call push_sub('calc_mode.calc_mode_is')
      calc_mode_is = (calc_mode_id == mode)

      call pop_sub()
    end function calc_mode_is

    integer function calc_mode_parallel_mask() result(par_mask)
      call push_sub('calc_mode.calc_mode_parallel_mask')

      par_mask = 0

      par_mask = ibset(par_mask, P_STRATEGY_DOMAINS - 1) ! all modes are parallel in domains
      par_mask = ibset(par_mask, P_STRATEGY_KPOINTS - 1)

      select case(calc_mode_id)
      case(CM_TD, CM_GS)
        par_mask = ibset(par_mask, P_STRATEGY_STATES - 1)
      case(CM_CASIDA)
        par_mask = ibset(par_mask, P_STRATEGY_OTHER - 1)
      end select

      call pop_sub()
    end function calc_mode_parallel_mask

    ! This function returns the default modes used for a calculation,
    ! that might be different from the modes available.
    integer function calc_mode_default_parallel_mask() result(par_mask)
      call push_sub('calc_mode.calc_mode_default_parallel_mask')

      par_mask = 0

      par_mask = ibset(par_mask, P_STRATEGY_DOMAINS - 1) ! all modes are parallel in domains

      select case(calc_mode_id)
      case(CM_TD)
        par_mask = ibset(par_mask, P_STRATEGY_STATES - 1)
      case(CM_CASIDA)
        par_mask = ibset(par_mask, P_STRATEGY_OTHER - 1)
      end select

      call pop_sub()
    end function calc_mode_default_parallel_mask

end module calc_mode_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
