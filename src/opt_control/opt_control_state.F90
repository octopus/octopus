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
!! $Id: opt_control_state.F90 6827 2010-07-19 16:13:27Z acastro $

#include "global.h"


!> This module holds the "oct_state_t" datatype, which contains all the information
!! about the "state" that is propagated, i.e. for the moment being just the electronic
!! state, that is a state_t datatype. However, it is intended to include the ionic positions
!! and velocities, also, so that one can propagate (and optimize) simultaneously the quantum
!! and classical subsystems.
module opt_control_state_m
  use states_m

  implicit none

  private
  public :: opt_control_state_t,    &
            opt_control_state_init, &
            opt_control_state_end

  !> This is the datatype that contains the objects that are propagated: in principle this
  !! could be both the quantum and the classical subsystems, but for the moment it is only
  !! the quantum susbsytesm. So this data type is merely a wrapper around the states_t data type.
  type opt_control_state_t
    private
    type(states_t) :: psi
  end type opt_control_state_t

  contains

  subroutine opt_control_state_init(ocs, initial_state)
    type(opt_control_state_t), intent(inout) :: ocs
    type(states_t), intent(in)               :: initial_state

    call states_copy(ocs%psi, initial_state)

  end subroutine opt_control_state_init


  subroutine opt_control_state_end(ocs)
    type(opt_control_state_t), intent(inout) :: ocs

    call states_end(ocs%psi)

  end subroutine opt_control_state_end

end module opt_control_state_m
