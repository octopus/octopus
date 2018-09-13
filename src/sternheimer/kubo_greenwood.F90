!! Copyright (C) 2018 X. Andrade
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

module kubo_greenwood_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use messages_oct_m
  use system_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                       &
       kubo_greenwood_run

contains
  
  subroutine kubo_greenwood_run(sys, hm)
    type(system_t),       intent(inout) :: sys
    type(hamiltonian_t),  intent(inout) :: hm

    PUSH_SUB(kubo_greewood_run)

    call messages_experimental('Kubo Greenwood')
    
    POP_SUB(kubo_greewood_run)
  end subroutine kubo_greenwood_run

end module kubo_greenwood_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
