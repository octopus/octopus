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
  use restart_oct_m
  use states_oct_m
  use states_restart_oct_m
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

    type(restart_t) :: gs_restart
    integer :: ierr
    
    PUSH_SUB(kubo_greewood_run)

    call messages_write('Info: Starting Kubo-Greenwood linear-response calculation.')
    call messages_info()
    
    call messages_experimental('Kubo Greenwood')

    call restart_init(gs_restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh = sys%gr%mesh, exact = .true.)
    if(ierr == 0) then
      call states_look_and_load(gs_restart, sys%st, sys%gr)
      call restart_end(gs_restart)
    else
      call messages_write("Cannot find occupied and unoccupied states, a previous ground state calculation is required.")
      call messages_fatal()
    end if
    
    POP_SUB(kubo_greewood_run)
  end subroutine kubo_greenwood_run

end module kubo_greenwood_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
