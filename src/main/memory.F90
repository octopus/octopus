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
!! $Id: memory.F90 4881 2009-01-20 20:58:05Z dstrubbe $

#include "global.h"

module memory_m
  use global_m
  use hamiltonian_m
  use mesh_m
  use messages_m
  use states_m
  use system_m
  
  implicit none

  private
  public :: memory_run

contains

  ! ---------------------------------------------------------
  subroutine memory_run(sys, hm)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm

    real(8) :: mesh_global, mesh_local, wfns, ground_state

    mesh_global = mesh_global_memory(sys%gr%mesh)/CNST(1024.0)**2
    mesh_local  = mesh_local_memory(sys%gr%mesh)/CNST(1024.0)**2

    write(message(1), '(a)')        "Mesh"
    write(message(2), '(a,f10.1,a)') "  global  ", mesh_global, " [Mb] (global)"
    write(message(3), '(a,f10.1,a)') "  local   ", mesh_local,  " [Mb] (par_domains)"
    call messages_info(3)

    wfns = states_wfns_memory(sys%st, sys%gr%mesh)/CNST(1024.0)**2

    write(message(1), '(a)')        "States"
    write(message(2), '(a,f10.1,a)') "  real    ", wfns,       " [Mb] (par_kpoints + par_states + par_domains)"
    write(message(3), '(a,f10.1,a)') "  complex ", M_TWO*wfns, " [Mb] (par_kpoints + par_states + par_domains)"
    call messages_info(3)

    if(states_are_real(sys%st)) then
      ground_state = mesh_global + mesh_local + wfns
    else
      ground_state = mesh_global + mesh_local + M_TWO*wfns
    end if

    write(message(1), '(a)')        "Total"
    write(message(2), '(a,f10.1,a)') "  gs      ", ground_state, " [Mb]"
    write(message(3), '(a,f10.1,a)') "  td      ", mesh_global + mesh_local + M_TWO*wfns, " [Mb]"
    call messages_info(3)

  end subroutine memory_run

end module memory_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
