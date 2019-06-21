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

module floquet_run_oct_m
  use iso_c_binding
  use calc_mode_par_oct_m
  use comm_oct_m
  use density_oct_m
  use distributed_oct_m
  use floquet_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use hamiltonian_base_oct_m
  use output_oct_m
  use io_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
  use loct_math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use multicomm_oct_m
  use parser_oct_m
  use partial_charges_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_dim_oct_m
  use states_io_oct_m
  use states_restart_oct_m
  use system_oct_m
  use td_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use v_ks_oct_m

  implicit none

  private
  public ::                       &
       floquet_run_init,          &
       floquet_run

contains
  
  
  
  subroutine floquet_run_init()
    PUSH_SUB(floquet_run_init)

    call calc_mode_par_set_parallelization(P_STRATEGY_OTHER, default = .false.) ! enabled not by default

    POP_SUB(floquet_run_init)
  end subroutine floquet_run_init
  !------------------------------------------------------------------------

  subroutine floquet_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch
    
    integer :: mode, ierr
    type(td_t)  :: td
    logical :: td_fromScratch

    FLOAT :: x
    type(restart_t) :: restart
    type(states_t)  :: dressed_st
      
    PUSH_SUB(floquet_run)

    call io_mkdir(FLOQUET_DIR)

    call parse_variable('FloquetMode', FLOQUET_NON_INTERACTING, mode) 

    ! Needed to read info about laser, timestep, etc. 
    call td_init(td, sys, hm)    
    
    if(mode == FLOQUET_NON_INTERACTING .or. mode == FLOQUET_FROZEN_PHONON) then
      
      call states_allocate_wfns(sys%st, sys%gr%mesh, wfs_type = TYPE_CMPLX)
      
      call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)

      if(ierr == 0) call states_load(restart, sys%st, sys%gr, ierr, label = ": gs")
      if (ierr /= 0) then
        message(1) = 'Unable to read ground-state wavefunctions.'
        call messages_fatal(1)
      end if
      call restart_end(restart)
      
      call density_calc(sys%st, sys%gr, sys%st%rho)
      call v_ks_calc(sys%ks, hm, sys%st, sys%geo)
                  
      call floquet_init(sys,hm%F,sys%st%d%dim)
      if (hm%F%boson == OPTION__FLOQUETBOSON__CLASSICAL_FLOQUET) call floquet_hamiltonians_init(hm ,sys%gr, sys%st, sys)
      call floquet_hamiltonian_solve(hm,sys%gr,sys,sys%st, fromScratch)
      
    else
      
      call floquet_init(sys,hm%F,sys%st%d%dim)
      call floquet_hamiltonians_init(hm ,sys%gr, sys%st, sys)
      call floquet_load_td_hamiltonians(hm, sys, ierr)

      if(ierr /= 0) then
        write(message(1),'(a)') 'Failed to load td-Hamiltonians.'
        write(message(2),'(a)') 'Begin time propagation to sample them.'
        call messages_warning(2)        
    
        td_fromScratch = .true.
        call td_run(sys, hm, td_fromScratch)    
      end if

      call states_allocate_wfns(sys%st, sys%gr%mesh, wfs_type = TYPE_CMPLX)
      call floquet_hamiltonian_solve(hm,sys%gr,sys,sys%st, fromScratch)      
      
    end if   

    call states_deallocate_wfns(sys%st)

    POP_SUB(floquet_run)
  end subroutine floquet_run

  

end module floquet_run_oct_m
