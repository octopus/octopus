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

  implicit none

  private
  public ::                       &
       floquet_run_init,          &
       floquet_run

contains
  
  
  
  subroutine floquet_run_init()
    PUSH_SUB(floquet_run_init)

    call calc_mode_par_set_parallelization(P_STRATEGY_OTHER, default = .true.)

    POP_SUB(floquet_run_init)
  end subroutine floquet_run_init
  !------------------------------------------------------------------------

  subroutine floquet_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch
    
    integer :: mode, ierr
    type(td_t)  :: td

    FLOAT :: x
    type(restart_t) :: restart
      
    PUSH_SUB(floquet_run)

    call io_mkdir(FLOQUET_DIR)

    call parse_variable('TDFloquetMode', FLOQUET_NON_INTERACTING, mode) 
    
    if(mode == FLOQUET_NON_INTERACTING .or. mode == FLOQUET_FROZEN_PHONON) then
      call td_init(td, sys, hm)
      
      call states_allocate_wfns(sys%st, sys%gr%mesh, alloc_Left = hm%cmplxscl%space)
!       call scf_init(td%scf, sys%gr, sys%geo, sys%st, hm)

      call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)

      if(ierr == 0) call states_load(restart, sys%st, sys%gr, ierr, label = ": gs")
      if (ierr /= 0) then
        message(1) = 'Unable to read ground-state wavefunctions.'
        call messages_fatal(1)
      end if
      call restart_end(restart)
      
      
!       call density_calc(sys%st, sys%gr, sys%st%rho)
!       call v_ks_calc(sys%ks, hm, sys%st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
!       x = minval(sys%st%eigenval(sys%st%st_start, :))
! #ifdef HAVE_MPI
!       if(sys%st%parallel_in_states) then
!         call MPI_Bcast(x, 1, MPI_FLOAT, 0, sys%st%mpi_grp%comm, mpi_err)
!       end if
! #endif
!       call hamiltonian_span(hm, minval(sys%gr%mesh%spacing(1:sys%gr%mesh%sb%dim)), x)
!       ! initialize Fermi energy
!       call states_fermi(sys%st, sys%gr%mesh)
!       call energy_calc_total(hm, sys%gr, sys%st)
      
      
      
      call floquet_init(sys,hm%F,hm%geo,sys%st%d%dim)
      call floquet_hamiltonians_init(hm ,sys%gr, sys%st, sys)
      call floquet_hamiltonian_solve(hm,sys%gr,sys,sys%st)
      
      call states_deallocate_wfns(sys%st)
      
      
    else
      
      
      
      call td_run(sys, hm, fromScratch)
    end if

    POP_SUB(floquet_run)
  end subroutine floquet_run

  

end module floquet_run_oct_m
