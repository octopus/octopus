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
  use energy_calc_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use lcao_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use parser_oct_m
  use rdmft_oct_m
  use restart_oct_m
  use scf_oct_m
  use simul_box_oct_m
  use species_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_io_oct_m
  use states_restart_oct_m
  use system_oct_m
  use v_ks_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                 &
    ground_state_run_init,  &
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
  subroutine ground_state_run(sys, hm, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    logical,             intent(inout) :: fromScratch

    type(scf_t)  :: scfv
    type(restart_t) :: restart_load, restart_dump
    integer      :: ierr

    PUSH_SUB(ground_state_run)

    call messages_write('Info: Allocating ground state wave-functions')
    call messages_info()

    if(sys%st%parallel_in_states) then
      call messages_experimental('State parallelization for ground state calculations')
    end if
    
    call states_allocate_wfns(sys%st, sys%gr%mesh)

#ifdef HAVE_MPI
    ! sometimes a deadlock can occur here (if some nodes can allocate and other cannot)
    if(sys%st%dom_st_kpt_mpi_grp%comm > 0) call MPI_Barrier(sys%st%dom_st_kpt_mpi_grp%comm, mpi_err)
#endif
    call messages_write('Info: Ground-state allocation done.')
    call messages_info()

    if(.not. fromScratch) then
      ! load wavefunctions
      ! in RDMFT we need the full ground state
      call restart_init(restart_load, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, &
                        mesh=sys%gr%mesh, exact = (sys%ks%theory_level == RDMFT))
      if(ierr == 0) &
        call states_load(restart_load, sys%st, sys%gr, ierr)

      if(ierr /= 0) then
        call messages_write("Unable to read wavefunctions.")
        call messages_new_line()
        call messages_write("Starting from scratch!")
        call messages_warning()
        fromScratch = .true.
      end if
    end if

    call scf_init(scfv, sys%gr, sys%geo, sys%st, sys%mc, hm)

    if(fromScratch) then
      if(sys%ks%theory_level == RDMFT) then
        call messages_write("RDMFT calculations cannot be started FromScratch")
        call messages_new_line()
        call messages_write("Run a DFT calculation with XCFunctional = oep_x first")
        call messages_fatal()
      else
        call lcao_run(sys, hm, lmm_r = scfv%lmm_r)
      end if
    else
      ! setup Hamiltonian
      call messages_write('Info: Setting up Hamiltonian.')
      call messages_info()
      call system_h_setup(sys, hm, calc_eigenval = .false.)
    end if

    call restart_init(restart_dump, RESTART_GS, RESTART_TYPE_DUMP, sys%mc, ierr, mesh=sys%gr%mesh)

    ! run self-consistency
    if (states_are_real(sys%st)) then
      call messages_write('Info: SCF using real wavefunctions.')
    else
      call messages_write('Info: SCF using complex wavefunctions.')
    end if
    call messages_info()

    if(sys%st%d%pack_states) call states_pack(sys%st)
    
    ! self-consistency for occupation numbers in RDMFT
    if(sys%ks%theory_level == RDMFT) then 
      call scf_rdmft(sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp,scfv%max_iter)
    else
      if(.not. fromScratch) then
        call scf_run(scfv, sys%mc, sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp, &
                     restart_load=restart_load, restart_dump=restart_dump)
        call restart_end(restart_load)
      else
        call scf_run(scfv, sys%mc, sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp, restart_dump=restart_dump)
      end if
    end if

    call scf_end(scfv)
    call restart_end(restart_dump)

    if(sys%st%d%pack_states) call states_unpack(sys%st)

    ! clean up
    call states_deallocate_wfns(sys%st)

    POP_SUB(ground_state_run)
  end subroutine ground_state_run

end module ground_state_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
