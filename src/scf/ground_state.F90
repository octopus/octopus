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
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_function_oct_m
  use lcao_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use pcm_oct_m
  use rdmft_oct_m
  use restart_oct_m
  use scf_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use electrons_oct_m
  use v_ks_oct_m

  implicit none

  private
  public ::                       &
    ground_state_run_init,        &
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
  subroutine ground_state_run(sys, fromScratch)
    type(electrons_t),   intent(inout) :: sys
    logical,             intent(inout) :: fromScratch

    type(scf_t)     :: scfv
    type(restart_t) :: restart_load, restart_dump
    integer         :: ierr
    type(rdm_t)     :: rdm

    PUSH_SUB(ground_state_run)

    call messages_write('Info: Allocating ground state wave-functions')
    call messages_info()

    if(sys%st%parallel_in_states) then
      call messages_experimental('State parallelization for ground state calculations')
    end if

    if (sys%hm%pcm%run_pcm) then
      if (sys%hm%pcm%epsilon_infty /= sys%hm%pcm%epsilon_0 .and. sys%hm%pcm%tdlevel /= PCM_TD_EQ) then
        message(1) = 'Non-equilbrium PCM is not active in a time-independent run.'
        message(2) = 'You set epsilon_infty /= epsilon_0, but epsilon_infty is not relevant for CalculationMode = gs.'
        message(3) = 'By definition, the ground state is in equilibrium with the solvent.'
        message(4) = 'Therefore, the only relevant dielectric constant is the static one.'
        message(5) = 'Nevertheless, the dynamical PCM response matrix is evaluated for benchamarking purposes.'
        call messages_warning(5)
      end if
    end if

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh, packed=.true.)

#ifdef HAVE_MPI
    ! sometimes a deadlock can occur here (if some nodes can allocate and other cannot)
    if(sys%st%dom_st_kpt_mpi_grp%comm > 0) call MPI_Barrier(sys%st%dom_st_kpt_mpi_grp%comm, mpi_err)
#endif
    call messages_write('Info: Ground-state allocation done.')
    call messages_info()

    if(.not. fromScratch) then
      ! load wavefunctions
      ! in RDMFT we need the full ground state
      call restart_init(restart_load, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, &
                        mesh=sys%gr%mesh, exact = (sys%ks%theory_level == RDMFT))
      if(ierr == 0) &
        call states_elec_load(restart_load, sys%namespace, sys%st, sys%gr, ierr)

      if(ierr /= 0) then
        call messages_write("Unable to read wavefunctions.")
        call messages_new_line()
        call messages_write("Starting from scratch!")
        call messages_warning()
        fromScratch = .true.
      end if
    end if

    call write_canonicalized_xyz_file("exec", "initial_coordinates", sys%geo, sys%gr%mesh, sys%namespace)

    if(sys%ks%theory_level /= RDMFT) then
      call scf_init(scfv, sys%namespace, sys%gr, sys%geo, sys%st, sys%mc, sys%hm)
    end if

    if (fromScratch .and. sys%ks%theory_level /= RDMFT) then
      call lcao_run(sys, lmm_r = scfv%lmm_r)
    else
      ! setup Hamiltonian
      call messages_write('Info: Setting up Hamiltonian.')
      call messages_info()
      call v_ks_h_setup(sys%namespace, sys%gr, sys%geo, sys%st, sys%ks, sys%hm, calc_eigenval = .false., calc_current = .false.)
    end if

    call restart_init(restart_dump, sys%namespace, RESTART_GS, RESTART_TYPE_DUMP, sys%mc, ierr, mesh=sys%gr%mesh)

    ! run self-consistency
    call scf_state_info(sys%st)

    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call sys%st%pack()
    
    ! self-consistency for occupation numbers and natural orbitals in RDMFT
    if(sys%ks%theory_level == RDMFT) then 
      call rdmft_init(rdm, sys%namespace, sys%gr, sys%st, sys%geo, sys%mc, fromScratch)
      call scf_rdmft(rdm, sys%namespace, sys%gr, sys%geo, sys%st, sys%ks, sys%hm, sys%outp, restart_dump)
      call rdmft_end(rdm, sys%gr)
    else
      if(.not. fromScratch) then
        call scf_run(scfv, sys%namespace, sys%mc, sys%gr, sys%geo, sys%st, sys%ks, sys%hm, sys%outp, &
                     restart_load=restart_load, restart_dump=restart_dump)
        call restart_end(restart_load)
      else
        call scf_run(scfv, sys%namespace, sys%mc, sys%gr, sys%geo, sys%st, sys%ks, sys%hm, sys%outp, &
          restart_dump=restart_dump)
      end if

      call scf_end(scfv)
    end if

    call restart_end(restart_dump)

    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call sys%st%unpack()

    ! clean up
    call states_elec_deallocate_wfns(sys%st)

    POP_SUB(ground_state_run)
  end subroutine ground_state_run

end module ground_state_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
