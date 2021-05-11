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

module electrons_ground_state_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_function_oct_m
  use ions_oct_m
  use lcao_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use output_oct_m
  use pcm_oct_m
  use rdmft_oct_m
  use restart_oct_m
  use scf_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use v_ks_oct_m

  implicit none

  private
  public ::                       &
    electrons_ground_state_run

contains

  ! ---------------------------------------------------------
  subroutine electrons_ground_state_run(namespace, mc, gr, ions, st, ks, hm, outp, space, fromScratch)
    type(namespace_t),        intent(in)    :: namespace
    type(multicomm_t),        intent(in)    :: mc
    type(grid_t),             intent(inout) :: gr
    type(ions_t),             intent(inout) :: ions
    type(states_elec_t),      intent(inout) :: st
    type(v_ks_t),             intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(output_t),           intent(in)    :: outp
    type(space_t),            intent(in)    :: space
    logical,                  intent(inout) :: fromScratch

    type(scf_t)     :: scfv
    type(restart_t) :: restart_load, restart_dump
    integer         :: ierr
    type(rdm_t)     :: rdm

    PUSH_SUB(ground_state_run_legacy)

    call messages_write('Info: Allocating ground state wave-functions')
    call messages_info()

    if (st%parallel_in_states) then
      call messages_experimental('State parallelization for ground state calculations')
    end if

    if (hm%pcm%run_pcm) then
      if (hm%pcm%epsilon_infty /= hm%pcm%epsilon_0 .and. hm%pcm%tdlevel /= PCM_TD_EQ) then
        message(1) = 'Non-equilbrium PCM is not active in a time-independent run.'
        message(2) = 'You set epsilon_infty /= epsilon_0, but epsilon_infty is not relevant for CalculationMode = gs.'
        message(3) = 'By definition, the ground state is in equilibrium with the solvent.'
        message(4) = 'Therefore, the only relevant dielectric constant is the static one.'
        message(5) = 'Nevertheless, the dynamical PCM response matrix is evaluated for benchamarking purposes.'
        call messages_warning(5)
      end if
    end if

    call states_elec_allocate_wfns(st, gr%mesh, packed=.true.)

#ifdef HAVE_MPI
    ! sometimes a deadlock can occur here (if some nodes can allocate and other cannot)
    if (st%dom_st_kpt_mpi_grp%comm > 0) call MPI_Barrier(st%dom_st_kpt_mpi_grp%comm, mpi_err)
#endif
    call messages_write('Info: Ground-state allocation done.')
    call messages_info()

    if (.not. fromScratch) then
      ! load wavefunctions
      ! in RDMFT we need the full ground state
      call restart_init(restart_load, namespace, RESTART_GS, RESTART_TYPE_LOAD, mc, ierr, mesh=gr%mesh, &
        exact = (ks%theory_level == RDMFT))
      if (ierr == 0) then
        call states_elec_load(restart_load, namespace, space, st, gr%mesh, hm%kpoints, ierr)
      end if

      if (ierr /= 0) then
        call messages_write("Unable to read wavefunctions.")
        call messages_new_line()
        call messages_write("Starting from scratch!")
        call messages_warning()
        fromScratch = .true.
      end if
    end if

    call write_canonicalized_xyz_file("exec", "initial_coordinates", ions, gr%sb, namespace)

    if (ks%theory_level /= RDMFT) then
      call scf_init(scfv, namespace, gr, ions, st, mc, hm, ks, space)
    end if

    if (fromScratch .and. ks%theory_level /= RDMFT) then
      call lcao_run(namespace, space, gr, ions, st, ks, hm, lmm_r = scfv%lmm_r)
    else
      ! setup Hamiltonian
      call messages_write('Info: Setting up Hamiltonian.')
      call messages_info()
      call v_ks_h_setup(namespace, space, gr, ions, st, ks, hm, calc_eigenval = .false., calc_current = .false.)
    end if

    call restart_init(restart_dump, namespace, RESTART_GS, RESTART_TYPE_DUMP, mc, ierr, mesh=gr%mesh)

    ! run self-consistency
    call scf_state_info(st)

    if (st%d%pack_states .and. hamiltonian_elec_apply_packed(hm)) then
      call st%pack()
    end if

    ! self-consistency for occupation numbers and natural orbitals in RDMFT
    if (ks%theory_level == RDMFT) then
      call rdmft_init(rdm, namespace, gr, st, mc, space, fromScratch)
      call scf_rdmft(rdm, namespace, space, gr, ions, st, ks, hm, outp, restart_dump)
      call rdmft_end(rdm)
    else
      if(.not. fromScratch) then
        call scf_run(scfv, namespace, space, mc, gr, ions, st, ks, hm, outp, restart_load=restart_load, restart_dump=restart_dump)
        call restart_end(restart_load)
      else
        call scf_run(scfv, namespace, space, mc, gr, ions, st, ks, hm, outp, restart_dump=restart_dump)
      end if

      call scf_end(scfv)
    end if

    call restart_end(restart_dump)

    if (st%d%pack_states .and. hamiltonian_elec_apply_packed(hm)) then
      call st%unpack()
    end if

    ! clean up
    call states_elec_deallocate_wfns(st)

    POP_SUB(ground_state_run_legacy)
  end subroutine electrons_ground_state_run

end module electrons_ground_state_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
