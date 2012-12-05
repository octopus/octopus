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

module ground_state_m
  use calc_mode_m
  use datasets_m
  use energy_calc_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lcao_m
  use mesh_m
  use messages_m
  use mpi_m
  use multicomm_m
  use parser_m
  use rdmft_m
  use restart_m
  use scf_m
  use simul_box_m
  use species_m
  use species_pot_m
  use states_m
  use states_calc_m
  use states_io_m
  use system_m
  use v_ks_m
  use varinfo_m

  implicit none

  private
  public ::                 &
    ground_state_run_init,  &
    ground_state_run

contains

  subroutine ground_state_run_init()

    PUSH_SUB(ground_state_run_init)

#ifdef HAVE_SCALAPACK
    call calc_mode_set_parallelization(P_STRATEGY_STATES, default = .false.)
#endif
    call calc_mode_set_scalapack_compat()
    call calc_mode_set_extra_states(.true.)

    POP_SUB(ground_state_run_init)
  end subroutine ground_state_run_init

  ! ---------------------------------------------------------
  subroutine ground_state_run(sys, hm, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    logical,             intent(inout) :: fromScratch

    type(scf_t)  :: scfv
    integer      :: ierr

    PUSH_SUB(ground_state_run)

    call messages_write('Info: Allocating ground state wave-functions')
    call messages_info()

    call states_allocate_wfns(sys%st, sys%gr%mesh, alloc_zphi = sys%st%open_boundaries)

#ifdef HAVE_MPI
    ! sometimes a deadlock can occur here (if some nodes can allocate and other cannot)
    if(sys%st%dom_st_kpt_mpi_grp%comm > 0) call MPI_Barrier(sys%st%dom_st_kpt_mpi_grp%comm, mpi_err)
#endif
    call messages_write('      done.')
    call messages_info()

    ! Read free states for ground-state open-boundary calculation.
    if(sys%st%open_boundaries) then
      call read_free_states(sys%st, sys%gr)
      ! allocate self_energy and calculate
      call states_init_self_energy(sys%st, sys%gr, hm%d%nspin, hm%d%ispin, hm%lead)
    end if

    if(.not. fromScratch) then
      ! load wavefunctions
      if (sys%ks%theory_level == RDMFT) then ! in RDMFT we need the exact number of states from either HF or DFT
        call restart_read(trim(restart_dir)//GS_DIR, sys%st, sys%gr, ierr, rdmft = .true.)
      else
        call restart_read(trim(restart_dir)//GS_DIR, sys%st, sys%gr, ierr)
      endif      

      if(ierr .ne. 0) then
        call messages_write("Could not load wavefunctions from '"//trim(restart_dir)//GS_DIR//"'")
        call messages_new_line()
        call messages_write("Starting from scratch!")
        call messages_warning()
        fromScratch = .true.
      end if
    end if

    call scf_init(scfv, sys%gr, sys%geo, sys%st, hm)

    if(fromScratch) then
      if(sys%ks%theory_level == RDMFT) then
        call messages_write("RDMFT calculations cannot be started FromScratch")
        call messages_new_line()
        call messages_write("Run a DFT or HF calculation first")
        call messages_fatal()
      else
        call lcao_run(sys, hm, lmm_r = scfv%lmm_r)
      endif
    else
      ! setup Hamiltonian
      call messages_write('Info: Setting up Hamiltonian.')
      call messages_info()
      call system_h_setup(sys, hm, calc_eigenval = .false.)
    end if

    ! run self-consistency
    if (states_are_real(sys%st)) then
      call messages_write('Info: SCF using real wavefunctions.')
    else
      call messages_write('Info: SCF using complex wavefunctions.')
    end if
    call messages_info()

    ! self-consistency for occupation numbers in RDMFT
    if(sys%ks%theory_level == RDMFT) then 
      call scf_occ(sys%gr, hm, sys%st, sys)
    endif

    if(sys%st%d%pack_states) call states_pack(sys%st)

    call scf_run(scfv, sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp)
    call scf_end(scfv)

    if(sys%st%d%pack_states) call states_unpack(sys%st)

    ! clean up
    call states_deallocate_wfns(sys%st)

    POP_SUB(ground_state_run)
  end subroutine ground_state_run

end module ground_state_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
