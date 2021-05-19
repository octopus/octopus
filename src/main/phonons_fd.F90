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

module phonons_fd_oct_m
  use density_oct_m
  use energy_calc_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use mesh_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scf_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use electrons_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use v_ks_oct_m
  use vibrations_oct_m

  implicit none

  private
  public :: phonons_run

contains

  ! ---------------------------------------------------------
  subroutine phonons_run(system)
    class(*), intent(inout) :: system

    PUSH_SUB(phonons_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = vib_modes not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call phonons_run_legacy(system)
    end select

    POP_SUB(phonons_run)
  end subroutine phonons_run

  ! ---------------------------------------------------------
  subroutine phonons_run_legacy(sys)
    type(electrons_t),      intent(inout) :: sys

    type(vibrations_t) :: vib
    integer :: ierr
    type(restart_t) :: gs_restart

    PUSH_SUB(phonons_run_legacy)

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    ! Why not? The symmetries are computed only for the unperturbed geometry,
    ! and are not valid when the atoms are displaced.
    ! FIXME: implement instead use of symmetry over dynamical matrix to make things more efficient.
    if(sys%st%symmetrize_density .or. sys%kpoints%use_symmetries) then
      message(1) = "Cannot compute vibrational modes by finite differences when symmetry is being used."
      message(2) = "Set KPointsUseSymmetries = no and SymmetrizeDensity = no, for gs run and this run."
      call messages_fatal(2)
    end if
    
    call init_()

    ! load wavefunctions
    call restart_init(gs_restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
    if(ierr == 0) then
      call states_elec_load(gs_restart, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr)
    end if
    if (ierr /= 0) then
      message(1) = "Unable to read wavefunctions."
      call messages_fatal(1)
    end if
    call restart_end(gs_restart)

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call messages_info(1)
    call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)

    call vibrations_init(vib, sys%ions, "fd", sys%namespace)

    !%Variable Displacement
    !%Type float
    !%Default 0.01 a.u.
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% When calculating phonon properties by finite differences (<tt>CalculationMode = vib_modes, 
    !% ResponseMethod = finite_differences</tt>), 
    !% <tt>Displacement</tt> controls how much the atoms are to be moved in order to calculate the 
    !% dynamical matrix.
    !%End
    call parse_variable(sys%namespace, 'Displacement', CNST(0.01), vib%disp, units_inp%length)

    ! calculate dynamical matrix
    call get_dyn_matrix(sys%gr, sys%namespace, sys%mc, sys%ions, sys%st, sys%ks, sys%hm, sys%outp, vib, &
                        sys%space)

    call vibrations_output(vib)
    
    call vibrations_end(vib)

    call end_()
    POP_SUB(phonons_run_legacy)

  contains

    ! ---------------------------------------------------------
    subroutine init_()

      PUSH_SUB(phonons_run_legacy.init_)
      call states_elec_allocate_wfns(sys%st, sys%gr%mesh)

      POP_SUB(phonons_run_legacy.init_)
    end subroutine init_

    ! ---------------------------------------------------------
    subroutine end_()

      PUSH_SUB(phonons_run_legacy.end_)
      call states_elec_deallocate_wfns(sys%st)

      POP_SUB(phonons_run_legacy.end_)
    end subroutine end_

  end subroutine phonons_run_legacy


  ! ---------------------------------------------------------
  subroutine get_dyn_matrix(gr, namespace, mc, ions, st, ks, hm, outp, vib, space)
    type(grid_t),     target, intent(inout) :: gr
    type(namespace_t),        intent(in)    :: namespace
    type(multicomm_t),        intent(in)    :: mc
    type(ions_t),             intent(inout) :: ions
    type(states_elec_t),      intent(inout) :: st
    type(v_ks_t),             intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(output_t),           intent(in)    :: outp
    type(vibrations_t),       intent(inout) :: vib
    type(space_t),            intent(in)    :: space

    type(scf_t)               :: scf
    type(mesh_t),     pointer :: mesh
    integer :: iatom, jatom, alpha, beta, imat, jmat
    FLOAT, allocatable :: forces(:,:), forces0(:,:)

    PUSH_SUB(get_dyn_matrix)

    mesh => gr%mesh

    call scf_init(scf, namespace, gr, ions, st, mc, hm, ks, space)
    SAFE_ALLOCATE(forces0(1:space%dim, 1:ions%natoms))
    SAFE_ALLOCATE(forces (1:space%dim, 1:ions%natoms))
    forces = M_ZERO
    forces0 = M_ZERO

    ! FIXME: why displace in + and -? Could just do + and take difference from undisplaced.
    
    do iatom = 1, ions%natoms
      do alpha = 1, space%dim
        imat = vibrations_get_index(vib, iatom, alpha)

        write(message(1), '(a,i3,3a)') 'Info: Moving atom ', iatom, ' in the +', index2axis(alpha), '-direction.'
        call messages_info(1)

        ! move atom iatom in direction alpha by dist
        ions%pos(alpha, iatom) = ions%pos(alpha, iatom) + vib%disp

        ! first force
        call hamiltonian_elec_epot_generate(hm, namespace, space, gr, ions, st)
        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval=.true.)
        call energy_calc_total (namespace, space, hm, gr, st)
        call scf_mix_clear(scf)
        call scf_run(scf, namespace, space, mc, gr, ions, st, ks, hm, outp, gs_run=.false., verbosity = VERB_COMPACT)
        forces0 = ions%tot_force

        write(message(1), '(a,i3,3a)') 'Info: Moving atom ', iatom, ' in the -', index2axis(alpha), '-direction.'
        call messages_info(1)

        ions%pos(alpha, iatom) = ions%pos(alpha, iatom) - M_TWO*vib%disp

        ! second force
        call hamiltonian_elec_epot_generate(hm, namespace, space, gr, ions, st)
        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval=.true.)
        call energy_calc_total(namespace, space, hm, gr, st)
        call scf_mix_clear(scf)
        call scf_run(scf, namespace, space, mc, gr, ions, st, ks, hm, outp, gs_run=.false., verbosity = VERB_COMPACT)
        forces = ions%tot_force

        ions%pos(alpha, iatom) = ions%pos(alpha, iatom) + vib%disp

        do jatom = 1, ions%natoms
          do beta = 1, space%dim
            jmat = vibrations_get_index(vib, jatom, beta)
            vib%dyn_matrix(jmat, imat) = (forces0(beta, jatom) - forces(beta, jatom)) / (M_TWO*vib%disp) &
              * vibrations_norm_factor(vib, ions, iatom, jatom)
          end do
        end do
        call vibrations_out_dyn_matrix_row(vib, imat)

      end do
    end do
    SAFE_DEALLOCATE_A(forces0)
    SAFE_DEALLOCATE_A(forces)
    call scf_end(scf)

    call vibrations_symmetrize_dyn_matrix(vib)
    call vibrations_diag_dyn_matrix(vib)

    POP_SUB(get_dyn_matrix)
  end subroutine get_dyn_matrix

end module phonons_fd_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
