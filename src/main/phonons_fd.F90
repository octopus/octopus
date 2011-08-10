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

module phonons_fd_m
  use datasets_m
  use density_m
  use energy_m
  use epot_m
  use geometry_m
  use global_m
  use grid_m
  use output_m
  use hamiltonian_m
  use io_m
  use lalg_adv_m
  use mesh_m
  use messages_m
  use multicomm_m
  use parser_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use unit_m
  use unit_system_m
  use utils_m 
  use v_ks_m
  use vibrations_m

  implicit none

  private
  public :: phonons_run

contains

  ! ---------------------------------------------------------
  subroutine phonons_run(sys, hm)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm

    type(vibrations_t) :: vib
    integer :: ierr

    PUSH_SUB(phonons_run)

    call init_()

    ! load wavefunctions
    call restart_read(trim(restart_dir)//GS_DIR, sys%st, sys%gr, sys%geo, ierr, exact = .true.)

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call messages_info(1)
    call system_h_setup(sys, hm)

    call vibrations_init(vib, sys%geo, sys%gr%sb, "fd")

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
    call parse_float(datasets_check('Displacement'), units_from_atomic(units_inp%length, CNST(0.01)), vib%disp)
    vib%disp = units_to_atomic(units_inp%length, vib%disp)

    ! calculate dynamical matrix
    call get_dyn_matrix(sys%gr, sys%geo, sys%mc, sys%st, sys%ks, hm, sys%outp, vib)

    call vibrations_output(vib)
    
    call vibrations_end(vib)

    call end_()
    POP_SUB(phonons_run)

  contains

    ! ---------------------------------------------------------
    subroutine init_()

      PUSH_SUB(phonons_run.init_)
      call states_allocate_wfns(sys%st, sys%gr%mesh)

      POP_SUB(phonons_run.init_)
    end subroutine init_

    ! ---------------------------------------------------------
    subroutine end_()

      PUSH_SUB(phonons_run.end_)
      call states_deallocate_wfns(sys%st)

      POP_SUB(phonons_run.end_)
    end subroutine end_

  end subroutine phonons_run


  ! ---------------------------------------------------------
  subroutine get_dyn_matrix(gr, geo, mc, st, ks, hm, outp, vib)
    type(grid_t), target, intent(inout) :: gr
    type(geometry_t),     intent(inout) :: geo
    type(multicomm_t),    intent(in)    :: mc
    type(states_t),       intent(inout) :: st
    type(v_ks_t),         intent(inout) :: ks
    type(hamiltonian_t),  intent(inout) :: hm
    type(output_t),       intent(in)    :: outp
    type(vibrations_t),   intent(inout) :: vib

    type(scf_t)               :: scf
    type(mesh_t),     pointer :: mesh
    integer :: iatom, jatom, alpha, beta, imat, jmat
    FLOAT, allocatable :: forces(:,:), forces0(:,:)

    PUSH_SUB(get_dyn_matrix)

    mesh => gr%mesh

    call scf_init(scf, gr, geo, st, hm)
    SAFE_ALLOCATE(forces0(1:geo%natoms, 1:3))
    SAFE_ALLOCATE(forces (1:geo%natoms, 1:3))
    forces = M_ZERO
    forces0 = M_ZERO

    do iatom = 1, geo%natoms
      do alpha = 1, gr%mesh%sb%dim
        write(message(1), '(a,i3,3a)') 'Info: Moving atom ', iatom, ' in the ', index2axis(alpha), '-direction.'
        call messages_info(1)

        imat = vibrations_get_index(vib, iatom, alpha)

        ! move atom iatom in direction alpha by dist
        geo%atom(iatom)%x(alpha) = geo%atom(iatom)%x(alpha) + vib%disp

        ! first force
        call hamiltonian_epot_generate(hm, gr, geo, st)
        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, hm, st, calc_eigenval=.true.)
        call total_energy (hm, gr, st, -1)
        call scf_run(scf, gr, geo, st, ks, hm, outp, gs_run=.false., verbosity = VERB_COMPACT)
        do jatom = 1, geo%natoms
          forces0(jatom, :) = geo%atom(jatom)%f(:)
        end do

        geo%atom(iatom)%x(alpha) = geo%atom(iatom)%x(alpha) - M_TWO*vib%disp

        ! second force
        call hamiltonian_epot_generate(hm, gr, geo, st)
        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, hm, st, calc_eigenval=.true.)
        call total_energy(hm, gr, st, -1)
        call scf_run(scf, gr, geo, st, ks, hm, outp, gs_run=.false., verbosity = VERB_COMPACT)
        do jatom = 1, geo%natoms
          forces(jatom, :) = geo%atom(jatom)%f(:)
        end do

        geo%atom(iatom)%x(alpha) = geo%atom(iatom)%x(alpha) + vib%disp

        do jatom = 1, geo%natoms
          do beta = 1, gr%mesh%sb%dim
            jmat = vibrations_get_index(vib, jatom, beta)
            vib%dyn_matrix(imat, jmat) = &
              (forces0(jatom, beta) - forces(jatom, beta)) / (M_TWO*vib%disp ) &
              * vibrations_norm_factor(vib, geo, iatom, jatom)
            call vibrations_out_dyn_matrix(vib, imat, jmat)
          end do
        end do

      end do
    end do
    SAFE_DEALLOCATE_A(forces0)
    SAFE_DEALLOCATE_A(forces)
    call scf_end(scf)

    call vibrations_diag_dyn_matrix(vib)

    POP_SUB(get_dyn_matrix)
  end subroutine get_dyn_matrix

end module phonons_fd_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
