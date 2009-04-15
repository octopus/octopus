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
  use energy_m
  use external_pot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lalg_adv_m
  use loct_parser_m
  use mesh_m
  use messages_m
  use multicomm_m
  use h_sys_output_m
  use profiling_m
  use vibrations_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use units_m
  use v_ks_m

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

    call init_()

    ! load wave-functions
    call restart_read(trim(restart_dir)//'gs', sys%st, sys%gr, sys%geo, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not load wave-functions: Starting from scratch"
      call write_warning(1)
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call system_h_setup(sys, hm)

    call vibrations_init(vib, sys%geo, sys%gr%sb)

    !%Variable Displacement
    !%Type float
    !%Default 0.01 a.u.
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% When calculating phonon properties by finite differences (<tt>CalculationMode = phonons</tt>)
    !% <tt>Displacement</tt> controls how much the atoms are to be moved in order to calculate the 
    !% dynamical matrix.
    !%End
    call loct_parse_float(datasets_check('Displacement'), CNST(0.01)/units_inp%length%factor, vib%disp)
    vib%disp = vib%disp*units_inp%length%factor

    ! calculate dynamical matrix
    call get_dyn_matrix(sys%gr, sys%geo, sys%mc, sys%st, sys%ks, hm, sys%outp, vib)

    call vibrations_output(vib, "_fd")
    
    call vibrations_end(vib)

    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()

      call push_sub('phonons.phonons_run')
      call states_allocate_wfns(sys%st, sys%gr%mesh)

    end subroutine init_

    ! ---------------------------------------------------------
    subroutine end_()
      call states_deallocate_wfns(sys%st)

      call pop_sub()
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
    type(h_sys_output_t), intent(in)    :: outp
    type(vibrations_t),      intent(inout) :: vib

    type(scf_t)               :: scf
    type(mesh_t),     pointer :: m

    integer :: i, j, alpha, beta
    FLOAT, allocatable :: forces(:,:), forces0(:,:)

    m   => gr%mesh

    call scf_init(scf, gr, geo, st, hm)
    ALLOCATE(forces0(geo%natoms, 3), geo%natoms*3)
    ALLOCATE(forces (geo%natoms, 3), geo%natoms*3)
    forces = M_ZERO
    forces0 = M_ZERO

    do i = 1, geo%natoms
      do alpha = 1, gr%mesh%sb%dim
        write(message(1), '(a,i3,a,i2)') 'Info: Moving atom ', i, ' in the direction ', alpha
        call write_info(1)

        ! move atom i in direction alpha by dist
        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) + vib%disp

        ! first force
        call hamiltonian_epot_generate(hm, gr, geo, st)
        call states_calc_dens(st, m%np)
        call v_ks_calc(gr, ks, hm, st, calc_eigenval=.true.)
        call total_energy (hm, gr, st, -1)
        call scf_run(scf, gr, geo, st, ks, hm, outp, gs_run=.false., verbosity = VERB_COMPACT)
        do j = 1, geo%natoms
          forces0(j, :) = geo%atom(j)%f(:)
        end do

        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) - M_TWO*vib%disp

        ! second force
        call hamiltonian_epot_generate(hm, gr, geo, st)
        call states_calc_dens(st, m%np)
        call v_ks_calc(gr, ks, hm, st, calc_eigenval=.true.)
        call total_energy(hm, gr, st, -1)
        call scf_run(scf, gr, geo, st, ks, hm, outp, gs_run=.false., verbosity = VERB_COMPACT)
        do j = 1, geo%natoms
          forces(j, :) = geo%atom(j)%f(:)
        end do

        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) + vib%disp

        do j = 1, geo%natoms
          do beta = 1, gr%mesh%sb%dim
            vib%dyn_matrix(vibrations_get_index(vib, i, alpha), vibrations_get_index(vib, j, beta)) = &
              (forces0(j, beta) - forces(j, beta)) / (M_TWO*vib%disp )
          end do
        end do

      end do
    end do
    SAFE_DEALLOCATE_A(forces0)
    SAFE_DEALLOCATE_A(forces)
    call scf_end(scf)

    call vibrations_normalize_dyn_matrix(vib, geo)
    call vibrations_diag_dyn_matrix(vib)

  end subroutine get_dyn_matrix

end module phonons_fd_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
