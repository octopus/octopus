!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module phonons_m
  use global_m
  use messages_m
  use datasets_m
  use units_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use lib_adv_alg_m
  use mesh_m
  use functions_m
  use output_m
  use external_pot_m
  use geometry_m
  use v_ks_m
  use hamiltonian_m
  use states_m
  use system_m
  use restart_m
  use scf_m
  use grid_m

  implicit none

  private
  public :: phonons_run

  type phonons_t
    integer :: dim
    FLOAT, pointer :: DM(:,:), freq(:)

    FLOAT :: disp
  end type phonons_t

contains

  ! ---------------------------------------------------------
  subroutine phonons_run(sys, h)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h

    type(phonons_t) :: ph
    integer :: i, j, iunit, ierr

    call init_()

    ! load wave-functions
    call restart_read(trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not load wave-functions: Starting from scratch"
      call write_warning(1)
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call system_h_setup(sys, h)

    ! create directory for output
    call io_mkdir('phonons')

    ph%dim = sys%gr%geo%natoms*sys%gr%sb%dim
    ALLOCATE(ph%DM(ph%dim, ph%dim), ph%dim*ph%dim)
    ALLOCATE(ph%freq(ph%dim), ph%dim)

    !%Variable Displacement
    !%Type float
    !%Default 0.01 a.u.
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% When calculating phonon properties by finite differences (<tt>CalculationMode = phonons</tt>)
    !% <tt>Displacement</tt> controls how much the atoms are to be moved in order to calculate the 
    !% dynamical matrix.
    !%End
    call loct_parse_float(check_inp('Displacement'), CNST(0.01)/units_inp%length%factor, ph%disp)
    ph%disp = ph%disp*units_inp%length%factor

    ! calculate dynamical matrix
    call get_DM(sys%gr, sys%st, sys%ks, h, sys%outp, ph)

    ! output phonon frequencies and eigenvectors
    iunit = io_open('phonons/freq', action='write')
    do i = 1, ph%dim
      write(iunit, *) i, sqrt(abs(ph%freq(i))) * 219474.63 ! output cm^-1
    end do
    call io_close(iunit)

    ! output phonon eigenvectors
    iunit = io_open('phonons/vec', action='write')
    do i = 1, ph%dim
      write(iunit, '(i6)', advance='no') i
      do j = 1, ph%dim
        write(iunit, '(es14.5)', advance='no') ph%DM(j, i)
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)

    deallocate(ph%DM, ph%freq)
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()

      call push_sub('phonons.phonons_run')
      call states_allocate_wfns(sys%st, sys%gr%m)

    end subroutine init_

    ! ---------------------------------------------------------
    subroutine end_()
      call states_deallocate_wfns(sys%st)

      call pop_sub()
    end subroutine end_

  end subroutine phonons_run


  ! ---------------------------------------------------------
  subroutine get_DM(gr, st, ks, h, outp, ph)
    type(grid_t), target, intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    type(v_ks_t),         intent(inout) :: ks
    type(hamiltonian_t),  intent(inout) :: h
    type(output_t),       intent(in)    :: outp
    type(phonons_t),      intent(inout) :: ph

    type(scf_t)               :: scf
    type(mesh_t),     pointer :: m
    type(geometry_t), pointer :: geo

    integer :: i, j, alpha, beta, n, iunit
    FLOAT, allocatable :: forces(:,:), forces0(:,:), tmpDM(:,:)

    m   => gr%m
    geo => gr%geo

    call scf_init(gr, scf, st, h)
    ALLOCATE(forces0(geo%natoms, 3), geo%natoms*3)
    ALLOCATE(forces (geo%natoms, 3), geo%natoms*3)
    forces = M_ZERO; forces0 = M_ZERO
    n = geo%natoms*NDIM

    call io_assign(iunit)
    iunit = io_open('phonons/DM', action='write')

    do i = 1, geo%natoms
      do alpha = 1, NDIM
        write(message(1), '(a,i3,a,i2)') 'Info: Moving atom ', i, ' in the direction ', alpha
        call write_info(1)

        ! move atom i in direction alpha by dist
        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) + ph%disp

        ! first force
        call epot_generate(h%ep, gr, st, h%reltype)
        call states_calc_dens(st, m%np, st%rho)
        call v_ks_calc(gr, ks, h, st, calc_eigenval=.true.)
        call hamiltonian_energy (h, gr, st, -1)
        call scf_run(scf, gr, st, ks, h, outp)
        do j = 1, geo%natoms
          forces0(j, :) = geo%atom(j)%f(:)
        end do

        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) - M_TWO*ph%disp

        ! second force
        call epot_generate(h%ep, gr, st, h%reltype)
        call states_calc_dens(st, m%np, st%rho)
        call v_ks_calc(gr, ks, h, st, calc_eigenval=.true.)
        call hamiltonian_energy(h, gr, st, -1)
        call scf_run(scf, gr, st, ks, h, outp)
        do j = 1, geo%natoms
          forces(j, :) = geo%atom(j)%f(:)
        end do

        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) + ph%disp

        do j = 1, geo%natoms
          do beta = 1, NDIM
            ph%DM(NDIM*(i-1) + alpha, NDIM*(j-1) + beta) = &
              (forces0(j, beta) - forces(j, beta)) / (M_TWO*ph%disp &
              * sqrt(geo%atom(i)%spec%weight*geo%atom(j)%spec%weight))
            write(iunit, '(es14.5)', advance='no') ph%DM(NDIM*(i-1) + alpha, NDIM*(j-1) + beta)
          end do
        end do
        write(iunit, '(1x)')

      end do
    end do
    deallocate(forces0, forces)
    call scf_end(scf)
    call io_close(iunit)

    !we need a temporary copy of DM, to avoid passing the same array twice
    ALLOCATE(tmpDM(ph%dim, ph%dim), ph%dim*ph%dim)
    
    tmpDM(1:ph%dim,1:ph%dim)=ph%DM(1:ph%dim,1:ph%dim)

    ! diagonalize DM
    call lalg_eigensolve(ph%dim, tmpDM, ph%DM, ph%freq)

    deallocate(tmpDM)
  end subroutine get_DM

end module phonons_m
