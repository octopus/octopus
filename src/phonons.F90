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

#include "global.h"

module phonons
  use global
  use lib_oct
  use lib_oct_parser
  use io
  use lib_adv_alg
  use external_pot
  use hamiltonian
  use states
  use system
  use scf

  implicit none
  
  private
  public :: phonons_run

  type phonons_type
    integer :: dim
    FLOAT, pointer :: DM(:,:), freq(:)

    FLOAT :: disp
  end type phonons_type

contains

  subroutine phonons_run(scf, sys, h)
    type(system_type), intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(scf_type), intent(inout) :: scf
    
    type(phonons_type) :: ph
    integer :: i, j, iunit

    ! create directory for output
    call loct_mkdir('phonons')

    ph%dim = sys%natoms*3
    allocate(ph%DM(ph%dim, ph%dim), ph%freq(ph%dim))

    call loct_parse_float("Displacement", CNST(0.01)/units_inp%length%factor, ph%disp)
    ph%disp = ph%disp*units_inp%length%factor

    ! calculate dynamical matrix
    call get_DM(scf, sys, h, ph)

    ! output phonon frequencies and eigenvectors
    call io_assign(iunit)
    open(iunit, file='phonons/freq', status='unknown')
    do i = 1, ph%dim
      write(iunit, *) i, sqrt(abs(ph%freq(i))) * 219474.63 ! output cm^-1
    end do
    call io_close(iunit)

    ! output phonon eigenvectors
    call io_assign(iunit)
    open(iunit, file='phonons/vec', status='unknown')
    do i = 1, ph%dim
      write(iunit, '(i6)', advance='no') i
      do j = 1, ph%dim
        write(iunit, '(es14.5)', advance='no') ph%DM(j, i)
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)

    deallocate(ph%DM, ph%freq)
  end subroutine phonons_run

  subroutine get_DM(scf, sys, h, ph)
    type(system_type), intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(scf_type), intent(inout) :: scf
    type(phonons_type), intent(inout) :: ph

    integer :: i, j, alpha, beta, n, iunit
    FLOAT :: energ
    FLOAT, allocatable :: forces(:,:), forces0(:,:)

    allocate(forces0(sys%natoms, 3), forces(sys%natoms, 3))
    n = sys%natoms*3

    call io_assign(iunit)
    open(iunit, file='phonons/DM', status='unknown')

    do i = 1, sys%natoms
      do alpha = 1, 3
        write(message(1), '(a,i3,a,i2)') 'Info: Moving atom ', i, ' in the direction ', alpha
        call write_info(1)

        ! move atom i in direction alpha by dist
        sys%atom(i)%x(alpha) = sys%atom(i)%x(alpha) + ph%disp

        ! first force
        call epot_generate(h%ep, sys%m, sys, h%Vpsl, h%reltype)
        call X(calcdens) (sys%st, sys%m%np, sys%st%rho)
        call X(h_calc_vhxc) (h, sys%m, sys%st, sys, calc_eigenval=.true.)
        call hamiltonian_energy (h, sys%st, sys%eii, -1)
        call scf_run(scf, sys, h)
        do j = 1, sys%natoms
          forces0(j, :) = sys%atom(j)%f(:)
        end do

        sys%atom(i)%x(alpha) = sys%atom(i)%x(alpha) - M_TWO*ph%disp

        ! second force
        call epot_generate(h%ep, sys%m, sys, h%Vpsl, h%reltype)
        call X(calcdens) (sys%st, sys%m%np, sys%st%rho)
        call X(h_calc_vhxc) (h, sys%m, sys%st, sys, calc_eigenval=.true.)
        call hamiltonian_energy(h, sys%st, sys%eii, -1)
        call scf_run(scf, sys, h)
        do j = 1, sys%natoms
          forces(j, :) = sys%atom(j)%f(:)
        end do

        sys%atom(i)%x(alpha) = sys%atom(i)%x(alpha) + ph%disp

        do j = 1, sys%natoms
          do beta = 1, 3
            ph%DM(3*(i-1) + alpha, 3*(j-1) + beta) = &
                 (forces0(j, beta) - forces(j, beta)) / (M_TWO*ph%disp &
                 * sqrt(sys%atom(i)%spec%weight*sys%atom(j)%spec%weight))
            write(iunit, '(es14.5)', advance='no') ph%DM(3*(i-1) + alpha, 3*(j-1) + beta)
          end do
        end do
        write(iunit, '(1x)')

      end do
    end do
    deallocate(forces0, forces)
    call io_close(iunit)

    ! diagonalize DM
    call lalg_eigensolve(ph%dim, ph%DM, ph%DM, ph%freq)
    
  end subroutine get_DM

end module phonons
