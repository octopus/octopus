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

#include "config_F90.h"

module phonons
  use liboct
  use system
  use hamiltonian
  use scf

  implicit none
  
  private
  public :: phonons_run

  type phonons_type
    integer :: dim
    real(r8), pointer :: DM(:,:), freq(:)

    real(r8) :: disp
  end type phonons_type

contains

  subroutine phonons_run(scf, sys, h)
    type(system_type), intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(scf_type), intent(inout) :: scf
    
    type(phonons_type) :: ph
    integer :: i, j, iunit

    ! create directory for output
    call oct_mkdir('phonons')

    ph%dim = sys%natoms*3
    allocate(ph%DM(ph%dim, ph%dim), ph%freq(ph%dim))

    call oct_parse_double("Displacement", 0.01_r8/units_inp%length%factor, ph%disp)
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
    real(r8) :: energ
    real(r8), allocatable :: forces(:,:), forces0(:,:)

#if defined(HAVE_LAPACK)
    real(r8), allocatable :: work(:)
    integer :: lwork, info
#endif

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
        call generate_external_pot(h, sys)
        call R_FUNC(calcdens) (sys%st, sys%m%np, sys%st%rho)
        call R_FUNC(hamiltonian_setup) (h, sys%m, sys%st, sys)
        call hamiltonian_energy (h, sys%st, sys%eii, -1)
        call scf_run(scf, sys, h)
        do j = 1, sys%natoms
          forces0(j, :) = sys%atom(j)%f(:)
        end do

        sys%atom(i)%x(alpha) = sys%atom(i)%x(alpha) - 2._r8*ph%disp

        ! second force
        call generate_external_pot(h, sys)
        call R_FUNC(calcdens) (sys%st, sys%m%np, sys%st%rho)
        call R_FUNC(hamiltonian_setup) (h, sys%m, sys%st, sys)
        call hamiltonian_energy(h, sys%st, sys%eii, -1)
        call scf_run(scf, sys, h)
        do j = 1, sys%natoms
          forces(j, :) = sys%atom(j)%f(:)
        end do

        sys%atom(i)%x(alpha) = sys%atom(i)%x(alpha) + ph%disp

        do j = 1, sys%natoms
          do beta = 1, 3
            ph%DM(3*(i-1) + alpha, 3*(j-1) + beta) = &
                 (forces0(j, beta) - forces(j, beta)) / (2._r8*ph%disp &
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
#if defined(HAVE_LAPACK)
    lwork = 3*ph%dim - 1
    allocate(work(lwork))
    call dsyev('v', 'u', ph%dim, ph%DM(1, 1), ph%dim, ph%freq(1), work(1), lwork, info)

    if(info.ne.0) then
      write(message(1),'(a,i5)') 'LAPACK "dsygv" returned error code ', info
      call write_fatal(1)
    endif
    deallocate(work)
#endif
    
  end subroutine get_DM

end module phonons
