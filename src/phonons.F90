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
  use units
  use lib_oct
  use lib_oct_parser
  use io
  use lib_adv_alg
  use mesh
  use functions
  use output
  use external_pot
  use geometry
  use hamiltonian
  use states
  use system
  use restart
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

  integer function phonons_run(sys, h) result(ierr)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h

    type(phonons_type) :: ph
    integer :: i, j, iunit, err

    ierr = 0
    call init_()

    ! load wave-functions
    call X(restart_read) ('tmp/restart_gs', sys%st, sys%m, err)
    if(err.ne.0) then
      message(1) = "Could not load wave-functions: Starting from scratch"
      call write_warning(1)
      
      ierr = 1
      call end_()
      return
    end if
    
    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call X(system_h_setup) (sys, h)

    ! create directory for output
    call io_mkdir('phonons')

    ph%dim = sys%geo%natoms*conf%dim
    allocate(ph%DM(ph%dim, ph%dim), ph%freq(ph%dim))

    call loct_parse_float("Displacement", CNST(0.01)/units_inp%length%factor, ph%disp)
    ph%disp = ph%disp*units_inp%length%factor

    ! calculate dynamical matrix
    call get_DM(sys%m, sys%f_der, sys%st, sys%geo, h, sys%outp, ph)

    ! output phonon frequencies and eigenvectors
    iunit = io_open('phonons/freq')
    do i = 1, ph%dim
      write(iunit, *) i, sqrt(abs(ph%freq(i))) * 219474.63 ! output cm^-1
    end do
    call io_close(iunit)

    ! output phonon eigenvectors
    iunit = io_open('phonons/vec')
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
    subroutine init_()
      call push_sub('phonons_run')

      ! allocate wfs
      allocate(sys%st%X(psi)(sys%m%np, sys%st%d%dim, sys%st%nst, sys%st%d%nik))
    end subroutine init_


    subroutine end_()
      deallocate(sys%st%X(psi))
      
      call pop_sub()
    end subroutine end_    

  end function phonons_run

  subroutine get_DM(m, f_der, st, geo, h, outp, ph)
    type(mesh_type),        intent(IN)    :: m
    type(f_der_type),       intent(inout) :: f_der
    type(states_type),      intent(inout) :: st
    type(geometry_type),    intent(inout) :: geo
    type(hamiltonian_type), intent(inout) :: h
    type(output_type),      intent(IN)    :: outp
    type(phonons_type),     intent(inout) :: ph

    type(scf_type) :: scf
    integer :: i, j, alpha, beta, n, iunit
    FLOAT, allocatable :: forces(:,:), forces0(:,:)

    call scf_init(scf, m, st, h)
    allocate(forces0(geo%natoms, 3), forces(geo%natoms, 3))
    forces = M_ZERO; forces0 = M_ZERO
    n = geo%natoms*conf%dim

    call io_assign(iunit)
    iunit = io_open('phonons/DM')

    do i = 1, geo%natoms
      do alpha = 1, conf%dim
        write(message(1), '(a,i3,a,i2)') 'Info: Moving atom ', i, ' in the direction ', alpha
        call write_info(1)

        ! move atom i in direction alpha by dist
        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) + ph%disp

        ! first force
        call epot_generate(h%ep, m, st, geo, h%reltype)
        call X(calcdens) (st, m%np, st%rho)
        call X(h_calc_vhxc) (h, m, f_der, st, calc_eigenval=.true.)
        call hamiltonian_energy (h, st, geo%eii, -1)
        call scf_run(scf, m, f_der, st, geo, h, outp)
        do j = 1, geo%natoms
          forces0(j, :) = geo%atom(j)%f(:)
        end do

        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) - M_TWO*ph%disp

        ! second force
        call epot_generate(h%ep, m, st, geo, h%reltype)
        call X(calcdens) (st, m%np, st%rho)
        call X(h_calc_vhxc) (h, m, f_der, st, calc_eigenval=.true.)
        call hamiltonian_energy(h, st, geo%eii, -1)
        call scf_run(scf, m, f_der, st, geo, h, outp)
        do j = 1, geo%natoms
          forces(j, :) = geo%atom(j)%f(:)
        end do

        geo%atom(i)%x(alpha) = geo%atom(i)%x(alpha) + ph%disp

        do j = 1, geo%natoms
          do beta = 1, conf%dim
            ph%DM(conf%dim*(i-1) + alpha, conf%dim*(j-1) + beta) = &
                 (forces0(j, beta) - forces(j, beta)) / (M_TWO*ph%disp &
                 * sqrt(geo%atom(i)%spec%weight*geo%atom(j)%spec%weight))
            write(iunit, '(es14.5)', advance='no') ph%DM(conf%dim*(i-1) + alpha, conf%dim*(j-1) + beta)
          end do
        end do
        write(iunit, '(1x)')

      end do
    end do
    deallocate(forces0, forces)
    call scf_end(scf)
    call io_close(iunit)

    ! diagonalize DM
    call lalg_eigensolve(ph%dim, ph%DM, ph%DM, ph%freq)
    
  end subroutine get_DM

end module phonons
