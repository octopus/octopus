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

module unocc
use global
use lib_oct_parser
use lib_oct
use mesh
use states
use system
use restart
use hamiltonian
use eigen_solver

implicit none

type unocc_type
  integer  :: max_iter ! maximum number of iterations
  FLOAT :: conv     ! convergence criterium for the eigenvectors

  type(states_type), pointer :: st
end type unocc_type

contains

logical function unocc_run(sys, h, fromScratch)
  type(system_type),      intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  logical,                intent(inout) :: fromScratch

  type(eigen_solver_type) :: eigens
  integer :: max_iter, iunit, ierr
  FLOAT   :: conv
  logical :: converged

  unocc_run = .true.
  call init_()

  if(.not.fromScratch) then
    call restart_load("tmp/restart_unocc", sys%st, sys%m, ierr)
    if(ierr > 0) then ! Fatal error are flagged by ierr > 0
      message(1) = "Could not load tmp/restart_unocc: Starting from scratch"
      call write_warning(1)

      fromScratch = .true.
    end if
  end if

  if(fromScratch) then
    call restart_load("tmp/restart_gs", sys%st, sys%m, ierr)
    if(ierr > 0) then ! Fatal error are flagged by ierr > 0
      message(1) = "Could not load tmp/restart_gs: Starting from scratch"
      call write_warning(1)

      unocc_run = .false.
      call end_()
      return
    end if
  end if

  ! setup hamiltonian (PROBABLY THIS SHOULD NOT BE DONE)
  call states_fermi(sys%st, sys%m)
  call X(calcdens)(sys%st, sys%m%np, sys%st%rho)
  
  ! setup Hamiltonian
  message(1) = 'Info: Setting up Hamiltonian.'
  call write_info(1)
  
  call X(h_calc_vhxc)(h, sys%m, sys%f_der, sys%st, calc_eigenval=.true.) ! get potentials
  call states_fermi(sys%st, sys%m)                            ! occupations
  call hamiltonian_energy(h, sys%st, sys%geo%eii, -1)         ! total energy
  
  call eigen_solver_run(eigens, sys%m, sys%f_der, sys%st, h, 1, converged)
  
  ! write output file
  call io_assign(iunit)
  call loct_mkdir("static")
  open(iunit, status='unknown', file='static/eigenvalues')
  if(converged) then
    write(iunit,'(a)') 'All unoccupied stated converged.'
  else
    write(iunit,'(a)') 'Some of the unoccupied states are not fully converged!'
  end if
  write(iunit,'(a, e17.6)') 'Criterium = ', eigens%final_tol
  write(iunit,'(1x)')
  call states_write_eigenvalues(iunit, sys%st%nst, sys%st, eigens%diff)
  call io_close(iunit)
  
  if (conf%periodic_dim>0 .and. sys%st%nik>sys%st%d%nspin) then
    call io_assign(iunit)
    open(iunit, status='unknown', file='static/bands.dat')
    call states_write_bands(iunit, sys%st%nst, sys%st)
    call io_close(iunit)
  end if
  
  ! write restart information.
  call restart_write("tmp/restart_unocc", sys%st, sys%m, ierr)
  
  ! output wave-functions
  call X(states_output) (sys%st, sys%m, sys%f_der, "static", sys%outp)

  call end_()
contains

  subroutine init_()
    integer :: max_iter, nus
    FLOAT :: conv

    call push_sub('unocc_run')

    call loct_parse_int("NumberUnoccStates", 5, nus)
    if(nus <= 0) then
      message(1) = "Input: NumberUnoccStates must be > 0"
      call write_fatal(1)
    end if

    ! fix states: THIS IS NOT OK
    sys%st%nst = sys%st%nst + nus
    sys%st%st_end = sys%st%nst

    deallocate(sys%st%eigenval, sys%st%occ)
    allocate(sys%st%X(psi) (sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))
    allocate(sys%st%eigenval(sys%st%nst, sys%st%nik), sys%st%occ(sys%st%nst, sys%st%nik))
    if(sys%st%d%ispin == SPINORS) then
      allocate(sys%st%mag(sys%st%nst, sys%st%d%nik, 2))
      sys%st%mag = M_ZERO
    end if
    sys%st%eigenval = huge(PRECISION)
    sys%st%occ      = M_ZERO

    ! now the eigen solver stuff
    call eigen_solver_init(eigens, sys%st, sys%m, 200)
  end subroutine init_

  subroutine end_()
    deallocate(sys%st%X(psi))
    call eigen_solver_end(eigens)

    call pop_sub()
  end subroutine end_
  
end function unocc_run

end module unocc
