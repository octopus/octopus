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

private
public :: unocc_type, &
          unocc_run

type unocc_type
  integer  :: max_iter ! maximum number of iterations
  FLOAT :: conv     ! convergence criterium for the eigenvectors

  type(states_type), pointer :: st
end type unocc_type

contains

! ---------------------------------------------------------
integer function unocc_run(sys, h, fromScratch) result(ierr)
  type(system_type),      intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  logical,                intent(inout) :: fromScratch

  type(eigen_solver_type) :: eigens
  integer :: iter, max_iter, iunit, err
  FLOAT   :: conv
  logical :: converged

  ierr = 0
  call init_()

  if(.not.fromScratch) then
    ! not all states will be read, that is the reason for the <0 instead of .ne.0
    call X(restart_read) ('tmp/restart_unocc', sys%st, sys%m, err)
    if(err < 0) then
      message(1) = "Could not load tmp/restart_unocc: Starting from scratch"
      call write_warning(1)

      fromScratch = .true.
    end if
  end if

  if(fromScratch) then
    call X(restart_read) ('tmp/restart_gs', sys%st, sys%m, err)
    if(err < 0) then
      message(1) = "Could not load tmp/restart_gs: Starting from scratch"
      call write_warning(1)

      ierr = 1
      call end_()
      return
    end if
    message(1) = "Loaded wave-functions from 'tmp/restart_gs'"
    call write_info(1)
  end if

  ! Setup Hamiltonian
  message(1) = 'Info: Setting up Hamiltonian.'
  call write_info(1)
  
  call X(states_calc_dens)(sys%st, sys%m%np, sys%st%rho)
  call X(h_calc_vhxc)(h, sys%m, sys%f_der, sys%st, calc_eigenval=.true.) ! get potentials
  call hamiltonian_energy(h, sys%st, sys%geo%eii, -1)         ! total energy
  

  message(1) = ''
  message(2) = stars
  message(3) = "Starting calculation of unoccupied states"
  call write_info(3)

  do iter = 1, max_iter
    call eigen_solver_run(eigens, sys%m, sys%f_der, sys%st, h, iter, converged)

    write(message(1), '(a,i4,a,g15.6)') "Iter = ", iter, ";  Error = ", maxval(sqrt(eigens%diff**2))
    call write_info(1)

    ! write restart information.
    call X(restart_write) ('tmp/restart_unocc', sys%st, sys%m, err)
    if(err.ne.0) then
      message(1) = 'Unsuccesfull write of "tmp/restart_unocc"'
      call write_fatal(1)
    end if

    if(maxval(abs(eigens%diff)) < eigens%final_tol) exit
  end do

  message(1) = stars
  message(2) = ''
  call write_info(2)

  ! write output file
  call io_mkdir('static')
  iunit = io_open('static/eigenvalues')

  if(converged) then
    write(iunit,'(a)') 'All unoccupied states converged.'
  else
    write(iunit,'(a)') 'Some of the unoccupied states are not fully converged!'
  end if
  write(iunit,'(a, e17.6)') 'Criterium = ', eigens%final_tol
  write(iunit,'(1x)')
  call states_write_eigenvalues(iunit, sys%st%nst, sys%st, eigens%diff)
  call io_close(iunit)
  
  if (conf%periodic_dim>0 .and. sys%st%d%nik>sys%st%d%nspin) then
    iunit = io_open('static/bands.dat')
    call states_write_bands(iunit, sys%st%nst, sys%st)
    call io_close(iunit)
  end if
  
  ! output wave-functions
  call X(states_output) (sys%st, sys%m, sys%f_der, "static", sys%outp)

  call end_()
contains

  ! ---------------------------------------------------------
  subroutine init_()
    integer :: nus
    FLOAT :: conv

    call push_sub('unocc_run')

    call loct_parse_int("NumberUnoccStates", 5, nus)
    if(nus <= 0) then
      message(1) = "Input: NumberUnoccStates must be > 0"
      call write_fatal(1)
    end if

    call loct_parse_int("MaximumIter", 20, max_iter)
    if(max_iter <= 0) max_iter = huge(max_iter)

    ! fix states: THIS IS NOT OK
    sys%st%nst = sys%st%nst + nus
    sys%st%st_end = sys%st%nst

    deallocate(sys%st%eigenval, sys%st%occ)
    allocate(sys%st%X(psi) (sys%m%np, sys%st%d%dim, sys%st%nst, sys%st%d%nik))
    allocate(sys%st%eigenval(sys%st%nst, sys%st%d%nik), sys%st%occ(sys%st%nst, sys%st%d%nik))
    if(sys%st%d%ispin == SPINORS) then
      allocate(sys%st%mag(sys%st%nst, sys%st%d%nik, 2))
      sys%st%mag = M_ZERO
    end if
    sys%st%eigenval = huge(PRECISION)
    sys%st%occ      = M_ZERO

    ! now the eigen solver stuff
    call eigen_solver_init(eigens, sys%st, sys%m, 50)

  end subroutine init_

  ! ---------------------------------------------------------
  subroutine end_()
    deallocate(sys%st%X(psi))
    call eigen_solver_end(eigens)

    call pop_sub()
  end subroutine end_
  
end function unocc_run

end module unocc
