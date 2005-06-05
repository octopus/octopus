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

module unocc
use global
use messages
use syslabels
use lib_oct_parser
use lib_oct
use mesh
use states
use system
use restart
use v_ks
use hamiltonian
use eigen_solver
use io

implicit none

private
public :: unocc_type, &
          unocc_run

type unocc_type
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
  integer :: iunit, err, ik, p
  R_TYPE, allocatable :: h_psi(:,:)
  logical :: converged

  ierr = 0
  call init_()

  if(.not.fromScratch) then
    ! not all states will be read, that is the reason for the <0 instead of .ne.0
    call X(restart_read) (trim(tmpdir)//'restart_unocc', sys%st, sys%m, err)
    if(err < 0) then
      message(1) = 'Could not load "'//trim(tmpdir)//'restart_unocc": Starting from scratch'
      call write_warning(1)

      fromScratch = .true.
    end if
  end if

  if(fromScratch) then
    call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, sys%m, err)
    if(err < 0) then
      message(1) = "Could not load '"//trim(tmpdir)//"restart_gs: Starting from scratch'"
      call write_warning(1)

      ierr = 1
      call end_()
      return
    end if
    message(1) = "Loaded wave-functions from '"//trim(tmpdir)//"restart_gs'"
    call write_info(1)
  end if

  ! Setup Hamiltonian
  message(1) = 'Info: Setting up Hamiltonian.'
  call write_info(1)
  
  call X(states_calc_dens)(sys%st, sys%m%np, sys%st%rho)
  call X(h_calc_vhxc)(sys%ks, h, sys%m, sys%f_der, sys%st, calc_eigenval=.true.) ! get potentials
  call hamiltonian_energy(h, sys%st, sys%geo%eii, -1)         ! total energy
  
  message(1) = "Info: Starting calculation of unoccupied states"
  call write_info(1)

  ! First, get the residues of the occupied states.
  ! These are assumed to be converged; otherwise one should do a SCF calculation.
  allocate(h_psi(sys%m%np, h%d%dim))
  do ik = 1, sys%st%d%nik
     do p = 1, eigens%converged
        call X(Hpsi)(h, sys%m, sys%f_der, sys%st%X(psi)(:,:, p, ik) , h_psi, ik)
        eigens%diff(p, ik) = X(states_residue)(sys%m, sys%st%d%dim, h_psi, sys%st%eigenval(p, ik), &
             sys%st%X(psi)(:, :, p, ik))
     enddo
  enddo
  deallocate(h_psi)

  call eigen_solver_run(eigens, sys%m, sys%f_der, sys%st, h, 1, converged, verbose = .true.)

  ! write restart information.
  call X(restart_write) (trim(tmpdir)//'restart_unocc', sys%st, sys%m, err)
  if(err.ne.0) then
    message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'restart_unocc"'
    call write_fatal(1)
  end if

  ! write output file
  call io_mkdir('static')
  iunit = io_open('static/eigenvalues', action='write')

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
    iunit = io_open('static/bands.dat', action='write')
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

    call push_sub('unocc_run')

    call loct_parse_int(check_inp('NumberUnoccStates'), 5, nus)
    if(nus <= 0) then
      message(1) = "Input: NumberUnoccStates must be > 0"
      call write_fatal(1)
    end if

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

    ! Having initial and final tolerance does not make sense in this case:
    eigens%init_tol       = eigens%final_tol
    eigens%final_tol_iter = 2
    eigens%converged      = sys%st%nst - nus

  end subroutine init_

  ! ---------------------------------------------------------
  subroutine end_()
    deallocate(sys%st%X(psi))
    call eigen_solver_end(eigens)

    call pop_sub()
  end subroutine end_
  
end function unocc_run

end module unocc
