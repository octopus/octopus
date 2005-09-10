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
use simul_box

implicit none

private
public :: unocc_run

contains

! ---------------------------------------------------------
integer function unocc_run(sys, h, fromScratch) result(ierr)
  type(system_type),      intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  logical,                intent(inout) :: fromScratch

  type(eigen_solver_type) :: eigens
  integer :: iunit, err, ik, p, occupied_states
  R_TYPE, allocatable :: h_psi(:,:)
  logical :: converged, l

  ierr = 0
  occupied_states = sys%st%nst
  call init_(sys%gr%m, sys%st)

  call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, sys%gr%m, err)
  if( (err .ne. 0)  .and.  (err < occupied_states) ) then
     message(1) = "Not all the occupied states could be read from '"//trim(tmpdir)//"restart_gs'"
     message(2) = "I will start a gs calculation."
     call write_warning(2)
     ierr = 1
     call end_()
     return
  else
     message(1) = "Loaded wave-functions from '"//trim(tmpdir)//"restart_gs'"
     call write_info(1)
  endif

  ! Setup Hamiltonian
  message(1) = 'Info: Setting up Hamiltonian.'
  call write_info(1)

  call X(states_calc_dens)(sys%st, sys%gr%m%np, sys%st%rho)
  call X(v_ks_calc)(sys%gr, sys%ks, h, sys%st, calc_eigenval=.true.) ! get potentials
  call hamiltonian_energy(h, sys%st, sys%gr%geo%eii, -1)         ! total energy

  message(1) = "Info: Starting calculation of unoccupied states"
  call write_info(1)

  ! First, get the residues of the occupied states.
  ! These are assumed to be converged; otherwise one should do a SCF calculation.
  allocate(h_psi(sys%gr%m%np, h%d%dim))
  do ik = 1, sys%st%d%nik
     do p = 1, eigens%converged
        call X(Hpsi)(h, sys%gr, sys%st%X(psi)(:,:, p, ik) , h_psi, ik)
        eigens%diff(p, ik) = X(states_residue)(sys%gr%m, sys%st%d%dim, h_psi, sys%st%eigenval(p, ik), &
             sys%st%X(psi)(:, :, p, ik))
     enddo
  enddo
  deallocate(h_psi)

  call eigen_solver_run(eigens, sys%gr, sys%st, h, 1, converged, verbose = .true.)

  ! write restart information.
  call X(restart_write) (trim(tmpdir)//'restart_gs', sys%st, sys%gr, err)
  if(err.ne.0) then
    message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'restart_gs"'
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
  call states_write_eigenvalues(iunit, sys%st%nst, sys%st, sys%gr%sb, eigens%diff)
  call io_close(iunit)

  if(simul_box_is_periodic(sys%gr%sb).and. sys%st%d%nik>sys%st%d%nspin) then
    iunit = io_open('static/bands.dat', action='write')
    call states_write_bands(iunit, sys%st%nst, sys%st, sys%gr%sb)
    call io_close(iunit)
  end if

  !%Variable WriteMatrixElements
  !%Type logical
  !%Section 9 Unoccupied States
  !%Description
  !% If true outputs the following matrix elements:
  !%   <i|T + V_ext|j>
  !%   <ij| 1/|r1-r2| |kl>
  !% in the directory ME
  !%End
  call loct_parse_logical(check_inp('WriteMatrixElements'), .false., l)
  if(l) call write_matrix_elements(sys, h)

  ! output wave-functions
  call X(states_output) (sys%st, sys%gr, "static", sys%outp)

  call end_()
contains


  ! ---------------------------------------------------------
  subroutine init_(m, st)
    type(mesh_type),   intent(in)    :: m
    type(states_type), intent(inout) :: st

    integer :: nus

    call push_sub('unocc.unocc_run')

    call loct_parse_int(check_inp('NumberUnoccStates'), 5, nus)
    if(nus <= 0) then
      message(1) = "Input: NumberUnoccStates must be > 0"
      call write_fatal(1)
    end if

    ! fix states: THIS IS NOT OK
    st%nst    = st%nst + nus
    st%st_end = st%nst

    deallocate(st%eigenval, st%occ)
    allocate(st%X(psi) (m%np, st%d%dim, st%nst, st%d%nik))
    allocate(st%eigenval(st%nst, st%d%nik), st%occ(st%nst, st%d%nik))
    if(st%d%ispin == SPINORS) then
      allocate(st%mag(st%nst, st%d%nik, 2))
      st%mag = M_ZERO
    end if
    st%eigenval = huge(PRECISION)
    st%occ      = M_ZERO

    ! now the eigen solver stuff
    call eigen_solver_init(sys%gr, eigens, st, 50)

    ! Having initial and final tolerance does not make sense in this case:
    eigens%init_tol       = eigens%final_tol
    eigens%final_tol_iter = 2
    eigens%converged      = st%nst - nus

  end subroutine init_


  ! ---------------------------------------------------------
  subroutine end_()
    deallocate(sys%st%X(psi))
    call eigen_solver_end(eigens)

    call pop_sub()
  end subroutine end_

end function unocc_run

! warning: only works for spin-unpolarized and 1 k-point
subroutine write_matrix_elements(sys, h)
  use mesh_function
  use poisson

  type(system_type), target, intent(inout) :: sys
  type(hamiltonian_type),    intent(in)    :: h

  type(states_type), pointer :: st
  type(mesh_type),   pointer :: m

  call io_mkdir("ME")

  st => sys%st
  m  => sys%gr%m

  message(1) = "Computing Matrix Elements"
  call write_info(1)

  message(1) = "  :: one-body"
  call write_info(1)
  call one_body()

  message(1) = "  :: two-body"
  call write_info(1)
  call two_body()

contains
  subroutine one_body()
    integer i, j, iunit
    R_TYPE :: me

    call io_assign(iunit)
    iunit = io_open('ME/1-body', action='write')

    do i = 1, st%nst
      do j = 1, st%nst
        if(j > i) cycle
      
        me = st%eigenval(i,1) - X(mf_integrate) (m, R_CONJ(st%X(psi) (:, 1, i, 1)) * &
            h%Vhxc(:, 1) * st%X(psi) (:, 1, j, 1))

        write(iunit, *) i, j, me
      end do
    end do
    
    call io_close(iunit)

  end subroutine one_body


  subroutine two_body()
    integer i, j, k, l, iunit
    R_TYPE :: me
    R_TYPE, allocatable :: n(:), v(:)

    call io_assign(iunit)
    iunit = io_open('ME/2-body', action='write')

    allocate(n(1:m%np), v(1:m%np))

    do i = 1, st%nst
      do j = 1, st%nst
        if(j > i) cycle

        n(:) = R_CONJ(st%X(psi) (:, 1, i, 1)) * st%X(psi) (:, 1, j, 1)
        call X(poisson_solve) (sys%gr, v, n)

        do k = 1, st%nst
          if(k > i) cycle
          do l = 1, st%nst
            if(l > k) cycle
            if(l > j) cycle

            me = X(mf_integrate) (m, v(:) * &
                st%X(psi) (:, 1, k, 1) * R_CONJ(st%X(psi) (:, 1, l, 1)))

            write(iunit, *) i, j, k, l, me
          end do
        end do
      end do
    end do

    call io_close(iunit)

  end subroutine two_body

end subroutine write_matrix_elements


end module unocc
