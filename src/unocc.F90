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
use hamiltonian
use eigen_solver

implicit none

type unocc_type
  integer  :: max_iter ! maximum number of iterations
  FLOAT :: conv     ! convergence criterium for the eigenvectors

  type(states_type), pointer :: st
end type unocc_type

contains

subroutine unocc_init(u, m, st, val_charge)
  type(unocc_type),  intent(out) :: u
  type(mesh_type),   intent(IN)  :: m
  type(states_type), intent(IN)  :: st
  FLOAT,             intent(in)  :: val_charge

  call push_sub('unocc_init')

  call loct_parse_int("UnoccMaximumIter", 200, u%max_iter)
  call loct_parse_float("UnoccConv", CNST(1e-4), u%conv)
  if(u%max_iter <= 0 .and. u%conv <= M_ZERO) then
    message(1) = "Input: Not all occ convergence criteria can be <= 0"
    message(2) = "Please set one of the following:"
    message(3) = "UnoccMaximumIter | UnoccConv"
    call write_fatal(3)
  end if
  if(u%max_iter <= 0) u%max_iter = huge(u%max_iter)

  ! allocate states structure
  allocate(u%st)
  call states_init(u%st, m, val_charge)

  call loct_parse_int("UnoccNumberStates", 5, u%st%nst)
  if(u%st%nst <= 0) then
    message(1) = "Input: UnoccNumberStates must be > 0"
    call write_fatal(1)
  end if

  ! setup variables
  u%st%nst = u%st%nst + st%nst
  u%st%st_end = u%st%nst
  allocate(u%st%X(psi) (m%np, u%st%dim, u%st%nst, u%st%nik))
  allocate(u%st%eigenval(u%st%nst, u%st%nik), u%st%occ(u%st%nst, u%st%nik))
  if(u%st%d%ispin == SPINORS) then
    allocate(u%st%mag(u%st%nst, u%st%d%nik, 2))
    u%st%mag = M_ZERO
  end if
  u%st%eigenval = M_ZERO
  u%st%occ      = M_ZERO

  call pop_sub()
end subroutine unocc_init

subroutine unocc_end(u)
  type(unocc_type), intent(inout) :: u
  
  if(associated(u%st)) then
    call states_end(u%st)
    deallocate(u%st)
    nullify   (u%st)
  end if

end subroutine unocc_end

subroutine unocc_run(u, m, st, h, outp)
  type(unocc_type),       intent(inout) :: u
  type(mesh_type),        intent(IN)    :: m
  type(states_type),      intent(IN)    :: st
  type(hamiltonian_type), intent(inout) :: h
  type(output_type),      intent(IN)    :: outp

  type(eigen_solver_type) :: eigens
  integer :: iunit
  logical :: converged

  call push_sub('unocc_run')

  ! Initialize eigens (not necessary to call eigen_init)
  eigens%es_type        = RS_CG
  eigens%init_tol       = u%conv 
  eigens%final_tol      = u%conv
  eigens%final_tol_iter = 1
  eigens%es_maxiter     = u%max_iter
  allocate(eigens%diff(u%st%nst, u%st%nik))

  call eigen_solver_run(eigens, m, u%st, h, 1, converged)

  ! write output file
  call io_assign(iunit)
  call loct_mkdir("static")
  open(iunit, status='unknown', file='static/eigenvalues')
  if(converged) then
    write(iunit,'(a)') 'All unoccupied stated converged.'
  else
    write(iunit,'(a)') 'Some of the unoccupied states are not fully converged!'
  end if
  write(iunit,'(a, e17.6)') 'Criterium = ', u%conv
  write(iunit,'(1x)')
  call states_write_eigenvalues(iunit, u%st%nst, u%st, eigens%diff)
  call io_close(iunit)
  
  if (conf%periodic_dim>0 .and. st%nik>st%d%nspin) then
    call io_assign(iunit)
    open(iunit, status='unknown', file='static/bands.dat')
    call states_write_bands(iunit, u%st%nst, u%st)
    call io_close(iunit)
  end if

  ! write restart information.
  call X(states_write_restart)("tmp/restart.occ", m, u%st) 

  ! output wave-functions
  call X(states_output) (u%st, m, "static", outp)

  ! Deallocate eigens..
  deallocate(eigens%diff); nullify(eigens%diff)

  call pop_sub()
end subroutine unocc_run

end module unocc
