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

module unocc
use liboct
use io
use states
use scf

implicit none

type unocc_type
  integer  :: max_iter ! maximum number of iterations
  real(r8) :: conv     ! convergence criterium for the eigenvectors

  type(states_type), pointer :: st
end type unocc_type

contains

subroutine unocc_init(u, m, st)
  type(unocc_type), intent(out) :: u
  type(mesh_type), intent(IN)   :: m
  type(states_type), intent(IN) :: st

  sub_name = 'unocc_init'; call push_sub()

  call oct_parse_int(C_string("UnoccMaximumIter"), 200, u%max_iter)
  call oct_parse_double(C_string("UnoccConv"), 1e-4_r8, u%conv)
  if(u%max_iter <= 0 .and. u%conv <= 0.0_r8) then
    message(1) = "Input: Not all occ convergence criteria can be <= 0"
    message(2) = "Please set one of the following:"
    message(3) = "UnoccMaximumIter | UnoccConv"
    call write_fatal(3)
  end if
  if(u%max_iter <= 0) u%max_iter = huge(u%max_iter)

  ! allocate states structure
  allocate(u%st)
  call states_dup(st, u%st)

  call oct_parse_int(C_string("UnoccNumberStates"), 5, u%st%nst)
  if(u%st%nst <= 0) then
    message(1) = "Input: UnoccNumberStates must be > 0"
    call write_fatal(1)
  end if

  ! setup variables
  u%st%nst = u%st%nst + st%nst
  u%st%st_end = u%st%nst
  nullify(u%st%R_FUNC(psi), u%st%eigenval, u%st%occ)
  allocate(u%st%R_FUNC(psi) (0:m%np, st%dim, u%st%nst, st%nik))
  allocate(u%st%eigenval(u%st%nst, st%nik), u%st%occ(u%st%nst, st%nik))
  u%st%eigenval = 0._r8
  u%st%occ      = 0._r8

  call pop_sub()
end subroutine unocc_init

subroutine unocc_end(u)
  type(unocc_type), intent(out) :: u
  
  if(associated(u%st)) then
    if(associated(u%st%R_FUNC(psi))) then
      deallocate(u%st%R_FUNC(psi), u%st%eigenval, u%st%occ)
      nullify   (u%st%R_FUNC(psi), u%st%eigenval, u%st%occ)
    end if
    deallocate(u%st)
    nullify   (u%st)
  end if

end subroutine unocc_end

subroutine unocc_run(u, scf, sys, h)
  type(unocc_type), intent(inout) :: u
  type(system_type), intent(inout) :: sys
  type(scf_type), intent(inout) :: scf
  type(hamiltonian_type), intent(inout) :: h

  integer :: iter, is, j, iunit
  logical :: finish, file_exists
  real(r8) :: tol, diff(u%st%nst, u%st%nik)
  type(states_type), pointer :: tmp_st

  sub_name = 'unocc_scf'; call push_sub()

  ! let us first set the sys%st structure
  tmp_st => sys%st
  sys%st => u%st

  iter_loop: do iter = 1, u%max_iter
    if(clean_stop()) exit

    call eigen_solver_run(scf%eigens, sys, h, iter, diff)
    tol = maxval(diff)

    ! compute eigenvalues
    call R_FUNC(hamiltonian_eigenval) (h, sys, 1, sys%st%nst) ! eigenvalues

    write(message(1), '(a,i5,a,f12.6)') 'Info: iter = ', iter, ' tol = ', tol
    call write_info(1)
    !call states_write_eigenvalues(stdout, sys%st%nst, sys%st, diff)
    
    ! save restart information
    call R_FUNC(states_write_restart)("restart.occ", sys%m, sys%st)

    finish = (u%conv > 0) .and. (tol <= u%conv)
    if(finish) then
      write(message(1), '(a, i4, a)') &
           'Info: Occupation analysis SCF converged in ', iter, ' iterations.'
      call write_info(1)
      exit iter_loop
    end if

  end do iter_loop

  ! write output file
  call io_assign(iunit)
  call oct_mkdir(C_string("static"))
  open(iunit, status='unknown', file='eigenvalues')

  if(.not.finish) then
    write(iunit, '(a, i4, a)') &
         'Occupation analysis SCF converged in ', iter, ' iterations.'
  else
    write(iunit,'(a)') 'Occupational analysis did *not* converge!'
  end if
  write(iunit,'(a, e17.6)') 'Criterium = ', u%conv
  write(iunit,'(1x)')
  call states_write_eigenvalues(iunit, sys%st%nst, sys%st, diff)
  call io_close(iunit)

  ! output wave-functions
  call states_output(sys%st, sys%m, "static", sys%outp)

  ! we now put this back
  sys%st => tmp_st

  call pop_sub()
end subroutine unocc_run


end module unocc
