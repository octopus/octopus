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

module static_pol
use global
use units
use lib_oct_parser
use lib_oct
use io
use mesh_function
use mesh
use system
use hamiltonian
use states
use geometry
use restart
use scf

implicit none

private
public :: static_pol_run

contains

integer function static_pol_run(sys, h, fromScratch) result(ierr)
  type(system_type),      intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  logical,                intent(inout) :: fromScratch

  type(scf_type)             :: scfv
  type(mesh_type),   pointer :: m    ! shortcuts
  type(states_type), pointer :: st

  integer :: iunit, ios, i_start, i, j, is, k, err
  FLOAT :: e_field
  FLOAT, allocatable :: Vpsl_save(:), trrho(:), dipole(:, :, :)
  logical :: resume, out_pol
 
  ierr = 0
  call init_()
 
  ! load wave-functions
  call X(restart_read) ("tmp/restart_gs", sys%st, sys%m, err)
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

  ! Allocate the dipole...
  allocate(dipole(conf%dim, conf%dim, 2))
  dipole = M_ZERO

  if(.not.fromScratch) then
    ! Finds out wether there is an existent restart.pol file.
    inquire(file='tmp/restart.pol', exist=resume)

    if(resume) then
      ! Finds out how many dipoles have already been written.
      call io_assign(iunit)
      open(unit = iunit, file='tmp/restart.pol', status='old', action='read')
      rewind(iunit)
      i_start = 1
      do i = 1, 3
        read(iunit, fmt=*, iostat = ios) ((dipole(i, j, k), j = 1, conf%dim), k = 1, 2)
        if(ios.ne.0) exit
        i_start = i_start + 1
      end do
      call io_close(iunit)
    else
      fromScratch = .true.
    end if
  end if

  if(fromScratch) then
    call io_assign(iunit)
    open(unit = iunit, file='tmp/restart.pol', status='unknown', action='write')
    call io_close(iunit)
    i_start = 1
  end if

  ! Save local pseudopotential
  allocate(Vpsl_save(m%np))
  Vpsl_save = h%ep%Vpsl

  ! Allocate the trrho to the contain the trace of the density.
  allocate(trrho(m%np))
  trrho = M_ZERO
  
  call scf_init(scfv, sys%m, sys%st, h)
  do i = i_start, conf%dim
    do k = 1, 2
      write(message(1), '(a)')
      write(message(2), '(a,i1,a,i1)')'Info: Calculating dipole moment for field ', i, ', #',k
      call write_info(2)
      
      h%ep%vpsl = vpsl_save + (-1)**k*m%x(:,i)*e_field
      
      call scf_run(scfv, m, sys%f_der, st, sys%geo, h, sys%outp)
      
      trrho = M_ZERO
      do is = 1, st%d%spin_channels
        trrho(:) = trrho(:) + st%rho(:, is)
      end do
      
      ! calculate dipole
      do j = 1, conf%dim
        dipole(i, j, k) = dmf_moment(m, trrho, j, 1)
      end do
      
    end do

    ! Writes down the dipole to the file
    call io_assign(iunit)
    open(unit=iunit, file='tmp/restart.pol', action='write', status='old', position='append')
    write(iunit, fmt='(6e20.12)') ((dipole(i, j, k), j = 1, conf%dim), k = 1, 2)
    call io_close(iunit)
    
    if(i == conf%dim) then 
      out_pol = .true.
    end if
  end do
  call scf_end(scfv)

  ! Closes the restart file
  call io_close(iunit)

  if(out_pol) then ! output pol file
    call loct_mkdir("linear")
    call io_assign(iunit)
    open(iunit, file='linear/polarizability', status='unknown')
    write(iunit, '(2a)', advance='no') '# Static polarizability tensor [', &
       trim(units_out%length%abbrev)
    if(conf%dim.ne.1) write(iunit, '(a,i1)', advance='no') '^', conf%dim
    write(iunit, '(a)') ']'
         
    do j = 1, conf%dim
      write(iunit, '(3f12.6)') (dipole(j, 1:conf%dim, 1) - dipole(j, 1:conf%dim, 2))/(M_TWO*e_field) &
         / units_out%length%factor**conf%dim
    end do
    call io_close(iunit)
  end if

  deallocate(Vpsl_save, trrho, dipole)
  
  call end_()

contains

  subroutine init_()
    call push_sub('static_pol_run')

    ! allocate wfs
    allocate(sys%st%X(psi)(sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))

    ! shortcuts
    m  => sys%m
    st => sys%st
 
    ! read in e_field value
    call loct_parse_float('POLStaticField', CNST(0.01)/units_inp%energy%factor*units_inp%length%factor, e_field)
    e_field = e_field * units_inp%energy%factor / units_inp%length%factor
    if (e_field <= M_ZERO) then
      write(message(1), '(a,e14.6,a)') "Input: '", e_field, "' is not a valid POLStaticField"
      message(2) = '(0 < POLStaticField)'
      call write_fatal(2)
    end if

  end subroutine init_

  subroutine end_()
    deallocate(sys%st%X(psi))
    
    call pop_sub()
  end subroutine end_

end function static_pol_run

end module static_pol
