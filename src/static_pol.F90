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

module static_pol
use liboct
use io
use units
use system
use hamiltonian
use scf

implicit none

contains

subroutine static_pol_run(scf, sys, h)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  type(scf_type), intent(inout) :: scf
  
  integer :: iunit, ios, i_start, i, j, is, k
  real(r8) :: e_field
  real(r8), allocatable :: Vpsl_save(:), trrho(:), dipole(:, :, :)
  logical :: resume, out_pol

  call push_sub('static_pol_run')

  ! read in e_field value
  call oct_parse_double('POLStaticField', 0.01_r8/units_inp%energy%factor*units_inp%length%factor, e_field)
  e_field = e_field * units_inp%energy%factor / units_inp%length%factor
  if (e_field <= 0._r8) then
    write(message(1), '(a,e14.6,a)') "Input: '", e_field, "' is not a valid POLStaticField"
    message(2) = '(0 < POLStaticField)'
    call write_fatal(2)
  end if

  ! Allocate the dipole...
  allocate(dipole(conf%dim, conf%dim, 2))
    dipole = M_ZERO

  ! Finds out wether there is an existent restart.pol file.
  inquire(file='tmp/restart.pol', exist=resume)

  if(resume) then
   ! Finds out how many dipoles have already been written.
   call io_assign(iunit)
   open(unit = iunit, file='tmp/restart.pol', status='old', action='read')
   rewind(iunit)
   i_start = 1
   do i = 1, 3
      read(iunit, fmt=*, iostat = ios) ((dipole(i, j, k), j = 1, 3), k = 1, 2)
      if(ios.ne.0) exit
      i_start = i_start + 1
   enddo
   call io_close(iunit)
  else
   call io_assign(iunit)
   open(unit = iunit, file='tmp/restart.pol', status='new', action='write')
   call io_close(iunit)
   i_start = 1
  endif

  ! Save local pseudopotential
  allocate(Vpsl_save(sys%m%np))
  Vpsl_save = h%Vpsl

  ! Allocate the trrho to the contain the trace of the density.
  allocate(trrho(sys%m%np))
    trrho = M_ZERO

  do i = i_start, conf%dim
     do k = 1, 2
        write(message(1), '(/,a,i1,a,i1)')'Info: Calculating dipole moment for field ', i, ', #',k
        call write_info(1)

        h%vpsl = vpsl_save + (-1)**k*sys%m%lxyz(i, :)*sys%m%h(i)*e_field

        call scf_run(scf, sys, h)

        trrho = M_ZERO
        do is = 1, sys%st%spin_channels
           trrho(:) = trrho(:) + sys%st%rho(:, is)
        enddo

        ! calculate dipole
        do j = 1, conf%dim
           dipole(i, j, k) = R_FUNC(mesh_moment)(sys%m, trrho, j, 1)
        enddo

     enddo

     ! Writes down the dipole to the file
     call io_assign(iunit)
     open(unit=iunit, file='tmp/restart.pol', action='write', status='old', position='append')
     write(iunit, fmt='(6e20.12)') ((dipole(i, j, k), j = 1, 3), k = 1, 2)
     call io_close(iunit)

     if(i == conf%dim) then 
       out_pol = .true.
     end if
  end do

  ! Closes the restart file
  call io_close(iunit)

  if(out_pol) then ! output pol file
    call oct_mkdir("linear")
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

  call pop_sub()
  return
end subroutine static_pol_run

end module static_pol
