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

subroutine nbo2xyz(sampling)
  use global
  use liboct
  use atom

  implicit none

  integer :: sampling

  character(len=80) :: sysname, str, nbofile, xyzfile
  integer :: ierr, natoms, ncatoms, nspecies, i, nbo_unit, xyz_unit, iter, j
  real(r8) :: dump
  type(atom_type), pointer :: atm(:)
  type(atom_classical_type), pointer :: catm(:)
  type(specie_type), pointer :: spec(:)

  write(*,*) 'step1'
  ! Opens and parses the inp file through the liboct
  ierr = oct_parse_init(C_string('inp'), C_string('out.oct'))
  if(ierr .ne. 0) then
    message(1) = "Error initializing liboct"
    call write_fatal(1)
  end if
  
  write(*,*) 'step2'
  ! Finds out if we want to be verbose or not.
  call oct_parse_int(C_string('verbose'), 30, conf%verbose)
  if(conf%verbose >= 999 .and. mpiv%node == 0) then
    message(1) = 'Entering DEBUG mode'
    call write_warning(1)
  end if
  
  ! Sets the dimensionaliy of the problem.
  conf%dim=3

  write(*,*) 'step3'
  ! Fixes the units.
  call units_init()

  ! Gets the system name
  call oct_parse_str('SystemName', 'system', sysname)

  write(*,*)'step4'
  ! Sets the filenames
  nbofile = trim(sysname)//'.nbo'
  xyzfile = trim(sysname)//'-movie.xyz'

  ! how many do we have?
  str = C_string("Species")
  nspecies = oct_parse_block_n(str)
  if (nspecies < 1) then
    message(1) = "Input: Species block not specified"
    message(2) = '% Species'
    message(3) = '   specie <params>'
    message(4) = '%'
    call write_fatal(4)    
  end if
  allocate(spec(nspecies))

  write(*,*) 'step5'
  do i = 1, nspecies
    call oct_parse_block_str(str, i-1, 0, spec(i)%label)
    call oct_parse_block_double(str, i-1, 1, spec(i)%weight)
  end do

  write(*,*) 'step6'
  ! Initializes the atom
  call atom_init(natoms, atm, ncatoms, catm, nspecies, spec)

  ! Opens the nbo file
  call io_assign(nbo_unit)
  open(unit = nbo_unit, file = nbofile, action = 'read')

  ! Opens the xyz file
  call io_assign(xyz_unit)
  open(unit = xyz_unit, file = xyzfile, action = 'write')

  ! Skips the header
  rewind(unit = nbo_unit)
  read(unit = nbo_unit, fmt = *); read(unit = nbo_unit, fmt = *)

  ierr = 0
  do while(ierr == 0)
     read(unit = nbo_unit, iostat = ierr, fmt = *) iter, dump, dump, dump, dump, &
         ((atm(i)%x(j), j = 1, 3), i = 1, natoms)
     if(mod(iter, sampling) == 0) then
       call write_xyz
     endif
  enddo

  call io_close(nbo_unit); call io_close(xyz_unit)
  stop

contains


subroutine write_xyz

  integer :: i
  
  ! xyz format, for easy plot in rasmol
    write(xyz_unit, '(i4)') natoms
    write(xyz_unit, '(1x)')
    do i = 1, natoms
      write(xyz_unit, '(6x,a,2x,3f12.6)') atm(i)%spec%label, atm(i)%x(:)/units_out%length%factor
    end do

  return
end subroutine write_xyz

end subroutine nbo2xyz
