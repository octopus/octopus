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

program xyzanim
  use global
  use lib_oct
  use lib_oct_parser
  use io
  use units
  use atom
  use geometry

  implicit none

  character(len=80) :: str, nbofile, xyzfile
  character(len=5000) :: line
  integer :: ierr, sampling, natoms, ncatoms, nspecies, i, nbo_unit, xyz_unit, iter, j
  FLOAT :: dump

  type(geometry_type) :: geo

  ! Initialize stuff
  call global_init()
  call units_init()

  if(units_out%length%name .ne. 'Angstrom') then
    message(1) = 'Output in atomic units. Set UnitsOutput="eVA" if you want it in angstroms'
    call write_warning(1)
  endif
  if(conf%verbose<999) conf%verbose = -1

  ! Sets the filenames
  write(nbofile, '(a)') 'td.general/coordinates'
  write(xyzfile, '(a)') 'movie.xyz'

  ! how many do we have?
  str = "Species"
  geo%nspecies = loct_parse_block_n(str)
  if (nspecies < 1) then
    message(1) = "Input: Species block not specified"
    message(2) = '% Species'
    message(3) = '   specie <params>'
    message(4) = '%'
    call write_fatal(4)    
  end if
  allocate(geo%specie(geo%nspecies))

  ! how often do we sample?
  call loct_parse_int('AnimationSampling', 100, sampling)
  if(sampling < 1) then
    message(1) = 'Sampling rate (AnimationSampling) should be bigger than 0'
    call write_fatal(1)
  end if

  do i = 1, geo%nspecies
    call loct_parse_block_string(str, i-1, 0, geo%specie(i)%label)
    call loct_parse_block_float (str, i-1, 1, geo%specie(i)%weight)
  end do

  ! Initializes the atom system
  call geometry_init(geo, no_species_init=.true.)

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
         ((geo%atom(i)%x(j), j = 1, 3), i = 1, geo%natoms)
     if(mod(iter, sampling) == 0) then
       call write_xyz
     endif
  enddo

  call io_close(nbo_unit); call io_close(xyz_unit)

contains


subroutine write_xyz

  integer :: i
  
  ! xyz format, for easy plot in rasmol
    write(xyz_unit, '(i4)') geo%natoms
    write(xyz_unit, '(i10)') iter
    do i = 1, natoms
      write(xyz_unit, '(6x,a,2x,3f12.6)') geo%atom(i)%spec%label, geo%atom(i)%x(:)/units_out%length%factor
    end do

  return
end subroutine write_xyz

end program xyzanim
