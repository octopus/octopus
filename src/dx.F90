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
module dx

use global
use mesh
use liboct

implicit none

contains

subroutine create_dx_file(data_unit, datafile, dxheader_unit, m, f)
  integer, intent(in)          :: data_unit, dxheader_unit
  type(mesh_type), intent(in)  :: m
  real(r8), intent(in)         :: f(1:m%np)
  character(len=*), intent(in) :: datafile

  character(len=80) :: pwd, line
  real(r8), allocatable :: fcube(:, :, :)
  integer :: ix, iy, iz

  call clear_str(pwd)
  call oct_getcwd(pwd)

  allocate(fcube(m%fft_n(1), m%fft_n(2), m%fft_n(3)))
  fcube = 0.0_r8
  call R_FUNC(mesh_to_cube) (m, f, fcube)
    do ix = 1, m%fft_n(1)
       do iy = 1, m%fft_n(2)
          do iz = 1, m%fft_n(3)
             write(data_unit, *) fcube(ix, iy, iz)
          end do
       end do
    end do

  write(line,'(a)')              'file = '//trim(pwd(1:len_trim(pwd)-1))//'/'//trim(datafile)
     write(dxheader_unit,'(a)') trim(line)
  write(line,'(a,i4,a,i4,a,i4)') 'grid = ', m%fft_n(1), ' x ', m%fft_n(2), ' x ', m%fft_n(3);
     write(dxheader_unit,'(a)') trim(line)
  write(line,'(a)')              'format = ascii'
     write(dxheader_unit,'(a)') trim(line)
  write(line,'(a)')              'interleaving = record'
     write(dxheader_unit,'(a)') trim(line)
  write(line,'(a)')              'majority = row'
     write(dxheader_unit,'(a)') trim(line)
  write(line,'(a)')              'field = field0'
     write(dxheader_unit,'(a)') trim(line)
  write(line,'(a)')              'type = float'
     write(dxheader_unit,'(a)') trim(line)
  write(line,'(a)')              'dependency = positions'
     write(dxheader_unit,'(a)') trim(line)
  write(line,'(a)')              'positions = regular, regular, regular, 0, 1, 0, 1, 0, 1'
     write(dxheader_unit,'(a)') trim(line)

end subroutine create_dx_file

end module dx


