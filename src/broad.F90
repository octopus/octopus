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

program broad
  use io
  use units

  implicit none

  type broad_type
    FLOAT :: b, energy_step, min_energy, max_energy
  end type broad_type

  integer :: ierr
  logical :: l
  type(broad_type) :: b

  ! Initialize stuff
  call global_init()
  call units_init()

  ! broadening to use
  call oct_parse_float("LinBroadnening", CNST(0.02)/units_inp%energy%factor, b%b)
  call oct_parse_float("LinEnergyStep",  CNST(0.01)/units_inp%energy%factor, b%energy_step)
  call oct_parse_float("LinMinEnergy",   M_ZERO, b%min_energy)
  call oct_parse_float("LinMaxEnergy",   M_ONE/units_inp%energy%factor, b%max_energy)

  b%b = b%b * units_inp%energy%factor
  b%energy_step = b%energy_step * units_inp%energy%factor
  b%min_energy  = b%min_energy  * units_inp%energy%factor
  b%max_energy  = b%max_energy  * units_inp%energy%factor

  ! calculate resonances
  message(1) = "Info: Broadening spectra"
  call write_info(1)

  call oct_parse_logical("LinEigenvalues", .true., l)
  if(l) then
    message(1) = "      Eigenvalues"
    call write_info(1)
    call calc_broad(b, 'linear', 'eps-diff', .true.)
  end if

  call oct_parse_logical("LinPetersilka", .true., l)
  if(l) then
    message(1) = "      Petersilka"
    call write_info(1)
    call calc_broad(b, 'linear', 'petersilka', .true.)
  end if

  call oct_parse_logical("LinCasida", .true., l)
  if(l) then
    message(1) = "      Casida"
    call write_info(1)
    call calc_broad(b, 'linear', 'casida', .false.)
  end if

contains
  subroutine calc_broad(b, dir, fname, extracols)
    type(broad_type), intent(in) :: b
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: fname
    logical, intent(in) :: extracols
    
    FLOAT, allocatable :: s(:,:)
    FLOAT :: w, e, f(4)
    integer :: istat, n, iunit, j1, j2

    n = (b%max_energy - b%min_energy) / b%energy_step
    allocate(s(4,n))
    s = M_ZERO

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+fname, status='old', iostat=istat)
    if(istat.ne.0) then
      write(message(1),'(3a)') "Could not open multipoles file!"
      call write_fatal(1)
    end if
    read(iunit, *) ! skip header

    do
      if(extracols) then
        read(iunit, *, end=100) j1, j2, e, f
      else
        read(iunit, *, end=100) e, f
      end if
      e = e * units_inp%energy%factor

      do j1 = 1, n
        w = b%min_energy + real(j1-1, PRECISION)*b%energy_step
        s(:, j1) = s(:, j1) + f(:)*b%b/((w-e)**2 + b%b**2)/M_PI ! lorentzian
      end do
    end do
100   continue
    call io_close(iunit)
      
    ! print spectra
    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/spectrum."+fname, status='unknown')
    do j1 = 1, n
      w = b%min_energy + real(j1-1, PRECISION)*b%energy_step
      write(iunit, '(5es14.6)') w/units_inp%energy%factor, s(:, j1)*units_inp%energy%factor
    end do
    call io_close(iunit)

    deallocate(s)
  end subroutine calc_broad

end program broad
