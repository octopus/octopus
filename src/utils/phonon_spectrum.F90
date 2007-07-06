!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

program phonon_spectrum
  use global_m
  use geometry_m
  use messages_m
  use datasets_m
  use io_m
  use lib_oct_parser_m
  use units_m

  implicit none

  integer :: iunit, ierr, ii, jj, iter, max_iter, ini_iter
  FLOAT :: tt
  FLOAT, allocatable :: vini(:,:), vaf(:), time(:), vafft(:)
  type(geometry_t) :: geo 

  ! Initialize stuff
  call global_init()
  call parser_init()
  call datasets_init(1)
  call io_init()
  if(in_debug_mode) then
    call io_mkdir('debug')
  end if
  call units_init()

  call geometry_init(geo)

  ini_iter = 0
  call loct_parse_int(check_inp('TDMaximumIter'), 1500, max_iter)

  ALLOCATE(vini(1:3, geo%natoms), 3*geo%natoms)
  ALLOCATE(vaf(ini_iter:max_iter),  max_iter-ini_iter)
  ALLOCATE(time(ini_iter:max_iter), max_iter-ini_iter + 2)
  
  ! Opens the coordinates files.
  iunit = io_open('td.general/coordinates', action='read')

  call io_skip_header(iunit)

  do while(.true.)

    read(unit = iunit, iostat = ierr, fmt = *) iter, tt, &
         ((geo%atom(ii)%x(jj), jj = 1, 3), ii = 1, geo%natoms),&
         ((geo%atom(ii)%v(jj), jj = 1, 3), ii = 1, geo%natoms)

    if (ierr /= 0) exit

    !get the initial velocities
    if (iter == ini_iter) then
      do ii = 1, geo%natoms
        vini(1:3, ii) = geo%atom(ii)%v(1:3)
      end do
    end if

    !calculate the vaf
    if (iter >= ini_iter) then
      time(iter) = tt
      vaf(iter) = M_ZERO
      do ii = 1, geo%natoms
        vaf(iter) = sum(geo%atom(ii)%v(1:3) * vini(1:3, ii))
      end do
    end if

  end do

  ! Close the file
  call io_close(iunit)

  !get the actual maximum time step
  max_iter = iter


  call fourier

  call geometry_end(geo)
  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()
  
contains
  
  subroutine fourier
    FLOAT :: ww, av
    FLOAT, parameter :: dw = 0.5/hartree_to_cm_inv
    integer, parameter :: max_freq = 20000
    integer :: ifreq, count

    ALLOCATE(vafft(1:max_freq), max_freq)

    !remove the dc component
    av = M_ZERO
    count = 0
    do jj = ini_iter+1, max_iter-1
      av = av + vaf(jj) * M_HALF*(time(jj+1)-time(jj-1))
      count = count + 1
    end do
    
    do jj = ini_iter+1, max_iter-1
      vaf(jj) = vaf(jj) - av/(M_HALF*(time(jj+1)-time(jj-1))*count)
    end do

    !now calculate the FT
    !$omp parallel do private(ww, jj)
    do ifreq = 1, max_freq
      ww = dw * ifreq
      vafft(ifreq) = M_ZERO
      do jj = ini_iter+1, max_iter-1
        vafft(ifreq) = vafft(ifreq) + cos(ww*time(jj))*vaf(jj)*M_HALF*(time(jj+1)-time(jj-1))
      end do
    end do

    iunit = io_open('td.general/vibrational_spectrum', action='write')

    do ifreq = 1, max_freq
      ww = dw * ifreq
      write(unit = iunit, iostat = ierr, fmt = *) ww*hartree_to_cm_inv, vafft(ifreq) / units_out%time%factor
    end do

    call io_close(iunit)

  end subroutine fourier

end program phonon_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
