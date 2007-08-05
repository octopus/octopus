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

  integer :: iunit, ierr, ii, jj, iter, read_iter, max_iter, ini_iter, end_iter
  FLOAT :: start_time, end_time
  FLOAT, allocatable :: vini(:,:), vaf(:), time(:)
  CMPLX, allocatable :: vafft(:)
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

  !These variables are documented somewhere else
  call loct_parse_int(check_inp('TDMaximumIter'), 1500, max_iter)
  call loct_parse_float(check_inp('SpecStartTime'),  M_ZERO, start_time)
  call loct_parse_float(check_inp('SpecEndTime'),  -M_ONE, end_time)

  if (end_time < M_ZERO) end_time = huge(REAL_PRECISION)

  ALLOCATE(vini(1:3, geo%natoms), 3*geo%natoms)
  ALLOCATE(vaf(0:max_iter),  max_iter-ini_iter)
  ALLOCATE(time(0:max_iter), max_iter-ini_iter + 2)
  
  ! Opens the coordinates files.
  iunit = io_open('td.general/coordinates', action='read')

  call io_skip_header(iunit)

  ini_iter = -1
  iter = 0
  do while(.true.)

    read(unit = iunit, iostat = ierr, fmt = *) read_iter, time(iter), &
         ((geo%atom(ii)%x(jj), jj = 1, 3), ii = 1, geo%natoms),&
         ((geo%atom(ii)%v(jj), jj = 1, 3), ii = 1, geo%natoms)

    if (ierr /= 0) then 
      iter = iter - 1 !last iteration is not valid
      exit
    end if

    ASSERT(iter == read_iter)

    if (time(iter) >= end_time) exit

    if (time(iter) >= start_time) then

      if(ini_iter == -1) then 
        ini_iter = iter
        do ii = 1, geo%natoms
          vini(1:3, ii) = geo%atom(ii)%v(1:3)
        end do
      end if

      !calculate the vaf
      vaf(iter) = M_ZERO
      do ii = 1, geo%natoms
        vaf(iter) = vaf(iter) + sum(geo%atom(ii)%v(1:3) * vini(1:3, ii))
      end do
      
    end if

    iter = iter + 1
  end do

  call io_close(iunit)

  if(ini_iter == 0 ) ini_iter = 1
  end_iter = iter - 1

  write (message(1), '(a)') "Read velocities from 'td.general/coordinates'"
  call write_info(1)

  call fourier

  call geometry_end(geo)
  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()
  
contains
  
  subroutine fourier
    FLOAT :: ww, av
    FLOAT, parameter :: dw = 1.0/hartree_to_cm_inv
    integer, parameter :: max_freq = 10000
    integer :: ifreq, count

    ALLOCATE(vafft(1:max_freq), max_freq)

    av = maxval(abs(vaf(ini_iter:end_iter)))

    if( av < CNST(1e-12)) then 
      write (message(1), '(a)') "Error: Velocity autocorrelation function is zero"
      call write_fatal(1)
    end if

    vaf(ini_iter:end_iter) = vaf(ini_iter:end_iter)/av

    !apply an envelope
    do jj = ini_iter, end_iter
      vaf(jj) = vaf(jj) * sin((time(jj)-time(ini_iter+1))*M_PI/(time(end_iter)-time(ini_iter)))
    end do

    !remove the dc component
    av = M_ZERO
    count = 0
    do jj = ini_iter, end_iter
      av = av + vaf(jj) * M_HALF*(time(jj+1)-time(jj-1))
      count = count + 1
    end do
    
    do jj = ini_iter, end_iter
      vaf(jj) = vaf(jj) - av/(M_HALF*(time(jj+1)-time(jj-1))*count)
    end do

    iunit = io_open('td.general/velocity_autocorrelation', action='write')

    do jj = ini_iter, end_iter
      write(unit = iunit, iostat = ierr, fmt = *) time(jj), vaf(jj)
    end do

    !print again to see the matching
    do jj = ini_iter, end_iter
      write(unit = iunit, iostat = ierr, fmt = *) (time(end_iter)-time(ini_iter)) + time(jj), vaf(jj)
    end do

    call io_close(iunit)

    write (message(1), '(a)') "Taking the fourier transform of the velocity autocorrelation function"
    call write_info(1)

    !now calculate the FT
    !$omp parallel do private(ww, jj)
    do ifreq = 1, max_freq
      ww = dw * ifreq
      vafft(ifreq) = M_ZERO
      do jj = ini_iter, end_iter
        vafft(ifreq) = vafft(ifreq) + &
             exp(M_zI * ww * time(jj) * units_out%time%factor) * vaf(jj)*M_HALF*(time(jj+1)-time(jj-1))
      end do
    end do

    write (message(1), '(a)') "Done."
    call write_info(1)

    iunit = io_open('td.general/vibrational_spectrum', action='write')
    
    do ifreq = 1, max_freq
      ww = dw * ifreq
      write(unit = iunit, iostat = ierr, fmt = '(4e20.10)') &
           ww*hartree_to_cm_inv, abs(vafft(ifreq)), real(vafft(ifreq)), aimag(vafft(ifreq))
    end do

    call io_close(iunit)

  end subroutine fourier

end program phonon_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
