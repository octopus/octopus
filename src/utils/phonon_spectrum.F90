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
  use math_m
  use messages_m
  use datasets_m
  use io_m
  use loct_parser_m
  use profiling_m
  use units_m
  use varinfo_m

  implicit none

  integer, parameter :: SPEC_VIBRATIONAL = 1, SPEC_INFRARED = 2
  
  integer :: mode

  integer :: iunit, ierr, ii, jj, iter, read_iter, max_iter, ini_iter, end_iter
  FLOAT :: start_time, end_time
  FLOAT, allocatable :: vaf(:), time(:), dipole(:,:)
  CMPLX, allocatable :: ftvaf(:), ftdipole(:,:)
  type(geometry_t) :: geo 
  
  FLOAT :: ww, av, irtotal
  FLOAT, parameter :: dw = M_ONE/hartree_to_cm_inv
  integer :: ifreq, idir
  integer, parameter :: max_freq = 10000
  
  ! Initialize stuff
  call global_init()
  call parser_init()
  call datasets_init(1)
  call io_init()
  if(in_debug_mode) then
    call io_mkdir('debug')
  end if
  call units_init()

  !These variables are documented in src/td/spectrum.F90
  call loct_parse_int(datasets_check('TDMaximumIter'), 1500, max_iter)
  call loct_parse_float(datasets_check('SpecStartTime'),  M_ZERO, start_time)
  call loct_parse_float(datasets_check('SpecEndTime'),  -M_ONE, end_time)


  !%Variable SpecVibrational
  !%Type integer
  !%Default vibrational
  !%Section Utilities::oct-vibrational
  !%Description
  !% This variable select the kind of spectrum that will be
  !% calculated from a molecular dynamics run.
  !%Option vibrational 1
  !% Vibrational spectrum from the velocity autocorrelation function.
  !%Option infrared    2
  !% Infrared spectrum obtained from the dipole moment.
  !%End
  call loct_parse_int  (datasets_check('SpecVibrational'), SPEC_VIBRATIONAL, mode)
  if(.not.varinfo_valid_option('SpecVibrational', mode)) call input_error('SpecVibrational')

  if (end_time < M_ZERO) end_time = huge(end_time)

  SAFE_ALLOCATE(time(0:max_iter+1))
  
  select case(mode)
    case(SPEC_VIBRATIONAL) 

      call geometry_init(geo)

      SAFE_ALLOCATE(vaf(0:max_iter))

      call read_vaf(vaf)

      av = maxval(abs(vaf(ini_iter:end_iter)))
      
      if( av < CNST(1e-12)) then 
        write (message(1), '(a)') "Error: Velocity autocorrelaion function is zero"
        call write_fatal(1)
      end if
      
      vaf(ini_iter:end_iter) = vaf(ini_iter:end_iter)/av
      
      SAFE_ALLOCATE(ftvaf(1:max_freq))
    
      call fourier(vaf, ftvaf)

      !print the vaf
      iunit = io_open('td.general/velocity_autocorrelation', action='write')
      
      do jj = ini_iter, end_iter
        write(unit = iunit, iostat = ierr, fmt = *) time(jj), vaf(jj)
      end do
      
      ! print again to see the matching
      do jj = ini_iter, end_iter
        write(unit = iunit, iostat = ierr, fmt = *) (time(end_iter)-time(ini_iter)) + time(jj), vaf(jj)
      end do
      
      call io_close(iunit)
      
      !and print the spectrum
      iunit = io_open('td.general/vibrational_spectrum', action='write')
      
      do ifreq = 1, max_freq
        ww = dw * ifreq
        write(unit = iunit, iostat = ierr, fmt = '(4e20.10)') &
             ww*hartree_to_cm_inv, abs(ftvaf(ifreq)), real(ftvaf(ifreq)), aimag(ftvaf(ifreq))
      end do
      
      call io_close(iunit)
      
      SAFE_DEALLOCATE_A(vaf)
      SAFE_DEALLOCATE_A(ftvaf)

      call geometry_end(geo)

    case(SPEC_INFRARED)

      SAFE_ALLOCATE(dipole(0:max_iter+1, 1:3))

      call read_dipole(dipole)

      SAFE_ALLOCATE(ftdipole(1:max_freq, 1:3))

      do idir = 1, 3
        call fourier(dipole(:, idir), ftdipole(:, idir))
      end do

      !and print the spectrum
      iunit = io_open('td.general/infrared', action='write')
      
      do ifreq = 1, max_freq
        ww = dw * ifreq
        irtotal = sqrt(sum( abs(ftdipole(ifreq, 1:3))**2 ))
        write(unit = iunit, iostat = ierr, fmt = '(5e20.10)') &
             ww*hartree_to_cm_inv, ww*irtotal, ww*abs(ftdipole(ifreq, 1:3))
      end do
      call io_close(iunit)

  end select
  
  SAFE_DEALLOCATE_A(time)

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()

contains

  subroutine read_vaf(vaf)
    FLOAT, intent(out) :: vaf(0:)

    FLOAT, allocatable :: vini(:,:)

    SAFE_ALLOCATE(vini(1:3, 1:geo%natoms))

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

    write (message(1), '(a)') "Read velocities from '"// &
      trim(io_workpath('td.general/coordinates'))//"'"
    call write_info(1)
    
    SAFE_DEALLOCATE_A(vini)

  end subroutine read_vaf

  subroutine read_dipole(dipole)
    FLOAT,   intent(out)   :: dipole(0:, :)

    FLOAT :: charge

    ! Opens the coordinates files.
    iunit = io_open('td.general/multipoles', action='read')

    call io_skip_header(iunit)

    ini_iter = -1
    iter = 0
    do while(.true.)

      read(unit = iunit, iostat = ierr, fmt = *) read_iter, time(iter), &
           charge, dipole(iter, 1), dipole(iter, 2), dipole(iter, 3)

      if (ierr /= 0) then 
        iter = iter - 1 !last iteration is not valid
        exit
      end if

      ASSERT(iter == read_iter)

      if (time(iter) >= end_time) exit

      if (time(iter) >= start_time .and. ini_iter == -1) then 
          ini_iter = iter
      end if

      iter = iter + 1
    end do

    call io_close(iunit)

    if(ini_iter == 0 ) ini_iter = 1
    end_iter = iter - 1

    write (message(1), '(a)') "Read dipole moment from '"// &
      trim(io_workpath('td.general/multipoles'))//"'."
    call write_info(1)
  end subroutine read_dipole
  
  subroutine fourier(fi, ftfi)
    FLOAT, intent(inout)  :: fi(:)
    CMPLX, intent(out)    :: ftfi(:)

    FLOAT :: ww, av
    integer :: ifreq, count

    !apply an envelope
    do jj = ini_iter, end_iter
      fi(jj) = fi(jj) * sin((time(jj)-time(ini_iter+1))*M_PI/(time(end_iter)-time(ini_iter)))
    end do

    !remove the dc component
    av = M_ZERO
    count = 0
    do jj = ini_iter, end_iter
      av = av + fi(jj) * M_HALF*(time(jj+1)-time(jj-1))
      count = count + 1
    end do
    
    do jj = ini_iter, end_iter
      fi(jj) = fi(jj) - av/(M_HALF*(time(jj+1)-time(jj-1))*count)
    end do

    write (message(1), '(a)') "Taking the fourier transform."
    call write_info(1)

    !now calculate the FT
    !$omp parallel do private(ww, jj)
    do ifreq = 1, max_freq
      ww = dw * ifreq
      ftfi(ifreq) = M_ZERO
      do jj = ini_iter, end_iter
        ftfi(ifreq) = ftfi(ifreq) + &
             exp(M_zI * ww * time(jj) * units_out%time%factor) * fi(jj)*M_HALF*(time(jj+1)-time(jj-1))
      end do
    end do
    !$omp end parallel do

    write (message(1), '(a)') "Done."
    call write_info(1)

  end subroutine fourier

end program phonon_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
