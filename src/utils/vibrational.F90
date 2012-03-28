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

  program vibrational
    use batch_m
    use command_line_m
    use datasets_m
    use geometry_m
    use global_m
    use io_m
    use math_m
    use messages_m
    use parser_m
    use profiling_m
    use simul_box_m
    use space_m
    use spectrum_m
    use unit_m
    use unit_system_m
    use varinfo_m

    implicit none

    integer, parameter :: SPEC_VIBRATIONAL = 1, SPEC_INFRARED = 2

    integer :: mode

    integer :: iunit, ierr, ii, jj, kk, iter, read_iter, max_iter, ini_iter, end_iter, ntime, nvaf, nvel, ivel
    FLOAT :: start_time, end_time
    FLOAT, allocatable :: vaf(:), time(:), dipole(:,:), velocities(:, :), ftvaf(:)
    CMPLX, allocatable :: ftdipole(:,:)
    type(geometry_t)  :: geo 
    type(space_t)     :: space
    type(simul_box_t) :: sb
    type(spec_t) :: spectrum
    type(batch_t) :: vafb, ftvafb
    FLOAT :: ww, av, irtotal
    FLOAT :: dw, max_energy, curtime
    integer :: ifreq, idir
    integer, parameter :: max_freq = 10000

    ! Initialize stuff
    call global_init()		 

    call getopt_init(ierr)
    mode = SPEC_VIBRATIONAL
    if(ierr.eq.0) call getopt_vibrational(mode)
    call getopt_end()

    call parser_init()
    call messages_init()

    call datasets_init(1)
    call io_init()

    call unit_system_init()

    !These variables are documented in src/td/spectrum.F90
    call parse_integer(datasets_check('TDMaximumIter'), 1500, max_iter)
    call parse_float(datasets_check('PropagationSpectrumStartTime'),  M_ZERO, start_time, units_inp%time)
    call parse_float(datasets_check('PropagationSpectrumEndTime'),  -M_ONE, end_time, units_inp%time)
    call parse_float(datasets_check('PropagationSpectrumMaxEnergy'), &
      units_from_atomic(units_inp%energy, units_to_atomic(unit_invcm, CNST(10000.0))), max_energy, units_inp%energy)

    dw = max_energy/(max_freq-M_ONE) !Initializes the wavevector step dw

    if (end_time < M_ZERO) end_time = huge(end_time)

    call space_init(space)
    call geometry_init(geo, space)
    call simul_box_init(sb, geo, space)

    select case(mode)
    case(SPEC_VIBRATIONAL) 

      ! Opens the coordinates files.
      iunit = io_open('td.general/coordinates', action='read')

      call io_skip_header(iunit)

      ntime = 1
      iter = 1

      ! check the number of time steps we will read
      do
        read(unit = iunit, iostat = ierr, fmt = *) read_iter, curtime, &
          & ((geo%atom(ii)%x(jj), jj = 1, 3), ii = 1, geo%natoms), &
          & ((geo%atom(ii)%v(jj), jj = 1, 3), ii = 1, geo%natoms)

        curtime = units_to_atomic(units_out%time, curtime)

        if(ierr.ne.0) then
          iter = iter - 1	! last iteration is not valid
          ntime = ntime - 1
          exit
        end if

        ASSERT(iter == read_iter + 1)

        if (curtime >= end_time) exit

        if (curtime >= start_time) ntime = ntime + 1 ! ntime counts how many steps are gonna be used

        iter = iter + 1 !counts number of timesteps (with time larger than zero up to SpecEndTime)
      end do

      call io_close(iunit)

      nvel = geo%natoms*space%dim

      SAFE_ALLOCATE(time(1:ntime))
      SAFE_ALLOCATE(velocities(1:nvel, 1:ntime))

      ! Opens the coordinates files.
      iunit = io_open('td.general/coordinates', action='read')

      call io_skip_header(iunit)

      ntime = 1
      iter = 1

      do
        read(unit = iunit, iostat = ierr, fmt = *) read_iter, curtime, &
          & ((geo%atom(ii)%x(jj), jj = 1, 3), ii = 1, geo%natoms), &
          & ((geo%atom(ii)%v(jj), jj = 1, 3), ii = 1, geo%natoms)

        curtime = units_to_atomic(units_out%time, curtime)

        if(ierr.ne.0) then
          iter = iter - 1	! last iteration is not valid
          ntime = ntime - 1
          exit
        end if

        ASSERT(iter == read_iter + 1)

        if (curtime >= end_time) exit

        if (curtime >= start_time) then

          time(ntime) = curtime
          ivel = 1
          do ii = 1, geo%natoms
            do jj = 1, space%dim
              velocities(ivel, ntime) = units_to_atomic(units_out%velocity, geo%atom(ii)%v(jj))
              ivel = ivel + 1
            end do
          end do

          ntime = ntime + 1 ! ntime counts how many steps are gonna be used
        end if

        iter = iter + 1 !counts number of timesteps (with time larger than zero up to SpecEndTime)
      end do

      call io_close(iunit)

      nvaf = ntime

      SAFE_ALLOCATE(vaf(1:ntime))

      call read_vaf(vaf)

      ini_iter=1
      end_iter=nvaf

      !print the vaf
      iunit = io_open('td.general/velocity_autocorrelation', action='write')

800   FORMAT(80('#'))      
      write(unit = iunit, iostat = ierr, fmt = 800) 
      write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
      write(unit = iunit, iostat = ierr, fmt = '(a1,4x,a6,a7,a1,10x,a10)') '#', &
        &  'time [',units_out%time%abbrev,']', 'VAF [a.u.]'
      write(unit = iunit, iostat = ierr, fmt = 800) 

      do jj = ini_iter, end_iter
        write(unit = iunit, iostat = ierr, fmt = *) units_from_atomic(units_out%time, time(jj)), vaf(jj)
      end do

      !print again to see the matching
      do jj = ini_iter, end_iter
        write(unit = iunit, iostat = ierr, fmt = *) &
          units_from_atomic(units_out%time, (time(end_iter)-time(ini_iter)) + time(jj)), vaf(jj)
      end do

      call io_close(iunit)


      SAFE_ALLOCATE(ftvaf(1:max_freq))
      
      call spectrum_init(spectrum)

      call batch_init(vafb, 1)
      call batch_add_state(vafb, vaf)

      call signal_damp(spectrum%damp, spectrum%damp_factor, 1, ntime, time(2) - time(1), vafb)

      call batch_init(ftvafb, 1)
      call batch_add_state(ftvafb, ftvaf)

      call fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
        1, ntime, M_ZERO, time(2) - time(1), vafb, 1, max_freq, dw, ftvafb)

      call batch_end(vafb)
      call batch_end(ftvafb)


      !print the vaf
      iunit = io_open('td.general/velocity_autocorrelation', action='write')

      write(unit = iunit, iostat = ierr, fmt = "(80('#'))") 
      write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
      write(unit = iunit, iostat = ierr, fmt = '(a1,4x,a6,a7,a1,10x,a10)') '#', &
        &  'time [',units_out%time%abbrev,']', 'VAF [a.u.]'
      write(unit = iunit, iostat = ierr, fmt = 800) 

      do jj = ini_iter, end_iter
        write(unit = iunit, iostat = ierr, fmt = *) units_from_atomic(units_out%time, time(jj)), vaf(jj)
      end do

      !print again to see the matching
      do jj = ini_iter, end_iter
        write(unit = iunit, iostat = ierr, fmt = *) &
          units_from_atomic(units_out%time, (time(end_iter)-time(ini_iter)) + time(jj)), vaf(jj)
      end do

      call io_close(iunit)

      !print the vaf
      iunit = io_open('td.general/velocity_autocorrelation', action='write')
      !and print the spectrum
      iunit = io_open('td.general/vibrational_spectrum', action='write')

      write(unit = iunit, iostat = ierr, fmt = 800) 
      write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
      write(unit = iunit, iostat = ierr, fmt = '(a17,8x,a15,5x,a13,5x,a13)') &
        '#   Energy [1/cm]', 'Spectrum [a.u.]', 'Re(FT of vaf)', 'Im(FT of vaf)'      
      write(unit = iunit, iostat = ierr, fmt = 800 ) 

      do ifreq = 1, max_freq
        ww = dw * ifreq
        write(unit = iunit, iostat = ierr, fmt = '(2e20.10)') units_from_atomic(unit_invcm, ww), ftvaf(ifreq)
      end do

      call io_close(iunit)

      SAFE_DEALLOCATE_A(vaf)

      SAFE_DEALLOCATE_A(ftvaf)

    case(SPEC_INFRARED)

      SAFE_ALLOCATE(dipole(0:max_iter+1, 1:3))

      call read_dipole(dipole)

      SAFE_ALLOCATE(ftdipole(1:max_freq, 1:3))

      do idir = 1, 3
        call fourier(dipole(:, idir), ftdipole(:, idir))
      end do

      !and print the spectrum
      iunit = io_open('td.general/infrared', action='write')

100   FORMAT(100('#'))

      write(unit = iunit, iostat = ierr, fmt = 100)
      write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
      write(unit = iunit, iostat = ierr, fmt = '(a25,a1,a1,a10)') &
        & '# all absorptions are in ',units_out%length%abbrev,'*',units_out%time%abbrev
      write(unit = iunit, iostat = ierr, fmt = '(a1)') '#'
      write(unit = iunit, iostat = ierr, fmt = '(a19,41x,a17)') '#        Energy    ', 'absorption'
      write(unit = iunit, iostat = ierr, fmt = '(a15,13x,a5,15x,a7,13x,a7,13x,a7)') &
        & '#        [1/cm]','total','FT(<x>)','FT(<y>)','FT(<z>)'
      write(unit = iunit, iostat = ierr, fmt = 0100) 

      do ifreq = 1, max_freq
        ww = dw * ifreq
        irtotal = sqrt(sum( abs(ftdipole(ifreq, 1:3))**2 ))
        write(unit = iunit, iostat = ierr, fmt = '(5e20.10)') &
          units_from_atomic(unit_invcm, ww), units_from_atomic(units_out%length*units_out%time, ww*irtotal), &
          (units_from_atomic(units_out%length*units_out%time,  ww*abs(ftdipole(ifreq, idir))), idir = 1, 3)
      end do
      call io_close(iunit)

    end select

    call simul_box_end(sb)
    call geometry_end(geo)
    call space_end(space)

    SAFE_DEALLOCATE_A(time)

    call io_end()
    call datasets_end()
    call messages_end()
    call parser_end()
    call global_end()

  contains

    subroutine read_vaf(vaf)
      FLOAT, intent(out) :: vaf(:)
      FLOAT, allocatable :: norm(:)
      FLOAT:: summand
      integer :: itm, itn, numiter

      write (message(1), '(a)') "Read velocities from '"// &
        trim(io_workpath('td.general/coordinates'))//"'"
      call messages_info(1)

      !calculating the vaf, formula from
      !
      ! http://www.timteatro.net/2010/09/29/velocity-autocorrelation-and-vibrational-spectrum-calculation/

      vaf = M_ZERO

      do itm = 1, ntime
        vaf(itm) = M_ZERO

        do itn = 1, ntime - itm + 1

          do ivel = 1, nvel
            vaf(itm) = vaf(itm) + velocities(ivel, itm + itn - 1)*velocities(ivel, itn)
          end do
        end do

        vaf(itm) = vaf(itm)/real(ntime - itm + 1, REAL_PRECISION)

      end do

      do itm = 2, ntime
        vaf(itm) = vaf(itm)/vaf(1)
      end do
      vaf(1) = M_ONE

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
          & charge, dipole(iter, 1), dipole(iter, 2), dipole(iter, 3)

        time(iter) =  units_to_atomic(units_out%time, time(iter))

        !dipole moment has unit of charge*length, charge has the same unit in both systems
        forall(ii= 1: 3) dipole(iter, ii) = units_to_atomic(units_out%length, dipole(iter, ii)) 


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
      call messages_info(1)
    end subroutine read_dipole

    subroutine fourier(fi, ftfi)
      implicit none
      FLOAT, intent(inout)  :: fi(:)
      CMPLX, intent(out)    :: ftfi(:)

      FLOAT :: ww, av, dt
      integer :: ifreq, count

      !apply an envelope
      do jj = ini_iter, end_iter
        fi(jj) = fi(jj) * sin((time(jj)-time(ini_iter+1))*M_PI/(time(end_iter)-time(ini_iter))) 
      end do

      dt = time(2) - time(1)

      !remove the DC component
      av = M_ZERO
      count = 0
      do jj = ini_iter, end_iter
        av = av + fi(jj)*dt 
        count = count + 1
      end do

      do jj = ini_iter, end_iter
        fi(jj) = fi(jj) - av/(dt*count) 
      end do

      write (message(1), '(a)') "Taking the fourier transform."
      call messages_info(1)

      !now calculate the FT
      !$omp parallel do private(ww, jj)
      do ifreq = 1, max_freq
        ww = dw * ifreq
        ftfi(ifreq) = M_ZERO
        do jj = ini_iter, end_iter
          ftfi(ifreq) = ftfi(ifreq) + &
            exp(M_zI * ww * time(jj))*fi(jj)*dt
        end do
      end do
      !$omp end parallel do

      write (message(1), '(a)') "Done."
      call messages_info(1)

    end subroutine fourier

  end program vibrational

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
