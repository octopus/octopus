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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

  program vibrational
    use batch_oct_m
    use command_line_oct_m
    use global_oct_m
    use io_oct_m
    use ions_oct_m
    use messages_oct_m
    use namespace_oct_m
    use parser_oct_m
    use profiling_oct_m
    use spectrum_oct_m
    use unit_oct_m
    use unit_system_oct_m

    implicit none

    integer :: iunit, ierr, ii, jj, iter, read_iter, ntime, nvaf, nvel, ivel
    FLOAT, allocatable :: vaf(:), time(:), velocities(:, :), ftvaf(:)
    type(ions_t),     pointer :: ions 
    type(spectrum_t) :: spectrum
    type(batch_t) :: vafb, ftvafb
    FLOAT :: ww, curtime, vaftime, deltat
    integer :: ifreq, max_freq
    integer :: skip
    
    ! Initialize stuff
    call global_init(is_serial = .true.)		 

    call getopt_init(ierr)
    call getopt_end()

    call parser_init()
    
    call messages_init()

    call io_init()

    call unit_system_init(global_namespace)

    call spectrum_init(spectrum, global_namespace, &
      default_energy_step = units_to_atomic(unit_invcm, CNST(0.2)), &
      default_max_energy  = units_to_atomic(unit_invcm, CNST(5000.0)))
 
    !%Variable VibrationalSpectrumTimeStepFactor
    !%Type integer
    !%Default 10
    !%Section Utilities::oct-vibrational_spectrum
    !%Description
    !% In the calculation of the vibrational spectrum, it is not necessary
    !% to read the velocity at every time step. This variable controls
    !% the integer factor between the simulation time step and the
    !% time step used to calculate the vibrational spectrum.
    !%End

    call messages_obsolete_variable(global_namespace, 'PropagationSpectrumTimeStepFactor', 'VibrationalSpectrumTimeStepFactor')
    call parse_variable(global_namespace, 'VibrationalSpectrumTimeStepFactor', 10, skip)
    if(skip <= 0) call messages_input_error(global_namespace, 'VibrationalSpectrumTimeStepFactor')

    max_freq = spectrum_nenergy_steps(spectrum)

    if (spectrum%end_time < M_ZERO) spectrum%end_time = huge(spectrum%end_time)

    ions => ions_t(global_namespace)

    ! Opens the coordinates files.
    iunit = io_open('td.general/coordinates', global_namespace, action='read')

    call io_skip_header(iunit)

    ntime = 1
    iter = 1

    ! check the number of time steps we will read
    do
      read(unit = iunit, iostat = ierr, fmt = *) read_iter, curtime, &
        ((ions%pos(jj, ii), jj = 1, 3), ii = 1, ions%natoms), &
        ((ions%vel(jj, ii), jj = 1, 3), ii = 1, ions%natoms)

      curtime = units_to_atomic(units_out%time, curtime)

      if(ierr /= 0 .or. curtime >= spectrum%end_time) then
        iter = iter - 1	! last iteration is not valid
        ntime = ntime - 1
        exit
      end if

      if(iter /= read_iter + 1) then
        call messages_write("Error while reading file 'td.general/coordinates',", new_line = .true.)
        call messages_write('expected iteration ')
        call messages_write(iter - 1)
        call messages_write(', got iteration ')
        call messages_write(read_iter)
        call messages_write('.')
        call messages_fatal()
      end if
      
      ! ntime counts how many steps are gonna be used
      if (curtime >= spectrum%start_time .and. mod(iter, skip) == 0) ntime = ntime + 1 

      iter = iter + 1 !counts number of timesteps (with time larger than zero up to SpecEndTime)
    end do

    call io_close(iunit)

    nvel = ions%natoms*ions%space%dim

    SAFE_ALLOCATE(time(1:ntime))
    SAFE_ALLOCATE(velocities(1:nvel, 1:ntime))

    ! Opens the coordinates files.
    iunit = io_open('td.general/coordinates', global_namespace, action='read')

    call io_skip_header(iunit)

    ntime = 1
    iter = 1

    do
      read(unit = iunit, iostat = ierr, fmt = *) read_iter, curtime, &
        ((ions%pos(jj, ii), jj = 1, 3), ii = 1, ions%natoms), &
        ((ions%vel(jj, ii), jj = 1, 3), ii = 1, ions%natoms)

      curtime = units_to_atomic(units_out%time, curtime)

      if(ierr /= 0 .or. curtime >= spectrum%end_time) then
        iter = iter - 1	! last iteration is not valid
        ntime = ntime - 1
        exit
      end if

      ASSERT(iter == read_iter + 1)
      
      if (curtime >= spectrum%start_time .and. mod(iter, skip) == 0) then

        time(ntime) = curtime
        ivel = 1
        do ii = 1, ions%natoms
          do jj = 1, ions%space%dim
            velocities(ivel, ntime) = units_to_atomic(units_out%velocity, ions%vel(jj, ii))
            ivel = ivel + 1
          end do
        end do

        ntime = ntime + 1 
      end if

      iter = iter + 1 !counts number of timesteps (with time larger than zero up to SpecEndTime)
    end do

    call io_close(iunit)

    deltat = time(2) - time(1)

    !%Variable VibrationalSpectrumTime
    !%Type integer
    !%Section Utilities::oct-vibrational_spectrum
    !%Description
    !% This variable controls the maximum time for the calculation of
    !% the velocity autocorrelation function. The default is the total
    !% propagation time.
    !%End
    call parse_variable(global_namespace, 'VibrationalSpectrumTime', ntime*deltat, vaftime)

    nvaf = int(vaftime/deltat)

    SAFE_ALLOCATE(vaf(1:ntime))

    call messages_write('Time step = ')
    call messages_write(deltat, units = units_out%time)
    call messages_info()
    
    call messages_new_line()
    call messages_write('Calculating the velocity autocorrelation function')
    call messages_info()

    call calculate_vaf(vaf)

   !print the vaf
    iunit = io_open('td.general/velocity_autocorrelation', global_namespace, action='write')

800 FORMAT(80('#'))      
    write(unit = iunit, iostat = ierr, fmt = 800) 
    write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
    write(unit = iunit, iostat = ierr, fmt = '(a,4x,a6,a7,a1,10x,a10)') &
      '#       Iter', 'time [',units_out%time%abbrev,']', 'VAF [a.u.]'
    write(unit = iunit, iostat = ierr, fmt = 800) 

    do jj = 1, nvaf
      write(unit = iunit, iostat = ierr, fmt = *) jj, units_from_atomic(units_out%time, (jj - 1)*deltat), vaf(jj)
    end do

    call io_close(iunit)


    SAFE_ALLOCATE(ftvaf(1:max_freq))

    ftvaf = M_ONE

    call batch_init(vafb, vaf)

    call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, 1, nvaf, M_ZERO, deltat, vafb)

    call batch_init(ftvafb, ftvaf)

    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
      1, nvaf, M_ZERO, deltat, vafb, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, ftvafb)

    call vafb%end()
    call ftvafb%end()


    !and print the spectrum
    iunit = io_open('td.general/vibrational_spectrum', global_namespace, action='write')

    write(unit = iunit, iostat = ierr, fmt = 800) 
    write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
    write(unit = iunit, iostat = ierr, fmt = '(a17,6x,a15)') '#   Energy [1/cm]', 'Spectrum [a.u.]'
    write(unit = iunit, iostat = ierr, fmt = 800 ) 

    do ifreq = 1, max_freq
      ww = spectrum%energy_step*(ifreq - 1) + spectrum%min_energy
      write(unit = iunit, iostat = ierr, fmt = '(2e20.10)') units_from_atomic(unit_invcm, ww), ftvaf(ifreq)
    end do

    call io_close(iunit)

    SAFE_DEALLOCATE_A(vaf)

    SAFE_DEALLOCATE_A(ftvaf)

    SAFE_DEALLOCATE_P(ions)

    SAFE_DEALLOCATE_A(time)

    call io_end()
    call messages_end()

    call parser_end()
    call global_end()

  contains

    subroutine calculate_vaf(vaf)
      FLOAT, intent(out) :: vaf(:)
      integer :: itm, itn

      PUSH_SUB(calculate_vaf)

      write (message(1), '(a)') "Read velocities from '"// &
        trim(io_workpath('td.general/coordinates', global_namespace))//"'"
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

        vaf(itm) = vaf(itm)/TOFLOAT(ntime - itm + 1)

      end do

      do itm = 2, ntime
        vaf(itm) = vaf(itm)/vaf(1)
      end do
      vaf(1) = M_ONE

      POP_SUB(calculate_vaf)
    end subroutine calculate_vaf

  end program vibrational

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
