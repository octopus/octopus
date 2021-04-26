!! Copyright (C) 2015 X. Andrade
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

  program conductivity
    use batch_oct_m
    use command_line_oct_m
    use geometry_oct_m
    use global_oct_m
    use grid_oct_m
    use io_oct_m
    use messages_oct_m
    use namespace_oct_m
    use parser_oct_m
    use profiling_oct_m
    use space_oct_m
    use spectrum_oct_m
    use species_oct_m
    use states_elec_oct_m
    use unit_oct_m
    use unit_system_oct_m

    implicit none

    integer :: iunit, ierr, ii, jj, iter, read_iter, ntime, nvel, ivel
    integer :: istart, iend, energy_steps, out_file
    FLOAT, allocatable :: time(:), velocities(:, :)
    FLOAT, allocatable :: total_current(:, :), ftcurr(:, :, :), curr(:, :)
    FLOAT, allocatable :: heat_current(:,:), ftheatcurr(:,:,:), heatcurr(:,:)
    CMPLX, allocatable :: invdielectric(:, :)
    type(geometry_t), pointer :: geo 
    type(space_t)     :: space
    type(spectrum_t) :: spectrum
    type(block_t)     :: blk
    type(batch_t) :: currb, ftcurrb, heatcurrb, ftheatcurrb
    FLOAT :: ww, curtime, deltat, velcm(1:MAX_DIM), vel0(1:MAX_DIM), current(1:MAX_DIM), integral(1:2), v0
    character(len=MAX_PATH_LEN) :: ref_filename
    integer :: ref_file, time_steps_ref, kk
    FLOAT, allocatable :: current_ref(:, :)
    FLOAT :: dt_ref, tt, start_time
    integer :: ifreq, max_freq
    integer :: skip
    logical :: from_forces
    character(len=120) :: header
    FLOAT :: excess_charge, qtot

    ! Initialize stuff
    call global_init(is_serial = .true.) 

    call getopt_init(ierr)
    call getopt_end()

    call parser_init()

    call messages_init()

    call messages_experimental('oct-conductivity')

    call io_init()

    call unit_system_init(global_namespace)

    call spectrum_init(spectrum, global_namespace, default_energy_step = CNST(0.0001), default_max_energy  = CNST(1.0))
 
    !%Variable ConductivitySpectrumTimeStepFactor
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-conductivity_spectrum
    !%Description
    !% In the calculation of the conductivity, it is not necessary
    !% to read the velocity at every time step. This variable controls
    !% the integer factor between the simulation time step and the
    !% time step used to calculate the conductivity.
    !%End

    call messages_obsolete_variable(global_namespace, 'PropagationSpectrumTimeStepFactor', 'ConductivitySpectrumTimeStepFactor')
    call parse_variable(global_namespace, 'ConductivitySpectrumTimeStepFactor', 1, skip)
    if(skip <= 0) call messages_input_error(global_namespace, 'ConductivitySpectrumTimeStepFactor')

    !%Variable ConductivityFromForces
    !%Type logical
    !%Default no
    !%Section Utilities::oct-conductivity_spectrum
    !%Description
    !% (Experimental) If enabled, Octopus will attempt to calculate the conductivity from the forces instead of the current. 
    !%End
    call parse_variable(global_namespace, 'ConductivityFromForces', .false., from_forces)
    if(from_forces) call messages_experimental('ConductivityFromForces')
    
    max_freq = spectrum_nenergy_steps(spectrum)
    
    if (spectrum%end_time < M_ZERO) spectrum%end_time = huge(spectrum%end_time)

    call space_init(space, global_namespace)
    geo => geometry_t(global_namespace, space)

    !We need the total charge
    call parse_variable(global_namespace, 'ExcessCharge', M_ZERO, excess_charge)
    qtot = -(geometry_val_charge(geo) + excess_charge)

    if(from_forces) then

      call messages_write('Info: Reading coordinates from td.general/coordinates')
      call messages_info()

      ! Opens the coordinates files.
      iunit = io_open('td.general/coordinates', global_namespace, action='read')

      call io_skip_header(iunit)
      call spectrum_count_time_steps(global_namespace, iunit, ntime, deltat)
      call io_close(iunit)

      nvel = geo%natoms*space%dim

      SAFE_ALLOCATE(time(1:ntime))
      SAFE_ALLOCATE(velocities(1:nvel, 1:ntime))

      ! Opens the coordinates files.
      iunit = io_open('td.general/coordinates', global_namespace, action='read', status='old', die=.false.)

      call io_skip_header(iunit)

      ntime = 1
      iter = 1

      do
        read(unit = iunit, iostat = ierr, fmt = *) read_iter, curtime, &
          ((geo%atom(ii)%x(jj), jj = 1, 3), ii = 1, geo%natoms), &
          ((geo%atom(ii)%v(jj), jj = 1, 3), ii = 1, geo%natoms)

        curtime = units_to_atomic(units_out%time, curtime)

        if(ierr /= 0 .or. curtime >= spectrum%end_time) then
          iter = iter - 1 ! last iteration is not valid
          ntime = ntime - 1
          exit
        end if

        ASSERT(iter == read_iter + 1)

        if (curtime >= spectrum%start_time .and. mod(iter, skip) == 0) then

          time(ntime) = curtime
          ivel = 1
          do ii = 1, geo%natoms
            do jj = 1, space%dim
              velocities(ivel, ntime) = units_to_atomic(units_out%velocity, geo%atom(ii)%v(jj))
              ivel = ivel + 1
            end do
          end do

          ntime = ntime + 1 
        end if

        iter = iter + 1 !counts number of timesteps (with time larger than zero up to SpecEndTime)
      end do

      call io_close(iunit)

      call messages_write('      done.')
      call messages_info()

    else !from_forces

      vel0(:) = M_ONE

      call messages_write('Info: Reading total current from td.general/total_current')
      call messages_info()

      iunit = io_open('td.general/total_current', global_namespace, action='read')
      
      if(iunit > 0) then
        
        call io_skip_header(iunit)
        call spectrum_count_time_steps(global_namespace, iunit, ntime, deltat)
        ntime = ntime + 1

        call io_skip_header(iunit)

        SAFE_ALLOCATE(total_current(1:space%dim, 1:ntime))
        SAFE_ALLOCATE(heat_current(1:space%dim, 1:ntime))
        SAFE_ALLOCATE(time(1:ntime))
        
        do iter = 1, ntime
          read(iunit, *) read_iter, time(iter), (total_current(ii, iter), ii = 1, space%dim)
        end do
        
        call io_close(iunit)

        if(parse_is_defined(global_namespace, 'TransientAbsorptionReference')) then
          call parse_variable(global_namespace, 'TransientAbsorptionReference', '.', ref_filename)
          ref_file = io_open(trim(ref_filename)//'/total_current', action='read', status='old', die=.false.)
          if(ref_file < 0) then
            message(1) = "Cannot open reference file '"//trim(io_workpath(trim(ref_filename)//'/total_current'))//"'"
            call messages_fatal(1)
          end if
          call io_skip_header(ref_file)
          call spectrum_count_time_steps(global_namespace, ref_file, time_steps_ref, dt_ref)
          if(time_steps_ref < ntime) then
            message(1) = "The reference calculation does not contain enought time steps"
            call messages_fatal(1)
          end if

          if(dt_ref /= time(2)-time(1)) then
            message(1) = "The time step of the reference calculation is different from the current calculation"
            call messages_fatal(1)
          end if

          !We remove the reference
          time_steps_ref = time_steps_ref + 1
          SAFE_ALLOCATE(current_ref(1:space%dim, 1:time_steps_ref))
          call io_skip_header(ref_file)
          do ii = 1, time_steps_ref
            read(ref_file, *) jj, tt, (current_ref(kk, ii), kk = 1, space%dim)
          end do
          call io_close(ref_file)
          do iter = 1, ntime
            do kk = 1, space%dim
              total_current(kk, iter) = total_current(kk, iter) - current_ref(kk, iter)
            end do
          end do
          SAFE_DEALLOCATE_A(current_ref)

          start_time = spectrum%start_time
          call parse_variable(global_namespace, 'GaugeFieldDelay', start_time, spectrum%start_time )
        end if

        
      else
        
        call messages_write("Cannot find the 'td.general/total_current' file.")
        call messages_write("Conductivity will only be calculated from the forces")
        call messages_warning()

        ntime = 1
        SAFE_ALLOCATE(total_current(1:space%dim, 1:ntime))
        SAFE_ALLOCATE(heat_current(1:space%dim, 1:ntime))
        SAFE_ALLOCATE(time(1:ntime))
        
        total_current(1:space%dim, 1:ntime) = M_ZERO
        
      end if

      iunit = io_open('td.general/total_heat_current', global_namespace, action='read', status='old', die=.false.)
      
      if(iunit > 0) then

        if(ntime == 1) then
          call io_skip_header(iunit)
          call spectrum_count_time_steps(global_namespace, iunit, ntime, deltat)
          ntime = ntime + 1  
        end if
        
        call io_skip_header(iunit)
        
        do iter = 1, ntime
          read(iunit, *) read_iter, time(iter), (heat_current(ii, iter), ii = 1, space%dim)
        end do
        
       call io_close(iunit)
       
      else
        
        call messages_write("Cannot find the 'td.general/heat_current' file.")
        call messages_write("Thermal conductivity will only be calculated from the forces")
        call messages_warning()
        
        heat_current(1:space%dim, 1:ntime) = M_ZERO
        
      end if
   end if
   

   SAFE_ALLOCATE(curr(ntime, 1:space%dim))
   SAFE_ALLOCATE(heatcurr(ntime, 1:space%dim))
   integral = M_ZERO

   if(from_forces) iunit = io_open('td.general/current_from_forces', global_namespace, action='write')

   do iter = 1, ntime

     if(from_forces) then
       if(iter == 1) then
         vel0 = M_ZERO
         ivel = 1
         do ii = 1, geo%natoms
           do jj = 1, space%dim
             vel0(jj) = vel0(jj) + velocities(ivel, iter)/TOFLOAT(geo%natoms)
             ivel = ivel + 1
           end do
         end do
       end if
        
       velcm = M_ZERO
       current = M_ZERO
        
       ivel = 1
       do ii = 1, geo%natoms
         do jj = 1, space%dim
           velcm(jj) = velcm(jj) + velocities(ivel, iter)/TOFLOAT(geo%natoms)
           current(jj) = current(jj) + species_mass(geo%atom(ii)%species)/geo%latt%rcell_volume*(velocities(ivel, iter) - vel0(jj))
           ivel = ivel + 1
         end do
       end do
        
       integral(1) = integral(1) + deltat/vel0(1)*(vel0(1)*qtot/geo%latt%rcell_volume + current(1))
       integral(2) = integral(2) + &
         deltat/vel0(1)*(vel0(1)*qtot/geo%latt%rcell_volume - total_current(1, iter)/geo%latt%rcell_volume)

       curr(iter, 1:space%dim) = vel0(1:space%dim)*qtot/geo%latt%rcell_volume + current(1:space%dim)
       
     else
       curr(iter, 1:space%dim)    = total_current(1:space%dim, iter)/geo%latt%rcell_volume
       heatcurr(iter,1:space%dim) = heat_current(1:space%dim, iter)/geo%latt%rcell_volume
     end if
        
     if(from_forces) write(iunit,*) iter, iter*deltat,  curr(iter, 1:space%dim)
      
   end do

   if(from_forces) call io_close(iunit)

   ! Find out the iteration numbers corresponding to the time limits.
   call spectrum_fix_time_limits(spectrum, ntime, deltat, istart, iend, max_freq)
   istart = max(1, istart)
   energy_steps = spectrum_nenergy_steps(spectrum)

   SAFE_ALLOCATE(ftcurr(1:energy_steps, 1:space%dim, 1:2))

   call batch_init(currb, 1, 1, space%dim, curr)
   call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart, iend, spectrum%start_time, deltat, currb)

   call batch_init(ftcurrb, 1, 1, space%dim, ftcurr(:, :, 1))
   call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
      istart, iend, spectrum%start_time, deltat, currb, spectrum%min_energy, spectrum%max_energy, &
      spectrum%energy_step, ftcurrb)
   call ftcurrb%end()

   call batch_init(ftcurrb, 1, 1, space%dim, ftcurr(:, :, 2))
   call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
     istart, iend, spectrum%start_time, deltat, currb, spectrum%min_energy, spectrum%max_energy, &
     spectrum%energy_step, ftcurrb)
      
   call ftcurrb%end()
   call currb%end()


   !and print the spectrum
   iunit = io_open('td.general/conductivity', global_namespace, action='write')

    
   write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'
    write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
    write(unit = iunit, iostat = ierr, fmt = '(a,a,a)') &
      '#  Energy [', trim(units_abbrev(units_out%energy)), '] Conductivity [a.u.] ReX ImX ReY ImY ReZ ImZ'
    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'

    v0 = sqrt(sum(vel0(1:space%dim)**2))
    if(.not. from_forces .or. v0 < epsilon(v0)) v0 = CNST(1.0)

    if( .not. from_forces .and. parse_is_defined(global_namespace, 'GaugeVectorField')) then
      if(parse_block(global_namespace, 'GaugeVectorField', blk) == 0) then
        do ii = 1, space%dim
          call parse_block_float(blk, 0, ii - 1, vel0(ii))
        end do
        call parse_block_end(blk)
      end if
      vel0(1:space%dim) = vel0(1:space%dim) / P_C
      v0 = sqrt(sum(vel0(1:space%dim)**2))
    end if


    do ifreq = 1, energy_steps
      ww = spectrum%energy_step*(ifreq - 1) + spectrum%min_energy
      write(unit = iunit, iostat = ierr, fmt = '(7e20.10)') units_from_atomic(units_out%energy, ww), &
           transpose(ftcurr(ifreq, 1:space%dim, 1:2)/v0)
    end do
    
    call io_close(iunit)
    

    !Compute the inverse dielectric function from the conductivity
    ! We have \chi = -i \sigma / \omega
    ! and \epsilon^-1 = 1 + 4 \pi \chi
    SAFE_ALLOCATE(invdielectric(1:space%dim, 1:energy_steps))
    do ifreq = 1, energy_steps
      ww = max((ifreq-1)*spectrum%energy_step + spectrum%min_energy, M_EPSILON)

      invdielectric(1:space%dim, ifreq) = (vel0(1:space%dim) + M_FOUR * M_PI * &
                  TOCMPLX(ftcurr(ifreq, 1:space%dim, 2),-ftcurr(ifreq, 1:space%dim, 1)) / ww )/v0

    end do

    out_file = io_open('td.general/inverse_dielectric_function_from_current', global_namespace, action='write')
    select case(space%dim)
    case(1)
      write(header, '(3a15)') '#        energy', 'Re x', 'Im x'
    case(2)
      write(header, '(5a15)') '#        energy', 'Re x', 'Im x', 'Re y', 'Im y'
    case(3)
      write(header, '(7a15)') '#        energy', 'Re x', 'Im x', 'Re y', 'Im y', 'Re z', 'Im z'
    end select
    write(out_file,'(a)') trim(header)
    do ifreq = 1, energy_steps
      ww = (ifreq-1)*spectrum%energy_step + spectrum%min_energy
      select case(space%dim)
      case(1)
        write(out_file, '(3e15.6)') ww,                                        &
           TOFLOAT(invdielectric(1, ifreq)), aimag(invdielectric(1, ifreq))
      case(2)
        write(out_file, '(5e15.6)') ww,                                        &
           TOFLOAT(invdielectric(1, ifreq)), aimag(invdielectric(1, ifreq)), &
           TOFLOAT(invdielectric(2, ifreq)), aimag(invdielectric(2, ifreq))
      case(3)
        write(out_file, '(7e15.6)') ww,                                        &
           TOFLOAT(invdielectric(1, ifreq)), aimag(invdielectric(1, ifreq)), &
           TOFLOAT(invdielectric(2, ifreq)), aimag(invdielectric(2, ifreq)), &
           TOFLOAT(invdielectric(3, ifreq)), aimag(invdielectric(3, ifreq))
      end select
    end do
    call io_close(out_file)

    SAFE_DEALLOCATE_A(ftcurr)
    SAFE_DEALLOCATE_A(invdielectric)

    
!!!!!!!!!!!!!!!!!!!!!
    SAFE_ALLOCATE(ftheatcurr(1:energy_steps, 1:space%dim, 1:2))

    ftheatcurr = M_ONE

    call batch_init(heatcurrb, 1, 1, space%dim, heatcurr)
    
    call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, 1, ntime, M_ZERO, deltat, heatcurrb)

    call batch_init(ftheatcurrb, 1, 1, space%dim, ftheatcurr(:, :, 1))
    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
      1, ntime, M_ZERO, deltat, heatcurrb, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, ftheatcurrb)
    call ftheatcurrb%end()

    call batch_init(ftheatcurrb, 1, 1, space%dim, ftheatcurr(:, :, 2))
    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
      1, ntime, M_ZERO, deltat, heatcurrb, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, ftheatcurrb)
    call ftheatcurrb%end()
    
    call heatcurrb%end()


    !and print the spectrum
    iunit = io_open('td.general/heat_conductivity', global_namespace, action='write')

    
    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'
    write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
    write(unit = iunit, iostat = ierr, fmt = '(a,a,a)') &
      '#  Energy [', trim(units_abbrev(units_out%energy)), '] Conductivity [a.u.] ReX ImX ReY ImY ReZ ImZ'
    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'

    v0 = sqrt(sum(vel0(1:space%dim)**2))
    if(.not. from_forces .or. v0 < epsilon(v0)) v0 = CNST(1.0)

    do ifreq = 1, energy_steps
      ww = spectrum%energy_step*(ifreq - 1) + spectrum%min_energy
      write(unit = iunit, iostat = ierr, fmt = '(7e20.10)') units_from_atomic(units_out%energy, ww), &
        transpose(ftheatcurr(ifreq, 1:space%dim, 1:2)/v0)
      !print *, ifreq, ftheatcurr(ifreq, 1:3, 1:2)
   end do
    
    call io_close(iunit)
    
    SAFE_DEALLOCATE_A(ftheatcurr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    SAFE_DEALLOCATE_P(geo)

    SAFE_DEALLOCATE_A(time)

    call io_end()
    call messages_end()

    call parser_end()
    call global_end()

  end program conductivity

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
