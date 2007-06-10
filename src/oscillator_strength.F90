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
!! $Id: cross_section.F90 2781 2007-03-23 10:58:32Z lorenzen $

#include "global.h"

program oscillator_strength
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use spectrum_m
  use lib_oct_m
  use getopt_m

  implicit none

  integer :: in_file(3), out_file(3), i, eq_axis, nspin, lmax, time_steps
  logical :: calculate_tensor
  type(spec_t) :: s
  integer :: ierr
  FLOAT :: omega, power, leftbound, rightbound, search_interval, total_time, dt
  character(len=100) :: filename


  ! Reads the information passed through the command line options (if available).
  call getopt_init(ierr)
  if(ierr.ne.0) then
    message(1) = "Your fortran compiler doesn't support command line arguments,"
    message(2) = "the oct-oscillator-strength command is not available."
    call write_fatal(2)
  end if
  filename = "multipoles"
  omega = CNST(0.1)
  search_interval = CNST(0.01)
  call getopt_oscillator_strength(omega, filename, search_interval)

  ! Initialize stuff
  call global_init()
  call parser_init()
  in_debug_mode = .false.
  call datasets_init(1)
  call io_init()
  if(in_debug_mode) then
     call io_mkdir('debug')
  end if

  s%damp        = SPECTRUM_DAMP_NONE
  s%transform   = SPECTRUM_TRANSFORM_SIN
  s%start_time  = M_ZERO
  s%end_time    = -M_ONE
  s%energy_step = CNST(0.001)
  s%max_energy  = CNST(1.0)
  s%damp_factor = M_ZERO

  call read_files(trim(filename))

  omega = omega * units_inp%energy%factor
  search_interval = search_interval * units_inp%energy%factor
  leftbound = omega - search_interval
  rightbound = omega + search_interval

  call loct_1dminimize(leftbound, rightbound, omega, func, ierr)
  if(ierr.ne.0) then
    write(0,'(a)') 'Could not find a maximum.'
    stop
  end if

  call ft(in_file(1), s, omega, power)
  total_time = dt*time_steps
  power = power / (M_ONE - sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))

  write(*, '(a,f14.8,a)') 'omega    = ', omega / units_inp%energy%factor, &
                          ' '//trim(units_inp%energy%abbrev)
  write(*, '(a,f14.8,a)') 'C(omega) = ', power / units_out%length%factor**2, &
                          ' '//trim(units_inp%length%abbrev)//'^2'
  write(*, '(a,f14.8,a)') '<0|P|I>  = ', sqrt(abs(power)) / units_inp%length%factor, &
                          ' '//trim(units_inp%length%abbrev)
  write(*, '(a,f14.8)')   'f[O->I]  = ', M_TWO*omega*power

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()

  contains

    subroutine func(x, fx)
      FLOAT, intent(in) :: x
      FLOAT, intent(out) :: fx
      call ft(in_file(1), s, x, fx)
      fx = -fx
    end subroutine func

    subroutine read_files(fname)
      character(len=*), intent(in) :: fname
      type(kick_t) :: kick
      
      in_file(1) = io_open(trim(fname), action='read', status='old', die=.false.)
      if(in_file(1) < 0) in_file(1) = io_open('td.general/'//trim(fname), action='read', status='old', die=.false.)
      if(in_file(1) >= 0) then
        ! OK, so we only have one file. Now we have to see what it has inside.
        call spectrum_mult_info(in_file(1), nspin, kick, time_steps, dt, units_inp, lmax=lmax)
        ! This should not be like this.
        call units_get(units_out, 'eVA')
      else
        write(message(1),'(a)') 'File '//trim(fname)//' not found.'
        call write_fatal(1)
      end if

    end subroutine read_files

    subroutine ft(in_file, s, omega, power)
      integer, intent(in) :: in_file
      type(spec_t),  intent(inout) :: s
      FLOAT, intent(in)   :: omega
      FLOAT, intent(out)  :: power
  
      integer :: nspin, lmax, time_steps, is, ie, ntiter, i, j, jj, isp
      type(kick_t) :: kick
      FLOAT :: dt, dump, w, x
      FLOAT, allocatable :: dipole(:, :, :)
      type(unit_system_t) :: file_units
  
      ! This function gives us back the unit connected to the "multipoles" file, the header information,
      ! the number of time steps, and the time step.
      call spectrum_mult_info(in_file, nspin, kick, time_steps, dt, file_units, lmax=lmax)
  
      ! Find out the iteration numbers corresponding to the time limits.
      is = 0
      ie = time_steps
      ntiter = ie - is + 1
  
      ! Read the dipole.
      call io_skip_header(in_file)
      ALLOCATE(dipole(3, 0:time_steps, nspin), 3*(time_steps+1)*nspin)
      do i = 0, time_steps
        select case(nspin)
        case(1)
          read(in_file, *) j, dump, dump, dipole(1:3, i, 1)
        case(2)
          read(in_file, *) j, dump, dump, dipole(1:3, i, 1), dump, dipole(1:3, i, 2)
        case(4)
          read(in_file, *) j, dump, dump, dipole(1:3, i, 1), dump, dipole(1:3, i, 2), &
            dump, dipole(1:3, i, 3), dump, dipole(1:3, i, 4)
        end select
  
        dipole(1:3, i, :) = dipole(1:3, i, :) * units_out%length%factor
      end do
      ! Now substract the initial dipole.
      do i = time_steps, 0, -1
        dipole(:, i, :) = dipole(:, i, :) - dipole(:, 0, :)
      end do
  
      w = omega
  
      power = M_ZERO
  
      do j = is, ie
        jj = j - is
  
        select case(s%transform)
        case(SPECTRUM_TRANSFORM_SIN)
          x = sin(w*jj*dt)
        case(SPECTRUM_TRANSFORM_COS)
          x = cos(w*jj*dt)
        case(SPECTRUM_TRANSFORM_EXP)
          x = exp(-w*jj*dt)
        end select
  
        do isp = 1, nspin
          power = power + x*sum(dipole(1:3, j, isp) * kick%pol(1:3,kick%pol_dir))
        end do
  
      end do
  
      power = ( power*dt / (-kick%delta_strength) ) / (dt*ntiter)
      call pop_sub()
    end subroutine ft
    
end program oscillator_strength

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
