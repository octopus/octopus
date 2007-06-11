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

module oscillator_strength_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use spectrum_m
  use lib_oct_m

  implicit none

  integer, parameter ::         &
    FIRST_ORDER_ONE_FILE  = 0,  &
    FIRST_ORDER_TWO_FILES = 1,  &
    SECOND_ORDER          = 2


  FLOAT, allocatable :: ot(:)
  type(kick_t) :: kick
  integer :: time_steps, ntiter, mode
  FLOAT :: dt

  contains

    subroutine ft(omega, power)
      FLOAT, intent(in)   :: omega
      FLOAT, intent(out)  :: power
      integer :: i, j
      FLOAT :: x
      power = M_ZERO

      select case(mode)
      case(FIRST_ORDER_ONE_FILE, FIRST_ORDER_TWO_FILES)
        do j = 0, time_steps
          x = sin(omega*j*dt)
          power = power + x*ot(j)
        end do
        power = - ( power*dt / (-kick%delta_strength) ) / (dt*ntiter)
      case(SECOND_ORDER)
        do j = 0, time_steps
          x = cos(omega*j*dt)
          power = power + x*ot(j)
        end do
        power =  - ( power*dt / (-kick%delta_strength**2) ) / (dt*ntiter)
      end select

    end subroutine ft

end module oscillator_strength_m


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
  use command_line_m
  use oscillator_strength_m

  implicit none

  integer :: in_file(3), i, lmax, ierr, nspin, order
  FLOAT :: omega, power, leftbound, rightbound, search_interval, total_time
  character(len=100) :: file1, file2

  ! Reads the information passed through the command line options (if available).
  call getopt_init(ierr)
  if(ierr.ne.0) then
    message(1) = "Your fortran compiler doesn't support command line arguments,"
    message(2) = "the oct-oscillator-strength command is not available."
    call write_fatal(2)
  end if
  file1 = "multipoles"
  file2 = ""
  omega = CNST(0.1)
  search_interval = CNST(0.01)
  order = 1
  call getopt_oscillator_strength(omega, file1, file2, search_interval, order)

  if(trim(file2).ne."") then
    if(order.eq.1) then
      mode = FIRST_ORDER_TWO_FILES
    else
      mode = SECOND_ORDER
    end if
  elseif(order.eq.2) then
    message(1) = 'In order to get the second order response, the program needs two files.'
    call write_fatal(1)
  end if

  ! Initialize stuff
  call global_init()
  in_debug_mode = .false.
  call datasets_init(1)
  call io_init(defaults = .true.)
  if(in_debug_mode) then
     call io_mkdir('debug')
  end if

  call read_files()

  omega = omega * units_inp%energy%factor
  search_interval = search_interval * units_inp%energy%factor
  leftbound = omega - search_interval
  rightbound = omega + search_interval

  call loct_1dminimize(leftbound, rightbound, omega, ft, ierr)
  if(ierr.ne.0) then
    write(0,'(a)') 'Could not find a maximum.'
    stop
  end if

  call ft(omega, power)
  power = -power
  total_time = dt*time_steps
  power = power / (M_ONE - sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))

  write(*, '(a,f14.8,a)') 'omega    = ', omega / units_inp%energy%factor, &
                          ' '//trim(units_inp%energy%abbrev)
  write(*, '(a,f14.8,a)') 'C(omega) = ', power / units_out%length%factor**2, &
                          ' '//trim(units_inp%length%abbrev)//'^2'
  write(*, '(a,f14.8,a)') '<0|P|I>  = ', sqrt(abs(power)) / units_inp%length%factor, &
                          ' '//trim(units_inp%length%abbrev)
  write(*, '(a,f14.8)')   'f[O->I]  = ', M_TWO*omega*power

  deallocate(ot)
  call io_end()
  call datasets_end()
  call global_end()

  contains

    subroutine read_files
      FLOAT :: dump
      integer :: j, isp
      FLOAT, allocatable :: dipole(:, :, :)
      
      in_file(1) = io_open(trim(file1), action='read', status='old', die=.false.)
      if(in_file(1) < 0) in_file(1) = io_open('td.general/'//trim(file1), action='read', status='old', die=.false.)
      if(in_file(1) >= 0) then
        call spectrum_mult_info(in_file(1), nspin, kick, time_steps, dt, units_inp, lmax=lmax)
        ! This should not be like this.
        call units_get(units_out, 'eVA')
      else
        write(message(1),'(a)') 'File '//trim(file1)//' not found.'
        call write_fatal(1)
      end if

      if( (mode.eq.FIRST_ORDER_TWO_FILES) .or. (mode.eq.SECOND_ORDER) ) then
        in_file(2) = io_open(trim(file2), action='read', status='old', die=.false.)
        if(in_file(2) < 0) in_file(2) = io_open('td.general/'//trim(file2), action='read', status='old', die=.false.)
        if(in_file(2) >= 0) then
        else
          write(message(1),'(a)') 'File '//trim(file2)//' not found.'
          call write_fatal(1)
        end if
      end if

      ! Find out the iteration numbers corresponding to the time limits.
      ntiter = time_steps + 1

      ! Read the dipole.
      ALLOCATE(dipole(3, 0:time_steps, nspin), 3*(time_steps+1)*nspin)
      ALLOCATE(ot(0:time_steps), time_steps+1)

      call io_skip_header(in_file(1))
      do i = 0, time_steps
        select case(nspin)
        case(1)
          read(in_file(1), *) j, dump, dump, dipole(1:3, i, 1)
        case(2)
          read(in_file(1), *) j, dump, dump, dipole(1:3, i, 1), dump, dipole(1:3, i, 2)
        case(4)
          read(in_file(1), *) j, dump, dump, dipole(1:3, i, 1), dump, dipole(1:3, i, 2), &
            dump, dipole(1:3, i, 3), dump, dipole(1:3, i, 4)
        end select
  
        dipole(1:3, i, :) = dipole(1:3, i, :) * units_out%length%factor
      end do
      ! Now substract the initial dipole.
      do i = time_steps, 0, -1
        dipole(:, i, :) = dipole(:, i, :) - dipole(:, 0, :)
      end do
      do i = 0, time_steps
        ot(i) = M_ZERO
        do isp = 1, nspin
          ot(i) = ot(i) + sum(dipole(1:3, i, isp) * kick%pol(1:3,kick%pol_dir))
        end do
      end do

      if( (mode.eq.FIRST_ORDER_TWO_FILES) .or. (mode.eq.SECOND_ORDER) ) then
        call io_skip_header(in_file(2))
        do i = 0, time_steps
          select case(nspin)
          case(1)
            read(in_file(2), *) j, dump, dump, dipole(1:3, i, 1)
          case(2)
            read(in_file(2), *) j, dump, dump, dipole(1:3, i, 1), dump, dipole(1:3, i, 2)
          case(4)
            read(in_file(2), *) j, dump, dump, dipole(1:3, i, 1), dump, dipole(1:3, i, 2), &
              dump, dipole(1:3, i, 3), dump, dipole(1:3, i, 4)
          end select
    
          dipole(1:3, i, :) = dipole(1:3, i, :) * units_out%length%factor
        end do
        ! Now substract the initial dipole.
        do i = time_steps, 0, -1
          dipole(:, i, :) = dipole(:, i, :) - dipole(:, 0, :)
        end do
        do i = 0, time_steps
          do isp = 1, nspin
            if(mode.eq.FIRST_ORDER_TWO_FILES) then
              ot(i) = ot(i) - sum(dipole(1:3, i, isp) * kick%pol(1:3,kick%pol_dir))
            else
              ot(i) = ot(i) + sum(dipole(1:3, i, isp) * kick%pol(1:3,kick%pol_dir))
            end if
          end do
        end do
        ot = M_HALF*ot

      end if

      deallocate(dipole)
    end subroutine read_files

end program oscillator_strength

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
