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
  use string_m
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
  integer :: time_steps, mode
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
        power = - ( power*dt / (-kick%delta_strength) ) / (dt*time_steps)
      case(SECOND_ORDER)
        do j = 0, time_steps
          x = cos(omega*j*dt)
          power = power + x*ot(j)
        end do
        power =  - ( power*dt / (-kick%delta_strength**2) ) / (dt*time_steps)
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

  integer :: in_file(3), i, lmax, ierr, nspin, order, nfrequencies, print_omega_file, run_mode
  FLOAT :: omega, power, leftbound, rightbound, search_interval, &
           total_time, dw, w, aw, min_aw, min_w, omega_orig, final_time
  FLOAT, allocatable :: warray(:), tarray(:)
  character(len=100) :: file1, file2

  integer, parameter :: ANALYZE_SIGNAL           = 1, &
                        GENERATE_NTHORDER_SIGNAL = 2

  ! Reads the information passed through the command line options (if available).
  call getopt_init(ierr)
  if(ierr.ne.0) then
    message(1) = "Your fortran compiler doesn't support command line arguments,"
    message(2) = "the oct-oscillator-strength command is not available."
    call write_fatal(2)
  end if

  ! Set the default values
  run_mode = ANALYZE_SIGNAL
  file1 = "ot"
  file2 = ""
  omega = - M_ONE
  search_interval = - M_ONE
  order = 1
  nfrequencies = 1000
  final_time = - M_ONE
  print_omega_file = 0

  ! Get the parameters from the command line.
  call getopt_oscillator_strength(run_mode, omega, file1, file2, search_interval, &
                                  order, nfrequencies, final_time, print_omega_file)

  ! Initialize stuff
  call global_init()
  in_debug_mode = .false.
  call datasets_init(1)
  call io_init(defaults = .true.)
  if(in_debug_mode) then
     call io_mkdir('debug')
  end if


  select case(run_mode)
  case(GENERATE_NTHORDER_SIGNAL)
    call generate_signal(order)
    stop
  case(ANALYZE_SIGNAL)
  case default
  end select


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

  if(trim(file1).eq.'ot') then
    call read_ot()
  else
    call read_multipoles()
  end if


  if(omega > M_ZERO) then
    omega = omega * units_inp%energy%factor
  else
    omega = M_HALF
  end if
  if(search_interval > M_ZERO) then
    search_interval = search_interval * units_inp%energy%factor
  else
    search_interval = M_HALF
  end if
  if(final_time > M_ZERO) then
    total_time = final_time * units_inp%time%factor
    if(total_time > dt*time_steps) then
      total_time = dt*time_steps
      write(0, '(a)')        '* WARNING: The reqeusted total time to process is larger than the time'
      write(0, '(a)')        '*          available in the input file.'
      write(0, '(a,f8.4,a)') '           The time has been adjusted to ', total_time / units_inp%time%factor, &
                           units_inp%time%abbrev
    end if
    time_steps = int(total_time / dt)
    total_time = time_steps * dt
  else
    total_time = dt*time_steps
  end if
  leftbound = omega - search_interval
  rightbound = omega + search_interval

  ALLOCATE(warray(nfrequencies), nfrequencies)
  ALLOCATE(tarray(nfrequencies), nfrequencies)
  dw = (rightbound-leftbound) / (nfrequencies - 1)
  min_w = omega
  min_aw = M_ZERO
  do i = 1, nfrequencies
    w = leftbound + (i-1)*dw
    warray(i) = w
    call ft(w, aw)
    tarray(i) = -aw
    if(aw < min_aw) then
      min_aw = aw
      min_w = w
    end if
  end do

  if(print_omega_file .eq. 1) call print_file()

  omega_orig = omega
  omega = min_w
  call loct_1dminimize(min_w - 2*dw, min_w + 2*dw, omega, ft, ierr)
  if(ierr.ne.0) then
    write(0,'(a)') 'Could not find a maximum.'
    write(0,'(a)')
    write(0, '(a,f12.8,a,f12.8,a)') '   Search interval = [', leftbound / units_inp%energy%factor, ',', &
                                                              rightbound / units_inp%energy%factor, ']'
    write(0, '(a,f12.4,a)')         '   Search discretization = ', dw / units_inp%energy%factor, &
                                    ' '//trim(units_inp%energy%abbrev)
    stop
  end if

  call ft(omega, power)
  power = -power
  power = power / (M_ONE - sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))

!!$  call write_ot()
  call write_ot(nspin, time_steps, dt, kick, units_inp, ot)

  select case(order)
  case(1)
    write(*, '(a,f12.8,a,f12.8,a)') 'omega    = ', omega / units_inp%energy%factor, &
                                    ' '//trim(units_inp%energy%abbrev)//' = ',      &
                                    omega, ' Ha'
    write(*, '(a,f12.8,a,f12.8,a)') 'C(omega) = ', power / units_inp%length%factor**2, &
                                    ' '//trim(units_inp%length%abbrev)//'^2 =',        &
                                    power, ' b^2'
    write(*, '(a,f12.8,a,f12.8,a)') '<0|P|I>  = ', sqrt(abs(power)) / units_inp%length%factor, &
                                    ' '//trim(units_inp%length%abbrev)//' = ',                 &
                                    sqrt(abs(power)),' b'
    write(*, '(a,f12.8)')           'f[O->I]  = ', M_TWO*omega*power
    write(*, '(a)')
    write(*, '(a,f12.8,a,f12.8,a)') '   Search interval = [', leftbound / units_inp%energy%factor, ',', &
                                                              rightbound / units_inp%energy%factor, ']'
    write(*, '(a,f12.4,a)')         '   Search discretization = ', dw / units_inp%energy%factor, &
                                    ' '//trim(units_inp%energy%abbrev)
  case(2)
    write(*, '(a,f12.8,a,f12.8,a)') 'omega    = ', omega / units_inp%energy%factor, &
                                    ' '//trim(units_inp%energy%abbrev)//' = ',      &
                                    omega, ' Ha'
    write(*, '(a,f12.8,a,f12.8,a)') 'C(omega) = ', power / units_inp%length%factor**3, &
                                    ' '//trim(units_inp%length%abbrev)//'^3 = ',       &
                                    power, ' b^3'
    write(*, '(a)')
    write(*, '(a,f12.8,a,f12.8,a)') '   Search interval = [', leftbound / units_inp%energy%factor, ',', &
                                                              rightbound / units_inp%energy%factor, ']'
  end select

  deallocate(ot)
  call io_end()
  call datasets_end()
  call global_end()

  contains


    ! ---------------------------------------------------------
    subroutine read_multipoles
      FLOAT :: dump
      integer :: j, isp
      FLOAT, allocatable :: dipole(:, :, :)

      in_file(1) = io_open(trim(file1), action='read', status='old', die=.false.)
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
        if(in_file(2) >= 0) then
        else
          write(message(1),'(a)') 'File '//trim(file2)//' not found.'
          call write_fatal(1)
        end if
      end if

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
  
        dipole(1:3, i, :) = dipole(1:3, i, :) * units_inp%length%factor
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
    
          dipole(1:3, i, :) = dipole(1:3, i, :) * units_inp%length%factor
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
    end subroutine read_multipoles


    ! ---------------------------------------------------------
    subroutine read_ot
      integer :: iunit, i
      character(len=100) :: line
      FLOAT :: dummy, t1, t2

      iunit = io_open('ot', action='read', status='old')
      if(iunit .eq. 0) then
        write(message(1),'(a)') 'A file called ot should be present and was not found.'
        call write_fatal(1)
      end if

      read(iunit, '(15x,i2)') nspin
      call kick_read(kick, iunit)
      read(iunit, '(a)')  line
      read(iunit, '(a)')  line
      call io_skip_header(iunit)

      ! Figure out about the units of the file
      i = index(line,'eV')
      if(i.ne.0) then
        call units_get(units_inp, 'eVA')
      else
        call units_get(units_inp, 'a.u')
      end if
      ! This should not be like this.
      call units_get(units_out, 'eVA')

      ! count number of time_steps
      time_steps = 0
      do
        read(iunit, *, end=100) dummy
        time_steps = time_steps + 1
        if(time_steps == 1) t1 = dummy
        if(time_steps == 2) t2 = dummy
      end do
      100 continue
      dt = (t2 - t1) * units_inp%time%factor ! units_out is OK

      call io_skip_header(iunit)

      ALLOCATE(ot(0:time_steps), time_steps+1)

      do i = 0, time_steps-1
        read(iunit, *) dummy, ot(i)
        ot(i) = ot(i) * units_inp%length%factor
      end do

    end subroutine read_ot


    ! ---------------------------------------------------------
    subroutine print_file
      integer :: iunit, i
      character(len=20) :: header_string

      iunit = io_open('omega', action='write', status='replace')

      write(iunit, '(a15,i2)')      '# nspin        ', nspin
      call kick_write(kick, iunit)
      write(iunit, '(a)') '#%'
      write(iunit, '(a1,a20)', advance = 'no') '#', str_center("omega", 20)
      write(header_string,'(a)') 'F(omega)'
      write(iunit, '(a20)', advance = 'yes') str_center(trim(header_string), 20)
      write(iunit, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_inp%energy%abbrev) // ']', 20)
      ! Here we should print the units of the transform.
      write(iunit, '(a)', advance = 'yes')

      do i = 1, nfrequencies
        write(iunit,'(2e20.8)') warray(i) / units_inp%energy%factor, &
                                tarray(i)
      end do

      call io_close(iunit)
    end subroutine print_file

end program oscillator_strength


subroutine generate_signal(order)
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use spectrum_m
  use lib_adv_alg_m

  implicit none

  integer, intent(in) :: order

  logical :: file_exists
  integer :: i, j, nspin, time_steps, lmax, nfiles, k, isp
  integer, allocatable :: iunit(:)
  FLOAT :: dt, lambda, det, dump, o0
  FLOAT, allocatable :: q(:), mu(:), qq(:, :), c(:)
  character(len=20) :: filename
  type(kick_t) :: kick
  type(unit_system_t) :: units
  FLOAT, allocatable :: multipole(:, :, :), ot(:)

  ! Find out how many files do we have
  nfiles = 0
  do
    write(filename,'(a11,i1)') 'multipoles.', nfiles+1
    inquire(file=trim(filename), exist  = file_exists)
    if(.not.file_exists) exit
    nfiles = nfiles + 1
  end do

  ! WARNING: Check that order is smaller or equal to nfiles
  if(nfiles == 0) then
    write(message(1),'(a)') 'No multipoles.x file was found'
    call write_fatal(1)
  endif
  if(order > nfiles) then
    write(message(1),'(a)') 'The order that you ask for is higher than the number'
    write(message(2),'(a)') 'of multipoles.x file that you supply.'
    call write_fatal(2)
  end if

  ! Open the files.
  allocate(iunit(nfiles))
  do j = 1, nfiles
    write(filename,'(a11,i1)') 'multipoles.', j
    iunit(j) = io_open(trim(filename), action='read', status='old', die=.false.)
  end do

  ALLOCATE(q(nfiles), nfiles)
  ALLOCATE(mu(nfiles), nfiles)
  ALLOCATE(qq(nfiles, nfiles), nfiles*nfiles)
  ALLOCATE(c(nfiles), nfiles)

  c        = M_ZERO
  c(order) = M_ONE

  ! Get the basic info from the first file
  call spectrum_mult_info(iunit(1), nspin, kick, time_steps, dt, units, lmax=lmax)
  lambda = kick%delta_strength
  q(1) = 1.0_8

  do j = 2, nfiles
    call spectrum_mult_info(iunit(j), nspin, kick, time_steps, dt, units, lmax=lmax)
    q(j) = kick%delta_strength / lambda
  end do

  do i = 1, nfiles
   do j = 1, nfiles
     qq(i, j) = q(j)**i
   end do
  end do

  det = lalg_inverter(nfiles, qq, invert = .true.)

  mu = matmul(qq, c)

  ALLOCATE(multipole(3, 0:time_steps, nspin), 2*(time_steps+1)*nspin)
  ALLOCATE(ot(0:time_steps), time_steps+1)
  multipole = M_ZERO
  ot = M_ZERO

  do j = 1, nfiles
    call io_skip_header(iunit(j))

    do i = 0, time_steps
      select case(nspin)
      case(1)
        read(iunit(j), *) k, dump, dump, multipole(1:3, i, 1)
      case(2)
        read(iunit(j), *) k, dump, dump, multipole(1:3, i, 1), dump, multipole(1:3, i, 2)
      case(4)
        read(iunit(j), *) k, dump, dump, multipole(1:3, i, 1), dump, multipole(1:3, i, 2), &
          dump, multipole(1:3, i, 3), dump, multipole(1:3, i, 4)
      end select
      multipole(1:3, i, :) = multipole(1:3, i, :) * units%length%factor

      dump = M_ZERO
      do isp = 1, nspin
        dump = dump + sum(multipole(1:3, i, isp) * kick%pol(1:3,kick%pol_dir))
      end do
      if(i == 0) o0 = dump

      ot(i) = ot(i) + mu(j)*(dump - o0)
    end do

  end do

  call write_ot(nspin, time_steps, dt, kick, units, ot(0:time_steps))

  ! Close files and exit.
  do j = 1, nfiles
    call io_close(iunit(j))
  end do
  deallocate(iunit, q, mu, qq, c, ot, multipole)
end subroutine generate_signal


! ---------------------------------------------------------
subroutine write_ot(nspin, time_steps, dt, kick, units, ot)
  use global_m
  use io_m
  use units_m
  use spectrum_m
  use string_m

  implicit none

  integer,             intent(in) :: nspin, time_steps
  FLOAT,               intent(in) :: dt
  type(kick_t),        intent(in) :: kick
  type(unit_system_t), intent(in) :: units
  FLOAT,               intent(in) :: ot(0:time_steps)

  integer :: iunit, i
  character(len=20) :: header_string

  iunit = io_open('ot', action='write', status='replace')

  write(iunit, '(a15,i2)')      '# nspin        ', nspin
  call kick_write(kick, iunit)

  ! Units
  write(iunit,'(a1)', advance = 'no') '#'
  write(iunit,'(a20)', advance = 'no') str_center('t', 19)
  write(iunit,'(a20)', advance = 'yes') str_center('<O>(t)', 20)
  write(iunit,'(a1)', advance = 'no') '#'
  write(header_string, '(a)') '['//trim(units%time%abbrev)//']'
  write(iunit,'(a20)', advance = 'no')  str_center(trim(header_string), 19)
  write(header_string, '(a)') '['//trim(units%length%abbrev)//']'
  write(iunit,'(a20)', advance = 'yes')  str_center(trim(header_string), 20)

  do i = 0, time_steps
    write(iunit, '(2e20.8)') i*dt / units%time%factor, ot(i) / units%length%factor
  end do

  call io_close(iunit)

end subroutine write_ot


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
