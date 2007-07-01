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

        do j = 0, time_steps
          x = sin(omega*j*dt)
          power = power + x*ot(j)
        end do
        power = power*dt / (dt*time_steps)

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

  integer :: run_mode, order, nfrequencies, iprint_omega_file, nspin, ierr
  FLOAT :: omega, search_interval, final_time

  integer, parameter :: ANALYZE_NTHORDER_SIGNAL           = 1, &
                        GENERATE_NTHORDER_SIGNAL = 2

  ! Reads the information passed through the command line options (if available).
  call getopt_init(ierr)
  if(ierr.ne.0) then
    message(1) = "Your fortran compiler doesn't support command line arguments,"
    message(2) = "the oct-oscillator-strength command is not available."
    call write_fatal(2)
  end if

  ! Set the default values
  run_mode = ANALYZE_NTHORDER_SIGNAL
  omega = - M_ONE
  search_interval = - M_ONE
  order = 1
  nfrequencies = 1000
  final_time = - M_ONE
  iprint_omega_file = 0

  ! Get the parameters from the command line.
  call getopt_oscillator_strength(run_mode, omega, search_interval, &
                                  order, nfrequencies, final_time, iprint_omega_file)

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
  case(ANALYZE_NTHORDER_SIGNAL)
    call analyze_signal(order, omega, search_interval, final_time, nfrequencies, iprint_omega_file)
  case default
  end select

  call io_end()
  call datasets_end()
  call global_end()
end program oscillator_strength


! ---------------------------------------------------------
subroutine analyze_signal(order, omega, search_interval, final_time, nfrequencies, iprint_omega_file)
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use spectrum_m
  use lib_adv_alg_m
  use oscillator_strength_m

  implicit none

  integer, intent(in)    :: order
  FLOAT,   intent(inout) :: omega
  FLOAT,   intent(inout) :: search_interval
  FLOAT,   intent(inout) :: final_time
  integer, intent(inout) :: nfrequencies
  integer, intent(inout) :: iprint_omega_file

  FLOAT :: total_time, leftbound, rightbound, dw, min_w, min_aw, w, aw, omega_orig, power
  FLOAT, allocatable :: warray(:), tarray(:)
  integer :: nspin, i, ierr
  type(unit_system_t) :: units

  call read_ot(nspin, units)

  if(omega > M_ZERO) then
    omega = omega * units%energy%factor
  else
    omega = M_HALF
  end if

  if(search_interval > M_ZERO) then
    search_interval = search_interval * units%energy%factor
  else
    search_interval = M_HALF
  end if

  if(final_time > M_ZERO) then
    total_time = final_time * units%time%factor
    if(total_time > dt*time_steps) then
      total_time = dt*time_steps
      write(0, '(a)')        '* WARNING: The reqeusted total time to process is larger than the time'
      write(0, '(a)')        '*          available in the input file.'
      write(0, '(a,f8.4,a)') '           The time has been adjusted to ', total_time / units%time%factor, &
                           units%time%abbrev
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
  warray = M_ZERO; tarray = M_ZERO
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

  if(iprint_omega_file .eq. 1) call print_omega_file(units, nspin, nfrequencies, warray, tarray)

  omega_orig = omega
  omega = min_w
  call loct_1dminimize(min_w - 2*dw, min_w + 2*dw, omega, ft, ierr)
  if(ierr.ne.0) then
    write(0,'(a)') 'Could not find a maximum.'
    write(0,'(a)')
    write(0, '(a,f12.8,a,f12.8,a)') '   Search interval = [', leftbound / units%energy%factor, ',', &
                                                              rightbound / units%energy%factor, ']'
    write(0, '(a,f12.4,a)')         '   Search discretization = ', dw / units%energy%factor, &
                                    ' '//trim(units%energy%abbrev)
    stop
  end if

  call ft(omega, power)
  power = -power
  power = power / (M_ONE - sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))

  select case(order)
  case(1)
    write(*, '(a,f12.8,a,f12.8,a)') 'omega    = ', omega / units%energy%factor, &
                                    ' '//trim(units%energy%abbrev)//' = ',      &
                                    omega, ' Ha'
    write(*, '(a,f12.8,a,f12.8,a)') 'C(omega) = ', power / units%length%factor**2, &
                                    ' '//trim(units%length%abbrev)//'^2 =',        &
                                    power, ' b^2'
    write(*, '(a,f12.8,a,f12.8,a)') '<0|P|I>  = ', sqrt(abs(power)) / units%length%factor, &
                                    ' '//trim(units%length%abbrev)//' = ',                 &
                                    sqrt(abs(power)),' b'
    write(*, '(a,f12.8)')           'f[O->I]  = ', M_TWO*omega*power
    write(*, '(a)')
    write(*, '(a,f12.8,a,f12.8,a)') '   Search interval = [', leftbound / units%energy%factor, ',', &
                                                              rightbound / units%energy%factor, ']'
    write(*, '(a,f12.4,a)')         '   Search discretization = ', dw / units%energy%factor, &
                                    ' '//trim(units%energy%abbrev)
  case(2)
    write(*, '(a,f12.8,a,f12.8,a)') 'omega    = ', omega / units%energy%factor, &
                                    ' '//trim(units%energy%abbrev)//' = ',      &
                                    omega, ' Ha'
    write(*, '(a,f12.8,a,f12.8,a)') 'C(omega) = ', power / units%length%factor**3, &
                                    ' '//trim(units%length%abbrev)//'^3 = ',       &
                                    power, ' b^3'
    write(*, '(a)')
    write(*, '(a,f12.8,a,f12.8,a)') '   Search interval = [', leftbound / units%energy%factor, ',', &
                                                              rightbound / units%energy%factor, ']'
  end select

  call modify_ot(nspin, time_steps, dt, kick, units, ot, omega, sqrt(abs(power)))

  deallocate(ot)

end subroutine analyze_signal


! ---------------------------------------------------------
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
  lambda = abs(kick%delta_strength)
  q(1) = kick%delta_strength / lambda

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

  ot = ot / lambda**order

  call write_ot(nspin, time_steps, dt, kick, units, ot(0:time_steps))

  ! Close files and exit.
  do j = 1, nfiles
    call io_close(iunit(j))
  end do
  deallocate(iunit, q, mu, qq, c, ot, multipole)
end subroutine generate_signal


! ---------------------------------------------------------
subroutine modify_ot(nspin, time_steps, dt, kick, units, ot, omega, cij)
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
  FLOAT,               intent(in) :: omega
  FLOAT,               intent(in) :: cij

  integer :: iunit, i
  character(len=20) :: header_string
  FLOAT :: modot

  iunit = io_open('ot.mod', action='write', status='replace')

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

  write(0, *) 'A', omega, (abs(cij)**2), kick%delta_strength

  do i = 0, time_steps
    modot = ot(i) + M_TWO * abs(cij)**2 * sin(omega*i*dt)
    write(iunit, '(3e20.8)') i*dt / units%time%factor, modot / units%length%factor, ot(i) / units%length%factor
  end do

  call io_close(iunit)
end subroutine modify_ot


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


! ---------------------------------------------------------
subroutine read_ot(nspin, units)!, ot)
  use global_m
  use io_m
  use messages_m
  use datasets_m
  use units_m
  use spectrum_m
  use string_m
  use oscillator_strength_m

  implicit none

  integer, intent(out) :: nspin

  type(unit_system_t), intent(inout) :: units

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
  call units_get(units, 'eVA')

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
  100 continue
  dt = (t2 - t1) * units%time%factor ! units_out is OK

  call io_skip_header(iunit)

  ALLOCATE(ot(0:time_steps), time_steps+1)

  do i = 0, time_steps-1
    read(iunit, *) dummy, ot(i)
    ot(i) = ot(i) * units%length%factor
  end do

end subroutine read_ot


! ---------------------------------------------------------
subroutine print_omega_file(units, nspin, nfrequencies, warray, tarray)
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use spectrum_m
  use lib_adv_alg_m
  use oscillator_strength_m

  implicit none

  type(unit_system_t), intent(in) :: units
  integer, intent(in)             :: nspin
  integer, intent(in)             :: nfrequencies
  FLOAT, intent(in)               :: warray(nfrequencies), tarray(nfrequencies)

  integer :: iunit, i
  character(len=20) :: header_string

  iunit = io_open('omega', action='write', status='replace')

  write(iunit, '(a15,i2)')      '# nspin        ', nspin
  call kick_write(kick, iunit)
  write(iunit, '(a)') '#%'
  write(iunit, '(a1,a20)', advance = 'no') '#', str_center("omega", 20)
  write(header_string,'(a)') 'F(omega)'
  write(iunit, '(a20)', advance = 'yes') str_center(trim(header_string), 20)
  write(iunit, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units%energy%abbrev) // ']', 20)
  ! Here we should print the units of the transform.
  write(iunit, '(a)', advance = 'yes')

  write(0, *) size(warray), size(tarray)
  do i = 1, nfrequencies
    write(iunit,'(2e20.8)') warray(i) / units%energy%factor, &
                            tarray(i)
  end do

  call io_close(iunit)
end subroutine print_omega_file


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
