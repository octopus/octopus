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

  integer, parameter :: SINE_TRANSFORM   = 1, &
                        COSINE_TRANSFORM = 2

  integer             :: observable(2)
  FLOAT               :: conversion_factor
  type(unit_system_t) :: units
  FLOAT, allocatable  :: ot(:)
  type(kick_t)        :: kick
  integer             :: time_steps
  FLOAT               :: total_time
  integer             :: mode
  FLOAT               :: dt

  contains

    subroutine ft(omega, power)
      FLOAT, intent(in)   :: omega
      FLOAT, intent(out)  :: power
      integer :: i, j
      FLOAT :: x
      power = M_ZERO

      select case(mode)

      case(SINE_TRANSFORM)

        do j = 0, time_steps
          x = sin(omega*j*dt)
          power = power + x*ot(j)
        end do
        power = power*dt / (dt*time_steps)

      case(COSINE_TRANSFORM)

        do j = 0, time_steps
          x = cos(omega*j*dt)
          power = power + x*ot(j)
        end do
        power = power*dt / (dt*time_steps)

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

  integer :: run_mode, order, nfrequencies, nspin, ierr, nresonances, l, m
  FLOAT :: omega, search_interval, final_time
  integer, parameter :: ANALYZE_NTHORDER_SIGNAL           = 1, &
                        GENERATE_NTHORDER_SIGNAL          = 2, &
                        READ_RESONANCES_FROM_FILE         = 3, &
                        GENERATE_OMEGA_FILE               = 4
  character(len=100) :: ffile

  ! Reads the information passed through the command line options (if available).
  call getopt_init(ierr)
  if(ierr.ne.0) then
    message(1) = "Your fortran compiler doesn't support command line arguments,"
    message(2) = "the oct-oscillator-strength command is not available."
    call write_fatal(2)
  end if

  ! Set the default values
  run_mode        = ANALYZE_NTHORDER_SIGNAL
  omega           = - M_ONE
  search_interval = - M_ONE
  order           = 1
  nfrequencies    = 1000
  final_time      = - M_ONE
  nresonances     = 1
  observable(1)   = -1
  observable(2)   = 0
  ffile           = ""

  ! Get the parameters from the command line.
  call getopt_oscillator_strength(run_mode, omega, search_interval,             &
                                  order, nresonances, nfrequencies, final_time, &
                                  observable(1), observable(2), ffile)

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
    call generate_signal(order, observable)
  case(ANALYZE_NTHORDER_SIGNAL)
    call analyze_signal(order, omega, search_interval, final_time, nresonances, nfrequencies)
  case(READ_RESONANCES_FROM_FILE)
    call read_resonances_file(order, ffile, search_interval, final_time, nfrequencies)
  case(GENERATE_OMEGA_FILE)
    call print_omega_file(omega, search_interval, final_time, nfrequencies)
  case default
  end select

  call io_end()
  call datasets_end()
  call global_end()
end program oscillator_strength
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine read_resonances_file(order, ffile, search_interval, final_time, nfrequencies)
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

  integer, intent(inout)       :: order
  character(len=*), intent(in) :: ffile
  FLOAT,   intent(inout)       :: search_interval
  FLOAT,   intent(in)          :: final_time
  integer, intent(in)          :: nfrequencies

  FLOAT :: dummy, leftbound, rightbound, w, power, dw
  integer :: iunit, nresonances, ios, i, j, k, npairs, nspin, order_in_file, nw_substracted, ierr
  logical :: file_exists
  FLOAT, allocatable :: wij(:), omega(:), cij(:)

  if(order.ne.2) then
    write(message(1),'(a)') 'The run mode #3 is only compatible with the analysis of the'
    write(message(2),'(a)') 'second order response.'
    call write_fatal(2)
  end if

  ! First, let us check that the file "ot" exists.
  inquire(file="ot", exist  = file_exists)
  if(.not.file_exists) then
    write(message(1),'(a)') "Could not find 'ot' file."
    call write_fatal(1)
  end if

  ! Now, we should find out which units the file "ot" has.
  call units_from_file(units, "ot", ierr)
  if(ierr.ne.0) then
    write(message(1),'(a)') "Could not retrieve units in the 'ot' file."
    call write_fatal(1)
  end if

  mode = COSINE_TRANSFORM

  iunit = io_open(trim(ffile), action='read', status='old', die=.false.)
  if(iunit.eq.0) then
    write(message(1),'(a)') 'Could not open '//trim(ffile)//' file.'
    call write_fatal(1)
  end if

  call io_skip_header(iunit)
  ! Count the number of resonances
  nresonances = 0
  do
    read(iunit, *, iostat = ios) dummy, dummy
    if(ios.ne.0) exit
    nresonances = nresonances + 1
  end do

  npairs = (nresonances*(nresonances-1))/2

  ALLOCATE(omega(nresonances), nresonances)
  ALLOCATE(cij(nresonances), nresonances)
  ALLOCATE(wij(npairs), npairs)

  call io_skip_header(iunit)
  do i = 1, nresonances
    read(iunit, *) omega(i), cij(i)
  end do

  k = 1
  do i = 1, nresonances
    do j = i + 1, nresonances
      wij(k) = omega(j) - omega(i)
      k = k + 1
    end do
  end do

  if(search_interval > M_ZERO) then
    search_interval = search_interval * units%energy%factor
  else
    search_interval = M_HALF
  end if

  call read_ot(nspin, order_in_file, nw_substracted)

  if(order_in_file.ne.order) then
    write(message(1), '(a)') 'Internal error in analyze_signal'
    call write_fatal(1)
  end if

  if(final_time > M_ZERO) then
    total_time = final_time * units%time%factor
    if(total_time > dt*time_steps) then
      total_time = dt*time_steps
      write(0, '(a)')        '* WARNING: The requested total time to process is larger than the time'
      write(0, '(a)')        '*          available in the input file.'
      write(0, '(a,f8.4,a)') '           The time has been adjusted to ', total_time / units%time%factor, &
                             units%time%abbrev
    end if
    time_steps = int(total_time / dt)
    total_time = time_steps * dt
  else
    total_time = dt*time_steps
  end if

  dw = (rightbound-leftbound) / (nfrequencies - 1)

  ! First, substract zero resonance...
  w = M_ZERO
  call resonance_second_order(w, power, nw_substracted, dw, leftbound, rightbound)
  call modify_ot(time_steps, dt, order, ot, omega, sqrt(abs(power)))
  nw_substracted = nw_substracted + 1

  ! Then, get all the others...
  do k = 1, npairs
    leftbound = wij(k) - search_interval
    rightbound = wij(k) + search_interval
    call find_resonance(wij(k), leftbound, rightbound, nfrequencies)
    call resonance_second_order(wij(k), power, nw_substracted, dw, leftbound, rightbound)
    call modify_ot(time_steps, dt, order, ot, wij(k), sqrt(abs(power)))
  end do
 
  deallocate(wij, cij, omega)
  call io_close(iunit)
end subroutine read_resonances_file
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine analyze_signal(order, omega, search_interval, final_time, nresonances, nfrequencies)
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

  integer, intent(inout) :: order
  FLOAT,   intent(inout) :: omega
  FLOAT,   intent(inout) :: search_interval
  FLOAT,   intent(inout) :: final_time
  integer, intent(inout) :: nresonances
  integer, intent(inout) :: nfrequencies

  FLOAT :: leftbound, rightbound, dw, power
  integer :: nspin, i, ierr, order_in_file, nw_substracted
  logical :: file_exists

  ! First, let us check that the file "ot" exists.
  inquire(file="ot", exist  = file_exists)
  if(.not.file_exists) then
    write(message(1),'(a)') "Could not find 'ot' file."
    call write_fatal(1)
  end if

  ! Now, we should find out which units the file "ot" has.
  call units_from_file(units, "ot", ierr)
  if(ierr.ne.0) then
    write(message(1),'(a)') "Could not retrieve units in the 'ot' file."
    call write_fatal(1)
  end if

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

  leftbound = omega - search_interval
  rightbound = omega + search_interval

  call read_ot(nspin, order_in_file, nw_substracted)

  if(order_in_file.ne.order) then
    write(message(1), '(a)') 'Internal error in analyze_signal'
    call write_fatal(1)
  end if

  if(mod(order, 2).eq.1) then
    mode = SINE_TRANSFORM
  else
    mode = COSINE_TRANSFORM
  end if

  if(final_time > M_ZERO) then
    total_time = final_time * units%time%factor
    if(total_time > dt*time_steps) then
      total_time = dt*time_steps
      write(0, '(a)')        '* WARNING: The requested total time to process is larger than the time'
      write(0, '(a)')        '*          available in the input file.'
      write(0, '(a,f8.4,a)') '           The time has been adjusted to ', total_time / units%time%factor, &
                         units%time%abbrev
    end if
    time_steps = int(total_time / dt)
    total_time = time_steps * dt
  else
    total_time = dt*time_steps
  end if

  dw = (rightbound-leftbound) / (nfrequencies - 1)

  do
    if(nw_substracted >= nresonances) exit

    if(mode == COSINE_TRANSFORM .and. nw_substracted == 0) then
      omega = M_ZERO
    else
      call find_resonance(omega, leftbound, rightbound, nfrequencies)
    end if

    select case(order)
    case(1)
      call resonance_first_order(omega, power, nw_substracted, dw, leftbound, rightbound)
    case(2)
      call resonance_second_order(omega, power, nw_substracted, dw, leftbound, rightbound)
    end select

    call modify_ot(time_steps, dt, order, ot, omega, sqrt(abs(power)))

    nw_substracted = nw_substracted + 1
  end do

  deallocate(ot)
end subroutine analyze_signal
! ---------------------------------------------------------


! ---------------------------------------------------------
! TODO This subroutine should be simplified.
subroutine find_resonance(omega, leftbound, rightbound, nfrequencies)
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use spectrum_m
  use oscillator_strength_m

  implicit none

  FLOAT, intent(inout) :: omega
  FLOAT, intent(in)    :: leftbound, rightbound
  integer, intent(in)  :: nfrequencies

  integer :: i, ierr
  FLOAT :: dw, w, aw, min_aw, min_w, omega_orig
  FLOAT, allocatable :: warray(:), tarray(:)

  ALLOCATE(warray(nfrequencies), nfrequencies)
  ALLOCATE(tarray(nfrequencies), nfrequencies)

  warray = M_ZERO; tarray = M_ZERO
  dw = (rightbound-leftbound) / (nfrequencies - 1)
  do i = 1, nfrequencies
    w = leftbound + (i-1)*dw
    warray(i) = w
    call ft(w, aw)
    tarray(i) = -aw
  end do

  min_w = omega
  min_aw = M_ZERO
  do i = 1, nfrequencies
    w = leftbound + (i-1)*dw
    if(-tarray(i) < min_aw) then
      min_aw = -tarray(i)
      min_w = w
    end if
  end do

  omega_orig = omega
  omega = min_w
  call loct_1dminimize(min_w - 2*dw, min_w + 2*dw, omega, ft, ierr)
  if(ierr.ne.0) then
    write(message(1),'(a)') 'Could not find a maximum.'
    write(message(2),'(a)')
    write(message(3), '(a,f12.8,a,f12.8,a)') '   Search interval = [', leftbound / units%energy%factor, ',', &
                                                          rightbound / units%energy%factor, ']'
    write(message(4), '(a,f12.4,a)')         '   Search discretization = ', dw / units%energy%factor, &
                                ' '//trim(units%energy%abbrev)
    call write_fatal(4)
  end if

  deallocate(warray, tarray)
end subroutine find_resonance
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine resonance_first_order(omega, power, nw_substracted, dw, leftbound, rightbound)
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use spectrum_m
  use oscillator_strength_m

  implicit none

  FLOAT, intent(in)               :: omega
  FLOAT, intent(out)              :: power
  integer, intent(in)             :: nw_substracted
  FLOAT, intent(in)               :: dw, leftbound, rightbound


  call ft(omega, power)
  power = -power
  select case(mode)
  case(SINE_TRANSFORM)
    power = power / (M_ONE - sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))
  case(COSINE_TRANSFORM)
    ! WARNING: Something should go here.
  end select

  write(*, '(a)')                 '******************************************************************'
  write(*, '(a,i3)')              'Resonance #', nw_substracted + 1
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
  write(*, '(a)')                 '******************************************************************'
  write(*, '(a)')


end subroutine resonance_first_order
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine resonance_second_order(omega, power, nw_substracted, dw, leftbound, rightbound)
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use spectrum_m
  use oscillator_strength_m

  implicit none

  FLOAT, intent(in)               :: omega
  FLOAT, intent(out)              :: power
  integer, intent(in)             :: nw_substracted
  FLOAT, intent(in)               :: dw, leftbound, rightbound

  call ft(omega, power)
  power = -power
  select case(mode)
  case(SINE_TRANSFORM)
    power = power / (M_ONE - sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))
  case(COSINE_TRANSFORM)
    ! WARNING: Something should go here.
  end select

  write(*, '(a)')                 '******************************************************************'
  write(*, '(a,i3)')              'Resonance #', nw_substracted + 1
  write(*, '(a,f12.8,a,f12.8,a)') 'omega    = ', omega / units%energy%factor, &
                                  ' '//trim(units%energy%abbrev)//' = ',      &
                                  omega, ' Ha'
  write(*, '(a,f12.8,a,f12.8,a)') 'C(omega) = ', power / units%length%factor**3, &
                                  ' '//trim(units%length%abbrev)//'^3 = ',       &
                                  power, ' b^3'
  write(*, '(a)')
  write(*, '(a,f12.8,a,f12.8,a)') '   Search interval = [', leftbound / units%energy%factor, ',', &
                                                            rightbound / units%energy%factor, ']'
  write(*, '(a)')                 '******************************************************************'
  write(*, '(a)')

end subroutine resonance_second_order
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine generate_signal(order, observable)
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
  integer, intent(in) :: observable(2)

  logical :: file_exists
  integer :: i, j, nspin, time_steps, lmax, nfiles, k, isp, add_lm, l, m, max_add_lm
  integer, allocatable :: iunit(:)
  FLOAT :: dt, lambda, det, dump, o0, conversion_factor, lequalonefactor
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

  if(kick%l > 0) then
    max_add_lm = (kick%l+1)**2-1
    ALLOCATE(multipole(max_add_lm, 0:time_steps, nspin), 2*(time_steps+1)*nspin)
    conversion_factor = units%length%factor ** kick%l
  else
    max_add_lm = 3
    ALLOCATE(multipole(3, 0:time_steps, nspin), 2*(time_steps+1)*nspin)
    conversion_factor = units%length%factor
  end if
  ALLOCATE(ot(0:time_steps), time_steps+1)
  multipole = M_ZERO
  ot = M_ZERO

  do j = 1, nfiles
    call io_skip_header(iunit(j))

    do i = 0, time_steps
      select case(nspin)
      case(1)
        read(iunit(j), *) k, dump, dump, (multipole(add_lm, i, 1), add_lm = 1, max_add_lm)
      case(2)
        read(iunit(j), *) k, dump, dump, (multipole(add_lm, i, 1), add_lm = 1, max_add_lm), &
                          dump, (multipole(add_lm, i, 2), add_lm = 1, max_add_lm)
      case(4)
        read(iunit(j), *) k, dump, dump, (multipole(add_lm, i, 1), add_lm = 1, max_add_lm), &
                          dump, (multipole(add_lm, i, 2), add_lm = 1, max_add_lm), &
                          dump, (multipole(add_lm, i, 3), add_lm = 1, max_add_lm), &
                          dump, (multipole(add_lm, i, 4), add_lm = 1, max_add_lm)
      end select
      multipole(1:3, i, :) = multipole(1:3, i, :) * conversion_factor

      select case(observable(1))
      case(-1)
        dump = M_ZERO
        if(kick%l > 0) then
          add_lm = 1
          lcycle1: do l = 1, kick%l
            mcycle1: do m = -l, l
              if(kick%l.eq.l .and. kick%m.eq.m) exit lcycle1
              add_lm = add_lm + 1
            end do mcycle1
          end do lcycle1
          ! The multipoles file treats differently the l=1 case.
          if(kick%l == 1) then
            lequalonefactor = sqrt(CNST(3.0)/(M_FOUR*M_Pi))
            select case(kick%m)
            case(-1); add_lm = 1
            case(0);  add_lm = 3
            case(1);  add_lm = 2
            end select
          else
            lequalonefactor = M_ONE
          end if
          do isp = 1, nspin
            dump = dump + lequalonefactor * multipole(add_lm, i, isp)
          end do
        else
          do isp = 1, nspin
            dump = dump + sum(multipole(1:3, i, isp) * kick%pol(1:3,kick%pol_dir))
          end do
        end if
      case(0)
        dump = M_ZERO
        do isp = 1, nspin
          dump = dump + multipole(observable(2), i, isp) 
        end do
      case default
        add_lm = 1
        lcycle2: do l = 1, observable(1)
          mcycle2: do m = -l, l
            add_lm = add_lm + 1
            if(observable(1).eq.l .and. observable(2).eq.2) exit lcycle2
          end do mcycle2
        end do lcycle2
        dump = M_ZERO
        do isp = 1, nspin
          dump = dump + multipole(add_lm, i, isp)
        end do
      end select

      if(i == 0) o0 = dump

      ot(i) = ot(i) + mu(j)*(dump - o0)

    end do

  end do

  ot = ot / lambda**order

  call write_ot(nspin, time_steps, dt, kick, units, order, ot(0:time_steps), observable)

  ! Close files and exit.
  do j = 1, nfiles
    call io_close(iunit(j))
  end do
  deallocate(iunit, q, mu, qq, c, ot, multipole)
end subroutine generate_signal
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine modify_ot(time_steps, dt, order, ot, omega, cij)
  use global_m
  use io_m
  use units_m
  use spectrum_m
  use string_m

  implicit none

  integer,             intent(in)    :: time_steps
  FLOAT,               intent(in)    :: dt
  integer,             intent(in)    :: order
  FLOAT,               intent(inout) :: ot(0:time_steps)
  FLOAT,               intent(in)    :: omega
  FLOAT,               intent(in)    :: cij

  integer :: i

  select case(mod(order, 2))
  case(1)
    do i = 0, time_steps
      ot(i) = ot(i) + M_TWO * abs(cij)**2 * sin(omega*i*dt)
    end do
  case(0)
    do i = 0, time_steps
      ot(i) = ot(i) - abs(cij)**2 * cos(omega*i*dt)
    end do
  end select

end subroutine modify_ot
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine write_ot(nspin, time_steps, dt, kick, units, order, ot, observable)
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
  integer,             intent(in) :: order
  FLOAT,               intent(in) :: ot(0:time_steps)
  integer,             intent(in) :: observable(2)

  integer :: iunit, i
  character(len=20) :: header_string
  FLOAT :: conversion_factor

  iunit = io_open('ot', action='write', status='replace')

  write(iunit, '(a15,i2)')      '# nspin        ', nspin
  write(iunit, '(a15,i2)')      '# Order        ', order
  write(iunit, '(a28,i2)')      '# Frequencies substracted = ', 0
  select case(observable(1))
  case(-1)
    write(iunit,'(a)') '# Observable operator = kick operator'
    if(kick%l > 0 ) then
      conversion_factor = units%length%factor ** kick%l
    else
      conversion_factor = units%length%factor
    end if
  case(0)
    select case(observable(2))
    case(1); write(iunit,'(a)') '# O = x'
    case(2); write(iunit,'(a)') '# O = y'
    case(3); write(iunit,'(a)') '# O = z'
    end select
    conversion_factor = units%length%factor
  case default
    conversion_factor = units%length%factor ** observable(1)
    write(iunit, '(a12,i1,a1,i2,a1)') '# (l, m) = (', observable(1),',',observable(2),')'
  end select
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
    write(iunit, '(2e20.8)') i*dt / units%time%factor, ot(i) / conversion_factor
  end do

  call io_close(iunit)
end subroutine write_ot


! ---------------------------------------------------------
subroutine read_ot(nspin, order, nw_substracted)
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
  integer, intent(out) :: order
  integer, intent(out) :: nw_substracted

  integer :: iunit, i, ierr
  character(len=100) :: line
  character(len=12)  :: dummychar
  FLOAT :: dummy, t1, t2

  iunit = io_open('ot', action='read', status='old')
  if(iunit .eq. 0) then
    write(message(1),'(a)') 'A file called ot should be present and was not found.'
    call write_fatal(1)
  end if

  read(iunit, '(15x,i2)') nspin
  read(iunit, '(15x,i2)') order
  read(iunit, '(28x,i2)') nw_substracted
  read(iunit, '(a)')      line

  i = index(line, 'Observable')
  if(index(line, 'Observable').ne.0) then
    observable(1) = -1
  elseif(index(line, '# O =').ne.0) then
    observable(1) = 0
    if(index(line,'x').ne.0) then
      observable(2) = 1
    elseif(index(line,'y').ne.0) then
      observable(2) = 2
    elseif(index(line,'z').ne.0) then
      observable(2) = 3
    end if
  elseif(index(line, '# (l, m) = ').ne.0) then
    read(line,'(a12,i1,a1,i2,a1)') dummychar(1:12), observable(1), dummychar(1:1), observable(2), dummychar(1:1)
  else
    write(message(1),'(a)') 'Problem reading ot file: could not figure out the shape'
    write(message(2),'(a)') 'of the observation operator.'
    call write_fatal(2)
  end if

  call kick_read(kick, iunit)
  read(iunit, '(a)')  line
  read(iunit, '(a)')  line
  call io_skip_header(iunit)

  ! Figure out about the units of the file
  call units_from_file(units, "ot", ierr)
  if(ierr.ne.0) then
    write(message(1), '(a)') 'Could not figure out the units in file "ot".'
    call write_fatal(1)
  end if

  select case(observable(1))
  case(-1)
    if(kick%l > 0) then
      conversion_factor = units%length%factor**kick%l
    else
      conversion_factor = units%length%factor
    end if
  case(0)
    conversion_factor = units%length%factor
  case default
    conversion_factor = units%length%factor**observable(1)
  end select


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
subroutine print_omega_file(omega, search_interval, final_time, nfrequencies)
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

  FLOAT,   intent(inout) :: omega
  FLOAT,   intent(inout) :: search_interval
  FLOAT,   intent(inout) :: final_time
  integer, intent(inout) :: nfrequencies

  integer :: iunit, i, ierr, nspin, order, nw_substracted
  logical :: file_exists
  character(len=20) :: header_string
  FLOAT, allocatable :: warray(:), tarray(:)
  FLOAT :: leftbound, rightbound, dw, power, w, aw

  ! First, let us check that the file "ot" exists.
  inquire(file="ot", exist  = file_exists)
  if(.not.file_exists) then
    write(message(1),'(a)') "Could not find 'ot' file."
    call write_fatal(1)
  end if

  ! Now, we should find out which units the file "ot" has.
  call units_from_file(units, "ot", ierr)
  if(ierr.ne.0) then
    write(message(1),'(a)') "Could not retrieve units in the 'ot' file."
    call write_fatal(1)
  end if

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

  leftbound = omega - search_interval
  rightbound = omega + search_interval

  ALLOCATE(warray(nfrequencies), nfrequencies)
  ALLOCATE(tarray(nfrequencies), nfrequencies)

  call read_ot(nspin, order, nw_substracted)

  if(mod(order, 2).eq.1) then
    mode = SINE_TRANSFORM
  else
    mode = COSINE_TRANSFORM
  end if

  if(final_time > M_ZERO) then
    total_time = final_time * units%time%factor
    if(total_time > dt*time_steps) then
      total_time = dt*time_steps
      write(0, '(a)')        '* WARNING: The requested total time to process is larger than the time'
      write(0, '(a)')        '*          available in the input file.'
      write(0, '(a,f8.4,a)') '           The time has been adjusted to ', total_time / units%time%factor, &
                         units%time%abbrev
    end if
    time_steps = int(total_time / dt)
    total_time = time_steps * dt
  else
    total_time = dt*time_steps
  end if

  warray = M_ZERO; tarray = M_ZERO
  dw = (rightbound-leftbound) / (nfrequencies - 1)
  do i = 1, nfrequencies
    w = leftbound + (i-1)*dw
    warray(i) = w
    call ft(w, aw)
    tarray(i) = -aw
  end do

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

  do i = 1, nfrequencies
    write(iunit,'(2e20.8)') warray(i) / units%energy%factor, &
                            tarray(i)
  end do

  deallocate(warray, tarray)
end subroutine print_omega_file
! ---------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
