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
  use datasets_m
  use global_m
  use io_m
  use kick_m
  use loct_m
  use loct_math_m
  use parser_m
  use messages_m
  use profiling_m
  use spectrum_m
  use string_m
  use unit_m
  use unit_system_m

  implicit none

  integer, parameter :: SINE_TRANSFORM   = 1, &
                        COSINE_TRANSFORM = 2

  type local_operator_t
    integer :: n_multipoles
    integer, pointer :: l(:), m(:)
    FLOAT,   pointer :: weight(:)
  end type local_operator_t

  integer             :: observable(2)
  type(unit_system_t) :: units
  FLOAT, allocatable  :: ot(:)
  type(kick_t)        :: kick
  integer             :: time_steps
  FLOAT               :: total_time
  integer             :: mode
  FLOAT               :: dt


  contains


    ! ---------------------------------------------------------
    subroutine ft(omega, power)
      FLOAT, intent(in)   :: omega
      FLOAT, intent(out)  :: power
      integer :: j
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


    ! ---------------------------------------------------------
    subroutine ft2(omega, power)
      FLOAT, intent(in)   :: omega
      FLOAT, intent(out)  :: power
      integer :: j
      FLOAT :: x
      power = M_ZERO

      select case(mode)

      case(SINE_TRANSFORM)

        do j = 0, time_steps
          x = sin(omega*j*dt)
          power = power + x*ot(j)
        end do
        ! The function should be *minus* the sine Fourier transform, since this
        ! is the function to be minimized.
        power = - (power*dt / (dt*time_steps))**2

      case(COSINE_TRANSFORM)

        do j = 0, time_steps
          x = cos(omega*j*dt)
          power = power + x*ot(j)
        end do
        power = - (power*dt / (dt*time_steps))**2

      end select

    end subroutine ft2


    ! ---------------------------------------------------------
    subroutine local_operator_copy(o, i)
      type(local_operator_t), intent(inout) :: o
      type(local_operator_t), intent(inout) :: i
      integer :: j

      o%n_multipoles = i%n_multipoles
      SAFE_ALLOCATE(     o%l(1:o%n_multipoles))
      SAFE_ALLOCATE(     o%m(1:o%n_multipoles))
      SAFE_ALLOCATE(o%weight(1:o%n_multipoles))

      do j = 1, o%n_multipoles
        o%l(j) = i%l(j)
        o%m(j) = i%m(j)
        o%weight(j) = i%weight(j)
      end do

    end subroutine local_operator_copy

end module oscillator_strength_m



! ---------------------------------------------------------
! ---------------------------------------------------------
! ---------------------------------------------------------
program oscillator_strength
  use global_m
  use messages_m
  use datasets_m
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use loct_m
  use profiling_m
  use command_line_m
  use oscillator_strength_m

  implicit none

  integer :: run_mode, order, nfrequencies, ierr, nresonances
  FLOAT :: omega, search_interval, final_time, damping
  integer, parameter :: ANALYZE_NTHORDER_SIGNAL           = 1, &
                        GENERATE_NTHORDER_SIGNAL          = 2, &
                        READ_RESONANCES_FROM_FILE         = 3, &
                        GENERATE_OMEGA_FILE               = 4
  character(len=100) :: ffile

  ! Reads the information passed through the command line options (if available).
  call getopt_init(ierr)
  if(ierr.ne.0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-oscillator-strength command is not available."
    call messages_fatal(2)
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
  damping         = CNST(0.1)/CNST(27.2114) ! This is the usual damping factor used in the literature.

  ! Get the parameters from the command line.
  call getopt_oscillator_strength(run_mode, omega, search_interval,             &
                                  order, nresonances, nfrequencies, final_time, &
                                  observable(1), observable(2), damping, ffile)
  call getopt_end()

  ! Initialize stuff
  call global_init()
  in_debug_mode = .false.
  call io_init(defaults = .true.)
  call datasets_init(1)
  if(in_debug_mode) then
     call io_mkdir('debug')
  end if

  select case(run_mode)
  case(GENERATE_NTHORDER_SIGNAL)
    call generate_signal(order, observable)
  case(ANALYZE_NTHORDER_SIGNAL)
    call analyze_signal(order, omega, search_interval, final_time, nresonances, nfrequencies, damping)
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
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use lalg_adv_m
  use oscillator_strength_m
  use profiling_m

  implicit none

  integer, intent(inout)       :: order
  character(len=*), intent(in) :: ffile
  FLOAT,   intent(inout)       :: search_interval
  FLOAT,   intent(in)          :: final_time
  integer, intent(in)          :: nfrequencies

  FLOAT :: dummy, leftbound, rightbound, w, power, dw
  integer :: iunit, nresonances, ios, i, j, k, npairs, nspin, order_in_file, nw_subtracted, ierr
  logical :: file_exists
  FLOAT, allocatable :: wij(:), omega(:), c0i(:)

  if(order.ne.2) then
    write(message(1),'(a)') 'The run mode #3 is only compatible with the analysis of the'
    write(message(2),'(a)') 'second-order response.'
    call messages_fatal(2)
  end if

  ! First, let us check that the file "ot" exists.
  inquire(file="ot", exist  = file_exists)
  if(.not.file_exists) then
    write(message(1),'(a)') "Could not find 'ot' file."
    call messages_fatal(1)
  end if

  ! Now, we should find out which units the file "ot" has.
  call unit_system_from_file(units, "ot", ierr)
  if(ierr.ne.0) then
    write(message(1),'(a)') "Could not retrieve units in the 'ot' file."
    call messages_fatal(1)
  end if

  mode = COSINE_TRANSFORM

  iunit = io_open(trim(ffile), action='read', status='old', die=.false.)
  if(iunit.eq.0) then
    write(message(1),'(a)') 'Could not open '//trim(ffile)//' file.'
    call messages_fatal(1)
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

  SAFE_ALLOCATE(omega(1:nresonances))
  SAFE_ALLOCATE(  c0i(1:nresonances))
  SAFE_ALLOCATE(wij(1:npairs))

  call io_skip_header(iunit)
  do i = 1, nresonances
    read(iunit, *) omega(i), c0i(i)
  end do

  k = 1
  do i = 1, nresonances
    do j = i + 1, nresonances
      wij(k) = omega(j) - omega(i)
      k = k + 1
    end do
  end do

  if(search_interval > M_ZERO) then
    search_interval = units_to_atomic(units%energy, search_interval)
  else
    search_interval = M_HALF
  end if

  call read_ot(nspin, order_in_file, nw_subtracted)

  if(order_in_file.ne.order) then
    write(message(1), '(a)') 'Error: The ot file should contain the second-order response'
    write(message(2), '(a)') '       in this run mode.'
    call messages_fatal(2)
  end if

  if(final_time > M_ZERO) then
    total_time = units_to_atomic(units%time, final_time)
    if(total_time > dt*time_steps) then
      total_time = dt*time_steps
      write(0, '(a)')        '* WARNING: The requested total time to process is larger than the time'
      write(0, '(a)')        '*          available in the input file.'
      write(0, '(a,f8.4,a)') '           The time has been adjusted to ', &
        units_from_atomic(units%time, total_time), units_abbrev(units%time)
    end if
    time_steps = int(total_time / dt)
    total_time = time_steps * dt
  else
    total_time = dt*time_steps
  end if

  dw = (rightbound-leftbound) / (nfrequencies - 1)

  ! First, subtract zero resonance...
  w = M_ZERO
  call resonance_second_order(w, power, nw_subtracted, leftbound, rightbound, M_ZERO, M_ZERO)
  call modify_ot(time_steps, dt, order, ot, w, power)
  nw_subtracted = nw_subtracted + 1

  ! Then, get all the others...
  k = 1
  do i = 1, nresonances
    do j = i + 1, nresonances
      leftbound = wij(k) - search_interval
      rightbound = wij(k) + search_interval
      call find_resonance(wij(k), leftbound, rightbound, nfrequencies)
      call resonance_second_order(wij(k), power, nw_subtracted, leftbound, rightbound, c0i(i), c0i(j))
      call modify_ot(time_steps, dt, order, ot, wij(k), power)
      nw_subtracted = nw_subtracted + 1
      k = k + 1
    end do
  end do
 
  SAFE_DEALLOCATE_A(wij)
  SAFE_DEALLOCATE_A(c0i)
  SAFE_DEALLOCATE_A(omega)
  call io_close(iunit)
end subroutine read_resonances_file
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine analyze_signal(order, omega, search_interval, final_time, nresonances, nfrequencies, damping)
  use global_m
  use messages_m
  use datasets_m
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use lalg_adv_m
  use oscillator_strength_m
  use profiling_m

  implicit none

  integer, intent(inout) :: order
  FLOAT,   intent(inout) :: omega
  FLOAT,   intent(inout) :: search_interval
  FLOAT,   intent(inout) :: final_time
  integer, intent(inout) :: nresonances
  integer, intent(inout) :: nfrequencies
  FLOAT,   intent(in)    :: damping

  FLOAT :: leftbound, rightbound, dw, power
  FLOAT, allocatable :: w(:), c0I2(:)
  integer :: nspin, i, ierr, order_in_file, nw_subtracted
  logical :: file_exists

  ! First, let us check that the file "ot" exists.
  inquire(file="ot", exist  = file_exists)
  if(.not.file_exists) then
    write(message(1),'(a)') "Could not find 'ot' file."
    call messages_fatal(1)
  end if

  ! Now, we should find out which units the file "ot" has.
  call unit_system_from_file(units, "ot", ierr)
  if(ierr.ne.0) then
    write(message(1),'(a)') "Could not retrieve units in the 'ot' file."
    call messages_fatal(1)
  end if

  if(omega > M_ZERO) then
    omega = units_to_atomic(units%energy, omega)
  else
    omega = M_HALF
  end if

  if(search_interval > M_ZERO) then
    search_interval = units_to_atomic(units%energy, search_interval)
  else
    search_interval = M_HALF
  end if

  leftbound = omega - search_interval
  rightbound = omega + search_interval

  call read_ot(nspin, order_in_file, nw_subtracted)

  if(order_in_file.ne.order) then
    write(message(1), '(a)') 'Internal error in analyze_signal'
    call messages_fatal(1)
  end if

  if(mod(order, 2).eq.1) then
    mode = SINE_TRANSFORM
  else
    mode = COSINE_TRANSFORM
  end if

  if(final_time > M_ZERO) then
    total_time = units_to_atomic(units%time, final_time)
    if(total_time > dt*time_steps) then
      total_time = dt*time_steps
      write(0, '(a)')        '* WARNING: The requested total time to process is larger than the time'
      write(0, '(a)')        '*          available in the input file.'
      write(0, '(a,f8.4,a)') '           The time has been adjusted to ', &
        units_from_atomic(units%time, total_time), units_abbrev(units%time)
    end if
    time_steps = int(total_time / dt)
    total_time = time_steps * dt
  else
    total_time = dt*time_steps
  end if

  dw = (rightbound-leftbound) / (nfrequencies - 1)

  SAFE_ALLOCATE(   w(1:nresonances))
  SAFE_ALLOCATE(c0I2(1:nresonances))
  w = M_ZERO
  c0I2 = M_ZERO

  i = 1
  do
    if(nw_subtracted >= nresonances) exit

    if(mode == COSINE_TRANSFORM .and. nw_subtracted == 0) then
      omega = M_ZERO
    else
      call find_resonance(omega, leftbound, rightbound, nfrequencies)
    end if

    select case(order)
    case(1)
      call resonance_first_order(omega, power, nw_subtracted, dw, leftbound, rightbound)
    case(2)
      call resonance_second_order(omega, power, nw_subtracted, leftbound, rightbound, M_ZERO, M_ZERO)
    end select

    w(i) = omega
    c0I2(i) = power

    call modify_ot(time_steps, dt, order, ot, omega, power)

    nw_subtracted = nw_subtracted + 1
    i = i + 1
  end do

  select case(order)
    case(1)
      call write_polarizability(nfrequencies, nresonances, dw, w, c0I2, damping)
  end select

  SAFE_DEALLOCATE_A(ot)
  SAFE_DEALLOCATE_A(w)
  SAFE_DEALLOCATE_A(c0I2)
end subroutine analyze_signal
! ---------------------------------------------------------


! ---------------------------------------------------------
! Implements the SOS formula of the polarizability, and writes
! down to the "polarizability" file the real and imaginary part
! of the dynamical polarizability.
! ---------------------------------------------------------
subroutine write_polarizability(nfrequencies, nresonances, dw, w, c0I2, gamma)
  use global_m
  use messages_m
  use datasets_m
  use io_m
  use profiling_m

  implicit none

  integer, intent(in) :: nfrequencies, nresonances
  FLOAT, intent(in)   :: dw
  FLOAT, intent(in)   :: w(nresonances), c0I2(nresonances)
  FLOAT, intent(in)   :: gamma

  integer :: iunit, i, j
  FLOAT :: e
  CMPLX :: pol

  iunit = io_open('polarizability', status='replace', action = 'write', die=.false.)
  write(iunit, '(a)') '# Polarizability file. Generated using the SOS formula with the following data:'
  write(iunit, '(a)') '#'

  do i = 1, nresonances
    write(iunit, '(a1,3e20.12)') '#', w(i), sqrt(abs(c0I2(i))), c0I2(i)
  end do

  write(iunit, '(a10,f12.6)') '# Gamma = ', gamma
  write(iunit, '(a)')         '#'

  do i = 1, nfrequencies
    e = (i-1)*dw
    pol = M_z0
    do j = 1, nresonances
      pol = pol + c0I2(j) * ( M_ONE/(w(j)- e - M_zI*gamma) + M_ONE/(w(j) + e + M_zI*gamma )  )
    end do
    write(iunit, '(3e20.12)') e, pol
  end do 

  call io_close(iunit)
end subroutine write_polarizability
! ---------------------------------------------------------


! ---------------------------------------------------------
! \todo This subroutine should be simplified.
subroutine find_resonance(omega, leftbound, rightbound, nfrequencies)
  use global_m
  use messages_m
  use datasets_m
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use oscillator_strength_m
  use profiling_m

  implicit none

  FLOAT, intent(inout) :: omega
  FLOAT, intent(in)    :: leftbound, rightbound
  integer, intent(in)  :: nfrequencies

  integer :: i, ierr
  FLOAT :: dw, w, aw, min_aw, min_w, omega_orig
  FLOAT, allocatable :: warray(:), tarray(:)

  SAFE_ALLOCATE(warray(1:nfrequencies))
  SAFE_ALLOCATE(tarray(1:nfrequencies))

  warray = M_ZERO; tarray = M_ZERO
  dw = (rightbound-leftbound) / (nfrequencies - 1)
  do i = 1, nfrequencies
    w = leftbound + (i-1)*dw
    warray(i) = w
    call ft2(w, aw)
    tarray(i) = aw
  end do

  min_w = omega
  min_aw = M_ZERO
  do i = 1, nfrequencies
    w = leftbound + (i-1)*dw
    if(tarray(i) < min_aw) then
      min_aw = tarray(i)
      min_w = w
    end if
  end do

  omega_orig = omega
  omega = min_w
#ifndef SINGLE_PRECISION
  call loct_1dminimize(min_w - 2*dw, min_w + 2*dw, omega, ft2, ierr)
#else
  stop "FIXME: cannot work in single-precision."
#endif
  if(ierr.ne.0) then
    write(message(1),'(a)') 'Could not find a maximum.'
    write(message(2),'(a)')
    write(message(3), '(a,f12.8,a,f12.8,a)') '   Search interval = [', &
      units_from_atomic(units%energy, leftbound), ',', units_from_atomic(units%energy, rightbound), ']'
    write(message(4), '(a,f12.4,a)')         '   Search discretization = ', &
      units_from_atomic(units%energy, dw), ' '//trim(units_abbrev(units%energy))
    call messages_fatal(4)
  end if

  SAFE_DEALLOCATE_A(warray)
  SAFE_DEALLOCATE_A(tarray)
end subroutine find_resonance
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine resonance_first_order(omega, power, nw_subtracted, dw, leftbound, rightbound)
  use global_m
  use messages_m
  use datasets_m
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use oscillator_strength_m
  use profiling_m

  implicit none

  FLOAT, intent(in)               :: omega
  FLOAT, intent(out)              :: power
  integer, intent(in)             :: nw_subtracted
  FLOAT, intent(in)               :: dw, leftbound, rightbound

  call ft(omega, power)

  select case(mode)
  case(SINE_TRANSFORM)
    power = power / (M_ONE - sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))
  case(COSINE_TRANSFORM)
    ! WARNING: Something should go here.
  end select

  write(*, '(a)')                 '******************************************************************'
  write(*, '(a,i3)')              'Resonance #', nw_subtracted + 1
  write(*, '(a,f12.8,a,f12.8,a)') 'omega    = ', units_from_atomic(units_out%energy, omega), &
                                  ' '//trim(units_abbrev(units_out%energy))//' = ',      &
                                  omega, ' Ha'
  write(*, '(a,f12.8,a,f12.8,a)') 'C(omega) = ', units_from_atomic(units_out%length**2, power), &
                                  ' '//trim(units_abbrev(units_out%length**2))//' =',        &
                                  power, ' b^2'
  write(*, '(a,f12.8,a,f12.8,a)') '<0|P|I>  = ', units_from_atomic(units_out%length, sqrt(abs(power))), &
                                  ' '//trim(units_abbrev(units_out%length))//' = ',                 &
                                  sqrt(abs(power)),' b'
  write(*, '(a,f12.8)')           'f[O->I]  = ', M_TWO*omega*power
  write(*, '(a)')
  write(*, '(a,f12.8,a,f12.8,a)') '   Search interval = [', units_from_atomic(units_out%energy, leftbound), ',', &
                                                            units_from_atomic(units_out%energy, rightbound), ']'
  write(*, '(a,f12.4,a)')         '   Search discretization = ', units_from_atomic(units_out%energy, dw), &
                                  ' '//trim(units_abbrev(units_out%energy))
  write(*, '(a)')                 '******************************************************************'
  write(*, '(a)')


end subroutine resonance_first_order
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine resonance_second_order(omega, power, nw_subtracted, leftbound, rightbound, c01, c02)
  use global_m
  use messages_m
  use datasets_m
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use oscillator_strength_m
  use profiling_m

  implicit none

  FLOAT, intent(in)               :: omega
  FLOAT, intent(out)              :: power
  integer, intent(in)             :: nw_subtracted
  FLOAT, intent(in)               :: leftbound, rightbound
  FLOAT, intent(in)               :: c01, c02

  call ft(omega, power)
  select case(mode)
  case(SINE_TRANSFORM)
    power = power / (M_ONE - sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))
  case(COSINE_TRANSFORM)
    ! WARNING: there is some difference between the omega=0 case and the rest.
    if(omega.ne.M_ZERO) then
      power = power / (M_ONE + sin(M_TWO*omega*total_time)/(M_TWO*omega*total_time))
    else
      power = power / M_TWO
    end if
  end select

  write(*, '(a)')                 '******************************************************************'
  write(*, '(a,i3)')              'Resonance #', nw_subtracted + 1
  write(*, '(a,f12.8,a,f12.8,a)') 'omega    = ', units_from_atomic(units_out%energy, omega), &
                                  ' '//trim(units_abbrev(units_out%energy))//' = ',      &
                                  omega, ' Ha'
  write(*, '(a,f12.8,a,f12.8,a)') 'C(omega) = ', units_from_atomic(units_out%length**3, power), &
                                  ' '//trim(units_abbrev(units_out%length**3))//' = ',       &
                                  power, ' b^3'

  if(c01*c02 .ne. M_ZERO) then
    write(*, '(a,f12.8)')         '    C(omega)/(C0i*C0j) = ', power / (c01 * c02)
  end if

  write(*, '(a)')
  write(*, '(a,f12.8,a,f12.8,a)') '   Search interval = [', units_from_atomic(units_out%energy, leftbound), ',', &
                                                            units_from_atomic(units_out%energy, rightbound), ']'
  write(*, '(a)')                 '******************************************************************'
  write(*, '(a)')

end subroutine resonance_second_order
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine generate_signal(order, observable)
  use global_m
  use messages_m
  use datasets_m
  use kick_m
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use lalg_adv_m

  use oscillator_strength_m, only : local_operator_t, &
                                    local_operator_copy
  use profiling_m

  implicit none

  integer, intent(in) :: order
  integer, intent(in) :: observable(2)

  logical :: file_exists
  integer :: i, j, nspin, time_steps, lmax, nfiles, k, add_lm, l, m, max_add_lm
  integer, allocatable :: iunit(:)
  FLOAT :: dt, lambda, det, dump, o0
  type(unit_t) :: mp_unit
  FLOAT, allocatable :: q(:), mu(:), qq(:, :), c(:)
  character(len=20) :: filename
  type(kick_t) :: kick
  type(unit_system_t) :: units
  FLOAT, allocatable :: multipole(:, :, :), ot(:), dipole(:, :)
  type(local_operator_t) :: kick_operator
  type(local_operator_t) :: obs

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
    call messages_fatal(1)
  endif
  if(order > nfiles) then
    write(message(1),'(a)') 'The order that you ask for is higher than the number'
    write(message(2),'(a)') 'of multipoles.x file that you supply.'
    call messages_fatal(2)
  end if

  ! Open the files.
  SAFE_ALLOCATE(iunit(1:nfiles))
  do j = 1, nfiles
    write(filename,'(a11,i1)') 'multipoles.', j
    iunit(j) = io_open(trim(filename), action='read', status='old', die=.false.)
  end do

  SAFE_ALLOCATE( q(1:nfiles))
  SAFE_ALLOCATE(mu(1:nfiles))
  SAFE_ALLOCATE(qq(1:nfiles, 1:nfiles))
  SAFE_ALLOCATE( c(1:nfiles))

  c        = M_ZERO
  c(order) = M_ONE

  ! Get the basic info from the first file
  call spectrum_mult_info(iunit(1), nspin, kick, time_steps, dt, units, lmax=lmax)

  ! Sets the kick operator...
  if(kick%n_multipoles > 0) then
    kick_operator%n_multipoles = kick%n_multipoles
    SAFE_ALLOCATE(     kick_operator%l(1:kick_operator%n_multipoles))
    SAFE_ALLOCATE(     kick_operator%m(1:kick_operator%n_multipoles))
    SAFE_ALLOCATE(kick_operator%weight(1:kick_operator%n_multipoles))
    do i = 1, kick_operator%n_multipoles
      kick_operator%l(i) = kick%l(i)
      kick_operator%m(i) = kick%m(i)
      kick_operator%weight(i) = kick%weight(i)
    end do
  else
    kick_operator%n_multipoles = 3
    SAFE_ALLOCATE(     kick_operator%l(1:kick_operator%n_multipoles))
    SAFE_ALLOCATE(     kick_operator%m(1:kick_operator%n_multipoles))
    SAFE_ALLOCATE(kick_operator%weight(1:kick_operator%n_multipoles))
    kick_operator%l(1:3) = 1
    kick_operator%m(1) = -1
    kick_operator%m(2) =  0
    kick_operator%m(3) =  1
    ! WARNING: not sure if m = -1 => y, and m = 1 => x. What is the convention?
    kick_operator%weight(1) = -sqrt((M_FOUR*M_PI)/M_THREE) * kick%pol(2, kick%pol_dir)
    kick_operator%weight(2) =  sqrt((M_FOUR*M_PI)/M_THREE) * kick%pol(3, kick%pol_dir)
    kick_operator%weight(3) = -sqrt((M_FOUR*M_PI)/M_THREE) * kick%pol(1, kick%pol_dir)
  end if

  ! Sets the observation operator
  select case(observable(1))
    case(-1)
      ! This means that the "observation operator" should be equal 
      ! to the "perturbation operator", i.e., the kick.
      call local_operator_copy(obs, kick_operator)
    case(0)
      ! This means that the observable is the dipole operator; observable(2) determines
      ! if it is x, y or z.
      obs%n_multipoles = 1
      SAFE_ALLOCATE(obs%l(1:1))
      SAFE_ALLOCATE(obs%m(1:1))
      SAFE_ALLOCATE(obs%weight(1:1))
      obs%l(1) = 1
      select case(observable(2))
        case(1)
          obs%m(1) = -1
          obs%weight(1) = -sqrt((M_FOUR*M_PI)/M_THREE)
        case(2)
          obs%m(1) =  1
          obs%weight(1) = sqrt((M_FOUR*M_PI)/M_THREE)
        case(3)
          obs%m(1) =  0
          obs%weight(1) = -sqrt((M_FOUR*M_PI)/M_THREE)
      end select
    case default
      ! This means that the observation operator is (l,m) = (observable(1), observable(2))
      obs%n_multipoles = 1
      SAFE_ALLOCATE(obs%l(1:1))
      SAFE_ALLOCATE(obs%m(1:1))
      SAFE_ALLOCATE(obs%weight(1:1))
      obs%weight(1) = M_ONE
      obs%l(1) = observable(1)
      obs%m(1) = observable(2)
  end select

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

  if(kick%n_multipoles > 0) then
    lmax = maxval(kick%l(1:obs%n_multipoles))
    max_add_lm = (lmax+1)**2-1
    SAFE_ALLOCATE(multipole(1:max_add_lm, 0:time_steps, 1:nspin))
    ! The units have nothing to do with the perturbing kick??
    mp_unit = units%length**kick%l(1)
  else
    max_add_lm = 3
    SAFE_ALLOCATE(multipole(1:3, 0:time_steps, 1:nspin))
    mp_unit = units%length
  end if
  SAFE_ALLOCATE(ot(0:time_steps))
  multipole = M_ZERO
  ot = M_ZERO

  SAFE_ALLOCATE(dipole(1:3, 1:nspin))

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
      multipole(1:max_add_lm, i, :) = units_to_atomic(mp_unit, multipole(1:max_add_lm, i, :))

      ! The dipole is treated differently in the multipoles file: first of all, 
      ! the program should have printed *minus* the dipole operator.
      dipole(1:3, 1:nspin) = - multipole(1:3, i, 1:nspin)
      ! And then it contains the "cartesian" dipole, opposed to the spherical dipole:
      multipole(1, i, 1:nspin) = -sqrt(M_THREE/(M_FOUR*M_PI)) * dipole(2, 1:nspin)
      multipole(2, i, 1:nspin) =  sqrt(M_THREE/(M_FOUR*M_PI)) * dipole(3, 1:nspin)
      multipole(3, i, 1:nspin) = -sqrt(M_THREE/(M_FOUR*M_PI)) * dipole(1, 1:nspin)

      dump = M_ZERO
      do k = 1, obs%n_multipoles
        add_lm = 1; l = 1
        lcycle: do
          do m = -l, l
            if(l == obs%l(k) .and. m == obs%m(k)) exit lcycle
            add_lm = add_lm + 1
          end do
          l = l + 1
        end do lcycle
        ! Warning: it should not add the nspin components?
        dump = dump + obs%weight(k) * sum(multipole(add_lm, i, 1:nspin))
      end do

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
  SAFE_DEALLOCATE_A(iunit)
  SAFE_DEALLOCATE_A(q)
  SAFE_DEALLOCATE_A(mu)
  SAFE_DEALLOCATE_A(qq)
  SAFE_DEALLOCATE_A(c)
  SAFE_DEALLOCATE_A(ot)
  SAFE_DEALLOCATE_A(multipole)
  SAFE_DEALLOCATE_A(dipole)
end subroutine generate_signal
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine modify_ot(time_steps, dt, order, ot, omega, power)
  use global_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use string_m
  use profiling_m

  implicit none

  integer,             intent(in)    :: time_steps
  FLOAT,               intent(in)    :: dt
  integer,             intent(in)    :: order
  FLOAT,               intent(inout) :: ot(0:time_steps)
  FLOAT,               intent(in)    :: omega
  FLOAT,               intent(in)    :: power

  integer :: i

  select case(mod(order, 2))
  case(1)
    do i = 0, time_steps
      ot(i) = ot(i) - M_TWO * power * sin(omega*i*dt)
    end do
  case(0)
    if(omega .eq. M_ZERO) then
      do i = 0, time_steps
        ot(i) = ot(i) - power * cos(omega*i*dt)
      end do
    else
      do i = 0, time_steps
        ot(i) = ot(i) - M_TWO * power * cos(omega*i*dt)
      end do
    end if
  end select

end subroutine modify_ot
! ---------------------------------------------------------


! ---------------------------------------------------------
subroutine write_ot(nspin, time_steps, dt, kick, units, order, ot, observable)
  use global_m
  use io_m
  use kick_m
  use profiling_m
  use spectrum_m
  use string_m
  use unit_m
  use unit_system_m

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
  type(unit_t) :: ot_unit

  iunit = io_open('ot', action='write', status='replace')

  write(iunit, '(a15,i2)')      '# nspin        ', nspin
  write(iunit, '(a15,i2)')      '# Order        ', order
  write(iunit, '(a28,i2)')      '# Frequencies subtracted = ', 0
  select case(observable(1))
  case(-1)
    write(iunit,'(a)') '# Observable operator = kick operator'
    if(kick%n_multipoles > 0 ) then
      ot_unit = units_out%length**kick%l(1)
    else
      ot_unit = units_out%length
    end if
  case(0)
    select case(observable(2))
    case(1); write(iunit,'(a)') '# O = x'
    case(2); write(iunit,'(a)') '# O = y'
    case(3); write(iunit,'(a)') '# O = z'
    end select
    ot_unit = units_out%length
  case default
    ot_unit = units_out%length**observable(1)
    write(iunit, '(a12,i1,a1,i2,a1)') '# (l, m) = (', observable(1),',',observable(2),')'
  end select
  call kick_write(kick, iunit)

  ! Units
  write(iunit,'(a1)', advance = 'no') '#'
  write(iunit,'(a20)', advance = 'no') str_center('t', 19)
  write(iunit,'(a20)', advance = 'yes') str_center('<O>(t)', 20)
  write(iunit,'(a1)', advance = 'no') '#'
  write(header_string, '(a)') '['//trim(units_abbrev(units_out%time))//']'
  write(iunit,'(a20)', advance = 'no')  str_center(trim(header_string), 19)
  write(header_string, '(a)') '['//trim(units_abbrev(units_out%length))//']'
  write(iunit,'(a20)', advance = 'yes')  str_center(trim(header_string), 20)

  do i = 0, time_steps
    write(iunit, '(2e20.8)') units_from_atomic(units_out%time, i*dt), units_from_atomic(ot_unit, ot(i))
  end do

  call io_close(iunit)
end subroutine write_ot


! ---------------------------------------------------------
subroutine read_ot(nspin, order, nw_subtracted)
  use global_m
  use io_m
  use messages_m
  use datasets_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use string_m
  use oscillator_strength_m
  use profiling_m

  implicit none

  integer, intent(out) :: nspin
  integer, intent(out) :: order
  integer, intent(out) :: nw_subtracted

  integer :: iunit, i, ierr
  character(len=100) :: line
  character(len=12)  :: dummychar
  FLOAT :: dummy, t1, t2
  type(unit_t) :: ot_unit

  iunit = io_open('ot', action='read', status='old')
  if(iunit .eq. 0) then
    write(message(1),'(a)') 'A file called ot should be present and was not found.'
    call messages_fatal(1)
  end if

  read(iunit, '(15x,i2)') nspin
  read(iunit, '(15x,i2)') order
  read(iunit, '(28x,i2)') nw_subtracted
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
    write(message(1),'(a)') 'Problem reading "ot" file: could not figure out the shape'
    write(message(2),'(a)') 'of the observation operator.'
    call messages_fatal(2)
  end if

  call kick_read(kick, iunit)
  read(iunit, '(a)')  line
  read(iunit, '(a)')  line
  call io_skip_header(iunit)

  ! Figure out about the units of the file
  call unit_system_from_file(units, "ot", ierr)
  if(ierr.ne.0) then
    write(message(1), '(a)') 'Could not figure out the units in file "ot".'
    call messages_fatal(1)
  end if

  select case(observable(1))
  case(-1)
    if(kick%n_multipoles > 0) then
      ot_unit = units_out%length**kick%l(1)
    else
      ot_unit = units_out%length
    end if
  case(0)
    ot_unit = units_out%length
  case default
    ot_unit = units_out%length**observable(1)
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
  dt = units_to_atomic(units%time, (t2 - t1)) ! units_out is OK

  call io_skip_header(iunit)

  SAFE_ALLOCATE(ot(0:time_steps))

  do i = 0, time_steps-1
    read(iunit, *) dummy, ot(i)
    ot(i) = units_to_atomic(ot_unit, ot(i))
  end do

end subroutine read_ot


! ---------------------------------------------------------
subroutine print_omega_file(omega, search_interval, final_time, nfrequencies)
  use global_m
  use messages_m
  use datasets_m
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use spectrum_m
  use lalg_adv_m
  use oscillator_strength_m
  use profiling_m

  implicit none

  FLOAT,   intent(inout) :: omega
  FLOAT,   intent(inout) :: search_interval
  FLOAT,   intent(inout) :: final_time
  integer, intent(inout) :: nfrequencies

  integer :: iunit, i, ierr, nspin, order, nw_subtracted
  logical :: file_exists
  character(len=20) :: header_string
  FLOAT, allocatable :: warray(:), tarray(:)
  FLOAT :: leftbound, rightbound, dw, w, aw

  ! First, let us check that the file "ot" exists.
  inquire(file="ot", exist  = file_exists)
  if(.not.file_exists) then
    write(message(1),'(a)') "Could not find 'ot' file."
    call messages_fatal(1)
  end if

  ! Now, we should find out which units the file "ot" has.
  call unit_system_from_file(units, "ot", ierr)
  if(ierr.ne.0) then
    write(message(1),'(a)') "Could not retrieve units in the 'ot' file."
    call messages_fatal(1)
  end if

  if(omega > M_ZERO) then
    omega = units_to_atomic(units%energy, omega)
  else
    omega = M_HALF
  end if

  if(search_interval > M_ZERO) then
    search_interval = units_to_atomic(units%energy, search_interval)
  else
    search_interval = M_HALF
  end if

  leftbound = omega - search_interval
  rightbound = omega + search_interval

  SAFE_ALLOCATE(warray(1:nfrequencies))
  SAFE_ALLOCATE(tarray(1:nfrequencies))

  call read_ot(nspin, order, nw_subtracted)

  if(mod(order, 2).eq.1) then
    mode = SINE_TRANSFORM
  else
    mode = COSINE_TRANSFORM
  end if

  if(final_time > M_ZERO) then
    total_time = units_to_atomic(units%time, final_time)
    if(total_time > dt*time_steps) then
      total_time = dt*time_steps
      write(0, '(a)')        '* WARNING: The requested total time to process is larger than the time'
      write(0, '(a)')        '*          available in the input file.'
      write(0, '(a,f8.4,a)') '           The time has been adjusted to ', &
           units_from_atomic(units_out%time, total_time), trim(units_abbrev(units_out%time))
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
    tarray(i) = aw
  end do

  iunit = io_open('omega', action='write', status='replace')
  write(iunit, '(a15,i2)')      '# nspin        ', nspin
  call kick_write(kick, iunit)
  write(iunit, '(a)') '#%'
  write(iunit, '(a1,a20)', advance = 'no') '#', str_center("omega", 20)
  write(header_string,'(a)') 'F(omega)'
  write(iunit, '(a20)', advance = 'yes') str_center(trim(header_string), 20)
  write(iunit, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20)
  ! Here we should print the units of the transform.
  write(iunit, '(a)', advance = 'yes')

  do i = 1, nfrequencies
    write(iunit,'(2e20.8)') units_from_atomic(units_out%energy, warray(i)), &
                            tarray(i)
  end do

  SAFE_DEALLOCATE_A(warray)
  SAFE_DEALLOCATE_A(tarray)
end subroutine print_omega_file
! ---------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
