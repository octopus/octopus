!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module spectrum
use global
use messages
use lib_oct_parser
use io
use units

implicit none

integer, parameter :: SPECTRUM_DAMP_NONE       = 0, &
                      SPECTRUM_DAMP_LORENTZIAN = 1, &
                      SPECTRUM_DAMP_POLYNOMIAL = 2, &
                      SPECTRUM_DAMP_GAUSSIAN   = 3

type spec_type
  FLOAT :: start_time  ! start time for the transform
  FLOAT :: end_time    ! when to stop the transform
  FLOAT :: energy_step ! step in energy mesh
  FLOAT :: min_energy  ! maximum of energy mesh
  FLOAT :: max_energy  ! maximum of energy mesh
end type spec_type

! For the strength function
type spec_sf
  ! input
  integer  :: transform        ! The Fourier transform to perform (sin, cos)
  integer  :: damp             ! Damp type (none, exp or pol)
  FLOAT :: damp_factor      ! factor used in damping
  FLOAT :: delta_strength   ! strength of the delta perturbation

  ! output
  integer :: no_e, nspin ! dimensions of sp
  FLOAT, pointer :: sp(:,:) ! do not forget to deallocate this
  FLOAT :: ewsum  ! electronic sum rule
  FLOAT :: alpha  ! Polariz. (sum rule)
  FLOAT :: alpha2 ! Polariz. (F.T.)
end type spec_sf

! For the rotational strength function
type spec_rsf
  ! input
  integer  :: damp             ! Damp type (none, exp or pol)
  FLOAT :: damp_factor      ! factor used in damping
  FLOAT :: delta_strength   ! strength of the delta perturbation
  FLOAT :: pol(3)           ! Direction of the perturbation

  ! output
  integer :: no_e
  CMPLX, pointer :: sp(:) ! do not forget to deallocate this
end type spec_rsf

type spec_sh
  ! input
  character :: pol

  ! output
  integer :: no_e ! dimensions of sp
  FLOAT, pointer :: sp(:) ! do not forget to deallocate this
end type spec_sh

contains

subroutine spectrum_strength_function(out_file, s, sf, print_info)
  character(len=*), intent(in) :: out_file
  type(spec_type), intent(inout) :: s
  type(spec_sf), intent(inout) :: sf
  logical, intent(in) :: print_info

  integer :: iunit, i, is, ie, &
      ntiter, j, jj, k, isp, time_steps
  FLOAT :: dump, dt, x
  FLOAT, allocatable :: dumpa(:)
  FLOAT, allocatable :: dipole(:,:)

  call spectrum_mult_info(iunit, sf%nspin, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

  ! load dipole from file
  allocate(dipole(0:time_steps, sf%nspin))
  do i = 0, time_steps
    read(iunit, *) j, dump, dipole(i,:)
    dipole(i,:) = dipole(i,:) * units_out%length%factor
  end do
  call io_close(iunit)

  ! subtract static dipole
  do i = 1, sf%nspin
     dipole(:, i) = dipole(:, i) - dipole(0, i)
  end do

  sf%no_e = (s%max_energy - s%min_energy) / s%energy_step
  allocate(sf%sp(0:sf%no_e, sf%nspin))
  sf%sp = M_ZERO
  sf%alpha = M_ZERO; sf%alpha2 = M_ZERO; sf%ewsum = M_ZERO

  ! Gets the damping function (here because otherwise it is awfully slow in "pol" mode...)
  allocate(dumpa(is:ie))
  do j = is, ie
     jj = j - is
      select case(sf%damp)
      case(SPECTRUM_DAMP_NONE)
        dumpa(j) = M_ONE
      case(SPECTRUM_DAMP_LORENTZIAN)
        dumpa(j)= exp(-jj*dt*sf%damp_factor)
      case(SPECTRUM_DAMP_POLYNOMIAL)
        dumpa(j) = M_ONE - M_THREE*(real(jj)/ntiter)**2                          &
            + M_TWO*(real(jj)/ntiter)**3
      case(SPECTRUM_DAMP_GAUSSIAN)
        dumpa(j)= exp(-(jj*dt)**2*sf%damp_factor**2)
      end select
  enddo

  do k = 0, sf%no_e
    do j = is, ie
      
      jj = j - is
      
      select case(sf%transform)
      case(1)
        x = sin((k*s%energy_step + s%min_energy)*jj*dt)
      case(2)
        x = cos((k*s%energy_step + s%min_energy)*jj*dt)
      end select

      do isp = 1, sf%nspin
        sf%sp(k, isp) = sf%sp(k, isp) + x*dumpa(j)*dipole(j, isp)

        ! polarizability sum rule
        if(k == 0) then
          sf%alpha2 = sf%alpha2 + dumpa(j)*dipole(j, isp)
        end if
      end do

    end do
    sf%sp(k, :) = sf%sp(k, :)*dt

    ! calculate strength function
    sf%sp(k, :) = (sf%sp(k, :)*s%energy_step*k*M_TWO)/(M_pi*sf%delta_strength)

    do isp = 1, sf%nspin
      if(k.ne.0) then
        sf%alpha = sf%alpha + (sf%sp(k, isp)/(k*s%energy_step)**2)*s%energy_step
      endif
      sf%ewsum = sf%ewsum + sf%sp(k, isp)*s%energy_step
    end do
  end do

  sf%alpha2 = sf%alpha2 * dt / sf%delta_strength
  deallocate(dipole, dumpa)

  ! output
  if(trim(out_file) .ne. '-') then
    iunit = io_open(out_file, action='write')

    ! should output units, etc...
    do i = 0, sf%no_e
      write(iunit,'(5e15.6)') (i*s%energy_step + s%min_energy) / units_out%energy%factor, &
           sf%sp(i, :) * units_out%energy%factor
    end do
    call io_close(iunit)
  end if

  ! print some info
  if(print_info) then
    write(message(1), '(a,i8)')    'Number of spin       = ', sf%nspin
    write(message(2), '(a,i8)')    'Number of time steps = ', ntiter
    write(message(3), '(a,i4)')    'SpecTransformMode    = ', sf%transform
    write(message(4), '(a,i4)')    'SpecDampMode         = ', sf%damp
    write(message(5), '(a,f10.4)') 'SpecDampFactor       = ', sf%damp_factor * units_out%time%factor
    write(message(6), '(a,f10.4)') 'SpecStartTime        = ', s%start_time   / units_out%time%factor
    write(message(7), '(a,f10.4)') 'SpecEndTime          = ', s%end_time     / units_out%time%factor
    write(message(8), '(a,f10.4)') 'SpecMinEnergy        = ', s%min_energy   / units_inp%energy%factor
    write(message(9), '(a,f10.4)') 'SpecMaxEnergy        = ', s%max_energy   / units_inp%energy%factor
    write(message(10),'(a,f10.4)') 'SpecEnergyStep       = ', s%energy_step  / units_inp%energy%factor
    call write_info(10)

    message(1) = ""
    write(message(2),'(a,f16.6)') 'Electronic sum rule  = ', sf%ewsum
    write(message(3),'(a,f16.6)') 'Polariz. (sum rule)  = ', sf%alpha  / units_inp%length%factor**3
    write(message(4),'(a,f16.6)') 'Polariz. (F.T.)      = ', sf%alpha2 / units_inp%length%factor**3
    call write_info(4)
  end if

  return
end subroutine spectrum_strength_function

subroutine spectrum_rotatory_strength(out_file, s, rsf, print_info)
  character(len=*), intent(in) :: out_file
  type(spec_type), intent(inout) :: s
  type(spec_rsf), intent(inout) :: rsf
  logical, intent(in) :: print_info

  integer :: iunit, i, is, ie, &
      ntiter, j, jj, k, time_steps
  FLOAT :: dump, dt
  CMPLX :: z
  FLOAT, allocatable :: dumpa(:)
  FLOAT, allocatable :: angular(:, :)

  call spectrum_file_info("td.general/angular", iunit, time_steps, dt, i)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

  ! load dipole from file
  allocate(angular(0:time_steps, 3))
  do i = 0, time_steps
    read(iunit, *) j, dump, angular(i, 1:3)
    !dipole(i,:) = dipole(i,:) * units_out%length%factor
  end do
  call io_close(iunit)

  ! subtract static dipole
  do i = 1, 3
     angular(:, i) = angular(:, i) - angular(0, i)
  end do

  rsf%no_e = (s%max_energy - s%min_energy) / s%energy_step
  allocate(rsf%sp(0:rsf%no_e))
  rsf%sp = M_z0

  ! Gets the damping function (here because otherwise it is awfully slow in "pol" mode...)
  allocate(dumpa(is:ie))
  do j = is, ie
     jj = j - is
      select case(rsf%damp)
      case(SPECTRUM_DAMP_NONE)
        dumpa(j) = M_ONE
      case(SPECTRUM_DAMP_LORENTZIAN)
        dumpa(j)= exp(-jj*dt*rsf%damp_factor)
      case(SPECTRUM_DAMP_POLYNOMIAL)
        dumpa(j) = 1.0 - 3.0*(real(jj)/ntiter)**2                          &
            + 2.0*(real(jj)/ntiter)**3
      case(SPECTRUM_DAMP_GAUSSIAN)
        dumpa(j)= exp(-(jj*dt)**2*rsf%damp_factor**2)
      end select
  enddo

  do k = 0, rsf%no_e
    do j = is, ie
      
      jj = j - is

      z = exp(M_zI*(k*s%energy_step + s%min_energy)*jj*dt)
      rsf%sp(k) = rsf%sp(k) + z*dumpa(j)*sum(angular(j, :)*rsf%pol(:))

    end do
    rsf%sp(k) = rsf%sp(k)*dt

  end do

  deallocate(angular, dumpa)

  ! output
  if(trim(out_file) .ne. '-') then
    iunit = io_open(out_file, action='write')

    ! should output units, etc...
    do i = 0, rsf%no_e
      write(iunit,'(5e15.6)') (i*s%energy_step + s%min_energy) / units_out%energy%factor, &
           rsf%sp(i) * (units_out%length%factor)**3
    end do
    call io_close(iunit)
  end if

  ! print some info
  if(print_info) then
    write(message(1), '(a,i8)')    'Number of time steps = ', ntiter
    write(message(2), '(a,i4)')    'SpecDampMode         = ', rsf%damp
    write(message(3), '(a,f10.4)') 'SpecDampFactor       = ', rsf%damp_factor * units_out%time%factor
    write(message(4), '(a,f10.4)') 'SpecStartTime        = ', s%start_time   / units_out%time%factor
    write(message(5), '(a,f10.4)') 'SpecEndTime          = ', s%end_time     / units_out%time%factor
    write(message(6), '(a,f10.4)') 'SpecMinEnergy        = ', s%min_energy   / units_inp%energy%factor
    write(message(7), '(a,f10.4)') 'SpecMaxEnergy        = ', s%max_energy   / units_inp%energy%factor
    write(message(8),'(a,f10.4)')  'SpecEnergyStep       = ', s%energy_step  / units_inp%energy%factor
    call write_info(8)
  end if

  return
end subroutine spectrum_rotatory_strength

subroutine spectrum_hs_from_mult(out_file, s, sh, print_info)
  character(len=*), intent(in) :: out_file
  type(spec_type), intent(inout) :: s
  type(spec_sh), intent(inout) :: sh
  logical, intent(in) :: print_info

  integer :: i, j, iunit, nspin, time_steps, is, ie, ntiter
  FLOAT :: dt, dump
  FLOAT, allocatable :: d(:,:)
  CMPLX :: c
  CMPLX, allocatable :: dipole(:), ddipole(:)

  call spectrum_mult_info(iunit, nspin, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)
  
  ! load dipole from file
  allocate(dipole(0:time_steps))
  allocate(d(3, nspin))
  do i = 1, time_steps
    read(iunit, *) j, dump, dump, dump, d
    select case(sh%pol)
    case('x')
      dipole(i) = -sum(d(3, :))
    case('y')
      dipole(i) = -sum(d(1, :))
    case('z')
      dipole(i) =  sum(d(2, :))
    case('+')
      dipole(i) = -sum(d(3, :) + M_zI*d(1, :)) / sqrt(M_TWO)
    case('-')
      dipole(i) = -sum(d(3, :) - M_zI*d(1, :)) / sqrt(M_TWO)
    end select
    dipole(i) = dipole(i) * units_out%length%factor * sqrt(M_FOUR*M_PI/M_THREE)
  end do
  deallocate(d)

  ! we now calculate the first time derivative
  allocate(ddipole(0:time_steps))
  ddipole(0) = (dipole(1) - dipole(0))/dt
  do i = 1, time_steps - 1
    ddipole(i) = (dipole(i + 1) - dipole(i - 1))/(M_TWO*dt)
  end do
  ddipole(time_steps) = (dipole(time_steps) - dipole(time_steps - 1))/dt

  ! and the second time derivative
  dipole(0) = (ddipole(1) - ddipole(0))/dt
  do i = 1, time_steps - 1
    dipole(i) = (ddipole(i + 1) - ddipole(i - 1))/(M_TWO*dt)
  end do
  dipole(time_steps) = (ddipole(time_steps) - ddipole(time_steps - 1))/dt
  deallocate(ddipole)

  ! now we Fourier transform
  sh%no_e = (s%max_energy - s%min_energy) / s%energy_step
  allocate(sh%sp(0:sh%no_e))
  sh%sp = M_ZERO
  
  do i = 0, sh%no_e
    c = M_z0
    do j = is, ie
      c = c + exp(M_zI*(i*s%energy_step + s%min_energy)*j*dt)*dipole(j)
    end do
    sh%sp(i) = abs(c)**2
  end do
  deallocate(dipole)

  ! output
  if(trim(out_file) .ne. '-') then
    iunit = io_open(trim(out_file) // "." // trim(sh%pol), action='write')

    ! should output units, etc...
    do i = 0, sh%no_e
      write(iunit,'(5e15.6)') (i*s%energy_step + s%min_energy) / units_out%energy%factor, &
           sh%sp(i) * units_out%energy%factor
    end do
    call io_close(iunit)
  end if

end subroutine spectrum_hs_from_mult

subroutine spectrum_hs_from_acc(out_file, s, sh, print_info)
  character(len=*), intent(in) :: out_file
  type(spec_type), intent(inout) :: s
  type(spec_sh), intent(inout) :: sh
  logical, intent(in) :: print_info

  integer :: i, j, iunit, time_steps, is, ie, ntiter
  FLOAT :: dt, dummy, a(3)
  CMPLX, allocatable :: acc(:)
  CMPLX :: c
  
  call spectrum_acc_info(iunit, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)
  
  ! load dipole from file
  allocate(acc(0:time_steps))
  acc = M_ZERO
  do i = 1, time_steps
    read(iunit, *) j, dummy, a
    select case(sh%pol)
    case('x')
      acc(i) = a(1)
    case('y')
      acc(i) = a(2)
    case('z')
      acc(i) = a(3)
    case('+')
      acc(i) = (a(1) + M_zI*a(2)) / sqrt(M_TWO)
    case('-')
      acc(i) = (a(1) - M_zI*a(2)) / sqrt(M_TWO)
    end select
    acc(i) = acc(i) * units_out%acceleration%factor
  end do

  ! now we Fourier transform
  sh%no_e = (s%max_energy - s%min_energy) / s%energy_step
  allocate(sh%sp(0:sh%no_e))
  sh%sp = M_ZERO
  
  do i = 0, sh%no_e
    c = M_z0
    do j = is, ie
      c = c + exp(M_zI*(i*s%energy_step + s%min_energy)*j*dt)*acc(j)
    end do
    sh%sp(i) = abs(c)**2
  end do
  deallocate(acc)

  ! output
  if(trim(out_file) .ne. '-') then
    iunit = io_open(trim(out_file) // "." // trim(sh%pol), action='write')

    ! should output units, etc...
    do i = 0, sh%no_e
      write(iunit,'(5e15.6)') (i*s%energy_step + s%min_energy) / units_out%energy%factor, &
           sh%sp(i) * units_out%energy%factor
    end do
    call io_close(iunit)
  end if

end subroutine spectrum_hs_from_acc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finds out some information about file "file", supposed to contain time
! dependent output information (e.g., multipoles, angular, coordinates, acc,...
! It opens the file for reading, assigning it unit "iunit".
! In principle, it should substitute both spectrum_mult_info and
! spectrum_acc_info.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectrum_file_info(file, iunit, time_steps, dt, n)
  character(len=*), intent(in) :: file
  integer, intent(out) :: iunit, n, time_steps
  FLOAT, intent(out) :: dt

  integer :: i, j
  FLOAT :: t1, t2, dummy

  iunit = io_open(file, action='read', status='old')

  ! First line may contain some informative integer n (or not...)
  read(iunit,'(10x,i2)', iostat = i) n
  if(i.ne.0) n = 0

  ! Skip header.
  read(iunit,*); read(iunit,*)

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
100 continue
  dt = (t2 - t1) * units_out%time%factor ! units_out is OK
  time_steps = time_steps - 1
  
  if(time_steps < 3) then
    write(message(1),'(a,a,a)') "Empty ", trim(adjustl(file)),"?"
    call write_fatal(1)
  end if

  rewind(iunit)
  read(iunit, *); read(iunit, *); read(iunit, *) ! skip header
end subroutine spectrum_file_info

subroutine spectrum_mult_info(iunit, nspin, time_steps, dt)
  integer, intent(out) :: iunit, nspin, time_steps
  FLOAT, intent(out) :: dt

  integer :: i, j
  FLOAT :: t1, t2, dummy

  ! open files
  call io_assign(iunit)
  iunit = io_open('multipoles', action='read', status='old', die=.false.)
  if(iunit < 0) then
    iunit = io_open('td.general/multipoles', action='read', status='old')
  end if
  
  ! read in dipole
  read(iunit, '(10x,i2)') nspin
  read(iunit, *); read(iunit, *) ! skip header

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
100 continue
  dt = (t2 - t1) * units_out%time%factor ! units_out is OK
  time_steps = time_steps - 1
  
  if(time_steps < 3) then
    message(1) = "Empty multipole file?"
    call write_fatal(1)
  end if

  rewind(iunit)
  read(iunit, *); read(iunit, *); read(iunit, *) ! skip header

end subroutine spectrum_mult_info

subroutine spectrum_acc_info(iunit, time_steps, dt)
  integer, intent(out) :: iunit, time_steps
  FLOAT, intent(out) :: dt

  integer :: i, j
  FLOAT :: t1, t2, dummy

  ! open files
  iunit = io_open('acceleration', action='read', status='old', die=.false.)
  if(iunit < 0) then
    iunit = io_open('td.general/acceleration', action='read', status='old')
  endif
  
  ! read in dipole
  read(iunit, *); read(iunit, *) ! skip header

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
100 continue
  dt = (t2 - t1) * units_out%time%factor ! units_out is OK
  time_steps = time_steps - 1
  
  if(time_steps < 3) then
    message(1) = "Empty multipole file?"
    call write_fatal(1)
  end if

  rewind(iunit)
  read(iunit, *); read(iunit, *) ! skip header

end subroutine spectrum_acc_info

subroutine spectrum_fix_time_limits(time_steps, dt, start_time, end_time, is, ie, ntiter)
  integer, intent(in) :: time_steps
  FLOAT, intent(in) :: dt
  FLOAT, intent(inout) :: start_time, end_time
  integer, intent(out) :: is, ie, ntiter

  FLOAT :: ts, te, dummy

  ts = M_ZERO; te = time_steps*dt
  if(start_time < ts) start_time = ts
  if(start_time > te) start_time = te
  if(end_time   > te .or. end_time <= M_ZERO) end_time   = te
  if(end_time   < ts) end_time   = ts

  if(end_time < start_time) then
    dummy = end_time ! swap
    end_time = start_time
    start_time = dummy
  end if
  is = int(start_time/dt)
  ie = int(end_time/dt)
  ntiter = ie - is + 1

end subroutine spectrum_fix_time_limits

end module spectrum
