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

#include "global.h"

module lasers
use mesh

implicit none

type laser_type
  complex(r8) :: pol(3) ! the polarization of the laser
  real(r8) :: A0     ! the initial amplitude of the laser
  real(r8) :: omega0 ! the average frequency of the laser

  integer  :: envelope

  real(r8) :: t0     ! the maximum of the pulse
  real(r8) :: tau0   ! the width of the pulse
  real(r8) :: tau1   ! for the ramped shape, the length of the "ramping" intervals 
                     ! (tau0 will be the length of the constant interval)

  real(r8), pointer :: numerical(:,:) ! when the laser is numerical
  real(r8) :: dt

  type(spline_type) :: phase      ! when reading a laser from a file
  type(spline_type) :: amplitude

  ! For the velocity gauge
  ! cosine envelope
  complex(r8) :: field_1
  real(r8) :: w1, a1, r1
end type laser_type

contains

subroutine laser_init(m, no_l, l)
  type(mesh_type), intent(IN) :: m
  integer, intent(out) :: no_l
  type(laser_type), pointer :: l(:)

  integer :: i, k, iunit
  real(r8) :: x, y, z
  character(len=80) :: str

  call push_sub('laser_init')

  str = 'TDLasers'
  no_l = oct_parse_block_n(str)
  if(no_l > 0) then
    allocate(l(no_l))

    do i = 1, no_l
      call oct_parse_block_complex(str, i-1, 0, l(i)%pol(1))
      call oct_parse_block_complex(str, i-1, 1, l(i)%pol(2))
      call oct_parse_block_complex(str, i-1, 2, l(i)%pol(3))
      call oct_parse_block_double (str, i-1, 3, l(i)%A0)
      call oct_parse_block_double (str, i-1, 4, l(i)%omega0)
      call oct_parse_block_int    (str, i-1, 5, l(i)%envelope)

      if(l(i)%envelope > 0.and.l(i)%envelope < 4) then
        call oct_parse_block_double (str, i-1, 6, l(i)%tau0)
        call oct_parse_block_double (str, i-1, 7, l(i)%t0)
        l(i)%tau0 = l(i)%tau0 * units_inp%time%factor
        l(i)%t0   = l(i)%t0   * units_inp%time%factor

        if(l(i)%envelope == 3) then
          call oct_parse_block_double (str, i-1, 8, l(i)%tau1)
          l(i)%tau1 = l(i)%tau1 * units_inp%time%factor ! adjust units
        end if
      else if(l(i)%envelope == 10) then ! read from file
        call from_file(i, l(i))
      end if

      ! Adjust units of common variables
      l(i)%A0     = l(i)%A0 * units_inp%energy%factor / units_inp%length%factor
      l(i)%omega0 = l(i)%omega0 * units_inp%energy%factor

      ! normalize polarization
      l(i)%pol(:) = l(i)%pol(:)/sqrt(sum(abs(l(i)%pol(:))**2))

      ! velocity gauge
      if(l(i)%envelope == 2) then
        ! setup some constants
        l(i)%w1 = M_PI/ (2._r8*l(i)%tau0)
        l(i)%a1 = -M_PI/2._r8 * (2._r8 + l(i)%t0/l(i)%tau0)
        l(i)%r1 = 1/((l(i)%w1 - l(i)%omega0)*(l(i)%w1 + l(i)%omega0))
        l(i)%field_1 = -M_zI * l(i)%omega0 * cos(l(i)%a1) &
             - l(i)%w1 * sin(l(i)%a1)
      end if

    end do
  else
    no_l = 0
  end if

  call pop_sub()

contains
  subroutine from_file(i, l)
    integer, intent(in) :: i
    type(laser_type), intent(out) :: l

    integer :: iunit, lines, j
    character(len=50) :: filename
    real(r8) :: dummy
    real(r8), allocatable :: t(:), am(:), ph(:)

    call oct_parse_block_double(str, i-1, 6, l%tau0)
    l%tau0   = l%tau0 * units_inp%time%factor

    ! open file
    call oct_parse_block_string(str, i-1, 7, filename)
    call io_assign(iunit)
    open(iunit, file=trim(filename), status='old', iostat=j)
    if(j.ne.0) then
      write(message(1),'(3a)') "Could not open file '", trim(filename), "'"
      call write_fatal(1)
    endif

    ! count lines in file
    lines = 0
    do
      read(iunit, *, err=100, end=100) dummy, dummy, dummy
      lines = lines + 1
    end do
100 continue
    rewind(iunit)

    ! allocate and read info
    allocate(t(lines), am(lines), ph(lines))
    do j = 1, lines
      read(iunit, *, err=100, end=100) t(j), am(j), ph(j)
    end do
    call io_close(iunit)

    ! build splines
    t = t*l%tau0
    call spline_init(l%amplitude)
    call spline_fit(lines, t, am, l%amplitude)

    call spline_init(l%phase)
    call spline_fit(lines, t, ph, l%phase)

    ! clean
    deallocate(t, am, ph)

  end subroutine from_file
end subroutine laser_init

subroutine laser_end(no_l, l)
  integer, intent(in) :: no_l
  type(laser_type), pointer :: l(:)

  integer :: i

  if(no_l <= 0) return

  do i = 1, no_l
    select case(l(i)%envelope)
      case(10)
        call spline_end(l(i)%amplitude)
        call spline_end(l(i)%phase)
      case(99)
        if(associated(l(i)%numerical)) then
          deallocate(l(i)%numerical); nullify(l(i)%numerical)
        end if
      end select
  end do
  
  deallocate(l); nullify(l)

end subroutine laser_end

subroutine laser_write_info(no_l, l, iunit)
  integer, intent(in) :: no_l, iunit
  type(laser_type), intent(in) :: l(no_l)

  integer :: i

  do i = 1, no_l
    write(iunit,'(i2,a)') i,':'
    write(iunit,'(3x,a,3(a1,f7.4,a1,f7.4,a1))') 'Polarization: ', &
         '(', real(l(i)%pol(1)), ',', aimag(l(i)%pol(1)), '), ', &
         '(', real(l(i)%pol(2)), ',', aimag(l(i)%pol(2)), '), ', &
         '(', real(l(i)%pol(3)), ',', aimag(l(i)%pol(3)), ')'
    write(iunit,'(3x,a,i2)')    'Envelope:     ', l(i)%envelope
    write(iunit,'(3x,a,f10.4,3a)') 'Frequency: ', l(i)%omega0/units_inp%energy%factor, &
         ' [', trim(units_inp%energy%abbrev), ']'
    write(iunit,'(3x,a,f10.4,5a)') 'Amplitude: ', &
         l(i)%A0/units_inp%energy%factor*units_inp%length%factor, &
         ' [', trim(units_inp%energy%abbrev), '/', trim(units_inp%length%abbrev), ']'
    write(iunit, '(3x,a,es14.4,a)') 'Intensity: ', &
         (l(i)%A0*5.14225e9_r8)**2*1.3272e-3_r8, " [W/cm^2]"
    write(iunit,'(3x,a,f10.4,3a)') 'Width:     ', l(i)%tau0/units_inp%time%factor, &
         ' [', trim(units_inp%time%abbrev), ']'
    if(l(i)%envelope > 0.and.l(i)%envelope < 4) then
      write(iunit,'(3x,a,f10.4,3a)') 'Middle t:  ', l(i)%t0/units_inp%time%factor, &
           ' [', trim(units_inp%time%abbrev), ']'
    else if(l(i)%envelope == 10) then
      write(iunit, '(3x,a)') 'Amplitude and phase information read from file'
    end if
    if(l(i)%envelope == 3) then
      write(iunit,'(3x,a,f10.4,3a)') 'Ramp time: ', l(i)%tau1/units_inp%time%factor, &
         ' [', trim(units_inp%time%abbrev), ']'
    end if
  end do

end subroutine laser_write_info

subroutine laser_field(no_l, l, t, field)
  integer,          intent(in)  :: no_l
  type(laser_type), intent(in)  :: l(no_l)
  real(r8),         intent(in)  :: t
  real(r8),         intent(out) :: field(conf%dim)

  complex(r8), allocatable :: amp(:)
  integer :: i

  if(no_l .eq. 0) then
    field(1:conf%dim) = M_ZERO
    return
  end if

  if(no_l == 1 .and. l(1)%envelope == 99) then
    i = int(abs(M_TWO*t/l(1)%dt) + M_HALF) ! steps of dt/2
    field(1:conf%dim) = l(1)%numerical(1:conf%dim, i)
  else
    allocate(amp(no_l));
    call laser_amplitude(no_l, l, t, amp)
    field = 0._r8
    do i = 1, no_l
      field(1:conf%dim) = field(1:conf%dim) + real(amp(i)*l(i)%pol(1:conf%dim))
    end do
    deallocate(amp)
  end if

  return
end subroutine laser_field

!returns the laser vector field A
subroutine laser_vector_field(no_l, l, t, field)
  integer, intent(in) :: no_l
  type(laser_type), intent(IN) :: l(no_l)
  real(r8), intent(in) :: t
  real(r8), intent(out) :: field(conf%dim)

  real(r8) :: r
  complex(r8), allocatable :: amp(:)
  integer :: i

  allocate(amp(no_l));
  call laser_vector_amplitude(no_l, l, t, amp)
  field = 0._r8
  do i = 1, no_l
    field(1:conf%dim) = field(1:conf%dim) + real(amp(i)*l(i)%pol(1:conf%dim))
  end do
  deallocate(amp)

  return
end subroutine laser_vector_field

! The following routines have to be changed in order to add more
! laser envelopes

! returns the amplitude of the electric field
subroutine laser_amplitude(no_l, l, t, amp)
  integer, intent(in) :: no_l
  type(laser_type), intent(in) :: l(no_l)
  real(r8), intent(in) :: t
  complex(r8), intent(out) :: amp(no_l)

  real(r8) :: r, ph
  integer :: i

  do i = 1, no_l
    ph = 0._r8; r = 0._r8
    select case (l(i)%envelope)
    case(1) ! Gaussian
      r = exp(-(t - l(i)%t0)**2 / (2.0_r8*l(i)%tau0**2))
    case(2) ! cosinoidal
      if(abs(t - l(i)%t0) <= l(i)%tau0) then
        r = cos( (M_Pi/2)*((t - 2*l(i)%tau0 - l(i)%t0)/l(i)%tau0) )
      end if
    case(3) ! ramped
      if(t > l(i)%t0-l(i)%tau0/2.0_r8-l(i)%tau1 .and. t <= l(i)%t0-l(i)%tau0/2.0_r8) then
        r = (t - (l(i)%t0 - l(i)%tau0/2.0 - l(i)%tau1))/l(i)%tau1
      elseif(t>l(i)%t0-l(i)%tau0/2.0_r8 .and. t <=l(i)%t0+l(i)%tau0/2.0_r8) then
        r = 1._r8
      elseif(t>l(i)%t0+l(i)%tau0/2.0_r8 .and. t <=l(i)%t0+l(i)%tau0/2.0_r8+l(i)%tau1) then
        r = (l(i)%t0 + l(i)%tau0/2.0_r8 + l(i)%tau1 - t)/l(i)%tau1
      endif
    case(10)
      r  = splint(l(i)%amplitude, t)
      ph = splint(l(i)%phase, t)
    end select

    amp(i) = l(i)%A0 * r * exp(M_zI * (l(i)%omega0*t + ph))
  end do

  return
end subroutine laser_amplitude

! returns the amplitude of the vector field A
subroutine laser_vector_amplitude(no_l, l, t, amp)
  integer, intent(in) :: no_l
  type(laser_type), intent(in) :: l(no_l)
  real(r8), intent(in) :: t
  complex(r8), intent(out) :: amp(no_l)

  complex(r8) :: r, tt
  integer :: i

  do i = 1, no_l
    select case (l(i)%envelope)
    case(1) ! Gaussian
      ! not yet implemented
      r = 0._r8
    case(2) ! cosinoidal
      if(t < l(i)%t0 - l(i)%tau0) then
        tt = l(i)%t0 - l(i)%tau0
      else if(abs(t-l(i)%t0) <= l(i)%tau0) then
        tt = t
      else
        tt = l(i)%t0 + l(i)%tau0
      end if

      r = - l(i)%r1 * (l(i)%field_1 + exp(M_zI * l(i)%omega0*tt) * &
           (M_zI*l(i)%omega0*cos(l(i)%a1 + l(i)%w1*tt) + l(i)%w1*sin(l(i)%a1 + l(i)%w1*tt)))
    case(3) ! ramp
      ! not yet implemented
      r = 0._r8
    case default
      r = 0._r8
    end select
    amp(i) = l(i)%A0 * r
  end do

  return
end subroutine laser_vector_amplitude

end module lasers
