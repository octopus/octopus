#include "config.h"

module lasers
use fdf
use global
use units
use mesh

implicit none

type laser_type
  complex(r8) :: pol(3) ! the polarization of the laser
  real(r8) :: A0     ! the initial amplitude of the laser
  real(r8) :: tau0   ! the width of the (gaussian) pulse
  real(r8) :: t0     ! the maximum of the pulse
  real(r8) :: tau1   ! for the ramped shape, the length of the "ramping" intervals 
                     ! (tau0 will be the length of the constant interval)
  real(r8) :: omega0 ! the average frequency of the laser

  integer  :: envelope

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
  character(len=256) :: c

  sub_name = 'laser_init'; call push_sub()

  no_l = fdf_integer('TDNumberLasers', 0)
  if(no_l > 0) then
    allocate(l(no_l))

    if (.not.fdf_block('TDLasers', iunit)) then ! read in info
      message(1) = 'Expecting TDLasers block'
      call write_fatal(1)
    end if

    do i = 1, no_l
      ! This read *sucks*, becose FORTRAN *sucks*
      ! I give up, if you know a better way, let me know.
      read(iunit,'(a)') c
      read(c, *) l(i)%pol(:), l(i)%A0, l(i)%envelope

      if(l(i)%envelope == 3) then
        read(c, *) l(i)%pol(:), l(i)%A0, l(i)%envelope, &
             l(i)%tau0, l(i)%t0, l(i)%omega0, l(i)%tau1
        l(i)%tau1 = l(i)%tau1 * units_inp%time%factor ! adjust units
      else
        read(c, *) l(i)%pol(:), l(i)%A0, l(i)%envelope,&
             l(i)%tau0, l(i)%t0, l(i)%omega0
      end if

      ! Adjust units of common variables
      l(i)%A0     = l(i)%A0 * units_inp%energy%factor / units_inp%length%factor
      l(i)%tau0   = l(i)%tau0   * units_inp%time%factor
      l(i)%t0     = l(i)%t0     * units_inp%time%factor
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
end subroutine laser_init

subroutine laser_end(no_l, l)
  integer, intent(in) :: no_l
  type(laser_type), pointer :: l(:)

  integer :: i

  if(no_l <= 0) return
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
    write(iunit,'(3x,a,f10.4,5a)') 'Amplitude: ', &
         l(i)%A0/units_inp%energy%factor*units_inp%length%factor, &
         ' [', trim(units_inp%energy%abbrev), '/', trim(units_inp%length%abbrev), ']'
    write(iunit,'(3x,a,f10.4,3a)') 'Width:     ', l(i)%tau0/units_inp%time%factor, &
         ' [', trim(units_inp%time%abbrev), ']'
    write(iunit,'(3x,a,f10.4,3a)') 'Middle t:  ', l(i)%t0/units_inp%time%factor, &
         ' [', trim(units_inp%time%abbrev), ']'
    if(l(i)%envelope == 3) then
      write(iunit,'(3x,a,f10.4,3a)') 'Ramp time: ', l(i)%tau1/units_inp%time%factor, &
         ' [', trim(units_inp%time%abbrev), ']'
    end if
    write(iunit,'(3x,a,f10.4,3a)') 'Frequency: ', l(i)%omega0/units_inp%energy%factor, &
         ' [', trim(units_inp%energy%abbrev), ']'
  end do

end subroutine laser_write_info

subroutine laser_field(no_l, l, t, field)
  integer, intent(in) :: no_l
  type(laser_type), intent(IN) :: l(no_l)
  real(r8), intent(in) :: t
  real(r8), intent(out) :: field(3)

  real(r8) :: r
  complex(r8), allocatable :: amp(:)
  integer :: i

  if(no_l .eq. 0) then
    field = 0._r8
    return
  end if

  allocate(amp(no_l));
  call laser_amplitude(no_l, l, t, amp)
  field = 0._r8
  do i = 1, no_l
    field(1:3) = field(1:3) + real(amp(i)*l(i)%pol(1:3))
  end do
  deallocate(amp)

  return
end subroutine laser_field

!returns the laser vector field A
subroutine laser_vector_field(no_l, l, t, field)
  integer, intent(in) :: no_l
  type(laser_type), intent(IN) :: l(no_l)
  real(r8), intent(in) :: t
  real(r8), intent(out) :: field(3)

  real(r8) :: r
  complex(r8), allocatable :: amp(:)
  integer :: i

  allocate(amp(no_l));
  call laser_vector_amplitude(no_l, l, t, amp)
  field = 0._r8
  do i = 1, no_l
    field(1:3) = field(1:3) + real(amp(i)*l(i)%pol(1:3))
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

  complex(r8) :: r
  integer :: i

  do i = 1, no_l
    select case (l(i)%envelope)
    case(1) ! Gaussian
      r = exp(-(t - l(i)%t0)**2 / (2.0_r8*l(i)%tau0**2))
    case(2) ! cosinoidal
      if(abs(t - l(i)%t0) <= l(i)%tau0) then
        r = cos( (M_Pi/2)*((t - 2*l(i)%tau0 - l(i)%t0)/l(i)%tau0) )
      else
        r = 0.0_r8
      endif
    case(3) ! ramped
      r = 0.0_r8
      if(t > l(i)%t0-l(i)%tau0/2.0_r8-l(i)%tau1 .and. t <= l(i)%t0-l(i)%tau0/2.0_r8) then
        r = (t - (l(i)%t0 - l(i)%tau0/2.0 - l(i)%tau1))/l(i)%tau1
      elseif(t>l(i)%t0-l(i)%tau0/2.0_r8 .and. t <=l(i)%t0+l(i)%tau0/2.0_r8) then
        r = 1._r8
      elseif(t>l(i)%t0+l(i)%tau0/2.0_r8 .and. t <=l(i)%t0+l(i)%tau0/2.0_r8+l(i)%tau1) then
        r = (l(i)%t0 + l(i)%tau0/2.0_r8 + l(i)%tau1 - t)/l(i)%tau1
      endif
    case default
      r = 0.0_r8
    end select

    amp(i) = l(i)%A0 * r * exp(M_zI * l(i)%omega0*t)
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
