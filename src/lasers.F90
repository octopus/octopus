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

module lasers_m
  use global_m
  use messages_m
  use datasets_m
  use io_m
  use lib_oct_parser_m
  use lib_oct_gsl_spline_m
  use units_m
  use mesh_m
  use simul_box_m
  use output_m

  implicit none

  private
  public ::                       &
    laser_t,                      &
    laser_init,                   &
    laser_end,                    &
    laser_write_info,             &
    laser_field,                  &
    laser_vector_field,           &
    laser_potential

  integer, parameter ::           &
    ENVELOPE_CONSTANT      =  0,  &
    ENVELOPE_GAUSSIAN      =  1,  &
    ENVELOPE_COSINOIDAL    =  2,  &
    ENVELOPE_TRAPEZOIDAL   =  3,  &
    ENVELOPE_FROM_FILE     = 10,  &
    ENVELOPE_NUMERICAL     = 99
 
  integer, parameter ::           &
    SPATIAL_PART_DIPOLE    =  1,  &
    SPATIAL_PART_FROM_FILE =  2

  type laser_t
    CMPLX :: pol(MAX_DIM) ! the polarization of the laser
    FLOAT :: A0           ! the initial amplitude of the laser
    FLOAT :: omega0       ! the average frequency of the laser

    integer  :: envelope
    integer  :: spatial_part

    FLOAT :: t0     ! the maximum of the pulse
    FLOAT :: tau0   ! the width of the pulse
    FLOAT :: tau1   ! for the ramped shape, the length of the "ramping" intervals
                    ! (tau0 will be the length of the constant interval)

    FLOAT, pointer :: numerical(:,:)    ! when the laser is numerical
    FLOAT :: dt

    type(loct_spline_t) :: phase     ! when reading a laser from a file
    type(loct_spline_t) :: amplitude1
    type(loct_spline_t) :: amplitude2
    type(loct_spline_t) :: amplitude3

    FLOAT, pointer :: v(:)              ! the spatial part.

    ! For the velocity gauge
    ! cosine envelope
    CMPLX :: field_1
    FLOAT :: w1, a1, r1
  end type laser_t

contains


  ! ---------------------------------------------------------
  subroutine laser_init(no_l, l, m)
    integer,     intent(out) :: no_l
    type(laser_t),   pointer :: l(:)
    type(mesh_t), intent(in) :: m

    integer(POINTER_SIZE) :: blk
    integer               :: i, no_c
    character(len=50)     :: filename1, filename2

    call push_sub('lasers.laser_init')

    !%Variable TDLasers
    !%Type block
    !%Section Time Dependent
    !%Description
    !% The block TDLasers describe the type and shape of time-dependent external perturbations
    !% that are applied to the system.
    !% Each line of the block describes a laser field; this way you can actually have more
    !% than one laser (e.g. a "pump" and a "probe").
    !% The syntax of each line is, then:
    !%
    !% <tt>%TDLasers
    !% <br>&nbsp;&nbsp;nx | ny | nz | amplitude | omega | envelope | tau0 | t0 | tau1 | filename1 | filename2
    !% <br>%</tt>
    !%
    !% The first three (possibly complex) numbers mark the polarization direction of the
    !% field. The "amplitude" is obviously the amplitude of the field. The "omega" is the
    !% frequency. The "envelope" decides the shape of the enveloping function.
    !% "tau0", "t0" and "tau1" are three paramenters that decide on the
    !% temporal shape of the pulse -- the exact details depend on the particular envelope.
    !%
    !% If the envelope is given in a file, this will be "filename1". If the spatial part
    !% of the field is given in a file, this will be "filename2".
    !%
    !% The last three columns ("tau1", "filename1" and "filename2") are optional; they will
    !% only be searched if needed.
    !%
    !% In order to give the spatial shape of the field in a file, the component "filename2"
    !% has to be present. If it is not present, then the laser field will be a dipolar field
    !% (which is the usual case).
    !%
    !% The possible options for the envelope are:
    !%
    !%Option envelope_constant 0
    !% The envelope is just the unit function. The laser is not a pulse, but a continuous wave (cw).
    !%Option envelope_gaussian 1
    !% The envelope is a Gaussian function. To fully determine its shape, you need to specify the width (tau0)
    !% and the peak time (t0).
    !%Option envelope_cosinusoidal 2
    !% The envelope is a half-cycle of a cosine function. To fully determine its shape, you need to specify the width (tau0)
    !% and the peak time (t0).
    !%Option envelope_trapezoidal 3
    !% The envelope looks like a trapezoid: it starts at zero, ramps linearly until one during tau1 units of time, stays
    !% at that maximum for tau0 units of time, and then decays linearly to zero during tau1 units of time.
    !%Option envelope_fromfile 10
    !% The shape of the laser pulse is read from a file, indicated in the field "filename1".
    !%Option envelope_numerical 99
    !% The laser shape is generated internally. This option is used if the code is run in the optimal-control mode.
    !%End

    no_l = 0
    if(loct_parse_block(check_inp('TDLasers'), blk) == 0) then
      no_l = loct_parse_block_n(blk)
      ALLOCATE(l(no_l), no_l)

      ! The structure of the block is:
      ! nx | ny | nz | amplitude | omega | envelope | tau0 | t0 | tau1 | filename1 | filename2
      do i = 1, no_l
        no_c = loct_parse_block_cols(blk, i-1)
        call loct_parse_block_cmplx(blk, i-1, 0, l(i)%pol(1))
        call loct_parse_block_cmplx(blk, i-1, 1, l(i)%pol(2))
        call loct_parse_block_cmplx(blk, i-1, 2, l(i)%pol(3))
        call loct_parse_block_float(blk, i-1, 3, l(i)%A0)
        call loct_parse_block_float(blk, i-1, 4, l(i)%omega0)
        call loct_parse_block_int  (blk, i-1, 5, l(i)%envelope)
        call loct_parse_block_float(blk, i-1, 6, l(i)%tau0)
        call loct_parse_block_float(blk, i-1, 7, l(i)%t0)
        if(no_c>8)  call loct_parse_block_float(blk, i-1, 8, l(i)%tau1)
        if(no_c>9) call loct_parse_block_string(blk, i-1, 9, filename1)
        if(no_c>10) call loct_parse_block_string(blk, i-1, 10, filename2)

        l(i)%A0     = l(i)%A0 * units_inp%energy%factor / units_inp%length%factor
        l(i)%omega0 = l(i)%omega0 * units_inp%energy%factor
        l(i)%tau0 = l(i)%tau0 * units_inp%time%factor
        l(i)%t0   = l(i)%t0   * units_inp%time%factor
        l(i)%tau1 = l(i)%tau1 * units_inp%time%factor

        if(l(i)%envelope == ENVELOPE_TRAPEZOIDAL .and. no_c < 9)  &
             call input_error('TDLasers')
        if(l(i)%envelope == ENVELOPE_FROM_FILE   .and. no_c < 10) &
             call input_error('TDLasers')
        if(no_c>10) l(i)%spatial_part = SPATIAL_PART_FROM_FILE

        if(l(i)%envelope == ENVELOPE_FROM_FILE)   &
             call get_envelope_from_file(l(i), trim(filename1))
        if(l(i)%spatial_part == SPATIAL_PART_FROM_FILE) &
             call get_spatial_part_from_file(l(i), trim(filename2))

        ! normalize polarization
        l(i)%pol(:) = l(i)%pol(:)/sqrt(sum(abs(l(i)%pol(:))**2))

        ! velocity gauge
        if(l(i)%envelope == ENVELOPE_COSINOIDAL) then
          ! setup some constants
          l(i)%w1 = M_PI/ (M_TWO*l(i)%tau0)
          l(i)%a1 = -M_PI/M_TWO * (M_TWO + l(i)%t0/l(i)%tau0)
          l(i)%r1 = 1/((l(i)%w1 - l(i)%omega0)*(l(i)%w1 + l(i)%omega0))
          l(i)%field_1 = -M_zI * l(i)%omega0 * cos(l(i)%a1) &
            - l(i)%w1 * sin(l(i)%a1)
        end if

      end do

      call loct_parse_block_end(blk)
    end if

    call pop_sub()

  contains


    ! ---------------------------------------------------------
    subroutine get_spatial_part_from_file(l, filename)
      type(laser_t), intent(inout) :: l
      character(len=*), intent(in) :: filename

      integer :: ierr

      ALLOCATE(l%v(m%np), m%np)
      call dinput_function(filename, m, l%v(:), ierr)
      if(ierr < 0) then
        write(message(1),'(a)') "Could not read the external field spatial part."
        write(message(2),'(a)') "Expected in file '"//trim(filename)//"'"
        call write_fatal(2)
      end if

    end subroutine get_spatial_part_from_file


    ! ---------------------------------------------------------
    subroutine get_envelope_from_file(l, filename)
      type(laser_t), intent(inout) :: l
      character(len=*), intent(in) :: filename

      integer :: iunit, lines, j
      FLOAT :: dummy
      FLOAT, allocatable :: t(:), am1(:), am2(:), am3(:)

      iunit = io_open(filename, action='read', status='old')

      ! count lines in file
      lines = 0
      do
        read(iunit, *, err=100, end=100) dummy, dummy, dummy
        lines = lines + 1
      end do
100   continue
      rewind(iunit)

      ! allocate and read info
      ALLOCATE( t(lines), lines)
      ALLOCATE(am1(lines), lines)
      ALLOCATE(am2(lines), lines)
      ALLOCATE(am3(lines), lines)
      am1(:) = M_ZERO
      am2(:) = M_ZERO
      am3(:) = M_ZERO
      
      SELECT CASE(m%sb%dim)
      CASE(1)
         do j = 1, lines
            read(iunit,*, err=100, end=100) t(j), am1(j)
         end do
      CASE(2)
         do j = 1, lines
            read(iunit,*, err=100, end=100) t(j), am1(j), am2(j)
         end do
      CASE(3)
         do j = 1, lines
            read(iunit,*, err=100, end=100) t(j), am1(j), am2(j), am3(j)
         end do
      END SELECT
      call io_close(iunit)

      ! build splines
      !t = t*l%tau0 ! to scale for different units(?)
      call loct_spline_init(l%amplitude1)
      call loct_spline_fit(lines, t, am1, l%amplitude1)

      call loct_spline_init(l%amplitude2)
      call loct_spline_fit(lines, t, am2, l%amplitude2)

      call loct_spline_init(l%amplitude3)
      call loct_spline_fit(lines, t, am3, l%amplitude3)
      ! clean
      deallocate(t, am1, am2, am3)

    end subroutine get_envelope_from_file
  end subroutine laser_init


  ! ---------------------------------------------------------
  subroutine laser_end(no_l, l)
    type(laser_t), pointer :: l(:)
    integer,    intent(in) :: no_l

    integer :: i

    if(no_l <= 0) return

    do i = 1, no_l
      select case(l(i)%envelope)
      case(ENVELOPE_FROM_FILE)
        call loct_spline_end(l(i)%amplitude1)
        call loct_spline_end(l(i)%amplitude2)
        call loct_spline_end(l(i)%amplitude3)
      case(ENVELOPE_NUMERICAL)
        if(associated(l(i)%numerical)) then
          deallocate(l(i)%numerical); nullify(l(i)%numerical)
        end if
      end select
      if(l(i)%spatial_part == SPATIAL_PART_FROM_FILE) then
        deallocate(l(i)%v); nullify(l(i)%v)
      end if
    end do

    deallocate(l); nullify(l)

  end subroutine laser_end


  ! ---------------------------------------------------------
  subroutine laser_write_info(no_l, l, iunit)
    integer,       intent(in) :: no_l, iunit
    type(laser_t), intent(in) :: l(no_l)

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
        (l(i)%A0*CNST(5.14225e9))**2*CNST(1.3272e-3), " [W/cm^2]"
      if(l(i)%envelope == ENVELOPE_FROM_FILE) then
        write(iunit, '(3x,a)') 'Amplitudes are read from file.'
      else
        write(iunit,'(3x,a,f10.4,3a)') 'Width:     ', l(i)%tau0/units_inp%time%factor, &
          ' [', trim(units_inp%time%abbrev), ']'
        write(iunit,'(3x,a,f10.4,3a)') 'Middle t:  ', l(i)%t0/units_inp%time%factor, &
          ' [', trim(units_inp%time%abbrev), ']'
        if(l(i)%envelope == ENVELOPE_TRAPEZOIDAL) then
          write(iunit,'(3x,a,f10.4,3a)') 'Ramp time: ', l(i)%tau1/units_inp%time%factor, &
            ' [', trim(units_inp%time%abbrev), ']'
        end if
      end if
      if(l(i)%spatial_part == SPATIAL_PART_FROM_FILE) then
        write(iunit,'(3x,a)') 'The spatial part of the field is read from a file.'
      end if
    end do

  end subroutine laser_write_info


  ! ---------------------------------------------------------
  subroutine laser_potential(sb, no_l, l, t, m, pot)
    type(simul_box_t), intent(in) :: sb
    integer,           intent(in) :: no_l
    type(laser_t),     intent(in) :: l(no_l)
    FLOAT,             intent(in) :: t
    type(mesh_t),      intent(in) :: m
    FLOAT,            intent(out) :: pot(m%np)

    CMPLX :: amp
    integer :: i, j, jj
    FLOAT :: field(MAX_DIM)

    if(no_l == 0) then
      pot(1:m%np) = M_ZERO
      return
    end if

    pot = M_ZERO
    do i = 1, no_l
      ! FIXME: need to check if calculation of jj is correct (time-grid mapping)
      ! either M_HALF or M_ONE
      if(l(i)%envelope == ENVELOPE_NUMERICAL) then
         jj = int(abs(M_TWO*t/l(1)%dt) + M_HALF) ! steps of dt/2
         field(1:sb%dim) = l(i)%numerical(1:sb%dim, jj)
      else if(l(i)%envelope == ENVELOPE_FROM_FILE) then
         field(1) = loct_splint(l(i)%amplitude1, t)
         field(2) = loct_splint(l(i)%amplitude2, t)
         field(3) = loct_splint(l(i)%amplitude3, t)
      else
        call laser_amplitude(l(i), t, amp)
        field(1:sb%dim) = real(amp*l(i)%pol(1:sb%dim))
      endif

      if(l(i)%spatial_part == SPATIAL_PART_FROM_FILE) then
        do j = 1, m%np
          pot(j) = pot(j) + real(amp*l(i)%v(j))
        end do
      else

        do j = 1, m%np
          pot(j) = pot(j) + sum(field(1:sb%dim)*m%x(j,1:sb%dim))
        end do
      end if
    end do

  end subroutine laser_potential


  ! ---------------------------------------------------------
  subroutine laser_field(sb, no_l, l, t, field)
    type(simul_box_t), intent(in) :: sb
    integer,          intent(in)  :: no_l
    type(laser_t),    intent(in)  :: l(no_l)
    FLOAT,            intent(in)  :: t
    FLOAT,            intent(out) :: field(sb%dim)

    CMPLX :: amp
    integer :: i

    if(no_l == 0) then
      field(1:sb%dim) = M_ZERO
      return
    end if

   ! FIXME: need to check if calculation of i is correct (time-grid mapping)
   ! either M_HALF or M_ONE
   ! FIXME: need to extend this to the possibility of various laser fields.
    if(no_l == 1 .and. l(1)%envelope == ENVELOPE_NUMERICAL) then
       i = int(abs(M_TWO*t/l(1)%dt) + M_HALF) ! steps of dt/2
       field(1:sb%dim) = l(1)%numerical(1:sb%dim, i)
    else if(l(1)%envelope == ENVELOPE_FROM_FILE) then
       field(1) = loct_splint(l(1)%amplitude1, t)
       field(2) = loct_splint(l(1)%amplitude2, t)
       field(3) = loct_splint(l(1)%amplitude3, t)
    else
       field = M_ZERO
       do i = 1, no_l
          call laser_amplitude(l(i), t, amp)
          field(1:sb%dim) = field(1:sb%dim) + real(amp*l(i)%pol(1:sb%dim))
       end do
    end if
    
  end subroutine laser_field


  ! ---------------------------------------------------------
  !returns the laser vector field A
  subroutine laser_vector_field(sb, no_l, l, t, field)
    type(simul_box_t), intent(in) :: sb
    integer,          intent(in)  :: no_l
    type(laser_t),    intent(in)  :: l(no_l)
    FLOAT,            intent(in)  :: t
    FLOAT,            intent(out) :: field(sb%dim)

    CMPLX :: amp
    integer :: i

    field = M_ZERO
    do i = 1, no_l
      call laser_vector_amplitude(l(i), t, amp)
      field(1:sb%dim) = field(1:sb%dim) + real(amp*l(i)%pol(1:sb%dim))
    end do

  end subroutine laser_vector_field


  ! The following routines have to be changed in order to add more
  ! laser envelopes

  ! ---------------------------------------------------------
  ! returns the amplitude of the electric field
  subroutine laser_amplitude(l, t, amp)
    type(laser_t), intent(in) :: l
    FLOAT,         intent(in) :: t
    CMPLX,        intent(out) :: amp

    FLOAT :: r, ph

    ph = M_ZERO; r = M_ZERO
    select case (l%envelope)
    case(ENVELOPE_CONSTANT)
      r = M_ONE
    case(ENVELOPE_GAUSSIAN)
      r = exp(-(t - l%t0)**2 / (M_TWO*l%tau0**2))
    case(ENVELOPE_COSINOIDAL)
      if(abs(t - l%t0) <= l%tau0) then
        r = cos( (M_Pi/2)*((t - 2*l%tau0 - l%t0)/l%tau0) )
      end if
    case(ENVELOPE_TRAPEZOIDAL)
      if(t > l%t0-l%tau0/M_TWO-l%tau1 .and. t <= l%t0-l%tau0/M_TWO) then
        r = (t - (l%t0 - l%tau0/M_TWO - l%tau1))/l%tau1
      elseif(t>l%t0-l%tau0/M_TWO .and. t <=l%t0+l%tau0/M_TWO) then
        r = M_ONE
      elseif(t>l%t0+l%tau0/M_TWO .and. t <=l%t0+l%tau0/M_TWO+l%tau1) then
        r = (l%t0 + l%tau0/M_TWO + l%tau1 - t)/l%tau1
      end if
    case(ENVELOPE_FROM_FILE)
      r  = M_ONE  !loct_splint(l%amplitude, t)
      ph = M_ZERO !loct_splint(l%phase, t)
    end select

    amp = l%A0 * r * exp(M_zI * (l%omega0*t + ph))

    return
  end subroutine laser_amplitude


  ! ---------------------------------------------------------
  ! returns the amplitude of the vector field A
  subroutine laser_vector_amplitude(l, t, amp)
    type(laser_t), intent(in) :: l
    FLOAT,         intent(in) :: t
    CMPLX,        intent(out) :: amp

    CMPLX :: r, tt

    select case (l%envelope)
    case(ENVELOPE_GAUSSIAN) ! Gaussian
      ! not yet implemented
      r = M_ZERO
    case(ENVELOPE_COSINOIDAL) ! cosinoidal
      if(t < l%t0 - l%tau0) then
        tt = l%t0 - l%tau0
      else if(abs(t-l%t0) <= l%tau0) then
        tt = t
      else
        tt = l%t0 + l%tau0
      end if

      r = - l%r1 * (l%field_1 + exp(M_zI * l%omega0*tt) * &
        (M_zI*l%omega0*cos(l%a1 + l%w1*tt) + l%w1*sin(l%a1 + l%w1*tt)))
    case(ENVELOPE_TRAPEZOIDAL) ! ramp
      ! not yet implemented
      r = M_ZERO
    case default
      r = M_ZERO
    end select
    amp = l%A0 * r

    return
  end subroutine laser_vector_amplitude

end module lasers_m
