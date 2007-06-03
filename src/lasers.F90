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
!! $Id$

#include "global.h"

module lasers_m
  use datasets_m
  use global_m
  use io_m
  use lib_oct_gsl_spline_m
  use lib_oct_parser_m
  use mesh_m
  use messages_m
  use output_m
  use simul_box_m
  use units_m
  use tdf_m

  implicit none

  private
  public ::                       &
    laser_t,                      &
    laser_init,                   &
    laser_end,                    &
    laser_write_info,             &
    laser_field,                  &
    laser_potential,              &
    laser_to_numerical

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

  integer, public, parameter ::           &
    E_FIELD_ELECTRIC         =  1,  &
    E_FIELD_MAGNETIC         =  2,  &
    E_FIELD_VECTOR_POTENTIAL =  3

  type laser_t
    integer :: field
    CMPLX :: pol(MAX_DIM) ! the polarization of the laser
    type(tdf_t) :: f
    integer  :: spatial_part
    FLOAT, pointer :: v(:) ! the spatial part.
  end type laser_t

contains


  ! ---------------------------------------------------------
  subroutine laser_to_numerical(l, dt, max_iter)
    type(laser_t) :: l
    FLOAT, intent(in) :: dt
    integer, intent(in) :: max_iter

    call push_sub('lasers.lasers_to_numerical')

    call tdf_to_numerical(l%f, max_iter, dt)

    call pop_sub()
  end subroutine laser_to_numerical


  ! ---------------------------------------------------------
  subroutine laser_init(no_l, l, m)
    integer,     intent(out) :: no_l
    type(laser_t),   pointer :: l(:)
    type(mesh_t), intent(in) :: m

    C_POINTER         :: blk
    integer           :: i, no_c, ierr
    character(len=50) :: filename1, filename2

    integer :: envelope
    FLOAT :: a0, omega0, t0, tau0, tau1

    call push_sub('lasers.laser_init')

    call obsolete_variable("TDLasers", "TDExternalFields")
    !%Variable TDExternalFields
    !%Type block
    !%Section Time Dependent
    !%Description
    !% The block TDExternalFields describe the type and shape of time-dependent external perturbations
    !% that are applied to the system.
    !%
    !% Each line of the block describes an external field; this way you can actually have more
    !% than one laser (e.g. a "pump" and a "probe").
    !%
    !% The syntax of each line is, then:
    !%
    !% <tt>%TDLasers
    !% <br>&nbsp;&nbsp; type | nx | ny | nz | amplitude | omega | envelope | tau0 | t0 | tau1 | filename1 | filename2
    !% <br>%</tt>
    !%
    !% The first number of each line describes which kind of external field is described by the line:
    !% an electric field ("electric_field"), a magnetic field ("magnetic field"), or a
    !% vector potential "vector_potential"). This last option, in the current version, is
    !% constant field in space, which permits to describe an electrical perturbation in the
    !% velocity gauge.
    !%
    !% The next three (possibly complex) numbers mark the polarization direction of the
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
    !%Option envelope_constant 0
    !% The envelope is just the unit function. The laser is not a pulse, but a continuous wave (cw).
    !%Option envelope_gaussian 1
    !% The envelope is a Gaussian function. To fully determine its shape, you need to 
    !% specify the width (tau0) and the peak time (t0).
    !%Option envelope_cosinusoidal 2
    !% The envelope is a half-cycle of a cosine function. To fully determine its shape, you need to 
    !% specify the width (tau0)
    !% and the peak time (t0).
    !%Option envelope_trapezoidal 3
    !% The envelope looks like a trapezoid: it starts at zero, ramps linearly until one during 
    !% tau1 units of time, stays at that maximum for tau0 units of time, and then decays 
    !% linearly to zero during tau1 units of time.
    !%Option envelope_fromfile 10
    !% The shape of the laser pulse is read from a file, indicated in the field "filename1".
    !%Option envelope_numerical 99
    !% The laser shape is generated internally. This option is used if the code is run in the 
    !% optimal-control mode.
    !%Option electric_field 1
    !% The external field is an electric field; the usual case when we want to describe a 
    !% laser in the length gauge.
    !%Option magnetic_field 2
    !% The external field is a (homogeneous) time-dependent magnetic field.
    !%Option vector_potential 3
    !% The external field is a time-depedent homogeneous vector potential, which may describe
    !% a laser field in the velocity gauge.
    !%End

    no_l = 0
    if(loct_parse_block(check_inp('TDExternalFields'), blk) == 0) then
      no_l = loct_parse_block_n(blk)
      ALLOCATE(l(no_l), no_l)

      ! The structure of the block is:
      ! nx | ny | nz | amplitude | omega | envelope | tau0 | t0 | tau1 | filename1 | filename2
      do i = 1, no_l
        no_c = loct_parse_block_cols(blk, i-1)

        call loct_parse_block_int(blk, i-1, 0, l(i)%field)
        call loct_parse_block_cmplx(blk, i-1, 1, l(i)%pol(1))
        call loct_parse_block_cmplx(blk, i-1, 2, l(i)%pol(2))
        call loct_parse_block_cmplx(blk, i-1, 3, l(i)%pol(3))
        call loct_parse_block_float(blk, i-1, 4, a0)
        call loct_parse_block_float(blk, i-1, 5, omega0)
        call loct_parse_block_int  (blk, i-1, 6, envelope)
        call loct_parse_block_float(blk, i-1, 7, tau0)
        call loct_parse_block_float(blk, i-1, 8, t0)
        if(no_c>9)  then
          call loct_parse_block_float(blk, i-1, 9, tau1)
        else
          tau1 = M_ZERO
        end if
        if(no_c>10) call loct_parse_block_string(blk, i-1, 10, filename1)
        if(no_c>11) call loct_parse_block_string(blk, i-1, 11, filename2)

        a0     = a0 * units_inp%energy%factor / units_inp%length%factor
        omega0 = omega0 * units_inp%energy%factor
        tau0 = tau0 * units_inp%time%factor
        t0   = t0   * units_inp%time%factor
        tau1 = tau1 * units_inp%time%factor


        if(envelope == ENVELOPE_TRAPEZOIDAL .and. no_c < 9)  &
             call input_error('TDExternalFields')
        if(envelope == ENVELOPE_FROM_FILE   .and. no_c < 10) &
             call input_error('TDExternalFields')

        select case(envelope)
        case(ENVELOPE_CONSTANT)
          call tdf_init_cw(l(i)%f, a0, omega0)
        case(ENVELOPE_GAUSSIAN)
          call tdf_init_gaussian(l(i)%f, a0, omega0, t0, tau0)
        case(ENVELOPE_COSINOIDAL)
          call tdf_init_cosinoidal(l(i)%f, a0, omega0, t0, tau0)
        case(ENVELOPE_TRAPEZOIDAL)
          call tdf_init_trapezoidal(l(i)%f, a0, omega0, t0, tau0, tau1)
        case(ENVELOPE_FROM_FILE)
          call tdf_init_fromfile(l(i)%f, trim(filename1), ierr)
        case(ENVELOPE_NUMERICAL)
          write(message(1),'(a)') 'Internal Error at laser_init'
          call write_fatal(1)
        end select

        if(no_c>10) l(i)%spatial_part = SPATIAL_PART_FROM_FILE

        if(l(i)%spatial_part == SPATIAL_PART_FROM_FILE) &
             call get_spatial_part_from_file(l(i), trim(filename2))

        ! normalize polarization
        l(i)%pol(:) = l(i)%pol(:)/sqrt(sum(abs(l(i)%pol(:))**2))

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


  end subroutine laser_init


  ! ---------------------------------------------------------
  subroutine laser_end(no_l, l)
    type(laser_t), pointer :: l(:)
    integer,    intent(in) :: no_l

    integer :: i

    call push_sub('laser.laser_end')

    if(no_l > 0) then
      do i = 1, no_l
        call tdf_end(l(i)%f)
        if(l(i)%spatial_part == SPATIAL_PART_FROM_FILE) then
          deallocate(l(i)%v); nullify(l(i)%v)
        end if
      end do
      deallocate(l); nullify(l)
    end if

    call pop_sub()
  end subroutine laser_end


  ! ---------------------------------------------------------
  subroutine laser_write_info(no_l, l, iunit)
    integer,       intent(in) :: no_l, iunit
    type(laser_t), intent(in) :: l(no_l)

    integer :: i

    do i = 1, no_l
      write(iunit,'(i2,a)') i,':'
      select case(l(i)%field)
      case(E_FIELD_ELECTRIC);   write(iunit, '(a)') '   Electric Field.'
      case(E_FIELD_MAGNETIC);   write(iunit, '(a)') '   Magnetic Field.'
      case(E_FIELD_VECTOR_POTENTIAL); write(iunit, '(a)') '   Potential vector.'
      end select
      write(iunit,'(3x,a,3(a1,f7.4,a1,f7.4,a1))') 'Polarization: ', &
        '(', real(l(i)%pol(1)), ',', aimag(l(i)%pol(1)), '), ', &
        '(', real(l(i)%pol(2)), ',', aimag(l(i)%pol(2)), '), ', &
        '(', real(l(i)%pol(3)), ',', aimag(l(i)%pol(3)), ')'
      call tdf_write(l(i)%f, iunit)
      if(l(i)%spatial_part == SPATIAL_PART_FROM_FILE) then
        write(iunit,'(3x,a)') 'The spatial part of the field is read from a file.'
      end if
    end do

  end subroutine laser_write_info


  ! ---------------------------------------------------------
  subroutine laser_potential(sb, l, t, m, pot)
    type(simul_box_t), intent(in) :: sb
    type(laser_t),     intent(in) :: l
    FLOAT,             intent(in) :: t
    type(mesh_t),      intent(in) :: m
    FLOAT,            intent(out) :: pot(:)

    CMPLX :: amp
    integer :: j, jj
    FLOAT :: field(MAX_DIM)

    call push_sub('lasers.laser_potential')


    pot = M_ZERO
    amp = tdf(l%f, t)
    field(1:sb%dim) = real(amp*l%pol(1:sb%dim))

    if(l%spatial_part == SPATIAL_PART_FROM_FILE) then
      do j = 1, m%np
        pot(j) = pot(j) + real(amp*l%v(j))
      end do
    else
      do j = 1, m%np
        pot(j) = pot(j) + sum(field(1:sb%dim)*m%x(j, 1:sb%dim))
      end do
    end if

    call pop_sub()
  end subroutine laser_potential


  ! ---------------------------------------------------------
  subroutine laser_field(sb, l, t, field)
    type(simul_box_t), intent(in) :: sb
    type(laser_t),    intent(in)  :: l
    FLOAT,            intent(in)  :: t
    FLOAT,            intent(out) :: field(sb%dim)

    CMPLX :: amp
    integer :: i

    field = M_ZERO
    amp = tdf(l%f, t)
    field(1:sb%dim) = field(1:sb%dim) + real(amp*l%pol(1:sb%dim))

  end subroutine laser_field



end module lasers_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
