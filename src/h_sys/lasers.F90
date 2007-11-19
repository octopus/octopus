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
    laser_vector_potential,       &
    laser_to_numerical

  integer, parameter ::           &
    ENVELOPE_CONSTANT      =  0,  &
    ENVELOPE_GAUSSIAN      =  1,  &
    ENVELOPE_COSINOIDAL    =  2,  &
    ENVELOPE_TRAPEZOIDAL   =  3,  &
    ENVELOPE_FROM_FILE     = 10,  &
    ENVELOPE_FROM_EXPR     = 99
 
  integer, public, parameter ::           &
    E_FIELD_ELECTRIC         =  1,  &
    E_FIELD_MAGNETIC         =  2,  &
    E_FIELD_VECTOR_POTENTIAL =  3

  integer, public, parameter ::     &
    E_FIELD_SCALAR_POTENTIAL =  4

  type laser_t
    integer :: field
    CMPLX :: pol(MAX_DIM) ! the polarization of the laser
    type(tdf_t) :: f

    FLOAT, pointer :: v(:)
    FLOAT, pointer :: a(:, :)
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
    integer           :: i, j, ierr, envelope
    character(len=50) :: filename1
    character(len=200) :: scalar_pot_expression
    character(len=200) :: envelope_expression
    FLOAT :: a0, omega0, t0, tau0, tau1, r, pot_re, pot_im
    FLOAT, allocatable :: x(:)

    call push_sub('lasers.laser_init')

    call obsolete_variable("TDLasers", "TDExternalFields")
    !%Variable TDExternalFields
    !%Type block
    !%Section Time Dependent
    !%Description
    !% The block TDExternalFields describe the type and shape of time-dependent 
    !% external perturbationss that are applied to the system.
    !%
    !% Each line of the block describes an external field; this way you can actually have more
    !% than one laser (e.g. a "pump" and a "probe").
    !%
    !% The syntax of each line is, then:
    !%
    !% <tt>%TDExternalField
    !% <br>&nbsp;&nbsp; type | nx | ny | nz | envelope | amplitude | omega | tau0 | t0 | tau1 | filename1
    !% <br>%</tt>
    !%
    !% The first number of each line describes which kind of external field is described
    !% by the line: an electric field ("electric_field"), a magnetic field ("magnetic field"), 
    !% or a vector potential "vector_potential"). This last option, in the current version, is
    !% constant field in space, which permits to describe an electrical perturbation in the
    !% velocity gauge.
    !%
    !% The next three (possibly complex) numbers mark the polarization direction of the
    !% field. The "amplitude" is obviously the amplitude of the field. The "omega" is the
    !% frequency. The "envelope" decides the shape of the enveloping function.
    !% "tau0", "t0" and "tau1" are three paramenters that decide on the
    !% temporal shape of the pulse -- the exact details depend on the particular envelope.
    !%
    !% If the envelope is given in a file, this will be "filename1".
    !%
    !% The last two columns ("tau1", "filename1") are optional; they will
    !% only be searched if needed.
    !%
    !%Option envelope_constant 0
    !% The envelope is just the unit function. The laser is not a pulse, but a 
    !% continuous wave (cw).
    !%Option envelope_gaussian 1
    !% The envelope is a Gaussian function. To fully determine its shape, you need to 
    !% specify the width (tau0) and the peak time (t0).
    !%Option envelope_cosinusoidal 2
    !% The envelope is a half-cycle of a cosine function. To fully determine its shape, 
    !% you need to specify the width (tau0) and the peak time (t0).
    !%Option envelope_trapezoidal 3
    !% The envelope looks like a trapezoid: it starts at zero, ramps linearly until one during 
    !% tau1 units of time, stays at that maximum for tau0 units of time, and then decays 
    !% linearly to zero during tau1 units of time.
    !%Option envelope_fromfile 10
    !% The shape of the laser pulse is read from a file, indicated in the field "filename1".
    !%Option envelope_fromexpr 99
    !% The shape of the laser pulse is parsed from an expression.
    !%Option electric_field 1
    !% The external field is an electric field; the usual case when we want to describe a 
    !% laser in the length gauge.
    !%Option magnetic_field 2
    !% The external field is a (homogeneous) time-dependent magnetic field.
    !%Option vector_potential 3
    !% The external field is a time-depedent homogeneous vector potential, which may describe
    !% a laser field in the velocity gauge.
    !%Option scalar_potential 4
    !% The external field is an arbitrary scalar potential, which may describe an
    !% inhomogeneous electrical field
    !%End

    no_l = 0
    if(loct_parse_block(check_inp('TDExternalFields'), blk) == 0) then
      no_l = loct_parse_block_n(blk)
      ALLOCATE(l(no_l), no_l)

      ! The structure of the block is:
      ! nx | ny | nz | amplitude | omega | envelope | tau0 | t0 | tau1 | filename1 | filename2
      do i = 1, no_l

        a0 = M_ZERO; omega0 = M_ZERO; tau0 = M_ZERO; t0 = M_ZERO; tau1 = M_ZERO

        call loct_parse_block_int(blk, i-1, 0, l(i)%field)

        select case(l(i)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          call loct_parse_block_string(blk, i-1, 1, scalar_pot_expression)
          j = 1
          l(i)%pol = M_z1
        case default
          call loct_parse_block_cmplx(blk, i-1, 1, l(i)%pol(1))
          call loct_parse_block_cmplx(blk, i-1, 2, l(i)%pol(2))
          call loct_parse_block_cmplx(blk, i-1, 3, l(i)%pol(3))
          j = 3
        end select

        call loct_parse_block_int  (blk, i-1, j+1, envelope)

        select case(envelope)
        case(ENVELOPE_CONSTANT)
          call loct_parse_block_float(blk, i-1, j+2, a0)
          call loct_parse_block_float(blk, i-1, j+3, omega0)
        case(ENVELOPE_GAUSSIAN)
          call loct_parse_block_float(blk, i-1, j+2, a0)
          call loct_parse_block_float(blk, i-1, j+3, omega0)
          call loct_parse_block_float(blk, i-1, j+4, tau0)
          call loct_parse_block_float(blk, i-1, j+5, t0)
        case(ENVELOPE_COSINOIDAL)
          call loct_parse_block_float(blk, i-1, j+2, a0)
          call loct_parse_block_float(blk, i-1, j+3, omega0)
          call loct_parse_block_float(blk, i-1, j+4, tau0)
          call loct_parse_block_float(blk, i-1, j+5, t0)
        case(ENVELOPE_TRAPEZOIDAL)
          call loct_parse_block_float(blk, i-1, j+2, a0)
          call loct_parse_block_float(blk, i-1, j+3, omega0)
          call loct_parse_block_float(blk, i-1, j+4, tau0)
          call loct_parse_block_float(blk, i-1, j+5, t0)
          call loct_parse_block_float(blk, i-1, j+6, tau1)
        case(ENVELOPE_FROM_FILE)
          call loct_parse_block_string(blk, i-1, j+2, filename1)
        case(ENVELOPE_FROM_EXPR) ! This should be envelope_from_formula
          call loct_parse_block_string(blk, i-1, j+2, envelope_expression)
        case default
          call input_error('TDExternalFieldsd')
        end select

        a0     = a0 * units_inp%energy%factor / units_inp%length%factor
        omega0 = omega0 * units_inp%energy%factor
        tau0 = tau0 * units_inp%time%factor
        t0   = t0   * units_inp%time%factor
        tau1 = tau1 * units_inp%time%factor

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
        case(ENVELOPE_FROM_EXPR)
          call tdf_init_fromexpr(l(i)%f, trim(envelope_expression))
        end select

        l(i)%pol(:) = l(i)%pol(:)/sqrt(sum(abs(l(i)%pol(:))**2))

        select case(l(i)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          ALLOCATE(l(i)%v(m%np_part), m%np_part)
          l(i)%v = M_ZERO
          ALLOCATE(x(MAX_DIM), MAX_DIM)
          do j = 1, m%np
            call mesh_r(m, j, r, x = x)
            call loct_parse_expression(pot_re, pot_im, x(1), x(2), x(3), &
              r, M_ZERO, trim(scalar_pot_expression))
            l(i)%v(j) = pot_re
          end do
          deallocate(x)

        case(E_FIELD_MAGNETIC)
          ! WARNING: note that for the moment we are ignoring the possibility of a complex
          ! polarizability vector for the td magnetic field case.
          ALLOCATE(l(i)%a(m%np_part, m%sb%dim), m%np_part*m%sb%dim)
          l(i)%a = M_ZERO
          ALLOCATE(x(m%sb%dim), m%sb%dim)
          do j = 1, m%np
            x(1:m%sb%dim) = m%x(j, 1:m%sb%dim)
            select case(m%sb%dim)
            case(2)
              l(i)%a(j, :) = (/x(2), -x(1)/) * sign(CNST(1.0), real(l(i)%pol(3), REAL_PRECISION))
            case(3)
              l(i)%a(j, :) = (/ x(2)*l(i)%pol(3) - x(3)*l(i)%pol(2), &
                                x(3)*l(i)%pol(1) - x(1)*l(i)%pol(3), &
                                x(1)*l(i)%pol(2) - x(2)*l(i)%pol(1)  /)
            end select
          enddo
          l(i)%a = -M_HALF * l(i)%a 
          deallocate(x)

        end select

      end do

      call loct_parse_block_end(blk)
    end if

    call pop_sub()

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
        select case(l(i)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          deallocate(l(i)%v); nullify(l(i)%v)
        case(E_FIELD_MAGNETIC)
          deallocate(l(i)%a); nullify(l(i)%a)
        end select
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
      case(E_FIELD_SCALAR_POTENTIAL); write(iunit, '(a)') '   Scalar Potential.'
      end select
      if(l(i)%field .ne. E_FIELD_SCALAR_POTENTIAL) then
        write(iunit,'(3x,a,3(a1,f7.4,a1,f7.4,a1))') 'Polarization: ', &
          '(', real(l(i)%pol(1)), ',', aimag(l(i)%pol(1)), '), ', &
          '(', real(l(i)%pol(2)), ',', aimag(l(i)%pol(2)), '), ', &
          '(', real(l(i)%pol(3)), ',', aimag(l(i)%pol(3)), ')'
      end if
      call tdf_write(l(i)%f, iunit)
    end do

  end subroutine laser_write_info


  ! ---------------------------------------------------------
  subroutine laser_potential(sb, l, m, pot, t)
    type(simul_box_t), intent(in) :: sb
    type(laser_t),     intent(in) :: l
    type(mesh_t),      intent(in) :: m
    FLOAT,            intent(out) :: pot(:)
    FLOAT, optional,   intent(in) :: t

    CMPLX :: amp
    integer :: j
    FLOAT :: field(MAX_DIM)

    call push_sub('lasers.laser_potential')

    pot = M_ZERO
    if(present(t)) then
      amp = tdf(l%f, t)
    else
      amp = M_z1
    end if

    select case(l%field)
    case(E_FIELD_SCALAR_POTENTIAL)
      pot(1:m%np) = pot(1:m%np) + real(amp)*l%v(1:m%np)
    case default
      field(1:sb%dim) = real(amp*l%pol(1:sb%dim))
      do j = 1, m%np
        pot(j) = pot(j) + sum(field(1:sb%dim)*m%x(j, 1:sb%dim))
      end do
    end select

    call pop_sub()
  end subroutine laser_potential


  ! ---------------------------------------------------------
  subroutine laser_vector_potential(l, m, a, t)
    type(laser_t),     intent(in) :: l
    type(mesh_t),      intent(in) :: m
    FLOAT,            intent(out) :: a(:, :)
    FLOAT, optional,   intent(in) :: t

    call push_sub('lasers.laser_vector_potential')

    if(present(t)) then
      a = l%a * tdf(l%f, t)
    else
      a = l%a
    end if

    call pop_sub()
  end subroutine laser_vector_potential


  ! ---------------------------------------------------------
  subroutine laser_field(sb, l, field, t)
    type(simul_box_t), intent(in) :: sb
    type(laser_t),    intent(in)  :: l
    FLOAT,            intent(out) :: field(sb%dim)
    FLOAT, optional,  intent(in)  :: t


    CMPLX :: amp

    field = M_ZERO
    if(l%field .eq. E_FIELD_SCALAR_POTENTIAL) return
    if(present(t)) then
      amp = tdf(l%f, t)
    else
      amp = M_z1
    end if
    field(1:sb%dim) = field(1:sb%dim) + real(amp*l%pol(1:sb%dim))

  end subroutine laser_field



end module lasers_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
