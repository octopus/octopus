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
  use derivatives_m
  use global_m
  use mpi_m
  use grid_m
  use io_m
  use parser_m
  use mesh_m
  use messages_m
  use profiling_m
  use simul_box_m
  use states_dim_m
  use unit_m
  use unit_system_m
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
    laser_to_numerical,           &
    laser_to_numerical_all,       &
    laser_kind,                   &
    laser_polarization,           &
    laser_get_f,                  &
    laser_set_f,                  &
    laser_get_phi,                &
    laser_set_phi,                &
    laser_set_f_value,            &
    laser_carrier_frequency,      &
    zvlaser_operator_linear,      &
    zvlaser_operator_quadratic

  integer, public, parameter ::     &
    E_FIELD_NONE             =  0,  &
    E_FIELD_ELECTRIC         =  1,  &
    E_FIELD_MAGNETIC         =  2,  &
    E_FIELD_VECTOR_POTENTIAL =  3,  &
    E_FIELD_SCALAR_POTENTIAL =  4

  type laser_t
    private
    integer :: field      = E_FIELD_NONE  ! which kind of external field it is (electric, magnetic...)
    CMPLX :: pol(MAX_DIM) = M_z0          ! the polarization of the laser.
    type(tdf_t) :: f                      ! The envelope.
    type(tdf_t) :: phi                    ! The phase
    FLOAT :: omega        = M_ZERO        ! The main, "carrier", frequency.

    FLOAT, pointer :: v(:)    => NULL()
    FLOAT, pointer :: a(:, :) => NULL()
  end type laser_t

contains


  ! ---------------------------------------------------------
  FLOAT function laser_carrier_frequency(laser) result(w0)
    type(laser_t), intent(in) :: laser

    PUSH_SUB(laser_carrier_frequency)
    w0 = laser%omega

    POP_SUB(laser_carrier_frequency)
  end function laser_carrier_frequency
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function laser_kind(laser)
    type(laser_t), intent(in) :: laser

    ! no push_sub allowed in pure function
    laser_kind = laser%field

  end function laser_kind
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  function laser_polarization(laser) result(pol)
    type(laser_t), intent(in) :: laser
    CMPLX :: pol(MAX_DIM)

    PUSH_SUB(laser_polarization)
    pol(1:MAX_DIM) = laser%pol(1:MAX_DIM)

    POP_SUB(laser_polarization)
  end function laser_polarization
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_get_f(laser, ff)
    type(laser_t), intent(in)    :: laser
    type(tdf_t),   intent(inout) :: ff

    PUSH_SUB(laser_get_f)
    call tdf_copy(ff, laser%f)

    POP_SUB(laser_get_f)
  end subroutine laser_get_f
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_set_f(laser, ff)
    type(laser_t), intent(inout) :: laser
    type(tdf_t),   intent(inout) :: ff

    PUSH_SUB(laser_set_f)

    call tdf_end(laser%f)
    call tdf_copy(laser%f, ff)

    POP_SUB(laser_set_f)
  end subroutine laser_set_f
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_get_phi(laser, phi)
    type(laser_t), intent(in)    :: laser
    type(tdf_t),   intent(inout) :: phi

    PUSH_SUB(laser_get_phi)
    call tdf_copy(phi, laser%phi)

    POP_SUB(laser_get_phi)
  end subroutine laser_get_phi
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_set_phi(laser, phi)
    type(laser_t), intent(inout) :: laser
    type(tdf_t),   intent(inout) :: phi

    PUSH_SUB(laser_set_phi)

    call tdf_end(laser%phi)
    call tdf_copy(laser%phi, phi)

    POP_SUB(laser_set_phi)
  end subroutine laser_set_phi
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_set_f_value(laser, ii, xx)
    type(laser_t), intent(inout) :: laser
    integer,       intent(in)    :: ii
    FLOAT,         intent(in)    :: xx

    PUSH_SUB(laser_set_f_value)
    call tdf_set_numerical(laser%f, ii, xx)

    POP_SUB(laser_set_f_value)
  end subroutine laser_set_f_value
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! The td functions that describe the laser field are transformed to a 
  ! "numerical" representation (i.e. time grid, values at this time grid).
  ! The possible phase and carrier frequency are evaluated and put together with
  ! the envelope, so that the envelope describes the full function (zero phase,
  ! zero carrier frequency).
  ! ---------------------------------------------------------
  subroutine laser_to_numerical_all(laser, dt, max_iter, omegamax)
    type(laser_t), intent(inout)  :: laser
    FLOAT,         intent(in)     :: dt
    integer,       intent(in)     :: max_iter
    FLOAT,         intent(in)     :: omegamax

    integer :: iter
    FLOAT   :: tt, fj, phi

    PUSH_SUB(lasers_to_numerical_all)

    call tdf_to_numerical(laser%f, max_iter, dt, omegamax)
    do iter = 1, max_iter + 1
      tt = (iter-1)*dt
      fj = tdf(laser%f, iter)
      phi = tdf(laser%phi, tt)
      call tdf_set_numerical(laser%f, iter, fj*cos(laser%omega*tt+phi))
    end do
    call tdf_end(laser%phi)
    call tdf_init_cw(laser%phi, M_ZERO, M_ZERO)
    laser%omega = M_ZERO

    POP_SUB(lasers_to_numerical_all)
  end subroutine laser_to_numerical_all
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! The td functions that describe the laser field are transformed to a 
  ! "numerical" representation (i.e. time grid, values at this time grid).
  ! ---------------------------------------------------------
  subroutine laser_to_numerical(laser, dt, max_iter, omegamax)
    type(laser_t), intent(inout)  :: laser
    FLOAT,         intent(in)     :: dt
    integer,       intent(in)     :: max_iter
    FLOAT,         intent(in)     :: omegamax
    PUSH_SUB(lasers_to_numerical)

    call tdf_to_numerical(laser%f, max_iter, dt, omegamax)
    call tdf_to_numerical(laser%phi, max_iter, dt, omegamax)

    POP_SUB(lasers_to_numerical)
  end subroutine laser_to_numerical
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_init(no_l, lasers, mesh)
    integer,       intent(out) :: no_l
    type(laser_t), pointer     :: lasers(:)
    type(mesh_t),  intent(in)  :: mesh

    type(block_t)     :: blk
    integer           :: il, ip, jj, ierr
    character(len=200) :: scalar_pot_expression
    character(len=200) :: envelope_expression
    FLOAT :: omega0, rr, pot_re, pot_im
    FLOAT, allocatable :: xx(:)

    PUSH_SUB(laser_init)

    call messages_obsolete_variable("TDLasers", "TDExternalFields")
    !%Variable TDExternalFields
    !%Type block
    !%Section Time-Dependent
    !%Description
    !% The block <tt>TDExternalFields</tt> describes the type and shape of time-dependent 
    !% external perturbations that are applied to the system.
    !%
    !% Each line of the block describes an external field; this way you can actually have more
    !% than one laser (<i>e.g.</i> a "pump" and a "probe").
    !%
    !% The syntax of each line is:
    !%
    !% <tt>%TDExternalField
    !% <br>&nbsp;&nbsp; type | ...other descriptors...
    !% <br>%</tt>
    !%
    !% The first element of each line describes which kind of external field is described
    !% by the line: (i) an electric field (<tt>electric_field</tt>); (ii) a magnetic field 
    !% (<tt>magnetic_field</tt>); (iii) a vector potential (<tt>vector_potential</tt>) -- this option, 
    !% in the current version, is a field constant in space, which permits us to describe 
    !% an electric perturbation in the velocity gauge; (iv) an arbitrary scalar potential
    !% (<tt>scalar_potential</tt>).
    !%
    !% The "other descriptors" depend on which kind of external field has been indicated in 
    !% the first column.
    !%
    !% (A) type = <tt>electric field, magnetic field, vector_potential</tt>
    !%
    !% For these cases, the syntax is:
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; type | nx | ny | nz | omega | envelope_function_name
    !% <br>%</tt>
    !%
    !% The three (possibly complex) numbers (<i>nx</i>, <i>ny</i>, <i>nz</i>) mark the polarization 
    !% direction of the field. The float <tt>omega</tt> will be the carrier frequency of the
    !% pulse. The envelope of the field is a time-dependent function whose definition
    !% must be given in a <tt>TDFunctions</tt> block. <tt>envelope_function_name</tt> is a string (and therefore
    !% it must be surrounded by quotation marks) that must match one of the function names
    !% given in the first column of the <tt>TDFunctions</tt> block.
    !%
    !% (B) type = <tt>scalar_potential</tt>
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; scalar_potential | "scalar_expression" | freq | envelope_function_name
    !% <br>%</tt>
    !%
    !% The scalar potential is not just a dipole, but any expression given by the string
    !% "scalar_expression". The temporal shape is determined by the envelope function
    !% defined by <tt>envelope_function_name</tt>.
    !%
    !% A NOTE ON UNITS:
    !%
    !% It is very common to describe the strength of a laser field by its intensity, rather
    !% than using the electric-field amplitude. In atomic units (or, more precisely, in any
    !% Gaussian system of units), the relationship between instantaneous electric field
    !% and intensity is:
    !% <math> I(t) = \frac{c}{8\pi} E^2(t) </math>.
    !%
    !% It is common to read intensities in W/cm^2. The dimensions of intensities are
    !% [W]/(L^2T), where [W] are the dimensions of energy. The relevant conversion factors
    !% are:
    !%
    !% <math> Hartree / (a_0^2 atomic_time) = 6.4364086e+15 W / cm^2 </math>
    !% 
    !% <math> eV / ( angstrom^2 (hbar/eV) ) = 2.4341348e+12 W / cm^2 </math>
    !%
    !% If, in atomic units, we set the electric-field amplitude to <math>E_0</math>, 
    !% then the intensity is:
    !%
    !% <math> I_0 = 3.51 10^16 W/cm^2 (E_0^2) </math>
    !%
    !% If, working with "Units = ev_angstrom", we set <math>E_0</math>, then the intensity is:
    !%
    !% <math> I_0 = 1.327 10^13 (E_0^2) W/cm^2 </math>
    !%
    !%Option electric_field 1
    !% The external field is an electric field, the usual case when we want to describe a 
    !% laser in the length gauge.
    !%Option magnetic_field 2
    !% The external field is a (homogeneous) time-dependent magnetic field.
    !%Option vector_potential 3
    !% The external field is a time-dependent homogeneous vector potential, which may describe
    !% a laser field in the velocity gauge.
    !%Option scalar_potential 4
    !% The external field is an arbitrary scalar potential, which may describe an
    !% inhomogeneous electrical field.
    !%End

    no_l = 0
    if(parse_block(datasets_check('TDExternalFields'), blk) == 0) then
      no_l = parse_block_n(blk)
      SAFE_ALLOCATE(lasers(1:no_l))

      do il = 1, no_l

        call parse_block_integer(blk, il-1, 0, lasers(il)%field)

        select case(lasers(il)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          call parse_block_string(blk, il-1, 1, scalar_pot_expression)
          jj = 1
          lasers(il)%pol = M_z1
        case default
          call parse_block_cmplx(blk, il-1, 1, lasers(il)%pol(1))
          call parse_block_cmplx(blk, il-1, 2, lasers(il)%pol(2))
          call parse_block_cmplx(blk, il-1, 3, lasers(il)%pol(3))
          jj = 3
        end select

        call parse_block_float(blk, il-1, jj+1, omega0)
        omega0 = units_to_atomic(units_inp%energy, omega0)

        lasers(il)%omega = omega0
     
        call parse_block_string(blk, il-1, jj+2, envelope_expression)

        ! For some reason, one cannot open a block if another one is already open.
        ! This is why I close blk before calling tdf_read, and then open it again.
        ! This should be fixed at the parser level.
        call parse_block_end(blk)
        call tdf_read(lasers(il)%f, trim(envelope_expression), ierr)
        ierr = parse_block(datasets_check('TDExternalFields'), blk)

        ! Check if there is a phase.
        if(parse_block_cols(blk, il-1) > jj+3) then
          call parse_block_string(blk, il-1, jj+3, envelope_expression)
          call parse_block_end(blk)
          call tdf_read(lasers(il)%phi, trim(envelope_expression), ierr)
          ierr = parse_block(datasets_check('TDExternalFields'), blk)
        else
          call tdf_init(lasers(il)%phi)
        end if

        lasers(il)%pol(:) = lasers(il)%pol(:)/sqrt(sum(abs(lasers(il)%pol(:))**2))

        select case(lasers(il)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          SAFE_ALLOCATE(lasers(il)%v(1:mesh%np_part))
          lasers(il)%v = M_ZERO
          SAFE_ALLOCATE(xx(1:MAX_DIM))
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, coords = xx)
            call parse_expression(pot_re, pot_im, mesh%sb%dim, xx, rr, M_ZERO, trim(scalar_pot_expression))
            lasers(il)%v(ip) = pot_re
          end do
          SAFE_DEALLOCATE_A(xx)

        case(E_FIELD_MAGNETIC)
          ! \warning: note that for the moment we are ignoring the possibility of a complex
          ! polarizability vector for the td magnetic-field case.
          SAFE_ALLOCATE(lasers(il)%a(1:mesh%np_part, 1:mesh%sb%dim))
          lasers(il)%a = M_ZERO
          SAFE_ALLOCATE(xx(1:mesh%sb%dim))
          do ip = 1, mesh%np
            xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
            select case(mesh%sb%dim)
            case(2)
              lasers(il)%a(ip, :) = (/xx(2), -xx(1)/) * sign(CNST(1.0), real(lasers(il)%pol(3), REAL_PRECISION))
            case(3)
              lasers(il)%a(ip, :) = (/ xx(2)*lasers(il)%pol(3) - xx(3)*lasers(il)%pol(2), &
                                xx(3)*lasers(il)%pol(1) - xx(1)*lasers(il)%pol(3), &
                                xx(1)*lasers(il)%pol(2) - xx(2)*lasers(il)%pol(1)  /)
            end select
          enddo
          lasers(il)%a = -M_HALF * lasers(il)%a 
          SAFE_DEALLOCATE_A(xx)

        end select

      end do

      call parse_block_end(blk)
    end if

    POP_SUB(laser_init)

  end subroutine laser_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_end(no_l, lasers)
    type(laser_t), pointer :: lasers(:)
    integer,    intent(in) :: no_l

    integer :: il

    PUSH_SUB(laser_end)

    if(no_l > 0) then
      do il = 1, no_l
        call tdf_end(lasers(il)%f)
        call tdf_end(lasers(il)%phi)
        select case(lasers(il)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          SAFE_DEALLOCATE_P(lasers(il)%v)
        case(E_FIELD_MAGNETIC)
          SAFE_DEALLOCATE_P(lasers(il)%a)
        end select
      end do
      SAFE_DEALLOCATE_P(lasers)
    end if

    POP_SUB(laser_end)
  end subroutine laser_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_write_info(lasers, iunit, dt, max_iter)
    type(laser_t),     intent(in) :: lasers(:)
    integer,           intent(in) :: iunit
    FLOAT,   optional, intent(in) :: dt 
    integer, optional, intent(in) :: max_iter 

    FLOAT :: tt, fluence, max_intensity, intensity, dt_
    CMPLX :: amp, val
    integer :: il, iter, idir, no_l, max_iter_

    if(.not.mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(laser_write_info)

    no_l = size(lasers)
    
    do il = 1, no_l

      if(present(dt)) then
        dt_ = dt
      else
        dt_ = tdf_dt(lasers(il)%f)
      end if
      if(present(max_iter)) then
        max_iter_ = max_iter
      else
        max_iter_ = tdf_niter(lasers(il)%f)
      end if

      write(iunit,'(i2,a)') il, ':'
      select case(lasers(il)%field)
      case(E_FIELD_ELECTRIC);   write(iunit, '(a)') '   Electric Field.'
      case(E_FIELD_MAGNETIC);   write(iunit, '(a)') '   Magnetic Field.'
      case(E_FIELD_VECTOR_POTENTIAL); write(iunit, '(a)') '   Vector Potential.'
      case(E_FIELD_SCALAR_POTENTIAL); write(iunit, '(a)') '   Scalar Potential.'
      end select
      if(lasers(il)%field .ne. E_FIELD_SCALAR_POTENTIAL) then
        write(iunit,'(3x,a,3(a1,f7.4,a1,f7.4,a1))') 'Polarization: ', &
          '(', real(lasers(il)%pol(1)), ',', aimag(lasers(il)%pol(1)), '), ', &
          '(', real(lasers(il)%pol(2)), ',', aimag(lasers(il)%pol(2)), '), ', &
          '(', real(lasers(il)%pol(3)), ',', aimag(lasers(il)%pol(3)), ')'
      end if
      write(iunit,'(3x,a,f14.8,3a)') 'Carrier frequency = ', &
        units_from_atomic(units_out%energy, lasers(il)%omega), &
        ' [', trim(units_abbrev(units_out%energy)), ']'
      write(iunit,'(3x,a)')       'Envelope: ' 
      call tdf_write(lasers(il)%f, iunit)

      if(.not.tdf_is_empty(lasers(il)%phi)) then
        write(iunit,'(3x,a)')       'Phase: ' 
        call tdf_write(lasers(il)%phi, iunit)
      end if

      ! 1 atomic unit of intensity = 3.5094448e+16 W / cm^2
      ! In a Gaussian system of units,
      ! I(t) = (1/(8\pi)) * c * E(t)^2
      ! (1/(8\pi)) * c = 5.4525289841210 a.u.
      if(lasers(il)%field .eq. E_FIELD_ELECTRIC) then
        fluence = M_ZERO

        max_intensity = M_ZERO
        do iter = 1, max_iter_
          tt = iter * dt_
          val = tdf(lasers(il)%f, tt)
          amp = val*exp(M_zI*(lasers(il)%omega*tt + tdf(lasers(il)%phi, tt)))
          intensity = M_ZERO
          do idir = 1, MAX_DIM
            intensity = intensity + CNST(5.4525289841210) * real(amp*lasers(il)%pol(idir))**2
          end do
          fluence = fluence + intensity
          if(intensity > max_intensity) max_intensity = intensity
        end do
        fluence = fluence * dt_
        write(iunit,'(a,es12.6,3a)') '   Peak intensity = ', max_intensity, ' [a.u]'
        write(iunit,'(a,es12.6,3a)') '                  = ', &
          max_intensity * 6.4364086e+15, ' [W/cm^2]'
        write(iunit,'(a,es12.6,a)')  '   Int. intensity = ', fluence, ' [a.u]'
        write(iunit,'(a,es12.6,a)')  '   Fluence        = ', &
          fluence / CNST(5.4525289841210) , ' [a.u]'
      end if

    end do

    POP_SUB(laser_write_info)
  end subroutine laser_write_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_potential(laser, mesh, pot, time)
    type(laser_t),     intent(in)    :: laser
    type(mesh_t),      intent(in)    :: mesh
    FLOAT,             intent(inout) :: pot(:)
    FLOAT, optional,   intent(in)    :: time

    CMPLX :: amp
    integer :: ip
    FLOAT :: field(MAX_DIM)

    PUSH_SUB(laser_potential)

    if(present(time)) then
      amp = tdf(laser%f, time) * exp(M_zI * ( laser%omega*time + tdf(laser%phi, time) ) )
    else
      amp = M_z1
    end if

    select case(laser%field)
    case(E_FIELD_SCALAR_POTENTIAL)
      pot(1:mesh%np) = pot(1:mesh%np) + real(amp)*laser%v(1:mesh%np)
    case default
      field(1:mesh%sb%dim) = real(amp*laser%pol(1:mesh%sb%dim))
      do ip = 1, mesh%np
        pot(ip) = pot(ip) + sum(field(1:mesh%sb%dim)*mesh%x(ip, 1:mesh%sb%dim))
      end do
    end select

    POP_SUB(laser_potential)
  end subroutine laser_potential
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_vector_potential(laser, mesh, aa, time)
    type(laser_t),   intent(in)    :: laser
    type(mesh_t),    intent(in)    :: mesh
    FLOAT,           intent(inout) :: aa(:, :)
    FLOAT, optional, intent(in)    :: time
    
    FLOAT   :: amp 
    integer :: ip, idir

    PUSH_SUB(laser_vector_potential)

    if(present(time)) then
      amp = real(tdf(laser%f, time)*exp(M_zI*(laser%omega*time + tdf(laser%phi, time))))
      forall(idir = 1:mesh%sb%dim, ip = 1:mesh%np) aa(ip, idir) = aa(ip, idir) + amp*laser%a(ip, idir)
    else
      forall(idir = 1:mesh%sb%dim, ip = 1:mesh%np) aa(ip, idir) = aa(ip, idir) + laser%a(ip, idir)
    end if

    POP_SUB(laser_vector_potential)
  end subroutine laser_vector_potential
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_field(laser, sb, field, time)
    type(laser_t),     intent(in)    :: laser
    type(simul_box_t), intent(in)    :: sb
    FLOAT,             intent(inout) :: field(:)
    FLOAT, optional,   intent(in)    :: time

    CMPLX :: amp

    !no PUSH SUB, called too often

    if(laser%field .eq. E_FIELD_SCALAR_POTENTIAL) then
      return
    endif
    if(present(time)) then
      amp = tdf(laser%f, time) * exp(M_zI * ( laser%omega * time + tdf(laser%phi, time) ) )
    else
      amp = M_z1
    end if
    field(1:sb%dim) = field(1:sb%dim) + real(amp*laser%pol(1:sb%dim))

  end subroutine laser_field

#include "undef.F90"
#include "real.F90"
#include "lasers_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lasers_inc.F90"

end module lasers_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
