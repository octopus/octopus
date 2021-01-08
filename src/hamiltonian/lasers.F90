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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module lasers_oct_m
  use derivatives_oct_m
  use global_oct_m
  use mpi_oct_m
  use parser_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use states_elec_dim_oct_m
  use symmetries_oct_m
  use symm_op_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use tdfunction_oct_m

  implicit none

  private
  public ::                       &
    laser_t,                      &
    laser_init,                   &
    laser_end,                    &
    laser_write_info,             &
    laser_field,                  &
    laser_electric_field,         & 
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
    laser_set_empty_phi,          &
    laser_set_f_value,            &
    laser_carrier_frequency,      &
    laser_set_frequency,          &
    laser_set_polarization,       &
    vlaser_operator_linear,       &
    vlaser_operator_quadratic


  integer, public, parameter ::     &
    E_FIELD_NONE             =  0,  &
    E_FIELD_ELECTRIC         =  1,  &
    E_FIELD_MAGNETIC         =  2,  &
    E_FIELD_VECTOR_POTENTIAL =  3,  &
    E_FIELD_SCALAR_POTENTIAL =  4

  type laser_t
    private
    integer :: field      = E_FIELD_NONE  !< which kind of external field it is (electric, magnetic...)
    CMPLX :: pol(MAX_DIM) = M_z0          !< the polarization of the laser.
    type(tdf_t) :: f                      !< The envelope.
    type(tdf_t) :: phi                    !< The phase
    FLOAT :: omega        = M_ZERO        !< The main, "carrier", frequency.

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
  integer pure elemental function laser_kind(laser)
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
  subroutine laser_set_empty_phi(laser)
    type(laser_t),          intent(inout) :: laser

    PUSH_SUB(laser_set_empty_phi)

    call tdf_init(laser%phi)

    POP_SUB(laser_set_empty_phi)
  end subroutine laser_set_empty_phi
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
  subroutine laser_set_frequency(laser, omega)
    type(laser_t), intent(inout) :: laser
    FLOAT,         intent(in)    :: omega

    PUSH_SUB(laser_carrier_frequency)
    laser%omega = omega

    POP_SUB(laser_set_frequency)
  end subroutine laser_set_frequency
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_set_polarization(laser, pol)
    type(laser_t), intent(inout) :: laser
    CMPLX,         intent(in)    :: pol(MAX_DIM)

    PUSH_SUB(laser_set_polarization)
    laser%pol(1:MAX_DIM) = pol(1:MAX_DIM)

    POP_SUB(laser_set_polarization)
  end subroutine laser_set_polarization
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> The td functions that describe the laser field are transformed to a 
  !! "numerical" representation (i.e. time grid, values at this time grid).
  !!
  !! The possible phase and carrier frequency are evaluated and put together with
  !! the envelope, so that the envelope describes the full function (zero phase,
  !! zero carrier frequency).
  ! ---------------------------------------------------------
  subroutine laser_to_numerical_all(laser, dt, max_iter, omegamax)
    type(laser_t),   intent(inout)  :: laser
    FLOAT,           intent(in)     :: dt
    integer,         intent(in)     :: max_iter
    FLOAT,           intent(in)     :: omegamax

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
  !> The td functions that describe the laser field are transformed to a 
  !! "numerical" representation (i.e. time grid, values at this time grid).
  ! ---------------------------------------------------------
  subroutine laser_to_numerical(laser, dt, max_iter, omegamax)
    type(laser_t),   intent(inout) :: laser
    FLOAT,           intent(in)    :: dt
    integer,         intent(in)    :: max_iter
    FLOAT,           intent(in)    :: omegamax
    
    PUSH_SUB(lasers_to_numerical)

    call tdf_to_numerical(laser%f, max_iter, dt, omegamax)
    call tdf_to_numerical(laser%phi,  max_iter, dt, omegamax)

    POP_SUB(lasers_to_numerical)
  end subroutine laser_to_numerical
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_init(lasers, namespace, no_l, mesh)
    type(laser_t),       pointer     :: lasers(:)
    type(namespace_t),   intent(in)  :: namespace
    integer,             intent(out) :: no_l
    type(mesh_t),        intent(in)  :: mesh

    type(block_t)     :: blk
    integer           :: il, ip, jj, ierr
    character(len=200) :: scalar_pot_expression
    character(len=200) :: envelope_expression, phase_expression
    FLOAT :: omega0, rr, pot_re, pot_im, xx(MAX_DIM)
    integer :: iop 

    PUSH_SUB(laser_init)

    call messages_obsolete_variable(namespace, "TDLasers", "TDExternalFields")
    !%Variable TDExternalFields
    !%Type block
    !%Section Time-Dependent
    !%Description
    !% The block <tt>TDExternalFields</tt> describes the type and shape of time-dependent 
    !% external perturbations that are applied to the system, in the form
    !% <math>f(x,y,z) \cos(\omega t + \phi (t)) g(t)</math>, where <math>f(x,y,z)</math> is defined by
    !% by a field type and polarization or a scalar potential, as below; <math>\omega</math>
    !% is defined by <tt>omega</tt>; <math>g(t)</math> is defined by
    !% <tt>envelope_function_name</tt>; and <math>\phi(t)</math> is the (time-dependent) phase from <tt>phase</tt>.
    !%
    !% These perturbations are only applied for time-dependent runs. If
    !% you want the value of the perturbation at time zero to be
    !% applied for time-independent runs, use <tt>TimeZero = yes</tt>.
    !%
    !% Each line of the block describes an external field; this way you can actually have more
    !% than one laser (<i>e.g.</i> a "pump" and a "probe").
    !%
    !% There are two ways to specify <math>f(x,y,z)</math> but both use the same <tt>omega | envelope_function_name [| phase]</tt>
    !% for the time-dependence.
    !% The float <tt>omega</tt> will be the carrier frequency of the
    !% pulse (in energy units). The envelope of the field is a time-dependent function whose definition
    !% must be given in a <tt>TDFunctions</tt> block. <tt>envelope_function_name</tt> is a string (and therefore
    !% it must be surrounded by quotation marks) that must match one of the function names
    !% given in the first column of the <tt>TDFunctions</tt> block.
    !% <tt>phase</tt> is optional and is taken to be zero if not provided, and is also a string specifying
    !% a time-dependent function.
    !%
    !% (A) type = <tt>electric field, magnetic field, vector_potential</tt>
    !%
    !% For these cases, the syntax is:
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; type | nx | ny | nz | omega | envelope_function_name | phase
    !% <br>%</tt>
    !%
    !% The <tt>vector_potential</tt> option (constant in space) permits us to describe
    !% an electric perturbation in the velocity gauge.
    !% The three (possibly complex) numbers (<tt>nx</tt>, <tt>ny</tt>, <tt>nz</tt>) mark the polarization
    !% direction of the field.
    !%
    !% (B) type = <tt>scalar_potential</tt>
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; scalar_potential | "spatial_expression" | omega | envelope_function_name | phase
    !% <br>%</tt>
    !%
    !% The scalar potential is any expression of the spatial coordinates given by the string
    !% "spatial_expression", allowing a field beyond the dipole approximation.
    !%
    !% For DFTB runs, only fields of type type = <tt>electric field</tt> are allowed for the moment, and the
    !% <tt>type</tt> keyword is omitted.
    !%
    !% A NOTE ON UNITS:
    !%
    !% It is very common to describe the strength of a laser field by its intensity, rather
    !% than using the electric-field amplitude. In atomic units (or, more precisely, in any
    !% Gaussian system of units), the relationship between instantaneous electric field
    !% and intensity is:
    !% <math> I(t) = \frac{c}{8\pi} E^2(t) </math>.
    !%
    !% It is common to read intensities in W/cm<math>^2</math>. The dimensions of intensities are
    !% [W]/(L<math>^2</math>T), where [W] are the dimensions of energy. The relevant conversion factors
    !% are:
    !%
    !% Hartree / (<math>a_0^2</math> atomic_time) = <math>6.4364086 \times 10^{15} \mathrm{W/cm}^2</math>
    !% 
    !% eV / ( &Aring;<math>^2 (\hbar</math>/eV) ) = <math>2.4341348 \times 10^{12} \mathrm{W/cm}^2</math>
    !%
    !% If, in atomic units, we set the electric-field amplitude to <math>E_0</math>,
    !% then the intensity is:
    !%
    !% <math> I_0 = 3.51 \times 10^{16} \mathrm{W/cm}^2 (E_0^2) </math>
    !%
    !% If, working with <tt>Units = ev_angstrom</tt>, we set <math>E_0</math>, then the intensity is:
    !%
    !% <math> I_0 = 1.327 \times 10^{13} (E_0^2) \mathrm{W/cm}^2 </math>
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
    if(parse_block(namespace, 'TDExternalFields', blk) == 0) then
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
        call tdf_read(lasers(il)%f, namespace, trim(envelope_expression), ierr)

        ! Check if there is a phase.
        if(parse_block_cols(blk, il-1) > jj+3) then
          call parse_block_string(blk, il-1, jj+3, phase_expression)
          call tdf_read(lasers(il)%phi, namespace, trim(phase_expression), ierr)
          if (ierr /= 0) then            
            write(message(1),'(3A)') 'Error in the "', trim(envelope_expression), '" field defined in the TDExternalFields block:'
            write(message(2),'(3A)') 'Time-dependent phase function "', trim(phase_expression), '" not found.'
            call messages_warning(2, namespace=namespace)
          end if
        else
          call tdf_init(lasers(il)%phi)
        end if

        lasers(il)%pol(:) = lasers(il)%pol(:)/sqrt(sum(abs(lasers(il)%pol(:))**2))

        select case(lasers(il)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          SAFE_ALLOCATE(lasers(il)%v(1:mesh%np_part))
          lasers(il)%v = M_ZERO
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, coords = xx)
            call parse_expression(pot_re, pot_im, mesh%sb%dim, xx, rr, M_ZERO, trim(scalar_pot_expression))
            lasers(il)%v(ip) = pot_re
          end do

        case(E_FIELD_MAGNETIC)
          ! \warning: note that for the moment we are ignoring the possibility of a complex
          ! polarizability vector for the td magnetic-field case.
          SAFE_ALLOCATE(lasers(il)%a(1:mesh%np_part, 1:mesh%sb%dim))
          lasers(il)%a = M_ZERO
          do ip = 1, mesh%np
            xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
            select case(mesh%sb%dim)
            case(2)
              lasers(il)%a(ip, :) = (/xx(2), -xx(1)/) * sign(CNST(1.0), real(lasers(il)%pol(3)))
            case(3)
              lasers(il)%a(ip, :) = (/ xx(2)*real(lasers(il)%pol(3)) - xx(3)*real(lasers(il)%pol(2)), &
                                xx(3)*real(lasers(il)%pol(1)) - xx(1)*real(lasers(il)%pol(3)), &
                                xx(1)*real(lasers(il)%pol(2)) - xx(2)*real(lasers(il)%pol(1))  /)
            end select
          end do
          lasers(il)%a = -M_HALF * lasers(il)%a 

        end select

      end do

      call parse_block_end(blk)
    end if

    if(mesh%sb%kpoints%use_symmetries) then
      do iop = 1, symmetries_number(mesh%sb%symm)
        if(iop == symmetries_identity_index(mesh%sb%symm)) cycle
        do il = 1, no_l
          if(.not. symm_op_invariant_cart(mesh%sb%symm%ops(iop), lasers(il)%pol(:), SYMPREC)) then
            message(1) = "The lasers break (at least) one of the symmetries used to reduce the k-points."
            message(2) = "Set SymmetryBreakDir accordingly to your laser fields."
            call messages_fatal(2, namespace=namespace)
          end if
        end do
      end do
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

    FLOAT :: tt, fluence, max_intensity, intensity, dt_, field(MAX_DIM), Up, maxfield,tmp
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
      if(lasers(il)%field /= E_FIELD_SCALAR_POTENTIAL) then
        write(iunit,'(3x,a,3(a1,f7.4,a1,f7.4,a1))') 'Polarization: ', &
          '(', TOFLOAT(lasers(il)%pol(1)), ',', aimag(lasers(il)%pol(1)), '), ', &
          '(', TOFLOAT(lasers(il)%pol(2)), ',', aimag(lasers(il)%pol(2)), '), ', &
          '(', TOFLOAT(lasers(il)%pol(3)), ',', aimag(lasers(il)%pol(3)), ')'
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
      if(lasers(il)%field  ==  E_FIELD_ELECTRIC .or. lasers(il)%field  ==  E_FIELD_VECTOR_POTENTIAL) then
        fluence = M_ZERO
        max_intensity = M_ZERO
        maxfield= M_ZERO   
        do iter = 1, max_iter_
          tt = iter * dt_
          field = M_ZERO
          call laser_electric_field(lasers(il), field, tt, dt_)
          intensity = M_ZERO
          do idir = 1, MAX_DIM
            intensity = intensity + CNST(5.4525289841210) * field(idir)**2
          end do
          fluence = fluence + intensity
          if(intensity > max_intensity) max_intensity = intensity

          tmp = sum(field(:)**2)
          if( tmp > maxfield) maxfield = tmp
        end do
        fluence = fluence * dt_

        write(iunit,'(a,es17.6,3a)') '   Peak intensity       = ', max_intensity, ' [a.u]'
        write(iunit,'(a,es17.6,3a)') '                        = ', &
          max_intensity * CNST(6.4364086e+15), ' [W/cm^2]'
        write(iunit,'(a,es17.6,a)')  '   Int. intensity       = ', fluence, ' [a.u]'
        write(iunit,'(a,es17.6,a)')  '   Fluence              = ', &
          fluence / CNST(5.4525289841210) , ' [a.u]'

        if(abs(lasers(il)%omega) > M_EPSILON)then
          ! Ponderomotive Energy is the cycle-averaged kinetic energy of 
          ! a free electron quivering in the field 
          ! Up = E^2/(4*\omega^2)
          !
          ! subroutine laser_to_numerical_all sets lasers%omega to zero
          !
          Up = maxfield/(4*lasers(il)%omega**2)

          write(iunit,'(a,es17.6,3a)') '   Ponderomotive energy = ', &
            units_from_atomic(units_out%energy, Up) ,& 
            ' [', trim(units_abbrev(units_out%energy)), ']'
        end if
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
      pot(1:mesh%np) = pot(1:mesh%np) + TOFLOAT(amp)*laser%v(1:mesh%np)
    case default
      field(1:mesh%sb%dim) = TOFLOAT(amp*laser%pol(1:mesh%sb%dim))
      do ip = 1, mesh%np
        ! The -1 sign is missing here. Check epot.F90 for the explanation.
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
      amp = TOFLOAT(tdf(laser%f, time)*exp(M_zI*(laser%omega*time + tdf(laser%phi, time))))
      do idir = 1, mesh%sb%dim
        do ip = 1, mesh%np
          aa(ip, idir) = aa(ip, idir) + amp*laser%a(ip, idir)
        end do
      end do
    else
      do idir = 1, mesh%sb%dim
        do ip = 1, mesh%np
          aa(ip, idir) = aa(ip, idir) + laser%a(ip, idir)
        end do
      end do
    end if

    POP_SUB(laser_vector_potential)
  end subroutine laser_vector_potential
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Retrieves the value of either the electric or the magnetic
  !! field. If the laser is given by a scalar potential, the field
  !! should be a function of space (the gradient of the scalar potential
  !! times the temporal dependence), but in that case the subroutine
  !! just returns the temporal function.
  subroutine laser_field(laser, field, time)
    type(laser_t),     intent(in)    :: laser
    FLOAT,             intent(inout) :: field(:)
    FLOAT, optional,   intent(in)    :: time

    integer :: dim
    CMPLX :: amp

    !no PUSH SUB, called too often

    dim = size(field)    

    if(present(time)) then
      amp = tdf(laser%f, time) * exp(M_zI * ( laser%omega * time + tdf(laser%phi, time) ) )
    else
      amp = M_z1
    end if
    if(laser%field  ==  E_FIELD_SCALAR_POTENTIAL) then
      ! In this case we will just return the value of the time function. The "field", in fact, 
      ! should be a function of the position in space (thus, a true "field"), given by the 
      ! gradient of the scalar potential.
      field(1) = field(1) + TOFLOAT(amp)
    else
      field(1:dim) = field(1:dim) + TOFLOAT(amp*laser%pol(1:dim))
    end if

  end subroutine laser_field


  ! ---------------------------------------------------------
  !> Returns a vector with the electric field, no matter whether the laser is described directly as
  !! an electric field, or with a vector potential in the velocity gauge.
  subroutine laser_electric_field(laser, field, time, dt)
    type(laser_t),     intent(in)    :: laser
    FLOAT,             intent(inout) :: field(:)
    FLOAT,             intent(in)    :: time
    FLOAT,             intent(in)    :: dt

    integer :: dim
    FLOAT, allocatable :: field1(:), field2(:)

    !no PUSH SUB, called too often

    dim = size(field)

    select case(laser%field)
    case(E_FIELD_ELECTRIC)
      field = M_ZERO
      call laser_field(laser, field(1:dim), time)
    case(E_FIELD_VECTOR_POTENTIAL)
      SAFE_ALLOCATE(field1(1:dim))
      SAFE_ALLOCATE(field2(1:dim))
      field1 = M_ZERO
      field2 = M_ZERO
      call laser_field(laser, field1(1:dim), time - dt)
      call laser_field(laser, field2(1:dim), time + dt)
      field = - (field2 - field1) / (M_TWO * P_C * dt)
      SAFE_DEALLOCATE_A(field1)
      SAFE_DEALLOCATE_A(field2)
    case default
      field = M_ZERO
    end select

  end subroutine laser_electric_field
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vlaser_operator_quadratic(laser, der, psi, hpsi)
    type(laser_t),       intent(in)    :: laser
    type(derivatives_t), intent(in)    :: der
    CMPLX,               intent(inout) :: psi(:,:)  !< psi(der%mesh%np_part, h%d%dim)
    CMPLX,               intent(inout) :: hpsi(:,:) !< hpsi(der%mesh%np_part, h%d%dim)

    integer :: ip
    logical :: vector_potential, magnetic_field

    FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), bb(1:MAX_DIM), b_prime(1:MAX_DIM)
    FLOAT, allocatable :: aa(:, :), a_prime(:, :)

    PUSH_SUB(vlaser_operator_quadratic)

    a_field = M_ZERO

    vector_potential = .false.
    magnetic_field = .false.

    select case (laser_kind(laser))
    case (E_FIELD_ELECTRIC) ! do nothing
    case (E_FIELD_MAGNETIC)
      if (.not. allocated(aa)) then 
        SAFE_ALLOCATE(aa(1:der%mesh%np_part, 1:der%mesh%sb%dim))
        aa = M_ZERO
        SAFE_ALLOCATE(a_prime(1:der%mesh%np_part, 1:der%mesh%sb%dim))
        a_prime = M_ZERO
      end if
      a_prime = M_ZERO
      call laser_vector_potential(laser, der%mesh, a_prime)
      aa = aa + a_prime
      b_prime = M_ZERO
      call laser_field(laser, b_prime(1:der%mesh%sb%dim))
      bb = bb + b_prime
      magnetic_field = .true.
    case (E_FIELD_VECTOR_POTENTIAL)
      a_field_prime = M_ZERO
      call laser_field(laser, a_field_prime(1:der%mesh%sb%dim))
      a_field = a_field + a_field_prime
      vector_potential = .true.
    end select

    if (magnetic_field) then
      do ip = 1, der%mesh%np
        hpsi(ip, :) = hpsi(ip, :) + M_HALF * &
          dot_product(aa(ip, 1:der%mesh%sb%dim), aa(ip, 1:der%mesh%sb%dim)) * psi(ip, :) / P_c**2
      end do
      SAFE_DEALLOCATE_A(aa)
      SAFE_DEALLOCATE_A(a_prime)
    end if
    if (vector_potential) then
      do ip = 1, der%mesh%np
        hpsi(ip, :) = hpsi(ip, :) + M_HALF * &
          dot_product(a_field(1:der%mesh%sb%dim), a_field(1:der%mesh%sb%dim))*psi(ip, :) / P_c**2
      end do
    end if

    POP_SUB(vlaser_operator_quadratic)
  end subroutine vlaser_operator_quadratic

  ! ---------------------------------------------------------
  subroutine vlaser_operator_linear(laser, der, std, psi, hpsi, ik, gyromagnetic_ratio, a_static)
    type(laser_t),           intent(in)      :: laser
    type(derivatives_t),     intent(in)      :: der
    type(states_elec_dim_t), intent(in) :: std
    CMPLX,                   intent(inout)   :: psi(:,:) 
    CMPLX,                   intent(inout)   :: hpsi(:,:)
    integer,                 intent(in)      :: ik
    FLOAT,                   intent(in)      :: gyromagnetic_ratio
    FLOAT, optional,         intent(in)      :: a_static(:,:)

    integer :: ip, idim
    logical :: electric_field, vector_potential, magnetic_field
    CMPLX, allocatable :: grad(:, :, :), lhpsi(:, :)

    FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), bb(1:MAX_DIM), b_prime(1:MAX_DIM)

    FLOAT, allocatable :: vv(:), pot(:), aa(:, :), a_prime(:, :)

    PUSH_SUB(vlaser_operator_linear)

    a_field = M_ZERO

    electric_field = .false.
    vector_potential = .false.
    magnetic_field = .false.

    select case (laser_kind(laser))
    case (E_FIELD_SCALAR_POTENTIAL)
      if (.not. allocated(vv)) then 
        SAFE_ALLOCATE(vv(1:der%mesh%np))
      end if
      vv = M_ZERO
      call laser_potential(laser, der%mesh, vv)
      electric_field = .true.

    case (E_FIELD_ELECTRIC)
      if (.not. allocated(vv)) then 
        SAFE_ALLOCATE(vv(1:der%mesh%np))
        vv = M_ZERO
        SAFE_ALLOCATE(pot(1:der%mesh%np))
      end if
      pot = M_ZERO
      call laser_potential(laser, der%mesh, pot)
      vv = vv + pot
      electric_field = .true.
      SAFE_DEALLOCATE_A(pot)

    case (E_FIELD_MAGNETIC)
      if (.not. allocated(aa)) then 
        SAFE_ALLOCATE(aa(1:der%mesh%np_part, 1:der%mesh%sb%dim))
        aa = M_ZERO
        SAFE_ALLOCATE(a_prime(1:der%mesh%np_part, 1:der%mesh%sb%dim))
        a_prime = M_ZERO
      end if
      a_prime = M_ZERO
      call laser_vector_potential(laser, der%mesh, a_prime)
      aa = aa + a_prime
      b_prime = M_ZERO
      call laser_field(laser, b_prime(1:der%mesh%sb%dim))
      bb = bb + b_prime
      magnetic_field = .true.
    case (E_FIELD_VECTOR_POTENTIAL)
      a_field_prime = M_ZERO
      call laser_field(laser, a_field_prime(1:der%mesh%sb%dim))
      a_field = a_field + a_field_prime
      vector_potential = .true.
    end select

    if (electric_field) then
      do idim = 1, std%dim
        hpsi(1:der%mesh%np, idim)= hpsi(1:der%mesh%np, idim) + vv(1:der%mesh%np) * psi(1:der%mesh%np, idim)
      end do
      SAFE_DEALLOCATE_A(vv)
    end if

    if (magnetic_field) then
      SAFE_ALLOCATE(grad(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:std%dim))
 
      do idim = 1, std%dim
        call zderivatives_grad(der, psi(:, idim), grad(:, :, idim))
      end do

      ! If there is a static magnetic field, its associated vector potential is coupled with
      ! the time-dependent one defined as a "laser" (ideally one should just add them all and
      ! do the calculation only once...). Note that h%ep%a_static already has been divided
      ! by P_c, and therefore here we only divide by P_c, and not P_c**2.
      !
      ! We put a minus sign, since for the moment vector potential for
      ! lasers and for the static magnetic field use a different
      ! convetion.
      if (present(a_static)) then
        do ip = 1, der%mesh%np
          hpsi(ip, :) = hpsi(ip, :) - dot_product(aa(ip, 1:der%mesh%sb%dim), a_static(ip, 1:der%mesh%sb%dim)) * psi(ip, :) / P_c
        end do
      end if

      select case (std%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        do ip = 1, der%mesh%np
          hpsi(ip, 1) = hpsi(ip, 1) - M_zI * dot_product(aa(ip, 1:der%mesh%sb%dim), grad(ip, 1:der%mesh%sb%dim, 1)) / P_c
        end do
      case (SPINORS)
        do ip = 1, der%mesh%np
          do idim = 1, std%dim
            hpsi(ip, idim) = hpsi(ip, idim) - M_zI * &
              dot_product(aa(ip, 1:der%mesh%sb%dim), grad(ip, 1:der%mesh%sb%dim, idim)) / P_c
          end do
        end do
      end select


      select case (std%ispin)
      case (SPIN_POLARIZED)
        SAFE_ALLOCATE(lhpsi(1:der%mesh%np, 1:std%dim))
        if(modulo(ik+1, 2) == 0) then ! we have a spin down
          lhpsi(1:der%mesh%np, 1) = - M_HALF / P_c * sqrt(dot_product(bb, bb)) * psi(1:der%mesh%np, 1)
        else
          lhpsi(1:der%mesh%np, 1) = + M_HALF / P_c * sqrt(dot_product(bb, bb)) * psi(1:der%mesh%np, 1)
        end if
        hpsi(1:der%mesh%np, :) = hpsi(1:der%mesh%np, :) + (gyromagnetic_ratio * M_HALF) * lhpsi(1:der%mesh%np, :)
        SAFE_DEALLOCATE_A(lhpsi)

      case (SPINORS)
        SAFE_ALLOCATE(lhpsi(1:der%mesh%np, 1:std%dim))
        lhpsi(1:der%mesh%np, 1) = M_HALF / P_c * (bb(3) * psi(1:der%mesh%np, 1) &
             + (bb(1) - M_zI * bb(2)) * psi(1:der%mesh%np, 2))
        lhpsi(1:der%mesh%np, 2) = M_HALF / P_c * (-bb(3) * psi(1:der%mesh%np, 2) &
             + (bb(1) + M_zI * bb(2)) * psi(1:der%mesh%np, 1))
        hpsi(1:der%mesh%np, :) = hpsi(1:der%mesh%np, :) + (gyromagnetic_ratio * M_HALF) * lhpsi(1:der%mesh%np, :)
        SAFE_DEALLOCATE_A(lhpsi)
      end select

      SAFE_DEALLOCATE_A(grad)
      SAFE_DEALLOCATE_A(aa)
      SAFE_DEALLOCATE_A(a_prime)
    end if

    if (vector_potential) then
      SAFE_ALLOCATE(grad(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:std%dim))

      do idim = 1, std%dim
        call zderivatives_grad(der, psi(:, idim), grad(:, :, idim))
      end do

      select case(std%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        do ip = 1, der%mesh%np
          hpsi(ip, 1) = hpsi(ip, 1) - M_zI * dot_product(a_field(1:der%mesh%sb%dim), grad(ip, 1:der%mesh%sb%dim, 1)) / P_c
        end do
      case (SPINORS)
        do ip = 1, der%mesh%np
          do idim = 1, std%dim
            hpsi(ip, idim) = hpsi(ip, idim) - M_zI * &
              dot_product(a_field(1:der%mesh%sb%dim), grad(ip, 1:der%mesh%sb%dim, idim)) / P_c
          end do
        end do
      end select
      SAFE_DEALLOCATE_A(grad)
    end if

    POP_SUB(vlaser_operator_linear)
  end subroutine vlaser_operator_linear

end module lasers_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
