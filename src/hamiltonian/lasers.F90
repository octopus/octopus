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
  use em_field_m
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
    laser_requires_gradient,      &
    dvlasers,                     &
    zvlasers,                     &
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
  FLOAT function laser_carrier_frequency(l) result(w0)
    type(laser_t), intent(in) :: l

    call push_sub('lasers.laser_carrier_frequency')
    w0 = l%omega

    call pop_sub()
  end function laser_carrier_frequency
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function laser_kind(l)
    type(laser_t), intent(in) :: l

    ! no push_sub allowed in pure function
    laser_kind = l%field

  end function laser_kind
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  function laser_polarization(l) result(pol)
    type(laser_t), intent(in) :: l
    CMPLX :: pol(MAX_DIM)

    call push_sub('lasers.laser_polarization')
    pol(1:MAX_DIM) = l%pol(1:MAX_DIM)

    call pop_sub()
  end function laser_polarization
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_get_f(l, f)
    type(laser_t), intent(in)    :: l
    type(tdf_t),   intent(inout) :: f

    call push_sub('lasers.laser_get')
    call tdf_copy(f, l%f)

    call pop_sub()
  end subroutine laser_get_f
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_set_f(l, f)
    type(laser_t), intent(inout) :: l
    type(tdf_t),   intent(inout) :: f

    call push_sub('lasers.laser_set_f')

    call tdf_end(l%f)
    call tdf_copy(l%f, f)

    call pop_sub()
  end subroutine laser_set_f
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_get_phi(l, phi)
    type(laser_t), intent(in)    :: l
    type(tdf_t),   intent(inout) :: phi

    call push_sub('lasers.laser_get_phi')
    call tdf_copy(phi, l%phi)

    call pop_sub()
  end subroutine laser_get_phi
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_set_phi(l, phi)
    type(laser_t), intent(inout) :: l
    type(tdf_t),   intent(inout) :: phi

    call push_sub('lasers.laser_set_phi')

    call tdf_end(l%phi)
    call tdf_copy(l%phi, phi)

    call pop_sub()
  end subroutine laser_set_phi
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_set_f_value(l, i, x)
    type(laser_t), intent(inout) :: l
    integer,       intent(in)    :: i
    FLOAT,         intent(in)    :: x

    call push_sub('lasers.laser_set_f_value')
    call tdf_set_numerical(l%f, i, x)

    call pop_sub()
  end subroutine laser_set_f_value
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! The td functions that describe the laser field are transformed to a 
  ! "numerical" representation (i.e. time grid, values at this time grid).
  ! The possible phase and carrier frequency are evaluated and put together with
  ! the envelope, so that the envelope describes the full function (zero phase,
  ! zero carrier frequency.
  ! ---------------------------------------------------------
  subroutine laser_to_numerical_all(l, dt, max_iter, omegamax)
    type(laser_t), intent(inout)  :: l
    FLOAT,         intent(in)     :: dt
    integer,       intent(in)     :: max_iter
    FLOAT,         intent(in)     :: omegamax

    integer :: j
    FLOAT   :: t, fj, phi

    call push_sub('lasers.lasers_to_numerical_all')

    call tdf_to_numerical(l%f, max_iter, dt, omegamax)
    do j = 1, max_iter + 1
      t = (j-1)*dt
      fj = tdf(l%f, j)
      phi = tdf(l%phi, t)
      call tdf_set_numerical(l%f, j, fj*cos(l%omega*t+phi))
    end do
    call tdf_end(l%phi)
    call tdf_init_cw(l%phi, M_ZERO, M_ZERO)
    l%omega = M_ZERO

    call pop_sub()
  end subroutine laser_to_numerical_all
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! The td functions that describe the laser field are transformed to a 
  ! "numerical" representation (i.e. time grid, values at this time grid).
  ! ---------------------------------------------------------
  subroutine laser_to_numerical(l, dt, max_iter, omegamax)
    type(laser_t), intent(inout)  :: l
    FLOAT,         intent(in)     :: dt
    integer,       intent(in)     :: max_iter
    FLOAT,         intent(in)     :: omegamax
    call push_sub('lasers.lasers_to_numerical')

    call tdf_to_numerical(l%f, max_iter, dt, omegamax)
    call tdf_to_numerical(l%phi, max_iter, dt, omegamax)

    call pop_sub()
  end subroutine laser_to_numerical
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_init(no_l, l, m)
    integer,     intent(out) :: no_l
    type(laser_t),   pointer :: l(:)
    type(mesh_t), intent(in) :: m

    type(block_t)     :: blk
    integer           :: i, j, ierr
    character(len=200) :: scalar_pot_expression
    character(len=200) :: envelope_expression
    FLOAT :: omega0, r, pot_re, pot_im
    FLOAT, allocatable :: x(:)

    call push_sub('lasers.laser_init')

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
    !%
    !% (A) type = <tt>electric field, magnetic field, vector_potential</tt>
    !%
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
    !%
    !% (B) type = <tt>scalar_potential</tt>
    !%
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; scalar_potential | "scalar_expression" | freq | envelope_function_name
    !% <br>%</tt>
    !%
    !% The scalar potential is not just a dipole, but any expression given by the string
    !% "scalar_expression". The temporal shape is determined by the envelope function
    !% defined by <tt>envelope_function_name</tt>.
    !% 
    !%
    !% A NOTE ON UNITS:
    !%
    !% It is very common to describe the strength of a laser field by its intensity, rather
    !% than using the electric-field amplitude. In atomic units (or, more precisely, in any
    !% Gaussian system of units), the relationship between instantaneous electric field
    !% and intensity is:
    !%
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
      SAFE_ALLOCATE(l(1:no_l))

      do i = 1, no_l

        call parse_block_integer(blk, i-1, 0, l(i)%field)

        select case(l(i)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          call parse_block_string(blk, i-1, 1, scalar_pot_expression)
          j = 1
          l(i)%pol = M_z1
        case default
          call parse_block_cmplx(blk, i-1, 1, l(i)%pol(1))
          call parse_block_cmplx(blk, i-1, 2, l(i)%pol(2))
          call parse_block_cmplx(blk, i-1, 3, l(i)%pol(3))
          j = 3
        end select

        call parse_block_float(blk, i-1, j+1, omega0)
        omega0 = units_to_atomic(units_inp%energy, omega0)

        l(i)%omega = omega0
     
        call parse_block_string(blk, i-1, j+2, envelope_expression)

        ! For some reason, one cannot open a block if another one is already open.
        ! This is why I close blk before calling tdf_read, and then open it again.
        ! This should be fixed at the parser level.
        call parse_block_end(blk)
        call tdf_read(l(i)%f, trim(envelope_expression), ierr)
        ierr = parse_block(datasets_check('TDExternalFields'), blk)

        ! Check if there is a phase.
        if(parse_block_cols(blk, i-1) > j+3) then
          call parse_block_string(blk, i-1, j+3, envelope_expression)
          call parse_block_end(blk)
          call tdf_read(l(i)%phi, trim(envelope_expression), ierr)
          ierr = parse_block(datasets_check('TDExternalFields'), blk)
        else
          call tdf_init(l(i)%phi)
        end if

        l(i)%pol(:) = l(i)%pol(:)/sqrt(sum(abs(l(i)%pol(:))**2))

        select case(l(i)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          SAFE_ALLOCATE(l(i)%v(1:m%np_part))
          l(i)%v = M_ZERO
          SAFE_ALLOCATE(x(1:MAX_DIM))
          do j = 1, m%np
            call mesh_r(m, j, r, x = x)
            call parse_expression(pot_re, pot_im, m%sb%dim, x, r, M_ZERO, trim(scalar_pot_expression))
            l(i)%v(j) = pot_re
          end do
          SAFE_DEALLOCATE_A(x)

        case(E_FIELD_MAGNETIC)
          ! WARNING: note that for the moment we are ignoring the possibility of a complex
          ! polarizability vector for the td magnetic field case.
          SAFE_ALLOCATE(l(i)%a(1:m%np_part, 1:m%sb%dim))
          l(i)%a = M_ZERO
          SAFE_ALLOCATE(x(1:m%sb%dim))
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
          SAFE_DEALLOCATE_A(x)

        end select

      end do

      call parse_block_end(blk)
    end if

    call pop_sub()

  end subroutine laser_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_end(no_l, l)
    type(laser_t), pointer :: l(:)
    integer,    intent(in) :: no_l

    integer :: i

    call push_sub('laser.laser_end')

    if(no_l > 0) then
      do i = 1, no_l
        call tdf_end(l(i)%f)
        call tdf_end(l(i)%phi)
        select case(l(i)%field)
        case(E_FIELD_SCALAR_POTENTIAL)
          SAFE_DEALLOCATE_P(l(i)%v); nullify(l(i)%v)
        case(E_FIELD_MAGNETIC)
          SAFE_DEALLOCATE_P(l(i)%a); nullify(l(i)%a)
        end select
      end do
      SAFE_DEALLOCATE_P(l); nullify(l)
    end if

    call pop_sub()
  end subroutine laser_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_write_info(l, iunit, dt, max_iter)
    type(laser_t), intent(in) :: l(:)
    integer,       intent(in) :: iunit
    FLOAT, optional, intent(in) :: dt 
    integer, optional, intent(in) :: max_iter 

    FLOAT :: t, fluence, max_intensity, intensity, dt_
    CMPLX :: amp, val
    integer :: i, j, k, no_l, max_iter_

    if(.not.mpi_grp_is_root(mpi_world)) return

    call push_sub('lasers.laser_write_info')

    no_l = size(l)
    
    do i = 1, no_l

      if(present(dt)) then
        dt_ = dt
      else
        dt_ = tdf_dt(l(i)%f)
      end if
      if(present(max_iter)) then
        max_iter_ = max_iter
      else
        max_iter_ = tdf_niter(l(i)%f)
      end if

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
      write(iunit,'(3x,a,f14.8,3a)') 'Carrier frequency = ', &
        units_from_atomic(units_out%energy, l(i)%omega), &
        ' [', trim(units_abbrev(units_out%energy)), ']'
      write(iunit,'(3x,a)')       'Envelope: ' 
      call tdf_write(l(i)%f, iunit)

      if(.not.tdf_is_empty(l(i)%phi)) then
        write(iunit,'(3x,a)')       'Phase: ' 
        call tdf_write(l(i)%phi, iunit)
      end if

      ! 1 atomic unit of intensity = 3.5094448e+16 W / cm^2
      ! In a Gaussian system of units,
      ! I(t) = (1/(8\pi)) * c * E(t)^2
      ! (1/(8\pi)) * c = 5.4525289841210 a.u.
      if(l(i)%field .eq. E_FIELD_ELECTRIC) then
        fluence = M_ZERO

        max_intensity = M_ZERO
        do j = 1, max_iter_
          t = j * dt_
          val = tdf(l(i)%f, t)
          amp = val*exp(M_zI*(l(i)%omega*t + tdf(l(i)%phi, t)))
          intensity = M_ZERO
          do k = 1, MAX_DIM
            intensity = intensity + CNST(5.4525289841210) * real(amp*l(i)%pol(k))**2
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

    call pop_sub()
  end subroutine laser_write_info
  ! ---------------------------------------------------------


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
      amp = tdf(l%f, t) * exp(M_zI * ( l%omega*t + tdf(l%phi, t) ) )
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


  ! ---------------------------------------------------------
  subroutine laser_vector_potential(l, a, t)
    type(laser_t),     intent(in) :: l
    FLOAT,            intent(out) :: a(:, :)
    FLOAT, optional,   intent(in) :: t

    call push_sub('lasers.laser_vector_potential')

    if(present(t)) then
      a = l%a * real( tdf(l%f, t) * exp(M_zI* ( l%omega *t + tdf(l%phi, t)  ) ) )
    else
      a = l%a
    end if

    call pop_sub()
  end subroutine laser_vector_potential
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine laser_field(sb, l, field, t)
    type(simul_box_t), intent(in)  :: sb
    type(laser_t),     intent(in)  :: l
    FLOAT,             intent(out) :: field(:)
    FLOAT, optional,   intent(in)  :: t

    CMPLX :: amp

    call push_sub('lasers.laser_field')

    field = M_ZERO
    if(l%field .eq. E_FIELD_SCALAR_POTENTIAL) then
      call pop_sub(); return
    endif
    if(present(t)) then
      amp = tdf(l%f, t) * exp(M_zI * ( l%omega * t + tdf(l%phi, t) ) )
    else
      amp = M_z1
    end if
    field(1:sb%dim) = field(1:sb%dim) + real(amp*l%pol(1:sb%dim))

    call pop_sub()
  end subroutine laser_field
  ! ---------------------------------------------------------

  logical elemental function laser_requires_gradient(this) result(req)
    type(laser_t),  intent(in)  :: this

    ! no push_sub allowed in elemental function
    req = (laser_kind(this) == E_FIELD_MAGNETIC .or. laser_kind(this) == E_FIELD_VECTOR_POTENTIAL)
    
  end function laser_requires_gradient

  ! ---------------------------------------------------------
  subroutine lasers_get_potentials(laser, mesh, time, em_field)
    type(laser_t),       intent(in)    :: laser
    type(mesh_t),        intent(inout) :: mesh
    FLOAT,               intent(in)    :: time
    type(em_field_t),    intent(inout) :: em_field
    
    FLOAT, allocatable :: pot(:, :)
    FLOAT :: mag(1:MAX_DIM)
    integer :: ip, idir

    call push_sub('lasers.lasers_get_potential')
    
    select case(laser_kind(laser))

    case(E_FIELD_SCALAR_POTENTIAL, E_FIELD_ELECTRIC)
      if(.not. associated(em_field%potential)) then 
        SAFE_ALLOCATE(em_field%potential(1:mesh%np))
        em_field%potential = M_ZERO
      end if
      
      SAFE_ALLOCATE(pot(1:mesh%np, 1:1))
      
      call laser_potential(mesh%sb, laser, mesh, pot(:, 1), time)
      
      forall(ip = 1:mesh%np) em_field%potential(ip) = em_field%potential(ip) + pot(ip, 1)

      SAFE_DEALLOCATE_A(pot)
      
    case(E_FIELD_MAGNETIC)
      if(.not. associated(em_field%vector_potential)) then 
        SAFE_ALLOCATE(em_field%vector_potential(1:mesh%np, 1:mesh%sb%dim))
        em_field%vector_potential = M_ZERO
      end if

      SAFE_ALLOCATE(pot(1:mesh%np, 1:mesh%sb%dim))

      call laser_vector_potential(laser, pot, time)
      call laser_field(mesh%sb, laser, mag, time)

      forall(idir = 1:mesh%sb%dim, ip = 1:mesh%np) 
        em_field%vector_potential(ip, idir) = em_field%vector_potential(ip, idir) + pot(ip, idir)
      end forall
      forall(idir = 1:mesh%sb%dim) em_field%uniform_magnetic_field(idir) = em_field%uniform_magnetic_field(idir) + mag(idir)

    case(E_FIELD_VECTOR_POTENTIAL)
      if(.not. associated(em_field%vector_potential)) then 
        SAFE_ALLOCATE(em_field%vector_potential(1:mesh%np, 1:MAX_DIM))
        em_field%vector_potential(mesh%np, 1:MAX_DIM) = M_ZERO
      end if
      
      call laser_field(mesh%sb, laser, mag, time)
      forall(idir = 1:mesh%sb%dim, ip = 1:mesh%np) 
        em_field%vector_potential(ip, idir) = em_field%vector_potential(ip, idir) + mag(idir)
      end forall
      
    end select
    
    call pop_sub()
  end subroutine lasers_get_potentials

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
