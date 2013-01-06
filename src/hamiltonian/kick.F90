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
!! $Id: kick.F90 $

#include "global.h"

module kick_m
  use c_pointer_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use ion_dynamics_m
  use loct_math_m
  use math_m
  use mesh_m
  use messages_m
  use parser_m
  use profiling_m
  use species_m
  use states_m
  use states_dim_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::               &
    kick_t,               &
    kick_init,            &
    kick_read,            &
    kick_write,           &
    kick_apply,           &
    kick_function_get


  integer, public, parameter ::        &
    KICK_FUNCTION_DIPOLE        = 0,   &
    KICK_FUNCTION_MULTIPOLE     = 1,   &
    KICK_FUNCTION_USER_DEFINED  = 2

  integer, public, parameter ::    &
    KICK_DENSITY_MODE        = 0,  &
    KICK_SPIN_MODE           = 1,  &
    KICK_SPIN_DENSITY_MODE   = 2

  integer, public, parameter ::    &
    QKICKMODE_NONE           = 0,  &
    QKICKMODE_EXP            = 1,  &
    QKICKMODE_COS            = 2,  &
    QKICKMODE_SIN            = 3,  &
    QKICKMODE_BESSEL         = 4


  type kick_t
    !> The time which the kick is applied (normally, this is zero)
    FLOAT             :: time
    !> The strength, and strength "mode".
    integer           :: delta_strength_mode
    FLOAT             :: delta_strength
    !> In case we use a normal dipole kick:
    FLOAT             :: pol(MAX_DIM, MAX_DIM)
    integer           :: pol_dir
    integer           :: pol_equiv_axes
    FLOAT             :: wprime(MAX_DIM)
    !> In case we have a general multipolar kick,
    !! the form of this "kick" will be (atomic units):
    !! \f[
    !! V(\vec{r}) = sum_{i=1}^{n\_multipoles} 
    !!                weight(i) * (e^2 / a_0^(l+1)) * r^l(i) * Y_{l(i),m(i)} (\vec{r})
    !! \f]
    !! which has units of energy; if we include the time-dependence (delta function):
    !! \f[    
    !! V(\vec{r}) = sum_{i=1}^{n\_multipoles} 
    !!                 weight(i) * (\hbar / a_0^l) * r^l(i) * Y_{l(i),m(i)} (\vec{r}) * \delta(t)
    !! \f]
    integer           :: n_multipoles
    integer, pointer  :: l(:), m(:)
    FLOAT, pointer    :: weight(:)
    FLOAT             :: qvector(MAX_DIM)
    FLOAT             :: qlength
    integer           :: qkick_mode
    integer           :: qbessel_l, qbessel_m
    !> In case we use a general function
    integer           :: function_mode
    character(len=200):: user_defined_function
  end type kick_t

contains

  ! ---------------------------------------------------------
  subroutine kick_init(kick, nspin, dim)
    type(kick_t), intent(out) :: kick
    integer,      intent(in)  :: nspin
    integer,      intent(in)  :: dim

    type(block_t) :: blk
    integer :: n_rows, irow, idir

    PUSH_SUB(kick_init)

    !%Variable TDDeltaKickTime
    !%Type float
    !%Default 0.0
    !%Section Time-Dependent::Response
    !%Description
    !% The delta-perturbation that can be applied by making use of the <tt>TDDeltaStrength</tt> variable,
    !% can be applied at the time given by this variable. Usually, this time is zero, since one wants
    !% to apply the delta pertubation or "kick" at a system at equilibrium, and no other time-dependent
    !% external potential is used. However, one may want to apply a kick on top of a laser field,
    !% for example.
    !%End
    call parse_float(datasets_check('TDDeltaKickTime'), M_ZERO, kick%time, units_inp%time)
    if(kick%time > M_ZERO) then
      call messages_experimental('TDDeltaKickTime > 0')
    end if

    !%Variable TDDeltaStrength
    !%Type float
    !%Default 0
    !%Section Time-Dependent::Response
    !%Description
    !% When no laser is applied, a delta (in time) perturbation with
    !% strength <tt>TDDeltaStrength</tt> can be applied. This is used to 
    !% calculate, <i>e.g.</i>, the linear optical spectra. If the ions are
    !% allowed to move, the kick will affect them also.
    !% The electric field is -(\hbar <i>k</i> / <i>e</i>) delta(<i>t</i>) for a dipole with
    !% zero wavevector, where <i>k</i> = <tt>TDDeltaStrength</tt>, which causes
    !% the wavefunctions instantaneously to acquire a phase exp(<i>ikx</i>).
    !% The unit is inverse length.
    !%End
    call parse_float(datasets_check('TDDeltaStrength'), M_ZERO, kick%delta_strength, units_inp%length**(-1))

    if(abs(kick%delta_strength) == M_ZERO) then
      kick%delta_strength_mode = 0
      kick%pol_equiv_axes = 0
      kick%pol(1:3, 1) = (/ M_ONE, M_ZERO, M_ZERO /)
      kick%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
      kick%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
      kick%pol_dir = 0
      kick%wprime = M_ZERO
      kick%n_multipoles = 0
      kick%qkick_mode = QKICKMODE_NONE
      nullify(kick%l)
      nullify(kick%m)
      nullify(kick%weight)
      POP_SUB(kick_init)
      return
    end if

    !%Variable TDDeltaStrengthMode
    !%Type integer
    !%Default kick_density
    !%Section Time-Dependent::Response
    !%Description
    !% When calculating the density response via real-time propagation,
    !% one needs to perform an initial kick on the KS system, at
    !% time zero. Depending on what kind of response property one wants to obtain,
    !% this kick may be done in several modes.
    !%Option kick_density 0
    !% The total density of the system is perturbed.
    !%Option kick_spin 1
    !% The individual spin densities are perturbed differently. Note that this mode
    !% is only possible if the run is done in spin-polarized mode, or with spinors.
    !%Option kick_spin_and_density 2
    !% A combination of the two above. Note that this mode
    !% is only possible if the run is done in spin-polarized mode, or with spinors.
    !%End
    call parse_integer(datasets_check('TDDeltaStrengthMode'), KICK_DENSITY_MODE, kick%delta_strength_mode)
    select case (kick%delta_strength_mode)
    case (KICK_DENSITY_MODE)
    case (KICK_SPIN_MODE, KICK_SPIN_DENSITY_MODE)
      if (nspin == 1) call input_error('TDDeltaStrengthMode')
    case default
      call input_error('TDDeltaStrengthMode')
    end select
    call messages_print_var_option(stdout, 'TDDeltaStrengthMode', kick%delta_strength_mode)

    nullify(kick%l)
    nullify(kick%m)
    nullify(kick%weight)

    kick%function_mode = KICK_FUNCTION_DIPOLE

    if(parse_isdef(datasets_check('TDDeltaUserDefined')).ne.0) then

      kick%function_mode = KICK_FUNCTION_USER_DEFINED
      kick%n_multipoles = 0

      !%Variable TDDeltaUserDefined
      !%Type string
      !%Section Time-Dependent::Response
      !%Description
      !% By default, the kick function will be a dipole. This will change if (1) the variable
      !% TDDeltaUserDefined is present in the inp file, or (2) if the block TDKickFunction 
      !% is present in the inp file. If both are present in the inp file, the TDKickFunction
      !% block will be ignored. The value of TDDeltaUserDefined should be a string describing
      !% the function that is going to be used as delta perturbation.
      !%End
      call parse_string(datasets_check('TDDeltaUserDefined'), "0", kick%user_defined_function)

      !%Variable TDKickFunction
      !%Type block
      !%Section Time-Dependent::Response
      !%Description
      !% If the block <tt>TDKickFunction</tt> is present in the input file, and the variable
      !% "TDDeltaUserDefined" is not present in the input file, the kick function to
      !% be applied at time zero of the time-propagation will not be a "dipole" function
      !% (<i>i.e.</i> phi => exp(i*k*z) phi), but a general multipole in the form r^l * Y_{lm}(r).
      !%
      !% Each line has two columns of integers: the (<i>l</i>,<i>m</i>) pair that defines the
      !% multipole. Any number of lines may be given, and the kick will be the sum of those
      !% multipoles.
      !%
      !% This feature allows calculation of quadrupole, octupole, etc., response functions.
      !%End
    else if(parse_block(datasets_check('TDKickFunction'), blk) == 0) then

      kick%function_mode = KICK_FUNCTION_MULTIPOLE
      n_rows = parse_block_n(blk)
      kick%n_multipoles = n_rows
      SAFE_ALLOCATE(     kick%l(1:n_rows))
      SAFE_ALLOCATE(     kick%m(1:n_rows))
      SAFE_ALLOCATE(kick%weight(1:n_rows))
      do irow = 1, n_rows
        call parse_block_integer(blk, irow - 1, 0, kick%l(irow))
        call parse_block_integer(blk, irow - 1, 1, kick%m(irow))
        call parse_block_float(blk, irow - 1, 2, kick%weight(irow))
        if( (kick%l(irow) < 0) .or. (abs(kick%m(irow)) > abs(kick%l(irow))) ) call input_error('TDkickFunction')
      end do

    else

      kick%function_mode = KICK_FUNCTION_DIPOLE
      kick%n_multipoles = 0

    end if
    

    ! Find out how many equivalent axes we have...
    !%Variable TDPolarizationEquivAxes
    !%Type integer
    !%Default 0
    !%Section Time-Dependent::Response
    !%Description
    !% Defines how many of the <tt>TDPolarization</tt> axes are equivalent. This information is stored in a file and then
    !% used by <tt>oct-propagation_spectrum</tt> to rebuild the full polarizability tensor from just the
    !% first <tt>TDPolarizationEquivAxes</tt> directions. This variable is also used by <tt>CalculationMode = vdw</tt>.
    !%End
    call parse_integer(datasets_check('TDPolarizationEquivAxes'), 0, kick%pol_equiv_axes)

    
    !%Variable TDPolarizationDirection
    !%Type integer
    !%Default 1
    !%Section Time-Dependent::Response
    !%Description
    !% When a delta potential is included in a time-dependent run, this
    !% variable defines in which direction the field will be applied
    !% by selecting one of the lines of <tt>TDPolarization</tt>. In a
    !% typical run (without using symmetry), the <tt>TDPolarization</tt> block
    !% would contain the three Cartesian unit vectors (the default
    !% value), and one would make 3 runs varying
    !% <tt>TDPolarization</tt> from 1 to 3.
    !% If one is using symmetry,  <tt>TDPolarization</tt> should run only from 1
    !% to <tt>TDPolarizationEquivAxes</tt>.
    !%End

    call parse_integer(datasets_check('TDPolarizationDirection'), 0, kick%pol_dir)

    if(kick%pol_dir < 1 .or. kick%pol_dir > dim) call input_error('TDPolarizationDirection')

    !%Variable TDPolarization
    !%Type block
    !%Section Time-Dependent::Response
    !%Description
    !% The (real) polarization of the delta electric field. Normally
    !% one needs three perpendicular polarization directions to calculate a
    !% spectrum (unless symmetry is used).
    !% The format of the block is:
    !%
    !% <tt>%TDPolarization
    !% <br>&nbsp;&nbsp;pol1x | pol1y | pol1z
    !% <br>&nbsp;&nbsp;pol2x | pol2y | pol2z
    !% <br>&nbsp;&nbsp;pol3x | pol3y | pol3z
    !% <br>%</tt>
    !%
    !% <tt>Octopus</tt> uses both this block and the variable
    !% <tt>TDPolarizationDirection</tt> to determine the polarization
    !% vector for the run. For example, if
    !% <tt>TDPolarizationDirection=2</tt> the polarization <tt>(pol2x,
    !% pol2y, pol2z)</tt> would be used.
    !%
    !% The default value for <tt>TDPolarization</tt> is the three
    !% Cartesian unit vectors (1,0,0), (0,1,0), and (0,0,1).
    !%
    !% Note that the directions do not necessarily need to be perpendicular
    !% when symmetries are used.
    !%
    !% WARNING: If you want to obtain the cross-section tensor, the
    !% <tt>TDPolarization</tt> block must be exactly the same for the run in
    !% each direction. The direction must be selected by the
    !% <tt>TDPolarizationDirection</tt> variable.
    !%
    !%End

    kick%pol(:, :) = M_ZERO
    if(parse_block(datasets_check('TDPolarization'), blk)==0) then
      n_rows = parse_block_n(blk)

      if(n_rows < dim) call input_error('TDPolarization')

      do irow = 1, n_rows
        do idir = 1, 3
          call parse_block_float(blk, irow - 1, idir - 1, kick%pol(idir, irow))
        end do
      end do
      if(n_rows < 3) kick%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
      if(n_rows < 2) kick%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
      call parse_block_end(blk)
    else
      ! Here the symmetry of the system should be analyzed, and the polarization
      ! basis built accordingly.
      kick%pol(1:3, 1) = (/ M_ONE,  M_ZERO, M_ZERO /)
      kick%pol(1:3, 2) = (/ M_ZERO, M_ONE,  M_ZERO /)
      kick%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE  /)
    end if

    ! Normalize
    do idir = 1, 3
      kick%pol(1:3, idir) = kick%pol(1:3, idir) / sqrt(sum(kick%pol(1:3, idir)**2))
    end do

    !%Variable TDPolarizationWprime
    !%Type block
    !%Section Time-Dependent::Response
    !%Description 
    !% This block is needed only when
    !% <tt>TDPolarizationEquivAxes</tt> is set to 3.  In such a case,
    !% the three directions (<i>pol1</i>, <i>pol2</i>, and <i>pol3</i>) defined in
    !% the <tt>TDPolarization</tt> block should be related by symmetry
    !% operations. If <i>A</i> is the symmetry operation that takes you
    !% from <i>pol1</i> to <i>pol2</i>, then <tt>TDPolarizationWprime</tt> 
    !% should be set to the direction defined by <i>A</i>^{-1} <i>pol3</i>.  
    !% For more information see MJT Oliveira
    !% <i>et al.</i>, <i>J. Nanoscience and Nanotechnology</i> <b>8</b>,
    !% 3392 (2008).
    !%End
    if(parse_block(datasets_check('TDPolarizationWprime'), blk)==0) then
      do idir = 1, 3
        call parse_block_float(blk, 0, idir - 1, kick%wprime(idir))
      end do
      kick%wprime(1:3) = kick%wprime(1:3) / sqrt(sum(kick%wprime(1:3)**2))
      call parse_block_end(blk)
    else
      kick%wprime(1:3) = (/ M_ZERO, M_ZERO, M_ONE /)
    end if

    !%Variable TDMomentumTransfer
    !%Type block
    !%Section Time-Dependent::Response
    !%Description
    !% Momentum-transfer vector for the calculation of the dynamic structure factor.
    !% When this variable is set, a non-dipole field is applied, and an output file
    !% <tt>ftchd</tt> is created (it contains the Fourier transform of the charge density
    !% at each time). The type of the applied external field can be set by
    !% an optional last number. Possible options are <tt>qexp</tt> (default), <tt>qcos</tt>,
    !% <tt>qsin</tt>, or <tt>qcos+qsin</tt>.
    !%Option qexp 1
    !% External field is exp(<i>iq.r</i>).
    !%Option qcos 2
    !% External field is cos(<i>q.r</i>).
    !%Option qsin 3
    !% External field is sin(<i>q.r</i>).
    !%Option qbessel 4
    !% External field is j_l(qr)*Y_lm(r), where q is the length of the momentum-transfer vector.
    !% In this case the block has to include two extra values (l and m).
    !%End

    if(parse_block(datasets_check('TDMomentumTransfer'), blk)==0) then
      do idir = 1, MAX_DIM
        call parse_block_float(blk, 0, idir - 1, kick%qvector(idir))
        kick%qvector(idir) = units_to_atomic(unit_one / units_inp%length, kick%qvector(idir))
      end do

      ! Read the calculation mode (exp, cos, sin, or bessel)
      if(parse_block_cols(blk, 0).gt.MAX_DIM) then

        call parse_block_integer(blk, 0, idir - 1, kick%qkick_mode)

        ! Read l and m if bessel mode (j_l*Y_lm) is used
        if(kick%qkick_mode.eq.QKICKMODE_BESSEL .and. parse_block_cols(blk, 0).eq.MAX_DIM+3) then
          call parse_block_integer(blk, 0, idir + 0, kick%qbessel_l)
          call parse_block_integer(blk, 0, idir + 1, kick%qbessel_m)
        else
          kick%qbessel_l = 0
          kick%qbessel_m = 0
        end if

      else
        kick%qkick_mode = QKICKMODE_EXP
      end if

      call parse_block_end(blk)
    else
      kick%qkick_mode = QKICKMODE_NONE
      kick%qvector(:) = M_ZERO
    end if

    kick%qlength = sqrt(sum(kick%qvector(:)**2))

    POP_SUB(kick_init)
  end subroutine kick_init


  ! ---------------------------------------------------------
  subroutine kick_read(kick, iunit)
    type(kick_t), intent(inout) :: kick
    integer,      intent(in)    :: iunit

    integer :: im, ierr
    character(len=100) :: line

    PUSH_SUB(kick_read)

    read(iunit, '(15x,i2)')     kick%delta_strength_mode
    read(iunit, '(15x,f18.12)') kick%delta_strength
    read(iunit, '(a)') line
    if(index(line,'defined').ne.0) then
      kick%function_mode = KICK_FUNCTION_USER_DEFINED
      read(line,'("# User defined:",1x,a)') kick%user_defined_function
    elseif(index(line,'multipole').ne.0) then
      kick%function_mode = KICK_FUNCTION_MULTIPOLE
      read(line, '("# N multipoles ",i3)') kick%n_multipoles
      SAFE_ALLOCATE(     kick%l(1:kick%n_multipoles))
      SAFE_ALLOCATE(     kick%m(1:kick%n_multipoles))
      SAFE_ALLOCATE(kick%weight(1:kick%n_multipoles))
      do im = 1, kick%n_multipoles
        read(iunit, '("# multipole    ",2i3,f18.12)') kick%l(im), kick%m(im), kick%weight(im)
      end do
    else
      kick%function_mode = KICK_FUNCTION_DIPOLE
      kick%n_multipoles = 0
      nullify(kick%l)
      nullify(kick%m)
      backspace(iunit)

      read(iunit, '(15x,3f18.12)') kick%pol(1:3, 1)
      read(iunit, '(15x,3f18.12)') kick%pol(1:3, 2)
      read(iunit, '(15x,3f18.12)') kick%pol(1:3, 3)
      read(iunit, '(15x,i2)')      kick%pol_dir
      read(iunit, '(15x,i2)')      kick%pol_equiv_axes
      read(iunit, '(15x,3f18.12)') kick%wprime(1:3)
    end if
    read(iunit, '(15x,f18.12)', iostat = ierr) kick%time
    if(ierr.ne.0) then
      kick%time = M_ZERO
      backspace(iunit)
    end if


    POP_SUB(kick_read)
  end subroutine kick_read


  ! ---------------------------------------------------------
  subroutine kick_write(kick, iunit, out)
    type(kick_t),          intent(in) :: kick
    integer,    optional,  intent(in) :: iunit
    type(c_ptr), optional, intent(in) :: out

    integer :: im
    character(len=120) :: aux

    PUSH_SUB(kick_write)

    if(present(iunit)) then
      write(iunit, '(a15,i1)')      '# kick mode    ', kick%delta_strength_mode
      write(iunit, '(a15,f18.12)')  '# kick strength', kick%delta_strength
       ! if this were to be read by humans, we would want units_from_atomic(units_out%length**(-1))
      if(kick%function_mode .eq. KICK_FUNCTION_USER_DEFINED) then
        write(iunit,'(a15,1x,a)')     '# User defined:', trim(kick%user_defined_function)
      elseif(kick%n_multipoles > 0) then
        write(iunit, '(a15,i3)')    '# N multipoles ', kick%n_multipoles
        do im = 1, kick%n_multipoles
          write(iunit, '(a15,2i3,f18.12)') '# multipole    ', kick%l(im), kick%m(im), kick%weight(im)
        end do
      else
        write(iunit, '(a15,3f18.12)') '# pol(1)       ', kick%pol(1:3, 1)
        write(iunit, '(a15,3f18.12)') '# pol(2)       ', kick%pol(1:3, 2)
        write(iunit, '(a15,3f18.12)') '# pol(3)       ', kick%pol(1:3, 3)
        write(iunit, '(a15,i1)')      '# direction    ', kick%pol_dir
        write(iunit, '(a15,i1)')      '# Equiv. axes  ', kick%pol_equiv_axes
        write(iunit, '(a15,3f18.12)') '# wprime       ', kick%wprime(1:3)
      end if
      write(iunit, '(a15,f18.12)') "# kick time    ", kick%time

    else if(present(out)) then
      write(aux, '(a15,i2)')      '# kick mode    ', kick%delta_strength_mode
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      write(aux, '(a15,f18.12)')  '# kick strength', kick%delta_strength
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      if(kick%function_mode .eq. KICK_FUNCTION_USER_DEFINED) then
        write(aux,'(a15,1x,a)')     '# User defined:', trim(kick%user_defined_function)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
      elseif(kick%n_multipoles > 0) then
        write(aux, '(a15,i3)')      '# N multipoles ', kick%n_multipoles
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        do im = 1, kick%n_multipoles
          write(aux, '(a15,2i3,f18.12)') '# multipole    ', kick%l(im), kick%m(im), kick%weight(im)
          call write_iter_string(out, aux)
          call write_iter_nl(out)
        end do
      else
        write(aux, '(a15,3f18.12)') '# pol(1)       ', kick%pol(1:3, 1)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,3f18.12)') '# pol(2)       ', kick%pol(1:3, 2)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,3f18.12)') '# pol(3)       ', kick%pol(1:3, 3)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,i2)')      '# direction    ', kick%pol_dir
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,i2)')      '# Equiv. axes  ', kick%pol_equiv_axes
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,3f18.12)') '# wprime       ', kick%wprime(1:3)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
      end if
      write(aux, '(a15,f18.12)') "# kick time    ", kick%time
      call write_iter_string(out, aux)
      call write_iter_nl(out)

    end if

    POP_SUB(kick_write)
  end subroutine kick_write


  ! ---------------------------------------------------------
  ! 
  subroutine kick_function_get(gr, kick, kick_function)
    type(grid_t),         intent(in)    :: gr
    type(kick_t),         intent(in)    :: kick
    CMPLX,                intent(out)   :: kick_function(:)

    integer :: ip, im
    FLOAT   :: xx(MAX_DIM)
    FLOAT   :: rkick, ikick, gylm(1:MAX_DIM), rr, ylm

    PUSH_SUB(kick_function_get)

    if(abs(kick%qlength) > M_EPSILON) then ! q-vector is set

      select case (kick%qkick_mode)
        case (QKICKMODE_COS)
           write(message(1), '(a,3F9.5,a)') 'Info: Using cos(q.r) field with q = (', kick%qvector(:), ')'
        case (QKICKMODE_SIN)
           write(message(1), '(a,3F9.5,a)') 'Info: Using sin(q.r) field with q = (', kick%qvector(:), ')'
        case (QKICKMODE_SIN + QKICKMODE_COS)
           write(message(1), '(a,3F9.5,a)') 'Info: Using sin(q.r)+cos(q.r) field with q = (', kick%qvector(:), ')'
        case (QKICKMODE_EXP)
           write(message(1), '(a,3F9.5,a)') 'Info: Using exp(iq.r) field with q = (', kick%qvector(:), ')'
        case (QKICKMODE_BESSEL)
           write(message(1), '(a,I2,a,I2,a,F9.5)') 'Info: Using j_l(qr)*Y_lm(r) field with (l,m)= (', &
                                                    kick%qbessel_l, ",", kick%qbessel_m,') and q = ', kick%qlength
        case default
           write(message(1), '(a,3F9.6,a)') 'Info: Unknown field type!'
      end select
      call messages_info(1)

      kick_function = M_z0
      do ip = 1, gr%mesh%np
        call mesh_r(gr%mesh, ip, rr, coords = xx)
        select case (kick%qkick_mode)
          case (QKICKMODE_COS)
            kick_function(ip) = kick_function(ip) + cos(sum(kick%qvector(:) * xx(:)))
          case (QKICKMODE_SIN)
            kick_function(ip) = kick_function(ip) + sin(sum(kick%qvector(:) * xx(:)))
          case (QKICKMODE_SIN+QKICKMODE_COS)
            kick_function(ip) = kick_function(ip) + sin(sum(kick%qvector(:) * xx(:)))
          case (QKICKMODE_EXP)
            kick_function(ip) = kick_function(ip) + exp(M_zI * sum(kick%qvector(:) * xx(:)))
          case (QKICKMODE_BESSEL)
            call grylmr(gr%mesh%x(ip, 1), gr%mesh%x(ip, 2), gr%mesh%x(ip, 3), kick%qbessel_l, kick%qbessel_m, ylm, gylm)
              kick_function(ip) = kick_function(ip) + loct_sph_bessel(kick%qbessel_l, kick%qlength*sqrt(sum(xx(:)**2)))*ylm
        end select
      end do

    else
      if(kick%function_mode .eq. KICK_FUNCTION_USER_DEFINED) then

        kick_function = M_z0
        do ip = 1, gr%mesh%np
          call mesh_r(gr%mesh, ip, rr, coords = xx)
            rkick = M_ZERO; ikick = M_ZERO
          call parse_expression(rkick, ikick, gr%sb%dim, xx, rr, M_ZERO, trim(kick%user_defined_function))
            kick_function(ip) = rkick
        end do

      elseif(kick%n_multipoles > 0) then

        kick_function = M_z0
        do im = 1, kick%n_multipoles
          do ip = 1, gr%mesh%np
            call mesh_r(gr%mesh, ip, rr, coords = xx)
            call loct_ylm(1, xx(1), xx(2), xx(3), kick%l(im), kick%m(im), ylm)
              kick_function(ip) = kick_function(ip) + kick%weight(im) * (rr**kick%l(im)) * ylm
          end do
        end do
      else
        forall(ip = 1:gr%mesh%np)
          kick_function(ip) = sum(gr%mesh%x(ip, 1:gr%mesh%sb%dim) * &
            kick%pol(1:gr%mesh%sb%dim, kick%pol_dir))
        end forall
      end if
    end if

    POP_SUB(kick_function_get)
  end subroutine kick_function_get


  ! ---------------------------------------------------------
  !> Applies the delta-function electric field \f$ E(t) = E_0 \Delta(t) \f$
  !! where \f$ E_0 = \frac{- k \hbar}{e} \f$ k = kick\%delta_strength.
  subroutine kick_apply(gr, st, ions, geo, kick)
    type(grid_t),         intent(in)    :: gr
    type(states_t),       intent(inout) :: st
    type(ion_dynamics_t), intent(in)    :: ions
    type(geometry_t),     intent(inout) :: geo
    type(kick_t),         intent(in)    :: kick

    integer :: iqn, ist, idim, ip, ispin, iatom
    CMPLX   :: cc(2), kick_value
    CMPLX, allocatable :: kick_function(:), psi(:, :)

    PUSH_SUB(kick_apply)

    ! The wavefunctions at time delta t read
    ! psi(delta t) = psi(t) exp(i k x)
    delta_strength: if(kick%delta_strength .ne. M_ZERO) then

        SAFE_ALLOCATE(kick_function(1:gr%mesh%np))
        call kick_function_get(gr, kick, kick_function)

        write(message(1),'(a,f11.6)') 'Info: Applying delta kick: k = ', kick%delta_strength
        select case (kick%function_mode)
        case (KICK_FUNCTION_DIPOLE)
          message(2) = "Info: kick function: dipole."
        case (KICK_FUNCTION_MULTIPOLE)
          message(2) = "Info: kick function: multipoles."
        case (KICK_FUNCTION_USER_DEFINED)
          message(2) = "Info: kick function: user defined function."
        end select
        select case (kick%delta_strength_mode)
        case (KICK_DENSITY_MODE)
          message(3) = "Info: Delta kick mode: Density mode"
        case (KICK_SPIN_MODE)
          message(3) = "Info: Delta kick mode: Spin mode"
        case (KICK_SPIN_DENSITY_MODE)
          message(3) = "Info: Delta kick mode: Density + Spin modes"
        end select
        call messages_info(3)

        SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))

        do iqn = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end

            call states_get_state(st, gr%mesh, ist, iqn, psi)

            select case (kick%delta_strength_mode)
            case (KICK_DENSITY_MODE)
              forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np)
                psi(ip, idim) = exp(M_zI*kick%delta_strength*kick_function(ip))*psi(ip, idim)
              end forall

            case (KICK_SPIN_MODE)
              ispin = states_dim_get_spin_index(st%d, iqn)
              do ip = 1, gr%mesh%np
                kick_value = M_zI*kick%delta_strength*kick_function(ip)

                cc(1) = exp(kick_value)
                cc(2) = exp(-kick_value)

                select case (st%d%ispin)
                case (SPIN_POLARIZED)
                  psi(ip, 1) = cc(ispin)*psi(ip, 1)
                case (SPINORS)
                  psi(ip, 1) = cc(1)*psi(ip, 1)
                  psi(ip, 2) = cc(2)*psi(ip, 2)
                end select
              end do

            case (KICK_SPIN_DENSITY_MODE)
              do ip = 1, gr%mesh%np
                cc(1) = exp(M_TWO*kick_value)
                select case (st%d%ispin)
                case (SPIN_POLARIZED)
                  if(is_spin_up(iqn)) then
                    psi(ip, 1) = cc(1)*psi(ip, 1)
                  end if
                case (SPINORS)
                  psi(ip, 1) = cc(1)*psi(ip, 1)
                end select
              end do
            end select

            call states_set_state(st, gr%mesh, ist, iqn, psi)

          end do
        end do

        SAFE_DEALLOCATE_A(psi)

        ! The nuclear velocities will be changed by
        ! Delta v_z = ( Z*e*E_0 / M) = - ( Z*k*\hbar / M)
        ! where M and Z are the ionic mass and charge, respectively.
        if(ion_dynamics_ions_move(ions)  .and. kick%delta_strength .ne. M_ZERO) then
          do iatom = 1, geo%natoms
            geo%atom(iatom)%v(1:gr%mesh%sb%dim) = geo%atom(iatom)%v(1:gr%mesh%sb%dim) + &
              kick%delta_strength * kick%pol(1:gr%mesh%sb%dim, kick%pol_dir) * &
              P_PROTON_CHARGE * species_zval(geo%atom(iatom)%spec) / &
              species_weight(geo%atom(iatom)%spec)
          end do
        end if

      SAFE_DEALLOCATE_A(kick_function)
    end if delta_strength

    POP_SUB(kick_apply)
  end subroutine kick_apply

end module kick_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
