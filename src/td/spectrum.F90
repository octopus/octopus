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

module spectrum_m
  use c_pointer_m
  use datasets_m
  use global_m
  use io_m
  use lalg_adv_m
  use loct_math_m
  use parser_m
  use math_m
  use messages_m
  use profiling_m
  use string_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                        &
    spec_t,                        &
    kick_t,                        &
    spectrum_init,                 &
    spectrum_cross_section,        &
    spectrum_cross_section_tensor, &
    spectrum_rotatory_strength,    &
    spectrum_hs_from_mult,         &
    spectrum_hs_from_acc,          &
    spectrum_hsfunction_init,      &
    spectrum_hsfunction_end,       &
    spectrum_hsfunction_min,       &
    spectrum_mult_info,            &
    spectrum_fix_time_limits,      &
    count_time_steps,              &
    kick_init,                     &
    kick_write,                    &
    kick_read

  integer, public, parameter ::    &
    SPECTRUM_DAMP_NONE       = 0,  &
    SPECTRUM_DAMP_LORENTZIAN = 1,  &
    SPECTRUM_DAMP_POLYNOMIAL = 2,  &
    SPECTRUM_DAMP_GAUSSIAN   = 3

  integer, public, parameter ::    &
    SPECTRUM_TRANSFORM_EXP   = 1,  &
    SPECTRUM_TRANSFORM_SIN   = 2,  &
    SPECTRUM_TRANSFORM_COS   = 3

  integer, public, parameter ::    &
    KICK_DENSITY_MODE        = 0,  &
    KICK_SPIN_MODE           = 1,  &
    KICK_SPIN_DENSITY_MODE   = 2

  type spec_t
    FLOAT   :: start_time          ! start time for the transform
    FLOAT   :: end_time            ! when to stop the transform
    FLOAT   :: energy_step         ! step in energy mesh
    FLOAT   :: max_energy          ! maximum of energy mesh
    integer :: damp                ! Damp type (none, exp or pol)
    integer :: transform           ! sine, cosine, or exponential transform
    FLOAT   :: damp_factor         ! factor used in damping
  end type spec_t

  type kick_t
    integer           :: delta_strength_mode
    FLOAT             :: delta_strength
    ! In case we use a normal dipole kick:
    FLOAT             :: pol(3, 3)
    integer           :: pol_dir
    integer           :: pol_equiv_axes
    FLOAT             :: wprime(3)
    ! In case we have a general multipolar kick
    ! The form of this "kick" will be (atomic units):
    ! V(\vec{r}) = sum_{i=1}^{n_multipoles} 
    !                 weight(i) * (e^2 / a_0^(l+1)) * r^l(i) * Y_{l(i),m(i)} (\vec{r})
    ! which has units of energy; if we include the time-dependence (delta function):
    ! V(\vec{r}) = sum_{i=1}^{n_multipoles} 
    !                 weight(i) * (\hbar / a_0^l) * r^l(i) * Y_{l(i),m(i)} (\vec{r}) * \delta(t)
    integer           :: n_multipoles
    integer, pointer  :: l(:), m(:)
    FLOAT, pointer    :: weight(:)
  end type kick_t

  ! Module variables, necessary to compute the function hsfunction, called by
  ! the C function loct_1dminimize
  FLOAT :: time_step_
  CMPLX, allocatable :: func_(:)
  integer :: is_, ie_

contains

  ! ---------------------------------------------------------
  subroutine spectrum_init(s)
    type(spec_t), intent(inout) :: s

    call push_sub('spectrum.spectrum_init')

    !%Variable SpecDampMode
    !%Type integer
    !%Default polynomial
    !%Section Utilities::oct-cross-section
    !%Description
    !% Decides which damping/filtering is to be applied in order to calculate
    !% spectra by calculating a Fourier transform.
    !%Option no 0
    !% No filtering at all.
    !%Option exponential 1
    !% Exponential filtering, corresponding to a Lorentzian-shaped spectrum.
    !%Option polynomial 2
    !% Third-order polynomial damping.
    !%Option gaussian 3
    !% Gaussian damping.
    !%End
    call parse_integer  (datasets_check('SpecDampMode'), SPECTRUM_DAMP_POLYNOMIAL, s%damp)
    if(.not.varinfo_valid_option('SpecDampMode', s%damp)) call input_error('SpecDampMode')

    !%Variable SpecTransform
    !%Type integer
    !%Default sine
    !%Section Utilities::oct-cross-section
    !%Description
    !% Decides which transform to perform.
    !%Option sine 2
    !% Sine transform <math>\int dt \sin(wt) f(t)</math>
    !%Option cosine 3
    !% Cosine transform <math>\int dt \cos(wt) f(t)</math>
    !%Option exponential 1
    !% Exponential transform <math>\int dt \exp(-wt) f(t)</math>
    !%End
    call parse_integer  (datasets_check('SpecTransform'), SPECTRUM_TRANSFORM_SIN, s%transform)
    if(.not.varinfo_valid_option('SpecTransform', s%transform)) call input_error('SpecTransform')

    !%Variable SpecStartTime
    !%Type float
    !%Default 0.0
    !%Section Utilities::oct-cross-section
    !%Description
    !% Processing is done for the given function in a time-window that starts at the
    !% value of this variable.
    !%End
    call parse_float(datasets_check('SpecStartTime'),  M_ZERO, s%start_time, units_inp%time)

    !%Variable SpecEndTime
    !%Type float
    !%Default -1.0 au
    !%Section Utilities::oct-cross-section
    !%Description
    !% Processing is done for the given function in a time-window that ends at the
    !% value of this variable. If set to a negative value, the maximum value from 
    !% the corresponding multipole file will used.
    !%End
    call parse_float(datasets_check('SpecEndTime'), -M_ONE, s%end_time, units_inp%time)

    !%Variable SpecEnergyStep
    !%Type float
    !%Default 0.01 eV
    !%Section Utilities::oct-cross-section
    !%Description
    !% Sampling rate for the spectrum.
    !%End
    call parse_float(datasets_check('SpecEnergyStep'), CNST(0.01)/(M_TWO*P_Ry), s%energy_step, units_inp%energy)
    

    !%Variable SpecMaxEnergy
    !%Type float
    !%Default 20 eV
    !%Section Utilities::oct-cross-section
    !%Description
    !% The Fourier transform is calculated for energies smaller than this value.
    !%End
    call parse_float(datasets_check('SpecMaxEnergy'), CNST(20.0)/(M_TWO*P_Ry), s%max_energy, units_inp%energy)

    !%Variable SpecDampFactor
    !%Type float
    !%Default 0.15 au
    !%Section Utilities::oct-cross-section
    !%Description
    !% If <tt>SpecDampMode = exponential</tt>, the damping parameter of the exponential
    !% is fixed through this variable.
    !%End
    call parse_float(datasets_check('SpecDampFactor'),  CNST(0.15), s%damp_factor, units_inp%time**(-1))

    call pop_sub()
  end subroutine spectrum_init


  ! ---------------------------------------------------------
  subroutine kick_write(k, iunit, out)
    type(kick_t),         intent(in) :: k
    integer,   optional,  intent(in) :: iunit
    type(c_ptr), optional,  intent(in) :: out

    integer :: i
    character(len=120) :: aux

    call push_sub('spectrum.kick_write')

    if(present(iunit)) then
      write(iunit, '(a15,i1)')      '# kick mode    ', k%delta_strength_mode
      write(iunit, '(a15,f18.12)')  '# kick strength', k%delta_strength
      if(k%n_multipoles > 0) then
        write(iunit, '(a15,i3)')      '# N multipoles ', k%n_multipoles
        do i = 1, k%n_multipoles
          write(iunit, '(a15,2i3,f18.12)')     '# multipole    ', k%l(i), k%m(i), k%weight(i)
        end do
      else
        write(iunit, '(a15,3f18.12)') '# pol(1)       ', k%pol(1:3, 1)
        write(iunit, '(a15,3f18.12)') '# pol(2)       ', k%pol(1:3, 2)
        write(iunit, '(a15,3f18.12)') '# pol(3)       ', k%pol(1:3, 3)
        write(iunit, '(a15,i1)')      '# direction    ', k%pol_dir
        write(iunit, '(a15,i1)')      '# Equiv. axes  ', k%pol_equiv_axes
        write(iunit, '(a15,3f18.12)') '# wprime       ', k%wprime(1:3)
      end if

    else if(present(out)) then
      write(aux, '(a15,i2)')      '# kick mode    ', k%delta_strength_mode
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      write(aux, '(a15,f18.12)')  '# kick strength', k%delta_strength
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      if(k%n_multipoles > 0) then
        write(aux, '(a15,i3)')      '# N multipoles ', k%n_multipoles
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        do i = 1, k%n_multipoles
          write(aux, '(a15,2i3,f18.12)') '# multipole    ', k%l(i), k%m(i), k%weight(i)
          call write_iter_string(out, aux)
          call write_iter_nl(out)
        end do
      else
        write(aux, '(a15,3f18.12)') '# pol(1)       ', k%pol(1:3, 1)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,3f18.12)') '# pol(2)       ', k%pol(1:3, 2)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,3f18.12)') '# pol(3)       ', k%pol(1:3, 3)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,i2)')      '# direction    ', k%pol_dir
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,i2)')      '# Equiv. axes  ', k%pol_equiv_axes
        call write_iter_string(out, aux)
        call write_iter_nl(out)
        write(aux, '(a15,3f18.12)') '# wprime       ', k%wprime(1:3)
        call write_iter_string(out, aux)
        call write_iter_nl(out)
      end if
    end if

    call pop_sub()
  end subroutine kick_write


  ! ---------------------------------------------------------
  subroutine kick_read(k, iunit)
    type(kick_t), intent(inout) :: k
    integer,      intent(in) :: iunit

    integer :: i
    character(len=100) :: line

    call push_sub('spectrum.kick_read')

    read(iunit, '(15x,i2)')      k%delta_strength_mode
    read(iunit, '(15x,f18.12)')  k%delta_strength
    read(iunit, '(a)') line
    if(index(line,'multipole').ne.0) then
      read(line, '("# N multipoles ",i3)') k%n_multipoles
      SAFE_ALLOCATE(     k%l(1:k%n_multipoles))
      SAFE_ALLOCATE(     k%m(1:k%n_multipoles))
      SAFE_ALLOCATE(k%weight(1:k%n_multipoles))
      do i = 1, k%n_multipoles
        read(iunit, '("# multipole    ",2i3,f18.12)') k%l(i), k%m(i), k%weight(i)
      end do
    else
      k%n_multipoles = 0
      nullify(k%l)
      nullify(k%m)
      backspace(iunit)

      read(iunit, '(15x,3f18.12)') k%pol(1:3, 1)
      read(iunit, '(15x,3f18.12)') k%pol(1:3, 2)
      read(iunit, '(15x,3f18.12)') k%pol(1:3, 3)
      read(iunit, '(15x,i2)')      k%pol_dir
      read(iunit, '(15x,i2)')      k%pol_equiv_axes
      read(iunit, '(15x,3f18.12)') k%wprime(1:3)
    end if

    call pop_sub()
  end subroutine kick_read


  ! ---------------------------------------------------------
  subroutine kick_init(k, nspin, dim)
    type(kick_t), intent(out) :: k
    integer,      intent(in)  :: nspin
    integer,      intent(in)  :: dim

    type(block_t) :: blk
    integer :: n_rows, i, j

    call push_sub('spectrum.kick_init')

    !%Variable TDDeltaStrength
    !%Type float
    !%Default 0.0
    !%Section Time-Dependent::Linear Response
    !%Description
    !% When no laser is applied, a delta (in time) perturbation with
    !% strength <tt>TDDeltaStrength</tt> can be applied. This is used to 
    !% calculate, <i>e.g.</i>, the linear optical spectra.
    !%
    !% Note that the "strength" here described is non-dimensional.
    !%End
    call parse_float(datasets_check('TDDeltaStrength'), M_ZERO, k%delta_strength)

    if(abs(k%delta_strength) == M_ZERO) then
      k%delta_strength_mode = 0
      k%pol_equiv_axes = 0
      k%pol(1:3, 1) = (/ M_ONE, M_ZERO, M_ZERO /)
      k%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
      k%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
      k%pol_dir = 0
      k%wprime = M_ZERO
      k%n_multipoles = 0
      nullify(k%l)
      nullify(k%m)
      nullify(k%weight)
      call pop_sub(); return
    end if

    !%Variable TDDeltaStrengthMode
    !%Type integer
    !%Default kick_density
    !%Section Time-Dependent::Linear Response
    !%Description
    !% When calculating the linear response of the density via the propagation
    !% in real time, one needs to perform an initial kick on the KS system, at
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
    call parse_integer(datasets_check('TDDeltaStrengthMode'), KICK_DENSITY_MODE, k%delta_strength_mode)
    select case (k%delta_strength_mode)
    case (KICK_DENSITY_MODE)
    case (KICK_SPIN_MODE, KICK_SPIN_DENSITY_MODE)
      if (nspin == 1) call input_error('TDDeltaStrengthMode')
    case default
      call input_error('TDDeltaStrengthMode')
    end select
    call messages_print_var_option(stdout, 'TDDeltaStrengthMode', k%delta_strength_mode)


    !%Variable TDKickFunction
    !%Type block
    !%Section Time-Dependent::Linear Response
    !%Description
    !% If the block <tt>TDKickFunction</tt> is present in the input file, the kick function to
    !% be applied at time zero of the time-propagation will not be a "dipole" function
    !% (i.e. phi => exp(i*k*z)phi), but a general multipole in the form r^l * Y_{lm}(r).
    !%
    !% The block <tt>TDKickFunction</tt> shall only contain one line, with two columns that shall
    !% be of integer type: those two integers will be the (<i>l</i>,<i>m</i>) pair that defines the
    !% multipole.
    !%
    !% This feature allows calculation of quadrupole, octupole, etc., response functions.
    !%End
    if(parse_block(datasets_check('TDKickFunction'), blk)==0) then
      n_rows = parse_block_n(blk)
      k%n_multipoles = n_rows
      SAFE_ALLOCATE(     k%l(1:n_rows))
      SAFE_ALLOCATE(     k%m(1:n_rows))
      SAFE_ALLOCATE(k%weight(1:n_rows))
      do i = 1, n_rows
        call parse_block_integer(blk, i-1, 0, k%l(i))
        call parse_block_integer(blk, i-1, 1, k%m(i))
        call parse_block_float(blk, i-1, 2, k%weight(i))
        if( (k%l(i) < 0) .or. (abs(k%m(i)) > abs(k%l(i))) ) call input_error('TDkickFunction')
      end do
    else
      k%n_multipoles = 0
      nullify(k%l)
      nullify(k%m)
      nullify(k%weight)
    end if
    

    ! Find out how many equivalent axes we have...
    !%Variable TDPolarizationEquivAxes
    !%Type integer
    !%Default 0
    !%Section Time-Dependent::Linear Response
    !%Description
    !% Defines how many of the <tt>TDPolarization</tt> axes are equivalent. This information can then
    !% be used by <tt>oct-cross-section</tt> to rebuild the full polarizability tensor from just the
    !% first <tt>TDPolarizationEquivAxes</tt> directions.
    !%End
    call parse_integer(datasets_check('TDPolarizationEquivAxes'), 0, k%pol_equiv_axes)

    
    !%Variable TDPolarizationDirection
    !%Type integer
    !%Default 1
    !%Section Time-Dependent::Linear Response
    !%Description
    !%
    !% When a delta potential is included in a time-dependent run, this
    !% variable defines in which direction the field will be applied
    !% by selecting one of the lines of <tt>TDPolarization</tt>. In a
    !% typical run (without using symmetry), the <tt>TDPolarization</tt> block
    !% would contain the three Cartesian unit vectors (the default
    !% value), and one would make 3 runs varying
    !% <tt>TDPolarization</tt> from 1 to 3.
    !%
    !%End

    call parse_integer(datasets_check('TDPolarizationDirection'), 0, k%pol_dir)

    if(k%pol_dir < 1 .or. k%pol_dir > dim) call input_error('TDPolarizationDirection')

    !%Variable TDPolarization
    !%Type block
    !%Section Time-Dependent::Linear Response
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
    !% WARNING: If you want to obtain the cross-section tensor, the
    !% <tt>TDPolarization</tt> block must be exactly the same for the run in
    !% each direction. The direction must be selected by the
    !% <tt>TDPolarizationDirection</tt> variable.
    !%
    !%End

    k%pol(:, :) = M_ZERO
    if(parse_block(datasets_check('TDPolarization'), blk)==0) then
      n_rows = parse_block_n(blk)

      if(n_rows < dim) call input_error('TDPolarization')

      do j = 1, n_rows
        do i = 1, 3
          call parse_block_float(blk, j-1, i-1, k%pol(i, j))
        end do
      end do
      if(n_rows<3) k%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
      if(n_rows<2) k%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
      call parse_block_end(blk)
    else
      ! Here the symmetry of the system should be analyzed, and the polarization
      ! basis built accordingly.
      k%pol(1:3, 1) = (/ M_ONE, M_ZERO, M_ZERO /)
      k%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
      k%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
    end if

    ! Normalize
    do i = 1, 3
      k%pol(1:3, i) = k%pol(1:3, i)/sqrt(sum(k%pol(1:3, i)**2))
    end do

    !%Variable TDPolarizationWprime
    !%Type block
    !%Section Time-Dependent::Linear Response
    !%Description
    !% Say you have a first symmetry operation (A)
    !% that takes you the first axis (p1) to the second axis (p2), and then
    !% a second symmetry operation (B) that takes you the second axis (p2) to the
    !% third (p3). Then wprime = A^{-1} p3.
    !%End
    if(parse_block(datasets_check('TDPolarizationWprime'), blk)==0) then
      do i = 1, 3
        call parse_block_float(blk, 0, i-1, k%wprime(i))
      end do
      k%wprime(1:3) = k%wprime(1:3)/sqrt(sum(k%wprime(1:3)**2))
    else
      k%wprime(1:3) = (/ M_ZERO, M_ZERO, M_ONE /)
    end if

    call pop_sub()
  end subroutine kick_init


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section_tensor(s, out_file, in_file)
    type(spec_t), intent(inout) :: s
    integer,      intent(in)    :: out_file
    integer,      intent(in)    :: in_file(:)

    character(len=20) :: header_string
    integer :: nspin, energy_steps, i, is, j, equiv_axes, n_files, k
    FLOAT, allocatable :: sigma(:, :, :, :), sigmap(:, :, :, :), sigmau(:, :, :),  &
      sigmav(:, :, :), sigmaw(:, :, :), p(:, :), ip(:, :)
    FLOAT :: dw, dump, average, anisotropy
    type(kick_t) :: kick

    call push_sub('spectrum.spectrum_cross_section_tensor')

    n_files = size(in_file)
    equiv_axes = 3 - n_files + 1

    call spectrum_cross_section_info(in_file(1), nspin, kick, energy_steps, dw)
    call io_skip_header(in_file(1))

    SAFE_ALLOCATE(sigma (1:3, 1:3, 0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmap(1:3, 1:3, 0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmau(1:3,      0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmav(1:3,      0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmaw(1:3,      0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(     p(1:3, 1:3))
    SAFE_ALLOCATE(    ip(1:3, 1:3))

    select case(equiv_axes)

    case(3)

      do i = 0, energy_steps
        read(in_file(1), *) dump, sigmau(1:3, i, 1:nspin)
      end do

      ! The first row of sigma is the vector that we have just read, but properly projected...
      do is = 1, nspin
        do i = 0, energy_steps
          sigmap(1, 1, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 1))
          sigmap(1, 2, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 2))
          sigmap(1, 3, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 3))
        end do
      end do

      ! The diagonal parts are also equal:
      sigmap(2, 2, :, :) = sigmap(1, 1, :, :)
      sigmap(3, 3, :, :) = Sigmap(1, 1, :, :)

      ! The (2,1) term and (3,1) term are equal by symmetry:
      sigmap(2, 1, :, :) = sigmap(1, 2, :, :)
      sigmap(3, 1, :, :) = sigmap(1, 3, :, :)

      ! But for the (2,3) term we need the wprime vector....
      do is = 1, nspin
        do i = 0, energy_steps
          sigmap(2, 3, i, is) = sum(sigmau(1:3, i, is)*kick%wprime(1:3))
          sigmap(3, 2, i, is) = sigmap(2, 3, i, is)
        end do
      end do

    case(2)

      call spectrum_cross_section_info(in_file(2), i, kick, j, dump)
      call io_skip_header(in_file(2))

      do i = 0, energy_steps
        read(in_file(1), *) dump, sigmau(1:3, i, 1:nspin)
        read(in_file(2), *) dump, sigmaw(1:3, i, 1:nspin)
      end do

      ! The first row of sigma is the vector that we have just read, but properly projected...
      do is = 1, nspin
        do i = 0, energy_steps
          sigmap(1, 1, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 1))
          sigmap(1, 2, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 2))
          sigmap(1, 3, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 3))
        end do
      end do

      ! The third row of sigma is also the vector that we have just read, but properly projected...
      do is = 1, nspin
        do i = 0, energy_steps
          sigmap(3, 1, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 1))
          sigmap(3, 2, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 2))
          sigmap(3, 3, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 3))
        end do
      end do

      ! The diagonal (2,2) is equal by symmetry to the (1,1)
      sigmap(2, 2, :, :) = sigmap(1, 1, :, :)

      ! The (2,1) term and (1,2) term are equal; the (2,3) and (3,2), also.
      sigmap(2, 1, :, :) = sigmap(1, 2, :, :)
      sigmap(2, 3, :, :) = sigmap(3, 2, :, :)

    case default

      call spectrum_cross_section_info(in_file(2), i, kick, j, dump)
      call spectrum_cross_section_info(in_file(3), i, kick, j, dump)
      call io_skip_header(in_file(2))
      call io_skip_header(in_file(3))

      do i = 0, energy_steps
        read(in_file(1), *) dump, sigmau(1:3, i, 1:nspin)
        read(in_file(2), *) dump, sigmav(1:3, i, 1:nspin)
        read(in_file(3), *) dump, sigmaw(1:3, i, 1:nspin)
      end do

      do is = 1, nspin
        do i = 0, energy_steps
          sigmap(1, 1, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 1))
          sigmap(1, 2, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 2))
          sigmap(1, 3, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 3))
        end do
      end do
      do is = 1, nspin
        do i = 0, energy_steps
          sigmap(2, 1, i, is) = sum( sigmav(1:3, i, is)*kick%pol(1:3, 1))
          sigmap(2, 2, i, is) = sum( sigmav(1:3, i, is)*kick%pol(1:3, 2))
          sigmap(2, 3, i, is) = sum( sigmav(1:3, i, is)*kick%pol(1:3, 3))
        end do
      end do
      do is = 1, nspin
        do i = 0, energy_steps
          sigmap(3, 1, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 1))
          sigmap(3, 2, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 2))
          sigmap(3, 3, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 3))
        end do
      end do

    end select

    ! And now, perform the necessary transformation.
    ip(1:3, 1:3) = kick%pol(1:3, 1:3)
    dump = lalg_inverter(3, ip)
    do is = 1, nspin
      do i = 0, energy_steps
        sigma(:, :, i, is) = matmul( transpose(ip), matmul(sigmap(:, :, i, is), ip) )
      end do
    end do

    ! Finally, write down the result
    write(out_file, '(a15,i2)')      '# nspin        ', nspin
    call kick_write(kick, out_file)
    write(out_file, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
    write(out_file, '(a20)', advance = 'no') str_center("(1/3)*Tr[sigma]", 20)
    write(out_file, '(a20)', advance = 'no') str_center("Anisotropy[sigma]", 20)
    do j = 1, nspin
      do i = 1, 3
        do k = 1, 3
          write(header_string,'(a6,i1,a1,i1,a1,i1,a1)') 'sigma(',i,',',k,',',j,')'
          write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
        end do
      end do
    end do
    write(out_file, '(1x)')
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20)
    do i = 1, 2+nspin*9
      write(out_file, '(a20)', advance = 'no')  str_center('['//trim(units_abbrev(units_out%length)) //'^2]', 20)
    end do
    write(out_file, '(1x)')

    ! The anisotropy (Delta alpha) of a second rank symmetric tensor alpha (such 
    ! as the cross section) is defined as:
    ! 
    ! (Delta alpha)^2 = (1/3) * ( 3*Tr(alpha^2) - (Tr(a))^2 )
    !
    ! The reason for this definition is that it is identically equal to:
    !
    ! (Delta alpha)^2 = (1/3) * ( (alpha_1-alpha_2)^2 + (alpha_1-alpha_3)^2 + (alpha_2-alpha_3)^2 )
    !
    ! where {alpha_1, alpha_2, alpha_3} are the eigenvalues of alpha. An "isotropic" tensor
    ! is characterized by having three equal eigenvalues, which leads to zero anisotropy. The
    ! more different that the eigenvalues are, the larger the anisotropy is.
    do i = 0, energy_steps

      p = M_ZERO
      do j = 1, min(2, nspin) ! we add spin up with spin down
         p(:, :) = p(:, :) + sigma(:, :, i, j)
      end do
      average = M_THIRD * ( p(1, 1) + p(2, 2) + p(3, 3) )
      ip = matmul(p, p)
      anisotropy = M_THIRD * ( M_THREE * ( ip(1, 1) + ip(2, 2) + ip(3, 3) ) - (M_THREE*average)**2 )

      ! Note that the cross-section elements do not have to be transformed to the proper units, since
      ! they have been read from the "cross_section_vector.x", that are already in the proper units.
      write(out_file,'(3e20.8)', advance = 'no') units_from_atomic(units_out%energy, (i*s%energy_step)), &
        average , sqrt(max(anisotropy, M_ZERO)) 
      do j = 1, nspin
        write(out_file,'(9e20.8)', advance = 'no') sigma(1:3, 1:3, i, j)
      end do
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(sigma)
    SAFE_DEALLOCATE_A(sigmap)
    SAFE_DEALLOCATE_A(sigmau)
    SAFE_DEALLOCATE_A(sigmav)
    SAFE_DEALLOCATE_A(sigmaw)
    SAFE_DEALLOCATE_A(p)
    SAFE_DEALLOCATE_A(ip)
    call pop_sub()
  end subroutine spectrum_cross_section_tensor


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section(in_file, out_file, s)
    integer, intent(in) :: in_file
    integer, intent(in) :: out_file
    type(spec_t),  intent(inout) :: s

    character(len=20) :: header_string
    integer :: nspin, lmax, time_steps, is, ie, ntiter, i, j, jj, isp, no_e, k, idir
    FLOAT   :: dt, dump, x, w, ewsum, polsum
    type(kick_t) :: kick
    FLOAT, allocatable :: dipole(:, :, :), sigma(:, :, :), dumpa(:), sf(:, :)
    type(unit_system_t) :: file_units

    call push_sub('spectrum.spectrum_cross_section')

    ! This function gives us back the unit connected to the "multipoles" file, the header information,
    ! the number of time steps, and the time step.
    call spectrum_mult_info(in_file, nspin, kick, time_steps, dt, file_units, lmax=lmax)

    ! Now we cannot process files that do not contain the dipole, or that contain more than the dipole.
    if(lmax.ne.1) then
      message(1) = 'multipoles file should contain the dipole -- and only the dipole.'
      call write_fatal(1)
    end if

    ! Find out the iteration numbers corresponding to the time limits.
    call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

    ! Read the dipole.
    call io_skip_header(in_file)
    SAFE_ALLOCATE(dipole(1:3, 0:time_steps, 1:nspin))
    do i = 0, time_steps
      select case(nspin)
      case(1)
        read(in_file, *) j, dump, dump, dipole(1:3, i, 1)
      case(2)
        read(in_file, *) j, dump, dump, dipole(1:3, i, 1), dump, dipole(1:3, i, 2)
      case(4)
        read(in_file, *) j, dump, dump, dipole(1:3, i, 1), dump, dipole(1:3, i, 2), &
          dump, dipole(1:3, i, 3), dump, dipole(1:3, i, 4)
      end select

      dipole(1:3, i, :) = units_to_atomic(file_units%length, dipole(1:3, i, :))
      
    end do

    ! Now subtract the initial dipole.
    do i = time_steps, 0, -1
      dipole(:, i, :) = dipole(:, i, :) - dipole(:, 0, :)
    end do

    ! Get the number of energy steps.
    no_e = s%max_energy / s%energy_step
    SAFE_ALLOCATE(sigma(1:3, 0:no_e, 1:nspin))
    SAFE_ALLOCATE(   sf(     0:no_e, nspin))
    sigma = M_ZERO
    sf    = M_ZERO

    ! Gets the damping function (here because otherwise it is awfully slow in "pol" mode...)_m
    SAFE_ALLOCATE(dumpa(is:ie))
    do j = is, ie
      jj = j - is
      select case(s%damp)
      case(SPECTRUM_DAMP_NONE)
        dumpa(j) = M_ONE
      case(SPECTRUM_DAMP_LORENTZIAN)
        dumpa(j)= exp(-jj*dt*s%damp_factor)
      case(SPECTRUM_DAMP_POLYNOMIAL)
        dumpa(j) = M_ONE - M_THREE*(real(jj)/ntiter)**2                          &
          + M_TWO*(real(jj)/ntiter)**3
      case(SPECTRUM_DAMP_GAUSSIAN)
        dumpa(j)= exp(-(jj*dt)**2*s%damp_factor**2)
      end select
    end do

    do k = 0, no_e
      w = k*s%energy_step
      do j = is, ie
        jj = j - is

        select case(s%transform)
        case(SPECTRUM_TRANSFORM_SIN)
          x = sin(w*jj*dt)
        case(SPECTRUM_TRANSFORM_COS)
          x = cos(w*jj*dt)
        case(SPECTRUM_TRANSFORM_EXP)
          x = exp(-w*jj*dt)
        end select

        do isp = 1, nspin
          sigma(1:3, k, isp) = sigma(1:3, k, isp) + x*dumpa(j)*dipole(1:3, j, isp)
          sf(k, isp) = sf(k, isp) + x*dumpa(j)*sum(dipole(1:3, j, isp)*kick%pol(1:3,kick%pol_dir))
        end do
      end do
      sigma(1:3, k, 1:nspin) = sigma(1:3, k, 1:nspin)*dt
      sf(k, 1:nspin) = sf(k, 1:nspin)*dt
      sf(k, 1:nspin) = - sf(k, 1:nspin) * (w*M_TWO)/(M_Pi*kick%delta_strength)
      sigma(1:3, k, 1:nspin) = -sigma(1:3, k, 1:nspin)*(M_FOUR*M_PI*w/P_c)/kick%delta_strength
    end do

    ewsum = sum(sf(0, 1:nspin)); polsum = M_ZERO
    do k = 1, no_e
      w = k*s%energy_step
      ewsum = ewsum + sum(sf(k, 1:nspin))
      polsum = polsum + sum(sf(k, 1:nspin))/w**2
    end do
    ewsum = ewsum * s%energy_step; polsum = polsum * s%energy_step

    write(out_file, '(a15,i2)')      '# nspin        ', nspin
    call kick_write(kick, out_file)
    write(out_file, '(a)') '#%'
    write(out_file, '(a,i8)')    '# Number of time steps = ', time_steps
    write(out_file, '(a,i4)')    '# SpecDampMode              = ', s%damp
    write(out_file, '(a,f10.4)') '# SpecDampFactor            = ', units_from_atomic(units_out%time**(-1), s%damp_factor)
    write(out_file, '(a,f10.4)') '# SpecStartTime             = ', units_from_atomic(units_out%time, s%start_time)
    write(out_file, '(a,f10.4)') '# SpecEndTime               = ', units_from_atomic(units_out%time, s%end_time)
    write(out_file, '(a,f10.4)') '# SpecMaxEnergy             = ', units_from_atomic(units_out%energy, s%max_energy)
    write(out_file, '(a,f10.4)') '# SpecEnergyStep            = ', units_from_atomic(units_out%energy, s%energy_step)
    write(out_file, '(a)') '#%'
    write(out_file, '(a,f16.6)') '# Electronic sum rule       = ', ewsum
    write(out_file, '(a,f16.6)') '# Polarizability (sum rule) = ', units_from_atomic(units_out%length**3, polsum)
    write(out_file, '(a)') '#%'
    
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center("Energy", 20)
    do j = 1, nspin
      do i = 1, 3
        write(header_string,'(a6,i1,a8,i1,a1)') 'sigma(',i,', nspin=',j,')'
        write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
      end do
    end do
    do j = 1, nspin
      write(header_string,'(a18,i1,a1)') 'StrengthFunction(', j, ')'
      write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
    end do
    write(out_file, '(1x)')
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20)
    do i = 1, nspin*3
      write(out_file, '(a20)', advance = 'no') str_center('['//trim(units_abbrev(units_out%length)) //'^2]', 20)
    end do
    do i = 1, nspin
      write(out_file, '(a20)', advance = 'no') str_center('[1/'//trim(units_abbrev(units_out%energy)) //']',20)
    end do
    write(out_file, '(1x)')

    do i = 0, no_e
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, i*s%energy_step)
      do j = 1, nspin
        write(out_file,'(3e20.8)', advance = 'no') (units_from_atomic(units_out%length**2, sigma(idir, i, j)), &
                                                    idir = 1, 3)
      end do
      do j = 1, nspin
        write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, sf(i, j))
      end do
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(dipole)
    SAFE_DEALLOCATE_A(sigma)
    call pop_sub()
  end subroutine spectrum_cross_section


  ! ---------------------------------------------------------
  subroutine spectrum_rotatory_strength(in_file, out_file, s)
    integer, intent(in) :: in_file
    integer, intent(in) :: out_file
    type(spec_t), intent(inout) :: s

    integer :: i, is, ie, ntiter, j, jj, k, time_steps, no_e, nspin
    FLOAT :: dump, dt, w
    type(kick_t) :: kick
    CMPLX :: z, sum1, sum2
    CMPLX, pointer :: sp(:)
    FLOAT, allocatable :: dumpa(:)
    FLOAT, allocatable :: angular(:, :)
    type(unit_system_t) :: file_units

    call push_sub('spectrum.spectrum_rotatory_strength')

    call spectrum_mult_info(in_file, nspin, kick, time_steps, dt, file_units)
    call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

    ! load dipole from file
    SAFE_ALLOCATE(angular(0:time_steps, 1:3))
    call io_skip_header(in_file)
    do i = 0, time_steps
      read(in_file, *) j, dump, angular(i, 1:3)
    end do

    ! subtract static dipole
    do i = 1, 3
      angular(:, i) = angular(:, i) - angular(0, i)
    end do

    no_e = s%max_energy / s%energy_step
    SAFE_ALLOCATE(sp(0:no_e))
    sp = M_z0
    sum1 = M_z0
    sum2 = M_z0

    ! Gets the damping function (here because otherwise it is awfully slow in "pol" mode...)
    SAFE_ALLOCATE(dumpa(is:ie))
    do j = is, ie
      jj = j - is
      select case(s%damp)
      case(SPECTRUM_DAMP_NONE)
        dumpa(j) = M_ONE
      case(SPECTRUM_DAMP_LORENTZIAN)
        dumpa(j)= exp(-jj*dt*s%damp_factor)
      case(SPECTRUM_DAMP_POLYNOMIAL)
        dumpa(j) = 1.0 - 3.0*(real(jj)/ntiter)**2                          &
          + 2.0*(real(jj)/ntiter)**3
      case(SPECTRUM_DAMP_GAUSSIAN)
        dumpa(j)= exp(-(jj*dt)**2*s%damp_factor**2)
      end select
    end do

    do k = 0, no_e
      w = k*s%energy_step
      do j = is, ie

        jj = j - is

        z = exp(M_zI * w * jj *dt)
        sp(k) = sp(k) + z*dumpa(j)*sum(angular(j, :)*kick%pol(1:3, kick%pol_dir))

      end do
      sp(k) = M_zI/(M_TWO*P_c*kick%delta_strength)*sp(k)*dt

      sum1 = sum1 + sp(k)*s%energy_step
      sum2 = sum2 + (sp(k)*w**2)*s%energy_step

    end do

    SAFE_DEALLOCATE_A(angular)
    SAFE_DEALLOCATE_A(dumpa)

    ! print some info
    write(message(1), '(a,i8)')    'Number of time steps = ', ntiter
    write(message(2), '(a,i4)')    'SpecDampMode         = ', s%damp
    write(message(3), '(a,f10.4)') 'SpecDampFactor       = ', units_from_atomic(units_out%time**(-1), s%damp_factor)
    write(message(4), '(a,f10.4)') 'SpecStartTime        = ', units_from_atomic(units_out%time, s%start_time)
    write(message(5), '(a,f10.4)') 'SpecEndTime          = ', units_from_atomic(units_out%time, s%end_time)
    write(message(6), '(a,f10.4)') 'SpecMaxEnergy        = ', units_from_atomic(units_inp%energy, s%max_energy) 
    write(message(7),'(a,f10.4)')  'SpecEnergyStep       = ', units_from_atomic(units_inp%energy, s%energy_step)
    message(8) = ""
    write(message(9), '(a,5e15.6,5e15.6)') 'R(0) sum rule = ', sum1
    write(message(10),'(a,5e15.6,5e15.6)') 'R(2) sum rule = ', sum2
    call write_info(10)


    ! Output to file
    write(out_file, '(a15,i2)')      '# nspin        ', nspin
    call kick_write(kick, out_file)
    write(out_file, '(a1,a20,a20,a20)') '#', str_center("Energy", 20), str_center("R", 20), str_center("Re[beta]", 20)
    write(out_file, '(a1,a20,a20,a20)') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
         str_center('['//trim(units_abbrev(units_out%length)) //'^3]', 20), &
         str_center('['//trim(units_abbrev(units_out%length)) //'^4]', 20)
    do i = 0, no_e
      write(out_file,'(e20.8,e20.8,e20.8)') units_from_atomic(units_out%energy, i*s%energy_step), &
        units_from_atomic(units_out%length**3, aimag(sp(i))/M_PI), &
        units_from_atomic(units_out%length**4, real(sp(i))*P_C/(M_THREE*max(i,1)*s%energy_step))
    end do

    call pop_sub()
  end subroutine spectrum_rotatory_strength
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_init(dt, is, ie, niter, acc)
    FLOAT,   intent(in) :: dt
    integer, intent(in) :: is, ie, niter
    CMPLX,   intent(in) :: acc(:)

    call push_sub('spectrum.spectrum_hsfunction_init')

    is_ = is
    ie_ = ie
    time_step_ = dt
    SAFE_ALLOCATE(func_(0:niter))
    func_ = acc

    call pop_sub()
  end subroutine spectrum_hsfunction_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_end
    call push_sub('spectrum.spectrum_hsfunction_end')
    SAFE_DEALLOCATE_A(func_)
    call pop_sub()
  end subroutine spectrum_hsfunction_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_min(a, b, omega, omega_min, func_min)
    FLOAT, intent(in)  :: a, b, omega
    FLOAT, intent(out) :: omega_min, func_min

    integer :: nfreqs, ierr, i
    FLOAT :: x, hsval, minhsval, dw, w, hsa, hsb

    call push_sub('spectrum.spectrum_hsfunction_min')

    ! x should be an initial guess for the minimum. So we do a quick search
    ! that we refine later calling 1dminimize.
    x = omega
    call hsfunction(x, minhsval)
    nfreqs = 100
    dw = (b-a)/(nfreqs-1)
    do i = 1, nfreqs
      w = a + dw*(i-1)
      call hsfunction(w, hsval)
      if(i .eq. 1)      hsa = hsval
      if(i .eq. nfreqs) hsb = hsval
      if(hsval < minhsval) then
        minhsval = hsval
        x = w
      end if
    end do

    if( hsa == minhsval ) then
      omega_min = a
      func_min = hsval
      call pop_sub(); return
    end if
    if( hsb == minhsval ) then
      omega_min = b
      func_min = hsval
      call pop_sub(); return
    end if

    ! Around x, we call some GSL sophisticated search algorithm to find the minimum.
#ifndef SINGLE_PRECISION
    call loct_1dminimize(a, b, x, hsfunction, ierr)
#else
    stop "FIXME: cannot work in single-precision."
#endif

    if(ierr .ne. 0) then
      write(message(1),'(a)') 'Could not find a maximum.'      
      call write_fatal(1)
    end if
    call hsfunction(x, hsval)
    omega_min = x
    func_min  = hsval

    call pop_sub()
  end subroutine spectrum_hsfunction_min
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine hsfunction(omega, power)
    FLOAT, intent(in)   :: omega
    FLOAT, intent(out)  :: power

    CMPLX   :: c, ez1, ez, z
    integer :: j

    call push_sub('spectrum.hsfunction')

    c = M_z0
    z = M_zI * omega * time_step_
    ez1 = exp((is_-1)*z)
    ez  = exp(z)
    do j = is_, ie_
      ! This would be easier, but slower.
      !c = c + exp(M_zI * omega * j * time_step_)*func_(j)
      ez1 = ez1 * ez
      c = c + ez1*func_(j)
    end do
    power = -abs(c)**2*time_step_**2

    call pop_sub()
  end subroutine hsfunction
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_from_mult(out_file, s, pol, w0)
    character(len=*),   intent(in) :: out_file
    type(spec_t),    intent(inout) :: s
    character , intent(in) :: pol
    FLOAT, optional, intent(in) :: w0

    integer :: i, j, iunit, nspin, time_steps, is, ie, ntiter, lmax
    FLOAT :: dt, dump
    type(kick_t) :: kick
    FLOAT, allocatable :: d(:,:)
    CMPLX, allocatable :: dipole(:), ddipole(:)
    type(unit_system_t) :: file_units

    call push_sub('spectrum.spectrum_hs_from_mult')

    call io_assign(iunit)
    iunit = io_open('multipoles', action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/multipoles', action='read', status='old')
    end if
    call spectrum_mult_info(iunit, nspin, kick, time_steps, dt, file_units, lmax=lmax)
    call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

    call io_skip_header(iunit)

    ! load dipole from file
    SAFE_ALLOCATE(dipole(0:time_steps))
    SAFE_ALLOCATE(ddipole(0:time_steps))
    SAFE_ALLOCATE(d(1:3, 1:nspin))

    do i = 1, time_steps
      read(iunit, *) j, dump, dump, d
      select case(pol)
      case('x')
        dipole(i) = -sum(d(1, :))
      case('y')
        dipole(i) = -sum(d(2, :))
      case('z')
        dipole(i) =  sum(d(3, :))
      case('+')
        dipole(i) = -sum(d(1, :) + M_zI*d(2, :)) / sqrt(M_TWO)
      case('-')
        dipole(i) = -sum(d(1, :) - M_zI*d(2, :)) / sqrt(M_TWO)
      end select
      dipole(i) = units_to_atomic(units_out%length, dipole(i))
    end do
    SAFE_DEALLOCATE_A(d)
    dipole(0) = dipole(1)
    call io_close(iunit)

    ! we now calculate the acceleration.
    ddipole(0) = M_ZERO
    do i = 1, time_steps - 1
      ddipole(i) = (dipole(i-1)+dipole(i+1)-M_TWO*dipole(i))/dt**2
    end do
    call interpolate( dt*(/ -3, -2, -1 /),   &
                      ddipole(time_steps-3:time_steps-1), &
                      M_ZERO, &
                      ddipole(time_steps) )

    call spectrum_hsfunction_init(dt, is, ie, time_steps, ddipole)
    call spectrum_hs(out_file, s, pol, w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(dipole)
    SAFE_DEALLOCATE_A(ddipole)

    call pop_sub()
  end subroutine spectrum_hs_from_mult
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_from_acc(out_file, s, pol, w0)
    character(len=*),   intent(in) :: out_file
    type(spec_t),    intent(inout) :: s
    character, intent(in) :: pol
    FLOAT, optional, intent(in) :: w0

    integer :: i, j, iunit, time_steps, is, ie, ntiter, ierr
    FLOAT :: dt, a(MAX_DIM)
    CMPLX, allocatable :: acc(:)

    call push_sub('spectrum.spectrum_hs_from_acc')

    call spectrum_acc_info(iunit, time_steps, dt)
    call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

    ! load dipole from file
    SAFE_ALLOCATE(acc(0:time_steps))
    acc = M_ZERO
    call io_skip_header(iunit)
    do i = 1, time_steps
      a = M_ZERO
      read(iunit, '(28x,e20.12)', advance = 'no', iostat = ierr) a(1)
      j = 2
      do while( (ierr.eq.0) .and. (j <= MAX_DIM) )
        read(iunit, '(e20.12)', advance = 'no', iostat = ierr) a(j)
      end do
      select case(pol)
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
      acc(i) = units_to_atomic(units_out%acceleration, acc(i))
    end do
    close(iunit)

    call spectrum_hsfunction_init(dt, is, ie, time_steps, acc)
    call spectrum_hs(out_file, s, pol, w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(acc)
    call pop_sub()
  end subroutine spectrum_hs_from_acc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs(out_file, s, pol, w0)
    character(len=*),   intent(in) :: out_file
    type(spec_t),    intent(inout) :: s
    character, intent(in) :: pol
    FLOAT, optional, intent(in) :: w0

    integer :: iunit, no_e, i
    FLOAT   :: omega, hsval, x
    FLOAT, allocatable :: sp(:)

    call push_sub('spectrum.spectrum_hs')

    if(present(w0)) then

      iunit = io_open(trim(out_file) // "." // trim(pol), action='write')
      write(iunit, '(a1,a20,a20)') '#', str_center("w", 20), str_center("H(w)", 20)
      write(iunit, '(a1,a20,a20)') '#', &
        str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
        str_center('[('//trim(units_abbrev(units_out%length))//'/' &
        //trim(units_abbrev(units_out%time))//')^2]' , 20)

      ! output
      omega = w0
      do while(omega <= s%max_energy)
        call spectrum_hsfunction_min(omega-w0, omega+w0, omega, x, hsval)

        write(iunit, '(1x,2e20.8)') units_from_atomic(units_out%energy, x), &
          units_from_atomic((units_out%length / units_out%time)**2, -hsval)

        ! 2*w0 because we assume that there are only odd peaks.
        omega = omega + 2*w0
      end do
      call io_close(iunit)

    else

      no_e = s%max_energy / s%energy_step
      SAFE_ALLOCATE(sp(0:no_e))
      sp = M_ZERO

      do i = 0, no_e
        call hsfunction(i*s%energy_step, sp(i))
        sp(i) = -sp(i)
      end do

      ! output
      if(trim(out_file) .ne. '-') then
        iunit = io_open(trim(out_file) // "." // trim(pol), action='write')
        write(iunit, '(a1,a20,a20)') '#', str_center("w", 20), str_center("H(w)", 20)
        write(iunit, '(a1,a20,a20)') &
          '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
          str_center('[('//trim(units_abbrev(units_out%length))//'/' &
            //trim(units_abbrev(units_out%time))//')^2]' , 20)
        do i = 0, no_e
          write(iunit, '(2e15.6)') units_from_atomic(units_out%energy, i*s%energy_step), &
            units_from_atomic((units_out%length / units_out%time)**2, sp(i))
        end do
        call io_close(iunit)
      end if
      SAFE_DEALLOCATE_A(sp)

    end if

    call pop_sub()
  end subroutine spectrum_hs
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_mult_info(iunit, nspin, kick, time_steps, dt, file_units, lmax)
    integer,           intent(in)  :: iunit
    integer,           intent(out) :: nspin
    type(kick_t),      intent(out) :: kick
    integer,           intent(out) :: time_steps
    FLOAT,             intent(out) :: dt
    type(unit_system_t), intent(out) :: file_units
    integer, optional, intent(out) :: lmax

    integer :: i
    character(len=100) :: line

    call push_sub('spectrum.spectrum_mult_info')

    rewind(iunit); read(iunit,*); read(iunit,*)
    read(iunit, '(15x,i2)')      nspin
    if(present(lmax)) then
      read(iunit, '(15x,i2)')      lmax
    end if
    call kick_read(kick, iunit)
    read(iunit, '(a)')           line
    read(iunit, '(a)')           line
    call io_skip_header(iunit)

    ! Figure out about the units of the file
    i = index(line,'eV')
    if(i.ne.0) then
      call unit_system_get(file_units, UNITS_EVA)
    else
      call unit_system_get(file_units, UNITS_ATOMIC)
    end if

    call count_time_steps(iunit, time_steps, dt)
    dt = units_to_atomic(file_units%time, dt) ! units_out is OK

    call pop_sub()
  end subroutine spectrum_mult_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------  
  subroutine count_time_steps(iunit, time_steps, dt)
    integer, intent(in)  :: iunit
    integer, intent(out) :: time_steps
    FLOAT,   intent(out) :: dt

    FLOAT :: t1, t2, dummy
    integer :: jj

    call push_sub('spectrum.count_time_steps')

    ! count number of time_steps
    time_steps = 0
    do
      read(iunit, *, end=100) jj, dummy
      time_steps = time_steps + 1
      if(time_steps == 1) t1 = dummy
      if(time_steps == 2) t2 = dummy
    end do
100 continue
    dt = (t2 - t1)
    time_steps = time_steps - 1
    
    if(time_steps < 3) then
      message(1) = "Empty file?"
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine count_time_steps
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section_info(iunit, nspin, kick, energy_steps, dw)
    integer, intent(in)           :: iunit
    integer, intent(out)          :: nspin
    type(kick_t), intent(out)     :: kick
    integer, intent(out)          :: energy_steps
    FLOAT,   intent(out)          :: dw

    FLOAT :: dummy, e1, e2

    call push_sub('spectrum.spectrum_cross_section_info')

    ! read in number of spin components
    read(iunit, '(15x,i2)')      nspin
    call kick_read(kick, iunit)
    call io_skip_header(iunit)

    ! count number of time_steps
    energy_steps = 0
    do
      read(iunit, *, end=100) dummy
      energy_steps = energy_steps + 1
      if(energy_steps == 1) e1 = dummy
      if(energy_steps == 2) e2 = dummy
    end do
100 continue
    dw = units_to_atomic(units_out%energy, e2 - e1)
    energy_steps = energy_steps - 1

    if(energy_steps < 3) then
      message(1) = "Empty multipole file?"
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine spectrum_cross_section_info


  ! ---------------------------------------------------------
  subroutine spectrum_acc_info(iunit, time_steps, dt)
    integer, intent(out) :: iunit, time_steps
    FLOAT,   intent(out) :: dt

    integer :: j
    FLOAT :: t1, t2, dummy

    call push_sub('spectrum.spectrum_acc_info')

    ! open files
    iunit = io_open('acceleration', action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/acceleration', action='read', status='old')
    end if

    ! read in dipole
    call io_skip_header(iunit)

    ! count number of time_steps
    time_steps = 0
    do
      read(iunit, *, end=100) j, dummy
      time_steps = time_steps + 1
      if(time_steps == 1) t1 = dummy
      if(time_steps == 2) t2 = dummy
    end do
100 continue
    dt = units_to_atomic(units_out%time, t2 - t1) ! units_out is OK
    time_steps = time_steps - 1

    if(time_steps < 3) then
      message(1) = "Empty multipole file?"
      call write_fatal(1)
    end if

    rewind(iunit)
    call pop_sub()
  end subroutine spectrum_acc_info


  ! ---------------------------------------------------------
  subroutine spectrum_fix_time_limits(time_steps, dt, start_time, end_time, is, ie, ntiter)
    integer,  intent(in) :: time_steps
    FLOAT,    intent(in) :: dt
    FLOAT, intent(inout) :: start_time, end_time
    integer, intent(out) :: is, ie, ntiter

    FLOAT :: ts, te, dummy

    call push_sub('spectrum.spectrum_fix_time_limits')

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

    call pop_sub()
  end subroutine spectrum_fix_time_limits

end module spectrum_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
