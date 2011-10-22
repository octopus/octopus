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
  use batch_m
  use c_pointer_m
  use compressed_sensing_m
  use datasets_m
  use global_m
  use io_m
  use lalg_adv_m
  use loct_math_m
  use math_m
  use messages_m
  use parser_m
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
    spectrum_dyn_structure_factor, &
    spectrum_rotatory_strength,    &
    spectrum_hs_from_mult,         &
    spectrum_hs_ar_from_mult,      &
    spectrum_hs_ar_from_acc,       &
    spectrum_hs_from_acc,          &
    spectrum_hs_from_vel,          &
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

  integer, public, parameter ::    &
    QKICKMODE_NONE           = 0,  &
    QKICKMODE_EXP            = 1,  &
    QKICKMODE_COS            = 2,  &
    QKICKMODE_SIN            = 3,  &
    QKICKMODE_BESSEL         = 4

  integer, public, parameter ::    &
    SPECTRUM_ABSORPTION      = 1,  &
    SPECTRUM_ENERGYLOSS      = 2

  integer, public, parameter ::       &
    SPECTRUM_FOURIER            = 1,  &
    SPECTRUM_COMPRESSED_SENSING = 2

  type spec_t
    FLOAT   :: start_time          ! start time for the transform
    FLOAT   :: end_time            ! when to stop the transform
    FLOAT   :: energy_step         ! step in energy mesh
    FLOAT   :: max_energy          ! maximum of energy mesh
    integer :: damp                ! damping type (none, exp or pol)
    integer :: transform           ! sine, cosine, or exponential transform
    FLOAT   :: damp_factor         ! factor used in damping
    integer :: spectype            ! damping type (none, exp or pol)
    integer :: method              ! fourier transform or compressed sensing 
    FLOAT   :: noise               ! the level of noise that is assumed in the time series for compressed sensing 
  end type spec_t

  type kick_t
    integer           :: delta_strength_mode
    FLOAT             :: delta_strength
    ! In case we use a normal dipole kick:
    FLOAT             :: pol(MAX_DIM, MAX_DIM)
    integer           :: pol_dir
    integer           :: pol_equiv_axes
    FLOAT             :: wprime(MAX_DIM)
    ! In case we have a general multipolar kick,
    ! the form of this "kick" will be (atomic units):
    ! V(\vec{r}) = sum_{i=1}^{n_multipoles} 
    !                 weight(i) * (e^2 / a_0^(l+1)) * r^l(i) * Y_{l(i),m(i)} (\vec{r})
    ! which has units of energy; if we include the time-dependence (delta function):
    ! V(\vec{r}) = sum_{i=1}^{n_multipoles} 
    !                 weight(i) * (\hbar / a_0^l) * r^l(i) * Y_{l(i),m(i)} (\vec{r}) * \delta(t)
    integer           :: n_multipoles
    integer, pointer  :: l(:), m(:)
    FLOAT, pointer    :: weight(:)
    FLOAT             :: qvector(MAX_DIM)
    FLOAT             :: qlength
    integer           :: qkick_mode
    integer           :: qbessel_l, qbessel_m
  end type kick_t

  ! Module variables, necessary to compute the function hsfunction, called by
  ! the C function loct_1dminimize
  FLOAT :: time_step_
  CMPLX, allocatable :: func_(:),func_ar_(:,:),pos_(:,:),tret_(:)
  CMPLX :: vv_(MAX_DIM)
  integer :: is_, ie_, default
  logical :: from_vel_

contains

  ! ---------------------------------------------------------
  subroutine spectrum_init(spectrum)
    type(spec_t), intent(inout) :: spectrum

    PUSH_SUB(spectrum_init)

    !%Variable PropagationSpectrumType
    !%Type integer
    !%Default AbsorptionSpectrum
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Type of spectrum to calculate.
    !%Option AbsorptionSpectrum 1
    !% Photoabsorption spectrum.
    !%Option EnergyLossSpectrum 2
    !% Dynamic structure factor (also known as energy-loss function or spectrum).
    !%End

    call parse_integer  (datasets_check('PropagationSpectrumType'), SPECTRUM_ABSORPTION, spectrum%spectype)
    if(.not.varinfo_valid_option('PropagationSpectrumType', spectrum%spectype)) call input_error('PropagationSpectrumType')

    !%Variable SpectrumMethod
    !%Type integer
    !%Default fourier
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Decides which method is used to obtain the spectrum. The
    !% default is the fourier transform.
    !%Option fourier 1
    !% The standard fourier transform.
    !%Option compressed_sensing 2
    !% (Experimental) Uses the compressed sensing technique.
    !%End
    call parse_integer  (datasets_check('SpectrumMethod'), SPECTRUM_FOURIER, spectrum%method)
    if(.not.varinfo_valid_option('SpectrumMethod', spectrum%method)) then
      call input_error('SpectrumMethod')
    endif

    !%Variable SpectrumSignalNoise
    !%Type float
    !%Default 3e-4 au
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% For compressed sensing, the signal to process, the
    !% time-dependent dipole in this case, is assumed to have some
    !% noise that is given by this quantity. The default value is
    !% 3.0e-4. This value is always assumed to be in atomic units.
    !%End
    call parse_float(datasets_check('SpectrumSignalNoise'), CNST(3.0e-4), spectrum%noise)

    if(spectrum%method == SPECTRUM_COMPRESSED_SENSING) then
      call messages_experimental('compressed sensing')
    end if

    !%Variable PropagationSpectrumDampMode
    !%Type integer
    !%Default polynomial
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Decides which damping/filtering is to be applied in order to
    !% calculate spectra by calculating a Fourier transform. The
    !% default is polynomial damping, except when compressed sensing
    !% is used. In that case the default is none.
    !%Option none 0
    !% No filtering at all.
    !%Option exponential 1
    !% Exponential filtering, corresponding to a Lorentzian-shaped spectrum.
    !%Option polynomial 2
    !% Third-order polynomial damping.
    !%Option gaussian 3
    !% Gaussian damping.
    !%End
    default = SPECTRUM_DAMP_POLYNOMIAL
    if(spectrum%method == SPECTRUM_COMPRESSED_SENSING) default = SPECTRUM_DAMP_NONE

    call parse_integer  (datasets_check('PropagationSpectrumDampMode'), default, spectrum%damp)
    if(.not.varinfo_valid_option('PropagationSpectrumDampMode', spectrum%damp)) call input_error('PropagationSpectrumDampMode')

    if(spectrum%method == SPECTRUM_COMPRESSED_SENSING .and. spectrum%damp /= SPECTRUM_DAMP_NONE) then
      message(1) = 'Using damping with compressed sensing, this is not required'
      message(2) = 'and can introduce noise in the spectra.'
      call messages_warning(2)
    end if

    !%Variable PropagationSpectrumTransform
    !%Type integer
    !%Default sine
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Decides which transform to perform.
    !%Option sine 2
    !% Sine transform <math>\int dt \sin(wt) f(t)</math>
    !%Option cosine 3
    !% Cosine transform <math>\int dt \cos(wt) f(t)</math>
    !%Option exponential 1
    !% Exponential transform <math>\int dt \exp(-wt) f(t)</math>
    !%End
    call parse_integer  (datasets_check('PropagationSpectrumTransform'), SPECTRUM_TRANSFORM_SIN, spectrum%transform)
    if(.not.varinfo_valid_option('PropagationSpectrumTransform', spectrum%transform)) then
      call input_error('PropagationSpectrumTransform')
    endif

    !%Variable PropagationSpectrumStartTime
    !%Type float
    !%Default 0.0
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Processing is done for the given function in a time-window that starts at the
    !% value of this variable.
    !%End
    call parse_float(datasets_check('PropagationSpectrumStartTime'),  M_ZERO, spectrum%start_time, units_inp%time)

    !%Variable PropagationSpectrumEndTime
    !%Type float
    !%Default -1.0 au
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Processing is done for the given function in a time-window that ends at the
    !% value of this variable. If set to a negative value, the maximum value from 
    !% the corresponding multipole file will used.
    !%End
    call parse_float(datasets_check('PropagationSpectrumEndTime'), -M_ONE, spectrum%end_time, units_inp%time)

    !%Variable PropagationSpectrumEnergyStep
    !%Type float
    !%Default 0.01 eV
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Sampling rate for the spectrum. If you supply a number equal or smaller than zero, then
    !% the sampling rate will be (2 * pi / T), where T is the total propagation time.
    !%End
    call parse_float(datasets_check('PropagationSpectrumEnergyStep'), CNST(0.01) / (M_TWO * P_Ry), &
      spectrum%energy_step, units_inp%energy)
    

    !%Variable PropagationSpectrumMaxEnergy
    !%Type float
    !%Default 20 eV
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% The Fourier transform is calculated for energies smaller than this value.
    !%End
    call parse_float(datasets_check('PropagationSpectrumMaxEnergy'), CNST(20.0) / (M_TWO * P_Ry), &
      spectrum%max_energy, units_inp%energy)

    !%Variable PropagationSpectrumDampFactor
    !%Type float
    !%Default 0.15 au
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% If <tt>PropagationSpectrumDampMode = exponential</tt>, the damping parameter of the exponential
    !% is fixed through this variable.
    !%End
    call parse_float(datasets_check('PropagationSpectrumDampFactor'), CNST(0.15), &
      spectrum%damp_factor, units_inp%time**(-1))

    POP_SUB(spectrum_init)
  end subroutine spectrum_init


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
      if(kick%n_multipoles > 0) then
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

    else if(present(out)) then
      write(aux, '(a15,i2)')      '# kick mode    ', kick%delta_strength_mode
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      write(aux, '(a15,f18.12)')  '# kick strength', kick%delta_strength
      call write_iter_string(out, aux)
      call write_iter_nl(out)
      if(kick%n_multipoles > 0) then
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
    end if

    POP_SUB(kick_write)
  end subroutine kick_write


  ! ---------------------------------------------------------
  subroutine kick_read(kick, iunit)
    type(kick_t), intent(inout) :: kick
    integer,      intent(in)    :: iunit

    integer :: im
    character(len=100) :: line

    PUSH_SUB(kick_read)

    read(iunit, '(15x,i2)')     kick%delta_strength_mode
    read(iunit, '(15x,f18.12)') kick%delta_strength
    read(iunit, '(a)') line
    if(index(line,'multipole').ne.0) then
      read(line, '("# N multipoles ",i3)') kick%n_multipoles
      SAFE_ALLOCATE(     kick%l(1:kick%n_multipoles))
      SAFE_ALLOCATE(     kick%m(1:kick%n_multipoles))
      SAFE_ALLOCATE(kick%weight(1:kick%n_multipoles))
      do im = 1, kick%n_multipoles
        read(iunit, '("# multipole    ",2i3,f18.12)') kick%l(im), kick%m(im), kick%weight(im)
      end do
    else
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

    POP_SUB(kick_read)
  end subroutine kick_read


  ! ---------------------------------------------------------
  subroutine kick_init(kick, nspin, dim)
    type(kick_t), intent(out) :: kick
    integer,      intent(in)  :: nspin
    integer,      intent(in)  :: dim

    type(block_t) :: blk
    integer :: n_rows, irow, idir

    PUSH_SUB(kick_init)

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


    !%Variable TDKickFunction
    !%Type block
    !%Section Time-Dependent::Response
    !%Description
    !% If the block <tt>TDKickFunction</tt> is present in the input file, the kick function to
    !% be applied at time zero of the time-propagation will not be a "dipole" function
    !% (<i>i.e.</i> phi => exp(i*k*z) phi), but a general multipole in the form r^l * Y_{lm}(r).
    !%
    !% Each line has two columns of integers: the (<i>l</i>,<i>m</i>) pair that defines the
    !% multipole. Any number of lines may be given, and the kick will be the sum of those
    !% multipoles.
    !%
    !% This feature allows calculation of quadrupole, octupole, etc., response functions.
    !%End
    if(parse_block(datasets_check('TDKickFunction'), blk) == 0) then
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
      kick%n_multipoles = 0
      nullify(kick%l)
      nullify(kick%m)
      nullify(kick%weight)
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
    !% <i>et al.</i>, <i>J. Nanoscience and Nanotechnology<i> <b>8</b>,
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
    !% an optional last number. Possible options are <tt>qexp</tt> (default), </tt>qcos</tt>,
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
          kick%qbessel_l = M_ZERO
          kick%qbessel_m = M_ZERO
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
  subroutine spectrum_cross_section_tensor(spectrum, out_file, in_file)
    type(spec_t), intent(inout) :: spectrum
    integer,      intent(in)    :: out_file
    integer,      intent(in)    :: in_file(:)

    character(len=20) :: header_string
    integer :: nspin, energy_steps, ie, is, equiv_axes, n_files, idir, jdir, ii, trash
    FLOAT, allocatable :: sigma(:, :, :, :), sigmap(:, :, :, :), sigmau(:, :, :),  &
      sigmav(:, :, :), sigmaw(:, :, :), pp(:, :), ip(:, :)
    FLOAT :: dw, dump, average, anisotropy
    type(kick_t) :: kick

    PUSH_SUB(spectrum_cross_section_tensor)

    n_files = size(in_file)
    equiv_axes = 3 - n_files + 1

    call spectrum_cross_section_info(in_file(1), nspin, kick, energy_steps, dw)
    ! on subsequent calls, do not overwrite energy_steps and dw
    call io_skip_header(in_file(1))

    SAFE_ALLOCATE(sigma (1:3, 1:3, 0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmap(1:3, 1:3, 0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmau(1:3,      0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmav(1:3,      0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmaw(1:3,      0:energy_steps, 1:nspin))
    SAFE_ALLOCATE(    pp(1:3, 1:3))
    SAFE_ALLOCATE(    ip(1:3, 1:3))

    select case(equiv_axes)

    case(3)

      do ie = 0, energy_steps
        read(in_file(1), *) dump, sigmau(1:3, ie, 1:nspin)
      end do

      ! The first row of sigma is the vector that we have just read, but properly projected...
      do is = 1, nspin
        do ie = 0, energy_steps
          sigmap(1, 1, ie, is) = sum(sigmau(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(1, 2, ie, is) = sum(sigmau(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(1, 3, ie, is) = sum(sigmau(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do

      ! The diagonal parts are also equal:
      sigmap(2, 2, :, :) = sigmap(1, 1, :, :)
      sigmap(3, 3, :, :) = sigmap(1, 1, :, :)

      ! The (2,1) term and (3,1) term are equal by symmetry:
      sigmap(2, 1, :, :) = sigmap(1, 2, :, :)
      sigmap(3, 1, :, :) = sigmap(1, 3, :, :)

      ! But for the (2,3) term we need the wprime vector....
      do is = 1, nspin
        do ie = 0, energy_steps
          sigmap(2, 3, ie, is) = sum(sigmau(1:3, ie, is) * kick%wprime(1:3))
          sigmap(3, 2, ie, is) = sigmap(2, 3, ie, is)
        end do
      end do

    case(2)

      call spectrum_cross_section_info(in_file(2), ie, kick, trash, dump)
      call io_skip_header(in_file(2))

      do ie = 0, energy_steps
        read(in_file(1), *) dump, sigmau(1:3, ie, 1:nspin)
        read(in_file(2), *) dump, sigmaw(1:3, ie, 1:nspin)
      end do

      ! The first row of sigma is the vector that we have just read, but properly projected...
      do is = 1, nspin
        do ie = 0, energy_steps
          sigmap(1, 1, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(1, 2, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(1, 3, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do

      ! The third row of sigma is also the vector that we have just read, but properly projected...
      do is = 1, nspin
        do ie = 0, energy_steps
          sigmap(3, 1, ie, is) = sum( sigmaw(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(3, 2, ie, is) = sum( sigmaw(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(3, 3, ie, is) = sum( sigmaw(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do

      ! The diagonal (2,2) is equal by symmetry to (1,1)
      sigmap(2, 2, :, :) = sigmap(1, 1, :, :)

      ! The (2,1) term and (1,2) term are equal; (2,3) and (3,2), also.
      sigmap(2, 1, :, :) = sigmap(1, 2, :, :)
      sigmap(2, 3, :, :) = sigmap(3, 2, :, :)

    case default

      call spectrum_cross_section_info(in_file(2), ie, kick, trash, dump)
      call spectrum_cross_section_info(in_file(3), ie, kick, trash, dump)
      call io_skip_header(in_file(2))
      call io_skip_header(in_file(3))

      do ie = 0, energy_steps
        read(in_file(1), *) dump, sigmau(1:3, ie, 1:nspin)
        read(in_file(2), *) dump, sigmav(1:3, ie, 1:nspin)
        read(in_file(3), *) dump, sigmaw(1:3, ie, 1:nspin)
      end do

      do is = 1, nspin
        do ie = 0, energy_steps
          sigmap(1, 1, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(1, 2, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(1, 3, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do
      do is = 1, nspin
        do ie = 0, energy_steps
          sigmap(2, 1, ie, is) = sum( sigmav(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(2, 2, ie, is) = sum( sigmav(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(2, 3, ie, is) = sum( sigmav(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do
      do is = 1, nspin
        do ie = 0, energy_steps
          sigmap(3, 1, ie, is) = sum( sigmaw(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(3, 2, ie, is) = sum( sigmaw(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(3, 3, ie, is) = sum( sigmaw(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do

    end select

    ! And now, perform the necessary transformation.
    ip(1:3, 1:3) = kick%pol(1:3, 1:3)
    dump = lalg_inverter(3, ip)
    do is = 1, nspin
      do ie = 0, energy_steps
        sigma(:, :, ie, is) = matmul( transpose(ip), matmul(sigmap(:, :, ie, is), ip) )
      end do
    end do

    ! Finally, write down the result
    write(out_file, '(a15,i2)')      '# nspin        ', nspin
    call kick_write(kick, out_file)
    write(out_file, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
    write(out_file, '(a20)', advance = 'no') str_center("(1/3)*Tr[sigma]", 20)
    write(out_file, '(a20)', advance = 'no') str_center("Anisotropy[sigma]", 20)
    do is = 1, nspin
      do idir = 1, 3
        do jdir = 1, 3
          write(header_string,'(a6,i1,a1,i1,a1,i1,a1)') 'sigma(', idir, ',', jdir, ',', is, ')'
          write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
        end do
      end do
    end do
    write(out_file, '(1x)')
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center('[' // trim(units_abbrev(units_out%energy)) // ']', 20)
    do ii = 1, 2 + nspin * 9
      write(out_file, '(a20)', advance = 'no')  str_center('[' // trim(units_abbrev(units_out%length**2)) // ']', 20)
    end do
    write(out_file, '(1x)')

    ! The anisotropy (Delta alpha) of a second-rank symmetric tensor alpha (such 
    ! as the cross section) is defined as:
    ! 
    ! (Delta alpha)^2 = (1/3) * ( 3 * Tr(alpha^2) - (Tr(a))^2 )
    !
    ! The reason for this definition is that it is identically equal to:
    !
    ! (Delta alpha)^2 = (1/3) * ( (alpha_1-alpha_2)^2 + (alpha_1-alpha_3)^2 + (alpha_2-alpha_3)^2 )
    !
    ! where {alpha_1, alpha_2, alpha_3} are the eigenvalues of alpha. An "isotropic" tensor
    ! is characterized by having three equal eigenvalues, which leads to zero anisotropy. The
    ! more different that the eigenvalues are, the larger the anisotropy is.
    do ie = 0, energy_steps

      pp = M_ZERO
      do is = 1, min(2, nspin) ! we add spin up with spin down
        pp(:, :) = pp(:, :) + sigma(:, :, ie, is)
      end do
      average = M_THIRD * ( pp(1, 1) + pp(2, 2) + pp(3, 3) )
      ip = matmul(pp, pp)
      anisotropy = M_THIRD * ( M_THREE * ( ip(1, 1) + ip(2, 2) + ip(3, 3) ) - (M_THREE * average)**2 )

      ! Note that the cross-section elements do not have to be transformed to the proper units, since
      ! they have been read from the "cross_section_vector.x", where they are already in the proper units.
      write(out_file,'(3e20.8)', advance = 'no') units_from_atomic(units_out%energy, (ie * spectrum%energy_step)), &
        average , sqrt(max(anisotropy, M_ZERO)) 
      do is = 1, nspin
        write(out_file,'(9e20.8)', advance = 'no') sigma(1:3, 1:3, ie, is)
      end do
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(sigma)
    SAFE_DEALLOCATE_A(sigmap)
    SAFE_DEALLOCATE_A(sigmau)
    SAFE_DEALLOCATE_A(sigmav)
    SAFE_DEALLOCATE_A(sigmaw)
    SAFE_DEALLOCATE_A(pp)
    SAFE_DEALLOCATE_A(ip)
    POP_SUB(spectrum_cross_section_tensor)
  end subroutine spectrum_cross_section_tensor


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section(in_file, out_file, spectrum)
    integer,      intent(in)    :: in_file
    integer,      intent(in)    :: out_file
    type(spec_t), intent(inout) :: spectrum

    character(len=20) :: header_string
    integer :: nspin, lmax, time_steps, istart, iend, ntiter, it, ii, isp, no_e, ie, idir, trash
    FLOAT   :: dt, dump, energy, ewsum, polsum
    type(kick_t) :: kick
    FLOAT, allocatable :: dipole(:, :, :), sigma(:, :, :), sf(:, :)
    type(unit_system_t) :: file_units
    type(batch_t) :: dipoleb, sigmab

    PUSH_SUB(spectrum_cross_section)

    ! This function gives us back the unit connected to the "multipoles" file, the header information,
    ! the number of time steps, and the time step.
    call spectrum_mult_info(in_file, nspin, kick, time_steps, dt, file_units, lmax=lmax)

    ! Now we cannot process files that do not contain the dipole, or that contain more than the dipole.
    if(lmax.ne.1) then
      message(1) = 'Multipoles file should contain the dipole -- and only the dipole.'
      call messages_fatal(1)
    end if

    ! Find out the iteration numbers corresponding to the time limits.
    call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

    ! Read the dipole.
    call io_skip_header(in_file)

    SAFE_ALLOCATE(dipole(0:time_steps, 1:3, 1:nspin))

    do it = 0, time_steps
      select case(nspin)
      case(1)
        read(in_file, *) trash, dump, dump, dipole(it, 1:3, 1)
      case(2)
        read(in_file, *) trash, dump, dump, dipole(it, 1:3, 1), dump, dipole(it, 1:3, 2)
      case(4)
        read(in_file, *) &
          trash, dump, dump, dipole(it, 1:3, 1), dump, dipole(it, 1:3, 2), dump, dipole(it, 1:3, 3), dump, dipole(it, 1:3, 4)
      end select

      dipole(it, 1:3, :) = units_to_atomic(file_units%length, dipole(it, 1:3, :))
      
    end do

    ! Now subtract the initial dipole.
    do it = time_steps, 0, -1
      dipole(it, :, :) = dipole(it, :, :) - dipole(0, :, :)
    end do

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! Get the number of energy steps.
    no_e = spectrum%max_energy / spectrum%energy_step
    SAFE_ALLOCATE(sigma(0:no_e, 1:3, 1:nspin))

    call batch_init(dipoleb, 3, 1, nspin, dipole)
    call batch_init(sigmab, 3, 1, nspin, sigma)

    call signal_damp(spectrum%damp, spectrum%damp_factor, istart + 1, iend + 1, dt, dipoleb)
    call fourier_transform(spectrum%method, spectrum%transform, spectrum%noise, &
      istart + 1, iend + 1, dt, dipoleb, 1, no_e + 1, spectrum%energy_step, sigmab)

    call batch_end(dipoleb)
    call batch_end(sigmab)

    SAFE_DEALLOCATE_A(dipole)

    SAFE_ALLOCATE(sf(0:no_e, nspin))

    do ie = 0, no_e
      energy = ie * spectrum%energy_step
      forall(isp = 1:nspin) sf(ie, isp) = sum(sigma(ie, 1:3, isp)*kick%pol(1:3, kick%pol_dir))
      sf(ie, 1:nspin) = -sf(ie, 1:nspin) * (energy * M_TWO) / (M_PI * kick%delta_strength)
      sigma(ie, 1:3, 1:nspin) = -sigma(ie, 1:3, 1:nspin)*(M_FOUR*M_PI*energy/P_c)/kick%delta_strength
    end do
    
    ewsum = sum(sf(0, 1:nspin))
    polsum = M_ZERO

    do ie = 1, no_e
      energy = ie * spectrum%energy_step
      ewsum = ewsum + sum(sf(ie, 1:nspin))
      polsum = polsum + sum(sf(ie, 1:nspin)) / energy**2
    end do

    ewsum = ewsum * spectrum%energy_step
    polsum = polsum * spectrum%energy_step

    write(out_file, '(a15,i2)')      '# nspin        ', nspin
    call kick_write(kick, out_file)
    write(out_file, '(a)') '#%'
    write(out_file, '(a,i8)')    '# Number of time steps = ', time_steps
    write(out_file, '(a,i4)')    '# PropagationSpectrumDampMode   = ', spectrum%damp
    write(out_file, '(a,f10.4)') '# PropagationSpectrumDampFactor = ', units_from_atomic(units_out%time**(-1), &
                                                                       spectrum%damp_factor)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumStartTime  = ', units_from_atomic(units_out%time, spectrum%start_time)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumEndTime    = ', units_from_atomic(units_out%time, spectrum%end_time)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumMaxEnergy  = ', units_from_atomic(units_out%energy, spectrum%max_energy)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumEnergyStep = ', units_from_atomic(units_out%energy, spectrum%energy_step)
    write(out_file, '(a)') '#%'
    write(out_file, '(a,f16.6)') '# Electronic sum rule       = ', ewsum
    write(out_file, '(a,f16.6)') '# Polarizability (sum rule) = ', units_from_atomic(units_out%length**3, polsum)
    write(out_file, '(a)') '#%'
    
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center("Energy", 20)
    do isp = 1, nspin
      do idir = 1, 3
        write(header_string,'(a6,i1,a8,i1,a1)') 'sigma(', idir, ', nspin=', isp, ')'
        write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
      end do
    end do
    do isp = 1, nspin
      write(header_string,'(a18,i1,a1)') 'StrengthFunction(', isp, ')'
      write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
    end do
    write(out_file, '(1x)')
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20)
    do ii = 1, nspin * 3
      write(out_file, '(a20)', advance = 'no') str_center('[' // trim(units_abbrev(units_out%length**2)) // ']', 20)
    end do
    do isp = 1, nspin
      write(out_file, '(a20)', advance = 'no') str_center('[/' // trim(units_abbrev(units_out%energy)) // ']', 20)
    end do
    write(out_file, '(1x)')

    do ie = 0, no_e
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, ie * spectrum%energy_step)
      do isp = 1, nspin
        write(out_file,'(3e20.8)', advance = 'no') (units_from_atomic(units_out%length**2, sigma(ie, idir, isp)), &
                                                    idir = 1, 3)
      end do
      do isp = 1, nspin
        write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, sf(ie, isp))
      end do
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(sigma)
    POP_SUB(spectrum_cross_section)
  end subroutine spectrum_cross_section

  ! ---------------------------------------------------------
  subroutine spectrum_dyn_structure_factor(in_file_sin, in_file_cos, out_file, spectrum)
    integer,      intent(in)    :: in_file_sin, in_file_cos
    integer,      intent(in)    :: out_file
    type(spec_t), intent(inout) :: spectrum

    character(len=20) :: header_string
    integer :: time_steps, time_steps_sin, time_steps_cos
    integer :: istart, iend, ntiter, it, jj, ii, no_e, ie, trash
    FLOAT   :: dt, dt_sin, dt_cos
    FLOAT   :: dump, dummy1, dummy2, dummy3, dummy4, energy, fsum
    type(kick_t) :: kick
    CMPLX :: xx
    CMPLX, allocatable :: ftchd(:), chi(:), damp(:)
    type(unit_system_t) :: file_units
    character(len=100) :: line

    PUSH_SUB(spectrum_dyn_structure_factor)

    ! Read information from ftchds.sin file

    rewind(in_file_sin)

    ! skip two lines 
    read(in_file_sin, *)
    read(in_file_sin, *)
    read(in_file_sin, '(15x,i2)')     kick%qkick_mode
    read(in_file_sin, '(10x,3f9.5)')  kick%qvector
    read(in_file_sin, '(15x,f18.12)') kick%delta_strength

    ! skip two lines 
    read(in_file_sin, *)
    read(in_file_sin, '(a)') line
    call io_skip_header(in_file_sin)
    call io_skip_header(in_file_cos)

    ! Figure out the units of the file
    ii = index(line, 'eV')
    if(ii.ne.0) then
      call unit_system_get(file_units, UNITS_EVA)
    else
      call unit_system_get(file_units, UNITS_ATOMIC)
    end if

    ! get time_steps and dt, and make sure that dt is the same in the two files
    call count_time_steps(in_file_sin, time_steps_sin, dt_sin)
    call count_time_steps(in_file_cos, time_steps_cos, dt_cos)

    if(dt_sin.ne.dt_cos) then
      message(1) = "dt is different in ftchds.cos and ftchds.sin!"
      call messages_fatal(1)
    end if

    time_steps = min(time_steps_sin, time_steps_cos)
    dt = units_to_atomic(file_units%time, dt_cos) ! units_out is OK

    ! Find out the iteration numbers corresponding to the time limits.
    call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

    ! Read the f-transformed charge density.
    call io_skip_header(in_file_sin)
    call io_skip_header(in_file_cos)

    SAFE_ALLOCATE(ftchd(0:time_steps))
    do it = 0, time_steps
      read(in_file_sin, *) trash, dump, dummy1, dummy2
      read(in_file_cos, *) trash, dump, dummy3, dummy4
      ftchd(it) = cmplx(dummy3-dummy2, dummy4+dummy1)
    end do

    ! Now subtract the initial value.
    do it = time_steps, 0, -1
      ftchd(it) = ftchd(it) - ftchd(0)
    end do

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! Get the number of energy steps.
    no_e = spectrum%max_energy / spectrum%energy_step

    SAFE_ALLOCATE(chi(0:no_e))
    chi = M_ZERO

    ! Gets the damp function
    SAFE_ALLOCATE(damp(istart:iend))
    do it = istart, iend
      jj = it - istart
      select case(spectrum%damp)
      case(SPECTRUM_DAMP_NONE)
        damp(it) = M_ONE
      case(SPECTRUM_DAMP_LORENTZIAN)
        damp(it)= exp(-jj * dt * spectrum%damp_factor)
      case(SPECTRUM_DAMP_POLYNOMIAL)
        damp(it) = M_ONE - M_THREE * (real(jj) / ntiter)**2 &
          + M_TWO * (real(jj) / ntiter)**3
      case(SPECTRUM_DAMP_GAUSSIAN)
        damp(it)= exp(-(jj * dt)**2 * spectrum%damp_factor**2)
      end select
    end do

    ! Fourier transformation from time to frequency
    do ie = 0, no_e
      energy = ie * spectrum%energy_step
      do it = istart, iend
        jj = it - istart

        xx = exp(M_zI * energy * jj * dt)
        chi(ie) = chi(ie) + xx * damp(it) * ftchd(it)

      end do
      chi(ie) = chi(ie) * dt / kick%delta_strength / M_PI
    end do

    ! Test f-sum rule
    fsum = M_ZERO
    do ie = 0, no_e
      energy = ie * spectrum%energy_step
      fsum = fsum + energy * aimag(chi(ie))
    end do
    fsum = spectrum%energy_step * fsum * 2/sum(kick%qvector(:)**2)

    write(out_file, '(a)') '#%'
    write(out_file, '(a,i8)')    '# Number of time steps = ', time_steps
    write(out_file, '(a,i4)')    '# PropagationSpectrumDampMode   = ', spectrum%damp
    write(out_file, '(a,f10.4)') '# PropagationSpectrumDampFactor = ', units_from_atomic(units_out%time**(-1), &
                                                                       spectrum%damp_factor)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumStartTime  = ', units_from_atomic(units_out%time, spectrum%start_time)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumEndTime    = ', units_from_atomic(units_out%time, spectrum%end_time)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumMaxEnergy  = ', units_from_atomic(units_out%energy, spectrum%max_energy)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumEnergyStep = ', units_from_atomic(units_out%energy, spectrum%energy_step)
    write(out_file, '(a,3f9.5)') '# qvector    : ', kick%qvector
    write(out_file, '(a,f10.4)') '# F-sum rule : ', fsum
    write(out_file, '(a)') '#%'

    write(out_file, '(a1,a20)', advance = 'no') '#', str_center("Energy", 20)
    write(header_string,'(a3)') 'chi'
    write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
    write(out_file, '(1x)')
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20)
    write(out_file, '(a20)', advance = 'no') str_center('[' // trim(units_abbrev(units_out%energy)) // '**(-1)]', 20)
    write(out_file, '(1x)')

    do ie = 0, no_e
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, ie * spectrum%energy_step)
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy**(-1), aimag(chi(ie)))
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(ftchd)
    SAFE_DEALLOCATE_A(chi)
    POP_SUB(spectrum_dyn_structure_factor)

  end subroutine spectrum_dyn_structure_factor


  ! ---------------------------------------------------------
  subroutine spectrum_rotatory_strength(in_file, out_file, spectrum)
    integer,      intent(in)    :: in_file
    integer,      intent(in)    :: out_file
    type(spec_t), intent(inout) :: spectrum

    integer :: istart, iend, ntiter, ie, idir, time_steps, no_e, nspin, trash, it
    FLOAT :: dump, dt, energy
    type(kick_t) :: kick
    CMPLX :: sum1, sum2, sp
    FLOAT, allocatable :: angular(:, :), resp(:), imsp(:)
    type(batch_t) :: angularb, respb, imspb
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_rotatory_strength)

    call spectrum_mult_info(in_file, nspin, kick, time_steps, dt, file_units)
    call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

    ! load angular momentum from file
    SAFE_ALLOCATE(angular(0:time_steps, 1:3))
    call io_skip_header(in_file)
    do ie = 0, time_steps
      read(in_file, *) trash, dump, angular(ie, 1:3)
    end do

    ! subtract static dipole
    do idir = 1, 3
      angular(:, idir) = angular(:, idir) - angular(0, idir)
    end do

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    no_e = spectrum%max_energy / spectrum%energy_step

    do it = istart, iend
      angular(it, 1) = sum(angular(it, 1:3)*kick%pol(1:3, kick%pol_dir))
    end do

    SAFE_ALLOCATE(resp(0:no_e))
    SAFE_ALLOCATE(imsp(0:no_e))

    call batch_init(angularb, 1)
    call batch_init(respb, 1)
    call batch_init(imspb, 1)

    call batch_add_state(angularb, angular(:, 1))
    call batch_add_state(respb, resp)
    call batch_add_state(imspb, imsp)

    call signal_damp(spectrum%damp, spectrum%damp_factor, istart + 1, iend + 1, dt, angularb)

    call fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
      istart + 1, iend + 1, dt, angularb, 1, no_e + 1, spectrum%energy_step, respb)
    call fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
      istart + 1, iend + 1, dt, angularb, 1, no_e + 1, spectrum%energy_step, imspb)

    call batch_end(angularb)
    call batch_end(respb)
    call batch_end(imspb)
    
    sum1 = M_Z0
    sum2 = M_Z0
    do ie = 0, no_e
      energy = ie * spectrum%energy_step

      sp = cmplx(resp(ie), imsp(ie))

      sp = sp*M_ZI/(M_TWO*P_c*kick%delta_strength)

      sum1 = sum1 + spectrum%energy_step*sp
      sum2 = sum2 + spectrum%energy_step*sp*energy**2

      resp(ie) = real(sp)
      imsp(ie) = aimag(sp)
    end do

    SAFE_DEALLOCATE_A(angular)

    ! print some info
    write(message(1), '(a,i8)')    'Number of time steps = ', ntiter
    write(message(2), '(a,i4)')    'PropagationSpectrumDampMode   = ', spectrum%damp
    write(message(3), '(a,f10.4)') 'PropagationSpectrumDampFactor = ', units_from_atomic(units_out%time**(-1), spectrum%damp_factor)
    write(message(4), '(a,f10.4)') 'PropagationSpectrumStartTime  = ', units_from_atomic(units_out%time, spectrum%start_time)
    write(message(5), '(a,f10.4)') 'PropagationSpectrumEndTime    = ', units_from_atomic(units_out%time, spectrum%end_time)
    write(message(6), '(a,f10.4)') 'PropagationSpectrumMaxEnergy  = ', units_from_atomic(units_inp%energy, spectrum%max_energy) 
    write(message(7),'(a,f10.4)')  'PropagationSpectrumEnergyStep = ', units_from_atomic(units_inp%energy, spectrum%energy_step)
    message(8) = ""
    write(message(9), '(a,5e15.6,5e15.6)') 'R(0) sum rule = ', sum1
    write(message(10),'(a,5e15.6,5e15.6)') 'R(2) sum rule = ', sum2
    call messages_info(10)


    ! Output to file
    write(out_file, '(a15,i2)')      '# nspin        ', nspin
    call kick_write(kick, out_file)
    write(out_file, '(a1,a20,a20,a20)') '#', str_center("Energy", 20), str_center("R", 20), str_center("Re[beta]", 20)
    write(out_file, '(a1,a20,a20,a20)') '#', str_center('[' // trim(units_abbrev(units_out%energy)) // ']', 20), &
         str_center('[' // trim(units_abbrev(units_out%length**3)) // ']', 20), &
         str_center('[' // trim(units_abbrev(units_out%length**4)) // ']', 20)
    write(out_file, '(a,5e15.6,5e15.6)') '# R(0) sum rule = ', sum1
    write(out_file, '(a,5e15.6,5e15.6)') '# R(2) sum rule = ', sum2
    do ie = 0, no_e
      write(out_file,'(e20.8,e20.8,e20.8)') units_from_atomic(units_out%energy, ie*spectrum%energy_step), &
        units_from_atomic(units_out%length**3, imsp(ie)/M_PI), &
        units_from_atomic(units_out%length**4, resp(ie)*P_C/(M_THREE*max(ie, 1)*spectrum%energy_step))
    end do

    SAFE_DEALLOCATE_A(resp)
    SAFE_DEALLOCATE_A(imsp)

    POP_SUB(spectrum_rotatory_strength)
  end subroutine spectrum_rotatory_strength
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_init(dt, is, ie, niter, acc, from_vel)
    FLOAT,   intent(in)           :: dt
    integer, intent(in)           :: is, ie, niter
    CMPLX,   intent(in)           :: acc(:)
    logical, optional, intent(in) :: from_vel

    PUSH_SUB(spectrum_hsfunction_init)

    is_ = is
    ie_ = ie
    time_step_ = dt
    SAFE_ALLOCATE(func_(0:niter))
    func_ = acc
    
    if (present(from_vel)) then
      from_vel_ = from_vel
    else
      from_vel_ = .false.
    end if

    POP_SUB(spectrum_hsfunction_init)
  end subroutine spectrum_hsfunction_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_end
    PUSH_SUB(spectrum_hsfunction_end)
    SAFE_DEALLOCATE_A(func_)
    POP_SUB(spectrum_hsfunction_end)
  end subroutine spectrum_hsfunction_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_min(aa, bb, dw, omega, omega_min, func_min)
    FLOAT, intent(in)  :: aa, bb, omega, dw
    FLOAT, intent(out) :: omega_min, func_min

    integer :: ierr
    FLOAT :: xx, hsval, minhsval, ww

    PUSH_SUB(spectrum_hsfunction_min)

    ! xx should be an initial guess for the minimum. So we do a quick search
    ! that we refine later calling 1dminimize.
    xx = omega
    call hsfunction(xx, minhsval)
    ww = aa - dw
    do while(ww<bb)
      ww = ww + dw
      call hsfunction(ww, hsval)
      if(hsval < minhsval) then
        minhsval = hsval
        xx = ww
      end if
    end do 

    ! Around xx, we call some GSL sophisticated search algorithm to find the minimum.
#ifndef SINGLE_PRECISION
    call loct_1dminimize(max(xx-CNST(10.0)*dw, aa), min(xx+CNST(10.0)*dw,bb), xx, hsfunction, ierr)
#else
    stop "FIXME: cannot work in single-precision."
#endif

    if(ierr .ne. 0) then
      write(message(1),'(a,f14.6,a)') 'spectrum_hsfunction_min: The maximum at', xx,' was not properly converged.'
      write(message(2),'(a,i5)')      'Error code: ierr'
      call messages_warning(2)
    end if
    call hsfunction(xx, hsval)
    omega_min = xx
    func_min  = hsval

    POP_SUB(spectrum_hsfunction_min)
  end subroutine spectrum_hsfunction_min
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine hsfunction(omega, power)
    FLOAT, intent(in)             :: omega
    FLOAT, intent(out)            :: power

    CMPLX   :: cc, ez1, ez, zz
    integer :: jj,dir

    PUSH_SUB(hsfunction)

    cc = M_z0
    zz = M_zI * omega * time_step_
    ez1 = exp((is_ - 1) * zz)
    ez  = exp(zz)

    if (from_vel_) then
      do jj = is_, ie_
        ! This would be easier, but slower.
        !cc = cc + exp(M_zI * omega * jj * time_step_)*func_(jj)
        ez1 = ez1 * ez
        cc = cc + ez1 * func_(jj) * omega
      end do
    else
      do jj = is_, ie_
        ez1 = ez1 * ez
        cc = cc + ez1 * func_(jj)
      end do
    end if
    
    power = -abs(cc)**2 * time_step_**2

    if(allocated(func_)) then
      do jj = is_, ie_
        ! This would be easier, but slower.
        !cc = cc + exp(M_zI * omega * jj * time_step_)*func_(jj)
        ez1 = ez1 * ez
        cc = cc + ez1 * func_(jj)
      end do
      power = -abs(cc)**2 * time_step_**2
    end if
    
    if(allocated(func_ar_)) then 
      power=M_ZERO
      do dir=1, MAX_DIM
        ez1 = exp((is_ - 1) * zz)
        cc = M_z0
        do jj = is_, ie_
          ! This would be easier, but slower.
          !cc = cc + exp(M_zI * omega * jj * time_step_)*func_(jj)
          ez1 = ez1 * ez 
          cc = cc + ez1 * func_ar_(dir,jj) &
               *exp(-M_zI * omega * tret_(jj) ) !integrate over the retarded time
        end do
        power = power - abs(cc)**2 * time_step_**2
!        pp(dir) = cc * time_step_
      end do
!      power = - sum(abs(zcross_product(vv_, pp))**2)


    end if 


    POP_SUB(hsfunction)
  end subroutine hsfunction
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_ar_init(dt, is, ie, niter, acc, pos,tret)
    FLOAT,   intent(in) :: dt
    integer, intent(in) :: is, ie, niter
    CMPLX,   intent(in) :: acc(:,:),pos(:,:),tret(:)

    PUSH_SUB(spectrum_hsfunction_ar_init)

    is_ = is
    ie_ = ie
    time_step_ = dt
    SAFE_ALLOCATE(func_ar_(1:MAX_DIM,0:niter))
    SAFE_ALLOCATE(pos_(1:MAX_DIM,0:niter))
    SAFE_ALLOCATE(tret_(0:niter))
    func_ar_ = acc
    pos_= pos  
    tret_=tret   
 

    POP_SUB(spectrum_hsfunction_ar_init)
  end subroutine spectrum_hsfunction_ar_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_ar_end
    PUSH_SUB(spectrum_hsfunction_ar_end)
    SAFE_DEALLOCATE_A(func_ar_)
    SAFE_DEALLOCATE_A(pos_)
    SAFE_DEALLOCATE_A(tret_)
    POP_SUB(spectrum_hsfunction_ar_end)
  end subroutine spectrum_hsfunction_ar_end
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine spectrum_hs_ar_from_acc(out_file, spectrum, vec, w0)
    character(len=*), intent(in)    :: out_file
    type(spec_t),     intent(inout) :: spectrum
    FLOAT,            intent(in)    :: vec(:)
    FLOAT,  optional, intent(in)    :: w0

    integer :: istep, trash, iunit, nspin, time_steps, istart, iend, ntiter, lmax, ierr, jj
    FLOAT :: dt, dump,aa(MAX_DIM)
    type(kick_t) :: kick
    FLOAT, allocatable :: dd(:,:)
    CMPLX, allocatable :: acc(:,:),PP(:,:),pos(:,:),nn(:,:),tret(:)
    FLOAT :: vv(1:MAX_DIM)   
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_hs_ar_from_acc)

    call spectrum_acc_info(iunit, time_steps, dt)
    call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

    ! load dipole from file
    SAFE_ALLOCATE(acc(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(PP(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(pos(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(nn(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(tret(0:time_steps))

    acc = M_ZERO
    pos = M_ZERO
    PP = M_ZERO
    nn = M_ZERO
    tret = M_ZERO

    call io_skip_header(iunit)
    do istep = 0, time_steps-1
      aa = M_ZERO
      read(iunit, '(28x,e20.12)', advance = 'no', iostat = ierr) aa(1)
      ! What on earth is the point of this with jj??
      jj = 2
      do while( (ierr.eq.0) .and. (jj <= MAX_DIM) )
       read(iunit, '(e20.12)', advance = 'no', iostat = ierr) aa(jj)
       jj = jj+1
      end do

!      read(iunit, *) trash, dump, aa

      acc(:,istep) = units_to_atomic(units_out%acceleration, aa(:))
!      write (*,*) istep, Real(acc(:,istep))
    end do
    close(iunit)


    ! Try to get the trajectory from multipole file

    iunit = io_open('multipoles', action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/multipoles', action='read', status='old')
    end if
    if (.not.(iunit < 0)) then
      call spectrum_mult_info(iunit, nspin, kick, time_steps, dt, file_units, lmax=lmax)
      call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

      call io_skip_header(iunit)
!    write (*,*) 

      SAFE_ALLOCATE(dd(1:3, 1:nspin))
      do istep = 0, time_steps-1
        read(iunit, *) trash, dump, dump, dd
        pos(1:MAX_DIM,istep) = -sum(dd(1:MAX_DIM, :),2)
        pos(:,istep) = units_to_atomic(units_out%length, pos(:,istep))
!        write (*,*) istep, Real(pos(:,istep))
      end do
      SAFE_DEALLOCATE_A(dd)
      pos(:,0) = pos(:,1)
      call io_close(iunit)

    end if 

!    write (*,*)  

    ! normalize vector and set to global var
!    vv = vec / sqrt(sum(vec(:)**2))  
    vv = vec
    write (*,*) "vv",vv



    PP(:,0) = M_ZERO
    do istep = 0, time_steps - 1
       nn(:,istep) = vv(:)-pos(:,istep)
       nn(:,istep) = nn(:,istep)/sqrt(sum(nn(:,istep)**2 ))
       tret(istep) = ddot_product(vv(:),Real(pos(:,istep)) )/P_C 
       PP(:,istep) = zcross_product(nn, zcross_product(nn, acc(:,istep))) 
!       write (*,*) istep, Real(PP(:,istep)),"acc", Real (acc(:,istep))
    end do

    call spectrum_hsfunction_ar_init(dt, istart, iend, time_steps, PP, pos,tret)
    call spectrum_hs(out_file, spectrum, 'a', w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(acc)
    SAFE_DEALLOCATE_A(PP)
    SAFE_DEALLOCATE_A(pos)
    SAFE_DEALLOCATE_A(nn)
    SAFE_DEALLOCATE_A(tret)

    POP_SUB(spectrum_hs_ar_from_acc)
  end subroutine spectrum_hs_ar_from_acc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_ar_from_mult(out_file, spectrum, vec, w0)
    character(len=*), intent(in)    :: out_file
    type(spec_t),     intent(inout) :: spectrum
    FLOAT,            intent(in)    :: vec(:)
    FLOAT,  optional, intent(in)    :: w0

    integer :: istep, trash, iunit, nspin, time_steps, istart, iend, ntiter, lmax
    FLOAT :: dt, dump
    type(kick_t) :: kick
    FLOAT, allocatable :: dd(:,:)
    CMPLX, allocatable :: dipole(:,:), ddipole(:,:), PP(:,:), tret(:)
    CMPLX :: vv(MAX_DIM)   
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_hs_ar_from_mult)


    iunit = io_open('multipoles', action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/multipoles', action='read', status='old')
    end if
    call spectrum_mult_info(iunit, nspin, kick, time_steps, dt, file_units, lmax=lmax)
    call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

    call io_skip_header(iunit)

    ! load dipole from file
    SAFE_ALLOCATE(dipole(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(ddipole(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(PP(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(tret(0:time_steps))
    SAFE_ALLOCATE(dd(1:3, 1:nspin))
    
    tret= M_ZERO

    do istep = 1, time_steps
      read(iunit, *) trash, dump, dump, dd
      dipole(1:MAX_DIM,istep) = -sum(dd(1:MAX_DIM, :),2)
      dipole(:,istep) = units_to_atomic(units_out%length, dipole(:,istep))
    end do
    SAFE_DEALLOCATE_A(dd)
    dipole(:,0) = dipole(:,1)
    call io_close(iunit)

    ! we now calculate the acceleration.
    ddipole(:,0) = M_ZERO
    do istep = 1, time_steps - 1
      ddipole(:,istep) = (dipole(:,istep - 1) + dipole(:,istep + 1) - M_TWO * dipole(:,istep)) / dt**2
    end do
    call interpolate( dt*(/ -3, -2, -1 /),   &
                      ddipole(1,time_steps - 3:time_steps - 1), &
                      M_ZERO, &
                      ddipole(1,time_steps) )
    call interpolate( dt*(/ -3, -2, -1 /),   &
                      ddipole(2,time_steps - 3:time_steps - 1), &
                      M_ZERO, &
                      ddipole(2,time_steps) )
    call interpolate( dt*(/ -3, -2, -1 /),   &
                      ddipole(3,time_steps - 3:time_steps - 1), &
                      M_ZERO, &
                      ddipole(3,time_steps) )

   ! normalize vector and set to global var
    vv = vec / sqrt(sum(vec(:)**2))  

    PP(:,0) = M_ZERO
    do istep = 1, time_steps - 1
!      write (*,*) istep, istep*dt, Real(ddipole(1,istep)), Real(ddipole(2,istep))
       tret(istep) = zdot_product(vv(:),dipole(:,istep) )/P_C        
       PP(:,istep) = zcross_product(vv, zcross_product(vv, ddipole(:,istep - 1))) 
!      PP(istep) = sum(abs(dipole(:,istep))**2)
!      write(*,*) istep, PP(istep)
      
    end do


    call spectrum_hsfunction_ar_init(dt, istart, iend, time_steps, PP, dipole,tret)
    call spectrum_hs(out_file, spectrum, 'a', w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(dipole)
    SAFE_DEALLOCATE_A(ddipole)
    SAFE_DEALLOCATE_A(PP)
    SAFE_DEALLOCATE_A(tret)


    POP_SUB(spectrum_hs_ar_from_mult)
  end subroutine spectrum_hs_ar_from_mult
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_from_mult(out_file, spectrum, pol, vec, w0)
    character(len=*), intent(in)    :: out_file
    type(spec_t),     intent(inout) :: spectrum
    character,        intent(in)    :: pol
    FLOAT,            intent(in)    :: vec(:)
    FLOAT,  optional, intent(in)    :: w0

    integer :: istep, trash, iunit, nspin, time_steps, istart, iend, ntiter, lmax
    FLOAT :: dt, dump, vv(MAX_DIM)  
    type(kick_t) :: kick
    FLOAT, allocatable :: dd(:,:)
    CMPLX, allocatable :: dipole(:), ddipole(:)
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_hs_from_mult)

    iunit = io_open('multipoles', action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/multipoles', action='read', status='old')
    end if
    call spectrum_mult_info(iunit, nspin, kick, time_steps, dt, file_units, lmax=lmax)
    call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    call io_skip_header(iunit)

    ! load dipole from file
    SAFE_ALLOCATE(dipole(0:time_steps))
    SAFE_ALLOCATE(ddipole(0:time_steps))
    SAFE_ALLOCATE(dd(1:3, 1:nspin))

    vv = vec / sqrt(sum(vec(:)**2))  

    do istep = 1, time_steps
      read(iunit, *) trash, dump, dump, dd
      select case(pol)
      case('x')
        dipole(istep) = -sum(dd(1, :))
      case('y')
        dipole(istep) = -sum(dd(2, :))
      case('z')
        dipole(istep) =  sum(dd(3, :))
      case('+')
        dipole(istep) = -sum(dd(1, :) + M_zI * dd(2, :)) / sqrt(M_TWO)
      case('-')
        dipole(istep) = -sum(dd(1, :) - M_zI * dd(2, :)) / sqrt(M_TWO)
      case('v')
        dipole(istep) = -sum(vv(1)*dd(1, :) + vv(2)*dd(2, :) + vv(3)*dd(3, :))
      end select
      dipole(istep) = units_to_atomic(units_out%length, dipole(istep))
    end do
    SAFE_DEALLOCATE_A(dd)
    dipole(0) = dipole(1)
    call io_close(iunit)

    ! we now calculate the acceleration.
    ddipole(0) = M_ZERO
    do istep = 1, time_steps - 1
      ddipole(istep) = (dipole(istep - 1) + dipole(istep + 1) - M_TWO * dipole(istep)) / dt**2
    end do
    call interpolate( dt*(/ -3, -2, -1 /),   &
                      ddipole(time_steps - 3:time_steps - 1), &
                      M_ZERO, &
                      ddipole(time_steps) )

    call spectrum_hsfunction_init(dt, istart, iend, time_steps, ddipole)
    call spectrum_hs(out_file, spectrum, pol, w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(dipole)
    SAFE_DEALLOCATE_A(ddipole)

    POP_SUB(spectrum_hs_from_mult)
  end subroutine spectrum_hs_from_mult
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_from_acc(out_file, spectrum, pol, vec, w0)
    character(len=*), intent(in)    :: out_file
    type(spec_t),     intent(inout) :: spectrum
    character,        intent(in)    :: pol
    FLOAT,            intent(in)    :: vec(:)
    FLOAT,  optional, intent(in)    :: w0

    integer :: istep, jj, iunit, time_steps, istart, iend, ntiter, ierr
    FLOAT :: dt, aa(MAX_DIM),vv(MAX_DIM)
    CMPLX, allocatable :: acc(:)

    PUSH_SUB(spectrum_hs_from_acc)

    call spectrum_acc_info(iunit, time_steps, dt)
    call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! load dipole from file
    SAFE_ALLOCATE(acc(0:time_steps))
    acc = M_ZERO
    vv = vec / sqrt(sum(vec(:)**2))  
    call io_skip_header(iunit)
    do istep = 1, time_steps
      aa = M_ZERO
      read(iunit, '(28x,e20.12)', advance = 'no', iostat = ierr) aa(1)
      ! What on earth is the point of this with jj??
      jj = 2
      do while( (ierr.eq.0) .and. (jj <= MAX_DIM) )
        read(iunit, '(e20.12)', advance = 'no', iostat = ierr) aa(jj)
        jj = jj + 1 
      end do
      select case(pol)
      case('x')
        acc(istep) = aa(1)
      case('y')
        acc(istep) = aa(2)
      case('z')
        acc(istep) = aa(3)
      case('+')
        acc(istep) = (aa(1) + M_zI * aa(2)) / sqrt(M_TWO)
      case('-')
        acc(istep) = (aa(1) - M_zI * aa(2)) / sqrt(M_TWO)
      case('v')
        acc(istep) = vv(1)*aa(1) + vv(2)*aa(2) + vv(3)*aa(3)
      end select
      acc(istep) = units_to_atomic(units_out%acceleration, acc(istep))
    end do
    close(iunit)

    call spectrum_hsfunction_init(dt, istart, iend, time_steps, acc)
    call spectrum_hs(out_file, spectrum, pol, w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(acc)
    POP_SUB(spectrum_hs_from_acc)
  end subroutine spectrum_hs_from_acc

  ! ---------------------------------------------------------
  
  subroutine spectrum_hs_from_vel(out_file, spectrum, pol, w0)
    character(len=*), intent(in)    :: out_file
    type(spec_t),     intent(inout) :: spectrum
    character,        intent(in)    :: pol
    FLOAT,  optional, intent(in)    :: w0

    integer :: istep, jj, iunit, time_steps, istart, iend, ntiter, ierr
    FLOAT :: dt, aa(MAX_DIM)
    CMPLX, allocatable :: vel(:)
    logical :: from_vel

    PUSH_SUB(spectrum_hs_from_vel)

    call spectrum_vel_info(iunit, time_steps, dt)
    call spectrum_fix_time_limits(time_steps, dt, spectrum%start_time, spectrum%end_time, istart, iend, ntiter)

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! load dipole from file
    SAFE_ALLOCATE(vel(0:time_steps))
    vel = M_ZERO
    call io_skip_header(iunit)
    do istep = 1, time_steps
      aa = M_ZERO
      read(iunit, '(28x,e20.12)', advance = 'no', iostat = ierr) aa(1)
      ! What on earth is the point of this with jj??
      jj = 2
      do while( (ierr.eq.0) .and. (jj <= MAX_DIM) )
        read(iunit, '(e20.12)', advance = 'no', iostat = ierr) aa(jj)
      end do
      select case(pol)
      case('x')
        vel(istep) = aa(1)
      case('y')
        vel(istep) = aa(2)
      case('z')
        vel(istep) = aa(3)
      case('+')
        vel(istep) = (aa(1) + M_zI * aa(2)) / sqrt(M_TWO)
      case('-')
        vel(istep) = (aa(1) - M_zI * aa(2)) / sqrt(M_TWO)
      end select
      vel(istep) = units_to_atomic(units_out%velocity, vel(istep))
    end do
    close(iunit)

    from_vel = .true.
    call spectrum_hsfunction_init(dt, istart, iend, time_steps, vel, from_vel)
    call spectrum_hs(out_file, spectrum, pol, w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(vel)
    POP_SUB(spectrum_hs_from_vel)
  end subroutine spectrum_hs_from_vel
  ! ---------------------------------------------------------



  ! ---------------------------------------------------------
  subroutine spectrum_hs(out_file, spectrum, pol, w0)
    character(len=*), intent(in)    :: out_file
    type(spec_t),     intent(inout) :: spectrum
    character,        intent(in)    :: pol
    FLOAT,  optional, intent(in)    :: w0

    integer :: iunit, no_e, ie
    FLOAT   :: omega, hsval, xx
    FLOAT, allocatable :: sp(:)

    PUSH_SUB(spectrum_hs)

    if(present(w0)) then

      iunit = io_open(trim(out_file) // "." // trim(pol), action='write')
      write(iunit, '(a1,a20,a20)') '#', str_center("w", 20), str_center("H(w)", 20)
      write(iunit, '(a1,a20,a20)') '#', &
        str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
        str_center('[('//trim(units_abbrev(units_out%length))//'/' &
        //trim(units_abbrev(units_out%time**2)), 20)
      
      ! output
      omega = w0
      do while(omega <= spectrum%max_energy)
        call spectrum_hsfunction_min(omega - w0, omega + w0, spectrum%energy_step, omega, xx, hsval)

        write(iunit, '(1x,2e20.8)') units_from_atomic(units_out%energy, xx), &
          units_from_atomic((units_out%length / units_out%time)**2, -hsval)
          
        ! 2 * w0 because we assume that there are only odd peaks.
        omega = omega + 2 * w0
      end do
      call io_close(iunit)

    else
      no_e = spectrum%max_energy / spectrum%energy_step
      SAFE_ALLOCATE(sp(0:no_e))
      sp = M_ZERO

      do ie = 0, no_e
        call hsfunction(ie * spectrum%energy_step, sp(ie))
        sp(ie) = -sp(ie)
      end do

      ! output
      if(trim(out_file) .ne. '-') then
        iunit = io_open(trim(out_file) // "." // trim(pol), action='write')
        write(iunit, '(a1,a20,a20)') '#', str_center("w", 20), str_center("H(w)", 20)
        
        write(iunit, '(a1,a20,a20)') &
          '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
          str_center('[('//trim(units_abbrev(units_out%length))//'/' &
            //trim(units_abbrev(units_out%time**2)), 20)
        
        do ie = 0, no_e
          write(iunit, '(2e15.6)') units_from_atomic(units_out%energy, ie * spectrum%energy_step), &
            units_from_atomic((units_out%length / units_out%time)**2, sp(ie))
        end do
        
        call io_close(iunit)
      end if
      SAFE_DEALLOCATE_A(sp)

    end if

    POP_SUB(spectrum_hs)
  end subroutine spectrum_hs
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_mult_info(iunit, nspin, kick, time_steps, dt, file_units, lmax)
    integer,             intent(in)  :: iunit
    integer,             intent(out) :: nspin
    type(kick_t),        intent(out) :: kick
    integer,             intent(out) :: time_steps
    FLOAT,               intent(out) :: dt
    type(unit_system_t), intent(out) :: file_units
    integer,   optional, intent(out) :: lmax

    integer :: ii
    character(len=100) :: line

    PUSH_SUB(spectrum_mult_info)

    rewind(iunit)
    read(iunit,*)
    read(iunit,*)
    read(iunit, '(15x,i2)') nspin
    if(present(lmax)) then
      read(iunit, '(15x,i2)') lmax
    end if
    call kick_read(kick, iunit)
    read(iunit, '(a)') line
    read(iunit, '(a)') line
    call io_skip_header(iunit)

    ! Figure out the units of the file
    ii = index(line,'eV')
    if(ii.ne.0) then
      call unit_system_get(file_units, UNITS_EVA)
    else
      call unit_system_get(file_units, UNITS_ATOMIC)
    end if

    call count_time_steps(iunit, time_steps, dt)
    dt = units_to_atomic(file_units%time, dt) ! units_out is OK

    POP_SUB(spectrum_mult_info)
  end subroutine spectrum_mult_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------  
  subroutine count_time_steps(iunit, time_steps, dt)
    integer, intent(in)  :: iunit
    integer, intent(out) :: time_steps
    FLOAT,   intent(out) :: dt

    FLOAT :: t1, t2, dummy
    integer :: trash

    PUSH_SUB(count_time_steps)

    ! count number of time_steps
    time_steps = 0
    do
      read(iunit, *, end=100) trash, dummy
      time_steps = time_steps + 1
      if(time_steps == 1) t1 = dummy
      if(time_steps == 2) t2 = dummy
    end do
100 continue
    dt = (t2 - t1)
    time_steps = time_steps - 1
    
    if(time_steps < 3) then
      message(1) = "Empty file?"
      call messages_fatal(1)
    end if

    POP_SUB(count_time_steps)
  end subroutine count_time_steps
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section_info(iunit, nspin, kick, energy_steps, dw)
    integer, intent(in)           :: iunit
    integer, intent(out)          :: nspin
    type(kick_t), intent(out)     :: kick
    integer, intent(out)          :: energy_steps
    FLOAT,   intent(out)          :: dw            ! energy step

    FLOAT :: dummy, e1, e2

    PUSH_SUB(spectrum_cross_section_info)

    ! read in number of spin components
    read(iunit, '(15x,i2)') nspin
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
      call messages_fatal(1)
    end if

    POP_SUB(spectrum_cross_section_info)
  end subroutine spectrum_cross_section_info


  ! ---------------------------------------------------------
  subroutine spectrum_acc_info(iunit, time_steps, dt)
    integer, intent(out) :: iunit, time_steps
    FLOAT,   intent(out) :: dt

    integer :: trash
    FLOAT :: t1, t2, dummy

    PUSH_SUB(spectrum_acc_info)

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
      read(iunit, *, end=100) trash, dummy
      time_steps = time_steps + 1
      if(time_steps == 1) t1 = dummy
      if(time_steps == 2) t2 = dummy
    end do
100 continue
    dt = units_to_atomic(units_out%time, t2 - t1) ! units_out is OK
    time_steps = time_steps - 1

    if(time_steps < 3) then
      message(1) = "Empty multipole file?"
      call messages_fatal(1)
    end if

    rewind(iunit)
    POP_SUB(spectrum_acc_info)
  end subroutine spectrum_acc_info
  
    ! ---------------------------------------------------------
  subroutine spectrum_vel_info(iunit, time_steps, dt)
    integer, intent(out) :: iunit, time_steps
    FLOAT,   intent(out) :: dt

    integer :: trash
    FLOAT :: t1, t2, dummy

    PUSH_SUB(spectrum_vel_info)

    ! open files
    iunit = io_open('velocity', action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/velocity', action='read', status='old')
    end if

    ! read in dipole
    call io_skip_header(iunit)

    ! count number of time_steps
    time_steps = 0
    do
      read(iunit, *, end=100) trash, dummy
      time_steps = time_steps + 1
      if(time_steps == 1) t1 = dummy
      if(time_steps == 2) t2 = dummy
    end do
100 continue
    dt = units_to_atomic(units_out%time, t2 - t1) ! units_out is OK
    time_steps = time_steps - 1

    if(time_steps < 3) then
      message(1) = "Empty multipole file?"
      call messages_fatal(1)
    end if

    rewind(iunit)
    POP_SUB(spectrum_vel_info)
  end subroutine spectrum_vel_info



  ! ---------------------------------------------------------
  subroutine spectrum_fix_time_limits(time_steps, dt, start_time, end_time, istart, iend, ntiter)
    integer, intent(in)    :: time_steps
    FLOAT,   intent(in)    :: dt
    FLOAT,   intent(inout) :: start_time, end_time
    integer, intent(out)   :: istart, iend, ntiter

    FLOAT :: ts, te, dummy

    PUSH_SUB(spectrum_fix_time_limits)

    ts = M_ZERO
    te = time_steps * dt

    if(start_time < ts) start_time = ts
    if(start_time > te) start_time = te
    if(end_time   > te .or. end_time <= M_ZERO) end_time = te
    if(end_time   < ts) end_time   = ts

    if(end_time < start_time) then
      dummy = end_time ! swap
      end_time = start_time
      start_time = dummy
    end if
    istart = int(start_time / dt)
    iend = int(end_time / dt)
    ntiter = iend - istart + 1

    POP_SUB(spectrum_fix_time_limits)
  end subroutine spectrum_fix_time_limits

  ! -------------------------------------------------------

  subroutine signal_damp(damp_type, damp_factor, time_start, time_end, time_step, time_function)
    integer,         intent(in)    :: damp_type
    FLOAT,           intent(in)    :: damp_factor    
    integer,         intent(in)    :: time_start
    integer,         intent(in)    :: time_end
    FLOAT,           intent(in)    :: time_step
    type(batch_t),   intent(inout) :: time_function

    integer :: itime, ii
    FLOAT   :: total_time, time, weight

    PUSH_SUB(signal_damp)

    ASSERT(batch_is_ok(time_function))
    ASSERT(batch_status(time_function) == BATCH_NOT_PACKED)

    total_time = time_step*(time_end - time_start + 1)

    do itime = time_start, time_end
      time = time_step*(itime - time_start)

      ! Gets the damp function
      select case(damp_type)
      case(SPECTRUM_DAMP_NONE)
        weight = M_ONE
      case(SPECTRUM_DAMP_LORENTZIAN)
        weight = exp(-time*damp_factor)
      case(SPECTRUM_DAMP_POLYNOMIAL)
        weight = M_ONE - M_THREE*(time/total_time)**2 + M_TWO*(time/total_time)**3
      case(SPECTRUM_DAMP_GAUSSIAN)
        weight = exp(-time**2*damp_factor**2)
      end select
            
      do ii = 1, time_function%nst_linear
        time_function%states_linear(ii)%dpsi(itime) = weight*time_function%states_linear(ii)%dpsi(itime)
      end do
      
    end do

    POP_SUB(signal_damp)

  end subroutine signal_damp

  ! -------------------------------------------------------

  subroutine fourier_transform(method, transform, noise, time_start, time_end, time_step, time_function, &
    energy_start, energy_end, energy_step, energy_function)
    integer,         intent(in)    :: method
    integer,         intent(in)    :: transform
    FLOAT,           intent(in)    :: noise
    integer,         intent(in)    :: time_start
    integer,         intent(in)    :: time_end
    FLOAT,           intent(in)    :: time_step
    type(batch_t),   intent(in)    :: time_function
    integer,         intent(in)    :: energy_start
    integer,         intent(in)    :: energy_end
    FLOAT,           intent(in)    :: energy_step
    type(batch_t),   intent(inout) :: energy_function

    integer :: itime, ienergy, ii
    FLOAT   :: time, energy, kernel
    type(compressed_sensing_t) :: cs

    PUSH_SUB(fourier_transform)

    ASSERT(batch_is_ok(time_function))
    ASSERT(batch_is_ok(energy_function))
    ASSERT(time_function%nst_linear == energy_function%nst_linear)
    ASSERT(batch_status(time_function) == batch_status(energy_function))
    ASSERT(batch_status(time_function) == BATCH_NOT_PACKED)

    select case(method)
    case(SPECTRUM_FOURIER)

      do ienergy = energy_start, energy_end

        energy = energy_step*(ienergy - energy_start)

        forall(ii = 1:energy_function%nst_linear) energy_function%states_linear(ii)%dpsi(ienergy) = 0.0

        do itime = time_start, time_end

          time = time_step*(itime - time_start)

          select case(transform)
          case(SPECTRUM_TRANSFORM_SIN)
            kernel = sin(energy*time)
          case(SPECTRUM_TRANSFORM_COS)
            kernel = cos(energy*time)
          case(SPECTRUM_TRANSFORM_EXP)
            kernel = exp(-energy*time)
          end select

          kernel = kernel*time_step

          do ii = 1, time_function%nst_linear
            energy_function%states_linear(ii)%dpsi(ienergy) = &
              energy_function%states_linear(ii)%dpsi(ienergy) + kernel*time_function%states_linear(ii)%dpsi(itime) 
          end do

        end do
      end do

    case(SPECTRUM_COMPRESSED_SENSING)

      ASSERT(transform == SPECTRUM_TRANSFORM_SIN)

      call compressed_sensing_init(cs, time_end - time_start + 1, time_step, time_step*(time_start - 1), &
        energy_end - energy_start + 1, energy_step, energy_step*(energy_start - 1), noise)

      do ii = 1, time_function%nst_linear
        call compressed_sensing_spectral_analysis(cs, time_function%states_linear(ii)%dpsi, &
          energy_function%states_linear(ii)%dpsi)
      end do

      call compressed_sensing_end(cs)

    end select

    POP_SUB(fourier_transform)

  end subroutine fourier_transform

end module spectrum_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
