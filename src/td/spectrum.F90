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

module spectrum_oct_m
  use batch_oct_m
  use iso_c_binding
  use compressed_sensing_oct_m
  use fft_oct_m
  use global_oct_m
  use io_oct_m
  use kick_oct_m
  use lalg_adv_oct_m
  use math_oct_m
  use messages_oct_m
  use minimizer_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pcm_oct_m
  use profiling_oct_m
  use string_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                        &
    spectrum_t,                    &
    spectrum_init,                 &
    spectrum_cross_section,        &
    spectrum_cross_section_tensor, &
    spectrum_dipole_power,         &
    spectrum_dyn_structure_factor, &
    spectrum_rotatory_strength,    &
    spectrum_hs_from_mult,         &
    spectrum_hs_ar_from_mult,      &
    spectrum_hs_ar_from_acc,       &
    spectrum_hs_from_acc,          &
    spectrum_hsfunction_init,      &
    spectrum_hsfunction_end,       &
    spectrum_hsfunction_min,       &
    spectrum_mult_info,            &
    spectrum_fix_time_limits,      &
    spectrum_count_time_steps,     &
    spectrum_signal_damp,          &
    spectrum_fourier_transform,    &
    spectrum_hs_from_current,      &
    spectrum_nenergy_steps

  integer, public, parameter ::    &
    SPECTRUM_DAMP_NONE       = 0,  &
    SPECTRUM_DAMP_LORENTZIAN = 1,  &
    SPECTRUM_DAMP_POLYNOMIAL = 2,  &
    SPECTRUM_DAMP_GAUSSIAN   = 3

  integer, public, parameter ::    &
    SPECTRUM_TRANSFORM_LAPLACE = 1,  &
    SPECTRUM_TRANSFORM_SIN     = 2,  &
    SPECTRUM_TRANSFORM_COS     = 3

  integer, public, parameter :: &
    SPECTRUM_ABSORPTION    = 1, &
    SPECTRUM_ENERGYLOSS    = 2, &
    SPECTRUM_P_POWER       = 3, &
    SPECTRUM_ROTATORY      = 4

  integer, public, parameter ::       &
    SPECTRUM_FOURIER            = 1,  &
    SPECTRUM_COMPRESSED_SENSING = 2

  type spectrum_t
    FLOAT   :: start_time          !< start time for the transform
    FLOAT   :: end_time            !< when to stop the transform
    FLOAT   :: energy_step         !< step in energy mesh
    FLOAT   :: min_energy          !< minimum of energy mesh
    FLOAT   :: max_energy          !< maximum of energy mesh
    integer :: damp                !< damping type (none, exp or pol)
    integer :: transform           !< sine, cosine, or exponential transform
    FLOAT   :: damp_factor         !< factor used in damping
    integer :: spectype            !< spectrum type (absorption, energy loss, or dipole power)
    integer :: method              !< Fourier transform or compressed sensing 
    FLOAT   :: noise               !< the level of noise that is assumed in the time series for compressed sensing 
    logical, private :: sigma_diag          !< diagonalize sigma tensor
  end type spectrum_t

  !> Module variables, necessary to compute the function hsfunction, called by
  !! the C function loct_1dminimize
  integer :: niter_
  FLOAT :: time_step_, energy_step_
  CMPLX, allocatable :: func_(:),func_ar_(:,:),pos_(:,:),tret_(:), funcw_(:)
  type(fft_t), save :: fft_handler
  integer :: is_, ie_, default

contains

  ! ---------------------------------------------------------
  subroutine spectrum_init(spectrum, namespace, default_energy_step, default_max_energy)
    type(spectrum_t),           intent(inout) :: spectrum
    type(namespace_t),          intent(in)    :: namespace
    FLOAT,            optional, intent(in)    :: default_energy_step
    FLOAT,            optional, intent(in)    :: default_max_energy

    FLOAT :: fdefault

    PUSH_SUB(spectrum_init)
    
    call messages_print_stress(stdout, "Spectrum Options", namespace=namespace)

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
    !%Option DipolePower 3
    !% Power spectrum of the dipole moment.
    !%Option RotatoryStrength 4
    !% Rotatory strength spectrum.
    !%End

    call parse_variable(namespace, 'PropagationSpectrumType', SPECTRUM_ABSORPTION, spectrum%spectype)

    if(.not.varinfo_valid_option('PropagationSpectrumType', spectrum%spectype)) then
      call messages_input_error('PropagationSpectrumType')
    end if
    call messages_print_var_option(stdout, 'PropagationSpectrumType', spectrum%spectype)
      
    !%Variable SpectrumMethod
    !%Type integer
    !%Default fourier
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Decides which method is used to obtain the spectrum.
    !%Option fourier 1
    !% The standard Fourier transform. Further specified by <tt>PropagationSpectrumTransform</tt>.
    !%Option compressed_sensing 2
    !% (Experimental) Uses the compressed sensing technique.
    !%End
    call parse_variable(namespace, 'SpectrumMethod', SPECTRUM_FOURIER, spectrum%method)
    if(.not.varinfo_valid_option('SpectrumMethod', spectrum%method)) then
      call messages_input_error('SpectrumMethod')
    end if
    call messages_print_var_option(stdout, 'SpectrumMethod', spectrum%method)

    if(spectrum%method == SPECTRUM_COMPRESSED_SENSING) then
      call messages_experimental('compressed sensing')

      !%Variable SpectrumSignalNoise
      !%Type float
      !%Default 0.0
      !%Section Utilities::oct-propagation_spectrum
      !%Description
      !% For compressed sensing, the signal to process, the
      !% time-dependent dipole in this case, is assumed to have some
      !% noise that is given by this dimensionless quantity.
      !%End
      call parse_variable(namespace, 'SpectrumSignalNoise', CNST(0.0), spectrum%noise)
      call messages_print_var_value(stdout, 'SpectrumSignalNoise', spectrum%noise)
    end if
 

    !%Variable PropagationSpectrumDampMode
    !%Type integer
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Decides which damping/filtering is to be applied in order to
    !% calculate spectra by calculating a Fourier transform. The
    !% default is polynomial damping, except when <tt>SpectrumMethod = compressed_sensing</tt>.
    !% In that case the default is none.
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

    call parse_variable(namespace, 'PropagationSpectrumDampMode', default, spectrum%damp)

    if(.not.varinfo_valid_option('PropagationSpectrumDampMode', spectrum%damp)) then
      call messages_input_error('PropagationSpectrumDampMode')
    end if
    call messages_print_var_option(stdout, 'PropagationSpectrumDampMode', spectrum%damp)

    if(spectrum%method == SPECTRUM_COMPRESSED_SENSING .and. spectrum%damp /= SPECTRUM_DAMP_NONE) then
      message(1) = 'Using damping with compressed sensing, this is not required'
      message(2) = 'and can introduce noise in the spectra.'
      call messages_warning(2, namespace=namespace)
    end if

    !%Variable PropagationSpectrumTransform
    !%Type integer
    !%Default sine
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Decides which transform to perform, if <tt>SpectrumMethod = fourier</tt>.
    !%Option sine 2
    !% Sine transform: <math>\int dt \sin(wt) f(t)</math>. Produces the imaginary part of the polarizability.
    !%Option cosine 3
    !% Cosine transform: <math>\int dt \cos(wt) f(t)</math>. Produces the real part of the polarizability.
    !%Option laplace 1
    !% Real exponential transform: <math>\int dt e^{-wt} f(t)</math>. Produces the real part of the polarizability at imaginary
    !% frequencies, <i>e.g.</i> for Van der Waals <math>C_6</math> coefficients.
    !% This is the only allowed choice for complex scaling.
    !%End
    call parse_variable(namespace, 'PropagationSpectrumTransform', SPECTRUM_TRANSFORM_SIN, spectrum%transform)
    if(.not.varinfo_valid_option('PropagationSpectrumTransform', spectrum%transform)) then
      call messages_input_error('PropagationSpectrumTransform')
    end if
    call messages_print_var_option(stdout, 'PropagationSpectrumTransform', spectrum%transform)

    !%Variable PropagationSpectrumStartTime
    !%Type float
    !%Default 0.0
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Processing is done for the given function in a time-window that starts at the
    !% value of this variable.
    !%End
    call parse_variable(namespace, 'PropagationSpectrumStartTime',  M_ZERO, spectrum%start_time, units_inp%time)
    call messages_print_var_value(stdout, 'PropagationSpectrumStartTime', spectrum%start_time, unit = units_out%time)

    !%Variable PropagationSpectrumEndTime
    !%Type float
    !%Default -1.0 au
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Processing is done for the given function in a time-window that ends at the
    !% value of this variable. If set to a negative value, the maximum value from 
    !% the corresponding multipole file will used.
    !%End
    call parse_variable(namespace, 'PropagationSpectrumEndTime', -M_ONE, spectrum%end_time, units_inp%time)
    call messages_print_var_value(stdout, 'PropagationSpectrumEndTime', spectrum%end_time, unit = units_out%time)

    !%Variable PropagationSpectrumEnergyStep
    !%Type float
    !%Default 0.01 eV
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% Sampling rate for the spectrum. If you supply a number equal or smaller than zero, then
    !% the sampling rate will be <math>2 \pi / T</math>, where <math>T</math> is the total propagation time.
    !%End
    fdefault = CNST(0.01)/(M_TWO*P_Ry)
    if(present(default_energy_step)) fdefault = default_energy_step
    call parse_variable(namespace, 'PropagationSpectrumEnergyStep', fdefault, spectrum%energy_step, units_inp%energy)
    call messages_print_var_value(stdout, 'PropagationSpectrumEnergyStep', spectrum%energy_step, unit = units_out%energy)

    !%Variable PropagationSpectrumMinEnergy
    !%Type float
    !%Default 0 
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% The Fourier transform is calculated for energies larger than this value.
    !%End
    call parse_variable(namespace, 'PropagationSpectrumMinEnergy', M_ZERO, spectrum%min_energy, units_inp%energy)
    call messages_print_var_value(stdout, 'PropagationSpectrumMinEnergy', spectrum%min_energy, unit = units_out%energy)

    
    !%Variable PropagationSpectrumMaxEnergy
    !%Type float
    !%Default 20 eV
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% The Fourier transform is calculated for energies smaller than this value.
    !%End
    fdefault = CNST(20.0)/(M_TWO*P_Ry)
    if(present(default_max_energy)) fdefault = default_max_energy
    call parse_variable(namespace, 'PropagationSpectrumMaxEnergy', fdefault, spectrum%max_energy, units_inp%energy)
    call messages_print_var_value(stdout, 'PropagationSpectrumMaxEnergy', spectrum%max_energy, unit = units_out%energy)

    !%Variable PropagationSpectrumDampFactor
    !%Type float
    !%Default -1.0
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% If <tt>PropagationSpectrumDampMode = exponential, gaussian</tt>, the damping parameter of the exponential
    !% is fixed through this variable.
    !% Default value ensure that the damping function adquires a 0.0001 value at the end of the propagation time.
    !%End
    call parse_variable(namespace, 'PropagationSpectrumDampFactor', -M_ONE, spectrum%damp_factor, units_inp%time**(-1))

    call messages_print_var_value(stdout, 'PropagationSpectrumDampFactor', spectrum%damp_factor, unit = units_out%time**(-1))

    !%Variable PropagationSpectrumSigmaDiagonalization
    !%Type logical
    !%Default .false.
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% If <tt>PropagationSpectrumSigmaDiagonalization = yes</tt>, the polarizability tensor is diagonalizied.
    !% This variable is only used if the cross_section_tensor is computed. 
    !%End
    call parse_variable(namespace, 'PropagationSpectrumSigmaDiagonalization', .false., spectrum%sigma_diag)
    call messages_print_var_value(stdout, 'PropagationSpectrumSigmaDiagonalization', spectrum%sigma_diag)

    call messages_print_stress(stdout, namespace=namespace)

    POP_SUB(spectrum_init)
  end subroutine spectrum_init


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section_tensor(spectrum, namespace, out_file, in_file)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(in)    :: out_file
    integer,           intent(in)    :: in_file(:)

    integer :: nspin, energy_steps, ie, is, equiv_axes, n_files, trash
    FLOAT, allocatable :: sigma(:, :, :, :), sigmap(:, :, :, :), sigmau(:, :, :),  &
      sigmav(:, :, :), sigmaw(:, :, :), ip(:, :)
    FLOAT :: dw, dump
    type(kick_t) :: kick

    PUSH_SUB(spectrum_cross_section_tensor)

    n_files = size(in_file)
    equiv_axes = 3 - n_files + 1

    call spectrum_cross_section_info(namespace, in_file(1), nspin, kick, energy_steps, dw)
    ! on subsequent calls, do not overwrite energy_steps and dw
    call io_skip_header(in_file(1))

    SAFE_ALLOCATE(sigma (1:3, 1:3, 1:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmap(1:3, 1:3, 1:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmau(1:3,      1:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmav(1:3,      1:energy_steps, 1:nspin))
    SAFE_ALLOCATE(sigmaw(1:3,      1:energy_steps, 1:nspin))
    SAFE_ALLOCATE(    ip(1:3, 1:3))

    select case(equiv_axes)

    case(3)

      do ie = 1, energy_steps
        read(in_file(1), *) dump, (sigmau(1:3, ie, is), is = 1, nspin)
      end do

      ! The first row of sigma is the vector that we have just read, but properly projected...
      do is = 1, nspin
        do ie = 1, energy_steps
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
        do ie = 1, energy_steps
          sigmap(2, 3, ie, is) = sum(sigmau(1:3, ie, is) * kick%wprime(1:3))
          sigmap(3, 2, ie, is) = sigmap(2, 3, ie, is)
        end do
      end do

    case(2)

      call spectrum_cross_section_info(namespace, in_file(2), ie, kick, trash, dump)
      call io_skip_header(in_file(2))

      do ie = 1, energy_steps
        read(in_file(1), *) dump, (sigmau(1:3, ie, is), is = 1, nspin)
        read(in_file(2), *) dump, (sigmaw(1:3, ie, is), is = 1, nspin)
      end do

      ! The first row of sigma is the vector that we have just read, but properly projected...
      do is = 1, nspin
        do ie = 1, energy_steps
          sigmap(1, 1, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(1, 2, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(1, 3, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do

      ! The third row of sigma is also the vector that we have just read, but properly projected...
      do is = 1, nspin
        do ie = 1, energy_steps
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

      call spectrum_cross_section_info(namespace, in_file(2), ie, kick, trash, dump)
      call spectrum_cross_section_info(namespace, in_file(3), ie, kick, trash, dump)
      call io_skip_header(in_file(2))
      call io_skip_header(in_file(3))

      do ie = 1, energy_steps
        read(in_file(1), *) dump, (sigmau(1:3, ie, is), is = 1, nspin)
        read(in_file(2), *) dump, (sigmav(1:3, ie, is), is = 1, nspin)
        read(in_file(3), *) dump, (sigmaw(1:3, ie, is), is = 1, nspin)
      end do

      do is = 1, nspin
        do ie = 1, energy_steps
          sigmap(1, 1, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(1, 2, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(1, 3, ie, is) = sum( sigmau(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do
      do is = 1, nspin
        do ie = 1, energy_steps
          sigmap(2, 1, ie, is) = sum( sigmav(1:3, ie, is) * kick%pol(1:3, 1))
          sigmap(2, 2, ie, is) = sum( sigmav(1:3, ie, is) * kick%pol(1:3, 2))
          sigmap(2, 3, ie, is) = sum( sigmav(1:3, ie, is) * kick%pol(1:3, 3))
        end do
      end do
      do is = 1, nspin
        do ie = 1, energy_steps
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
      do ie = 1, energy_steps
        sigma(:, :, ie, is) = matmul( transpose(ip), matmul(sigmap(:, :, ie, is), ip) )
      end do
    end do

    ! Finally, write down the result
    call spectrum_cross_section_tensor_write(out_file, sigma, nspin, spectrum%energy_step, &
               spectrum%min_energy, energy_steps, kick)

    ! Diagonalize sigma tensor
    if (spectrum%sigma_diag) then
      call spectrum_sigma_diagonalize(namespace, sigma, nspin, spectrum%energy_step, spectrum%min_energy, energy_steps, kick)
    end if

    SAFE_DEALLOCATE_A(sigma)
    SAFE_DEALLOCATE_A(sigmap)
    SAFE_DEALLOCATE_A(sigmau)
    SAFE_DEALLOCATE_A(sigmav)
    SAFE_DEALLOCATE_A(sigmaw)
    SAFE_DEALLOCATE_A(ip)

    POP_SUB(spectrum_cross_section_tensor)
  end subroutine spectrum_cross_section_tensor


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section_tensor_write(out_file, sigma, nspin, energy_step, min_energy, energy_steps, kick)
    integer,                intent(in) :: out_file
    FLOAT,                  intent(in) :: sigma(:, :, :, :) !< (3, 3, energy_steps, nspin) already converted to units
    integer,                intent(in) :: nspin
    FLOAT,                  intent(in) :: energy_step, min_energy
    integer,                intent(in) :: energy_steps
    type(kick_t), optional, intent(in) :: kick !< if present, will write itself and nspin

    integer :: is, idir, jdir, ie, ii
    FLOAT :: average, anisotropy
    FLOAT, allocatable :: pp(:,:), pp2(:,:), ip(:,:)
    logical :: spins_singlet, spins_triplet
    character(len=20) :: header_string

    PUSH_SUB(spectrum_cross_section_tensor_write)

    spins_singlet = .true.
    spins_triplet = .false.
    if(present(kick)) then
      write(out_file, '(a15,i2)')      '# nspin        ', nspin
      call kick_write(kick, out_file)
      select case(kick%delta_strength_mode)
      case (KICK_SPIN_MODE)
        spins_triplet = .true.
        spins_singlet = .false.
      case (KICK_SPIN_DENSITY_MODE)
        spins_triplet = .true.
      end select
    end if

    write(out_file, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
    write(out_file, '(a20)', advance = 'no') str_center("(1/3)*Tr[sigma]", 20)
    write(out_file, '(a20)', advance = 'no') str_center("Anisotropy[sigma]", 20)
    if (spins_triplet .and. spins_singlet) then
      write(out_file, '(a20)', advance = 'no') str_center("(1/3)*Tr[sigma-]", 20)
    end if
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
    if (spins_triplet .and. spins_singlet) then
      write(out_file, '(a20)', advance = 'no')  str_center('[' // trim(units_abbrev(units_out%length**2)) // ']', 20)
    end if
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

    SAFE_ALLOCATE(pp(1:3, 1:3))
    if (spins_triplet .and. spins_singlet) SAFE_ALLOCATE(pp2(1:3, 1:3))
    SAFE_ALLOCATE(ip(1:3, 1:3))

    do ie = 1, energy_steps

      pp(:, :) = sigma(:, :, ie, 1)
      if (nspin >= 2) then
        if (spins_singlet .and. spins_triplet) then
          pp2(:, :) = pp(:, :) - sigma(:, :, ie, 2)
          pp(:, :)  = pp(:, :) + sigma(:, :, ie, 2)
        elseif (spins_triplet .and. .not.spins_singlet) then
          pp(:, :) = pp(:, :) - sigma(:, :, ie, 2)
        elseif (spins_singlet .and. .not.spins_triplet) then
          pp(:, :) = pp(:, :) + sigma(:, :, ie, 2)
        end if
      end if

      average = M_THIRD * ( pp(1, 1) + pp(2, 2) + pp(3, 3) )
      ip = matmul(pp, pp)
      anisotropy = M_THIRD * ( M_THREE * ( ip(1, 1) + ip(2, 2) + ip(3, 3) ) - (M_THREE * average)**2 )

      ! Note that the cross-section elements do not have to be transformed to the proper units, since
      ! they have been read from the "cross_section_vector.x", where they are already in the proper units.
      write(out_file,'(3e20.8)', advance = 'no') units_from_atomic(units_out%energy, ((ie-1) * energy_step + min_energy)), &
        average, sqrt(max(anisotropy, M_ZERO))

      if (spins_singlet .and. spins_triplet) then
        average =  M_THIRD * ( pp2(1, 1) + pp2(2, 2) + pp2(3, 3) )
        write(out_file,'(1e20.8)', advance = 'no') average
      end if

      do is = 1, nspin
        write(out_file,'(9e20.8)', advance = 'no') sigma(1:3, 1:3, ie, is)
      end do
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(pp)
    if (spins_triplet .and. spins_singlet) then 
      SAFE_DEALLOCATE_A(pp2)
    end if
    SAFE_DEALLOCATE_A(ip)
    POP_SUB(spectrum_cross_section_tensor_write)
  end subroutine spectrum_cross_section_tensor_write


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section(spectrum, namespace, in_file, out_file, ref_file)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(in)    :: in_file
    integer,           intent(in)    :: out_file
    integer, optional, intent(in)    :: ref_file

    character(len=20) :: header_string
    integer :: nspin, ref_nspin, lmax, ref_lmax, time_steps, &
               ref_time_steps, istart, iend, ntiter, it, ii, isp, no_e, ie, idir
    FLOAT   :: dt, ref_dt, energy, ewsum, polsum
    type(kick_t) :: kick, ref_kick
    FLOAT, allocatable :: dipole(:, :, :), ref_dipole(:, :, :), sigma(:, :, :), sf(:, :)
    type(unit_system_t) :: file_units, ref_file_units
    type(batch_t) :: dipoleb, sigmab

    type(pcm_min_t) :: pcm

    PUSH_SUB(spectrum_cross_section)

    ! This function gives us back the unit connected to the "multipoles" file, the header information,
    ! the number of time steps, and the time step.
    call spectrum_mult_info(namespace, in_file, nspin, kick, time_steps, dt, file_units, lmax=lmax)

    if(present(ref_file)) then
      call spectrum_mult_info(namespace, ref_file, ref_nspin, ref_kick, &
        ref_time_steps, ref_dt, ref_file_units, lmax = ref_lmax)
      if( (nspin /= ref_nspin)           .or. &
          (time_steps /= ref_time_steps) .or. &
          (.not.(dt .app. ref_dt))         .or. &
          (lmax /= ref_lmax) ) then
        write(message(1),'(a)') 'The multipoles and reference multipoles files do not match.'
        call messages_fatal(1, namespace=namespace)
      end if
    end if

    ! Now we cannot process files that do not contain the dipole, or that contain more than the dipole.
    if(lmax /= 1) then
      message(1) = 'Multipoles file should contain the dipole -- and only the dipole.'
      call messages_fatal(1, namespace=namespace)
    end if

    if(kick%function_mode /= KICK_FUNCTION_DIPOLE) then
      message(1) = "Kick function must have been dipole to run this utility."
      call messages_fatal(1, namespace=namespace)
    end if

    if(kick%pol_dir < 1) then
      message(1) = "Kick polarization direction is not set. Probably no kick was used."
      call messages_fatal(1, namespace=namespace)
    end if

    ! Find out the iteration numbers corresponding to the time limits.
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    SAFE_ALLOCATE(dipole(0:time_steps, 1:3, 1:nspin))
    call spectrum_read_dipole(namespace, in_file, dipole)

    if(present(ref_file)) then
      SAFE_ALLOCATE(ref_dipole(0:time_steps, 1:3, 1:nspin))
      call spectrum_read_dipole(namespace, ref_file, ref_dipole)
    end if

    ! parsing and re-printing to output useful PCM data
    call pcm_min_input_parsing_for_spectrum(pcm, namespace)

    ! adding the dipole generated by the PCM polarization charges due solute
    if(pcm%localf) &
      call spectrum_add_pcm_dipole(namespace, dipole, time_steps, dt, nspin)

    ! Now subtract the initial dipole.
    if(present(ref_file)) then
      dipole = dipole - ref_dipole
    else
      do it = 1, time_steps
        dipole(it, :, :) = dipole(it, :, :) - dipole(0, :, :)
      end do
      dipole(0, :, :) = M_ZERO
    end if

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! Get the number of energy steps.
    no_e = spectrum_nenergy_steps(spectrum)
    SAFE_ALLOCATE(sigma(1:no_e, 1:3, 1:nspin))

    if ( pcm%localf ) then

      ! for numerical reasons, we cannot add the difference (d(t)-d(t0)) of PCM dipoles here -- although it would look cleaner

      ! in the PCM local field case \sigma(\omega) \propto \Im{\alpha(\omega)\epsilon(\omega)} not just \propto \Im{\alpha(\omega)}
      ! since the dielectric function is complex as well, we need to compute both the real and imaginary part of the polarizability
      call spectrum_times_pcm_epsilon(spectrum, pcm, dipole, sigma, nspin, istart, iend, kick%time, dt, no_e)

      write(out_file,'(a57)') "Cross-section spectrum contains full local field effects."

    else

      call batch_init(dipoleb, 3, 1, nspin, dipole)
      call batch_init(sigmab, 3, 1, nspin, sigma)

      call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart + 1, iend + 1, kick%time, dt, dipoleb)
      call spectrum_fourier_transform(spectrum%method, spectrum%transform, spectrum%noise, &
       istart + 1, iend + 1, kick%time, dt, dipoleb, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, sigmab)

      call dipoleb%end()
      call sigmab%end()

    end if

    if( pcm%run_pcm ) &
      call spectrum_over_pcm_refraction_index(spectrum, pcm, sigma, nspin, no_e)



    SAFE_DEALLOCATE_A(dipole)
    if(present(ref_file)) then
      SAFE_DEALLOCATE_A(ref_dipole)
    end if

    SAFE_ALLOCATE(sf(1:no_e, nspin))

    if (abs(kick%delta_strength) < CNST(1e-12)) kick%delta_strength = M_ONE
    do ie = 1, no_e
      energy = (ie-1) * spectrum%energy_step + spectrum%min_energy
      forall(isp = 1:nspin) sf(ie, isp) = sum(sigma(ie, 1:3, isp)*kick%pol(1:3, kick%pol_dir))
      sf(ie, 1:nspin) = -sf(ie, 1:nspin) * (energy * M_TWO) / (M_PI * kick%delta_strength)
      sigma(ie, 1:3, 1:nspin) = -sigma(ie, 1:3, 1:nspin)*(M_FOUR*M_PI*energy/P_c)/kick%delta_strength
    end do
    
    ! The formulae below are only correct in this particular case.
    if(kick%delta_strength_mode == KICK_DENSITY_MODE .and. spectrum%transform == SPECTRUM_TRANSFORM_SIN) then
      ewsum = sum(sf(1, 1:nspin))
      polsum = M_ZERO

      do ie = 2, no_e
        energy = (ie-1) * spectrum%energy_step + spectrum%min_energy
        ewsum = ewsum + sum(sf(ie, 1:nspin))
        polsum = polsum + sum(sf(ie, 1:nspin)) / energy**2
      end do

      ewsum = ewsum * spectrum%energy_step
      polsum = polsum * spectrum%energy_step
    end if

    write(out_file, '(a15,i2)')      '# nspin        ', nspin
    call kick_write(kick, out_file)
    write(out_file, '(a)') '#%'
    write(out_file, '(a,i8)')    '# Number of time steps = ', time_steps
    call spectrum_write_info(spectrum, out_file)
    write(out_file, '(a)') '#%'
    if(kick%delta_strength_mode == KICK_DENSITY_MODE .and. spectrum%transform == SPECTRUM_TRANSFORM_SIN) then
      write(out_file, '(a,f16.6)') '# Electronic sum rule       = ', ewsum
      write(out_file, '(a,f16.6,1x,a)') '# Static polarizability (from sum rule) = ', &
        units_from_atomic(units_out%length**3, polsum), trim(units_abbrev(units_out%length))
      write(out_file, '(a)') '#%'
    end if
    
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
      write(out_file, '(a20)', advance = 'no') str_center('[' // trim(units_abbrev(unit_one/units_out%energy)) // ']', 20)
    end do
    write(out_file, '(1x)')

    do ie = 1, no_e
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, &
                                    (ie-1) * spectrum%energy_step + spectrum%min_energy)
      do isp = 1, nspin
        write(out_file,'(3e20.8)', advance = 'no') (units_from_atomic(units_out%length**2, sigma(ie, idir, isp)), &
                                                    idir = 1, 3)
      end do
      do isp = 1, nspin
        write(out_file,'(e20.8)', advance = 'no') units_from_atomic(unit_one/units_out%energy, sf(ie, isp))
      end do
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(sigma)
    POP_SUB(spectrum_cross_section)
  end subroutine spectrum_cross_section

  ! ---------------------------------------------------------

  subroutine spectrum_read_dipole(namespace, in_file, dipole)
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(in)    :: in_file
    FLOAT,             intent(out)   :: dipole(0:, :, :)

    integer :: nspin, lmax, time_steps, trash, it, idir, ispin
    FLOAT   :: dt,  dump
    type(kick_t) :: kick
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_read_dipole)

    ! This function gives us back the unit connected to the "multipoles" file, the header information,
    ! the number of time steps, and the time step.
    call spectrum_mult_info(namespace, in_file, nspin, kick, time_steps, dt, file_units, lmax = lmax)

    ! Read the dipole.
    call io_skip_header(in_file)

    do it = 0, time_steps
      read(in_file, *) trash, dump, (dump, (dipole(it, idir, ispin), idir = 1, 3), ispin = 1, nspin)
    end do
    dipole(:,:,:) = units_to_atomic(file_units%length, dipole(:,:,:))
    
    POP_SUB(spectrum_read_dipole)

  end subroutine spectrum_read_dipole

  ! ---------------------------------------------------------

  subroutine spectrum_add_pcm_dipole(namespace, dipole, time_steps, dt, nspin)
    type(namespace_t), intent(in)    :: namespace
    FLOAT,             intent(inout) :: dipole(0:, :, :)
    integer,           intent(in)    :: time_steps
    FLOAT,             intent(in)    :: dt
    integer,           intent(in)    :: nspin

    type(pcm_t) :: pcm
    FLOAT :: dipole_pcm(1:3)
    integer :: ia, it, ii

    ! unit io variables
    integer :: asc_unit_test
    integer :: cavity_unit
    integer :: asc_vs_t_unit, asc_vs_t_unit_check
    integer :: dipole_vs_t_unit_check, dipole_vs_t_unit_check1
    integer :: iocheck
    integer :: aux_int
    FLOAT :: aux_float, aux_float1, aux_vec(1:3)
    character(len=23) :: asc_vs_t_unit_format
    character(len=16) :: asc_vs_t_unit_format_tail

    PUSH_SUB(spectrum_add_pcm_dipole)

    ! reading PCM cavity from standard output file in two steps

    ! first step - counting tesserae
    asc_unit_test = io_open(PCM_DIR//'ASC_e.dat', namespace, action='read')
    pcm%n_tesserae = 0
    iocheck = 1
    do while( iocheck >= 0 )
      read(asc_unit_test,*,IOSTAT=iocheck) aux_vec(1:3), aux_float, aux_int
      if( iocheck >= 0 ) pcm%n_tesserae = pcm%n_tesserae + 1
    end do
    call io_close(asc_unit_test)

    ! intermezzo - allocating PCM tessellation and polarization charges arrays
    SAFE_ALLOCATE(pcm%tess(1:pcm%n_tesserae))
    SAFE_ALLOCATE(pcm%q_e(1:pcm%n_tesserae))
    SAFE_ALLOCATE(pcm%q_e_in(1:pcm%n_tesserae)) ! with auxiliary role

    ! second step - reading of PCM tessellation arrays from standard output file
    !               writing the cavity to debug-purpose file
    asc_unit_test = io_open(PCM_DIR//'ASC_e.dat', namespace, action='read')
    cavity_unit = io_open(PCM_DIR//'cavity_check.xyz', namespace, action='write')
    write(cavity_unit,'(I3)') pcm%n_tesserae
    write(cavity_unit,*)
    do ia = 1, pcm%n_tesserae
      read(asc_unit_test,*) pcm%tess(ia)%point(1:3), aux_float, aux_int
      write(cavity_unit,'(A1,3(1X,F14.8))') 'H', pcm%tess(ia)%point(1:3)
    end do
    call io_close(asc_unit_test)
    call io_close(cavity_unit)

    write (asc_vs_t_unit_format_tail,'(I5,A11)') pcm%n_tesserae,'(1X,F14.8))'
    write (asc_vs_t_unit_format,'(A)') '(F14.8,'//trim(adjustl(asc_vs_t_unit_format_tail))

    ! Now, summary: * read the time-dependent PCM polarization charges due to solute electrons, pcm%q_e
    !               * compute the real-time dipole generated by pcm%q_e, dipole_pcm
    !               * add it to the real-time molecular dipole
    !               * write the total dipole and its PCM contribution to debug-purpose files
    ! N.B. we assume nuclei fixed in time

    ! opening time-dependent PCM charges standard and debug-purpose file
    asc_vs_t_unit = io_open(PCM_DIR//'ASC_e_vs_t.dat', namespace, action='read', form='formatted')
    asc_vs_t_unit_check = io_open(PCM_DIR//'ASC_e_vs_t_check.dat', namespace, action='write', form='formatted')

    ! opening time-dependent PCM and total dipole debug-purpose files
    dipole_vs_t_unit_check = io_open(PCM_DIR//'dipole_e_vs_t_check.dat', namespace, action='write', form='formatted')
    dipole_vs_t_unit_check1 = io_open(PCM_DIR//'dipole_e_vs_t_check1.dat', namespace, action='write', form='formatted')

    ! reading PCM charges for the zeroth-iteration - not used - pcm%q_e_in is only auxiliary here
    read(asc_vs_t_unit,trim(adjustl(asc_vs_t_unit_format))) aux_float1, ( pcm%q_e_in(ia) , ia=1,pcm%n_tesserae )

    do it = 1, time_steps

      ! reading real-time PCM charges due to electrons per timestep
      read(asc_vs_t_unit,trim(adjustl(asc_vs_t_unit_format))) aux_float, ( pcm%q_e(ia) , ia=1,pcm%n_tesserae )

      ! computing real-time PCM dipole per timestep
      call pcm_dipole(dipole_pcm(1:3), -pcm%q_e(1:pcm%n_tesserae), pcm%tess, pcm%n_tesserae)

      ! adding PCM dipole to the molecular dipole per timestep
      dipole(it, 1, 1:nspin) = dipole(it, 1, 1:nspin) + dipole_pcm(1)
      dipole(it, 2, 1:nspin) = dipole(it, 2, 1:nspin) + dipole_pcm(2)
      dipole(it, 3, 1:nspin) = dipole(it, 3, 1:nspin) + dipole_pcm(3)

      ! since we always have a kick for the optical spectrum in Octopus
      ! the first-iteration dipole should be equal to the zeroth-iteration one
      ! in any case, made them equal by hand
      if ( it == 1 ) then
        dipole(0, 1, 1:nspin) = dipole(1, 1, 1:nspin)
        dipole(0, 2, 1:nspin) = dipole(1, 2, 1:nspin)
        dipole(0, 3, 1:nspin) = dipole(1, 3, 1:nspin)
      end if

      ! writing real-time PCM charges and dipole, and the total dipole for debug purposes
      write(asc_vs_t_unit_check,trim(adjustl(asc_vs_t_unit_format))) aux_float, ( pcm%q_e(ia) , ia=1,pcm%n_tesserae )
      write(dipole_vs_t_unit_check,'(F14.8,3(1X,F14.8))') aux_float, dipole_pcm
      write(dipole_vs_t_unit_check1,'(F14.8,3(1X,F14.8))') aux_float, dipole(it,:,1)

    end do

    ! closing PCM and debug files
    call io_close(asc_vs_t_unit)
    call io_close(asc_vs_t_unit_check)
    call io_close(dipole_vs_t_unit_check)
    call io_close(dipole_vs_t_unit_check1)

    ! deallocating PCM arrays
    SAFE_DEALLOCATE_A(pcm%tess)
    SAFE_DEALLOCATE_A(pcm%q_e)
    SAFE_DEALLOCATE_A(pcm%q_e_in)

    POP_SUB(spectrum_add_pcm_dipole)

  end subroutine spectrum_add_pcm_dipole

  ! ---------------------------------------------------------

  subroutine spectrum_times_pcm_epsilon(spectrum, pcm, dipole, sigma, nspin, istart, iend, kick_time, dt, no_e)
    type(spectrum_t),   intent(in)    :: spectrum
    type(pcm_min_t)   , intent(in)    :: pcm
    FLOAT, allocatable, intent(inout) :: sigma(:, :, :)
    FLOAT, allocatable, intent(in)    :: dipole(:, :, :)
    integer,            intent(in)    :: nspin
    FLOAT,              intent(in)    :: kick_time
    integer,            intent(in)    :: istart, iend
    FLOAT,              intent(in)    :: dt
    integer,            intent(in)    :: no_e

    FLOAT, allocatable :: sigmap(:, :, :)
    type(batch_t) :: dipoleb, sigmab

    integer :: ie

    CMPLX, allocatable :: eps(:)

    PUSH_SUB(spectrum_times_pcm_epsilon)

    ! imaginary part of the polarizability

    call batch_init(dipoleb, 3, 1, nspin, dipole)
    call batch_init(sigmab, 3, 1, nspin, sigma)

    call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart + 1, iend + 1, kick_time, dt, dipoleb)
    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
      istart + 1, iend + 1, kick_time, dt, dipoleb, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, sigmab)

    call dipoleb%end()
    call sigmab%end()

    ! real part of the polarizability

    SAFE_ALLOCATE(sigmap(1:no_e, 1:3, 1:nspin))

    call batch_init(dipoleb, 3, 1, nspin, dipole)
    call batch_init(sigmab, 3, 1, nspin, sigmap)

    call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart + 1, iend + 1, kick_time, dt, dipoleb)
    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
      istart + 1, iend + 1, kick_time, dt, dipoleb, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, sigmab)

    call dipoleb%end()
    call sigmab%end()

    SAFE_ALLOCATE(eps(1:no_e))

    ! multiplying by the dielectric function and taking the imaginary part of the product

    do ie = 1, no_e
      call pcm_eps(pcm, eps(ie), (ie-1)*spectrum%energy_step + spectrum%min_energy)
      sigma(ie, 1:3, 1:nspin) = sigma(ie, 1:3, 1:nspin) * REAL(eps(ie), REAL_PRECISION) + sigmap(ie, 1:3, 1:nspin) *AIMAG(eps(ie))
    end do

    SAFE_DEALLOCATE_A(sigmap)
    SAFE_DEALLOCATE_A(eps)

    POP_SUB(spectrum_times_pcm_epsilon)

  end subroutine spectrum_times_pcm_epsilon

  ! ---------------------------------------------------------

  subroutine spectrum_over_pcm_refraction_index(spectrum, pcm, sigma, nspin, no_e)
    type(spectrum_t),   intent(in)    :: spectrum
    type(pcm_min_t)   , intent(in)    :: pcm
    FLOAT, allocatable, intent(inout) :: sigma(:, :, :)
    integer,            intent(in)    :: nspin
    integer,            intent(in)    :: no_e

    integer :: ie

    CMPLX, allocatable :: eps(:)

    PUSH_SUB(spectrum_over_pcm_refraction_index)

    SAFE_ALLOCATE(eps(1:no_e))

    ! dividing by the refraction index - n(\omega)=\sqrt{\frac{|\epsilon(\omega)|+\Re[\epsilon(\omega)]}{2}}

    do ie = 1, no_e
      call pcm_eps(pcm, eps(ie), (ie-1)*spectrum%energy_step + spectrum%min_energy)
      sigma(ie, 1:3, 1:nspin) = sigma(ie, 1:3, 1:nspin) / sqrt( CNST(0.5) * ( ABS(eps(ie)) + REAL(eps(ie), REAL_PRECISION) ) )
    end do

    SAFE_DEALLOCATE_A(eps)

    POP_SUB(spectrum_over_pcm_refraction_index)

  end subroutine spectrum_over_pcm_refraction_index
  
  ! ---------------------------------------------------------
  subroutine spectrum_dipole_power(spectrum, namespace, in_file, out_file)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(in)    :: in_file
    integer,           intent(in)    :: out_file

    character(len=20) :: header_string
    integer :: nspin, lmax, time_steps, istart, iend, ntiter, it, ii, isp, no_e, ie, idir
    FLOAT   :: dt
    FLOAT, allocatable :: dipole(:, :, :), transform_cos(:, :, :), transform_sin(:, :, :), power(:, :, :)
    type(unit_system_t) :: file_units
    type(batch_t) :: dipoleb, transformb_cos, transformb_sin
    type(kick_t) :: kick

    PUSH_SUB(spectrum_dipole_power)

    ! This function gives us back the unit connected to the "multipoles" file, the header information,
    ! the number of time steps, and the time step.
    call spectrum_mult_info(namespace, in_file, nspin, kick, time_steps, dt, file_units, lmax=lmax)

    ! Now we cannot process files that do not contain the dipole, or that contain more than the dipole.
    if (lmax /= 1) then
      message(1) = 'Multipoles file should contain the dipole -- and only the dipole.'
      call messages_fatal(1, namespace=namespace)
    end if

    ! Find out the iteration numbers corresponding to the time limits.
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    SAFE_ALLOCATE(dipole(0:time_steps, 1:3, 1:nspin))
    call spectrum_read_dipole(namespace, in_file, dipole)

    ! Now subtract the initial dipole.
    do it = 1, time_steps
      dipole(it, :, :) = dipole(it, :, :) - dipole(0, :, :)
    end do
    dipole(0, :, :) = M_ZERO

    if (spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! Get the number of energy steps.
    no_e = spectrum_nenergy_steps(spectrum)
    SAFE_ALLOCATE(transform_cos(1:no_e, 1:3, 1:nspin))
    SAFE_ALLOCATE(transform_sin(1:no_e, 1:3, 1:nspin))
    SAFE_ALLOCATE(power(1:no_e, 1:3, 1:nspin))


    call batch_init(dipoleb, 3, 1, nspin, dipole)
    call batch_init(transformb_cos, 3, 1, nspin, transform_cos)
    call batch_init(transformb_sin, 3, 1, nspin, transform_sin)

    call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart + 1, iend + 1, spectrum%start_time, dt, dipoleb)

    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
         istart + 1, iend + 1, spectrum%start_time, dt, dipoleb, spectrum%min_energy,  &
         spectrum%max_energy, spectrum%energy_step, transformb_cos)
    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
         istart + 1, iend + 1, spectrum%start_time, dt, dipoleb, spectrum%min_energy, &
         spectrum%max_energy, spectrum%energy_step, transformb_sin)

    do ie = 1, no_e
      power(ie, :, :) = (transform_sin(ie, :, :)**2 + transform_cos(ie, :, :)**2)
    end do

    call dipoleb%end()
    call transformb_cos%end()
    call transformb_sin%end()

    SAFE_DEALLOCATE_A(dipole)
    SAFE_DEALLOCATE_A(transform_sin)
    SAFE_DEALLOCATE_A(transform_cos)

    write(out_file, '(a15,i2)')      '# nspin        ', nspin
    write(out_file, '(a)') '#%'
    write(out_file, '(a,i8)')    '# Number of time steps = ', time_steps
    call spectrum_write_info(spectrum, out_file)
    write(out_file, '(a)') '#%'
    
    write(out_file, '(a1,a20,1x)', advance = 'no') '#', str_center("Energy", 20)
    do isp = 1, nspin
      do idir = 1, 3
        write(header_string,'(a6,i1,a8,i1,a1)') 'power(', idir, ', nspin=', isp, ')'
        write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
      end do
    end do
    write(out_file, '(1x)')
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20)
    do ii = 1, nspin * 3
      write(out_file, '(a20)', advance = 'no') str_center('[' // trim(units_abbrev(units_out%length**2)) // ']', 20)
    end do
    write(out_file, '(1x)')

    do ie = 1, no_e
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, &
                                    (ie-1) * spectrum%energy_step + spectrum%min_energy)
      do isp = 1, nspin
        write(out_file,'(3e20.8)', advance = 'no') (units_from_atomic(units_out%length**2, power(ie, idir, isp)), &
                                                    idir = 1, 3)
      end do
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(power)

    POP_SUB(spectrum_dipole_power)
  end subroutine spectrum_dipole_power

  ! ---------------------------------------------------------
  subroutine spectrum_dyn_structure_factor(spectrum, namespace, in_file_sin, in_file_cos, out_file)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(in)    :: in_file_sin, in_file_cos
    integer,           intent(in)    :: out_file

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
    if(ii /= 0) then
      call unit_system_get(file_units, UNITS_EVA)
    else
      call unit_system_get(file_units, UNITS_ATOMIC)
    end if

    ! get time_steps and dt, and make sure that dt is the same in the two files
    call spectrum_count_time_steps(namespace, in_file_sin, time_steps_sin, dt_sin)
    call spectrum_count_time_steps(namespace, in_file_cos, time_steps_cos, dt_cos)

    if(dt_sin /= dt_cos) then
      message(1) = "dt is different in ftchds.cos and ftchds.sin!"
      call messages_fatal(1, namespace=namespace)
    end if

    time_steps = min(time_steps_sin, time_steps_cos)
    dt = units_to_atomic(file_units%time, dt_cos) ! units_out is OK

    ! Find out the iteration numbers corresponding to the time limits.
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    ! Read the f-transformed charge density.
    call io_skip_header(in_file_sin)
    call io_skip_header(in_file_cos)

    SAFE_ALLOCATE(ftchd(0:time_steps))
    do it = 0, time_steps
      read(in_file_sin, *) trash, dump, dummy1, dummy2
      read(in_file_cos, *) trash, dump, dummy3, dummy4
      ftchd(it) = TOCMPLX(dummy3-dummy2, dummy4+dummy1)
    end do

    ! Now subtract the initial value.
    do it = 1, time_steps
      ftchd(it) = ftchd(it) - ftchd(0)
    end do
    ftchd(0) = M_ZERO

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! Get the number of energy steps.
    no_e = spectrum_nenergy_steps(spectrum)

    SAFE_ALLOCATE(chi(1:no_e))
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
    if (abs(kick%delta_strength) < 1.d-12) kick%delta_strength = M_ONE
    do ie = 1, no_e
      energy = (ie-1) * spectrum%energy_step + spectrum%min_energy
      do it = istart, iend
        jj = it - istart

        xx = exp(M_zI * energy * jj * dt)
        chi(ie) = chi(ie) + xx * damp(it) * ftchd(it)

      end do
      chi(ie) = chi(ie) * dt / kick%delta_strength / M_PI
    end do

    ! Test f-sum rule
    fsum = M_ZERO
    do ie = 1, no_e
      energy = (ie-1) * spectrum%energy_step + spectrum%min_energy
      fsum = fsum + energy * aimag(chi(ie))
    end do
    fsum = spectrum%energy_step * fsum * 2/sum(kick%qvector(:,1)**2)

    write(out_file, '(a)') '#%'
    write(out_file, '(a,i8)')    '# Number of time steps = ', time_steps
    call spectrum_write_info(spectrum, out_file)
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

    do ie = 1, no_e
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, &
                                  (ie-1) * spectrum%energy_step + spectrum%min_energy)
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy**(-1), aimag(chi(ie)))
      write(out_file, '(1x)')
    end do

    SAFE_DEALLOCATE_A(ftchd)
    SAFE_DEALLOCATE_A(chi)
    POP_SUB(spectrum_dyn_structure_factor)

  end subroutine spectrum_dyn_structure_factor


  ! ---------------------------------------------------------
  subroutine spectrum_rotatory_strength(spectrum, namespace, in_file, out_file)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(in)    :: in_file
    integer,           intent(in)    :: out_file

    integer :: istart, iend, ntiter, ie, idir, time_steps, no_e, nspin, trash, it
    FLOAT :: dump, dt, energy
    type(kick_t) :: kick
    CMPLX :: sum1, sum2, sp
    FLOAT, allocatable :: angular(:, :), resp(:), imsp(:)
    type(batch_t) :: angularb, respb, imspb
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_rotatory_strength)

    call spectrum_mult_info(namespace, in_file, nspin, kick, time_steps, dt, file_units)
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    ! load angular momentum from file
    SAFE_ALLOCATE(angular(0:time_steps, 1:3))
    call io_skip_header(in_file)
    do ie = 0, time_steps
      read(in_file, *) trash, dump, (angular(ie, idir), idir = 1, 3)
    end do

    ! subtract static dipole
    do idir = 1, 3
      angular(:, idir) = angular(:, idir) - angular(0, idir)
    end do

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    no_e = spectrum_nenergy_steps(spectrum)

    do it = istart, iend
      angular(it, 1) = sum(angular(it, 1:3)*kick%pol(1:3, kick%pol_dir))
    end do

    SAFE_ALLOCATE(resp(1:no_e))
    SAFE_ALLOCATE(imsp(1:no_e))

    call batch_init(angularb, 1, 1)
    call batch_init(respb, 1, 1)
    call batch_init(imspb, 1, 1)

    call angularb%add_state(angular(:, 1))
    call respb%add_state(resp)
    call imspb%add_state(imsp)

    call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart + 1, iend + 1, kick%time, dt, angularb)

    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
      istart + 1, iend + 1, kick%time, dt, angularb, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, respb)
    call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
      istart + 1, iend + 1, kick%time, dt, angularb, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, imspb)

    call angularb%end()
    call respb%end()
    call imspb%end()
    
    sum1 = M_Z0
    sum2 = M_Z0
    if (abs(kick%delta_strength) < 1.d-12) kick%delta_strength = M_ONE
    do ie = 1, no_e
      energy = (ie-1) * spectrum%energy_step + spectrum%min_energy

      sp = TOCMPLX(resp(ie), imsp(ie))

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
    do ie = 1, no_e
      write(out_file,'(e20.8,e20.8,e20.8)') units_from_atomic(units_out%energy, (ie-1)*spectrum%energy_step+spectrum%min_energy), &
        units_from_atomic(units_out%length**3, imsp(ie)/M_PI), &
        units_from_atomic(units_out%length**4, resp(ie)*P_C/(M_THREE*max((ie-1),1)*spectrum%energy_step))
    end do

    SAFE_DEALLOCATE_A(resp)
    SAFE_DEALLOCATE_A(imsp)

    POP_SUB(spectrum_rotatory_strength)
  end subroutine spectrum_rotatory_strength
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_init(dt, is, ie, niter, acc)
    FLOAT,             intent(in)    :: dt
    integer,           intent(in)    :: is, ie, niter
    CMPLX,             intent(in)    :: acc(:)

    integer :: nn(3), j, optimize_parity(3)
    logical :: optimize(3)

    PUSH_SUB(spectrum_hsfunction_init)

    is_ = is
    ie_ = ie
    time_step_ = dt
    niter_ = niter
    energy_step_ = (M_TWO * M_PI) / (niter * time_step_)
    SAFE_ALLOCATE(func_(0:niter))
    SAFE_ALLOCATE(funcw_(0:niter))
    func_ = M_z0
    func_ = acc
    nn(1:3) = (/ niter, 1, 1 /)
    optimize(1:3) = .false.
    optimize_parity(1:3) = -1

    call fft_init(fft_handler, nn(1:3), 1, FFT_COMPLEX, FFTLIB_FFTW, optimize, optimize_parity)
    call zfft_forward1(fft_handler, func_(0:niter-1), funcw_(0:niter-1))
    do j = 0, niter - 1
      funcw_(j) = -abs(funcw_(j))**2 * dt**2
    end do

    POP_SUB(spectrum_hsfunction_init)
  end subroutine spectrum_hsfunction_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_end

    PUSH_SUB(spectrum_hsfunction_end)

    call fft_end(fft_handler)

    SAFE_DEALLOCATE_A(func_)
    SAFE_DEALLOCATE_A(funcw_)
    POP_SUB(spectrum_hsfunction_end)
  end subroutine spectrum_hsfunction_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hsfunction_min(namespace, aa, bb, omega_min, func_min)
    type(namespace_t), intent(in)  :: namespace
    FLOAT,             intent(in)  :: aa, bb
    FLOAT,             intent(out) :: omega_min, func_min

    integer :: ierr, ie
    FLOAT :: xx, hsval, minhsval, ww, xa, xb, hxa, hxb

    PUSH_SUB(spectrum_hsfunction_min)

    ! xx should be an initial guess for the minimum. So we do a quick search
    ! that we refine later calling 1dminimize.
    !xx = omega
    !call hsfunction(xx, minhsval)

    ierr = 0

    ie = int(aa/energy_step_)
    ww = ie * energy_step_
    if(ww < aa) then
      ie = ie + 1
      ww = ie * energy_step_
    end if
    xx = ie * energy_step_
    minhsval = real(funcw_(ie))
    do while(ww <= bb)
      hsval = real(funcw_(ie))
      if(hsval < minhsval) then
        minhsval = hsval
        xx = ww
      end if
      ie = ie + 1
      ww = ie * energy_step_
    end do

    ! Around xx, we call some GSL sophisticated search algorithm to find the minimum.
    ! First, we get the value of the function at the extremes of the interval
    xa = max(xx-energy_step_, aa)
    xb = min(xx+energy_step_, bb)
    call hsfunction(xa, hxa)
    call hsfunction(xb, hxb)

    if(hxa <= minhsval) then
      xx = xa
      minhsval = hxa
    elseif(hxb <= minhsval) then
      xx = xb
      minhsval = hxb
    else
      call loct_1dminimize(xa, xb, xx, hsfunction, ierr)
    end if

    if(ierr /= 0) then
      write(message(1),'(a,f14.6,a)') 'spectrum_hsfunction_min: The maximum at', xx,' was not properly converged.'
      write(message(2),'(a,i12)')      'Error code: ', ierr
      call messages_warning(2, namespace=namespace)
    end if
    call hsfunction(xx, hsval)
    omega_min = xx
    func_min  = hsval

    POP_SUB(spectrum_hsfunction_min)
  end subroutine spectrum_hsfunction_min
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine hsfunction(omega, power)
    FLOAT, intent(in)  :: omega
    FLOAT, intent(out) :: power

    CMPLX   :: cc, ez1, ez, zz
    integer :: jj,dir

    PUSH_SUB(hsfunction)

    cc = M_z0
    zz = M_zI * omega * time_step_
    ez1 = exp((is_ - 1) * zz)
    ez  = exp(zz)
    do jj = is_, ie_
      ez1 = ez1 * ez
      cc = cc + ez1 * func_(jj)
    end do
    power = -abs(cc)**2 * time_step_**2

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
  ! FIXME: why is this never called?
  subroutine spectrum_hsfunction_ar_end

    PUSH_SUB(spectrum_hsfunction_ar_end)

    SAFE_DEALLOCATE_A(func_ar_)
    SAFE_DEALLOCATE_A(pos_)
    SAFE_DEALLOCATE_A(tret_)

    POP_SUB(spectrum_hsfunction_ar_end)
  end subroutine spectrum_hsfunction_ar_end
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine spectrum_hs_ar_from_acc(spectrum, namespace, out_file, vec, w0)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: out_file
    FLOAT,             intent(in)    :: vec(:)
    FLOAT,  optional,  intent(in)    :: w0

    integer :: istep, trash, iunit, nspin, time_steps, istart, iend, ntiter, lmax, ierr, jj
    FLOAT :: dt, dump, aa(MAX_DIM)
    CMPLX :: nn(MAX_DIM)
    type(kick_t) :: kick
    FLOAT, allocatable :: dd(:,:)
    CMPLX, allocatable :: acc(:,:),PP(:,:),pos(:,:),tret(:)
    FLOAT :: vv(1:MAX_DIM)   
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_hs_ar_from_acc)

    call spectrum_tdfile_info(namespace, 'acceleration', iunit, time_steps, dt)
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    ! load dipole from file
    SAFE_ALLOCATE(acc(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(PP(1:MAX_DIM,0:time_steps))
    SAFE_ALLOCATE(pos(1:MAX_DIM,0:time_steps))
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
      ! FIXME: parsing of file depends on how code was compiled (MAX_DIM)!!!
      jj = 2
      do while( (ierr == 0) .and. (jj <= MAX_DIM) )
       read(iunit, '(e20.12)', advance = 'no', iostat = ierr) aa(jj)
       jj = jj+1
      end do

!      read(iunit, *) trash, dump, aa

      acc(:,istep) = units_to_atomic(units_out%acceleration, aa(:))
!      write (*,*) istep, Real(acc(:,istep))
    end do
    close(iunit)


    ! Try to get the trajectory from multipole file

    iunit = io_open('multipoles', namespace, action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/multipoles', namespace, action='read', status='old')
    end if
    if (.not.(iunit < 0)) then
      call spectrum_mult_info(namespace, iunit, nspin, kick, time_steps, dt, file_units, lmax=lmax)
      call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

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
    vv(1:3) = vec(1:3)

    PP(:,0) = M_ZERO
    do istep = 0, time_steps - 1
       nn(:) = vv(:)-pos(:,istep)
       nn(:) = nn(:)/sqrt(sum(nn(:)**2 ))
       tret(istep) = ddot_product(vv(:),real(pos(:,istep),REAL_PRECISION))/P_C 
       PP(:,istep) = zcross_product(nn, zcross_product(nn, acc(:,istep))) 
!       write (*,*) istep, Real(PP(:,istep)),"acc", Real (acc(:,istep))
    end do

    call spectrum_hsfunction_ar_init(dt, istart, iend, time_steps, PP, pos, tret)
    call spectrum_hs(spectrum, namespace, out_file, 'a', w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(acc)
    SAFE_DEALLOCATE_A(PP)
    SAFE_DEALLOCATE_A(pos)
    SAFE_DEALLOCATE_A(tret)

    POP_SUB(spectrum_hs_ar_from_acc)
  end subroutine spectrum_hs_ar_from_acc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_ar_from_mult(spectrum, namespace, out_file, vec, w0)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: out_file
    FLOAT,             intent(in)    :: vec(:)
    FLOAT,  optional,  intent(in)    :: w0

    integer :: istep, trash, iunit, nspin, time_steps, istart, iend, ntiter, lmax
    FLOAT :: dt, dump
    type(kick_t) :: kick
    FLOAT, allocatable :: dd(:,:)
    CMPLX, allocatable :: dipole(:,:), ddipole(:,:), PP(:,:), tret(:)
    CMPLX :: vv(MAX_DIM)   
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_hs_ar_from_mult)


    iunit = io_open('multipoles', namespace, action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/multipoles', namespace, action='read', status='old')
    end if
    call spectrum_mult_info(namespace, iunit, nspin, kick, time_steps, dt, file_units, lmax=lmax)
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

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
    vv(1:3) = vec(1:3) / sqrt(sum(vec(1:3)**2))  

    PP(:,0) = M_ZERO
    do istep = 1, time_steps - 1
!      write (*,*) istep, istep*dt, Real(ddipole(1,istep)), Real(ddipole(2,istep))
       tret(istep) = zdot_product(vv(:),dipole(:,istep) )/P_C        
       PP(:,istep) = zcross_product(vv, zcross_product(vv, ddipole(:,istep - 1))) 
!      PP(istep) = sum(abs(dipole(:,istep))**2)
!      write(*,*) istep, PP(istep)
      
    end do


    call spectrum_hsfunction_ar_init(dt, istart, iend, time_steps, PP, dipole,tret)
    call spectrum_hs(spectrum, namespace, out_file, 'a', w0)
    call spectrum_hsfunction_end()

    SAFE_DEALLOCATE_A(dipole)
    SAFE_DEALLOCATE_A(ddipole)
    SAFE_DEALLOCATE_A(PP)
    SAFE_DEALLOCATE_A(tret)


    POP_SUB(spectrum_hs_ar_from_mult)
  end subroutine spectrum_hs_ar_from_mult
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_from_mult(spectrum, namespace, out_file, pol, vec, w0)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: out_file
    character,         intent(in)    :: pol
    FLOAT,             intent(in)    :: vec(:)
    FLOAT,   optional, intent(in)    :: w0

    integer :: istep, trash, iunit, nspin, time_steps, istart, iend, ntiter, lmax, no_e, ie
    FLOAT :: dt, dump, vv(MAX_DIM)  
    type(kick_t) :: kick
    FLOAT, allocatable :: dd(:,:)
    FLOAT, allocatable :: sps(:), spc(:), racc(:)
    CMPLX, allocatable :: dipole(:), ddipole(:)
    type(batch_t) :: acc_batch, sps_batch, spc_batch
    type(unit_system_t) :: file_units

    PUSH_SUB(spectrum_hs_from_mult)

    iunit = io_open('multipoles', namespace, action='read', status='old', die=.false.)
    if(iunit < 0) then
      iunit = io_open('td.general/multipoles', namespace, action='read', status='old')
    end if
    call spectrum_mult_info(namespace, iunit, nspin, kick, time_steps, dt, file_units, lmax=lmax)
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    call io_skip_header(iunit)

    ! load dipole from file
    SAFE_ALLOCATE(dipole(0:time_steps))
    SAFE_ALLOCATE(ddipole(0:time_steps))
    SAFE_ALLOCATE(dd(1:3, 1:nspin))

    vv(1:3) = vec(1:3) / sqrt(sum(vec(1:3)**2))  

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

    if(present(w0)) then

      call spectrum_hsfunction_init(dt, istart, iend, time_steps, ddipole)
      call spectrum_hs(spectrum, namespace, out_file, pol, w0)
      call spectrum_hsfunction_end()

    else

      SAFE_ALLOCATE(racc(0:time_steps))
      racc = real(ddipole, REAL_PRECISION)

      no_e = spectrum_nenergy_steps(spectrum)
      SAFE_ALLOCATE(sps(1:no_e))
      SAFE_ALLOCATE(spc(1:no_e))
      sps = M_ZERO
      spc = M_ZERO

      call batch_init(acc_batch, 1, 1)
      call batch_init(sps_batch, 1, 1)
      call batch_init(spc_batch, 1, 1)

      call acc_batch%add_state(racc)
      call sps_batch%add_state(sps)
      call spc_batch%add_state(spc)

      call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, acc_batch, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, spc_batch)
      call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, acc_batch, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, sps_batch)

      do ie = 1, no_e
        sps(ie) = (sps(ie)**2 + spc(ie)**2)
      end do

      call spectrum_hs_output(spectrum, namespace, out_file, pol, no_e, sps)   

      call acc_batch%end()
      call sps_batch%end()
      call spc_batch%end()

      SAFE_DEALLOCATE_A(racc)

    end if

    SAFE_DEALLOCATE_A(dipole)
    SAFE_DEALLOCATE_A(ddipole)

    POP_SUB(spectrum_hs_from_mult)
  end subroutine spectrum_hs_from_mult
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_from_acc(spectrum, namespace, out_file, pol, vec, w0)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: out_file
    character,         intent(in)    :: pol
    FLOAT,             intent(in)    :: vec(:)
    FLOAT,   optional, intent(in)    :: w0

    integer :: istep, jj, iunit, time_steps, istart, iend, ntiter, ierr, no_e, ie
    FLOAT :: dt, aa(MAX_DIM),vv(MAX_DIM)
    CMPLX, allocatable :: acc(:)
    FLOAT, allocatable :: racc(:), sps(:), spc(:)
    type(batch_t) :: acc_batch, sps_batch, spc_batch

    PUSH_SUB(spectrum_hs_from_acc)

    call spectrum_tdfile_info(namespace, 'acceleration', iunit, time_steps, dt)
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! load dipole from file
    SAFE_ALLOCATE(acc(0:time_steps))
    acc = M_ZERO
    vv = vec / sqrt(sum(vec(:)**2))  
    call io_skip_header(iunit)

    do istep = 1, time_steps
      aa = M_ZERO
      read(iunit, '(28x,e20.12)', advance = 'no', iostat = ierr) aa(1)
      ! FIXME: parsing of file depends on how code was compiled (MAX_DIM)!!!
      jj = 2
      do while( (ierr == 0) .and. (jj <= MAX_DIM) )
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

    if(present(w0)) then

      call spectrum_hsfunction_init(dt, istart, iend, time_steps, acc)
      call spectrum_hs(spectrum, namespace, out_file, pol, w0)
      call spectrum_hsfunction_end()

    else

      SAFE_ALLOCATE(racc(0:time_steps))
      racc = real(acc, REAL_PRECISION)

      no_e = spectrum_nenergy_steps(spectrum)
      SAFE_ALLOCATE(sps(1:no_e))
      SAFE_ALLOCATE(spc(1:no_e))
      sps = M_ZERO
      spc = M_ZERO

      call batch_init(acc_batch, 1, 1)
      call batch_init(sps_batch, 1, 1)
      call batch_init(spc_batch, 1, 1)

      call acc_batch%add_state(racc)
      call sps_batch%add_state(sps)
      call spc_batch%add_state(spc)

      call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, acc_batch, spectrum%min_energy, &
        spectrum%max_energy, spectrum%energy_step, spc_batch)
      call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, acc_batch, spectrum%min_energy, &
        spectrum%max_energy, spectrum%energy_step, sps_batch)

      do ie = 1, no_e
        sps(ie) = (sps(ie)**2 + spc(ie)**2)
      end do

      call spectrum_hs_output(spectrum, namespace, out_file, pol, no_e, sps)   

      call acc_batch%end()
      call sps_batch%end()
      call spc_batch%end()

      SAFE_DEALLOCATE_A(racc)

    end if

    SAFE_DEALLOCATE_A(acc)
    POP_SUB(spectrum_hs_from_acc)
  end subroutine spectrum_hs_from_acc
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine spectrum_hs_from_current(spectrum, namespace, out_file, pol, vec, w0)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: out_file
    character,         intent(in)    :: pol
    FLOAT,             intent(in)    :: vec(:)
    FLOAT,   optional, intent(in)    :: w0

    integer :: istep, jj, iunit, time_steps, istart, iend, ntiter, ierr, no_e, ie
    FLOAT :: dt, cc(MAX_DIM),vv(MAX_DIM)
    CMPLX, allocatable :: cur(:)
    FLOAT, allocatable :: rcur(:), sps(:), spc(:)
    type(batch_t) :: cur_batch, sps_batch, spc_batch

    PUSH_SUB(spectrum_hs_from_current)

    call spectrum_tdfile_info(namespace, 'total_current', iunit, time_steps, dt)
    call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! load dipole from file
    SAFE_ALLOCATE(cur(0:time_steps))
    cur = M_ZERO
    vv = vec / sqrt(sum(vec(:)**2))  
    call io_skip_header(iunit)

    do istep = 1, time_steps
      cc = M_ZERO
      read(iunit, '(28x,e20.12)', advance = 'no', iostat = ierr) cc(1)
      ! FIXME: parsing of file depends on how code was compiled (MAX_DIM)!!!
      jj = 2
      do while( (ierr == 0) .and. (jj <= MAX_DIM) )
        read(iunit, '(e20.12)', advance = 'no', iostat = ierr) cc(jj)
        jj = jj + 1 
      end do
      select case(pol)
      case('x')
        cur(istep) = cc(1)
      case('y')
        cur(istep) = cc(2)
      case('z')
        cur(istep) = cc(3)
      case('+')
        cur(istep) = (cc(1) + M_zI * cc(2)) / sqrt(M_TWO)
      case('-')
        cur(istep) = (cc(1) - M_zI * cc(2)) / sqrt(M_TWO)
      case('v')
        cur(istep) = vv(1)*cc(1) + vv(2)*cc(2) + vv(3)*cc(3)
      end select
      cur(istep) = units_to_atomic(units_out%velocity, cur(istep))
    end do
    close(iunit)

    if(present(w0)) then

      call spectrum_hsfunction_init(dt, istart, iend, time_steps, cur)
      call spectrum_hs(spectrum, namespace, out_file, pol, w0)
      call spectrum_hsfunction_end()

    else

      SAFE_ALLOCATE(rcur(0:time_steps))
      rcur = real(cur, REAL_PRECISION)

      no_e = spectrum_nenergy_steps(spectrum)
      SAFE_ALLOCATE(sps(1:no_e))
      SAFE_ALLOCATE(spc(1:no_e))
      sps = M_ZERO
      spc = M_ZERO

      call batch_init(cur_batch, 1, 1)
      call batch_init(sps_batch, 1, 1)
      call batch_init(spc_batch, 1, 1)

      call cur_batch%add_state(rcur)
      call sps_batch%add_state(sps)
      call spc_batch%add_state(spc)

      call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, cur_batch, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, spc_batch)
      call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, cur_batch, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, sps_batch)

      do ie = 1, no_e
        sps(ie) = (sps(ie)**2 + spc(ie)**2) * ((ie-1) * spectrum%energy_step + spectrum%min_energy)**2
      end do

      call spectrum_hs_output(spectrum, namespace, out_file, pol, no_e, sps)   

      call cur_batch%end()
      call sps_batch%end()
      call spc_batch%end()

      SAFE_DEALLOCATE_A(rcur)

    end if

    SAFE_DEALLOCATE_A(cur)
    POP_SUB(spectrum_hs_from_current)
  end subroutine spectrum_hs_from_current


  ! ---------------------------------------------------------
  subroutine spectrum_hs(spectrum, namespace, out_file, pol, w0)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: out_file
    character,         intent(in)    :: pol
    FLOAT,  optional,  intent(in)    :: w0

    integer :: iunit, no_e, ie
    FLOAT   :: omega, hsval, xx
    FLOAT, allocatable :: sp(:)

    PUSH_SUB(spectrum_hs)

    if(present(w0)) then

      iunit = io_open(trim(out_file) // "." // trim(pol), namespace, action='write')
      write(iunit, '(a1,a20,a20)') '#', str_center("w", 20), str_center("H(w)", 20)
      write(iunit, '(a1,a20,a20)') '#', &
        str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
        str_center('[('//trim(units_abbrev(units_out%length))//'/' &
        //trim(units_abbrev(units_out%time**2)), 20)
      
      ! output
      omega = w0
      do while(omega <= spectrum%max_energy)
        call spectrum_hsfunction_min(namespace, omega - w0, omega + w0, xx, hsval)

        write(iunit, '(1x,2e20.8)') units_from_atomic(units_out%energy, xx), &
          units_from_atomic((units_out%length / units_out%time)**2, -hsval)
          
        ! 2 * w0 because we assume that there are only odd peaks.
        omega = omega + 2 * w0
      end do
      call io_close(iunit)

    else
      no_e = spectrum_nenergy_steps(spectrum)
      SAFE_ALLOCATE(sp(1:no_e))
      sp = M_ZERO

      do ie = 1, no_e
        call hsfunction((ie-1) * spectrum%energy_step + spectrum%min_energy, sp(ie))
        sp(ie) = -sp(ie)
      end do

      call spectrum_hs_output(spectrum, namespace, out_file, pol, no_e, sp)

      SAFE_DEALLOCATE_A(sp)

    end if

    POP_SUB(spectrum_hs)
  end subroutine spectrum_hs
  ! ---------------------------------------------------------


  subroutine spectrum_hs_output(spectrum, namespace, out_file, pol, no_e, sp)
    type(spectrum_t),  intent(inout) :: spectrum
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: out_file
    character,         intent(in)    :: pol
    integer,           intent(in)    :: no_e
    FLOAT,             intent(in)    :: sp(:)

    integer :: iunit, ie

    PUSH_SUB(spectrum_hs_output)

      ! output
    if(trim(out_file) /= '-') then
      iunit = io_open(trim(out_file) // "." // trim(pol), namespace, action='write')
      write(iunit, '(a1,a20,a20)') '#', str_center("w", 20), str_center("H(w)", 20)
       
      write(iunit, '(a1,a20,a20)') &
        '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
        str_center('[('//trim(units_abbrev(units_out%length))//'/' &
          //trim(units_abbrev(units_out%time**2)), 20)
        
      do ie = 1, no_e
        write(iunit, '(2e15.6)') units_from_atomic(units_out%energy, (ie-1) * spectrum%energy_step + spectrum%min_energy), &
          units_from_atomic((units_out%length / units_out%time)**2, sp(ie))
      end do
        
      call io_close(iunit)
    end if

    POP_SUB(spectrum_hs_output)
  end subroutine spectrum_hs_output


  ! ---------------------------------------------------------
  subroutine spectrum_mult_info(namespace, iunit, nspin, kick, time_steps, dt, file_units, lmax)
    type(namespace_t),   intent(in)  :: namespace
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
    call kick_read(kick, iunit, namespace)
    read(iunit, '(a)') line
    read(iunit, '(a)') line
    call io_skip_header(iunit)

    ! Figure out the units of the file
    ii = index(line,'eV')
    if(ii /= 0) then
      call unit_system_get(file_units, UNITS_EVA)
    else
      call unit_system_get(file_units, UNITS_ATOMIC)
    end if

    call spectrum_count_time_steps(namespace, iunit, time_steps, dt)
    dt = units_to_atomic(file_units%time, dt) ! units_out is OK

    POP_SUB(spectrum_mult_info)
  end subroutine spectrum_mult_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------  
  subroutine spectrum_count_time_steps(namespace, iunit, time_steps, dt)
    type(namespace_t), intent(in)  :: namespace
    integer,           intent(in)  :: iunit
    integer,           intent(out) :: time_steps
    FLOAT,             intent(out) :: dt

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
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(count_time_steps)
  end subroutine spectrum_count_time_steps
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_cross_section_info(namespace, iunit, nspin, kick, energy_steps, dw)
    type(namespace_t), intent(in)  :: namespace
    integer,           intent(in)  :: iunit
    integer,           intent(out) :: nspin
    type(kick_t),      intent(out) :: kick
    integer,           intent(out) :: energy_steps
    FLOAT,             intent(out) :: dw            !< energy step

    FLOAT :: dummy, e1, e2

    PUSH_SUB(spectrum_cross_section_info)

    ! read in number of spin components
    read(iunit, '(15x,i2)') nspin
    call kick_read(kick, iunit, namespace)
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

    if(energy_steps < 3) then
      message(1) = "Empty multipole file?"
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(spectrum_cross_section_info)
  end subroutine spectrum_cross_section_info


  ! ---------------------------------------------------------
  subroutine spectrum_tdfile_info(namespace, fname, iunit, time_steps, dt)
    type(namespace_t), intent(in)  :: namespace
    character(len=*),  intent(in)  :: fname
    integer,           intent(out) :: iunit, time_steps
    FLOAT,             intent(out) :: dt

    integer :: trash
    FLOAT :: t1, t2, dummy
    character(len=256) :: filename

    PUSH_SUB(spectrum_tdfile_info)

    ! open files
    filename = trim('td.general/')//trim(fname)
    iunit = io_open(filename, namespace, action='read', status='old', die=.false.)

    if(iunit < 0) then
      filename = trim('./')//trim(fname)
      iunit = io_open(filename, namespace, action='read', status='old')
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
      message(1) = "Empty file?"
      call messages_fatal(1, namespace=namespace)
    end if

    rewind(iunit)
    POP_SUB(spectrum_tdfile_info)
  end subroutine spectrum_tdfile_info

  
  ! ---------------------------------------------------------
  subroutine spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)
    type(spectrum_t), intent(inout) :: spectrum
    integer,          intent(in)    :: time_steps
    FLOAT,            intent(in)    :: dt
    integer,          intent(out)   :: istart, iend, ntiter

    FLOAT :: ts, te, dummy

    PUSH_SUB(spectrum_fix_time_limits)

    ts = M_ZERO
    te = time_steps * dt

    if(spectrum%start_time < ts) spectrum%start_time = ts
    if(spectrum%start_time > te) spectrum%start_time = te
    if(spectrum%end_time   > te .or. spectrum%end_time <= M_ZERO) spectrum%end_time = te
    if(spectrum%end_time   < ts) spectrum%end_time   = ts

    if(spectrum%end_time < spectrum%start_time) then
      dummy = spectrum%end_time ! swap
      spectrum%end_time = spectrum%start_time
      spectrum%start_time = dummy
    end if
    istart = int(spectrum%start_time / dt)
    iend = int(spectrum%end_time / dt)
    ntiter = iend - istart + 1

    ! Get default damp factor
    if (spectrum%damp /= SPECTRUM_DAMP_NONE .and. spectrum%damp /= SPECTRUM_DAMP_POLYNOMIAL &
         .and. spectrum%damp_factor == -M_ONE) then
      select case(spectrum%damp)
        case(SPECTRUM_DAMP_LORENTZIAN)
          spectrum%damp_factor =  -log(0.0001)/(spectrum%end_time-spectrum%start_time)
        case(SPECTRUM_DAMP_GAUSSIAN)
          spectrum%damp_factor =  sqrt(-log(0.0001)/(spectrum%end_time-spectrum%start_time)**2)
      end select
    end if


    POP_SUB(spectrum_fix_time_limits)
  end subroutine spectrum_fix_time_limits

  ! -------------------------------------------------------

  subroutine spectrum_signal_damp(damp_type, damp_factor, time_start, time_end, t0, time_step, time_function)
    integer,            intent(in)    :: damp_type
    FLOAT,              intent(in)    :: damp_factor    
    integer,            intent(in)    :: time_start
    integer,            intent(in)    :: time_end
    FLOAT,              intent(in)    :: t0
    FLOAT,              intent(in)    :: time_step
    type(batch_t),      intent(inout) :: time_function

    integer :: itime, ii
    FLOAT   :: time
    FLOAT, allocatable :: weight(:)

    PUSH_SUB(signal_damp)

    ASSERT(time_function%is_ok())
    ASSERT(time_function%status() == BATCH_NOT_PACKED)

    SAFE_ALLOCATE(weight(time_start:time_end))

    do itime = time_start, time_end
      time = time_step*(itime-1)

      ! Gets the damp function
      select case(damp_type)
      case(SPECTRUM_DAMP_NONE)
        weight(itime) = M_ONE
      case(SPECTRUM_DAMP_LORENTZIAN)
        if (time < t0) then
          weight(itime) = M_ONE
        else
          weight(itime) = exp(-(time - t0)*damp_factor)
        end if
      case(SPECTRUM_DAMP_POLYNOMIAL)
        if (time < t0) then
          weight(itime) = M_ONE
        else
          weight(itime) = M_ONE - M_THREE*( (time - t0)/(time_step*(time_end - 1) - t0) )**2 + &
               M_TWO*( (time - t0)/(time_step*(time_end - 1) - t0) )**3
        end if
      case(SPECTRUM_DAMP_GAUSSIAN)
        if (time < t0) then
          weight(itime) = M_ONE
        else
          weight(itime) = exp(-(time - t0)**2*damp_factor**2)
        end if
      end select
    end do
            
    if(time_function%type() == TYPE_CMPLX) then
      do ii = 1, time_function%nst_linear
        do itime = time_start, time_end
          time_function%states_linear(ii)%zpsi(itime) = weight(itime)*time_function%states_linear(ii)%zpsi(itime)
        end do      
      end do
    else     
      do ii = 1, time_function%nst_linear
        do itime = time_start, time_end
          time_function%states_linear(ii)%dpsi(itime) = weight(itime)*time_function%states_linear(ii)%dpsi(itime)
        end do
      end do
    end if      

    SAFE_DEALLOCATE_A(weight)

    POP_SUB(signal_damp)

  end subroutine spectrum_signal_damp
  ! -------------------------------------------------------

  ! -------------------------------------------------------
  !> Computes the sine, cosine, (or "exponential") Fourier transform of the real function given in the
  !! time_function batch. 
  !!
  !! The initial and final integration times are given by time_start and time_end,
  !! whereas the initial and final computed energies are given by energy_start and energy_end. The result
  !! is placed in the energy_function batch. The cosine Fourier function is computed by multiplying the
  !! real function by \f$ \cos(w*(t-t0)) \f$, the sine Fourier transform is computed by multiplying the real function
  !! by \f$ \sin(w*(t-t0)) \f$, and the "exponential" transform is computed by multiplying the real function by
  !! \f$ e(-I*w*t0)*e(-w*t) \f$.
  subroutine spectrum_fourier_transform(method, transform, noise, time_start, time_end, t0, time_step, time_function, &
    energy_start, energy_end, energy_step, energy_function)
    integer,                  intent(in)    :: method
    integer,                  intent(in)    :: transform
    FLOAT,                    intent(in)    :: noise
    integer,                  intent(in)    :: time_start
    integer,                  intent(in)    :: time_end
    FLOAT,                    intent(in)    :: t0
    FLOAT,                    intent(in)    :: time_step
    type(batch_t),            intent(in)    :: time_function
    FLOAT,                    intent(in)    :: energy_start
    FLOAT,                    intent(in)    :: energy_end
    FLOAT,                    intent(in)    :: energy_step
    type(batch_t),            intent(inout) :: energy_function

    integer :: itime, ienergy, ii, energy_steps
    FLOAT   :: energy, sinz, cosz
    CMPLX :: ez, eidt
    type(compressed_sensing_t) :: cs

    PUSH_SUB(fourier_transform)
    
    ASSERT(time_function%is_ok())
    ASSERT(energy_function%is_ok())
    ASSERT(time_function%nst_linear == energy_function%nst_linear)
    ASSERT(time_function%status() == energy_function%status())
    ASSERT(time_function%status() == BATCH_NOT_PACKED)
    ASSERT(time_function%type() == TYPE_FLOAT)
    ASSERT(energy_function%type() == TYPE_FLOAT)

    energy_steps = nint((energy_end-energy_start) / energy_step) + 1 

    select case(method)

    case(SPECTRUM_FOURIER)

      do ienergy = 1, energy_steps

        energy = energy_step*(ienergy - 1) + energy_start

        do ii = 1, energy_function%nst_linear
          energy_function%states_linear(ii)%dpsi(ienergy) = M_ZERO
        end do

        select case(transform)

        ! The sine and cosine transforms are computed as the real and imaginary part of the exponential.
        ! One can compute the exponential by successive multiplications, instead of calling the sine or
        ! cosine function at each time step.
        case(SPECTRUM_TRANSFORM_SIN)

          eidt = exp(M_zI * energy * time_step )
          ez = exp(M_zI * energy * ( (time_start-1)*time_step - t0) )
          sinz = aimag(ez)
          do itime = time_start, time_end
            do ii = 1, time_function%nst_linear
              energy_function%states_linear(ii)%dpsi(ienergy) = &
                energy_function%states_linear(ii)%dpsi(ienergy) + &
                  time_function%states_linear(ii)%dpsi(itime) * sinz
            end do
            ez = ez * eidt
            sinz = aimag(ez)
          end do

        case(SPECTRUM_TRANSFORM_COS)

          eidt = exp(M_zI * energy * time_step)
          ez = exp(M_zI * energy * ( (time_start-1)*time_step - t0) )
          cosz = real(ez, REAL_PRECISION)
          do itime = time_start, time_end
            do ii = 1, time_function%nst_linear
              energy_function%states_linear(ii)%dpsi(ienergy) = &
                energy_function%states_linear(ii)%dpsi(ienergy) + &
                  time_function%states_linear(ii)%dpsi(itime) * cosz
            end do
            ez = ez * eidt
            cosz = real(ez, REAL_PRECISION)
          end do

        case(SPECTRUM_TRANSFORM_LAPLACE)
        
          eidt = exp( -energy * time_step)
          ez = exp( -energy * ( (time_start-1)*time_step - t0) )
          do itime = time_start, time_end
            do ii = 1, time_function%nst_linear
              energy_function%states_linear(ii)%dpsi(ienergy) = &
                energy_function%states_linear(ii)%dpsi(ienergy) + &
                real( time_function%states_linear(ii)%dpsi(itime) * ez, REAL_PRECISION)
            end do
            ez = ez * eidt
          end do
        end select

        ! The total sum must be multiplied by time_step in order to get the integral.
        do ii = 1, time_function%nst_linear
            energy_function%states_linear(ii)%dpsi(ienergy) = &
              energy_function%states_linear(ii)%dpsi(ienergy) * time_step
        end do
        

      end do

    case(SPECTRUM_COMPRESSED_SENSING)

      call compressed_sensing_init(cs, transform, &
        time_end - time_start + 1, time_step, time_step*(time_start - 1) - t0, &
        energy_steps, energy_step, energy_start, noise)

      do ii = 1, time_function%nst_linear
        call compressed_sensing_spectral_analysis(cs, time_function%states_linear(ii)%dpsi, &
          energy_function%states_linear(ii)%dpsi)
      end do

      call compressed_sensing_end(cs)

    end select

    POP_SUB(fourier_transform)

  end subroutine spectrum_fourier_transform

  ! ---------------------------------------------------------
  subroutine spectrum_sigma_diagonalize(namespace, sigma, nspin, energy_step, min_energy, energy_steps, kick)
    type(namespace_t),      intent(in) :: namespace
    FLOAT,                  intent(in) :: sigma(:, :, :, :) !< (3, 3, energy_steps, nspin) already converted to units
    integer,                intent(in) :: nspin
    FLOAT,                  intent(in) :: energy_step, min_energy
    integer,                intent(in) :: energy_steps
    type(kick_t), optional, intent(in) :: kick !< if present, will write itself and nspin

    integer :: is, idir, jdir, ie, info, out_file, out_file_t
    FLOAT, allocatable :: work(:,:) 
    CMPLX, allocatable :: w(:)
    character(len=20) :: header_string
    logical :: spins_singlet, spins_triplet, symmetrize
    FLOAT, allocatable :: pp(:,:), pp2(:,:)

    PUSH_SUB(spectrum_sigma_diagonalize)

    
    !%Variable PropagationSpectrumSymmetrizeSigma
    !%Type logical
    !%Default .false.
    !%Section Utilities::oct-propagation_spectrum
    !%Description
    !% The polarizablity tensor has to be real and symmetric. Due to numerical accuracy, 
    !% that is not extricly conserved when computing it from different time-propations.
    !% If <tt>PropagationSpectrumSymmetrizeSigma = yes</tt>, the polarizability tensor is
    !% symmetrized before its diagonalizied.
    !% This variable is only used if the cross_section_tensor is computed. 
    !%End
    call parse_variable(namespace, 'PropagationSpectrumSymmetrizeSigma', .false., symmetrize)
    call messages_print_var_value(stdout, 'PropagationSpectrumSymmetrizeSigma', symmetrize)

    spins_singlet = .true.
    spins_triplet = .false.
    if(present(kick)) then
      select case(kick%delta_strength_mode)
      case (KICK_SPIN_MODE)
        spins_triplet = .true.
        spins_singlet = .false.
      case (KICK_SPIN_DENSITY_MODE)
        spins_triplet = .true.
      end select
    end if
    
    if (spins_singlet .and. spins_triplet) then
      out_file = io_open('cross_section_diagonal-sigma_s', namespace, action='write')
      out_file_t = io_open('cross_section_diagonal-sigma_t', namespace, action='write')
    else
      out_file = io_open('cross_section_diagonal-sigma', namespace, action='write')
    end if

    write(out_file, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
    do idir = 1, 3
      write(out_file, '(a20)', advance = 'no') str_center("Real part", 20)
      if (.not.symmetrize) write(out_file, '(a20)', advance = 'no') str_center("Imaginary part", 20)
      do jdir = 1, 3
        write(header_string,'(a7,i1,a1,i1,a1,i1,a1)') 'vector(', idir, ',', jdir, ',', is, ')'
        write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
      end do
    end do
    write(out_file, '(1x)')
    write(out_file, '(a1,a20)', advance = 'no') '#', str_center('[' // trim(units_abbrev(units_out%energy)) // ']', 20)

    do idir = 1, 3
      write(out_file, '(a20)', advance = 'no')  str_center('[' // trim(units_abbrev(units_out%length**2)) // ']', 20)
      if (.not.symmetrize) then
        write(out_file, '(a20)', advance = 'no')  str_center('[' // trim(units_abbrev(units_out%length**2)) // ']', 20)
      end if
      do jdir = 1, 3
        write(out_file, '(a20)', advance = 'no')  str_center('[ - ]', 20)
      end do
    end do
    write(out_file, '(1x)')

    if (spins_singlet .and. spins_triplet) then
      write(out_file_t, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
      do idir = 1, 3
        write(out_file_t, '(a20)', advance = 'no') str_center("Real part", 20)
        if (.not.symmetrize) write(out_file_t, '(a20)', advance = 'no') str_center("Imaginary part", 20)
        do jdir = 1, 3
          write(header_string,'(a7,i1,a1,i1,a1,i1,a1)') 'vector(', idir, ',', jdir, ',', is, ')'
          write(out_file_t, '(a20)', advance = 'no') str_center(trim(header_string), 20)
        end do
      end do
      write(out_file_t, '(1x)')
      write(out_file_t, '(a1,a20)', advance = 'no') '#', str_center('[' // trim(units_abbrev(units_out%energy)) // ']', 20)
     
      do idir = 1, 3
        write(out_file_t, '(a20)', advance = 'no')  str_center('[' // trim(units_abbrev(units_out%length**2)) // ']', 20)
        if(.not.symmetrize) then
          write(out_file_t, '(a20)', advance = 'no')  str_center('[' // trim(units_abbrev(units_out%length**2)) // ']', 20)
        end if
        do jdir = 1, 3
          write(out_file_t, '(a20)', advance = 'no')  str_center('[ - ]', 20)
        end do
      end do
      write(out_file_t, '(1x)')
    end if

    SAFE_ALLOCATE(pp(1:3, 1:3))
    if (spins_triplet .and. spins_singlet) SAFE_ALLOCATE(pp2(1:3, 1:3))
    SAFE_ALLOCATE(w(1:3))
    SAFE_ALLOCATE(work(1:3, 1:3))
    do ie = 1, energy_steps

      pp(:, :) = sigma(:, :, ie, 1)
      if (nspin >= 2) then
        if (spins_singlet .and. spins_triplet) then
          pp2(:, :) = pp(:, :) - sigma(:, :, ie, 2)
          pp(:, :)  = pp(:, :) + sigma(:, :, ie, 2)
        elseif (spins_triplet .and. .not.spins_singlet) then
          pp(:, :) = pp(:, :) - sigma(:, :, ie, 2)
        elseif (spins_singlet .and. .not.spins_triplet) then
          pp(:, :) = pp(:, :) + sigma(:, :, ie, 2)
        end if
      end if

      if (symmetrize) then
        do idir = 1, 3
          do jdir = idir + 1, 3
            pp(idir, jdir) = (pp(idir, jdir) + pp(jdir, idir) )/2.
            pp(jdir, idir) = pp(idir, jdir)
          end do 
        end do
      end if

      work(1:3, 1:3) = pp(1:3, 1:3)
      call lalg_eigensolve_nonh(3, work, w, err_code = info, sort_eigenvectors = .true.)
      ! Note that the cross-section elements do not have to be transformed to the proper units, since
      ! they have been read from the "cross_section_vector.x", where they are already in the proper units.

      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, ((ie-1) * energy_step + min_energy))
      do idir = 3, 1, -1
        if (symmetrize) then
          write(out_file,'(2e20.8)', advance = 'no') real(w(idir))
        else
          write(out_file,'(2e20.8)', advance = 'no') w(idir)
        end if
      
        do jdir = 1, 3
          write(out_file,'(e20.8)', advance = 'no') work(jdir, idir)
        end do
      end do 
      write(out_file, '(1x)')

      if (spins_singlet .and. spins_triplet) then
        if (symmetrize) then
          do idir = 1, 3
            do jdir = idir + 1, 3
              pp2(idir, jdir) = (pp2(idir, jdir) + pp2(jdir, idir) )/2.
              pp2(jdir, idir) = pp2(idir, jdir)
            end do 
          end do
        end if
        work(1:3, 1:3) = -pp2(1:3, 1:3)
        call lalg_eigensolve_nonh(3, work, w, err_code = info, sort_eigenvectors = .true.)
        ! Note that the cross-section elements do not have to be transformed to the proper units, since
        ! they have been read from the "cross_section_vector.x", where they are already in the proper units.
      
        write(out_file_t,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, (ie * energy_step + min_energy))
        do idir = 3, 1, -1
          if (symmetrize) then
            write(out_file_t,'(2e20.8)', advance = 'no') real(w(idir))
          else
            write(out_file_t,'(2e20.8)', advance = 'no') w(idir)
          end if
        
          do jdir = 1, 3
            write(out_file_t,'(e20.8)', advance = 'no') work(jdir, idir)
          end do
        end do 
        write(out_file_t, '(1x)')
      end if
    end do

    call io_close(out_file)

    SAFE_DEALLOCATE_A(pp)
    if (spins_triplet .and. spins_singlet) then 
      SAFE_DEALLOCATE_A(pp2)
      call io_close(out_file_t)
    end if
    SAFE_DEALLOCATE_A(w)
    SAFE_DEALLOCATE_A(work)

    POP_SUB(spectrum_sigma_diagonalize)
  end subroutine spectrum_sigma_diagonalize

  pure integer function spectrum_nenergy_steps(spectrum) result(no_e)
    type(spectrum_t), intent(in) :: spectrum

    no_e = nint((spectrum%max_energy-spectrum%min_energy) / spectrum%energy_step) + 1
  end function spectrum_nenergy_steps

  subroutine spectrum_write_info(spectrum, out_file)
    type(spectrum_t), intent(in) :: spectrum
    integer,          intent(in) :: out_file

    PUSH_SUB(spectrum_write_info)

    write(out_file, '(a,i4)')    '# PropagationSpectrumDampMode   = ', spectrum%damp
    write(out_file, '(a,f10.4)') '# PropagationSpectrumDampFactor = ', units_from_atomic(units_out%time**(-1), &
                                                                       spectrum%damp_factor)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumStartTime  = ', units_from_atomic(units_out%time, spectrum%start_time)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumEndTime    = ', units_from_atomic(units_out%time, spectrum%end_time)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumMinEnergy  = ', units_from_atomic(units_out%energy, spectrum%min_energy)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumMaxEnergy  = ', units_from_atomic(units_out%energy, spectrum%max_energy)
    write(out_file, '(a,f10.4)') '# PropagationSpectrumEnergyStep = ', units_from_atomic(units_out%energy, spectrum%energy_step)

    POP_SUB(spectrum_write_info)
  end subroutine spectrum_write_info

end module spectrum_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
