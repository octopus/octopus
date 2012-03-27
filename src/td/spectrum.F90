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
  use fft_m
  use global_m
  use io_m
  use kick_m
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
    spectrum_init,                 &
    spectrum_cross_section,        &
    spectrum_cross_section_tensor, &
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
    count_time_steps,              &
    signal_damp,                   &
    fourier_transform

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
    SPECTRUM_ABSORPTION      = 1,  &
    SPECTRUM_ENERGYLOSS      = 2

  integer, public, parameter ::       &
    SPECTRUM_FOURIER            = 1,  &
    SPECTRUM_COMPRESSED_SENSING = 2

  type spec_t
    FLOAT   :: start_time          !< start time for the transform
    FLOAT   :: end_time            !< when to stop the transform
    FLOAT   :: energy_step         !< step in energy mesh
    FLOAT   :: max_energy          !< maximum of energy mesh
    integer :: damp                !< damping type (none, exp or pol)
    integer :: transform           !< sine, cosine, or exponential transform
    FLOAT   :: damp_factor         !< factor used in damping
    integer :: spectype            !< damping type (none, exp or pol)
    integer :: method              !< fourier transform or compressed sensing 
    FLOAT   :: noise               !< the level of noise that is assumed in the time series for compressed sensing 
  end type spec_t

  !> Module variables, necessary to compute the function hsfunction, called by
  !! the C function loct_1dminimize
  integer :: niter_
  FLOAT :: time_step_, energy_step_
  CMPLX, allocatable :: func_(:),func_ar_(:,:),pos_(:,:),tret_(:), funcw_(:)
  type(fft_t) :: fft_handler
  CMPLX :: vv_(MAX_DIM)
  integer :: is_, ie_, default

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
  subroutine spectrum_cross_section(in_file, out_file, spectrum, ref_file)
    integer,           intent(in)    :: in_file
    integer,           intent(in)    :: out_file
    type(spec_t),      intent(inout) :: spectrum
    integer, optional, intent(in)    :: ref_file

    character(len=20) :: header_string
    integer :: nspin, ref_nspin, lmax, ref_lmax, time_steps, &
               ref_time_steps, istart, iend, ntiter, it, ii, isp, no_e, ie, idir, trash
    FLOAT   :: dt, ref_dt, dump, energy, ewsum, polsum
    type(kick_t) :: kick, ref_kick
    FLOAT, allocatable :: dipole(:, :, :), ref_dipole(:, :, :), sigma(:, :, :), sf(:, :)
    type(unit_system_t) :: file_units, ref_file_units
    type(batch_t) :: dipoleb, sigmab

    PUSH_SUB(spectrum_cross_section)

    ! This function gives us back the unit connected to the "multipoles" file, the header information,
    ! the number of time steps, and the time step.
    call spectrum_mult_info(in_file, nspin, kick, time_steps, dt, file_units, lmax=lmax)

    if(present(ref_file)) then
      call spectrum_mult_info(ref_file, ref_nspin, ref_kick, &
        ref_time_steps, ref_dt, ref_file_units, lmax = ref_lmax)
      if( (nspin .ne. ref_nspin)           .or. &
          (time_steps .ne. ref_time_steps) .or. &
          (.not.(dt .app. ref_dt))         .or. &
          (lmax .ne. ref_lmax) ) then
        write(message(1),'(a)') 'The multipoles and reference multipoles files do not match.'
        call messages_fatal(1)
      end if
    end if

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

    if(present(ref_file)) then
      call io_skip_header(ref_file)
      SAFE_ALLOCATE(ref_dipole(0:time_steps, 1:3, 1:nspin))
      do it = 0, time_steps
        select case(nspin)
        case(1)
          read(ref_file, *) trash, dump, dump, ref_dipole(it, 1:3, 1)
        case(2)
          read(ref_file, *) trash, dump, dump, ref_dipole(it, 1:3, 1), dump, ref_dipole(it, 1:3, 2)
        case(4)
          read(ref_file, *) &
            trash, dump, dump, ref_dipole(it, 1:3, 1), dump, ref_dipole(it, 1:3, 2), &
            dump, ref_dipole(it, 1:3, 3), dump, ref_dipole(it, 1:3, 4)
        end select
        ref_dipole(it, 1:3, :) = units_to_atomic(file_units%length, ref_dipole(it, 1:3, :))
      end do
    end if

    ! Now subtract the initial dipole.
    if(present(ref_file)) then
      dipole = dipole - ref_dipole
    else
      do it = time_steps, 0, -1
        dipole(it, :, :) = dipole(it, :, :) - dipole(0, :, :)
      end do
    end if

    if(spectrum%energy_step <= M_ZERO) spectrum%energy_step = M_TWO * M_PI / (dt*time_steps)

    ! Get the number of energy steps.
    no_e = spectrum%max_energy / spectrum%energy_step
    SAFE_ALLOCATE(sigma(0:no_e, 1:3, 1:nspin))

    call batch_init(dipoleb, 3, 1, nspin, dipole)
    call batch_init(sigmab, 3, 1, nspin, sigma)

    call signal_damp(spectrum%damp, spectrum%damp_factor, istart + 1, iend + 1, dt, dipoleb)
    call fourier_transform(spectrum%method, spectrum%transform, spectrum%noise, &
      istart + 1, iend + 1, kick%time, dt, dipoleb, 1, no_e + 1, spectrum%energy_step, sigmab)

    call batch_end(dipoleb)
    call batch_end(sigmab)

    SAFE_DEALLOCATE_A(dipole)
    if(present(ref_file)) then
      SAFE_DEALLOCATE_A(ref_dipole)
    end if

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
      write(out_file, '(a20)', advance = 'no') str_center('[' // trim(units_abbrev(unit_one/units_out%energy)) // ']', 20)
    end do
    write(out_file, '(1x)')

    do ie = 0, no_e
      write(out_file,'(e20.8)', advance = 'no') units_from_atomic(units_out%energy, ie * spectrum%energy_step)
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
      ftchd(it) = TOCMPLX(dummy3-dummy2, dummy4+dummy1)
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
      istart + 1, iend + 1, kick%time, dt, angularb, 1, no_e + 1, spectrum%energy_step, respb)
    call fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
      istart + 1, iend + 1, kick%time, dt, angularb, 1, no_e + 1, spectrum%energy_step, imspb)

    call batch_end(angularb)
    call batch_end(respb)
    call batch_end(imspb)
    
    sum1 = M_Z0
    sum2 = M_Z0
    do ie = 0, no_e
      energy = ie * spectrum%energy_step

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
  subroutine spectrum_hsfunction_init(dt, is, ie, niter, acc)
    FLOAT,   intent(in)           :: dt
    integer, intent(in)           :: is, ie, niter
    CMPLX,   intent(in)           :: acc(:)

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
  subroutine spectrum_hsfunction_min(aa, bb, omega, omega_min, func_min)
    FLOAT, intent(in)  :: aa, bb, omega
    FLOAT, intent(out) :: omega_min, func_min

    integer :: ierr, ie
    FLOAT :: xx, hsval, minhsval, ww, xa, xb, hxa, hxb

    PUSH_SUB(spectrum_hsfunction_min)

    ! xx should be an initial guess for the minimum. So we do a quick search
    ! that we refine later calling 1dminimize.
    !xx = omega
    !call hsfunction(xx, minhsval)

    ie = int(aa/energy_step_)
    ww = ie * energy_step_
    if(ww < aa) then
      ie = ie + 1
      ww = ie * energy_step_
    end if
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
#ifndef SINGLE_PRECISION
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

    integer :: istep, trash, iunit, nspin, time_steps, istart, iend, ntiter, lmax, no_e, ie
    FLOAT :: dt, dump, vv(MAX_DIM)  
    type(kick_t) :: kick
    FLOAT, allocatable :: dd(:,:)
    FLOAT, allocatable :: sps(:), spc(:), racc(:)
    CMPLX, allocatable :: dipole(:), ddipole(:)
    type(batch_t) :: acc_batch, sps_batch, spc_batch
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

    if(present(w0)) then

      call spectrum_hsfunction_init(dt, istart, iend, time_steps, ddipole)
      call spectrum_hs(out_file, spectrum, pol, w0)
      call spectrum_hsfunction_end()

    else

      SAFE_ALLOCATE(racc(0:time_steps))
      racc = ddipole

      no_e = spectrum%max_energy / spectrum%energy_step
      SAFE_ALLOCATE(sps(0:no_e))
      SAFE_ALLOCATE(spc(0:no_e))
      sps = M_ZERO
      spc = M_ZERO

      call batch_init(acc_batch, 1)
      call batch_init(sps_batch, 1)
      call batch_init(spc_batch, 1)

      call batch_add_state(acc_batch, racc)
      call batch_add_state(sps_batch, sps)
      call batch_add_state(spc_batch, spc)

      call fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, acc_batch, 1, no_e + 1, spectrum%energy_step, spc_batch)
      call fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, acc_batch, 1, no_e + 1, spectrum%energy_step, sps_batch)

      do ie = 0, no_e
        sps(ie) = (sps(ie)**2 + spc(ie)**2)
      end do

      call spectrum_hs_output(out_file, spectrum, pol, no_e, sps)   

      call batch_end(acc_batch)
      call batch_end(sps_batch)
      call batch_end(spc_batch)

      SAFE_DEALLOCATE_A(racc)

    end if

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

    integer :: istep, jj, iunit, time_steps, istart, iend, ntiter, ierr, no_e, ie
    FLOAT :: dt, aa(MAX_DIM),vv(MAX_DIM)
    CMPLX, allocatable :: acc(:)
    FLOAT, allocatable :: racc(:), sps(:), spc(:)
    type(batch_t) :: acc_batch, sps_batch, spc_batch

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

    if(present(w0)) then

      call spectrum_hsfunction_init(dt, istart, iend, time_steps, acc)
      call spectrum_hs(out_file, spectrum, pol, w0)
      call spectrum_hsfunction_end()

    else

      SAFE_ALLOCATE(racc(0:time_steps))
      racc = acc

      no_e = spectrum%max_energy / spectrum%energy_step
      SAFE_ALLOCATE(sps(0:no_e))
      SAFE_ALLOCATE(spc(0:no_e))
      sps = M_ZERO
      spc = M_ZERO

      call batch_init(acc_batch, 1)
      call batch_init(sps_batch, 1)
      call batch_init(spc_batch, 1)

      call batch_add_state(acc_batch, racc)
      call batch_add_state(sps_batch, sps)
      call batch_add_state(spc_batch, spc)

      call fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, acc_batch, 1, no_e + 1, spectrum%energy_step, spc_batch)
      call fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
        istart + 1, iend + 1, M_ZERO, dt, acc_batch, 1, no_e + 1, spectrum%energy_step, sps_batch)

      do ie = 0, no_e
        sps(ie) = (sps(ie)**2 + spc(ie)**2)
      end do

      call spectrum_hs_output(out_file, spectrum, pol, no_e, sps)   

      call batch_end(acc_batch)
      call batch_end(sps_batch)
      call batch_end(spc_batch)

      SAFE_DEALLOCATE_A(racc)

    end if

    SAFE_DEALLOCATE_A(acc)
    POP_SUB(spectrum_hs_from_acc)
  end subroutine spectrum_hs_from_acc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine spectrum_hs_fourier_transform(spectrum, no_e, time_steps, acc, sp)
    type(spec_t),     intent(inout) :: spectrum
    integer,          intent(in)    :: no_e, time_steps
    FLOAT,            intent(in)    :: acc(0:time_steps)
    FLOAT,            intent(inout) :: sp(0:no_e)

    PUSH_SUB(spectrum_hs_fourier_transform)

    POP_SUB(spectrum_hs_fourier_transform)
  end subroutine spectrum_hs_fourier_transform
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
        call spectrum_hsfunction_min(omega - w0, omega + w0, omega, xx, hsval)

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

      call spectrum_hs_output(out_file, spectrum, pol, no_e, sp)

      SAFE_DEALLOCATE_A(sp)

    end if

    POP_SUB(spectrum_hs)
  end subroutine spectrum_hs
  ! ---------------------------------------------------------


  subroutine spectrum_hs_output(out_file, spectrum, pol, no_e, sp)
    character(len=*), intent(in)    :: out_file
    type(spec_t),     intent(inout) :: spectrum
    character,        intent(in)    :: pol
    integer,          intent(in)    :: no_e
    FLOAT,            intent(in)    :: sp(0:no_e)

    integer :: iunit, ie

    PUSH_SUB(spectrum_hs_output)

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

    POP_SUB(spectrum_hs_output)
  end subroutine spectrum_hs_output


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
  subroutine fourier_transform(method, transform, noise, time_start, time_end, t0, time_step, time_function, &
    energy_start, energy_end, energy_step, energy_function)
    integer,         intent(in)    :: method
    integer,         intent(in)    :: transform
    FLOAT,           intent(in)    :: noise
    integer,         intent(in)    :: time_start
    integer,         intent(in)    :: time_end
    FLOAT,           intent(in)    :: t0
    FLOAT,           intent(in)    :: time_step
    type(batch_t),   intent(in)    :: time_function
    integer,         intent(in)    :: energy_start
    integer,         intent(in)    :: energy_end
    FLOAT,           intent(in)    :: energy_step
    type(batch_t),   intent(inout) :: energy_function

    integer :: itime, ienergy, ii
    FLOAT   :: time, energy!, kernel
    CMPLX :: ez, eidt
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

        select case(transform)

        ! The sine and cosine transforms are computed as the real and imaginary part of the exponential.
        ! One can compute the exponential by successive multiplications, instead of calling the sine or
        ! cosine function at each time step.
        case(SPECTRUM_TRANSFORM_SIN)

          eidt = exp(M_zI * energy * time_step )
          ez = exp(-M_zI * energy * t0)
          do itime = time_start, time_end
            time = time_step*(itime - time_start)
            do ii = 1, time_function%nst_linear
              energy_function%states_linear(ii)%dpsi(ienergy) = &
                energy_function%states_linear(ii)%dpsi(ienergy) + &
                aimag( time_function%states_linear(ii)%dpsi(itime) * ez)
            end do
            ez = ez * eidt
          end do

        case(SPECTRUM_TRANSFORM_COS)

          eidt = exp(M_zI * energy * time_step)
          ez = exp(-M_zI * energy * t0)
          do itime = time_start, time_end
            time = time_step*(itime - time_start)
            do ii = 1, time_function%nst_linear
              energy_function%states_linear(ii)%dpsi(ienergy) = &
                energy_function%states_linear(ii)%dpsi(ienergy) + &
                real( time_function%states_linear(ii)%dpsi(itime) * ez, REAL_PRECISION)
            end do
            ez = ez * eidt
          end do

        case(SPECTRUM_TRANSFORM_EXP)

          eidt = exp( -energy * time_step)
          ez = exp(-M_zI * energy * t0)
          do itime = time_start, time_end
            time = time_step*(itime - time_start)
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
        time_end -time_start + 1, time_step, time_step*(time_start - 1) - t0, &
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
