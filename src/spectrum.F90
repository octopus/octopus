
!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module spectrum
use string
use global
use messages
use lib_oct_parser
use syslabels
use io
use units
use lib_adv_alg

implicit none

private
public :: spec_type,                     &
          kick_type,                     &
          spec_sh,                       &
          spectrum_init,                 &
          spectrum_cross_section,        &
          spectrum_cross_section_tensor, &
          spectrum_rotatory_strength,    &
          spectrum_hs_from_mult,         &
          spectrum_hs_from_acc,          &
          spectrum_mult_info,            &
          kick_init

integer, public, parameter :: SPECTRUM_DAMP_NONE       = 0, &
                              SPECTRUM_DAMP_LORENTZIAN = 1, &
                              SPECTRUM_DAMP_POLYNOMIAL = 2, &
                              SPECTRUM_DAMP_GAUSSIAN   = 3

integer, public, parameter :: KICK_DENSITY_MODE      = 0, &
                              KICK_SPIN_MODE         = 1, &
                              KICK_SPIN_DENSITY_MODE = 2

type spec_type
  FLOAT   :: start_time  ! start time for the transform
  FLOAT   :: end_time    ! when to stop the transform
  FLOAT   :: energy_step ! step in energy mesh
  FLOAT   :: max_energy  ! maximum of energy mesh
  integer :: damp     ! Damp type (none, exp or pol)
  FLOAT   :: damp_factor ! factor used in damping
end type spec_type

type kick_type
  FLOAT             :: pol(3, 3)
  integer           :: pol_dir
  integer           :: delta_strength_mode
  FLOAT             :: delta_strength
  integer           :: pol_equiv_axis
  FLOAT             :: wprime(3)
end type kick_type

type spec_sh
  ! input
  character :: pol

  ! output
  integer :: no_e ! dimensions of sp
  FLOAT, pointer :: sp(:) ! do not forget to deallocate this
end type spec_sh

contains

subroutine spectrum_init(s)
  type(spec_type), intent(inout) :: s

  call push_sub('spectrum.spectrum_init')

  !%Variable SpecDampMode
  !%Type integer
  !%Section 13 SpectrumCalculations
  !%Description
  !% Decides which damping/filtering is to be applied in order to calculate
  !% spectra by calculating a Fourier transform
  !%Option no 0
  !% No filtering at all.
  !%Option exponential 1
  !% Exponential filtering, corresponding with a Lorentzian-shaped spectrum
  !%Option polynomial 2
  !% Third-order polynomial damping.
  !%Option gaussian 3
  !% Gaussian damping
  !%End

  call loct_parse_float(check_inp('SpecStartTime'),  M_ZERO,      s%start_time)
  call loct_parse_float(check_inp('SpecEndTime'),   -M_ONE,       s%end_time)
  call loct_parse_float(check_inp('SpecEnergyStep'), CNST(0.01),  s%energy_step)
  call loct_parse_float(check_inp('SpecMaxEnergy'),  CNST(20.0),  s%max_energy)
  call loct_parse_int  (check_inp('SpecDampMode'), SPECTRUM_DAMP_POLYNOMIAL, s%damp)
  call loct_parse_float(check_inp('SpecDampFactor'),  CNST(0.15), s%damp_factor)
  s%start_time      = s%start_time      * units_inp%time%factor
  s%end_time        = s%end_time        * units_inp%time%factor
  s%energy_step     = s%energy_step     * units_inp%energy%factor
  s%max_energy      = s%max_energy      * units_inp%energy%factor
  s%damp_factor     = s%damp_factor     / units_inp%time%factor

  call pop_sub()
end subroutine spectrum_init

subroutine kick_init(k, nspin)
  type(kick_type), intent(out) :: k
  integer,         intent(in)  :: nspin
  integer(POINTER_SIZE) :: blk
  integer :: n, i, j

  call push_sub('spectrum.kick_init')

  call loct_parse_float(check_inp('TDDeltaStrength'), M_ZERO, k%delta_strength)
  ! units are 1/length
  k%delta_strength = k%delta_strength / units_inp%length%factor

  if(k%delta_strength <= M_ZERO) then
     k%delta_strength_mode = 0
     k%pol_equiv_axis = 0
     k%pol(1:3, 1) = (/ M_ONE, M_ZERO, M_ZERO /)
     k%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
     k%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
     k%pol_dir = 0
     k%wprime = M_ZERO
     call pop_sub()
     return
  endif

  !%Variable TDDeltaStrengthMode
  !%Type integer
  !%Section 10 Time Dependent
  !%Description
  !% When calculating the linear response of the density via the propagation
  !% in real time, one needs to perfrom an initical kick on the KS system, at 
  !% time zero. Depending on what kind response property one wants to obtain,
  !% this kick may be done in several modes.
  !%Option kick_density 0
  !% The total density of the system is perturbed.
  !%Option kick_spin 1
  !% The individual spin densities are perturbed differently. Note that this mode
  !% is only possible if the run is done in spin polarized mode, or with spinors.
  !%Option kick_spin_and_density 2
  !% A combination of the two above. Note that this mode
  !% is only possible if the run is done in spin polarized mode, or with spinors.
  !%End
  call loct_parse_int(check_inp('TDDeltaStrengthMode'), KICK_DENSITY_MODE, k%delta_strength_mode)
  select case (k%delta_strength_mode)
    case (KICK_DENSITY_MODE)
    case (KICK_SPIN_MODE, KICK_SPIN_DENSITY_MODE)
      if (nspin == 1) call input_error('TDDeltaStrengthMode')
    case default
      call input_error('TDDeltaStrengthMode')
  end select

  ! Find out how many equivalent axis we have...
  ! WARNING: TODO: document this variable.
  call loct_parse_int(check_inp('TDPolarizationEquivAxis'), 0, k%pol_equiv_axis)
  call loct_parse_int(check_inp('TDPolarizationDirection'), 1, k%pol_dir)

  k%pol(:, :) = M_ZERO
  if(loct_parse_block(check_inp('TDPolarization'), blk)==0) then
        n = loct_parse_block_n(blk)
        do j = 1, n
           do i = 1, 3
              call loct_parse_block_float(blk, j-1, i-1, k%pol(i, j))
           end do
        enddo
        if(n<3) k%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
        if(n<2) k%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
        call loct_parse_block_end(blk)
  else
        ! Here the symmetry of the system should be analized, and the polarization
        ! basis, built accordingly.
        k%pol_dir = 1
        k%pol(1:3, 1) = (/ M_ONE, M_ZERO, M_ZERO /)
        k%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
        k%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
  endif
  ! Normalize:
  do i = 1, 3
     k%pol(1:3, i) = k%pol(1:3, i)/sqrt(sum(k%pol(1:3, i)**2))
  enddo

  if(loct_parse_block(check_inp('TDPolarizationWprime'), blk)==0) then
    do i = 1, 3
       call loct_parse_block_float(blk, 0, i-1, k%wprime(i))
    enddo
    k%wprime(1:3) = k%wprime(1:3)/sqrt(sum(k%wprime(1:3)**2))
  else
    k%wprime(1:3) = (/ M_ONE, M_ZERO, M_ZERO /)
  endif

  call pop_sub()
end subroutine kick_init

subroutine spectrum_skip_header(iunit)
  integer, intent(in) :: iunit
  character(len=1) :: a

  rewind(iunit)
  read(iunit,'(a)') a
  do while(a=='#')
     read(iunit,'(a)') a
  enddo
  backspace(iunit)

end subroutine spectrum_skip_header

subroutine spectrum_cross_section_tensor(s, out_file, in_file)
  type(spec_type),  intent(inout) :: s
  integer, intent(in)    :: out_file
  integer, intent(in)    :: in_file(:)

  character(len=20) :: header_string
  integer :: nspin, energy_steps, i, is, j, equiv_axis, n_files, k
  FLOAT, allocatable :: sigma(:, :, :, :), sigmap(:, :, :, :), sigmau(:, :, :),  &
                        sigmav(:, :, :), sigmaw(:, :, :), p(:, :), ip(:, :)
  FLOAT :: dw, dump, average, anisotropy
  type(kick_type) :: kick

  call push_sub('spectrum.spectrum_cross_section_tensor')

  n_files = size(in_file)
  select case(n_files)
    case(1); equiv_axis = 3
    case(2); equiv_axis = 2
    case(3); equiv_axis = 1
  end select

  call spectrum_cross_section_info(in_file(1), nspin, kick, energy_steps, dw)
  call spectrum_skip_header(in_file(1))

  allocate(sigma (3, 3, 0:energy_steps, nspin), &
           sigmap(3, 3, 0:energy_steps, nspin), &
           sigmau(3, 0:energy_steps, nspin), &
           sigmav(3, 0:energy_steps, nspin), &
           sigmaw(3, 0:energy_steps, nspin), &
           p(3, 3), ip(3, 3))

  select case(equiv_axis)

  case(3)

      do i = 0, energy_steps
         read(in_file(1), *) dump, sigmau(1:3, i, 1:nspin)
      enddo

      ! The first row of sigma is the vector that we have just read, but properly projected...
      do is = 1, nspin
         do i = 0, energy_steps
            sigmap(1, 1, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 1))
            sigmap(1, 2, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 2))
            sigmap(1, 3, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 3))
         enddo
      enddo

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
         enddo
      enddo

  case(2)

      call spectrum_cross_section_info(in_file(2), i, kick, j, dump)
      call spectrum_skip_header(in_file(2))

      do i = 0, energy_steps
         read(in_file(1), *) dump, sigmau(1:3, i, 1:nspin)
         read(in_file(2), *) dump, sigmaw(1:3, i, 1:nspin)
      enddo

      ! The first row of sigma is the vector that we have just read, but properly projected...
      do is = 1, nspin
         do i = 0, energy_steps
            sigmap(1, 1, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 1))
            sigmap(1, 2, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 2))
            sigmap(1, 3, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 3))
         enddo
      enddo

      ! The third row of sigma is also the vector that we have just read, but properly projected...
      do is = 1, nspin
         do i = 0, energy_steps
            sigmap(3, 1, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 1))
            sigmap(3, 2, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 2))
            sigmap(3, 3, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 3))
         enddo
      enddo

      ! The diagonal (2,2) is equal by symmetry to the (1,1)
      sigmap(2, 2, :, :) = sigmap(1, 1, :, :)

      ! The (2,1) term and (1,2) term are equal; the (2,3) and (3,2), also.
      sigmap(2, 1, :, :) = sigmap(1, 2, :, :)
      sigmap(2, 3, :, :) = sigmap(3, 2, :, :)

  case default

      call spectrum_cross_section_info(in_file(2), i, kick, j, dump)
      call spectrum_cross_section_info(in_file(3), i, kick, j, dump)
      call spectrum_skip_header(in_file(2))
      call spectrum_skip_header(in_file(3))

      do i = 0, energy_steps
         read(in_file(1), *) dump, sigmau(1:3, i, 1:nspin)
         read(in_file(2), *) dump, sigmav(1:3, i, 1:nspin)
         read(in_file(3), *) dump, sigmaw(1:3, i, 1:nspin)
      enddo

      do is = 1, nspin
         do i = 0, energy_steps
            sigmap(1, 1, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 1))
            sigmap(1, 2, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 2))
            sigmap(1, 3, i, is) = sum( sigmau(1:3, i, is)*kick%pol(1:3, 3))
         enddo
      enddo
      do is = 1, nspin
         do i = 0, energy_steps
            sigmap(2, 1, i, is) = sum( sigmav(1:3, i, is)*kick%pol(1:3, 1))
            sigmap(2, 2, i, is) = sum( sigmav(1:3, i, is)*kick%pol(1:3, 2))
            sigmap(2, 3, i, is) = sum( sigmav(1:3, i, is)*kick%pol(1:3, 3))
         enddo
      enddo
      do is = 1, nspin
         do i = 0, energy_steps
            sigmap(3, 1, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 1))
            sigmap(3, 2, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 2))
            sigmap(3, 3, i, is) = sum( sigmaw(1:3, i, is)*kick%pol(1:3, 3))
         enddo
      enddo

  end select

  ! And now, perform the necessary transformation.
  ip(1:3, 1:3) = kick%pol(1:3, 1:3)
  dump = lalg_inverter(3, ip)
  do is = 1, nspin
     do i = 0, energy_steps
        sigma(:, :, i, is) = matmul( transpose(ip), matmul(sigmap(:, :, i, is), ip) )
     enddo
  enddo

  ! Finally, write down the result
  write(out_file, '(a15,i2)')      '# nspin        ', nspin
  write(out_file, '(a15,3f18.12)') '# pol(1)       ', kick%pol(1:3, 1)
  write(out_file, '(a15,3f18.12)') '# pol(2)       ', kick%pol(1:3, 2)
  write(out_file, '(a15,3f18.12)') '# pol(3)       ', kick%pol(1:3, 3)
  write(out_file, '(a15,i1)')      '# direction    ', kick%pol_dir
  write(out_file, '(a15,i1)')      '# kick mode    ', kick%delta_strength_mode
  write(out_file, '(a15,f18.12)')  '# kick strength', kick%delta_strength
  write(out_file, '(a15,i1)')      '# Equiv. axis  ', kick%pol_equiv_axis
  write(out_file, '(a15,3f18.12)') '# wprime       ', kick%wprime(1:3)
  header_string = str_center("Energy", 20)
  write(out_file, '(a20)', advance = 'no') header_string
  header_string = str_center("(1/3)*Tr[sigma]", 20)
  write(out_file, '(a20)', advance = 'no') header_string
  header_string = str_center("Anisotropy[sigma]", 20)
  write(out_file, '(a20)', advance = 'no') header_string
  do j = 1, nspin
     do i = 1, 3
        do k = 1, 3
           write(header_string,'(a6,i1,a1,i1,a1,i1,a1)') 'sigma(',i,',',k,',',j,')'
           header_string = str_center(trim(header_string),20)
           write(out_file, '(a20)', advance = 'no') header_string         
        enddo
     enddo
  enddo
  write(out_file, *)
  header_string = str_center('['//trim(units_out%energy%abbrev) // ']',20)
  write(out_file, '(a20)', advance = 'no') header_string
  do i = 1, 2+nspin*9
     header_string = str_center('['//trim(units_out%length%abbrev) //'^2]',20)
     write(out_file, '(a20)', advance = 'no') header_string
  enddo
  write(out_file,*)

  do i = 0, energy_steps
      average = M_ZERO
      anisotropy = M_ZERO
      do j = 1, nspin
         average = average + M_THIRD* ( sigma(1, 1, i, 1) + sigma(2, 2, i, 1) + sigma(3, 3, i, 1) )
         sigmap(:, :, i, 1) = matmul(sigma(:, :, i, 1),sigma(:, :, i, 1))
         anisotropy = anisotropy + &
                      M_THIRD * ( M_THREE * (sigmap(1, 1, i, 1) + sigmap(2, 2, i, 1) + sigmap(3, 3, i, 1)) - &
                      (sigma(1, 1, i, 1) + sigma(2, 2, i, 1) + sigma(3, 3, i, 1))**2 )
      enddo
      write(out_file,'(3e20.8)', advance = 'no') (i*s%energy_step) / units_out%energy%factor, &
           average / (units_out%length%factor**2), sqrt(anisotropy) / (units_out%length%factor**2)
      do j = 1, nspin
         write(out_file,'(9e20.8)', advance = 'no') sigma(1:3, 1:3, i, j) / (units_out%length%factor**2)
      enddo
      write(out_file,'(a)', advance = 'yes')
  end do

  deallocate(sigma, sigmap, sigmau, sigmav, sigmaw, p, ip)
  call pop_sub()
end subroutine spectrum_cross_section_tensor


subroutine spectrum_cross_section(in_file, out_file, s)
  integer, intent(in)    :: in_file
  integer, intent(in)    :: out_file
  type(spec_type),  intent(inout) :: s

  character(len=20) :: header_string
  integer :: nspin, lmax, time_steps, is, ie, ntiter, i, j, jj, isp, no_e, k
  FLOAT   :: dt, dump, x, w, ewsum, polsum
  type(kick_type) :: kick
  FLOAT, allocatable :: dipole(:, :, :), sigma(:, :, :), dumpa(:), sf(:, :)

  call push_sub('spectrum.spectrum_cross_section')

  ! This function gives us back the unit connected to the "multipoles" file, the header information,
  ! the number of time steps, and the time step. 
  call spectrum_mult_info(in_file, nspin, lmax, kick, time_steps, dt)

  ! Now we cannot process files that do not contain the dipole, or that contain more than the dipole.
  if(lmax.ne.1) then
    message(1) = 'multipoles file should contain the dipole -- and only the dipole.'
    call write_fatal(1)
  endif

  ! Find out the iteration numbers corresponding to the time limits.
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

  ! Read the dipole.
  call spectrum_skip_header(in_file)
  allocate(dipole(3, 0:time_steps, nspin))
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

    dipole(1:3, i, :) = dipole(1:3, i, :) * units_out%length%factor
  end do
  ! Now substract the initial dipole.
  do i = time_steps, 0, -1
     dipole(:, i, :) = dipole(:, i, :) - dipole(:, 0, :)
  enddo

  ! Get the number of energy steps.
  no_e = s%max_energy / s%energy_step
  allocate(sigma(3, 0:no_e, nspin)); sigma = M_ZERO
  allocate(sf(0:no_e, nspin)); sf = M_ZERO

  ! Gets the damping function (here because otherwise it is awfully slow in "pol" mode...)
  allocate(dumpa(is:ie))
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
  enddo

  do k = 0, no_e
    w = k*s%energy_step
    do j = is, ie
      jj = j - is
      x = sin(w*jj*dt)
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
  enddo
  ewsum = ewsum * s%energy_step; polsum = polsum * s%energy_step

  write(message(1), '(a,i8)')    'Number of spin       = ', nspin
  write(message(2), '(a,i8)')    'Number of time steps = ', time_steps
  write(message(3), '(a,i4)')    'SpecDampMode         = ', s%damp
  write(message(4), '(a,f10.4)') 'SpecDampFactor       = ', s%damp_factor * units_out%time%factor
  write(message(5), '(a,f10.4)') 'SpecStartTime        = ', s%start_time   / units_out%time%factor
  write(message(6), '(a,f10.4)') 'SpecEndTime          = ', s%end_time     / units_out%time%factor
  write(message(7), '(a,f10.4)') 'SpecMaxEnergy        = ', s%max_energy   / units_inp%energy%factor
  write(message(8), '(a,f10.4)') 'SpecEnergyStep       = ', s%energy_step  / units_inp%energy%factor
  call write_info(8)

  message(1) = ""
  write(message(2),'(a,f16.6)')   'Electronic sum rule  = ', ewsum
    write(message(3),'(a,f16.6)') 'Polariz. (sum rule)  = ', polsum  / units_inp%length%factor**3
  call write_info(3)

  write(out_file, '(a15,i2)')      '# nspin        ', nspin
  write(out_file, '(a15,3f18.12)') '# pol(1)       ', kick%pol(1:3, 1)
  write(out_file, '(a15,3f18.12)') '# pol(2)       ', kick%pol(1:3, 2)
  write(out_file, '(a15,3f18.12)') '# pol(3)       ', kick%pol(1:3, 3)
  write(out_file, '(a15,i1)')      '# direction    ', kick%pol_dir
  write(out_file, '(a15,i1)')      '# kick mode    ', kick%delta_strength_mode
  write(out_file, '(a15,f18.12)')  '# kick strength', kick%delta_strength
  write(out_file, '(a15,i1)')      '# Equiv. axis  ', kick%pol_equiv_axis
  write(out_file, '(a15,3f18.12)') '# wprime       ', kick%wprime(1:3)
  header_string = str_center("Energy", 20)
  write(out_file, '(a20)', advance = 'no') header_string
  do j = 1, nspin
     do i = 1, 3
           write(header_string,'(a6,i1,a8,i1,a1)') 'sigma(',i,', nspin=',j,')'
           header_string = str_center(trim(header_string),20)
           write(out_file, '(a20)', advance = 'no') header_string         
     enddo
  enddo
  do j = 1, nspin
     write(header_string,'(a18,i1,a1)') 'StrengthFunction(',j,')'
     header_string = str_center(trim(header_string),20)
     write(out_file, '(a20)', advance = 'no') header_string
  enddo
  write(out_file, *)
  header_string = str_center('['//trim(units_out%energy%abbrev) // ']',20)
  write(out_file, '(a20)', advance = 'no') header_string
  do i = 1, nspin*3
     header_string = str_center('['//trim(units_out%length%abbrev) //'^2]',20)
     write(out_file, '(a20)', advance = 'no') header_string
  enddo
  do i = 1, nspin
     header_string = str_center('[1/'//trim(units_out%energy%abbrev) //']',20)
     write(out_file, '(a20)', advance = 'no') header_string
  enddo
  write(out_file,*)

  do i = 0, no_e
     write(out_file,'(e20.8)', advance = 'no') i*s%energy_step / units_out%energy%factor
     do j = 1, nspin
        write(out_file,'(3e20.8)', advance = 'no') sigma(1:3, i, j) / (units_out%length%factor**2) 
     enddo
     do j = 1, nspin
        write(out_file,'(e20.8)', advance = 'no') sf(i, j) * units_out%energy%factor
     enddo
     write(out_file, *) 
  end do

  deallocate(dipole, sigma)
  call pop_sub()
end subroutine spectrum_cross_section

subroutine spectrum_rotatory_strength(in_file, out_file, s)
  integer, intent(in) :: in_file
  integer, intent(in) :: out_file
  type(spec_type), intent(inout) :: s

  integer :: i, is, ie, ntiter, j, jj, k, time_steps, no_e, nspin
  FLOAT :: dump, dt
  type(kick_type) :: kick
  CMPLX :: z
  CMPLX, pointer :: sp(:)
  FLOAT, allocatable :: dumpa(:)
  FLOAT, allocatable :: angular(:, :)

  call push_sub('spectrum_rotatory_strength')

  call spectrum_angular_info(in_file, nspin, kick, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

  ! load dipole from file
  allocate(angular(0:time_steps, 3))
  call spectrum_skip_header(in_file)
  do i = 0, time_steps
    read(in_file, *) j, dump, angular(i, 1:3)
    !dipole(i,:) = dipole(i,:) * units_out%length%factor
  end do

  ! subtract static dipole
  do i = 1, 3
     angular(:, i) = angular(:, i) - angular(0, i)
  end do

  no_e = s%max_energy / s%energy_step
  allocate(sp(0:no_e))
  sp = M_z0

  ! Gets the damping function (here because otherwise it is awfully slow in "pol" mode...)
  allocate(dumpa(is:ie))
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
  enddo

  do k = 0, no_e
    do j = is, ie

      jj = j - is

      z = exp(M_zI * k * s%energy_step * jj *dt)
      sp(k) = sp(k) + z*dumpa(j)*sum(angular(j, :)*kick%pol(1:3, kick%pol_dir))

    end do
    sp(k) = sp(k)*dt

  end do

  deallocate(angular, dumpa)

  ! should output units, etc...
  do i = 0, no_e
      write(out_file,'(5e15.6)') i*s%energy_step / units_out%energy%factor, &
           sp(i) * (units_out%length%factor)**3
  end do

  ! print some info
  write(message(1), '(a,i8)')    'Number of time steps = ', ntiter
  write(message(2), '(a,i4)')    'SpecDampMode         = ', s%damp
  write(message(3), '(a,f10.4)') 'SpecDampFactor       = ', s%damp_factor * units_out%time%factor
  write(message(4), '(a,f10.4)') 'SpecStartTime        = ', s%start_time   / units_out%time%factor
  write(message(5), '(a,f10.4)') 'SpecEndTime          = ', s%end_time     / units_out%time%factor
  write(message(6), '(a,f10.4)') 'SpecMaxEnergy        = ', s%max_energy   / units_inp%energy%factor
  write(message(7),'(a,f10.4)')  'SpecEnergyStep       = ', s%energy_step  / units_inp%energy%factor
  call write_info(7)

  call pop_sub()
end subroutine spectrum_rotatory_strength

subroutine spectrum_hs_from_mult(out_file, s, sh)
  character(len=*), intent(in) :: out_file
  type(spec_type), intent(inout) :: s
  type(spec_sh), intent(inout) :: sh

  integer :: i, j, iunit, nspin, time_steps, is, ie, ntiter, lmax
  FLOAT :: dt, dump
  type(kick_type) :: kick
  FLOAT, allocatable :: d(:,:)
  CMPLX :: c
  CMPLX, allocatable :: dipole(:), ddipole(:)

  call io_assign(iunit)
  iunit = io_open('multipoles', action='read', status='old', die=.false.)
  if(iunit < 0) then
    iunit = io_open('td.general/multipoles', action='read', status='old')
  end if
  call spectrum_mult_info(iunit, nspin, lmax, kick, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

  call spectrum_skip_header(iunit)
  ! load dipole from file
  allocate(dipole(0:time_steps))
  allocate(d(3, nspin))
  do i = 1, time_steps
    read(iunit, *) j, dump, dump, dump, d
    select case(sh%pol)
    case('x')
      dipole(i) = -sum(d(3, :))
    case('y')
      dipole(i) = -sum(d(1, :))
    case('z')
      dipole(i) =  sum(d(2, :))
    case('+')
      dipole(i) = -sum(d(3, :) + M_zI*d(1, :)) / sqrt(M_TWO)
    case('-')
      dipole(i) = -sum(d(3, :) - M_zI*d(1, :)) / sqrt(M_TWO)
    end select
    dipole(i) = dipole(i) * units_out%length%factor * sqrt(M_FOUR*M_PI/M_THREE)
  end do
  deallocate(d)

  ! we now calculate the first time derivative
  allocate(ddipole(0:time_steps))
  ddipole(0) = (dipole(1) - dipole(0))/dt
  do i = 1, time_steps - 1
    ddipole(i) = (dipole(i + 1) - dipole(i - 1))/(M_TWO*dt)
  end do
  ddipole(time_steps) = (dipole(time_steps) - dipole(time_steps - 1))/dt

  ! and the second time derivative
  dipole(0) = (ddipole(1) - ddipole(0))/dt
  do i = 1, time_steps - 1
    dipole(i) = (ddipole(i + 1) - ddipole(i - 1))/(M_TWO*dt)
  end do
  dipole(time_steps) = (ddipole(time_steps) - ddipole(time_steps - 1))/dt
  deallocate(ddipole)

  ! now we Fourier transform
  sh%no_e = s%max_energy / s%energy_step
  allocate(sh%sp(0:sh%no_e))
  sh%sp = M_ZERO

  do i = 0, sh%no_e
    c = M_z0
    do j = is, ie
      c = c + exp(M_zI * i * s%energy_step * j * dt)*dipole(j)
    end do
    sh%sp(i) = abs(c)**2
  end do
  deallocate(dipole)

  ! output
  if(trim(out_file) .ne. '-') then
    iunit = io_open(trim(out_file) // "." // trim(sh%pol), action='write')

    ! should output units, etc...
    do i = 0, sh%no_e
      write(iunit,'(5e15.6)') i*s%energy_step / units_out%energy%factor, &
           sh%sp(i) * units_out%energy%factor
    end do
    call io_close(iunit)
  end if

end subroutine spectrum_hs_from_mult

subroutine spectrum_hs_from_acc(out_file, s, sh)
  character(len=*), intent(in) :: out_file
  type(spec_type), intent(inout) :: s
  type(spec_sh), intent(inout) :: sh

  integer :: i, j, iunit, time_steps, is, ie, ntiter
  FLOAT :: dt, dummy, a(3)
  CMPLX, allocatable :: acc(:)
  CMPLX :: c

  call spectrum_acc_info(iunit, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

  ! load dipole from file
  allocate(acc(0:time_steps))
  acc = M_ZERO
  call spectrum_skip_header(iunit)
  do i = 1, time_steps
    read(iunit, *) j, dummy, a
    select case(sh%pol)
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
    acc(i) = acc(i) * units_out%acceleration%factor
  end do

  ! now we Fourier transform
  sh%no_e = s%max_energy / s%energy_step
  allocate(sh%sp(0:sh%no_e))
  sh%sp = M_ZERO

  do i = 0, sh%no_e
    c = M_z0
    do j = is, ie
      c = c + exp(M_zI * i * s%energy_step * j * dt)*acc(j)
    end do
    sh%sp(i) = abs(c)**2
  end do
  deallocate(acc)

  ! output
  if(trim(out_file) .ne. '-') then
    iunit = io_open(trim(out_file) // "." // trim(sh%pol), action='write')

    ! should output units, etc...
    do i = 0, sh%no_e
      write(iunit,'(5e15.6)') i*s%energy_step / units_out%energy%factor, &
           sh%sp(i) * units_out%energy%factor
    end do
    call io_close(iunit)
  end if

end subroutine spectrum_hs_from_acc

subroutine spectrum_angular_info(iunit, nspin, kick, time_steps, dt)
  integer, intent(in)          :: iunit
  integer, intent(out)         :: nspin, time_steps
  type(kick_type), intent(out) :: kick
  FLOAT, intent(out)           :: dt

  integer :: j
  FLOAT :: t1, t2, dummy

  call push_sub('spectrum.spectrum_angular_info')

  rewind(iunit); read(iunit,*); read(iunit,*)
  read(iunit, '(15x,i2)')      nspin
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 1)
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 2)
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 3)
  read(iunit, '(15x,i2)')      kick%pol_dir
  read(iunit, '(15x,i2)')      kick%delta_strength_mode
  read(iunit, '(15x,f18.12)')  kick%delta_strength
  read(iunit, '(15x,i2)')      kick%pol_equiv_axis
  read(iunit, '(15x,3f18.12)') kick%wprime(1:3)
  call spectrum_skip_header(iunit)

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
100 continue
  dt = (t2 - t1) * units_out%time%factor ! units_out is OK
  time_steps = time_steps - 1

  if(time_steps < 3) then
    write(message(1),'(a)') "Empty file?"
    call write_fatal(1)
  end if

  call pop_sub()
end subroutine spectrum_angular_info

subroutine spectrum_mult_info(iunit, nspin, lmax, kick, time_steps, dt)
  integer, intent(in)  :: iunit
  integer, intent(out) :: nspin
  integer, intent(out) :: lmax
  type(kick_type), intent(out) :: kick
  integer, intent(out) :: time_steps
  FLOAT,   intent(out) :: dt

  integer :: j
  FLOAT :: t1, t2, dummy

  call push_sub('spectrum.spectrum_mult_info')

  rewind(iunit); read(iunit,*); read(iunit,*)
  read(iunit, '(15x,i2)')      nspin
  read(iunit, '(15x,i2)')      lmax
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 1)
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 2)
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 3)
  read(iunit, '(15x,i2)')      kick%pol_dir
  read(iunit, '(15x,i2)')      kick%delta_strength_mode
  read(iunit, '(15x,f18.12)')  kick%delta_strength
  read(iunit, '(15x,i2)')      kick%pol_equiv_axis
  read(iunit, '(15x,3f18.12)') kick%wprime(1:3)
  call spectrum_skip_header(iunit)

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
100 continue
  dt = (t2 - t1) * units_out%time%factor ! units_out is OK
  time_steps = time_steps - 1

  if(time_steps < 3) then
    message(1) = "Empty multipole file?"
    call write_fatal(1)
  end if

  call pop_sub()
end subroutine spectrum_mult_info

subroutine spectrum_cross_section_info(iunit, nspin, kick, energy_steps, dw)
  integer, intent(in)  :: iunit
  integer, intent(out) :: nspin
  type(kick_type), intent(out) :: kick
  integer, intent(out) :: energy_steps
  FLOAT,   intent(out) :: dw

  integer :: j
  FLOAT :: dummy, e1, e2

  call push_sub('spectrum.spectrum_cross_section_info')

  ! read in number of spin components
  read(iunit, '(15x,i2)')      nspin
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 1)
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 2)
  read(iunit, '(15x,3f18.12)') kick%pol(1:3, 3)
  read(iunit, '(15x,i2)')      kick%pol_dir
  read(iunit, '(15x,i2)')      kick%delta_strength_mode
  read(iunit, '(15x,f18.12)')  kick%delta_strength
  read(iunit, '(15x,i2)')      kick%pol_equiv_axis
  read(iunit, '(15x,3f18.12)') kick%wprime(1:3)
  read(iunit, *); read(iunit, *) ! skip header

  ! count number of time_steps
  energy_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    energy_steps = energy_steps + 1
    if(energy_steps == 1) e1 = dummy
    if(energy_steps == 2) e2 = dummy
  end do
100 continue
  dw = (e2 - e1) * units_out%energy%factor
  energy_steps = energy_steps - 1

  if(energy_steps < 3) then
    message(1) = "Empty multipole file?"
    call write_fatal(1)
  end if

  call pop_sub()
end subroutine spectrum_cross_section_info

subroutine spectrum_acc_info(iunit, time_steps, dt)
  integer, intent(out) :: iunit, time_steps
  FLOAT, intent(out) :: dt

  integer :: j
  FLOAT :: t1, t2, dummy

  ! open files
  iunit = io_open('acceleration', action='read', status='old', die=.false.)
  if(iunit < 0) then
    iunit = io_open('td.general/acceleration', action='read', status='old')
  endif

  ! read in dipole
  call spectrum_skip_header(iunit)

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
100 continue
  dt = (t2 - t1) * units_out%time%factor ! units_out is OK
  time_steps = time_steps - 1

  if(time_steps < 3) then
    message(1) = "Empty multipole file?"
    call write_fatal(1)
  end if

  rewind(iunit)
end subroutine spectrum_acc_info

subroutine spectrum_fix_time_limits(time_steps, dt, start_time, end_time, is, ie, ntiter)
  integer, intent(in) :: time_steps
  FLOAT, intent(in) :: dt
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

end module spectrum
