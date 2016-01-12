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
!! $Id$

#include "global.h"

module pes_spm_m
  use global_m
  use io_m
  use mesh_m
  use messages_m
  use mpi_m
  use comm_m
  use parser_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use unit_m
  use unit_system_m
  use hamiltonian_m
  use lasers_m
  use varinfo_m
  use mesh_interpolation_m

  implicit none

  private

  public ::                             &
    pes_spm_t,                           &
    pes_spm_init,                        &
    pes_spm_init_write,                  &
    pes_spm_output,                      &
    pes_spm_calc,                        &
    pes_spm_end,                         &
    pes_spm_dump,                        &
    pes_spm_load

  type pes_spm_t
    integer                    :: nspoints                  !< how many points we store the wf
    FLOAT, pointer             :: rcoords(:,:)              !< coordinates of the sample points
    CMPLX, pointer             :: wf(:,:,:,:,:)   => NULL() !< wavefunctions at sample points
    FLOAT, pointer             :: dq(:,:)         => NULL() !< part 1 of Volkov phase (recipe phase) 
    FLOAT, pointer             :: domega(:)       => NULL() !< part 2 of Volkov phase (recipe phase)
    integer                    :: recipe                    !< type of calculation (RAW/PHASE)
    CMPLX, pointer             :: wfft(:,:,:,:,:) => NULL() !< Fourier transform of wavefunction
    FLOAT                      :: omegamax                  !< maximum frequency of the spectrum
    FLOAT                      :: delomega                  !< frequency spacing of the spectrum
    integer                    :: nomega                    !< number of frequencies of the spectrum
    logical                    :: onfly                     !< spectrum is calculated on-the-fly when true
    integer                    :: save_iter                 !< output interval and size of arrays 
    logical                    :: sphgrid                   !< use a spherical grid (instead of sample points from input)
    integer                    :: nstepsphi, nstepstheta
    FLOAT                      :: thetamin
    type(mesh_interpolation_t) :: interp
  end type pes_spm_t

  integer, parameter :: &
    M_RAW   = 1,        &
    M_PHASE = 2

contains

  ! ---------------------------------------------------------
  subroutine pes_spm_init(this, mesh, st, save_iter)
    type(pes_spm_t), intent(out) :: this
    type(mesh_t),    intent(in)  :: mesh
    type(states_t),  intent(in)  :: st
    integer,         intent(in)  :: save_iter

    type(block_t) :: blk
    integer       :: sdim, mdim
    integer       :: isp
    integer       :: ith, iph
    FLOAT         :: theta, phi, radius
    FLOAT         :: xx(MAX_DIM)

    PUSH_SUB(pes_spm_init)

    sdim = st%d%dim
    mdim = mesh%sb%dim

    message(1) = 'Info: Calculating PES using sample point technique.'
    call messages_info(1)

    !%Variable PES_spm_points
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% List of points at which to calculate the photoelectron spectrum by the sample point 
    !% method. If no points are given, a spherical grid is generated automatically.
    !% The exact syntax is:
    !%
    !% <tt>%PES_spm_points
    !% <br>&nbsp;&nbsp;x1 | y1 | z1
    !% <br>%
    !% </tt>
    !%End
    call messages_obsolete_variable('PhotoElectronSpectrumPoints', 'PES_spm_points')
    this%sphgrid = .false.
    if (parse_block('PES_spm_points', blk) < 0) then
      this%sphgrid = .true.
    end if

    if(this%sphgrid) then
      message(1) = 'Info: Using spherical grid.'
    else
      message(1) = 'Info: Using sample points from block.'
    end if
    call messages_info(1)

    !%Variable PES_spm_recipe
    !%Type integer
    !%Default phase
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The type for calculating the photoelectron spectrum in the sample point method.
    !%Option raw 1
    !% Calculate the photoelectron spectrum according to A. Pohl, P.-G. Reinhard, and 
    !% E. Suraud, <i>Phys. Rev. Lett.</i> <b>84</b>, 5090 (2000).
    !%Option phase 2
    !% Calculate the photoelectron spectrum by including the Volkov phase (approximately), see
    !% P. M. Dinh, P. Romaniello, P.-G. Reinhard, and E. Suraud, <i>Phys. Rev. A.</i> <b>87</b>, 032514 (2013).
    !%End
    call parse_variable('PES_spm_recipe', M_RAW, this%recipe)
    if(.not.varinfo_valid_option('PES_spm_recipe', this%recipe, is_flag = .true.)) &
      call messages_input_error('PES_spm_recipe')
    call messages_print_var_option(stdout, "PES_spm_recipe", this%recipe)

    !%Variable PES_spm_OmegaMax
    !%Type float
    !%Default 0.0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% If PES_spm_OmegaMax > 0, the photoelectron spectrum is directly calculated during 
    !% time-propagation, evaluated by the PES_spm method. PES_spm_OmegaMax is then the maximum frequency 
    !% (approximate kinetic energy) and PES_spm_DeltaOmega the spacing in frequency domain of the spectrum.
    !%End
    call parse_variable('PES_spm_OmegaMax', units_to_atomic(units_inp%energy, M_ZERO), this%omegamax)
    this%onfly = .false.
    if(this%omegamax > M_ZERO) then
      this%onfly = .true.
      message(1) = 'Info: Calculating PES during time propagation.'
      call messages_info(1)
      call messages_print_var_value(stdout, "PES_spm_OmegaMax", this%omegamax)
    end if
 
    !%Variable PES_spm_DeltaOmega
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The spacing in frequency domain for the photoelectron spectrum (if PES_spm_OmegaMax > 0).
    !% By default is set to PES_spm_OmegaMax/500. 
    !%End
    call parse_variable('PES_spm_DeltaOmega', units_to_atomic(units_inp%energy, this%omegamax/CNST(500)), this%delomega)
    if(this%onfly) then
      if(this%delomega <= M_ZERO) call messages_input_error('PES_spm_DeltaOmega')
      call messages_print_var_value(stdout, "PES_spm_DeltaOmega", this%delomega)
    end if

    !%Variable PES_spm_StepsThetaR
    !%Type integer
    !%Default 45
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in theta (0 <= theta <= pi) for the spherical grid (if no 
    !% PES_spm_points are given).
    !%End
    call parse_variable('PES_spm_StepsThetaR', 45, this%nstepstheta)
    if(this%sphgrid .and. this%nstepstheta < 0) call messages_input_error('PES_spm_StepsThetaR')

    !%Variable PES_spm_StepsPhiR
    !%Type integer
    !%Default 90
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in phi (0 <= phi <= 2 pi) for the spherical grid (if no
    !% PES_spm_points are given).
    !%End
    call parse_variable('PES_spm_StepsPhiR', 90, this%nstepsphi)
    if(this%sphgrid) then
      if(this%nstepsphi < 0)  call messages_input_error('PES_spm_StepsPhiR')
      if(this%nstepsphi == 0) this%nstepsphi = 1
    end if

    !%Variable PES_spm_Radius
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The radius of the sphere for the spherical grid (if no PES_spm_points
    !% are given).
    !%End
    if(this%sphgrid) then
      if(parse_is_defined('PES_spm_Radius')) then
        call parse_variable('PES_spm_Radius', M_ZERO, radius)
        if(radius <= M_ZERO) call messages_input_error('PES_spm_Radius')
        call messages_print_var_value(stdout, "PES_spm_Radius", radius)
      else
        select case(mesh%sb%box_shape)
        case(PARALLELEPIPED)
          radius = minval(mesh%sb%lsize(1:mdim))
        case(SPHERE)
          radius = mesh%sb%rsize
        case default
          message(1) = "Spherical grid not implemented for this box shape."
          message(2) = "Specify sample points with block PES_spm_points."
          call messages_fatal(2)
        end select
      end if
    end if

    call mesh_interpolation_init(this%interp, mesh)

    if(.not. this%sphgrid) then
      this%nspoints = parse_block_n(blk)
    else
      ! setting values for spherical grid 
      select case(mdim)
      case(1)
        this%thetamin = M_PI / M_TWO
        this%nstepstheta   = 0
        this%nstepsphi     = 2
        this%nspoints  = this%nstepsphi

      case(2)
        this%thetamin = M_PI / M_TWO
        this%nstepstheta   = 0
        this%nspoints  = this%nstepsphi

        call messages_print_var_value(stdout, "PES_spm_StepsPhiR", this%nstepsphi)

      case(3)
        this%thetamin = M_ZERO
        if(this%nstepstheta <= 1) this%nstepsphi = 1
        this%nspoints  = this%nstepsphi * (this%nstepstheta - 1) + 2

        call messages_print_var_value(stdout, "PES_spm_StepsPhiR", this%nstepsphi)
        call messages_print_var_value(stdout, "PES_spm_StepsThetaR", this%nstepstheta)

      end select
    end if

    call messages_print_var_value(stdout, "Number of sample points", this%nspoints)

    SAFE_ALLOCATE(this%rcoords(1:mdim, 1:this%nspoints))

    if(.not. this%sphgrid) then

      message(1) = 'Info: Reading sample points.'
      call messages_info(1)

      ! read points from input file
      do isp = 1, this%nspoints
        call parse_block_float(blk, isp - 1, 0, xx(1), units_inp%length)
        call parse_block_float(blk, isp - 1, 1, xx(2), units_inp%length)
        call parse_block_float(blk, isp - 1, 2, xx(3), units_inp%length)
        this%rcoords(1:mdim, isp) = xx(1:mdim)
      end do
      call parse_block_end(blk)

    else ! this%sphgrid == .true.

      message(1) = 'Info: Initializing spherical grid.'
      call messages_info(1)

      ! initializing spherical grid
      isp = 0
      do ith = 0, this%nstepstheta
        if(ith == 0) then
          theta = this%thetamin     ! pi/2 for 1d and 2d. zero for 3d.
        else
          theta = ith * M_PI / this%nstepstheta
        end if
        do iph = 0, this%nstepsphi - 1
          if(iph == 0) then
            phi = M_ZERO
          else
            phi = iph * M_TWO * M_PI / this%nstepsphi
          end if
          isp = isp + 1
                        this%rcoords(1, isp) = radius * cos(phi) * sin(theta)
          if(mdim >= 2) this%rcoords(2, isp) = radius * sin(phi) * sin(theta)
          if(mdim == 3) this%rcoords(3, isp) = radius * cos(theta)
          if(theta == M_ZERO .or. theta == M_PI) exit
        end do
      end do

    end if

    SAFE_ALLOCATE(this%wf(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nspoints, 0:save_iter-1))

    if(this%recipe == M_PHASE) then
      SAFE_ALLOCATE(this%dq(1:this%nspoints, 0:save_iter-1))
      SAFE_ALLOCATE(this%domega(0:save_iter-1))
      this%dq = M_ZERO
      this%domega = M_ZERO
    end if

    if(this%onfly) then
      this%nomega = nint(this%omegamax/this%delomega)
      SAFE_ALLOCATE(this%wfft(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nspoints, 1:this%nomega))
      this%wfft = M_z0
    end if

    this%save_iter = save_iter

    POP_SUB(pes_spm_init)
  end subroutine pes_spm_init


  ! ---------------------------------------------------------
  subroutine pes_spm_end(this)
    type(pes_spm_t), intent(inout) :: this

    PUSH_SUB(pes_spm_end)

    SAFE_DEALLOCATE_P(this%wf)
    SAFE_DEALLOCATE_P(this%rcoords)

    SAFE_DEALLOCATE_P(this%wfft)

    SAFE_DEALLOCATE_P(this%dq)
    SAFE_DEALLOCATE_P(this%domega)

    POP_SUB(pes_spm_end)
  end subroutine pes_spm_end


  ! ---------------------------------------------------------
  subroutine pes_spm_calc(this, st, mesh, dt, iter, hm)
    type(pes_spm_t),     intent(inout) :: this
    type(states_t),      intent(in)    :: st
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter
    type(hamiltonian_t), intent(in)    :: hm

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim
    integer            :: ii

    CMPLX, allocatable :: psistate(:), wfftact(:,:,:,:,:)
    CMPLX              :: rawfac
    CMPLX, allocatable :: phasefac(:)
    integer            :: iom
    FLOAT              :: omega
    CMPLX, allocatable :: interp_values(:)

    PUSH_SUB(pes_spm_calc)

    ii = mod(iter-1, this%save_iter)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = mesh%sb%dim

    SAFE_ALLOCATE(psistate(1:mesh%np_part))
    SAFE_ALLOCATE(interp_values(1:this%nspoints))

    if(this%onfly) then
      SAFE_ALLOCATE(wfftact(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nspoints, 1:this%nomega))
      wfftact = M_z0
    endif

    if(this%recipe == M_PHASE) then
      SAFE_ALLOCATE(phasefac(1:this%nspoints))
    end if

    ! needed for allreduce, otherwise it will take values from previous cycle
    this%wf(:,:,:,:,ii) = M_z0

    do ik = kptst, kptend 
      do ist = stst, stend
        do isdim = 1, sdim
          call states_get_state(st, mesh, isdim, ist, ik, psistate(1:mesh%np_part))
          call mesh_interpolation_evaluate(this%interp, this%nspoints, psistate(1:mesh%np_part), &
            this%rcoords(1:mdim, 1:this%nspoints), interp_values(1:this%nspoints))
          this%wf(ist, isdim, ik, :, ii) = interp_values(:)
        end do
      end do
    end do
    if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
      ! interpolated values have already been communicated over domains
      call comm_allreduce(st%st_kpt_mpi_grp%comm, this%wf(:,:,:,:,ii))
#endif
    end if

    if(this%recipe == M_PHASE) then
      call pes_spm_calc_rcphase(this, mesh, iter, dt, hm, ii)
    end if

    if(this%onfly) then
      do iom = 1, this%nomega
        omega = iom*this%delomega
        rawfac = exp(M_zI * omega * iter * dt) * dt * sqrt(M_TWO * omega) / (M_TWO * M_PI)**(mdim/M_TWO)

        select case(this%recipe)
        case(M_RAW)
          wfftact(stst:stend, 1:sdim, kptst:kptend, :, iom) = this%wf(stst:stend, 1:sdim, kptst:kptend, :, ii) * rawfac
        case(M_PHASE)
          phasefac(:) = exp(M_zI * (this%domega(ii) - sqrt(M_TWO * omega) * this%dq(:, ii)))

          do ik = kptst, kptend
            do ist = stst, stend
              do isdim = 1, sdim
                wfftact(ist, isdim, ik, :, iom) = this%wf(ist, isdim, ik, :, ii) * phasefac(:) * rawfac
              end do
            end do
          end do
        end select
        if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
          call comm_allreduce(st%st_kpt_mpi_grp%comm, wfftact(:,:,:,:,iom))
#endif
        end if
      end do

      this%wfft = this%wfft + wfftact
    end if

    SAFE_DEALLOCATE_A(psistate)
    SAFE_DEALLOCATE_A(interp_values)
    SAFE_DEALLOCATE_A(phasefac)

    SAFE_DEALLOCATE_A(wfftact)

    POP_SUB(pes_spm_calc)
  end subroutine pes_spm_calc


  ! ---------------------------------------------------------
  subroutine pes_spm_output(this, mesh, st, iter, dt)
    type(pes_spm_t), intent(in) :: this
    type(mesh_t),    intent(in) :: mesh
    type(states_t),  intent(in) :: st
    integer,         intent(in) :: iter
    FLOAT,           intent(in) :: dt

    integer            :: ist, ik, isdim
    integer            :: ii, jj
    integer            :: isp, save_iter, isp_save
    integer            :: iom, ith, iph, iphi
    FLOAT              :: omega, theta, phi
    CMPLX              :: vfu
    FLOAT              :: wfu
    FLOAT, allocatable :: wffttot(:,:)
    FLOAT              :: spctrsum, weight
    character(len=4)   :: filenr
    integer            :: iunitone, iunittwo
    integer            :: mdim

    PUSH_SUB(pes_spm_output)

    mdim = mesh%sb%dim

    save_iter = this%save_iter

    if(this%onfly) then
      SAFE_ALLOCATE(wffttot(1:this%nomega, 1:this%nspoints))
      wffttot = M_ZERO

      ! calculate total spectrum
      do iom = 1, this%nomega
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            do isdim = 1, st%d%dim
              wffttot(iom, :) = abs(this%wfft(ist, isdim, ik, :, iom))**M_TWO + wffttot(iom, :)
            end do
          end do
        end do
      end do
    end if

    if(.not. this%sphgrid .or. debug%info) then   ! too much output for spherical grid
      do isp = 1, this%nspoints
        write(filenr, '(i4.4)') isp
   
        iunitone = io_open('td.general/'//'PES_spm.'//filenr//'.wavefunctions.out', action='write', position='append')
        if(this%recipe == M_PHASE) &
          iunittwo = io_open('td.general/'//'PES_spm.'//filenr//'.phase.out', action='write', position='append')
   
        do ii = 1, save_iter - mod(iter, save_iter)
          jj = iter - save_iter + ii + mod(save_iter - mod(iter, save_iter), save_iter)
          write(iunitone, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do isdim = 1, st%d%dim
                vfu = units_from_atomic(sqrt(units_out%length**(-3)), this%wf(ist, isdim, ik, isp, ii-1))
                write(iunitone, '(1x,e18.10E3,1x,e18.10E3)', advance='no') &
                  real(vfu),  aimag(vfu)
              end do
            end do
          end do
          write(iunitone, '(1x)', advance='yes')

          if(this%recipe == M_PHASE) then
            write(iunittwo, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
            write(iunittwo, '(1x,e18.10E3,1x,e18.10E3)', advance='no') &
              this%domega(ii-1), this%dq(isp, ii-1)
            write(iunittwo, '(1x)', advance='yes')
          end if
        end do
        call io_close(iunitone)
        if(this%recipe == M_PHASE) call io_close(iunittwo)

        if(this%onfly) then
          iunitone = io_open('td.general/'//'PES_spm.'//filenr//'.spectrum.out', action='write', position='rewind')
          write(iunitone, '(a44)') '# frequency, total spectrum, orbital spectra'
          do iom = 1, this%nomega 
            omega = iom*this%delomega
            write(iunitone, '(e17.10, 1x, e17.10)', advance='no') omega, wffttot(iom, isp)
            do ik = 1, st%d%nik
              do ist = 1, st%nst 
                do isdim = 1, st%d%dim
                  wfu = abs(this%wfft(ist, isdim, ik, isp, iom))**M_TWO
                  write(iunitone,'(1x,e18.10e3)', advance='no') wfu
                end do
              end do
            end do
          write(iunitone,'(1x)', advance='yes')
          end do
          call io_close(iunitone)
        end if
   
      end do
    end if

    if(this%onfly .and. this%sphgrid) then
      iunittwo = io_open('td.general/'//'PES_spm.distribution.out', action='write', position='rewind')
      iunitone = io_open('td.general/'//'PES_spm.power.sum', action='write', position='rewind')
      write(iunitone, '(a23)') '# omega, total spectrum'

      select case(mdim)
      case(1)
        write(iunittwo, '(a40)') '# omega, distribution (left/right point)'

        do iom = 1, this%nomega
          omega = iom * this%delomega
          write(iunittwo, '(5(1x,e18.10E3))') omega, wffttot(iom, 2), wffttot(iom, 1)
          write(iunitone, '(2(1x,e18.10E3))') omega, sum(wffttot(iom, :)) / M_TWO * sqrt(M_TWO * omega)
        end do

      case(2)
        write(iunittwo, '(a26)') '# omega, phi, distribution'

        do iom = 1, this%nomega
          omega = iom * this%delomega

          spctrsum = M_ZERO
          do iph = 0, this%nstepsphi - 1
            spctrsum = spctrsum + wffttot(iom, iph + 1) * this%nstepsphi / M_TWO / M_PI
            phi = iph * M_TWO * M_PI / this%nstepsphi
            write(iunittwo,'(5(1x,e18.10E3))') omega, phi, wffttot(iom, iph + 1)
          end do
          ! just repeat the result for output
          write(iunittwo,'(5(1x,e18.10E3))') omega, M_TWO * M_PI, wffttot(iom, 1)
          write(iunittwo, '(1x)', advance='yes')
          write(iunitone, '(2(1x,e18.10E3))') omega, spctrsum
        end do

      case(3)
        write(iunittwo, '(a33)') '# omega, theta, phi, distribution'

        do iom = 1, this%nomega
          omega = iom * this%delomega
          isp = 0
          spctrsum = M_ZERO

          do ith = 0, this%nstepstheta
            theta = ith * M_PI / this%nstepstheta

            if(ith == 0 .or. ith == this%nstepstheta) then
              weight = (M_ONE - cos(M_PI / this%nstepstheta / M_TWO)) * M_TWO * M_PI
            else
              weight = abs(cos(theta - M_PI / this%nstepstheta / M_TWO) - cos(theta + M_PI / this%nstepstheta / M_TWO)) &
                * M_TWO * M_PI / this%nstepsphi
            end if

            do iph = 0, this%nstepsphi - 1
              isp = isp + 1
              spctrsum = spctrsum + wffttot(iom, isp) * weight

              phi = iph * M_TWO * M_PI / this%nstepsphi
              if(iph == 0) isp_save = isp
              write(iunittwo,'(5(1x,e18.10E3))') omega, theta, phi, wffttot(iom, isp)
   
              ! just repeat the result for output
              if(iph == (this%nstepsphi - 1)) then
                write(iunittwo,'(5(1x,e18.10E3))') omega, theta, M_TWO * M_PI, wffttot(iom, isp_save)
              end if
               
              ! just repeat the result for output
              if(theta == M_ZERO .or. theta == M_PI) then
                do iphi = 1, this%nstepsphi
                  phi = iphi * M_TWO * M_PI / this%nstepsphi
                  write(iunittwo,'(5(1x,e18.10E3))') omega, theta, phi, wffttot(iom, isp)
                end do
                exit
              end if
            end do

            write(iunittwo, '(1x)', advance='yes')
          end do
          write(iunitone, '(2(1x,e18.10E3))') omega, spctrsum * this%delomega
        end do
      end select
      call io_close(iunittwo)
      call io_close(iunitone)

    end if

    SAFE_DEALLOCATE_A(wffttot)

    POP_SUB(pes_spm_output)
  end subroutine pes_spm_output

  ! ---------------------------------------------------------
  subroutine pes_spm_init_write(this, mesh, st)
    type(PES_spm_t), intent(in) :: this
    type(mesh_t),    intent(in) :: mesh
    type(states_t),  intent(in) :: st

    integer          :: ist, ik, isdim
    integer          :: isp
    FLOAT            :: xx(MAX_DIM)
    character(len=4) :: filenr
    integer          :: iunit

    PUSH_SUB(pes_spm_init_write)

    xx = M_ZERO
    if(mpi_grp_is_root(mpi_world)) then
      if(.not. this%sphgrid .or. debug%info) then   ! too much output for spherical grid
        do isp = 1, this%nspoints
          write(filenr, '(i4.4)') isp
   
          iunit = io_open('td.general/'//'PES_spm.'//filenr//'.wavefunctions.out', action='write')
          xx(1:mesh%sb%dim) = this%rcoords(1:mesh%sb%dim, isp)
          write(iunit,'(a1)') '#'
          write(iunit, '(a7,f17.6,a1,f17.6,a1,f17.6,5a)') &
            '# R = (',xx(1),' ,',xx(2),' ,',xx(3), &
            ' )  [', trim(units_abbrev(units_inp%length)), ']'
   
          write(iunit,'(a1)') '#'  
          write(iunit, '(a3,14x)', advance='no') '# t' 
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do isdim = 1, st%d%dim
                write(iunit, '(3x,a8,i3,a7,i3,a8,i3,3x)', advance='no') &
                  "ik = ", ik, " ist = ", ist, " idim = ", isdim
              end do
            end do
          end do
          write(iunit, '(1x)', advance='yes')
   
          call io_close(iunit)
   
          if(this%recipe == M_PHASE) then
            iunit = io_open('td.general/'//'PES_spm.'//filenr//'.phase.out', action='write')
            write(iunit,'(a24)') '# time, dq(t), dOmega(t)'
            call io_close(iunit)
          end if
        end do
      end if
    end if

    POP_SUB(pes_spm_init_write)
  end subroutine pes_spm_init_write

  ! ---------------------------------------------------------
  subroutine pes_spm_dump(restart, this, st, ierr)
    type(restart_t), intent(in)  :: restart    
    type(pes_spm_t), intent(in)  :: this
    type(states_t),  intent(in)  :: st
    integer,         intent(out) :: ierr
    
    integer :: err, iunit
    
    PUSH_SUB(pes_spm_dump)

    err = 0 
    ierr = 0
    
    if (restart_skip(restart)) then
      POP_SUB(pes_spm_dump)
      return
    end if
    
    if (debug%info) then
      message(1) = "Debug: Writing PES_spm restart."
      call messages_info(1)
    end if

    if(this%onfly) then 
      call zrestart_write_binary(restart, 'pesspm', this%nspoints*st%d%dim*st%nst*st%d%nik*this%nomega, &
        this%wfft, err) 
    end if

    if(this%recipe == M_PHASE) then
      iunit = restart_open(restart, 'rcphase')
      write(iunit, *)  this%domega(:), this%dq(:,:)
      call restart_close(restart, iunit)
    end if

    if (err /= 0) ierr = ierr + 1
    
    if (debug%info) then
      message(1) = "Debug: Writing PES_spm restart done."
      call messages_info(1)
    end if
    
    POP_SUB(pes_spm_dump)
  end subroutine pes_spm_dump

  ! ---------------------------------------------------------
  subroutine pes_spm_load(restart, this, st, ierr)
    type(restart_t), intent(in)    :: restart    
    type(pes_spm_t), intent(inout) :: this
    type(states_t),  intent(inout) :: st
    integer,         intent(out)   :: ierr
    
    integer :: err, iunit
    
    PUSH_SUB(pes_spm_load)
    
    err = 0
    ierr = 0
    
    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(pes_spm_load)
      return
    end if
    
    if (debug%info) then
      message(1) = "Debug: Reading PES_spm restart."
      call messages_info(1)
    end if

    if(this%onfly) then
      call zrestart_read_binary(restart, 'pesspm', this%nspoints*st%d%dim*st%nst*st%d%nik*this%nomega, &
        this%wfft, err)
    end if

    if(this%recipe == M_PHASE) then
      iunit = restart_open(restart, 'rcphase')
      read(iunit, *)  this%domega(:), this%dq(:,:)
      call restart_close(restart, iunit)
    end if

    if (err /= 0) ierr = ierr + 1
    
    if(debug%info) then
      message(1) = "Debug: Reading PES_spm restart done."
      call messages_info(1)
    end if
    
    POP_SUB(pes_spm_load)
  end subroutine pes_spm_load 

  ! ---------------------------------------------------------
  subroutine pes_spm_calc_rcphase(this, mesh, iter, dt, hm, ii)
    type(pes_spm_t),     intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: ii
    
    integer :: mdim
    integer :: isp
    integer :: il, iprev
    FLOAT   :: vp(1:MAX_DIM)
    FLOAT   :: rr

    PUSH_SUB(pes_spm_calc_rcphase)

    mdim = mesh%sb%dim

    vp = M_ZERO
    do il = 1, hm%ep%no_lasers
      call laser_field(hm%ep%lasers(il), vp(1:mdim), iter*dt)
    end do
    vp(1:mdim) = -vp(1:mdim)

    iprev = ii - 1
    if(ii == 0) iprev = this%save_iter - 1

    do isp = 1, this%nspoints
      rr = sqrt(dot_product(this%rcoords(1:mdim, isp), this%rcoords(1:mdim, isp)))
      this%dq(isp, ii) = this%dq(isp, iprev) &
        + dot_product(this%rcoords(1:mdim, isp), vp(1:mdim)) / (P_C * rr) * dt
    end do

    this%domega(ii) = this%domega(iprev) &
      + dot_product(vp(1:mdim), vp(1:mdim)) / (M_TWO * P_C**M_TWO) * dt

    POP_SUB(pes_spm_calc_rcphase)
  end subroutine pes_spm_calc_rcphase

end module pes_spm_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
