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

module pes_rc_m
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

  private

  public ::                             &
    pes_rc_t,                           &
    pes_rc_init,                        &
    pes_rc_init_write,                  &
    pes_rc_output,                      &
    pes_rc_calc,                        &
    pes_rc_end,                         &
    pes_rc_dump,                        &
    pes_rc_load

  type PES_rc_t
    integer                    :: npoints                   !< how many points we store the wf
    integer, pointer           :: points(:)                 !< which points to use (local index)
    integer, pointer           :: points_global(:)          !< global index of the points
    FLOAT, pointer             :: coords(:,:)               !< coordinates of the sample points
    CMPLX, pointer             :: wf(:,:,:,:,:)   => NULL() !< wavefunctions at sample points
    integer, pointer           :: rankmin(:)                !< partition of the mesh containing the points
    FLOAT, pointer             :: dq(:,:)         => NULL() !< part 1 of Volkov phase (recipe phase) 
    FLOAT, pointer             :: domega(:)       => NULL() !< part 2 of Volkov phase (recipe phase)
    integer                    :: recipe                    !< type of calculation (RAW/PHASE)
    CMPLX, pointer             :: wfft(:,:,:,:,:) => NULL() !< Fourier transform of wavefunction
    FLOAT                      :: omegamax                  !< maximum frequency of the spectrum
    FLOAT                      :: delomega                  !< frequency spacing of the spectrum
    integer                    :: nomega                    !< number of frequencies of the spectrum
    logical                    :: onfly                     !< spectrum is calculated on-the-fly when true
    integer                    :: save_iter                 !< output interval and size of arrays 
    logical                    :: interpolation             !< use a an interpolated scheme for sample points
    integer                    :: nstepsphi, nstepstheta
    FLOAT                      :: thetamin
    type(mesh_interpolation_t) :: interp
  end type PES_rc_t

  integer, parameter :: &
    M_RAW   = 1,        &
    M_PHASE = 2

contains

  ! ---------------------------------------------------------
  subroutine PES_rc_init(pesrc, mesh, st, save_iter)
    type(PES_rc_t), intent(out) :: pesrc
    type(mesh_t),   intent(in)  :: mesh
    type(states_t), intent(in)  :: st
    integer,        intent(in)  :: save_iter

    type(block_t) :: blk
    integer       :: ip
    FLOAT         :: xx(MAX_DIM)
    FLOAT         :: dmin
    integer       :: rankmin
    logical       :: fromblk
    FLOAT         :: phi, theta, radius, default_radius
    integer       :: iph, ith

    PUSH_SUB(PES_rc_init)

    message(1) = 'Info: Calculating PES using rc technique.'
    call messages_info(1)

    !%Variable PhotoElectronSpectrumPoints
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% List of points at which to calculate the photoelectron spectrum by Suraud method.
    !% Required when <tt>PhotoElectronSpectrum = pes_rc</tt>.
    !% The exact syntax is:
    !%
    !% <tt>%PhotoElectronSpectrumPoints
    !% <br>&nbsp;&nbsp;x1 | y1 | z1
    !% <br>%
    !% </tt>
    !%End
    call messages_obsolete_variable('PES_rc_points', 'PhotoElectronSpectrumPoints')
    fromblk = .true.
    if (parse_block('PhotoElectronSpectrumPoints', blk) < 0) then
      fromblk = .false.
    end if

    if(fromblk) then
      message(1) = 'Info: Using PhotoElectronSpectrumPoints from block.'
    else
      message(1) = 'Info: Using spherical grid with interpolation.'
    end if
    call messages_info(1)

    !%Variable PES_rc_recipe
    !%Type integer
    !%Default raw
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
    call parse_variable('PES_rc_recipe', M_RAW, pesrc%recipe)
    if(.not.varinfo_valid_option('PES_rc_recipe', pesrc%recipe, is_flag = .true.)) &
      call messages_input_error('PES_rc_recipe')
    call messages_print_var_option(stdout, "PES_rc_recipe", pesrc%recipe)

    !%Variable PES_rc_OmegaMax
    !%Type float
    !%Default 0.0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% If PES_rc_OmegaMax > 0, the photoelectron spectrum is directly calculated during 
    !% time-propagation, evaluated by the PES_rc method. PES_rc_OmegaMax is then the maximum frequency 
    !% (approximate kinetic energy) and PES_rc_DelOmega the spacing in frequency domain of the spectrum.
    !%End
    call parse_variable('PES_rc_OmegaMax', units_to_atomic(units_inp%energy, M_ZERO), pesrc%omegamax)
    pesrc%onfly = .false.
    if(pesrc%omegamax > M_ZERO) then
      pesrc%onfly = .true.
      message(1) = 'Info: Calculating PES during time propagation.'
      call messages_info(1)
      call messages_print_var_value(stdout, "PES_rc_OmegaMax", pesrc%omegamax)
    end if
 
    !%Variable PES_rc_DelOmega
    !%Type float
    !%Default 0.005
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The spacing in frequency domain for the photoelectron spectrum (if PES_rc_OmegaMax > 0).
    !%End
    call parse_variable('PES_rc_DelOmega', units_to_atomic(units_inp%energy, CNST(0.005)), pesrc%delomega)
    if(pesrc%onfly) then
      if(pesrc%delomega <= M_ZERO) call messages_input_error('PES_rc_DelOmega')
      call messages_print_var_value(stdout, "PES_rc_DelOmega", pesrc%delomega)
    end if

    !%Variable PES_rc_ThetaSteps
    !%Type integer
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in theta (0 <= theta <= pi) for the spherical grid (if no 
    !% PhotoElectronSpectrumPoints are given).
    !%End
    call parse_variable('PES_rc_ThetaSteps', 45, pesrc%nstepstheta)
    if(.not.fromblk .and. pesrc%nstepstheta < 0) call messages_input_error('PES_rc_ThetaSteps')

    !%Variable PES_rc_PhiSteps
    !%Type integer
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in phi (0 <= phi <= 2 pi) for the spherical grid (if no
    !% PhotoElectronSpectrumPoints are given).
    !%End
    call parse_variable('PES_rc_PhiSteps', 90, pesrc%nstepsphi)
    if(.not.fromblk) then
      if(pesrc%nstepsphi < 0)  call messages_input_error('PES_rc_PhiSteps')
      if(pesrc%nstepsphi == 0) pesrc%nstepsphi = 1
    end if

    !%Variable PES_rc_Radius
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The radius of the sphere for the interpolation (if no PhotoElectronSpectrumPoints
    !% are given).
    !%End
    if(.not.fromblk) then
      select case(mesh%sb%box_shape)
      case(PARALLELEPIPED)
        default_radius = minval(mesh%sb%lsize(1:mesh%sb%dim))
      case(SPHERE)
        default_radius = mesh%sb%rsize
      case default
        message(1) = "Spherical grid not implemented for this box shape."
        message(2) = "Specify sample points with block PhotoElectronSpectrumPoints."
        call messages_fatal(2)
      end select
    end if
    call parse_variable('PES_rc_Radius', default_radius, radius)
    if(.not.fromblk) then
      if(radius <= M_ZERO) call messages_input_error('PES_rc_Radius')
      call messages_print_var_value(stdout, "PES_rc_Radius", radius)
    end if

    if(fromblk) then
      pesrc%interpolation = .false.
      pesrc%npoints = parse_block_n(blk)
    else
      pesrc%interpolation = .true.
      call mesh_interpolation_init(pesrc%interp, mesh)

      ! setting values for spherical grid 
      select case(mesh%sb%dim)
      case(1)
        pesrc%thetamin = M_PI / M_TWO
        pesrc%nstepstheta   = 0
        pesrc%nstepsphi     = 2
        pesrc%npoints  = pesrc%nstepsphi

      case(2)
        pesrc%thetamin = M_PI / M_TWO
        pesrc%nstepstheta   = 0
        pesrc%npoints  = pesrc%nstepsphi

        call messages_print_var_value(stdout, "PES_rc_StepsPhi", pesrc%nstepsphi)

      case(3)
        pesrc%thetamin = M_ZERO
        if(pesrc%nstepstheta <= 1) pesrc%nstepsphi = 1
        pesrc%npoints  = pesrc%nstepsphi * (pesrc%nstepstheta - 1) + 2

        call messages_print_var_value(stdout, "PES_rc_StepsPhi", pesrc%nstepsphi)
        call messages_print_var_value(stdout, "PES_rc_StepsTheta", pesrc%nstepstheta)

      end select
    end if

    call messages_print_var_value(stdout, "Number of PhotoElectronSpectrumPoints", pesrc%npoints)

    SAFE_ALLOCATE(pesrc%coords(1:mesh%sb%dim, 1:pesrc%npoints))

    if(fromblk) then
      SAFE_ALLOCATE(pesrc%points(1:pesrc%npoints))
      SAFE_ALLOCATE(pesrc%points_global(1:pesrc%npoints))
      SAFE_ALLOCATE(pesrc%rankmin(1:pesrc%npoints))
      pesrc%points_global = 0

      ! read points from input file
      do ip = 1, pesrc%npoints
        call parse_block_float(blk, ip - 1, 0, xx(1), units_inp%length)
        call parse_block_float(blk, ip - 1, 1, xx(2), units_inp%length)
        call parse_block_float(blk, ip - 1, 2, xx(3), units_inp%length)
        pesrc%coords(1:mesh%sb%dim, ip) = xx(1:mesh%sb%dim)
      end do
      call parse_block_end(blk)

      message(1) = 'Info: Calculating nearest points.'
      call messages_info(1)

      do ip = 1, pesrc%npoints
        ! nearest point
        pesrc%points(ip)  = mesh_nearest_point(mesh, pesrc%coords(1:mesh%sb%dim, ip), dmin, rankmin)
        pesrc%rankmin(ip) = rankmin

        if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
          if(mesh%mpi_grp%rank == rankmin) &
            pesrc%points_global(ip) = mesh%vp%local(mesh%vp%xlocal + pesrc%points(ip) - 1)
          call comm_allreduce(mesh%mpi_grp%comm, pesrc%points_global)
#endif
        else
          pesrc%points_global(ip) = pesrc%points(ip)
        end if
      end do

    else ! fromblk == .false.

      message(1) = 'Info: Initializing spherical grid.'
      call messages_info(1)

      ! initializing spherical grid
      ip = 0
      do ith = 0, pesrc%nstepstheta
        if(ith == 0) then
          theta = pesrc%thetamin     ! pi/2 for 1d and 2d. zero for 3d.
        else
          theta = ith * M_PI / pesrc%nstepstheta
        end if
        do iph = 0, pesrc%nstepsphi - 1
          if(iph == 0) then
            phi = M_ZERO
          else
            phi = iph * M_TWO * M_PI / pesrc%nstepsphi
          end if
          ip = ip + 1
                               pesrc%coords(1, ip) = radius * cos(phi) * sin(theta)
          if(mesh%sb%dim >= 2) pesrc%coords(2, ip) = radius * sin(phi) * sin(theta)
          if(mesh%sb%dim == 3) pesrc%coords(3, ip) = radius * cos(theta)
          if(theta == M_ZERO .or. theta == M_PI) exit
        end do
      end do

    end if

    SAFE_ALLOCATE(pesrc%wf(1:st%nst, 1:st%d%dim, 1:st%d%nik, 1:pesrc%npoints, 0:save_iter-1))

    if(pesrc%recipe == M_PHASE) then
      SAFE_ALLOCATE(pesrc%dq(1:pesrc%npoints, 0:save_iter-1))
      SAFE_ALLOCATE(pesrc%domega(0:save_iter-1))
      pesrc%dq = M_ZERO
      pesrc%domega = M_ZERO
    end if

    if(pesrc%onfly) then
      pesrc%nomega = nint(pesrc%omegamax/pesrc%delomega)
      SAFE_ALLOCATE(pesrc%wfft(1:st%nst, 1:st%d%dim, 1:st%d%nik, 1:pesrc%npoints, 1:pesrc%nomega))
      pesrc%wfft = M_z0
    end if

    pesrc%save_iter = save_iter

    POP_SUB(PES_rc_init)
  end subroutine PES_rc_init


  ! ---------------------------------------------------------
  subroutine PES_rc_end(pesrc)
    type(PES_rc_t), intent(inout) :: pesrc

    PUSH_SUB(PES_rc_end)

    SAFE_DEALLOCATE_P(pesrc%points)
    SAFE_DEALLOCATE_P(pesrc%points_global)
    SAFE_DEALLOCATE_P(pesrc%wf)
    SAFE_DEALLOCATE_P(pesrc%rankmin)
    SAFE_DEALLOCATE_P(pesrc%coords)

    SAFE_DEALLOCATE_P(pesrc%wfft)

    SAFE_DEALLOCATE_P(pesrc%dq)
    SAFE_DEALLOCATE_P(pesrc%domega)

    POP_SUB(PES_rc_end)
  end subroutine PES_rc_end


  ! ---------------------------------------------------------
  subroutine PES_rc_calc(pesrc, st, mesh, dt, iter, hm)
    type(PES_rc_t),      intent(inout) :: pesrc
    type(states_t),      intent(in)    :: st
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter
    type(hamiltonian_t), intent(in)    :: hm

    CMPLX, allocatable :: psi(:,:,:,:), psistate(:), wfftact(:,:,:,:,:)
    integer            :: ip, ii
    integer            :: dim, stst, stend, kptst, kptend
    integer            :: ik, ist, idim
    logical            :: contains_ip
    CMPLX              :: rawfac
    CMPLX, allocatable :: phasefac(:)
    FLOAT              :: omega
    integer            :: iom
    CMPLX, allocatable :: interp_values(:)

    PUSH_SUB(PES_rc_calc)

    ii = mod(iter-1, pesrc%save_iter)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    dim    = st%d%dim

    if(pesrc%interpolation) then
      SAFE_ALLOCATE(psistate(1:mesh%np_part))
      SAFE_ALLOCATE(interp_values(1:pesrc%npoints))
    else
      SAFE_ALLOCATE(psi(1:st%nst, dim, 1:1, 1:st%d%nik))
    end if

    if(pesrc%onfly) then
      SAFE_ALLOCATE(wfftact(1:st%nst, dim, 1:st%d%nik, 1:pesrc%npoints, 1:pesrc%nomega))
      wfftact = M_z0
    endif

    if(pesrc%recipe == M_PHASE) then
      SAFE_ALLOCATE(phasefac(1:pesrc%npoints))
    end if

    ! needed for allreduce, otherwise it will take values from previous cycle
    pesrc%wf(:,:,:,:,ii) = M_z0

    if(pesrc%interpolation) then
      do ik = kptst, kptend 
        do ist = stst, stend
          do idim = 1, dim
            call states_get_state(st, mesh, idim, ist, ik, psistate(1:mesh%np_part))
            call mesh_interpolation_evaluate(pesrc%interp, pesrc%npoints, psistate(1:mesh%np_part), &
              pesrc%coords(1:mesh%sb%dim, 1:pesrc%npoints), interp_values(1:pesrc%npoints))
            pesrc%wf(ist, idim, ik, :, ii) = interp_values(:)
          end do
        end do
      end do
      if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
        ! interpolated values have already been communicated over domains
        call comm_allreduce(st%st_kpt_mpi_grp%comm, pesrc%wf(:,:,:,:,ii))
#endif
      end if
    else ! pesrc%interpolation == .false.

      contains_ip = .true.

      do ip = 1, pesrc%npoints

#if defined(HAVE_MPI)
        if(mesh%parallel_in_domains) then
          if(mesh%mpi_grp%rank  ==  pesrc%rankmin(ip)) then ! needed if mesh%parallel_in_domains is true
            contains_ip = .true.
          else
            contains_ip = .false.
          end if
        end if
#endif

        if(contains_ip) then
          call states_get_points(st, pesrc%points(ip), pesrc%points(ip), psi)
          pesrc%wf(stst:stend, dim, kptst:kptend, ip, ii) = psi(stst:stend, dim, 1, kptst:kptend)
        end if
      end do
#if defined(HAVE_MPI)
      call comm_allreduce(mpi_world%comm, pesrc%wf(:,:,:,:,ii))
#endif
    end if

    if(pesrc%recipe == M_PHASE) then
      call pes_rc_calc_rcphase(pesrc, mesh, iter, dt, hm, ii)
    end if

    if(pesrc%onfly) then
      do iom = 1, pesrc%nomega
        omega = iom*pesrc%delomega
        rawfac = exp(M_zI * omega * iter * dt) * dt

        select case(pesrc%recipe)
        case(M_RAW)
          wfftact(stst:stend, 1:dim, kptst:kptend, :, iom) = pesrc%wf(stst:stend, 1:dim, kptst:kptend, :, ii) * rawfac
        case(M_PHASE)
          phasefac(:) = exp(M_zI * (pesrc%domega(ii) - sqrt(M_TWO * omega) * pesrc%dq(:, ii)))

          do ik = kptst, kptend
            do ist = stst, stend
              do idim = 1, dim
                wfftact(ist, idim, ik, :, iom) = pesrc%wf(ist, idim, ik, :, ii) * phasefac(:) * rawfac
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

      pesrc%wfft = pesrc%wfft + wfftact
    end if

    SAFE_DEALLOCATE_A(psistate)
    SAFE_DEALLOCATE_A(interp_values)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(phasefac)

    SAFE_DEALLOCATE_A(wfftact)

    POP_SUB(PES_rc_calc)
  end subroutine PES_rc_calc


  ! ---------------------------------------------------------
  subroutine PES_rc_output(pesrc, st, iter, dt)
    type(PES_rc_t), intent(in) :: pesrc
    type(states_t), intent(in) :: st
    integer,        intent(in) :: iter
    FLOAT,          intent(in) :: dt

    integer            :: ip, ii, jj, ik, ist, idim, iom, iph, ith, iphi
    integer            :: iunit, iunitone, iunittwo
    CMPLX              :: vfu
    character(len=4)   :: filenr
    FLOAT              :: wfu, omega
    FLOAT, allocatable :: wffttot(:,:)
    FLOAT              :: wffttotsave
    integer            :: save_iter
    FLOAT              :: theta, phi

    PUSH_SUB(PES_rc_output)

    save_iter = pesrc%save_iter

    if(pesrc%onfly) then
      SAFE_ALLOCATE(wffttot(1:pesrc%nomega, 1:pesrc%npoints))
      wffttot = M_ZERO

      ! calculate total spectrum
      do iom = 1, pesrc%nomega
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            do idim = 1, st%d%dim
              wffttot(iom, :) = real(pesrc%wfft(ist, idim, ik, :, iom))**2 + aimag(pesrc%wfft(ist, idim, ik,:, iom))**2 &
                + wffttot(iom, :)
            end do
          end do
        end do
      end do
    end if

    if(.not.pesrc%interpolation .or. in_debug_mode) then   ! too much output for interpolation
      do ip = 1, pesrc%npoints
        write(filenr, '(i4.4)') ip
   
        iunitone = io_open('td.general/'//'PES_rc.'//filenr//'.wavefunctions.out', action='write', position='append')
        if(pesrc%recipe == M_PHASE) &
          iunittwo = io_open('td.general/'//'PES_rc.'//filenr//'.phase.out', action='write', position='append')
   
        do ii = 1, save_iter - mod(iter, save_iter)
          jj = iter - save_iter + ii + mod(save_iter - mod(iter, save_iter), save_iter)
          write(iunitone, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do idim = 1, st%d%dim
                vfu = units_from_atomic(sqrt(units_out%length**(-3)), pesrc%wf(ist, idim, ik, ip, ii-1))
                write(iunitone, '(1x,e18.10E3,1x,e18.10E3)', advance='no') &
                  real(vfu),  aimag(vfu) 
              end do
            end do
          end do
          write(iunitone, '(1x)', advance='yes')

          if(pesrc%recipe == M_PHASE) then
            write(iunittwo, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
            write(iunittwo, '(1x,e18.10E3,1x,e18.10E3)', advance='no') &
              pesrc%domega(ii-1), pesrc%dq(ip, ii-1)
            write(iunittwo, '(1x)', advance='yes')
          end if
        end do
        call io_close(iunitone)
        if(pesrc%recipe == M_PHASE) call io_close(iunittwo)

        if(pesrc%onfly) then
          iunit = io_open('td.general/'//'PES_rc.'//filenr//'.spectrum.out', action='write', position='rewind')
          write(iunit, '(a44)') '# frequency, total spectrum, orbital spectra'
          do iom = 1, pesrc%nomega 
            omega = iom*pesrc%delomega
            write(iunit, '(e17.10, 1x, e17.10)', advance='no') omega, wffttot(iom, ip)
            do ik = 1, st%d%nik
              do ist = 1, st%nst 
                do idim = 1, st%d%dim
                  wfu = real(pesrc%wfft(ist, idim, ik, ip, iom))**2 + aimag(pesrc%wfft(ist, idim, ik, ip, iom))**2
                  write(iunit,'(1x,e18.10e3)', advance='no') wfu
                end do
              end do
            end do
          write(iunit,'(1x)', advance='yes')
          end do
          call io_close(iunit)
        end if
   
      end do
    end if

    if(pesrc%onfly .and. pesrc%interpolation) then
      iunit = io_open('td.general/'//'PES_rc.distribution.out', action='write', position='rewind')
      write(iunit, '(a33)') '# omega, theta, phi, distribution'

      do iom = 1, pesrc%nomega
        omega = iom * pesrc%delomega
        ip = 0
        do ith = 0, pesrc%nstepstheta
          if(ith == 0) then
            theta = pesrc%thetamin     ! pi/2 for 1d and 2d. zero for 3d.
          else
            theta = ith * M_PI / pesrc%nstepstheta
          end if
          do iph = 0, pesrc%nstepsphi - 1
            ip = ip + 1
            if(iph == 0) then
              wffttotsave = wffttot(iom, ip)
              phi = M_ZERO
            else
              phi = iph * M_TWO * M_PI / pesrc%nstepsphi
            end if
            write(iunit,'(5(1x,e18.10E3))') omega, theta, phi, wffttot(iom, ip)

            ! just repeat the result for output
            if(iph == (pesrc%nstepsphi - 1)) then
              write(iunit,'(5(1x,e18.10E3))') omega, theta, M_TWO * M_PI, wffttotsave
            end if
             
            ! just repeat the result for output
            if(theta == M_ZERO .or. theta == M_PI) then
              do iphi = 1, pesrc%nstepsphi
                phi = iphi * M_TWO * M_PI / pesrc%nstepsphi
                write(iunit,'(5(1x,e18.10E3))') omega, theta, phi, wffttot(iom, ip)
              end do
              exit
            end if
          end do
          if(pesrc%nstepsphi > 0) write(iunit, '(1x)', advance='yes')
        end do
        if(pesrc%nstepsphi == 0) write(iunit, '(1x)', advance='yes')
      end do
      call io_close(iunit)

    end if

    SAFE_DEALLOCATE_A(wffttot)

    POP_SUB(PES_rc_output)
  end subroutine PES_rc_output

  ! ---------------------------------------------------------
  subroutine PES_rc_init_write(pesrc, mesh, st)
    type(PES_rc_t), intent(in) :: pesrc
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(in) :: st

    integer          :: ip,ik,ist,idim,iunit
    FLOAT            :: xx(MAX_DIM)
    character(len=4) :: filenr

    PUSH_SUB(PES_rc_init_write)

    xx = M_ZERO
    if(mpi_grp_is_root(mpi_world)) then
      if(.not.pesrc%interpolation .or. in_debug_mode) then   ! too much output for interpolation
        do ip = 1, pesrc%npoints
          write(filenr, '(i4.4)') ip
   
          iunit = io_open('td.general/'//'PES_rc.'//filenr//'.wavefunctions.out', action='write')
          xx(1:mesh%sb%dim) = pesrc%coords(1:mesh%sb%dim, ip)
          write(iunit,'(a1)') '#'
          write(iunit, '(a7,f17.6,a1,f17.6,a1,f17.6,5a)') &
            '# R = (',xx(1),' ,',xx(2),' ,',xx(3), &
            ' )  [', trim(units_abbrev(units_inp%length)), ']'
   
          write(iunit,'(a1)') '#'  
          write(iunit, '(a3,14x)', advance='no') '# t' 
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do idim = 1, st%d%dim
                write(iunit, '(3x,a8,i3,a7,i3,a8,i3,3x)', advance='no') &
                  "ik = ", ik, " ist = ", ist, " idim = ", idim
              end do
            end do
          end do
          write(iunit, '(1x)', advance='yes')
   
          call io_close(iunit)
   
          if(pesrc%recipe == M_PHASE) then
            iunit = io_open('td.general/'//'PES_rc.'//filenr//'.phase.out', action='write')
            write(iunit,'(a24)') '# time, dq(t), dOmega(t)'
            call io_close(iunit)
          end if
        end do
      end if
    end if

    POP_SUB(PES_rc_init_write)
  end subroutine PES_rc_init_write

  ! ---------------------------------------------------------
  subroutine pes_rc_dump(restart, pesrc, st, ierr)
    type(restart_t), intent(in)  :: restart    
    type(pes_rc_t),  intent(in)  :: pesrc
    type(states_t),  intent(in)  :: st
    integer,         intent(out) :: ierr
    
    integer :: err, iunit
    
    PUSH_SUB(pes_rc_dump)

    ierr = 0
    
    if (restart_skip(restart)) then
      POP_SUB(pes_rc_dump)
      return
    end if
    
    if (in_debug_mode) then
      message(1) = "Debug: Writing pes_rc restart."
      call messages_info(1)
    end if

    if(pesrc%onfly) then 
      call zrestart_write_binary(restart, 'pesrc', pesrc%npoints*st%d%dim*st%nst*st%d%nik*pesrc%nomega, &
        pesrc%wfft, err) 
    end if

    if(pesrc%recipe == M_PHASE) then
      iunit = restart_open(restart, 'rcphase')
      write(iunit, *)  pesrc%domega(:), pesrc%dq(:,:)
      call restart_close(restart, iunit)
    end if

    if (err /= 0) ierr = ierr + 1
    
    if (in_debug_mode) then
      message(1) = "Debug: Writing pes_rc restart done."
      call messages_info(1)
    end if
    
    POP_SUB(pes_rc_dump)
  end subroutine pes_rc_dump

  ! ---------------------------------------------------------
  subroutine pes_rc_load(restart, pesrc, st, ierr)
    type(restart_t), intent(in)    :: restart    
    type(pes_rc_t),  intent(inout) :: pesrc
    type(states_t),  intent(inout) :: st
    integer,         intent(out)   :: ierr
    
    integer :: err, iunit
    
    PUSH_SUB(pes_rc_load)
    
    ierr = 0
    
    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(pes_rc_load)
      return
    end if
    
    if (in_debug_mode) then
      message(1) = "Debug: Reading pes_rc restart."
      call messages_info(1)
    end if

    if(pesrc%onfly) then
      call zrestart_read_binary(restart, 'pesrc', pesrc%npoints*st%d%dim*st%nst*st%d%nik*pesrc%nomega, &
        pesrc%wfft, err)
    end if

    if(pesrc%recipe == M_PHASE) then
      iunit = restart_open(restart, 'rcphase')
      read(iunit, *)  pesrc%domega(:), pesrc%dq(:,:)
      call restart_close(restart, iunit)
    end if

    if (err /= 0) ierr = ierr + 1
    
    if(in_debug_mode) then
      message(1) = "Debug: Reading pes_rc restart done."
      call messages_info(1)
    end if
    
    POP_SUB(pes_rc_load)
  end subroutine pes_rc_load 

  ! ---------------------------------------------------------
  subroutine PES_rc_calc_rcphase(pesrc, mesh, iter, dt, hm, ii)
    type(PES_rc_t),      intent(inout) :: pesrc
    type(mesh_t),        intent(in)    :: mesh
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: ii
    
    integer :: dim, il, ip
    FLOAT   :: vp(1:MAX_DIM)
    FLOAT   :: xx(MAX_DIM), er(MAX_DIM)

    PUSH_SUB(PES_rc_calc_rcphase)

    dim = mesh%sb%dim

    vp = M_ZERO
    do il = 1, hm%ep%no_lasers
      call laser_field(hm%ep%lasers(il), vp(1:dim), iter*dt)
    end do
    vp(1:dim) = -vp(1:dim)

    do ip = 1, pesrc%npoints
      xx(1:dim) = pesrc%coords(1:dim, ip)
      er(1:dim) = xx(1:dim)/sqrt(dot_product(xx(1:dim), xx(1:dim)))
      pesrc%dq(ip, ii) = pesrc%dq(ip, ii) + dot_product(er(1:dim), vp(1:dim))/P_C * dt
    end do

    pesrc%domega(ii) = pesrc%domega(ii) + dot_product(vp(1:dim), vp(1:dim)) / (M_TWO * P_C**M_TWO) * dt

    POP_SUB(PES_rc_calc_rcphase)
  end subroutine PES_rc_calc_rcphase

end module pes_rc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
