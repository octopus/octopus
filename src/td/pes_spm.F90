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

module pes_spm_oct_m
  use box_parallelepiped_oct_m
  use box_sphere_oct_m
  use comm_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use lasers_oct_m
  use mesh_interpolation_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m


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
    private
    integer                    :: nspoints                  !< how many points we store the wf
    FLOAT, allocatable         :: rcoords(:,:)              !< coordinates of the sample points
    FLOAT, allocatable         :: rcoords_nrm(:,:)
    CMPLX, allocatable, public :: wf(:,:,:,:,:)             !< wavefunctions at sample points
    FLOAT, allocatable         :: dq(:,:)                   !< part 1 of Volkov phase (recipe phase) 
    FLOAT, allocatable         :: domega(:)                 !< part 2 of Volkov phase (recipe phase)
    integer                    :: recipe                    !< type of calculation (RAW/PHASE)
    CMPLX, allocatable         :: wfft(:,:,:,:,:)           !< Fourier transform of wavefunction
    FLOAT                      :: omegamax                  !< maximum frequency of the spectrum
    FLOAT                      :: delomega                  !< frequency spacing of the spectrum
    integer                    :: nomega                    !< number of frequencies of the spectrum
    logical                    :: onfly                     !< spectrum is calculated on-the-fly when true
    integer                    :: save_iter                 !< output interval and size of arrays 
    logical                    :: sphgrid                   !< use a spherical grid (instead of sample points from input)
    integer                    :: nstepsphir, nstepsthetar
    type(mesh_interpolation_t) :: interp
  end type pes_spm_t

  integer, parameter :: &
    M_RAW   = 1,        &
    M_PHASE = 2

contains

  ! ---------------------------------------------------------
  subroutine pes_spm_init(this, namespace, mesh, st, save_iter)
    type(pes_spm_t),      intent(out) :: this
    type(namespace_t),    intent(in)  :: namespace
    type(mesh_t),         intent(in)  :: mesh
    type(states_elec_t),  intent(in)  :: st
    integer,              intent(in)  :: save_iter

    type(block_t) :: blk
    integer       :: stst, stend, kptst, kptend, sdim, mdim
    integer       :: imdim
    integer       :: isp
    integer       :: ith, iph
    FLOAT         :: thetar, phir, radius
    FLOAT         :: xx(MAX_DIM)

    PUSH_SUB(pes_spm_init)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = mesh%sb%dim

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
    call messages_obsolete_variable(namespace, 'PhotoElectronSpectrumPoints', 'PES_spm_points')
    this%sphgrid = .false.
    if (parse_block(namespace, 'PES_spm_points', blk) < 0) then
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
    call parse_variable(namespace, 'PES_spm_recipe', M_PHASE, this%recipe)
    if(.not.varinfo_valid_option('PES_spm_recipe', this%recipe, is_flag = .true.)) then
      call messages_input_error(namespace, 'PES_spm_recipe')
    end if
    call messages_print_var_option(stdout, "PES_spm_recipe", this%recipe)

    !%Variable PES_spm_OmegaMax
    !%Type float
    !%Default 0.0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% If <tt>PES_spm_OmegaMax > 0</tt>, the photoelectron spectrum is directly calculated during
    !% time-propagation, evaluated by the PES_spm method. <tt>PES_spm_OmegaMax</tt> is then the maximum frequency
    !% (approximate kinetic energy) and <tt>PES_spm_DeltaOmega</tt> the spacing in frequency domain of the spectrum.
    !%End
    call parse_variable(namespace, 'PES_spm_OmegaMax', units_to_atomic(units_inp%energy, M_ZERO), this%omegamax)
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
    !% The spacing in frequency domain for the photoelectron spectrum (if <tt>PES_spm_OmegaMax > 0</tt>).
    !% The default is <tt>PES_spm_OmegaMax/500</tt>.
    !%End
    call parse_variable(namespace, 'PES_spm_DeltaOmega', units_to_atomic(units_inp%energy, this%omegamax/CNST(500)), this%delomega)
    if(this%onfly) then
      if(this%delomega <= M_ZERO) call messages_input_error(namespace, 'PES_spm_DeltaOmega')
      call messages_print_var_value(stdout, "PES_spm_DeltaOmega", this%delomega)
    end if

    !%Variable PES_spm_StepsThetaR
    !%Type integer
    !%Default 45
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in <math>\theta</math> (<math>0 \le \theta \le \pi</math>) for the spherical grid (if no
    !% <tt>PES_spm_points</tt> are given).
    !%End
    call parse_variable(namespace, 'PES_spm_StepsThetaR', 45, this%nstepsthetar)
    if(this%sphgrid .and. this%nstepsthetar < 0) call messages_input_error(namespace, 'PES_spm_StepsThetaR')

    !%Variable PES_spm_StepsPhiR
    !%Type integer
    !%Default 90
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in <math>\phi</math> (<math>0 \le \phi \le 2 \pi</math>) for the spherical grid (if no
    !% <tt>PES_spm_points</tt> are given).
    !%End
    call parse_variable(namespace, 'PES_spm_StepsPhiR', 90, this%nstepsphir)
    if(this%sphgrid) then
      if(this%nstepsphir < 0)  call messages_input_error(namespace, 'PES_spm_StepsPhiR')
      if(this%nstepsphir == 0) this%nstepsphir = 1
    end if

    !%Variable PES_spm_Radius
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The radius of the sphere for the spherical grid (if no <tt>PES_spm_points</tt>
    !% are given).
    !%End
    if(this%sphgrid) then
      if(parse_is_defined(namespace, 'PES_spm_Radius')) then
        call parse_variable(namespace, 'PES_spm_Radius', M_ZERO, radius)
        if(radius <= M_ZERO) call messages_input_error(namespace, 'PES_spm_Radius')
        call messages_print_var_value(stdout, "PES_spm_Radius", radius)
      else
        select type (box => mesh%sb%box)
        type is (box_sphere_t)
          radius = box%radius
        type is (box_parallelepiped_t)
          radius = minval(box%half_length(1:mdim))
        class default
          message(1) = "Spherical grid not implemented for this box shape."
          message(2) = "Specify sample points with block PES_spm_points."
          call messages_fatal(2, namespace=namespace)
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
        this%nstepsthetar   = 0
        this%nstepsphir     = 2
        this%nspoints  = this%nstepsphir

      case(2)
        this%nstepsthetar   = 0
        this%nspoints  = this%nstepsphir

      case(3)
        if(this%nstepsthetar <= 1) then
          this%nstepsphir = 1
          this%nstepsthetar = 1
        end if
        this%nspoints  = this%nstepsphir * (this%nstepsthetar - 1) + 2

      end select

      if(mdim == 3) call messages_print_var_value(stdout, "PES_spm_StepsThetaR", this%nstepsthetar)
      call messages_print_var_value(stdout, "PES_spm_StepsPhiR", this%nstepsphir)

    end if

    call messages_print_var_value(stdout, "Number of sample points", this%nspoints)

    SAFE_ALLOCATE(this%rcoords(1:mdim, 1:this%nspoints))
    SAFE_ALLOCATE(this%rcoords_nrm(1:mdim, 1:this%nspoints))

    if(.not. this%sphgrid) then

      message(1) = 'Info: Reading sample points from input.'
      call messages_info(1)

      ! read points from input file
      do isp = 1, this%nspoints
        do imdim = 1, mdim
          call parse_block_float(blk, isp - 1, imdim - 1, xx(imdim), units_inp%length)
        end do
        this%rcoords(1:mdim, isp) = xx(1:mdim)
        this%rcoords_nrm(1:mdim, isp) = this%rcoords(1:mdim, isp) / &
          sqrt(dot_product(this%rcoords(:, isp), this%rcoords(:, isp)))
      end do
      call parse_block_end(blk)

    else ! this%sphgrid == .true.

      message(1) = 'Info: Initializing spherical grid.'
      call messages_info(1)

      ! initializing spherical grid
      thetar = M_PI / M_TWO
      isp = 0
      do ith = 0, this%nstepsthetar
        if(mdim == 3) thetar = ith * M_PI / this%nstepsthetar
        do iph = 0, this%nstepsphir - 1
          isp = isp + 1
          phir = iph * M_TWO * M_PI / this%nstepsphir
                        this%rcoords_nrm(1, isp) = cos(phir) * sin(thetar)
          if(mdim >= 2) this%rcoords_nrm(2, isp) = sin(phir) * sin(thetar)
          if(mdim == 3) this%rcoords_nrm(3, isp) = cos(thetar)
          this%rcoords(1:mdim, isp) = radius * this%rcoords_nrm(1:mdim, isp)
          if(mdim == 3 .and. (ith == 0 .or. ith == this%nstepsthetar)) exit
        end do
      end do

    end if

    SAFE_ALLOCATE(this%wf(stst:stend, 1:sdim, kptst:kptend, 1:this%nspoints, 0:save_iter-1))

    if(this%recipe == M_PHASE) then
      SAFE_ALLOCATE(this%dq(1:this%nspoints, 0:save_iter-1))
      SAFE_ALLOCATE(this%domega(0:save_iter-1))
      this%dq = M_ZERO
      this%domega = M_ZERO
    end if

    if(this%onfly) then
      this%nomega = nint(this%omegamax/this%delomega)
      SAFE_ALLOCATE(this%wfft(stst:stend, 1:sdim, kptst:kptend, 1:this%nspoints, 1:this%nomega))
      this%wfft = M_z0
    end if

    this%save_iter = save_iter

    POP_SUB(pes_spm_init)
  end subroutine pes_spm_init


  ! ---------------------------------------------------------
  subroutine pes_spm_end(this)
    type(pes_spm_t), intent(inout) :: this

    PUSH_SUB(pes_spm_end)

    SAFE_DEALLOCATE_A(this%wf)
    SAFE_DEALLOCATE_A(this%rcoords)
    SAFE_DEALLOCATE_A(this%rcoords_nrm)

    SAFE_DEALLOCATE_A(this%wfft)

    SAFE_DEALLOCATE_A(this%dq)
    SAFE_DEALLOCATE_A(this%domega)

    POP_SUB(pes_spm_end)
  end subroutine pes_spm_end


  ! ---------------------------------------------------------
  subroutine pes_spm_calc(this, st, mesh, dt, iter, hm)
    type(pes_spm_t),     intent(inout) :: this
    type(states_elec_t), intent(in)    :: st
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter
    type(hamiltonian_elec_t), intent(in)    :: hm

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim
    integer            :: itstep

    CMPLX, allocatable :: psistate(:), wfftact(:,:,:,:,:)
    CMPLX              :: rawfac
    CMPLX, allocatable :: phasefac(:)
    integer            :: iom
    FLOAT              :: omega
    CMPLX, allocatable :: interp_values(:)

    PUSH_SUB(pes_spm_calc)

    itstep = mod(iter-1, this%save_iter)

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
    this%wf(:,:,:,:,itstep) = M_z0

    do ik = kptst, kptend 
      do ist = stst, stend
        do isdim = 1, sdim
          call states_elec_get_state(st, mesh, isdim, ist, ik, psistate(1:mesh%np_part))
          call mesh_interpolation_evaluate(this%interp, this%nspoints, psistate(1:mesh%np_part), &
            this%rcoords(1:mdim, 1:this%nspoints), interp_values(1:this%nspoints))
          this%wf(ist, isdim, ik, :, itstep) = st%occ(ist, ik) * interp_values(:)
        end do
      end do
    end do

    if(this%recipe == M_PHASE) then
      call pes_spm_calc_rcphase(this, mesh, iter, dt, hm, itstep)
    end if

    if(this%onfly) then
      do iom = 1, this%nomega
        omega = iom*this%delomega
        rawfac = exp(M_zI * omega * iter * dt) * dt / (M_TWO * M_PI)**(mdim/M_TWO)

        if(this%recipe == M_RAW) then
          wfftact(stst:stend, 1:sdim, kptst:kptend, :, iom) = &
            rawfac * sqrt(M_TWO * omega) * this%wf(stst:stend, 1:sdim, kptst:kptend, :, itstep)
        else
          phasefac(:) = rawfac * exp(M_zI * (this%domega(itstep) - sqrt(M_TWO * omega) * this%dq(:, itstep)))

          do ik = kptst, kptend
            do ist = stst, stend
              do isdim = 1, sdim
                wfftact(ist, isdim, ik, :, iom) = phasefac(:) * this%wf(ist, isdim, ik, :, itstep) * sqrt(M_TWO * omega)
              end do
            end do
          end do
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
  subroutine pes_spm_output(this, mesh, st, namespace, iter, dt)
    type(pes_spm_t),     intent(in) :: this
    type(mesh_t),        intent(in) :: mesh
    type(states_elec_t), intent(in) :: st
    type(namespace_t),   intent(in) :: namespace
    integer,             intent(in) :: iter
    FLOAT,               intent(in) :: dt

    integer            :: ist, ik, isdim
    integer            :: ii, jj
    integer            :: isp, save_iter, isp_save
    integer            :: iom, ith, iph, iphi, itot
    FLOAT              :: omega, thetar, phir
    CMPLX              :: vfu
    FLOAT              :: weight
    FLOAT, allocatable :: spctrsum(:,:,:,:), spctrout(:,:)
    character(len=80)  :: filenr
    integer            :: iunitone, iunittwo
    integer            :: stst, stend, kptst, kptend, sdim, mdim

    PUSH_SUB(pes_spm_output)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = mesh%sb%dim

    save_iter = this%save_iter

    if(this%onfly) then
      SAFE_ALLOCATE(spctrsum(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nomega))
      spctrsum = M_ZERO

      SAFE_ALLOCATE(spctrout(1:this%nspoints, 1:this%nomega))
      spctrout = M_ZERO

      do ik = kptst, kptend 
        do ist = stst, stend
          do isdim = 1, sdim

            if(this%sphgrid) then

              select case(mdim)
              case(1)
                weight = M_HALF
   
                do iom = 1, this%nomega
                  spctrsum(ist, isdim, ik, iom) = spctrsum(ist, isdim, ik, iom) + &
                    sum(abs(this%wfft(ist, isdim, ik, :, iom)**M_TWO),1) * weight
                end do
   
              case(2)
                weight = M_TWO * M_PI / this%nstepsphir
   
                do iom = 1, this%nomega
                  do iph = 0, this%nstepsphir - 1
                    spctrsum(ist, isdim, ik, iom) = spctrsum(ist, isdim, ik, iom) + &
                      abs(this%wfft(ist, isdim, ik, iph, iom))**M_TWO * weight
                  end do
                end do
   
              case(3)
                do iom = 1, this%nomega
                  isp = 0
       
                  do ith = 0, this%nstepsthetar
                    thetar = ith * M_PI / this%nstepsthetar
       
                    if(ith == 0 .or. ith == this%nstepsthetar) then
                      weight = (M_ONE - cos(M_PI / this%nstepsthetar / M_TWO)) * M_TWO * M_PI
                    else
                      weight = abs(cos(thetar - M_PI / this%nstepsthetar / M_TWO) - &
                        cos(thetar + M_PI / this%nstepsthetar / M_TWO)) * M_TWO * M_PI / this%nstepsphir
                    end if
       
                    do iph = 0, this%nstepsphir - 1
                      isp = isp + 1
                      spctrsum(ist, isdim, ik, iom) = spctrsum(ist, isdim, ik, iom) + &
                        abs(this%wfft(ist, isdim, ik, isp, iom))**M_TWO * weight
       
                      if(ith == 0 .or. ith == this%nstepsthetar) exit
                    end do
                  end do
                end do
              end select

              ! distribution
              spctrout(1:this%nspoints, 1:this%nomega) = spctrout(1:this%nspoints, 1:this%nomega) + &
                abs(this%wfft(ist, isdim, ik, 1:this%nspoints, 1:this%nomega))**M_TWO

            end if
          end do
        end do
      end do

      if(st%parallel_in_states .or. st%d%kpt%parallel) then
        ! total spectrum = sum over all states
        call comm_allreduce(st%st_kpt_mpi_grp, spctrout)

        ! orbital spectra
        call comm_allreduce(st%st_kpt_mpi_grp, spctrsum)
      end if

      ! -----------------------------------------------------------------
      ! OUTPUT FOR SPHERICAL GRID 
      ! -----------------------------------------------------------------
      if(this%sphgrid) then
        iunittwo = io_open('td.general/'//'PES_spm.distribution.out', namespace, action='write', position='rewind')
        iunitone = io_open('td.general/'//'PES_spm.power.sum', namespace, action='write', position='rewind')
        write(iunitone, '(a23)') '# omega, total spectrum'
   
        select case(mdim)
        case(1)
          write(iunittwo, '(a40)') '# omega, distribution (left/right point)'
   
          do iom = 1, this%nomega
            omega = iom * this%delomega
            write(iunittwo, '(5(1x,e18.10E3))') omega, spctrout(2, iom), spctrout(1, iom)
            write(iunitone, '(2(1x,e18.10E3))') omega, sum(sum(sum(spctrsum(:,:,:,iom),1),1),1) * sqrt(M_TWO * omega)
          end do
   
        case(2)
          write(iunittwo, '(a26)') '# omega, phi, distribution'
   
          do iom = 1, this%nomega
            omega = iom * this%delomega
   
            do iph = 0, this%nstepsphir - 1
              phir = iph * M_TWO * M_PI / this%nstepsphir
              write(iunittwo,'(5(1x,e18.10E3))') omega, phir, spctrout(iph + 1, iom)
            end do
            ! just repeat the result for output
            write(iunittwo,'(5(1x,e18.10E3))') omega, M_TWO * M_PI, spctrout(1, iom)
            write(iunittwo, '(1x)', advance='yes')
   
            write(iunitone, '(2(1x,e18.10E3))', advance='no') omega, sum(sum(sum(spctrsum(:,:,:,iom),1),1),1) * sqrt(M_TWO * omega)
   
            do ik = 1, st%d%nik
              do ist = 1, st%nst
                do isdim = 1, st%d%dim
                  write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, iom) * sqrt(M_TWO * omega)
                end do
              end do
            end do
            write(iunitone, '(1x)', advance='yes')
                                                                                                    
          end do
   
        case(3)
          write(iunittwo, '(a33)') '# omega, theta, phi, distribution'
   
          do iom = 1, this%nomega
            omega = iom * this%delomega
            isp = 0
   
            do ith = 0, this%nstepsthetar
              thetar = ith * M_PI / this%nstepsthetar
   
              do iph = 0, this%nstepsphir - 1
                isp = isp + 1
   
                phir = iph * M_TWO * M_PI / this%nstepsphir
                if(iph == 0) isp_save = isp
                write(iunittwo,'(5(1x,e18.10E3))') omega, thetar, phir, spctrout(isp, iom)
     
                ! just repeat the result for output
                if(this%nstepsphir > 1 .and. iph == (this%nstepsphir - 1)) then
                  write(iunittwo,'(5(1x,e18.10E3))') omega, thetar, M_TWO * M_PI, spctrout(isp_save, iom)
                end if
                 
                ! just repeat the result for output
                if(ith == 0 .or. ith == this%nstepsthetar) then
                  if(this%nstepsphir > 1) then
                    do iphi = 1, this%nstepsphir
                      phir = iphi * M_TWO * M_PI / this%nstepsphir
                      write(iunittwo,'(5(1x,e18.10E3))') omega, thetar, phir, spctrout(isp, iom)
                    end do
                  end if
                  exit
                end if
              end do
   
              if(this%nstepsphir > 1 .or. ith == this%nstepsthetar) write(iunittwo, '(1x)', advance='yes')
            end do

            ! write total and orbital spectra
            write(iunitone, '(2(1x,e18.10E3))', advance='no') omega, sum(sum(sum(spctrsum(:,:,:,iom),1),1),1) * sqrt(M_TWO * omega)
            do ik = 1, st%d%nik
              do ist = 1, st%nst
                do isdim = 1, st%d%dim
                  write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, iom) * sqrt(M_TWO * omega)
                end do
              end do
            end do
            write(iunitone, '(1x)', advance='yes')
   
          end do
        end select
        call io_close(iunittwo)
        call io_close(iunitone)

      end if
    end if

    ! -----------------------------------------------------------------
    ! DEBUG OUTPUT 
    ! -----------------------------------------------------------------
    if(.not. this%sphgrid .or. debug%info) then   ! too much output for spherical grid
      if(mpi_grp_is_root(mesh%mpi_grp)) then
        do ik = kptst, kptend
          do ist = stst, stend
            do isdim = 1, sdim
              itot = ist + (ik-1) * st%nst + (isdim-1) * st%nst*st%d%kpt%nglobal
              write(*,*) 'TEST', itot
              write(filenr, '(i10.10)') itot
   
              iunitone = io_open('td.general/'//'PES_spm.'//trim(filenr)//'.wavefunctions.out', &
                namespace, action='write', position='append')
     
              do ii = 1, save_iter - mod(iter, save_iter)
                jj = iter - save_iter + ii + mod(save_iter - mod(iter, save_iter), save_iter)
                write(iunitone, '(e17.10)', advance='yes') units_from_atomic(units_inp%time, jj * dt)
   
                do isp = 1, this%nspoints
                  vfu = units_from_atomic(sqrt(units_out%length**(-3)), this%wf(ist, isdim, ik, isp, ii-1))
                  write(iunitone, '(1x,e18.10E3,1x,e18.10E3)', advance='no') TOFLOAT(vfu), aimag(vfu)
                end do
                write(iunitone, '(1x)', advance='yes')
              end do
   
              call io_close(iunitone)
            end do
          end do
        end do

        if(this%onfly) then
          do ik = kptst, kptend
            do ist = stst, stend
              do isdim = 1, sdim
                itot = ist + (ik-1) * st%nst + (isdim-1) * st%nst*st%d%kpt%nglobal
                write(filenr, '(i10.10)') itot
   
                iunitone = io_open('td.general/'//'PES_spm.'//trim(filenr)//'.spectrum.out', &
                  namespace, action='write', position='rewind')
                write(iunitone, '(a48)') '# frequency, orbital spectrum at sampling points'
   
                do iom = 1, this%nomega 
                  omega = iom*this%delomega
                  write(iunitone, '(e17.10, 1x)', advance='no') omega
   
                  do isp = 1, this%nspoints
                    write(iunitone, '(e17.10, 1x)', advance='no') abs(this%wfft(ist, isdim, ik, isp, iom))**M_TWO
                  end do
   
                  write(iunitone, '(1x)', advance='yes')
                end do
   
                call io_close(iunitone)
              end do
            end do
          end do
        end if
      end if

      if(mpi_grp_is_root(mpi_world)) then
        if(this%recipe == M_PHASE) then
          do isp = 1, this%nspoints
            write(filenr, '(i10.10)') isp
            iunittwo = io_open('td.general/'//'PES_spm.'//trim(filenr)//'.phase.out', namespace, action='write', position='append')
      
            do ii = 1, save_iter - mod(iter, save_iter)
              jj = iter - save_iter + ii + mod(save_iter - mod(iter, save_iter), save_iter)
      
              write(iunittwo, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
              write(iunittwo, '(1x,e18.10E3,1x,e18.10E3)', advance='no') this%domega(ii-1), this%dq(isp, ii-1)
              write(iunittwo, '(1x)', advance='yes')
            end do
   
            call io_close(iunittwo)
          end do
        end if

        if(this%onfly) then
          iunitone = io_open('td.general/'//'PES_spm.total.out', namespace, action='write', position='rewind')
          write(iunitone, '(a46)') '# frequency, total spectrum at sampling points'
          do iom = 1, this%nomega
            omega = iom*this%delomega
  
            write(iunitone, '(e17.10, 1x)', advance='no') omega 
            do isp = 1, this%nspoints
              write(iunitone, '(e17.10, 1x)', advance='no') spctrout(isp, iom)
            end do
   
            write(iunitone, '(1x)', advance='yes')
          end do
          call io_close(iunitone)
        end if
      end if
    end if

    SAFE_DEALLOCATE_A(spctrout)
    SAFE_DEALLOCATE_A(spctrsum)

    POP_SUB(pes_spm_output)
  end subroutine pes_spm_output

  ! ---------------------------------------------------------
  subroutine pes_spm_init_write(this, mesh, st, namespace)
    type(PES_spm_t),     intent(in) :: this
    type(mesh_t),        intent(in) :: mesh
    type(states_elec_t), intent(in) :: st
    type(namespace_t),   intent(in) :: namespace

    integer           :: ist, ik, isdim
    integer           :: isp
    FLOAT             :: xx(MAX_DIM)
    character(len=80) :: filenr
    integer           :: iunit

    PUSH_SUB(pes_spm_init_write)

    xx = M_ZERO
    if(mpi_grp_is_root(mpi_world)) then
      if(.not. this%sphgrid .or. debug%info) then   ! too much output for spherical grid
        do isp = 1, this%nspoints
          write(filenr, '(i10.10)') isp
   
          iunit = io_open('td.general/'//'PES_spm.'//trim(filenr)//'.wavefunctions.out', namespace, action='write')
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
            iunit = io_open('td.general/'//'PES_spm.'//trim(filenr)//'.phase.out', namespace, action='write')
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
    type(restart_t),     intent(in)  :: restart    
    type(pes_spm_t),     intent(in)  :: this
    type(states_elec_t), intent(in)  :: st
    integer,             intent(out) :: ierr
    
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
    type(restart_t),     intent(in)    :: restart    
    type(pes_spm_t),     intent(inout) :: this
    type(states_elec_t), intent(inout) :: st
    integer,             intent(out)   :: ierr
    
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
    type(hamiltonian_elec_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: ii
    
    integer :: mdim
    integer :: isp
    integer :: il, iprev
    FLOAT   :: vp(1:MAX_DIM)
    FLOAT   :: rdota

    PUSH_SUB(pes_spm_calc_rcphase)

    mdim = mesh%sb%dim

    vp = M_ZERO
    do il = 1, hm%ext_lasers%no_lasers
      call laser_field(hm%ext_lasers%lasers(il), vp(1:mdim), iter*dt)
    end do
    vp(1:mdim) = -vp(1:mdim)

    iprev = ii - 1
    if(ii == 0) iprev = this%save_iter - 1

    do isp = 1, this%nspoints
      rdota = dot_product(this%rcoords_nrm(1:mdim, isp), vp(1:mdim))
      this%dq(isp, ii) = this%dq(isp, iprev) + rdota * dt / P_C
    end do

    this%domega(ii) = this%domega(iprev) &
      + dot_product(vp(1:mdim), vp(1:mdim)) / (M_TWO * P_C**M_TWO) * dt

    POP_SUB(pes_spm_calc_rcphase)
  end subroutine pes_spm_calc_rcphase

end module pes_spm_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
