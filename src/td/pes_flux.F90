!! Copyright (C) 2015 P. Wopperer
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

module pes_flux_m
  use global_m
  use messages_m
  use mesh_m
  use profiling_m
  use mpi_m
  use comm_m
  use loct_m
  use loct_math_m
  use math_m
  use parser_m
  use par_vec_m
  use states_m
  use grid_m
  use derivatives_m
  use hamiltonian_m
  use lasers_m
  use io_function_m
  use restart_m
  use io_binary_m
  use io_m
  use simul_box_m
  use varinfo_m
  use mesh_interpolation_m

  implicit none

  private

  public ::                    &
    pes_flux_t,                &
    pes_flux_init,             &
    pes_flux_end,              &
    pes_flux_save,             &
    pes_flux_output,           &
    pes_flux_load,             &
    pes_flux_dump

  type pes_flux_t
    integer          :: nkpnts                         !< total number of k-points
    integer          :: nkpnts_start, nkpnts_end       !< start/end of index for k-points on the current node
    integer          :: nk
    integer          :: nk_start, nk_end
    integer          :: nstepsthetak, nstepsphik       !< parameters for k-mesh
    integer          :: nstepsomegak
    integer          :: nstepsomegak_start, nstepsomegak_end
    FLOAT            :: dk, dthetak                    !< parameters for k-mesh
    FLOAT, pointer   :: kcoords_cub(:,:)               !< coordinates of k-points
    FLOAT, pointer   :: kcoords_sph(:,:,:)

    integer          :: shape                          !< shape of the surface (= cubic/spherical)
    integer          :: nsrfcpnts                      !< total number of surface points
    integer          :: nsrfcpnts_start, nsrfcpnts_end !< number of surface points on node
    FLOAT, pointer   :: srfcnrml(:,:)                  !< vectors normal to the surface (includes surface element)
    FLOAT, pointer   :: rcoords(:,:)                   !< coordinates of the surface points
    integer, pointer :: srfcpnt(:)                     !< for cubic surface: returns index of the surface points
    integer, pointer :: rankmin(:)                     !< for cubic surface: returns node which has the surface point
    integer          :: lmax                           !< for spherical surface
    CMPLX, pointer   :: ylm_r(:,:,:)                   !< for spherical surface
    CMPLX, pointer   :: ylm_k(:,:,:)                   !< for spherical surface
    FLOAT, pointer   :: j_l(:,:)                       !< for spherical surface
    FLOAT            :: radius

    integer          :: tdsteps                        !<
    integer          :: tdstepsinterval                !<

    CMPLX, pointer   :: wf(:,:,:,:,:)                  !< wavefunction
    CMPLX, pointer   :: gwf(:,:,:,:,:,:)               !< gradient of wavefunction
    FLOAT, pointer   :: veca(:,:)                      !< vector potential
    CMPLX, pointer   :: conjgphase_prev_cub(:)         !< Volkov phase for all k-points from previous time step
    CMPLX, pointer   :: conjgphase_prev_sph(:,:)
    CMPLX, pointer   :: conjgplanewf_cub(:,:)          !< plane wave factor
    CMPLX, pointer   :: spctramp_cub(:,:,:,:)          !< spectral amplitude
    CMPLX, pointer   :: spctramp_sph(:,:,:,:,:)

    logical          :: usememory                      !< whether conjgplanewf should be kept in memory
    type(mesh_interpolation_t) :: interp

  end type pes_flux_t

  integer, parameter ::   &
    M_CUBIC      = 1,     &
    M_SPHERICAL  = 2

contains


  ! ---------------------------------------------------------
  subroutine pes_flux_init(this, mesh, st, hm)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(in)    :: hm

    type(block_t)      :: blk
    FLOAT              :: border(MAX_DIM)       ! distance of surface from border
    FLOAT              :: offset(MAX_DIM)       ! offset for border
    integer            :: sdim, mdim
    integer            :: imdim
    integer            :: isp, ikp
    integer            :: il
    integer            :: ikk, ith, iph, iomk
    FLOAT              :: kmax, kact, thetak, phik

    FLOAT, allocatable :: k_dot_aux(:)
    integer            :: nstepsphir, nstepsthetar
#if defined(HAVE_MPI)
    integer            :: mpisize, mpirank
#endif
    integer            :: ll, mm

    PUSH_SUB(pes_flux_init)

    sdim   = st%d%dim
    mdim   = mesh%sb%dim

#if defined(HAVE_MPI)
    call mpi_comm_size(mpi_comm_world, mpisize, mpi_err)
    call mpi_comm_rank(mpi_comm_world, mpirank, mpi_err)
#endif

    call messages_experimental("PhotoElectronSpectrum with t-surff")

    do il = 1, hm%ep%no_lasers
      if(laser_kind(hm%ep%lasers(il)) /= E_FIELD_VECTOR_POTENTIAL) then
        message(1) = 't-surff only works in velocity gauge.'
        call messages_fatal(1)
      end if
    end do

    message(1) = 'Info: Calculating PES using t-surff technique.'
    call messages_info(1)

    ! -----------------------------------------------------------------
    ! Setting up r-mesh (the surface)
    ! -----------------------------------------------------------------
    !%Variable PES_Flux_Shape
    !%Type integer
    !%Default sph
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The shape of the surface.
    !%Option cub 1
    !% Uses a parallelepiped as surface. All surface points are grid points.
    !% Choose the location of the surface with variable PES_Flux_Lsize.
    !%Option sph 2
    !% Constructs a sphere with radius PES_Flux_Radius. Points on the sphere 
    !% are interpolated by trilinear interpolation.
    !%End
    call parse_variable('PES_Flux_Shape', M_SPHERICAL, this%shape)
    if(.not.varinfo_valid_option('PES_Flux_Shape', this%shape, is_flag = .true.)) &
      call messages_input_error('PES_Flux_Shape')
    call messages_print_var_option(stdout, 'PES_Flux_Shape', this%shape)

    !%Variable PES_Flux_Offset
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Shifts the surface for PES_Flux_Shape = cub. The syntax is:
    !%
    !% <tt>%PES_Flux_Offset
    !% <br>&nbsp;&nbsp;xshift | yshift | zshift
    !% <br>%
    !% </tt>
    !%End
    offset = M_ZERO
    if(parse_block('PES_Flux_Offset', blk) == 0) then
      call parse_block_float(blk, 0, 0, offset(1))
      call parse_block_float(blk, 0, 1, offset(2))
      call parse_block_float(blk, 0, 2, offset(3))
    end if

    !%Variable PES_Flux_Lmax
    !%Type integer
    !%Default 1
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Maximum order of the spherical harmonic to be integrated on an equidistant spherical 
    !% grid (to be changed to Gauss-Legendre quadrature).
    !%End
    if(this%shape == M_SPHERICAL) then
      call parse_variable('PES_Flux_Lmax', 1, this%lmax)
    end if

    !%Variable PES_Flux_Radius
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The radius of the sphere, if PES_Flux_Shape == sph.
    !%End

    !%Variable PES_Flux_Lsize
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Sets the location of the surface for PES_Flux_Shape = cub. The syntax is:
    !%
    !% <tt>%PES_Flux_Lsize
    !% <br>&nbsp;&nbsp;xsize | ysize | zsize
    !% <br>%
    !% </tt>
    !%
    !% where xsize, ysize, and zsize are with respect to the origin. The surface can 
    !% be shifted with PES_Flux_Offset.
    !%End
    if(this%shape == M_CUBIC) then
      if(parse_block('PES_Flux_Lsize', blk) == 0) then
        call parse_block_float(blk, 0, 0, border(1))
        call parse_block_float(blk, 0, 1, border(2))
        call parse_block_float(blk, 0, 2, border(3))
        border(1:mdim) = int(border(1:mdim)/mesh%spacing(1:mdim))*mesh%spacing(1:mdim)
      else
        select case(mesh%sb%box_shape)
        case(PARALLELEPIPED)
          border(1:mdim) = mesh%sb%lsize(1:mdim)
        case(SPHERE)
          border(1:mdim) = mesh%sb%rsize/sqrt(M_TWO)
        case default
          message(1) = 'PES_Flux_Lsize not specified. No default values available for this box shape.'
          message(2) = 'Specify the location of the parallelepiped with block PES_Flux_Lsize.'
          call messages_fatal(2)
        end select
        message(1) = 'PES_Flux_Lsize not specified. Using default values.'
        call messages_info(1)
        call messages_print_var_value(stdout, 'PES_Flux_Lsize', border(1:mdim))
      end if

    else

      if(parse_is_defined('PES_Flux_Radius')) then
        call parse_variable('PES_Flux_Radius', M_ZERO, this%radius)
        if(this%radius <= M_ZERO) call messages_input_error('PES_Flux_Radius')
        call messages_print_var_value(stdout, 'PES_Flux_Radius', this%radius)
      else
        select case(mesh%sb%box_shape)
        case(PARALLELEPIPED)
          this%radius = minval(mesh%sb%lsize(1:mdim))
        case(SPHERE)
          this%radius = mesh%sb%rsize
        case default
          message(1) = 'PES_Flux_Radius not specified. No default values available for this box shape.'
          message(2) = 'Specify the radius of the sphere with variable PES_Flux_Radius.'
          call messages_fatal(2)
        end select
        message(1) = 'PES_Flux_Radius not specified. Using default values.'
        call messages_info(1)
        call messages_print_var_value(stdout, 'PES_Flux_Radius', this%radius)
      end if

    end if

    !%Variable PES_Flux_StepsThetaR
    !%Type integer
    !%Default 45
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in theta (0 <= theta <= pi) for the spherical surface.
    !%End
    if(this%shape == M_SPHERICAL) then
      call parse_variable('PES_Flux_StepsThetaR', 45, nstepsthetar)
      if(nstepsthetar < 0) call messages_input_error('PES_Flux_StepsThetaR')
    end if

    !%Variable PES_Flux_StepsPhiR
    !%Type integer
    !%Default 90
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in phi (0 <= phi <= 2 pi) for the spherical surface.
    !%End
    if(this%shape == M_SPHERICAL) then
      call parse_variable('PES_Flux_StepsPhiR', 90, nstepsphir)
      if(nstepsphir < 0) call messages_input_error('PES_Flux_StepsPhiR')
    end if

    ! get the surface points
    if(this%shape == M_CUBIC) then
      call pes_flux_getcube(this, mesh, border, offset)
    else
      if(mdim /= 3) then
        message(1) = 'Spherical grid works only in 3d.'
        call messages_fatal(1)
      end if
      call mesh_interpolation_init(this%interp, mesh)
      call pes_flux_getsphere(this, mesh, nstepsthetar, nstepsphir, offset)
    end if

    ! distribute the surface points on nodes,
    ! since mesh domains may have different numbers of surface points.
    call pes_flux_distribute(1, this%nsrfcpnts, this%nsrfcpnts_start, this%nsrfcpnts_end)
#if defined(HAVE_MPI)
    call MPI_Barrier(mpi_world%comm, mpi_err)
    write(*,*) &
      'Number of surface points on node ', mpirank, ' : ', this%nsrfcpnts_end - this%nsrfcpnts_start + 1
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
    if(mpi_grp_is_root(mpi_world)) then
      write(*,*) 'Info: Number of surface points (total):', this%nsrfcpnts
      do isp = 1, this%nsrfcpnts
        write(223,*) isp, this%rcoords(:, isp)
      end do
      flush(223)
    end if

    ! get the values of the spherical harmonics for the surface points for M_SPHERICAL
    if(this%shape == M_SPHERICAL) then
      SAFE_ALLOCATE(this%ylm_r(0:this%lmax, -this%lmax:this%lmax, 1:this%nsrfcpnts))
      this%ylm_r = M_z0

      do isp = 1, this%nsrfcpnts
        do ll = 0, this%lmax
          do mm = -ll, ll
            call ylmr(this%rcoords(1, isp), this%rcoords(2, isp), this%rcoords(3, isp), ll, mm, this%ylm_r(ll, mm, isp))
            this%ylm_r(ll, mm, isp) = conjg(this%ylm_r(ll, mm, isp))
          end do
        end do
      end do
    end if

    ! -----------------------------------------------------------------
    ! Setting up k-mesh
    ! 1D = 2 points, 2D = polar coordinates, 3D = spherical coordinates
    ! -----------------------------------------------------------------
    !%Variable PES_Flux_Kmax
    !%Type float
    !%Default 1.0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The maximum value of |k|.
    !%End
    call parse_variable('PES_Flux_Kmax', M_ONE, kmax)
    call messages_print_var_value(stdout, "PES_Flux_Kmax", kmax)
    if(kmax <= M_ZERO) call messages_input_error('PES_Flux_Kmax')

    !%Variable PES_Flux_DeltaK
    !%Type float
    !%Default 0.002
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Spacing of the k-mesh in |k| (equidistant).
    !%End
    call parse_variable('PES_Flux_DeltaK', CNST(0.002), this%dk)
    if(this%dk <= M_ZERO) call messages_input_error('PES_Flux_DeltaK')
    call messages_print_var_value(stdout, "PES_Flux_DeltaK", this%dk)

    !%Variable PES_Flux_StepsThetaK
    !%Type integer
    !%Default 45
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in theta (0 <= theta <= pi) for the spherical grid in k.
    !%End
    call parse_variable('PES_Flux_StepsThetaK', 45, this%nstepsthetak)
    if(this%nstepsthetak < 0) call messages_input_error('PES_Flux_StepsThetaK')

    !%Variable PES_Flux_StepsPhiK
    !%Type integer
    !%Default 90
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in phi (0 <= phi <= 2 pi) for the spherical grid in k.
    !%End
    call parse_variable('PES_Flux_StepsPhiK', 90, this%nstepsphik)
    if(this%nstepsphik < 0) call messages_input_error('PES_Flux_StepsPhiK')

    if(this%nstepsphik == 0) this%nstepsphik = 1

    select case(mdim)
    case(1)
      this%nstepsthetak = 0
      this%nstepsphik   = 2
      this%nstepsomegak = this%nstepsphik

    case(2)
      this%nstepsthetak = 0
      this%nstepsomegak = this%nstepsphik

    case(3)
      if(this%nstepsthetak <= 1) then
        this%nstepsphik   = 1
        this%nstepsthetak = 1
      end if
      this%dthetak = M_PI / this%nstepsthetak
      this%nstepsomegak = this%nstepsphik * (this%nstepsthetak - 1) + 2

    end select

    if(mdim == 3) call messages_print_var_value(stdout, "PES_Flux_StepsThetaK", this%nstepsthetak)
    call messages_print_var_value(stdout, "PES_Flux_StepsPhiK", this%nstepsphik)

    this%nk     = nint(kmax/this%dk)
    this%nkpnts = this%nstepsomegak * this%nk

    if(this%shape == M_SPHERICAL) then
      ! we split the k-mesh in radial & angular part
      call pes_flux_distribute(1, this%nk, this%nk_start, this%nk_end)
      call pes_flux_distribute(1, this%nstepsomegak, this%nstepsomegak_start, this%nstepsomegak_end)

      SAFE_ALLOCATE(this%ylm_k(0:this%lmax, -this%lmax:this%lmax, 1:this%nstepsomegak))
      SAFE_ALLOCATE(this%j_l(0:this%lmax, 1:this%nk))
      SAFE_ALLOCATE(this%kcoords_sph(1:3, 1:this%nk, 1:this%nstepsomegak))

      this%ylm_k = M_z0
      this%kcoords_sph = M_ZERO
      this%j_l = M_ZERO

      ! spherical harmonics & kcoords_sph
      iomk = 0
      do ith = 0, this%nstepsthetak
        thetak = ith * this%dthetak
        do iph = 0, this%nstepsphik - 1
          phik = iph * M_TWO * M_PI / this%nstepsphik
          iomk = iomk + 1
          do ll = 0, this%lmax
            do mm = -ll, ll
              call ylmr(cos(phik) * sin(thetak), sin(phik) * sin(thetak), cos(thetak), ll, mm, this%ylm_k(ll, mm, iomk))
            end do
          end do
          this%kcoords_sph(1, this%nk_start:this%nk_end, iomk) = cos(phik) * sin(thetak)
          this%kcoords_sph(2, this%nk_start:this%nk_end, iomk) = sin(phik) * sin(thetak)
          this%kcoords_sph(3, this%nk_start:this%nk_end, iomk) = cos(thetak)
          if(ith == 0 .or. ith == this%nstepsthetak) exit
        end do
      end do

      ! Bessel functions & kcoords_sph
      do ikk = this%nk_start, this%nk_end
        kact = ikk * this%dk
        do ll = 0, this%lmax
          this%j_l(ll, ikk) = loct_sph_bessel(ll, kact * this%radius) * &
                              M_TWO * M_PI / (M_TWO * M_PI)**M_THREE/M_TWO
        end do
        this%kcoords_sph(:, ikk, :) = kact * this%kcoords_sph(:, ikk, :)
      end do
#if defined(HAVE_MPI)
      call comm_allreduce(mpi_world%comm, this%kcoords_sph)
      call comm_allreduce(mpi_world%comm, this%j_l)
#endif
    else
      ! we do not split the k-mesh
      call pes_flux_distribute(1, this%nkpnts, this%nkpnts_start, this%nkpnts_end)
#if defined(HAVE_MPI)
      call MPI_Barrier(mpi_world%comm, mpi_err)
      write(*,*) &
        'Number of k-points on node ', mpirank, ' : ', this%nkpnts_end - this%nkpnts_start + 1
      call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
      if(mpi_grp_is_root(mpi_world)) then
        write(*,*) 'Info: Number of k-points (total):', this%nkpnts
      end if

      SAFE_ALLOCATE(this%kcoords_cub(1:mdim, 1:this%nkpnts))
      this%kcoords_cub = M_ZERO

      thetak = M_PI / M_TWO
      ikp = 0
      do ikk = 1, this%nk
        do ith = 0, this%nstepsthetak
          if(mdim == 3) thetak = ith * this%dthetak
          do iph = 0, this%nstepsphik - 1
            ikp = ikp + 1
            phik = iph * M_TWO * M_PI / this%nstepsphik
            kact = ikk * this%dk
                          this%kcoords_cub(1, ikp) = kact * cos(phik) * sin(thetak)
            if(mdim >= 2) this%kcoords_cub(2, ikp) = kact * sin(phik) * sin(thetak)
            if(mdim == 3) this%kcoords_cub(3, ikp) = kact * cos(thetak)
            if(thetak == M_ZERO .or. thetak == M_PI) exit
          end do
        end do
      end do
    end if

    ! -----------------------------------------------------------------
    ! Options for time integration 
    ! -----------------------------------------------------------------
    !%Variable PES_Flux_TDStepsInterval
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End
    call parse_variable('PES_Flux_TDStepsInterval', 1, this%tdstepsinterval)

    !%Variable PES_Flux_TDSteps
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End
    if(this%shape == M_CUBIC) then
      call parse_variable('PES_Flux_TDSteps', 1, this%tdsteps)
    else
#if defined(HAVE_MPI)
      this%tdsteps = mpisize
#else
      this%tdsteps = M_ONE
#endif
    end if

    ! -----------------------------------------------------------------
    ! Other stuff
    ! -----------------------------------------------------------------
    !%Variable PES_Flux_UseMemory
    !%Type logical
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End
    if(this%shape == M_CUBIC) &
      call parse_variable('PES_Flux_UseMemory', .true., this%usememory)

    SAFE_ALLOCATE(this%wf(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nsrfcpnts, 1:this%tdsteps))
    this%wf = M_z0

    SAFE_ALLOCATE(this%gwf(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nsrfcpnts, 1:this%tdsteps, 1:mdim))
    this%gwf = M_z0

    SAFE_ALLOCATE(this%veca(1:mdim, 1:this%tdsteps))
    this%veca = M_ZERO

    if(this%shape == M_SPHERICAL) then
      SAFE_ALLOCATE(this%spctramp_sph(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nk, 1:this%nstepsomegak))
      this%spctramp_sph = M_z0

      SAFE_ALLOCATE(this%conjgphase_prev_sph(1:this%nk, 1:this%nstepsomegak))
      this%conjgphase_prev_sph = M_z1

    else
      SAFE_ALLOCATE(this%spctramp_cub(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nkpnts))
      this%spctramp_cub = M_z0

      SAFE_ALLOCATE(this%conjgphase_prev_cub(1:this%nkpnts))
      this%conjgphase_prev_cub = M_z1

      if(this%usememory) then
        SAFE_ALLOCATE(k_dot_aux(1:this%nkpnts))
        SAFE_ALLOCATE(this%conjgplanewf_cub(1:this%nkpnts, 1:this%nsrfcpnts))
        this%conjgplanewf_cub = M_z0

        do isp = this%nsrfcpnts_start, this%nsrfcpnts_end
          k_dot_aux = M_ZERO
          do imdim = 1, mdim
            k_dot_aux(:) = k_dot_aux(:) + this%kcoords_cub(imdim, :) * this%rcoords(imdim, isp)
          end do
          this%conjgplanewf_cub(:, isp) = exp(-M_zI * k_dot_aux(:)) / (M_TWO * M_PI)**(mdim/M_TWO)
        end do
        SAFE_DEALLOCATE_A(k_dot_aux)
#if defined(HAVE_MPI)
        call comm_allreduce(mpi_world%comm, this%conjgplanewf_cub)
#endif
      end if
    end if

    POP_SUB(pes_flux_init)
  end subroutine pes_flux_init

  ! ---------------------------------------------------------
  subroutine pes_flux_end(this)
    type(pes_flux_t), intent(inout) :: this

    PUSH_SUB(pes_flux_end)

    if(this%shape == M_SPHERICAL) then
      SAFE_DEALLOCATE_P(this%kcoords_sph)
      SAFE_DEALLOCATE_P(this%ylm_k)
      SAFE_DEALLOCATE_P(this%j_l)
      SAFE_DEALLOCATE_P(this%ylm_r)
      SAFE_DEALLOCATE_P(this%conjgphase_prev_sph)
      SAFE_DEALLOCATE_P(this%spctramp_sph)
    else
      SAFE_DEALLOCATE_P(this%kcoords_cub)
      SAFE_DEALLOCATE_P(this%conjgphase_prev_cub)
      SAFE_DEALLOCATE_P(this%spctramp_cub)
      if(this%usememory) then
        SAFE_DEALLOCATE_P(this%conjgplanewf_cub)
      end if

      SAFE_DEALLOCATE_P(this%srfcpnt)
      SAFE_DEALLOCATE_P(this%rankmin)
    end if

    SAFE_DEALLOCATE_P(this%srfcnrml)
    SAFE_DEALLOCATE_P(this%rcoords)

    SAFE_DEALLOCATE_P(this%wf)
    SAFE_DEALLOCATE_P(this%gwf)
    SAFE_DEALLOCATE_P(this%veca)

    POP_SUB(pes_flux_end)
  end subroutine pes_flux_end

  ! ---------------------------------------------------------
  subroutine pes_flux_save(this, mesh, st, gr, hm, iter, maxiter, dt)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    integer,             intent(in)    :: maxiter
    FLOAT,               intent(in)    :: dt

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim, imdim
    integer            :: isp, itstep
    integer            :: il, ip_local
    CMPLX, allocatable :: gpsi(:,:), psi(:)
    CMPLX, allocatable :: interp_values(:)
    logical            :: contains_isp

    PUSH_SUB(pes_flux_save)

    if(iter > 0 .and. mod(iter, this%tdstepsinterval) == 0) then

      stst   = st%st_start
      stend  = st%st_end
      kptst  = st%d%kpt%start
      kptend = st%d%kpt%end
      sdim   = st%d%dim
      mdim   = mesh%sb%dim

      SAFE_ALLOCATE(psi(1:mesh%np_part))
      SAFE_ALLOCATE(gpsi(1:mesh%np_part, 1:mdim))

      if(this%shape == M_SPHERICAL) then
        SAFE_ALLOCATE(interp_values(1:this%nsrfcpnts))
      end if

      if(iter > 0 .and. mod(int(iter/this%tdstepsinterval), this%tdsteps) == 0) then
        itstep = this%tdsteps
      else
        itstep = mod(int(iter/this%tdstepsinterval), this%tdsteps) 
      end if

      ! clean up fields when a new cycle begins
      if(itstep == 1) then
         this%wf  = M_z0
        this%gwf  = M_z0
        this%veca = M_ZERO
      end if

      ! get and save current laser field
      do il = 1, hm%ep%no_lasers
        call laser_field(hm%ep%lasers(il), this%veca(1:mdim, itstep), iter*dt)
      end do
      this%veca(:, itstep) = - this%veca(:, itstep)

      ! save wavefunctions & gradients
      do ik = kptst, kptend
        do ist = stst, stend
          do isdim = 1, sdim
            call states_get_state(st, mesh, isdim, ist, ik, psi)
            call zderivatives_grad(gr%der, psi, gpsi, .true.)

            if(this%shape == M_SPHERICAL) then
              call mesh_interpolation_evaluate(this%interp, this%nsrfcpnts, psi(1:mesh%np_part), &
                this%rcoords(1:mdim, 1:this%nsrfcpnts), interp_values(1:this%nsrfcpnts))
              this%wf(ist, isdim, ik, :, itstep) = st%occ(ist, ik) * interp_values(:)
              do imdim = 1, mdim
                call mesh_interpolation_evaluate(this%interp, this%nsrfcpnts, gpsi(1:mesh%np_part, imdim), &
                  this%rcoords(1:mdim, 1:this%nsrfcpnts), interp_values(1:this%nsrfcpnts))
                this%gwf(ist, isdim, ik, :, itstep, imdim) = st%occ(ist, ik) * interp_values(:)
              end do

            else

              contains_isp = .true.
              do isp = 1, this%nsrfcpnts
                if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
                  if(mesh%mpi_grp%rank == this%rankmin(isp)) then
                    contains_isp = .true.
                  else
                    contains_isp = .false.
                  end if
#endif
                end if

                if(contains_isp) then
                  ip_local = this%srfcpnt(isp)
                  this%wf(ist, isdim, ik, isp, itstep) = st%occ(ist, ik) * psi(ip_local)
                  this%gwf(ist, isdim, ik, isp, itstep, 1:mdim) = &
                    st%occ(ist, ik) * gpsi(ip_local, 1:mdim) ! * mesh%spacing(1:mdim) ?
                end if
              end do
            end if
          end do
        end do
      end do

      if(this%shape == M_SPHERICAL) then 

        if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
          do isp = 1, this%nsrfcpnts
            call comm_allreduce(st%st_kpt_mpi_grp%comm, this%wf(:,:,:, isp, itstep))
            do imdim = 1, mdim
              call comm_allreduce(st%st_kpt_mpi_grp%comm, this%gwf(:,:,:, isp, itstep, imdim))
            end do
          end do
#endif
        endif

      else
#if defined(HAVE_MPI)
        call comm_allreduce(mpi_world%comm, this%wf(:,:,:,:, itstep))
        do imdim = 1, mdim
          call comm_allreduce(mpi_world%comm, this%gwf(:,:,:,:, itstep, imdim))
        end do
#endif
      end if

      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(gpsi)

      if(iter > 0) then
        if(mod(iter, this%tdsteps * this%tdstepsinterval) == 0 .or. iter == maxiter) then
          if(this%shape == M_CUBIC) then
            call pes_flux_integrate_cub(this, mesh, st, dt)
          else
            call pes_flux_integrate_sph(this, mesh, st, dt)
          end if
        end if
      end if

    end if

    POP_SUB(pes_flux_save)
  end subroutine pes_flux_save

  ! ---------------------------------------------------------
  subroutine pes_flux_integrate_cub(this, mesh, st, dt)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st
    FLOAT,               intent(in)    :: dt

    integer            :: sdim, mdim
    integer            :: ist, ik, isdim, imdim
    integer            :: isp, ikp, itstep
    integer            :: idir
    CMPLX, allocatable :: Jk_cub(:,:,:,:), spctramp_cub(:,:,:,:)
    CMPLX, allocatable :: conjgplanewf_cub(:)
    CMPLX, allocatable :: conjgphase_cub(:,:)
    FLOAT, allocatable :: k_dot_aux(:)
    integer            :: ikp_start, ikp_end, isp_start, isp_end
    FLOAT              :: vec

    PUSH_SUB(pes_flux_integrate_cub)

    ! this routine is parallelized over surface points since the number of 
    ! states is in most cases less than the number of surface points

    sdim      = st%d%dim
    mdim      = mesh%sb%dim
    ikp_start = this%nkpnts_start
    ikp_end   = this%nkpnts_end
    isp_start = this%nsrfcpnts_start
    isp_end   = this%nsrfcpnts_end

    SAFE_ALLOCATE(k_dot_aux(1:this%nkpnts))

    SAFE_ALLOCATE(Jk_cub(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nkpnts))
    SAFE_ALLOCATE(spctramp_cub(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nkpnts))
    spctramp_cub = M_z0

    if(.not. this%usememory) then
      SAFE_ALLOCATE(conjgplanewf_cub(1:this%nkpnts))
    end if

    SAFE_ALLOCATE(conjgphase_cub(1:this%nkpnts, 0:this%tdsteps))
    conjgphase_cub = M_z0

    ! calculate Volkov phase using the previous time step
    conjgphase_cub(:, 0) = this%conjgphase_prev_cub(:)
    do ikp = ikp_start, ikp_end
      do itstep = 1, this%tdsteps
        vec = sum((this%kcoords_cub(1:mdim, ikp) - this%veca(1:mdim, itstep) / P_c)**2)
        conjgphase_cub(ikp, itstep) = conjgphase_cub(ikp, itstep - 1) * exp(M_zI * vec * dt * this%tdstepsinterval / M_TWO)
      end do
    end do
#if defined(HAVE_MPI)
    call comm_allreduce(mpi_world%comm, conjgphase_cub)
#endif
    this%conjgphase_prev_cub(:) = conjgphase_cub(:, this%tdsteps)

    ! integrate over time & surface (on node)
    do isp = isp_start, isp_end 
      do idir = 1, mdim
        ! calculate flux only along the surface normal
        if(this%srfcnrml(idir, isp) == M_ZERO) cycle

        Jk_cub = M_z0

        if(.not. this%usememory) then
          k_dot_aux(:) = M_ZERO
          do imdim = 1, mdim
            k_dot_aux(:) = k_dot_aux(:) + this%kcoords_cub(imdim, :) * this%rcoords(imdim, isp)
          end do
          conjgplanewf_cub(:) = exp(-M_zI * k_dot_aux(:)) / (M_TWO * M_PI)**(mdim/M_TWO)
        end if

        do ik = 1, st%d%nik 
          do ist = 1, st%nst
            do isdim = 1, sdim

              do itstep = 1, this%tdsteps
                Jk_cub(ist, isdim, ik, 1:this%nkpnts) = &
                  Jk_cub(ist, isdim, ik, 1:this%nkpnts) + conjgphase_cub(1:this%nkpnts, itstep) * &
                  (this%wf(ist, isdim, ik, isp, itstep) * &
                   (M_TWO * this%veca(idir, itstep) / P_c - this%kcoords_cub(idir, 1:this%nkpnts)) + &
                   this%gwf(ist, isdim, ik, isp, itstep, idir) * M_zI)
              end do

              if(this%usememory) then
                Jk_cub(ist, isdim, ik, 1:this%nkpnts) = &
                  Jk_cub(ist, isdim, ik, 1:this%nkpnts) * this%conjgplanewf_cub(1:this%nkpnts, isp)
              else
                Jk_cub(ist, isdim, ik, 1:this%nkpnts) = &
                  Jk_cub(ist, isdim, ik, 1:this%nkpnts) * conjgplanewf_cub(1:this%nkpnts)
              end if

            end do
          end do
        end do
        spctramp_cub(:,:,:,:) = spctramp_cub(:,:,:,:) + Jk_cub(:,:,:,:) * this%srfcnrml(idir, isp)
      end do
    end do

#if defined(HAVE_MPI)
    call comm_allreduce(mpi_world%comm, spctramp_cub)
#endif

    this%spctramp_cub = this%spctramp_cub + spctramp_cub

    SAFE_DEALLOCATE_A(k_dot_aux)
    SAFE_DEALLOCATE_A(Jk_cub)
    SAFE_DEALLOCATE_A(spctramp_cub)
    SAFE_DEALLOCATE_A(conjgplanewf_cub)
    SAFE_DEALLOCATE_A(conjgphase_cub)

    POP_SUB(pes_flux_integrate_cub)
  end subroutine pes_flux_integrate_cub

  ! ---------------------------------------------------------
  subroutine pes_flux_integrate_sph(this, mesh, st, dt)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st
    FLOAT,               intent(in)    :: dt

    integer            :: sdim
    integer            :: ist, ik, isdim, imdim
    integer            :: itstep, lmax
    integer            :: ll, mm 
    CMPLX, allocatable :: sigma1(:,:), sigma2(:,:), sigma3(:)
    CMPLX, allocatable :: integ1_s(:,:,:,:,:,:,:), integ2_s(:,:,:,:,:,:)
    CMPLX, allocatable :: integ11_t(:,:), integ12_t(:,:,:)
    CMPLX, allocatable :: integ21_t(:), integ22_t(:,:)
    CMPLX, allocatable :: spctramp_sph(:,:,:,:,:)
    integer            :: ikk, ikk_start, ikk_end
    integer            :: iomk
    CMPLX, allocatable :: conjgphase_sph(:,:,:)
    FLOAT              :: vec
    integer            :: tdsteps_start, tdsteps_end

    PUSH_SUB(pes_flux_integrate_sph)

    ! this routine is parallelized over time steps

    call pes_flux_distribute(1, this%tdsteps, tdsteps_start, tdsteps_end)

    sdim = st%d%dim

    lmax       = this%lmax
    ikk_start  = this%nk_start
    ikk_end    = this%nk_end

    SAFE_ALLOCATE(conjgphase_sph(1:this%nk, 1:this%nstepsomegak, 0:this%tdsteps))
    conjgphase_sph = M_z0

    ! calculate Volkov phase using the previous time step
    conjgphase_sph(:,:, 0) = this%conjgphase_prev_sph(:,:)
    do ikk = ikk_start, ikk_end
      do iomk = 1, this%nstepsomegak
        do itstep = 1, this%tdsteps
          vec = sum((this%kcoords_sph(1:3, ikk, iomk) - this%veca(1:3, itstep) / P_c)**2)
          conjgphase_sph(ikk, iomk, itstep) = conjgphase_sph(ikk, iomk, itstep - 1) * &
            exp(M_zI * vec * dt * this%tdstepsinterval / M_TWO)
        end do
      end do
    end do
#if defined(HAVE_MPI)
    call comm_allreduce(mpi_world%comm, conjgphase_sph)
#endif
    this%conjgphase_prev_sph(:,:) = conjgphase_sph(:,:, this%tdsteps)

    ! surface integral S_lm (for time range)
    SAFE_ALLOCATE(sigma1(1:3, 1:this%nsrfcpnts))
    SAFE_ALLOCATE(sigma2(1:3, 1:this%nsrfcpnts))
    SAFE_ALLOCATE(sigma3(1:3))

    SAFE_ALLOCATE(integ1_s(1:st%nst, 1:sdim, 1:st%d%nik, 0:lmax, -lmax:lmax, 1:3, tdsteps_start:tdsteps_end))
    SAFE_ALLOCATE(integ2_s(1:st%nst, 1:sdim, 1:st%d%nik, 0:lmax, -lmax:lmax, tdsteps_start:tdsteps_end))

    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do isdim = 1, sdim

          do itstep = tdsteps_start, tdsteps_end

            do ll = 0, lmax 
              do mm = -ll, ll
                do imdim = 1, 3 
                  sigma1(imdim, :) = &
                    this%ylm_r(ll, mm, :) * this%srfcnrml(imdim, :) * &
                    this%wf(ist, isdim, ik, :, itstep)

                  sigma2(imdim, :) = &
                    this%ylm_r(ll, mm, :) * this%srfcnrml(imdim, :) * &
                    this%gwf(ist, isdim, ik, :, itstep, imdim)

                  ! surface integral
                  integ1_s(ist, isdim, ik, ll, mm, imdim, itstep) = sum(sigma1(imdim, :))
                  sigma3(imdim) = sum(sigma2(imdim, :))
                end do

                integ2_s(ist, isdim, ik, ll, mm, itstep) = sum(sigma3)

              end do
            end do

          end do

        end do
      end do
    end do

    SAFE_DEALLOCATE_A(sigma1)
    SAFE_DEALLOCATE_A(sigma2)
    SAFE_DEALLOCATE_A(sigma3)

    ! spectral amplitude for time range
    SAFE_ALLOCATE(integ11_t(1:this%nstepsomegak, 1:3))
    SAFE_ALLOCATE(integ21_t(1:this%nstepsomegak))
    SAFE_ALLOCATE(integ12_t(1:this%nk, 1:this%nstepsomegak, 1:3))
    SAFE_ALLOCATE(integ22_t(1:this%nk, 1:this%nstepsomegak))

    SAFE_ALLOCATE(spctramp_sph(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nk, 1:this%nstepsomegak))
    spctramp_sph = M_z0

    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do isdim = 1, sdim
  
          do itstep = tdsteps_start, tdsteps_end
  
            integ12_t = M_z0
            integ22_t = M_z0

            do ll = 0, lmax

              integ11_t = M_z0
              integ21_t = M_z0
  
              do mm = -ll, ll
                ! multiply with spherical harmonics & sum over all mm
                do imdim = 1, 3
                  integ11_t(1:this%nstepsomegak, imdim) = integ11_t(1:this%nstepsomegak, imdim) + &
                    integ1_s(ist, isdim, ik, ll, mm, imdim, itstep) * this%ylm_k(ll, mm, 1:this%nstepsomegak)
                end do
                integ21_t(1:this%nstepsomegak) = integ21_t(1:this%nstepsomegak) + &
                  integ2_s(ist, isdim, ik, ll, mm, itstep) * this%ylm_k(ll, mm, 1:this%nstepsomegak)
              end do
              ! multiply with Bessel function & sum over all ll
              do ikk = 1, this%nk
                integ12_t(ikk, 1:this%nstepsomegak, 1:3) = integ12_t(ikk, 1:this%nstepsomegak, 1:3) + &
                  integ11_t(1:this%nstepsomegak, 1:3) * this%j_l(ll, ikk) * ( - M_zI)**ll
                integ22_t(ikk, 1:this%nstepsomegak) = integ22_t(ikk, 1:this%nstepsomegak) + &
                  integ21_t(1:this%nstepsomegak) * this%j_l(ll, ikk) * ( - M_zI)**ll
              end do
            end do
            ! integrate over time & sum over dimensions
            do imdim = 1, 3
              spctramp_sph(ist, isdim, ik, :,:) = spctramp_sph(ist, isdim, ik, :,:) + &
                conjgphase_sph(:,:, itstep) * (integ12_t(:,:, imdim) * &
                 (M_TWO * this%veca(imdim, itstep)  / P_c - this%kcoords_sph(imdim, :,:)))
            end do
              spctramp_sph(ist, isdim, ik, :, :) = spctramp_sph(ist, isdim, ik, :,:) + &
                conjgphase_sph(:,:, itstep) * integ22_t(:,:) * M_zI
          end do

        end do
      end do
    end do

    SAFE_DEALLOCATE_A(integ1_s)
    SAFE_DEALLOCATE_A(integ2_s)

    SAFE_DEALLOCATE_A(integ11_t)
    SAFE_DEALLOCATE_A(integ12_t)
    SAFE_DEALLOCATE_A(integ21_t)
    SAFE_DEALLOCATE_A(integ22_t)

    SAFE_DEALLOCATE_A(conjgphase_sph)

    ! sum over time
#if defined(HAVE_MPI)
    call comm_allreduce(mpi_world%comm, spctramp_sph)
#endif
    this%spctramp_sph = this%spctramp_sph + spctramp_sph

    SAFE_DEALLOCATE_A(spctramp_sph)

    POP_SUB(pes_flux_integrate_sph)
  end subroutine pes_flux_integrate_sph

  ! ---------------------------------------------------------
  subroutine pes_flux_getcube(this, mesh, border, offset)
    type(mesh_t),     intent(in)    :: mesh
    type(pes_flux_t), intent(inout) :: this
    FLOAT,            intent(in)    :: border(1:MAX_DIM)
    FLOAT,            intent(in)    :: offset(1:MAX_DIM)

    integer, allocatable  :: which_surface(:)
    FLOAT                 :: xx(MAX_DIM), dd
    integer               :: mdim, imdim, idir, isp
    integer               :: ip_global, ip_start, ip_end
    integer               :: rankmin, nsurfaces

    PUSH_SUB(pes_flux_getcube)

    ! this routine is parallelized over the mesh in any case

    mdim = mesh%sb%dim

    call pes_flux_distribute(1, mesh%np_global, ip_start, ip_end)

    SAFE_ALLOCATE(which_surface(1:mesh%np_global))
    which_surface = 0

    ! get the surface points
    this%nsrfcpnts = 0
    do ip_global = ip_start, ip_end
      nsurfaces = 0
      xx(1:MAX_DIM) = mesh_x_global(mesh, ip_global) - offset(1:MAX_DIM)

      ! check whether the point is inside the cube
      if(all(abs(xx(1:mdim)) <= border(1:mdim))) then
        ! check whether the point is close to any border
        do imdim = 1, mdim
          dd = border(imdim) - abs(xx(imdim))
          if(dd < mesh%spacing(imdim)/M_TWO) then
            nsurfaces = nsurfaces + 1
            which_surface(ip_global) = int(sign(M_ONE, xx(imdim))) * imdim  ! +-x=+-1, +-y=+-2, +-z=+-3
          end if
        end do
        
        ! check whether the point is close to one border only (not an edge)
        if(nsurfaces == 1) then
          this%nsrfcpnts = this%nsrfcpnts + 1
        else
          which_surface(ip_global) = 0
        end if
      end if
    end do

#if defined(HAVE_MPI)
    call comm_allreduce(mpi_world%comm, this%nsrfcpnts)
    call comm_allreduce(mpi_world%comm, which_surface)
#endif

    SAFE_ALLOCATE(this%srfcpnt(1:this%nsrfcpnts))
    SAFE_ALLOCATE(this%srfcnrml(1:mdim, 1:this%nsrfcpnts))
    SAFE_ALLOCATE(this%rcoords(1:mdim, 1:this%nsrfcpnts))
    SAFE_ALLOCATE(this%rankmin(1:this%nsrfcpnts))

    this%srfcnrml = M_ZERO

    isp = 0
    do ip_global = 1, mesh%np_global
      if(which_surface(ip_global) /= 0) then
        isp = isp + 1
        ! coordinate of surface point
        xx(1:MAX_DIM) = mesh_x_global(mesh, ip_global)
        this%rcoords(1:mdim, isp) = xx(1:mdim)
        ! local ip & node which has the surface point
        this%srfcpnt(isp) = mesh_nearest_point(mesh, this%rcoords(1:mdim, isp), dd, rankmin)
        this%rankmin(isp) = rankmin
        ! surface normal
        idir = abs(which_surface(ip_global))
        this%srfcnrml(idir, isp) = sign(1, which_surface(ip_global))
        ! add the surface element (of directions orthogonal to the normal vector)
        do imdim = 1, mdim
          if(imdim == idir) cycle
          this%srfcnrml(idir, isp) = this%srfcnrml(idir, isp) * mesh%spacing(imdim)
        end do
      end if
    end do

    SAFE_DEALLOCATE_A(which_surface)

    POP_SUB(pes_flux_getcube)
  end subroutine pes_flux_getcube

  ! ---------------------------------------------------------
  subroutine pes_flux_getsphere(this, mesh, nstepsthetar, nstepsphir, offset)
    type(pes_flux_t), intent(inout) :: this
    type(mesh_t),     intent(in)    :: mesh
    integer,          intent(inout) :: nstepsthetar, nstepsphir
    FLOAT,            intent(in)    :: offset(1:MAX_DIM)
    
    integer :: mdim
    integer :: isp, ith, iph
    FLOAT   :: thetar, phir
    FLOAT   :: weight, dthetar

    PUSH_SUB(pes_flux_getsphere)

    mdim = mesh%sb%dim

    if(nstepsphir == 0) nstepsphir = 1

    select case(mdim)
    case(1)
      nstepsthetar = 0
      nstepsphir   = 2
      this%nsrfcpnts = nstepsphir
  
    case(2)
      nstepsthetar = 0
      this%nsrfcpnts = nstepsphir
  
    case(3)
      if(nstepsthetar <= 1) then
        nstepsphir = 1
        nstepsthetar = 1
      end if
      dthetar = M_PI / nstepsthetar
      this%nsrfcpnts  = nstepsphir * (nstepsthetar - 1) + 2
  
    end select

    SAFE_ALLOCATE(this%srfcnrml(1:mdim, 1:this%nsrfcpnts))
    SAFE_ALLOCATE(this%rcoords(1:mdim, 1:this%nsrfcpnts))

    ! initializing spherical grid
    thetar = M_PI / M_TWO
    isp = 0
    do ith = 0, nstepsthetar

      select case(mdim)
      case(1)
        weight = M_ONE

      case(2)
        weight = M_TWO * M_PI / nstepsphir

      case(3)
        thetar  = ith * dthetar

        if(ith == 0 .or. ith == nstepsthetar) then
          weight = (M_ONE - cos(dthetar/M_TWO)) * M_TWO * M_PI
        else
          weight = abs(cos(thetar - dthetar/M_TWO) - cos(thetar + dthetar/M_TWO)) &
            * M_TWO * M_PI / nstepsphir
        end if

      end select

      do iph = 0, nstepsphir - 1     ! 2*pi is the same than zero
        isp = isp + 1
        phir = iph * M_TWO * M_PI / nstepsphir
                      this%srfcnrml(1, isp) = cos(phir) * sin(thetar)
        if(mdim >= 2) this%srfcnrml(2, isp) = sin(phir) * sin(thetar)
        if(mdim == 3) this%srfcnrml(3, isp) = cos(thetar)
        this%rcoords(1:mdim, isp) = this%radius * this%srfcnrml(1:mdim, isp)
        ! here we also include the surface elements
        this%srfcnrml(1:mdim, isp) = weight * this%srfcnrml(1:mdim, isp)
        if(thetar == M_ZERO .or. thetar == M_PI) exit
      end do
    end do

    POP_SUB(pes_flux_getsphere)
  end subroutine pes_flux_getsphere

  ! ---------------------------------------------------------
  subroutine pes_flux_output(this, mesh, sb, st, dt)
    type(pes_flux_t), intent(inout)    :: this
    type(mesh_t),        intent(in)    :: mesh
    type(simul_box_t),   intent(in)    :: sb
    type(states_t),      intent(in)    :: st
    FLOAT,               intent(in)    :: dt

    integer            :: sdim, mdim
    integer            :: ist, ik, isdim
    integer            :: ikp, iomk
    integer            :: ikk, ith, iph, iphi
    FLOAT              :: phik, thetak, kk, kact

    integer            :: iunit
    FLOAT, allocatable :: spctrout_cub(:), spctrout_sph(:,:)
    FLOAT              :: spctroutsave
    FLOAT              :: weight, spctrsum

    PUSH_SUB(pes_flux_output)

    sdim   = st%d%dim
    mdim   = mesh%sb%dim

    if(this%shape == M_SPHERICAL) then
      SAFE_ALLOCATE(spctrout_sph(1:this%nk, 1:this%nstepsomegak))
      spctrout_sph = M_ZERO
    else
      SAFE_ALLOCATE(spctrout_cub(1:this%nkpnts))
      spctrout_cub = M_ZERO
    end if

    ! calculate the total spectrum
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do isdim = 1, sdim
          if(this%shape == M_SPHERICAL) then
            spctrout_sph(1:this%nk, 1:this%nstepsomegak) = spctrout_sph(1:this%nk, 1:this%nstepsomegak) + &
              abs(this%spctramp_sph(ist, isdim, ik, 1:this%nk, 1:this%nstepsomegak))**M_TWO * (dt * this%tdstepsinterval)**M_TWO
          else
            spctrout_cub(1:this%nkpnts) = spctrout_cub(1:this%nkpnts) + &
              abs(this%spctramp_cub(ist, isdim, ik, 1:this%nkpnts))**M_TWO * (dt * this%tdstepsinterval)**M_TWO
          end if
        end do
      end do
    end do

    if(mpi_grp_is_root(mpi_world)) then
      iunit = io_open('td.general/PES_flux.distribution.out', action='write', position='rewind')
      write(iunit, '(a29)') '# k, theta, phi, distribution'

      if(this%shape == M_SPHERICAL) then
        do ikk = 1, this%nk 
          kact = ikk * this%dk
          iomk = 0
          do ith = 0, this%nstepsthetak
            thetak = ith * this%dthetak
            do iph = 0, this%nstepsphik - 1
              iomk = iomk + 1
              if(iph == 0) spctroutsave = spctrout_sph(ikk, iomk)
              phik = iph * M_TWO * M_PI / this%nstepsphik
              write(iunit,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_sph(ikk, iomk)

              ! just repeat the result for output
              if(iph == (this%nstepsphik - 1)) &
                write(iunit,'(5(1x,e18.10E3))') kact, thetak, M_TWO * M_PI, spctroutsave

              ! just repeat the result for output
              if(ith == 0 .or. ith == this%nstepsthetak) then
                do iphi = 1, this%nstepsphik
                  phik = iphi * M_TWO * M_PI / this%nstepsphik
                  write(iunit,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_sph(ikk, iomk)
                end do
                exit
              end if
            end do
            write(iunit, '(1x)', advance='yes')
          end do
          write(iunit, '(1x)', advance='yes')
        end do

      else
        thetak = M_PI / M_TWO
        ikp    = 0
        do ikk = 1, this%nk
          kk = ikk * this%dk
          do ith = 0, this%nstepsthetak
            if(mdim == 3) thetak = ith * this%dthetak
            do iph = 0, this%nstepsphik - 1
              ikp = ikp + 1
              if(iph == 0) spctroutsave = spctrout_cub(ikp)
              phik = iph * M_TWO * M_PI / this%nstepsphik
              write(iunit,'(5(1x,e18.10E3))') kk, thetak, phik, spctrout_cub(ikp)
   
              ! just repeat the result for output
              if(iph == (this%nstepsphik - 1)) &
                write(iunit,'(5(1x,e18.10E3))') kk, thetak, M_TWO * M_PI, spctroutsave
               
              ! just repeat the result for output
              if(thetak == M_ZERO .or. thetak == M_PI) then
                do iphi = 1, this%nstepsphik
                  phik = iphi * M_TWO * M_PI / this%nstepsphik
                  write(iunit,'(5(1x,e18.10E3))') kk, thetak, phik, spctrout_cub(ikp)
                end do
                exit
              end if
            end do
            if(this%nstepsphik > 0) write(iunit, '(1x)', advance='yes')
          end do
          if(this%nstepsphik == 0) write(iunit, '(1x)', advance='yes')
        end do
      end if
      call io_close(iunit)
    end if

    if(mpi_grp_is_root(mpi_world)) then
      iunit = io_open('td.general/PES_flux.power.sum', action='write', position='rewind')
      write(iunit, '(a29)') '# E,  distribution'

      if(this%shape == M_SPHERICAL) then
        do ikk = 1, this%nk
          kact = ikk * this%dk
          iomk = 0
          spctrsum = M_ZERO

          do ith = 0, this%nstepsthetak
            thetak = ith * this%dthetak

            if(ith == 0 .or. ith == this%nstepsthetak) then
              weight = (M_ONE - cos(this%dthetak/M_TWO)) * M_TWO * M_PI
            else
              weight = abs(cos(thetak - this%dthetak/M_TWO) - cos(thetak + this%dthetak/M_TWO)) &
                * M_TWO * M_PI / this%nstepsphik
            end if

            do iph = 0, this%nstepsphik - 1
              iomk = iomk + 1
              spctrsum = spctrsum + spctrout_sph(ikk, iomk) * weight

              if(thetak == M_ZERO .or. thetak == M_PI) exit

            end do
          end do

          write(iunit,'(2(1x,e18.10E3))') kact**M_TWO / M_TWO, spctrsum * kact
        end do
      end if
      call io_close(iunit)
    end if

    SAFE_DEALLOCATE_A(spctrout_cub)
    SAFE_DEALLOCATE_A(spctrout_sph)

    POP_SUB(pes_flux_output)
  end subroutine pes_flux_output

  ! ---------------------------------------------------------
  subroutine pes_flux_dump(restart, this, mesh, st, ierr)
    type(restart_t),  intent(in)  :: restart
    type(pes_flux_t), intent(in)  :: this
    type(mesh_t),     intent(in)  :: mesh
    type(states_t),   intent(in)  :: st
    integer,          intent(out) :: ierr

    integer :: err, mdim, sdim, rstdim

    PUSH_SUB(pes_flux_dump)

    mdim   = mesh%sb%dim
    sdim   = st%d%dim
    rstdim = st%nst * sdim * st%d%nik

    if(restart_skip(restart)) then
      POP_SUB(pes_flux_dump)
      return
    end if

    if(in_debug_mode) then
      message(1) = "Debug: Writing pes_flux restart."
      call messages_info(1)
    end if

!    call zrestart_write_binary(restart, 'pesflux1', rstdim * this%nkpnts, this%spctramp, err)
!    call zrestart_write_binary(restart, 'pesflux2', rstdim * this%nsrfcpnts * this%tdsteps, this%wf, err)
!    call zrestart_write_binary(restart, 'pesflux3', rstdim * this%nsrfcpnts * this%tdsteps * mdim, this%gwf, err)
!    call zrestart_write_binary(restart, 'pesflux4', (this%tdsteps + 1) * this%nkpnts, this%conjgphase, err)

    if(err /= 0) ierr = ierr + 1

    if(in_debug_mode) then
      message(1) = "Debug: Writing pes_flux restart done."
      call messages_info(1)
    end if

    POP_SUB(pes_flux_dump)
  end subroutine pes_flux_dump

  ! ---------------------------------------------------------
  subroutine pes_flux_load(restart, this, mesh, st, ierr)
    type(restart_t),     intent(in)    :: restart
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(in)    :: st
    integer,             intent(out)   :: ierr

    integer :: err, mdim, sdim, rstdim

    PUSH_SUB(pes_flux_load)

    mdim   = mesh%sb%dim
    sdim   = st%d%dim
    rstdim = st%nst * sdim * st%d%nik

    if(restart_skip(restart)) then
      ierr = -1
      POP_SUB(pes_flux_load)
      return
    end if

    if(in_debug_mode) then
      message(1) = "Debug: Reading pes_flux restart."
      call messages_info(1)
    end if

!    call zrestart_read_binary(restart, 'pesflux1', rstdim * this%nkpnts, this%spctramp, err)
!    call zrestart_read_binary(restart, 'pesflux2', rstdim * this%nsrfcpnts * this%tdsteps, this%wf, err)
!    call zrestart_read_binary(restart, 'pesflux3', rstdim * this%nsrfcpnts * this%tdsteps * mdim, this%gwf, err)
!    call zrestart_read_binary(restart, 'pesflux4', (this%tdsteps + 1) * this%nkpnts, this%conjgphase, err)

    if(err /= 0) ierr = ierr + 1
   
    if(in_debug_mode) then
      message(1) = "Debug: Reading pes_flux restart done."
      call messages_info(1)
    end if

    POP_SUB(pes_flux_load)
  end subroutine pes_flux_load

  ! ---------------------------------------------------------
  subroutine pes_flux_distribute(istart_global, iend_global, istart, iend)
    integer,          intent(in)    :: istart_global
    integer,          intent(in)    :: iend_global
    integer,          intent(inout) :: istart
    integer,          intent(inout) :: iend

#if defined(HAVE_MPI)
    integer, allocatable :: dimrank(:)
    integer              :: mpisize, mpirank, irest, irank
    integer              :: inumber
#endif

    PUSH_SUB(pes_flux_distribute)

#if defined(HAVE_MPI)
    call mpi_comm_size(mpi_comm_world, mpisize, mpi_err)
    call mpi_comm_rank(mpi_comm_world, mpirank, mpi_err)

    SAFE_ALLOCATE(dimrank(0:mpisize-1))

    inumber = iend_global - istart_global + 1
    irest = mod(inumber, mpisize)
    do irank = 0, mpisize - 1
      if(irest > 0) then
        dimrank(irank) = int(inumber / mpisize) + 1
        irest = irest - 1
      else
        dimrank(irank) = int(inumber / mpisize)
      end if
    end do

    iend   = istart_global + sum(dimrank(:mpirank)) - 1
    istart = iend - dimrank(mpirank) + 1

    SAFE_DEALLOCATE_A(dimrank)

#else
    istart = istart_global
    iend   = iend_global
#endif

    POP_SUB(pes_flux_distribute)
  end subroutine pes_flux_distribute

end module pes_flux_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
