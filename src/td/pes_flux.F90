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

    integer          :: shape                          !< shape of the surface (= cubic/spherical)
    integer          :: nsrfcpnts                      !< total number of surface points
    integer          :: nsrfcpnts_start, nsrfcpnts_end !< number of surface points on node
    FLOAT, pointer   :: srfcnrml(:,:)                  !< vectors normal to the surface (includes surface element)
    FLOAT, pointer   :: rcoords(:,:)                   !< coordinates of the surface points
    integer, pointer :: srfcpnt(:)                     !< for cubic surface: returns index of the surface points
    integer, pointer :: rankmin(:)                     !< for cubic surface: returns node which has the surface point

    integer          :: tdsteps                        !<
    integer          :: tdstepsinterval                !<

    CMPLX, pointer   :: wf(:,:,:,:,:)                  !< wavefunction
    CMPLX, pointer   :: gwf(:,:,:,:,:,:)               !< gradient of wavefunction
    FLOAT, pointer   :: veca(:,:)                      !< vector potential
    CMPLX, pointer   :: conjgphase_prev_cub(:)         !< Volkov phase for all k-points from previous time step
    CMPLX, pointer   :: conjgplanewf_cub(:,:)          !< plane wave factor
    CMPLX, pointer   :: spctramp_cub(:,:,:,:)          !< spectral amplitude

    logical          :: usememory                      !< whether conjgplanewf should be kept in memory
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
    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: imdim
    integer            :: isp, ikp
    integer            :: il
    integer            :: ikk, ith, iph
    FLOAT              :: kmax, thetak, phik

    FLOAT, allocatable :: k_dot_aux(:)
#if defined(HAVE_MPI)
    integer            :: mpisize, mpirank
#endif
    FLOAT              :: kact

    PUSH_SUB(pes_flux_init)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
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

    message(1) = 'Info: Calculation PES using t-surff technique.'
    call messages_info(1)

    ! -----------------------------------------------------------------
    ! Setting up r-mesh (the surface)
    ! -----------------------------------------------------------------
    !%Variable PES_Flux_Shape
    !%Type integer
    !%Default cub
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The shape of the surface.
    !%Option cub 1
    !% Uses a parallelepiped as surface. All surface points are grid points.
    !% Choose the location of the surface with variable PES_Flux_Lsize.
    !%Option sph 2
    !% Constructs a sphere with radius PES_Flux_Radius. Points on the sphere 
    !% are interpolated by triangular interpolation.
    !%End
    call parse_variable('PES_Flux_Shape', M_CUBIC, this%shape)
    if(.not.varinfo_valid_option('PES_Flux_Shape', this%shape, is_flag = .true.)) &
      call messages_input_error('PES_Flux_Shape')
    call messages_print_var_option(stdout, 'PES_Flux_Shape', this%shape)

    if(this%shape == M_SPHERICAL) then
      message(1) = 'Spherical surface not yet implemented.'
      call messages_fatal(1)
    end if

    !%Variable PES_Flux_Offset
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End
    offset = M_ZERO
    if(parse_block('PES_Flux_Offset', blk) == 0) then
      call parse_block_float(blk, 0, 0, offset(1))
      call parse_block_float(blk, 0, 1, offset(2))
      call parse_block_float(blk, 0, 2, offset(3))
    end if

    !%Variable PES_Flux_Lsize
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End
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

    ! get the surface points
    call pes_flux_getcube(this, mesh, border, offset)

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

    SAFE_ALLOCATE(this%kcoords_cub(1:this%nkpnts, 1:mdim))
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
                        this%kcoords_cub(ikp, 1) = kact * cos(phik) * sin(thetak)
          if(mdim >= 2) this%kcoords_cub(ikp, 2) = kact * sin(phik) * sin(thetak)
          if(mdim == 3) this%kcoords_cub(ikp, 3) = kact * cos(thetak)
          if(thetak == M_ZERO .or. thetak == M_PI) exit
        end do
      end do
    end do

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
    call parse_variable('PES_Flux_TDSteps', 1, this%tdsteps)

    ! -----------------------------------------------------------------
    ! Other stuff
    ! -----------------------------------------------------------------
    !%Variable PES_Flux_UseMemory
    !%Type logical
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End
    call parse_variable('PES_Flux_UseMemory', .true., this%usememory)

    SAFE_ALLOCATE(this%wf(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nsrfcpnts, 1:this%tdsteps))
    this%wf = M_z0

    SAFE_ALLOCATE(this%gwf(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nsrfcpnts, 1:this%tdsteps, 1:mdim))
    this%gwf = M_z0

    SAFE_ALLOCATE(this%veca(1:this%tdsteps, 1:mdim))
    this%veca = M_ZERO

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
          k_dot_aux(:) = k_dot_aux(:) + this%kcoords_cub(:, imdim) * this%rcoords(imdim, isp)
        end do
        this%conjgplanewf_cub(:, isp) = exp(-M_zI * k_dot_aux(:)) / (M_TWO * M_PI)**(mdim/M_TWO)
      end do
      SAFE_DEALLOCATE_A(k_dot_aux)
#if defined(HAVE_MPI)
      call comm_allreduce(mpi_world%comm, this%conjgplanewf_cub)
#endif
    end if

    POP_SUB(pes_flux_init)
  end subroutine pes_flux_init

  ! ---------------------------------------------------------
  subroutine pes_flux_end(this)
    type(pes_flux_t), intent(inout) :: this

    PUSH_SUB(pes_flux_end)

    SAFE_DEALLOCATE_P(this%kcoords_cub)
    SAFE_DEALLOCATE_P(this%conjgphase_prev_cub)
    SAFE_DEALLOCATE_P(this%spctramp_cub)
    if(this%usememory) then
      SAFE_DEALLOCATE_P(this%conjgplanewf_cub)
    end if

    SAFE_DEALLOCATE_P(this%srfcpnt)
    SAFE_DEALLOCATE_P(this%rankmin)

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
        call laser_field(hm%ep%lasers(il), this%veca(itstep, 1:mdim), iter*dt)
      end do
      this%veca(itstep, :) = - this%veca(itstep, :)

      ! save wavefunctions & gradients
      do ik = kptst, kptend
        do ist = stst, stend
          do isdim = 1, sdim
            call states_get_state(st, mesh, isdim, ist, ik, psi)
            call zderivatives_grad(gr%der, psi, gpsi, .true.)

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
                this%gwf(ist, isdim, ik, isp, itstep, 1:mdim) = st%occ(ist, ik) * gpsi(ip_local, 1:mdim) ! * mesh%spacing(1:mdim) ?
              end if
            end do
          end do
        end do
      end do

#if defined(HAVE_MPI)
      call comm_allreduce(mpi_world%comm, this%wf(:,:,:,:, itstep))
      do imdim = 1, mdim
        call comm_allreduce(mpi_world%comm, this%gwf(:,:,:,:, itstep, imdim))
      end do
#endif

      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(gpsi)
   
      if(iter > 0) then
        if(mod(iter, this%tdsteps * this%tdstepsinterval) == 0 .or. iter == maxiter) then
          call pes_flux_integrate_cub(this, mesh, st, dt)
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
        vec = sum((this%kcoords_cub(ikp, 1:mdim) - this%veca(itstep, 1:mdim) / P_c)**2)
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
        if(this%srfcnrml(isp, idir) == M_ZERO) cycle

        Jk_cub = M_z0

        if(.not. this%usememory) then
          k_dot_aux(:) = M_ZERO
          do imdim = 1, mdim
            k_dot_aux(:) = k_dot_aux(:) + this%kcoords_cub(:, imdim) * this%rcoords(imdim, isp)
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
                   (M_TWO * this%veca(itstep, idir) / P_c - this%kcoords_cub(1:this%nkpnts,idir)) + &
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
        spctramp_cub(:,:,:,:) = spctramp_cub(:,:,:,:) + Jk_cub(:,:,:,:) * this%srfcnrml(isp, idir)
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
    SAFE_ALLOCATE(this%srfcnrml(1:this%nsrfcpnts, 1:mdim))
    SAFE_ALLOCATE(this%rcoords(1:mdim, 1:this%nsrfcpnts))
    SAFE_ALLOCATE(this%rankmin(1:this%nsrfcpnts))

    this%srfcnrml = M_ZERO

    isp = 0
    do ip_global = 1, mesh%np_global
      if(which_surface(ip_global) /= 0) then
        isp = isp + 1
        ! coordinate of surface point
        this%rcoords(1:mdim, isp) = mesh_x_global(mesh, ip_global) 
        ! local ip & node which has the surface point
        this%srfcpnt(isp) = mesh_nearest_point(mesh, this%rcoords(1:mdim, isp), dd, rankmin)
        this%rankmin(isp) = rankmin
        ! surface normal
        idir = abs(which_surface(ip_global))
        this%srfcnrml(isp, idir) = sign(1, which_surface(ip_global))
        ! add the surface element (of directions orthogonal to the normal vector)
        do imdim = 1, mdim
          if(imdim == idir) cycle
          this%srfcnrml(isp, idir) = this%srfcnrml(isp, idir) * mesh%spacing(imdim)
        end do
      end if
    end do

    SAFE_DEALLOCATE_A(which_surface)

    POP_SUB(pes_flux_getcube)
  end subroutine pes_flux_getcube

  ! ---------------------------------------------------------
  subroutine pes_flux_output(this, mesh, sb, st, dt)
    type(pes_flux_t), intent(inout)    :: this
    type(mesh_t),        intent(in)    :: mesh
    type(simul_box_t),   intent(in)    :: sb
    type(states_t),      intent(in)    :: st
    FLOAT,               intent(in)    :: dt

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim
    integer            :: ikp
    integer            :: ikk, ith, iph, iphi
    FLOAT              :: phik, thetak, kk

    integer            :: iunit
    FLOAT, allocatable :: spctrout_cub(:)
    FLOAT              :: spctroutsave

    PUSH_SUB(pes_flux_output)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = mesh%sb%dim

    SAFE_ALLOCATE(spctrout_cub(1:this%nkpnts))
    spctrout_cub = M_ZERO

    ! calculate the total spectrum
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do isdim = 1, sdim
          spctrout_cub(1:this%nkpnts) = spctrout_cub(1:this%nkpnts) + &
            abs(this%spctramp_cub(ist, isdim, ik, 1:this%nkpnts))**M_TWO * (dt * this%tdstepsinterval)**M_TWO
        end do
      end do
    end do

    if(mpi_grp_is_root(mpi_world)) then
      iunit = io_open('td.general/PES_flux.distribution.out', action='write', position='rewind')
      write(iunit, '(a29)') '# k, theta, phi, distribution'

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

    SAFE_DEALLOCATE_A(spctrout_cub)

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
#endif
    integer              :: inumber

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
