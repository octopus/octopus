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
    integer           :: nkpnts_global         !< global number of k-points
    integer           :: nkpnts_start          !< start of index for k-points on the current node
    integer           :: nkpnts_end
    integer           :: nk, nphi, ntheta
    FLOAT             :: delk, phimin, delphi
    FLOAT             :: thetamin, deltheta
    integer           :: srfcshape             !< shape of the surface (= cubic/spherical)
    integer           :: nsrfcpnts             !< total number of points contructing surface
    integer           :: tdsteps
    integer           :: tdstepsinterval

    integer, pointer  :: nsrfcpnt_start(:)
    integer, pointer  :: nsrfcpnt_end(:)
    FLOAT, pointer    :: kpnt(:,:)             !< coordinates of all k-points
    integer, pointer  :: srfcpnt(:)            !< returns the (global) index of the points on the surface
    FLOAT, pointer    :: srfcnrml(:,:)         !< (unit) vectors normal to the surface
    CMPLX, pointer    :: spctramp(:)           !< spectral amplitude
    CMPLX, pointer    :: conjgphase(:,:)       !< current Volkov phase for all k-points
    CMPLX, pointer    :: conjgplanewf(:,:)
    CMPLX, pointer    ::  wf(:,:)
    CMPLX, pointer    :: gwf(:,:,:)
    integer, pointer  :: mpirank_srfcpnt(:)    !< partition node of current surface point (for par_domains)
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

    type(block_t)        :: blk
    FLOAT                :: border(MAX_DIM)       ! distance of surface from border
    FLOAT                :: offset(MAX_DIM)       ! offset for border
    integer              :: ikp, isp, il, ith, iph, ikk, dim, ii, irest
    FLOAT                :: phi, theta, kk, krr
    FLOAT                :: phimax, thetamax, kmax
    FLOAT                :: xx(MAX_DIM)
#if defined(HAVE_MPI)
    integer, allocatable :: dimrank(:)
    integer              :: mpisize, mpirank, ierr
#endif

    PUSH_SUB(pes_flux_init)

    dim    = mesh%sb%dim
    offset = M_ZERO

    call messages_experimental("PhotoElectronSpectrum with t-surff")

#if defined(HAVE_MPI)
    call mpi_comm_size(mpi_comm_world, mpisize, ierr)
    call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
#endif

    do il = 1, hm%ep%no_lasers
      if(laser_kind(hm%ep%lasers(il)) /= E_FIELD_VECTOR_POTENTIAL) then
        message(1) = 't-surff only works in velocity gauge.'
        call messages_fatal(1)
      end if
    end do

    message(1) = 'Info: Calculation PES using t-surff technique.'
    call messages_info(1)

    ! surface
    if(parse_block('PESSurface', blk) < 0) then
      select case(mesh%sb%box_shape)
      case(PARALLELEPIPED)
        border(1:dim) = mesh%sb%lsize(1:dim)
      case(SPHERE)
        border(1:dim) = mesh%sb%rsize/sqrt(M_TWO)
      case default
        message(1) = &
          "PESSurface not specified. No default values available for this box shape. Specify a surface with block PESSurface."
        call messages_fatal(1)
      end select
      message(1) = "PESSurface not specified. Using default values."
      call messages_info(1)
    else
      call parse_block_float(blk, 0, 0, border(1))
      call parse_block_float(blk, 0, 1, border(2))
      call parse_block_float(blk, 0, 2, border(3))
      border(1:dim) = int(border(1:dim)/mesh%spacing(1:dim))*mesh%spacing(1:dim)
    end if
    !%Variable PESSurface
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    if(parse_block('PESSurfaceOffset', blk) == 0) then
      call parse_block_float(blk, 0, 0, offset(1))
      call parse_block_float(blk, 0, 1, offset(2))
      call parse_block_float(blk, 0, 2, offset(3))
    end if
    !%Variable PESSurfaceOffset
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call pes_flux_getsrfc(this, mesh, border, offset)

    ! k-mesh in 1D (2 points), 2D (polar coordinates), & 3D (spherical
    ! coordinates)
    call parse_variable('PESSurfaceKmax', M_ONE, kmax)
    !%Variable PESSurfaceKmax
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call parse_variable('PESSurfaceDeltaK', CNST(0.002), this%delk)
    !%Variable PESSurfaceDeltaK
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call parse_variable('PESSurfacePhiMin', M_ZERO, this%phimin)
    !%Variable PESSurfacePhiMin
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call parse_variable('PESSurfacePhiMax', M_TWO * M_PI, phimax)
    !%Variable PESSurfacePhiMax
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call parse_variable('PESSurfaceDelPhi', CNST((M_TWO * M_PI)/360), this%delphi)
    !%Variable PESSurfaceDelPhi
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call parse_variable('PESSurfaceThetaMin', M_ZERO, this%thetamin)
    !%Variable PESSurfaceThetaMin
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call parse_variable('PESSurfaceThetaMax', M_PI, thetamax)
    !%Variable PESSurfaceThetaMax
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call parse_variable('PESSurfaceDelTheta', CNST(M_PI/180), this%deltheta)
    !%Variable PESSurfaceDelTheta
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    select case(dim)
    case(1)
      phimax = M_PI
      this%phimin = M_ZERO
      this%delphi = M_PI
      this%thetamin = M_PI / M_TWO
      thetamax = M_PI / M_TWO
    case(2)
      thetamax = M_PI / M_TWO
      this%thetamin = M_PI / M_TWO
      this%deltheta = M_ONE
    end select

    this%nphi   = nint((phimax - this%phimin)/this%delphi)
    this%ntheta = nint((thetamax - this%thetamin)/this%deltheta)
    this%nk     = nint(kmax/this%delk)
    this%nkpnts_global = (this%nphi + 1) * (this%ntheta + 1) * this%nk 

    write(*,*) "Info: number of k-points:", this%nkpnts_global

    ! k-points (global)
    SAFE_ALLOCATE(this%kpnt(1:this%nkpnts_global, 1:dim))
    this%kpnt = M_ZERO

    ikp = 0
    do ith = 0, this%ntheta
      theta = ith * this%deltheta + this%thetamin
      do iph = 0, this%nphi
        phi = iph * this%delphi + this%phimin
        do ikk = 1, this%nk
          kk = ikk * this%delk
          ikp = ikp + 1
                       this%kpnt(ikp, 1) = kk * cos(phi) * sin(theta)
          if(dim >= 2) this%kpnt(ikp, 2) = kk * sin(phi) * sin(theta)
          if(dim == 3) this%kpnt(ikp, 3) = kk * cos(theta)
        end do
      end do
    end do

    ! other stuff
    call parse_variable('PESSurfaceTDStepsInterval', 1, this%tdstepsinterval)
    !%Variable PESSurfaceTDStepsInterval
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    call parse_variable('PESSurfaceTDSteps', 1, this%tdsteps)
    !%Variable PESSurfaceTDSteps
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%End

    SAFE_ALLOCATE(this%wf(1:this%nsrfcpnts, 1:this%tdsteps))
    this%wf = M_z0

    SAFE_ALLOCATE(this%gwf(1:this%nsrfcpnts, 1:this%tdsteps, 1:dim))
    this%gwf = M_z0

    SAFE_ALLOCATE(this%spctramp(1:this%nkpnts_global))
    this%spctramp = M_z0

    ! distribute the k-points on the nodes
#if defined(HAVE_MPI)
    SAFE_ALLOCATE(dimrank(0:mpisize-1))
    irest = mod(this%nkpnts_global, mpisize)
    do ii = 0, mpisize-1
      if(irest > 0) then
        dimrank(ii) = int(this%nkpnts_global / mpisize) + 1
        irest = irest - 1
      else
        dimrank(ii) = int(this%nkpnts_global / mpisize)
      end if
    end do

    this%nkpnts_end   = sum(dimrank(:mpirank))
    this%nkpnts_start = this%nkpnts_end - dimrank(mpirank) + 1
    SAFE_DEALLOCATE_A(dimrank)

!    write(*,*) mpirank, st%d%kpt%start, st%d%kpt%end, st%st_start, st%st_end
!    write(*,*) mpirank, this%nkpnts_start, this%nkpnts_end
#else
    this%nkpnts_end   = this%nkpnts_global
    this%nkpnts_start = 1

!    write(*,*) st%d%kpt%start, st%d%kpt%end, st%st_start, st%st_end
#endif

    SAFE_ALLOCATE(this%conjgphase(1:this%nkpnts_global, 0:this%tdsteps))
    this%conjgphase(:,:) = M_z0
    this%conjgphase(this%nkpnts_start:this%nkpnts_end, :) = M_z1

    SAFE_ALLOCATE(this%conjgplanewf(this%nkpnts_start:this%nkpnts_end, 1:this%nsrfcpnts))
    do ikp = this%nkpnts_start, this%nkpnts_end
      do isp = 1, this%nsrfcpnts
        ! the sign of the phase here should be -1 of the one for the Volkov phase.
        xx = mesh_x_global(mesh, this%srfcpnt(isp))
        krr = dot_product(this%kpnt(ikp, 1:dim), xx(1:dim))
        this%conjgplanewf(ikp, isp) = exp(-M_zI * krr) / (M_TWO * M_PI)**(dim/M_TWO)
      end do
    end do

    POP_SUB(pes_flux_init)
  end subroutine pes_flux_init

  ! ---------------------------------------------------------
  subroutine pes_flux_end(this)
    type(pes_flux_t), intent(inout) :: this

    PUSH_SUB(pes_flux_end)

    SAFE_DEALLOCATE_P(this%kpnt)
    SAFE_DEALLOCATE_P(this%spctramp)
    SAFE_DEALLOCATE_P(this%conjgphase)
    SAFE_DEALLOCATE_P(this%conjgplanewf)
    
    SAFE_DEALLOCATE_P(this%srfcpnt)
    SAFE_DEALLOCATE_P(this%srfcnrml)

    SAFE_DEALLOCATE_P(this%nsrfcpnt_start)
    SAFE_DEALLOCATE_P(this%nsrfcpnt_end)

    SAFE_DEALLOCATE_P(this%wf)
    SAFE_DEALLOCATE_P(this%gwf)

    SAFE_DEALLOCATE_P(this%mpirank_srfcpnt)

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

    integer            :: dim, itstep, ik, ist, idim, isp, il, ikp, dir, mpidim
    CMPLX, allocatable :: gpsi(:,:), psi(:)
    CMPLX, allocatable :: wfact(:), gwfact(:,:)
    FLOAT              :: vp(1:MAX_DIM)
    integer            :: ierr
    FLOAT              :: vec
    integer            :: mpirank, ip

    PUSH_SUB(pes_flux_save)

    if(iter > 0 .and. mod(iter, this%tdstepsinterval) == 0) then
      dim  = mesh%sb%dim
   
      SAFE_ALLOCATE(psi(1:mesh%np_part))
      SAFE_ALLOCATE(gpsi(1:mesh%np, 1:dim))

      SAFE_ALLOCATE( wfact(1:this%nsrfcpnts))
      SAFE_ALLOCATE(gwfact(1:this%nsrfcpnts, 1:dim))

       wfact = M_z0
      gwfact = M_z0
   
      if(iter > 0 .and. mod(int(iter/this%tdstepsinterval), this%tdsteps) == 0) then
        itstep = this%tdsteps
      else
        itstep = mod(int(iter/this%tdstepsinterval), this%tdsteps) 
      end if

      ! clean up fields when a new cycle begins
      if(itstep == 1) then
         this%wf = M_z0
        this%gwf = M_z0
      end if

      ! get current laser field
      vp = M_ZERO
      do il = 1, hm%ep%no_lasers
        call laser_field(hm%ep%lasers(il), vp(1:dim), iter*dt) 
      end do
      vp(1:dim) = -vp(1:dim)
   
      ! save Volkov phase using the previous time step
      do ikp = this%nkpnts_start, this%nkpnts_end
        vec = sum((this%kpnt(ikp, 1:dim) - vp(1:dim) / P_c)**2)
        this%conjgphase(ikp, itstep) = this%conjgphase(ikp, itstep - 1) * exp(M_zI * vec * dt * this%tdstepsinterval / M_TWO)
      end do
      if(itstep == this%tdsteps) &
        this%conjgphase(this%nkpnts_start:this%nkpnts_end, 0) = this%conjgphase(this%nkpnts_start:this%nkpnts_end, itstep)

      ! save wavefunctions & gradients
      mpirank = 0
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            call states_get_state(st, mesh, idim, ist, ik, psi)
            call zderivatives_grad(gr%der, psi, gpsi)
            do isp = 1, this%nsrfcpnts
              ip = this%srfcpnt(isp)
#if defined(HAVE_MPI)
              if(mesh%parallel_in_domains) then
                call mpi_comm_rank(mpi_comm_world, mpirank, ierr)
                ip = vec_global2local(mesh%vp, this%srfcpnt(isp), mesh%vp%partno)
              end if
#endif              
              if(mpirank == this%mpirank_srfcpnt(isp)) then
                 wfact(isp)        =  wfact(isp)        + st%occ(ist, ik) * psi(ip)
                gwfact(isp, 1:dim) = gwfact(isp, 1:dim) + st%occ(ist, ik) * &
                  (psi(ip) * M_TWO * vp(1:dim) / P_c + gpsi(ip, 1:dim) * M_zI) ! * mesh%spacing(1:dim))
              else
                 wfact(isp) = M_z0
                gwfact(isp, 1:dim) = M_z0
              end if
            end do
          end do
        end do
      end do
   
      ! communicate wf, gwf
#if defined(HAVE_MPI)
      mpidim = this%nsrfcpnts
      call MPI_Allreduce( wfact(:),    this%wf(:,itstep),   mpidim,     MPI_CMPLX, mpi_sum, &
        mpi_comm_world, mpi_err)
      call MPI_Allreduce(gwfact(:,:), this%gwf(:,itstep,:), mpidim*dim, MPI_CMPLX, mpi_sum, &
        mpi_comm_world, mpi_err)
#else
      this%wf(1:this%nsrfcpnts,itstep)    =  wfact(1:this%nsrfcpnts)
      this%gwf(1:this%nsrfcpnts,itstep,1:dim) = gwfact(1:this%nsrfcpnts,1:dim)
#endif

      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(gpsi)
      SAFE_DEALLOCATE_A(wfact)
      SAFE_DEALLOCATE_A(gwfact)

      if(iter > 0) then
        if(mod(iter, this%tdsteps * this%tdstepsinterval) == 0 .or. iter == maxiter) &
          call pes_flux_integrate(this, mesh, st)
      end if

    end if

!      ! accumulate flux through surface
!      do ik = st%d%kpt%start, st%d%kpt%end
!        do ist = st%st_start, st%st_end
!          do idim = 1, st%d%dim
!            call states_get_state(st, mesh, idim, ist, ik, psi)
!            call zderivatives_grad(gr%der, psi, gpsi)
!            do isp = 1, this%nsrfcpnts
!               wf(isp) = st%occ(ist, ik) * psi(this%srfcpnt(isp))
!              gwf(isp, 1:dim) = st%occ(ist, ik) * gpsi(this%srfcpnt(isp), 1:dim) 
!            end do
!            do dir = 1, dim
!              do isp = start(dir), end(dir)
!                this%Jk(1:this%nkpnts, idim, ist, ik, isp, dir) =    &
!                  this%Jk(1:this%nkpnts, idim, ist, ik, isp, dir) +  &
!                  conjg(this%conjgplanewf(1:this%nkpnts, isp)) *             & 
!                  (    this%kpnt(1:this%nkpnts, dir) * wf(isp)       &
!                                   - M_zI * gwf(isp, dir)  &
!                  - M_TWO * vp(dir) / P_c *  wf(isp))
!              end do
!            end do
!          end do
!        end do
!      end do
!   
!    end if   ! this%tdstepsinterval

    POP_SUB(pes_flux_save)
  end subroutine pes_flux_save

  ! ---------------------------------------------------------
  subroutine pes_flux_integrate(this, mesh, st)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st

    integer            :: dim, ist, itstep, dir, isp
    integer            :: start(1:MAX_DIM), end(1:MAX_DIM), startk, endk
    CMPLX, allocatable :: Jk(:)
    FLOAT              :: etime1, etime2

    PUSH_SUB(pes_flux_integrate)

    dim = mesh%sb%dim
    start = 0
    end = 0
    start(1:dim) = this%nsrfcpnt_start(1:dim)
    end(1:dim)   = this%nsrfcpnt_end(1:dim)

    startk = this%nkpnts_start
    endk   = this%nkpnts_end
     
    SAFE_ALLOCATE(Jk(startk:endk))

    ! integrate over time & surface
    do dir = 1, dim
      do isp = start(dir), end(dir)
        etime1 = loct_clock()
        Jk = M_z0
        etime2 = loct_clock()
        do itstep = 1, this%tdsteps
          Jk(startk:endk) = Jk(startk:endk) + this%conjgphase(startk:endk, itstep) * &
            (this%gwf(isp, itstep, dir) - this%kpnt(startk:endk, dir) * this%wf(isp, itstep))
        end do
!        write(*,*) 'time for t-loop:', loct_clock() - etime2
        etime2 = loct_clock()
        Jk(startk:endk) = Jk(startk:endk) * this%conjgplanewf(startk:endk, isp)
        this%spctramp(startk:endk) = this%spctramp(startk:endk) + Jk(startk:endk) * this%srfcnrml(isp, dir)
!        write(*,*) 'time for k-slice:', loct_clock() - etime2
!        write(*,*) 'time for srfcpnt:', loct_clock() - etime1
      end do
    end do

    SAFE_DEALLOCATE_A(Jk)

    POP_SUB(pes_flux_integrate)
  end subroutine pes_flux_integrate

  ! ---------------------------------------------------------
  subroutine pes_flux_getsrfc(this, mesh, border, offset)
    type(mesh_t),     intent(in)    :: mesh
    type(pes_flux_t), intent(inout) :: this
    FLOAT,            intent(in)    :: border(1:MAX_DIM)
    FLOAT,            intent(in)    :: offset(1:MAX_DIM)

    integer, allocatable  :: which_surface(:), aux(:)
    FLOAT                 :: xx(MAX_DIM), dd, dmin
    integer               :: ip, ierr, idim, isp, dir, start, dim

    PUSH_SUB(pes_flux_getsrfc)

    dim = mesh%sb%dim

    SAFE_ALLOCATE(this%nsrfcpnt_start(1:dim))
    SAFE_ALLOCATE(this%nsrfcpnt_end(1:dim))

    SAFE_ALLOCATE(which_surface(1:mesh%np_global))
    which_surface = 0

    SAFE_ALLOCATE(aux(1:dim))
    aux = 0

    this%nsrfcpnts = 0 
    do ip = 1, mesh%np_global
      isp = 0
      xx = mesh_x_global(mesh, ip)
      xx(1:dim) = xx(1:dim) + offset(1:dim)
      do idim = 1, dim
        ! distance to a border
        dd = abs(xx(idim)) - border(idim)
        if(abs(dd) < mesh%spacing(idim)/M_TWO .and. dd <= M_ZERO) then 
          isp = isp + 1 
          which_surface(ip) = int(sign(M_ONE, xx(idim))) * idim     ! +-x=+-1, +-y=+-2
        end if
      end do
      if(isp > 1) then
        ! edges or corners are not counted
        which_surface(ip) = 0
      else if(isp == 1) then
        ! points in absorbing zone are not counted either
        do idim = 1, dim 
          dd = abs(xx(idim)) - border(idim)
          if(dd >= mesh%spacing(idim)) then 
            which_surface(ip) = 0
          end if
        end do
      end if
    end do

    this%nsrfcpnts = 0 
    do ip = 1, mesh%np_global
      if(which_surface(ip) /= 0) this%nsrfcpnts = this%nsrfcpnts + 1
    end do

    SAFE_ALLOCATE(this%srfcpnt(1:this%nsrfcpnts))
    this%srfcpnt = 0

    SAFE_ALLOCATE(this%srfcnrml(1:this%nsrfcpnts, 1:dim))
    this%srfcnrml = M_ZERO

    ! number of surface points with normal vector in +-x, +-y, and +-z direction separately
    do ip = 1, mesh%np_global
      if(which_surface(ip) /= 0) then
        dir = abs(which_surface(ip))
        aux(dir) = aux(dir) + 1
      end if
    end do

    start = 1
    do dir = 1, dim
      this%nsrfcpnt_start(dir) = start
      start = start + aux(dir)
      this%nsrfcpnt_end(dir)   = start - 1
    end do

    ! Fill up this%srfcpnt in correct ordering
    aux(1:dim) = this%nsrfcpnt_start(1:dim) - 1
    do ip = 1, mesh%np_global
      if(which_surface(ip) /= 0) then
        dir = abs(which_surface(ip))
        aux(dir) = aux(dir) + 1

        this%srfcpnt(aux(dir)) = ip
        ! surface normal should point to the inside? Does not make a difference.
        ! add the surface element !!!
        this%srfcnrml(aux(dir), dir) = sign(1, which_surface(ip))
      end if
    end do

!    call dio_function_output(C_OUTPUT_HOW_MESH_INDEX, ".", "surface", mesh, &
!      dble(which_surface), unit_one, ierr)

!    call dio_function_output(C_OUTPUT_HOW_PLANE_Z, ".", "surface", mesh, &
!      dble(which_surface), unit_one, ierr)

    SAFE_DEALLOCATE_A(which_surface)
    SAFE_DEALLOCATE_A(aux)

    write(*,*) 'Info: number of surface points:', this%nsrfcpnts
    write(*,*) 'Info: surface points divisions (start):', this%nsrfcpnt_start(1:dim)
    write(*,*) 'Info: surface points divisions (end):  ', this%nsrfcpnt_end(1:dim)

    SAFE_ALLOCATE(this%mpirank_srfcpnt(1:this%nsrfcpnts))
    do isp = 1, this%nsrfcpnts
      xx = mesh_x_global(mesh, this%srfcpnt(isp))
      ip = mesh_nearest_point(mesh, xx, dmin, this%mpirank_srfcpnt(isp))
    end do

!    write out the field is_on_srfc on plot it!
!    call dio_function_output(512, &
!      ".", "surface", mesh, this%is_on_surface, unit_one, ierr)

!    call dio_function_output(512, &
!      ".", "surface", mesh, dble(which_surface), unit_one, ierr)

    POP_SUB(pes_flux_getsrfc)
  end subroutine pes_flux_getsrfc

  ! ---------------------------------------------------------
  subroutine pes_flux_output(this, sb, dt)
    type(pes_flux_t), intent(inout)    :: this
    type(simul_box_t),   intent(in)    :: sb
    FLOAT,               intent(in)    :: dt

    integer            :: dim, ikp, ikk, iph, iunit, ith
    FLOAT              :: phi, theta, kk
    CMPLX, allocatable :: spctrout(:)

    PUSH_SUB(pes_flux_output)

    dim = sb%dim
    SAFE_ALLOCATE(spctrout(1:this%nkpnts_global))
    spctrout = M_z0

#if defined(HAVE_MPI)
    call MPI_Reduce(this%spctramp(1:this%nkpnts_global), spctrout(1:this%nkpnts_global), &
      this%nkpnts_global, MPI_CMPLX, mpi_sum, 0, mpi_comm_world, mpi_err)
#else
    spctrout(1:this%nkpnts_global) = this%spctramp(1:this%nkpnts_global)
#endif

    if(mpi_grp_is_root(mpi_world)) then
      iunit = io_open('td.general/PESflux_map.z=0', action='write', position='rewind')
      ikp = 0
      do ith = 0, this%ntheta
        theta = ith * this%deltheta + this%thetamin
        do iph = 0, this%nphi
          phi = iph * this%delphi + this%phimin
          do ikk = 1, this%nk
            kk = ikk * this%delk
            ikp = ikp + 1
            write(iunit,'(3f18.10, 1e18.10)') kk, theta, phi, &
              (real(spctrout(ikp))**2 + aimag(spctrout(ikp))**2) * (dt * this%tdstepsinterval)**2
          end do
          write(iunit,'(1x)')
          if(dim == 1) write(iunit,'(1x)')
        end do
      end do
      call io_close(iunit)
    end if

    SAFE_DEALLOCATE_A(spctrout)

    POP_SUB(pes_flux_output)
  end subroutine pes_flux_output

  ! ---------------------------------------------------------
  subroutine pes_flux_dump(restart, this, mesh, ierr)
    type(restart_t),  intent(in)  :: restart
    type(pes_flux_t), intent(in)  :: this
    type(mesh_t),     intent(in)  :: mesh
    integer,          intent(out) :: ierr

    integer    :: err, dim, dumpdim
    CMPLX, allocatable :: phasedump(:,:), spctrdump(:)

    PUSH_SUB(pes_flux_dump)

    dim = mesh%sb%dim
    dumpdim = this%nkpnts_global * (this%tdsteps + 1)

    if(restart_skip(restart)) then
      POP_SUB(pes_flux_dump)
      return
    end if

    SAFE_ALLOCATE(spctrdump(1:this%nkpnts_global))
    SAFE_ALLOCATE(phasedump(1:this%nkpnts_global, 0:this%tdsteps))

    if(in_debug_mode) then
      message(1) = "Debug: Writing pes_flux restart."
      call messages_info(1)
    end if

#if defined(HAVE_MPI)
    call MPI_Reduce(this%spctramp(1:this%nkpnts_global), spctrdump(1:this%nkpnts_global), &
      this%nkpnts_global, MPI_CMPLX, mpi_sum, 0, mpi_comm_world, mpi_err)
    call MPI_Reduce(this%conjgphase(1:this%nkpnts_global,0:this%tdsteps), phasedump(1:this%nkpnts_global,0:this%tdsteps), &
      dumpdim, MPI_CMPLX, mpi_sum, 0, mpi_comm_world, mpi_err)
#else
    spctrdump = this%spctramp
    phasedump = this%conjgphase
#endif

    call zrestart_write_binary(restart, 'pesflux1', this%nkpnts_global, spctrdump, err)
    call zrestart_write_binary(restart, 'pesflux2', dumpdim, phasedump, err)
    call zrestart_write_binary(restart, 'pesflux3', this%nsrfcpnts*this%tdsteps, this%wf, err)
    call zrestart_write_binary(restart, 'pesflux4', this%nsrfcpnts*this%tdsteps*dim, this%gwf, err)

    if(err /= 0) ierr = ierr + 1

    if(in_debug_mode) then
      message(1) = "Debug: Writing pes_flux restart done."
      call messages_info(1)
    end if

    SAFE_DEALLOCATE_A(phasedump)
    SAFE_DEALLOCATE_A(spctrdump)

    POP_SUB(pes_flux_dump)
  end subroutine pes_flux_dump

  ! ---------------------------------------------------------
  subroutine pes_flux_load(restart, this, mesh, ierr)
    type(restart_t),     intent(in)    :: restart
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    integer,             intent(out)   :: ierr

    integer   :: err, dim, loaddim
    CMPLX, allocatable :: phaseload(:,:), spctrload(:)

    PUSH_SUB(pes_flux_load)

    dim = mesh%sb%dim
    loaddim = this%nkpnts_global * (this%tdsteps + 1)
    ierr = 0

    if(restart_skip(restart)) then
      ierr = -1
      POP_SUB(pes_flux_load)
      return
    end if

    SAFE_ALLOCATE(spctrload(1:this%nkpnts_global))
    SAFE_ALLOCATE(phaseload(1:this%nkpnts_global, 0:this%tdsteps))

    if(in_debug_mode) then
      message(1) = "Debug: Reading pes_flux restart."
      call messages_info(1)
    end if

    call zrestart_read_binary(restart, 'pesflux1', this%nkpnts_global, spctrload, err)
    call zrestart_read_binary(restart, 'pesflux2', loaddim, phaseload, err)
    call zrestart_read_binary(restart, 'pesflux3', this%nsrfcpnts*this%tdsteps, this%wf, err)
    call zrestart_read_binary(restart, 'pesflux4', this%nsrfcpnts*this%tdsteps*dim, this%gwf, err)

    this%spctramp = M_z0
    this%spctramp(this%nkpnts_start:this%nkpnts_end) = spctrload(this%nkpnts_start:this%nkpnts_end)

    this%conjgphase = M_z0
    this%conjgphase(this%nkpnts_start:this%nkpnts_end,:) = phaseload(this%nkpnts_start:this%nkpnts_end,:)

    if(err /= 0) ierr = ierr + 1
   
    if(in_debug_mode) then
      message(1) = "Debug: Reading pes_flux restart done."
      call messages_info(1)
    end if

    SAFE_DEALLOCATE_A(spctrload)
    SAFE_DEALLOCATE_A(phaseload)

    POP_SUB(pes_flux_load)
  end subroutine pes_flux_load

end module pes_flux_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
