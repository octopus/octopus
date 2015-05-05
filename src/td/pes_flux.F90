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
  use parser_m
  use states_m
  use grid_m
  use derivatives_m
  use hamiltonian_m
  use lasers_m
  use io_m
  use simul_box_m

  implicit none

  private

  public ::                    &
    pes_flux_t,                &
    pes_flux_init,             &
    pes_flux_end,              &
    pes_flux_calc,             &
    pes_flux_output

  type pes_flux_t
    integer           :: nkpnts                !< total number of k-points
    integer           :: nk, nphi
    FLOAT             :: delk, phimin, delphi
    integer           :: srfcshape             !< shape of the surface (= cubic/spherical)
    integer           :: nsrfcpnts             !< total number of points contructing surface
    integer           :: tdstepsinterval

    integer, pointer  :: nsrfcpnt_start(:)
    integer, pointer  :: nsrfcpnt_end(:)
    FLOAT, pointer    :: kpnt(:,:)             !< coordinates of all k-points
    integer, pointer  :: srfcpnt(:)            !< returns the index of the points on the surface
    FLOAT, pointer    :: srfcnrml(:,:)         !< (unit) vectors normal to the surface
    FLOAT, pointer    :: vlkvphase(:)          !< current Volkov phase for all k-points
    CMPLX, pointer    :: spctramp(:,:,:,:)     !< spectral amplitude
    CMPLX, pointer    :: Jk(:,:,:,:,:,:)       !< current density operator   
    CMPLX, pointer    :: phik(:,:)             !< Volkov waves
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
    integer              :: id, ikp, isp, il, iph, ikk, dim, start, end, nn
    FLOAT                :: phi, kk, krr
    FLOAT                :: phimax, kmax

    PUSH_SUB(pes_flux_init)

    dim = mesh%sb%dim

    call messages_experimental("PhotoElectronSpectrum with t-surff")

    do il = 1, hm%ep%no_lasers
      if(laser_kind(hm%ep%lasers(il)) /= E_FIELD_VECTOR_POTENTIAL) then
        message(1) = 't-surff only works in velocity gauge.'
        call messages_fatal(1)
      end if
    end do

    if(dim > 1 .and. mesh%sb%box_shape /= PARALLELEPIPED) then
      message(1) = 'Surff only works with BoxShape = parallelepiped in dim > 1.'
      call messages_fatal(1)
    end if

    message(1) = 'Info: Calculation PES using t-surff technique.'
    call messages_info(1)

    ! surface
    if(parse_block('PESSurface', blk) < 0) then
      border(:) = hm%ab_width
      message(1) = "PESSurface not specified. Using default values."
      call messages_info(1)
    else
      ASSERT(dim == 3)
      call parse_block_float(blk, 0, 0, border(1))
      call parse_block_float(blk, 0, 1, border(2))
      call parse_block_float(blk, 0, 2, border(3))
    end if
    call pes_flux_getsrfc(this, mesh, border)

    ! k-mesh in 1D (2 points) & 2D (polar coordinates)
    call parse_variable('PESSurfaceKmax', M_ONE, kmax)
    call parse_variable('PESSurfaceDeltaK', CNST(0.002), this%delk)
    call parse_variable('PESSurfacePhiMin', M_ZERO, this%phimin)
    call parse_variable('PESSurfacePhiMax', M_TWO * M_PI, phimax)
    call parse_variable('PESSurfaceDelPhi', CNST((M_TWO * M_PI)/360), this%delphi)

    if(dim == 1) then
      phimax = M_PI
      this%phimin = M_ZERO
      this%delphi = M_PI
    end if
    this%nphi   = nint((phimax - this%phimin)/this%delphi)
    this%nk     = nint(kmax/this%delk)
    this%nkpnts = (this%nphi + 1) * this%nk 

    ! k-points
    SAFE_ALLOCATE(this%kpnt(1:this%nkpnts, 1:dim))
    this%kpnt = M_ZERO

    ikp = 0
    do iph = 0, this%nphi
      phi = iph * this%delphi + this%phimin
      do ikk = 1, this%nk
        kk = ikk * this%delk
        ikp = ikp + 1
                             this%kpnt(ikp, 1) = kk * cos(phi)
        if(dim == 2) this%kpnt(ikp, 2) = kk * sin(phi)
      end do
    end do

    ! other stuff
    call parse_variable('PESSurfaceTDStepsInterval', 1, this%tdstepsinterval)

    SAFE_ALLOCATE(this%vlkvphase(1:this%nkpnts))
    this%vlkvphase(:) = M_ZERO

    start = st%st_start
    end = st%st_end
    SAFE_ALLOCATE(this%spctramp(1:this%nkpnts, 1:st%d%dim, start:end, 1:st%d%nik))
    this%spctramp = M_z0

    nn = this%nsrfcpnts
    SAFE_ALLOCATE(this%Jk(1:this%nkpnts, 1:st%d%dim, start:end, 1:st%d%nik, 1:nn, 1:dim))
    this%Jk = M_z0

    SAFE_ALLOCATE(this%phik(1:this%nkpnts, 1:this%nsrfcpnts))
    this%phik = M_z0 
    
    do ikp = 1, this%nkpnts
      do isp = 1, this%nsrfcpnts
        krr = dot_product(this%kpnt(ikp,1:dim), mesh%x(this%srfcpnt(isp),1:dim))
        this%phik(ikp,isp) = exp(M_zI * krr) / (M_TWO * M_PI)**(dim/M_TWO) 
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
    SAFE_DEALLOCATE_P(this%Jk)
    SAFE_DEALLOCATE_P(this%phik)
    
    SAFE_DEALLOCATE_P(this%srfcpnt)
    SAFE_DEALLOCATE_P(this%srfcnrml)

    SAFE_DEALLOCATE_P(this%nsrfcpnt_start)
    SAFE_DEALLOCATE_P(this%nsrfcpnt_end)

    SAFE_DEALLOCATE_P(this%vlkvphase)

    POP_SUB(pes_flux_end)
  end subroutine pes_flux_end

  ! ---------------------------------------------------------
  subroutine pes_flux_calc(this, mesh, st, gr, hm, iter, dt)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt

    integer            :: ik, ist, idim, ikp, isp, rankmin, il, iph, ikk, dim,nsrf,nkp, dir
    FLOAT              :: dmin
    CMPLX, allocatable :: gpsi(:,:), psi(:)
    CMPLX, allocatable :: fluxx(:)
    FLOAT, allocatable :: vp(:)
    CMPLX, allocatable :: wf(:), gwf(:,:)
    FLOAT              :: phi, kk
    integer            :: start(1:MAX_DIM), end(1:MAX_DIM)
    FLOAT              :: krr, vec
    CMPLX              :: planewf   ! plane waves for each k-point on the surface

    PUSH_SUB(pes_flux_calc)

    if(iter > 0 .and. mod(iter, this%tdstepsinterval) == 0) then

      dim  = mesh%sb%dim
      nsrf = this%nsrfcpnts
      nkp  = this%nkpnts

      start(1:dim) = this%nsrfcpnt_start(1:dim)
      end(1:dim)   = this%nsrfcpnt_end(1:dim)

      SAFE_ALLOCATE(wf(1:this%nsrfcpnts))
      wf = M_z0

      SAFE_ALLOCATE(gwf(1:this%nsrfcpnts, 1:dim))
      gwf = M_z0

      SAFE_ALLOCATE(psi(1:mesh%np_part))
      psi = M_z0

      SAFE_ALLOCATE(gpsi(1:mesh%np, 1:dim))
      gpsi = M_z0

      SAFE_ALLOCATE(fluxx(1:dim))
      fluxx = M_z0

      SAFE_ALLOCATE(vp(1:dim))
      vp = M_ZERO

      ! calculate the vector potential
      do il = 1, hm%ep%no_lasers
        ! add current fields
        call laser_field(hm%ep%lasers(il), vp, iter*dt)
      end do

      ! update the Volkov phase adding the previous time steps
      do ikp = 1, this%nkpnts
        vec = sum((this%kpnt(ikp, 1:dim) - vp(1:dim) / P_c)**2)
        this%phik(ikp,:) = this%phik(ikp,:)* exp(-M_zI* vec * dt / M_TWO) 
      end do

      ! accumulate flux through surface
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            call states_get_state(st, mesh, idim, ist, ik, psi)
            call zderivatives_grad(gr%der, psi, gpsi)
            do isp = 1, this%nsrfcpnts
               wf(isp) = st%occ(ist, ik) * psi(this%srfcpnt(isp))
              gwf(isp, 1:dim) = st%occ(ist, ik) * gpsi(this%srfcpnt(isp), 1:dim) 
            end do
            do dir = 1, dim
              do isp = start(dir), end(dir)
                this%Jk(1:this%nkpnts, idim, ist, ik, isp, dir) =    &
                  this%Jk(1:this%nkpnts, idim, ist, ik, isp, dir) +  &
                  conjg(this%phik(1:this%nkpnts, isp)) *             & 
                  (    this%kpnt(1:this%nkpnts, dir) * wf(isp)       &
                                   - M_zI * gwf(isp, dir)  &
                  - M_TWO * vp(dir) / P_c *  wf(isp))
              end do
            end do
          end do
        end do
      end do
   
    end if   ! this%tdstepsinterval

    SAFE_DEALLOCATE_A(wf)
    SAFE_DEALLOCATE_A(gwf)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(gpsi)
    SAFE_DEALLOCATE_A(fluxx)
    SAFE_DEALLOCATE_A(vp)

    POP_SUB(pes_flux_calc)
  end subroutine pes_flux_calc

  ! ---------------------------------------------------------
  subroutine pes_flux_getsrfc(this, mesh, border)
    type(mesh_t),     intent(in)    :: mesh
    type(pes_flux_t), intent(inout) :: this
    FLOAT,            intent(in)    :: border(1:MAX_DIM)

    integer, allocatable  :: which_surface(:), aux(:)
    FLOAT                 :: xx(MAX_DIM), rr, dd
    integer               :: ip, ierr, idim, isp, dir, start, dim

    PUSH_SUB(pes_flux_getsrfc)

    dim = mesh%sb%dim

    SAFE_ALLOCATE(this%nsrfcpnt_start(1:dim))
    SAFE_ALLOCATE(this%nsrfcpnt_end(1:dim))

    SAFE_ALLOCATE(which_surface(1:mesh%np))
    which_surface = 0

    SAFE_ALLOCATE(aux(1:dim))
    aux = 0

    this%nsrfcpnts = 0 
    do ip = 1, mesh%np
      isp = 0
      call mesh_r(mesh, ip, rr, coords=xx)
      do idim = 1, dim
        ! distance to a border
        dd = abs(xx(idim)) - (mesh%sb%lsize(idim) - border(idim))
        if(abs(dd) < mesh%spacing(idim) .and. dd >= M_ZERO) then 
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
          dd = abs(xx(idim)) - (mesh%sb%lsize(idim) - border(idim))
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
    do ip = 1, mesh%np
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
    do ip = 1, mesh%np
      if(which_surface(ip) /= 0) then
        dir = abs(which_surface(ip))
        aux(dir) = aux(dir) + 1

        this%srfcpnt(aux(dir)) = ip
        ! surface normal should point to the inside? Does not make a difference.
        ! add the surface element !!!
        this%srfcnrml(aux(dir), dir) = sign(1, which_surface(ip))
      end if
    end do

    SAFE_DEALLOCATE_A(which_surface)
    SAFE_DEALLOCATE_A(aux)

    write(*,*) 'Info: number of surface points:', this%nsrfcpnts
    do isp = 1, this%nsrfcpnts
      write(*,*) isp, mesh%x(this%srfcpnt(isp),:)
    end do

!    write out the field is_on_srfc on plot it!
!    call dio_function_output(512, &
!      ".", "surface", mesh, this%is_on_surface, unit_one, ierr)

!    call dio_function_output(512, &
!      ".", "surface", mesh, dble(which_surface), unit_one, ierr)

    POP_SUB(pes_flux_getsrfc)
  end subroutine pes_flux_getsrfc

  ! ---------------------------------------------------------
  subroutine pes_flux_output(this, mesh, st)
    type(pes_flux_t), intent(inout)    :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(in)    :: st

    integer            :: ikp, ikk, iph, iunit, idim, ik, ist, isp
    FLOAT              :: phi, kk
    FLOAT, allocatable :: summ(:)

    PUSH_SUB(pes_flux_output)

    SAFE_ALLOCATE(summ(1:this%nkpnts))
    summ = M_ZERO

    this%spctramp(:, :, :, :) = M_z0
    do ikp = 1, this%nkpnts
      do ik = st%d%kpt%start, st%d%kpt%end
        do idim = 1, st%d%dim
          do ist = st%st_start, st%st_end
            do isp = 1, this%nsrfcpnts
              this%spctramp(ikp, idim, ist, ik) = this%spctramp(ikp, idim, ist, ik) + &
                               dot_product(this%Jk(ikp, idim, ist, ik, isp,:), this%srfcnrml(isp,:)) 
  !                              &
  !                              * dt * this%interval / M_TWO
            end do
          end do
        end do
      end do
    end do
    

    ! sum over all states, spins, etc.  & write out (put this in pes_flux_output & add mpi_communication)
    do ikp = 1, this%nkpnts
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            summ(ikp) = summ(ikp) + &
              real(this%spctramp(ikp, idim, ist, ik))**2 + aimag(this%spctramp(ikp, idim, ist, ik))**2
          end do
        end do
      end do
    end do

    iunit = io_open('td.general/PESflux_map.z=0', action='write', position='rewind')

    ikp = 0
    do iph = 0, this%nphi
      phi = iph * this%delphi + this%phimin
      do ikk = 1, this%nk
        kk = ikk * this%delk
        ikp = ikp + 1
        write(iunit,'(2f18.10, 1e18.10)') kk, phi, summ(ikp)
      end do
      write(iunit,'(1x)')
      if(mesh%sb%dim == 1) write(iunit,'(1x)')
    end do

    call io_close(iunit)
    SAFE_DEALLOCATE_A(summ)

    POP_SUB(pes_flux_output)
  end subroutine pes_flux_output

end module pes_flux_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
