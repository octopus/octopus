#include "global.h"

module pes_flux_m
  use global_m
  use messages_m
!  use datasets_m
  use mesh_m
  use profiling_m
  use varinfo_m
  use parser_m
  use states_m
  use grid_m
  use derivatives_m
  use hamiltonian_m
  use lasers_m
  use io_function_m
  use io_m
  use unit_m
  use unit_system_m
  use simul_box_m

  implicit none

  private

  public ::                    &
    pes_flux_t,                &
    pes_flux_init,             &
    pes_flux_end,              &
    pes_flux_calc,             &
    pes_flux_init_boundaries,  &
    pes_flux_apply_boundaries

  type pes_flux_t
    integer           :: nkpnts         ! total number of k-points
    integer           :: nk, nphi
    FLOAT             :: delk, phimin, delphi
    integer           :: srfcshape      ! shape of the surface (= cubic/spherical)
    integer           :: nsrfcpnts      ! total number of points contructing surface
    integer           :: abmethod       ! method for absorbing boundaries
    integer           :: output
    integer           :: interval

    FLOAT, pointer    :: kpnt(:,:)      ! coordinates of all k-points
    integer, pointer  :: srfcpnt(:)     ! returns the index of the points on the surface
    FLOAT, pointer    :: srfcnrml(:,:)  ! (unit) vectors normal to the surface
    FLOAT, pointer    :: vlkvphase(:)   ! current Volkov phase for all k-points
    CMPLX, pointer    :: spctramp(:,:,:,:) ! spectral amplitude
  end type pes_flux_t

  integer, parameter ::   &
    M_CUBIC      = 1,     &
    M_SPHERICAL  = 2

contains


  ! ---------------------------------------------------------
  subroutine pes_flux_init(flux, mesh, st, hm)
    type(pes_flux_t),    intent(inout) :: flux
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(in)    :: hm

    type(block_t)        :: blk
    FLOAT, allocatable   :: border(:)       ! distance of surface from border
    integer              :: id, ikp, isp, il, iph, ikk
    FLOAT                :: phi, kk
!    type(mesh_t)         :: kmesh
    FLOAT                :: phimax, kmax

    PUSH_SUB(pes_flux_init)

    do il = 1, hm%ep%no_lasers
      if(laser_kind(hm%ep%lasers(il)) /= E_FIELD_VECTOR_POTENTIAL) then
        message(1) = 'Surff only works in velocity gauge.'
        call messages_fatal(1)
      end if
    end do

    if(mesh%sb%dim > 1 .and. mesh%sb%box_shape /= PARALLELEPIPED) then
      message(1) = 'Surff only works with BoxShape = parallelepiped in dim > 1.'
      call messages_fatal(1)
    end if

    message(1) = 'Info: Calculation PES using t-surff technique.'
    call messages_info(1)
  ! read in all the needed input data
  ! - surface points: either cubic or sphere (first cubic; sphere needs interp.)
  ! - k points (automatically at first, later manually from block)
  ! - 

!    if(parse_block(datasets_check('PhotoElectronSpectrumPoints'), blk) < 0) then
!      message(1) = 'The PhotoElectronSpectrumPoints block is required when PhotoElectronSpectrum = pes_flux'
!      call messages_fatal(1)
!    end if

!    flux%nsrfcpnts = parse_block_n(blk)
!    if(flux%nsrfcpnts /= 2) then
!      message(1) = 'Surff only works with 2 surface points.'
!      call messages_fatal(1)
!    end if


    !%Variable SurffShape
    !%Type integer
    !%Default m_cubic
    !%Section Time-Dependent::Surff
    !%Description
    !% The surface shape.
    !%Option m_cubic 1
    !% cubic surface.
    !%End
!    call parse_integer(datasets_check('SurffShape'), M_CUBIC, flux%srfcshape)
!    if(.not.varinfo_valid_option('SurffShape', flux%srfcshape)) &
!      call input_error('SurffShape')
!    call messages_print_var_option(stdout, "SurffShape", flux%srfcshape)

    ! output
    call parse_integer('PESSurfaceOutput', 1, flux%output)

    ! surface
    SAFE_ALLOCATE(border(1:MAX_DIM))
    if(parse_block('PESSurface', blk) < 0) then
      border(:) = hm%ab_width
      message(1) = "PESSurface not specified. Using default values."
      call messages_info(1)
    else
      call parse_block_float(blk, 0, 0, border(1))
      call parse_block_float(blk, 0, 1, border(2))
      call parse_block_float(blk, 0, 2, border(3))
    end if
    call pes_flux_getsrfc(flux, mesh, border)
    SAFE_DEALLOCATE_A(border)

    ! k-mesh in 1D (2 points) & 2D (polar coordinates)
    call parse_float('PESSurfaceKmax', M_ONE, kmax)
    call parse_float('PESSurfaceDeltaK', CNST(0.002), flux%delk)
    call parse_float('PESSurfacePhiMin', M_ZERO, flux%phimin)
    call parse_float('PESSurfacePhiMax', M_TWO * M_PI, phimax)
    call parse_float('PESSurfaceDelPhi', CNST((M_TWO * M_PI)/360), flux%delphi)

    if(mesh%sb%dim == 1) then
      phimax = M_PI
      flux%phimin = M_ZERO
      flux%delphi = M_PI
    end if
    flux%nphi   = nint((phimax - flux%phimin)/flux%delphi)
    flux%nk     = nint(kmax/flux%delk)
    flux%nkpnts = (flux%nphi + 1) * flux%nk 

!    write(*,*) 'testb01'
    write(*,*) flux%nphi, flux%nk, flux%nkpnts, flux%nsrfcpnts

    ! k-points
    SAFE_ALLOCATE(flux%kpnt(1:flux%nkpnts, 1:mesh%sb%dim))
    flux%kpnt = M_ZERO

    ikp = 0
    do iph = 0, flux%nphi
      phi = iph * flux%delphi + flux%phimin
      do ikk = 1, flux%nk
        kk = ikk * flux%delk
        ikp = ikp + 1
                             flux%kpnt(ikp, 1) = kk * cos(phi)
        if(mesh%sb%dim == 2) flux%kpnt(ikp, 2) = kk * sin(phi)
      end do
    end do

!    write(*,*) 'testb02'

    ! other stuff
    call parse_integer('PESSurfaceInterval', 1, flux%interval)

    SAFE_ALLOCATE(flux%vlkvphase(1:flux%nkpnts))
    flux%vlkvphase(:) = M_ZERO

    SAFE_ALLOCATE(flux%spctramp(1:flux%nkpnts, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    flux%spctramp = M_z0

    write(*,*) 'end init.'

    POP_SUB(pes_flux_init)
  end subroutine pes_flux_init

  ! ---------------------------------------------------------
  subroutine pes_flux_end(flux)
    type(pes_flux_t), intent(inout) :: flux

    PUSH_SUB(pes_flux_end)

    SAFE_DEALLOCATE_P(flux%kpnt)
    SAFE_DEALLOCATE_P(flux%spctramp)

    SAFE_DEALLOCATE_P(flux%srfcpnt)
    SAFE_DEALLOCATE_P(flux%srfcnrml)

    SAFE_DEALLOCATE_P(flux%vlkvphase)

    POP_SUB(pes_flux_end)
  end subroutine pes_flux_end


  ! ---------------------------------------------------------
  subroutine pes_flux_init_boundaries()

  ! generate the boundaries conditions
  ! at the beginning only mask function supported
  ! should call an external module for boundary conditions

  end subroutine pes_flux_init_boundaries

  ! ---------------------------------------------------------
  subroutine pes_flux_apply_boundaries()

  ! apply the chosen absorption method every time step

  end subroutine pes_flux_apply_boundaries

  ! ---------------------------------------------------------
  subroutine pes_flux_volkov(flux)
    type(pes_flux_t),  intent(inout) :: flux


  
  ! - generate the Volkov state at given r-point for given k-point and 
  !   given vector potential (at time t):
  !   One probably just needs to apply a new factor like in 
  !   "pes_mask_volkov_time_evolution_wf"
  ! - do not generate the gradient of the Volkov state (this is just a 
  !   multiplication with i*k)

  end subroutine pes_flux_volkov

  ! ---------------------------------------------------------
  subroutine pes_flux_volkovphase


  end subroutine pes_flux_volkovphase

  ! ---------------------------------------------------------
  subroutine pes_flux_derpsi()

    ! calculate the gradient of psi at a given r-point on surface
    ! for sphere probably interpolation necessary


  end subroutine pes_flux_derpsi
  ! ---------------------------------------------------------
  subroutine pes_flux_calc(flux, mesh, st, gr, hm, iter, dt)
    type(pes_flux_t),    intent(inout) :: flux
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt

    integer            :: ik, ist, idim, ikp, isp, rankmin, il, iph, ikk
    FLOAT              :: dmin
    CMPLX, allocatable :: gzpsi(:,:)
    CMPLX, allocatable :: fluxx(:)
    FLOAT, allocatable :: summ(:)
    FLOAT, allocatable :: vp(:)
    CMPLX, allocatable :: wf(:,:,:,:), gwf(:,:,:,:,:)
    FLOAT              :: phi, kk
    FLOAT              :: krr, vec
    CMPLX              :: planewf   ! plane waves for each k-point on the surface

    PUSH_SUB(pes_flux_calc)

    SAFE_ALLOCATE(wf(1:flux%nsrfcpnts, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    wf = M_z0

    SAFE_ALLOCATE(gwf(1:flux%nsrfcpnts, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik, 1:mesh%sb%dim))
    gwf = M_z0

    SAFE_ALLOCATE(gzpsi(1:mesh%np, 1:mesh%sb%dim))
    gzpsi = M_z0

    SAFE_ALLOCATE(fluxx(1:mesh%sb%dim))
    fluxx = M_z0

    SAFE_ALLOCATE(summ(1:flux%nkpnts))
    summ = M_ZERO

    SAFE_ALLOCATE(vp(1:mesh%sb%dim))
    vp = M_ZERO

    ! calculate the vector potential
    do il = 1, hm%ep%no_lasers
      ! add current fields
      call laser_field(hm%ep%lasers(il), vp, iter*dt)
    end do

    ! update the Volkov phase adding the previous time steps
    do ikp = 1, flux%nkpnts
      vec = sum((flux%kpnt(ikp, :) - vp(:) / P_c)**2)
      flux%vlkvphase(ikp) = flux%vlkvphase(ikp) + vec * dt / M_TWO
    end do

    ! accumulate flux through surface every flux%interval (save time)
    if(flux%interval > 0 .and. mod(iter, flux%interval) == 0) then
!      write(*,*) 'test01'

      ! get states and derivatives
      ! change this to states_get_state...
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
!            write(*,*) 'before grad.'
            call zderivatives_grad(gr%der, st%zpsi(:, idim, ist, ik), gzpsi)
!            write(*,*) 'after grad.'
            do isp = 1, flux%nsrfcpnts
              ! here one should call an interpolation scheme, now
              ! we just take the wave function at the nearest grid point
              wf(isp, idim, ist, ik) = st%occ(ist, ik) * st%zpsi(flux%srfcpnt(isp), idim, ist, ik)
              gwf(isp, idim, ist, ik, :) = st%occ(ist, ik) * gzpsi(flux%srfcpnt(isp), :)
            end do
          end do
        end do
      end do
   
!      write(*,*) 'test02'
   
      ! accumulate flux through surface
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            do ikp = 1, flux%nkpnts
              do isp = 1, flux%nsrfcpnts
                krr = dot_product(flux%kpnt(ikp, :), mesh%x(flux%srfcpnt(isp), :))
                planewf = exp(M_zI * krr) / (M_TWO * M_PI)**M_TWOTHIRD
   
                ! conjg() or not?
                fluxx(:) = conjg(planewf * exp(-M_zI * flux%vlkvphase(ikp))) *          &
                     (flux%kpnt(ikp, :) *  wf(isp, idim, ist, ik)                       &
                                 - M_zI * gwf(isp, idim, ist, ik, :)                    &
                  - M_TWO * vp(:) / P_c *  wf(isp, idim, ist, ik))
                flux%spctramp(ikp, idim, ist, ik) = flux%spctramp(ikp, idim, ist, ik) + &
                  dot_product(fluxx(:), flux%srfcnrml(isp, :)) * dt * flux%interval / M_TWO
                ! times the surface element (integrate that into srfcnrml)
                ! Sign of planephase? Minus, obtained in tests in weak field regime. 
                ! Sign of vlkvphase ???
              end do
            end do
          end do
        end do
      end do
   
!      write(*,*) 'test03'
   
      ! sum over all states, spins, etc.  & write out (put this in pes_flux_output & add mpi_communication)
      do ikp = 1, flux%nkpnts
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end
            do idim = 1, st%d%dim
              summ(ikp) = summ(ikp) + real(flux%spctramp(ikp, idim, ist, ik))**2 + aimag(flux%spctramp(ikp, idim, ist, ik))**2
            end do
          end do
        end do
      end do

      if(mod(iter, flux%output) == 0) then
        open(unit=300,position='rewind')
        ikp = 0
        do iph = 0, flux%nphi
          phi = iph * flux%delphi + flux%phimin
          do ikk = 1, flux%nk
            kk = ikk * flux%delk
            ikp = ikp + 1
            write(300,'(2f18.10, 1e18.10)') kk, phi, summ(ikp)
          end do
          write(300,*)
          if(mesh%sb%dim == 1) write(300,*)
        end do
        flush(300)
        close(300)
      end if

    end if   ! flux%interval

!    write(*,*) 'test04'

    SAFE_DEALLOCATE_A(wf)
    SAFE_DEALLOCATE_A(gwf)
    SAFE_DEALLOCATE_A(gzpsi)
    SAFE_DEALLOCATE_A(fluxx)
    SAFE_DEALLOCATE_A(summ)
    SAFE_DEALLOCATE_A(vp)

!    write(*,*) 'end flux calc.'

    POP_SUB(pes_flux_calc)
  end subroutine pes_flux_calc
  ! ---------------------------------------------------------
  subroutine pes_flux_grad(flux)
    type(pes_flux_t), intent(inout) :: flux


  end subroutine pes_flux_grad

  ! ---------------------------------------------------------
  subroutine pes_flux_getsrfc(flux, mesh, border)
    type(mesh_t),     intent(in)    :: mesh
    type(pes_flux_t), intent(inout) :: flux
    FLOAT,            intent(in)    :: border(1:MAX_DIM)

    integer, allocatable  :: which_surface(:)
    FLOAT                 :: xx(MAX_DIM), rr, dd
    integer               :: ip, ierr, idim, isp

    PUSH_SUB(pes_flux_getsrfc)

    SAFE_ALLOCATE(which_surface(1:mesh%np))
    which_surface(:) = 0

    flux%nsrfcpnts = 0 
    do ip = 1, mesh%np
      isp = 0
      call mesh_r(mesh, ip, rr, coords=xx)
      do idim = 1, mesh%sb%dim
        ! distance to a border
        dd = abs(xx(idim)) - (mesh%sb%lsize(idim) - border(idim))
        if(abs(dd) < mesh%spacing(idim) .and. dd >= M_ZERO) then 
          isp = isp + 1 
          which_surface(ip) = int(sign(M_ONE, xx(idim))) * idim     ! +-x=+-1, +-y=+-2
          flux%nsrfcpnts = flux%nsrfcpnts + 1
        end if
      end do
      if(isp > 1) then
        which_surface(ip) = 0                                ! corners are not counted
        flux%nsrfcpnts = flux%nsrfcpnts - isp
      else if(isp == 1) then
        do idim = 1, mesh%sb%dim
          dd = abs(xx(idim)) - (mesh%sb%lsize(idim) - border(idim))
          if(dd >= mesh%spacing(idim)) then 
            which_surface(ip) = 0
            flux%nsrfcpnts = flux%nsrfcpnts - 1
          end if
        end do
      end if
      ! if(isp == M_ONE) flux%is_on_surface(ip) = M_ONE      ! corners are not counted
      ! flux%is_on_surface(ip) = isp                         ! corners are counted 2x
    end do

    SAFE_ALLOCATE(flux%srfcpnt(1:flux%nsrfcpnts))
    flux%srfcpnt(:) = 0

    SAFE_ALLOCATE(flux%srfcnrml(1:flux%nsrfcpnts, 1:mesh%sb%dim))
    flux%srfcnrml = M_ZERO

    isp = 0 
    do ip = 1, mesh%np
      if(which_surface(ip) /= 0) then
        isp = isp + 1

        flux%srfcpnt(isp)  = ip
        ! surface normal should point to the inside? Does not make a difference.
        ! add the surface element !!!
        flux%srfcnrml(isp, abs(which_surface(ip))) = sign(1, which_surface(ip))
      end if
    end do

    SAFE_DEALLOCATE_A(which_surface)

    write(*,*) 'Surface points:', flux%nsrfcpnts
    do isp = 1, flux%nsrfcpnts
      write(*,*) isp, mesh%x(flux%srfcpnt(isp),:)
    end do

!    write out the field is_on_srfc on plot it!
!    call dio_function_output(512, &
!      ".", "surface", mesh, flux%is_on_surface, unit_one, ierr)

!    call dio_function_output(512, &
!      ".", "surface", mesh, dble(which_surface), unit_one, ierr)

    POP_SUB(pes_flux_getsrfc)
  end subroutine pes_flux_getsrfc

end module pes_flux_m
