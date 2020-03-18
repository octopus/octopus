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

#include "global.h"

module pes_flux_oct_m
  use boundary_op_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use kpoints_oct_m
  use io_oct_m
  use lasers_oct_m
  use loct_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_interpolation_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private

  public ::                    &
    pes_flux_t,                &
    pes_flux_init,             &
    pes_flux_end,              &
    pes_flux_calc,             &
    pes_flux_output,           &
    pes_flux_load,             &
    pes_flux_dump,             &
    pes_flux_reciprocal_mesh_gen, &
    pes_flux_pmesh,            &
    pes_flux_map_from_states,  &
    pes_flux_out_energy

  type pes_flux_t
    private
    !< NOTE: On unfortunate choice of nomenclature. In this module we use the variable k to indicate 
    !< the momentum of photoelectrons. This has to be distinguished from the electron pseudomomentum 
    !< in periodic systems and its discretized values (AKA kpoints). 
    
    integer          :: nkpnts                         !< total number of k-points
    integer          :: nkpnts_start, nkpnts_end       !< start/end of index for k-points on the current node
    integer          :: nk
    integer          :: nk_start, nk_end
    integer          :: nstepsthetak, nstepsphik       !< parameters for k-mesh
    integer          :: nstepsomegak
    integer          :: nstepsomegak_start, nstepsomegak_end
    FLOAT            :: dk                             !< parameters for k-mesh
    FLOAT, pointer   :: krad(:)                        !< the grid for the radial part of the k-mesh
    FLOAT, pointer   :: kcoords_cub(:,:,:)             !< coordinates of k-points
    FLOAT, pointer   :: kcoords_sph(:,:,:)
    integer          :: kgrid                          !< how is the grid in k: polar/cartesian

    integer, public  :: surf_shape                     !< shape of the surface (= cube/sphere/planes)
    integer          :: nsrfcpnts                      !< total number of surface points
    integer          :: nsrfcpnts_start, nsrfcpnts_end !< for cubic surface: number of surface points on node
    FLOAT, pointer   :: srfcnrml(:,:)                  !< vectors normal to the surface (includes surface element)
    FLOAT, pointer   :: rcoords(:,:)                   !< coordinates of the surface points
    integer, pointer :: srfcpnt(:)                     !< for cubic surface: returns local index of the surface points
    integer, pointer :: rankmin(:)                     !< for cubic surface: returns node which has the surface point
    integer          :: lmax                           !< for spherical surface
    CMPLX, pointer   :: ylm_r(:,:,:)                   !< for spherical surface
    CMPLX, pointer   :: ylm_k(:,:,:)                   !< for spherical surface
    FLOAT, pointer   :: j_l(:,:)                       !< for spherical surface
    FLOAT            :: radius
    integer, pointer :: face_idx_range(:,:)            !< face_idx_start(nface,1:2) 1 (2) start (end) idx of face nface in rcoords(:,:)


    CMPLX, pointer   :: expg(:,:,:)                    !< Fourier basis on a face of the cube
    FLOAT, pointer   :: sinc(:,:,:,:)                  !< Sync function for Fourier components parallel to cube faces
    FLOAT, pointer   :: LLr(:,:)                       !< Coordinates of the face edges
    integer, pointer :: NN(:,:)                        !< Number of points on each face mapping coord

    integer          :: tdsteps                        !< = sys%outp%restart_write_interval (M_PLANES/M_CUBIC)
                                                       !< = mesh%mpi_grp%size (M_SPHERICAL)
    integer          :: max_iter                       !< td%max_iter
    integer          :: save_iter                      !< sys%outp%restart_write_interval
    integer          :: itstep                         !< 1 <= itstep <= tdsteps

    CMPLX, pointer   :: wf(:,:,:,:,:)                  !< wavefunction
    CMPLX, pointer   :: gwf(:,:,:,:,:,:)               !< gradient of wavefunction
    FLOAT, pointer   :: veca(:,:)                      !< vector potential
    CMPLX, pointer   :: conjgphase_prev_cub(:,:)       !< Volkov phase for all k-points from previous time step
    CMPLX, pointer   :: conjgphase_prev_sph(:,:)
    CMPLX, pointer   :: conjgplanewf_cub(:,:,:)        !< plane wave factor
    CMPLX, pointer   :: spctramp_cub(:,:,:,:)          !< spectral amplitude
    CMPLX, pointer   :: spctramp_sph(:,:,:,:,:)

    integer, public  :: ll(3)                          !< the dimensions of a cubic mesh containing the momentum-space
                                                       !< mesh. Used when working with semi-periodic systems 
    integer, public  :: ngpt                           !< Number of free Gpoints use to increase resoltion                        

    logical          :: usememory                      !< whether conjgplanewf should be kept in memory
    logical          :: avoid_ab
    type(mesh_interpolation_t) :: interp
      
    logical          :: parallel_in_momentum           !< whether we are parallelizing over the k-mesh  
    logical          :: arpes_grid
    logical          :: surf_interp                    !< interpolate points on the surface

    integer          :: dim                            !< dimensions
    integer          :: pdim                           !< periodi dimensions


  end type pes_flux_t

  integer, parameter ::   &
    M_CUBIC      = 1,     &
    M_SPHERICAL  = 2,     &
    M_PLANES     = 3

  integer, parameter ::   &
    M_POLAR      = 1,     &
    M_CARTESIAN  = 2
    

contains


  ! ---------------------------------------------------------
  subroutine pes_flux_init(this, namespace, mesh, st, hm, save_iter, max_iter)
    type(pes_flux_t),    intent(inout) :: this
    type(namespace_t),   intent(in)    :: namespace
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: st
    type(hamiltonian_elec_t), intent(in)    :: hm
    integer,             intent(in)    :: save_iter
    integer,             intent(in)    :: max_iter

    type(block_t)      :: blk
    FLOAT              :: border(MAX_DIM)       ! distance of surface from border
    FLOAT              :: offset(MAX_DIM)       ! offset for border
    integer            :: stst, stend, kptst, kptend, sdim, mdim, pdim
    integer            :: imdim
    integer            :: isp, ikp
    integer            :: il, ik
    integer            :: ikk, ith, iph, iomk
    FLOAT              :: kmax, kact, thetak, phik, kpoint(1:3)

    FLOAT, allocatable :: k_dot_aux(:)
    integer            :: nstepsphir, nstepsthetar
    integer            :: ll, mm
    integer            :: default_shape
    FLOAT              :: fc_ptdens        !< density of points on face 
    

    PUSH_SUB(pes_flux_init)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = mesh%sb%dim
    pdim   = mesh%sb%periodic_dim

    this%surf_interp = .false.

    do il = 1, hm%ep%no_lasers
      if(laser_kind(hm%ep%lasers(il)) /= E_FIELD_VECTOR_POTENTIAL) then
        message(1) = 't-surff only works in velocity gauge.'
        call messages_fatal(1, namespace=namespace)
      end if
    end do

    message(1) = 'Info: Calculating PES using t-surff technique.'
    call messages_info(1)

    ! -----------------------------------------------------------------
    ! Setting up r-mesh (the surface)
    ! -----------------------------------------------------------------
    !%Variable PES_Flux_Shape
    !%Type integer
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The shape of the surface.
    !%Option cub 1
    !% Uses a parallelepiped as surface. All surface points are grid points.
    !% Choose the location of the surface with variable <tt>PES_Flux_Lsize</tt>
    !% (default for 1D and 2D).
    !%Option sph 2
    !% Constructs a sphere with radius <tt>PES_Flux_Radius</tt>. Points on the sphere
    !% are interpolated by trilinear interpolation (default for 3D).
    !%Option pln 3
    !% This option is for periodic systems. 
    !% Constructs planes perpendicular to the non-periodic dimension 
    !% symmetrically placed at positive and negative values of <tt>PES_Flux_Lsize</tt>.
    !%End
    default_shape = M_SPHERICAL
    if(mdim <= 2) default_shape = M_CUBIC
    if(simul_box_is_periodic(mesh%sb)) default_shape = M_PLANES
    
    call parse_variable(namespace, 'PES_Flux_Shape', default_shape, this%surf_shape)
    if(.not.varinfo_valid_option('PES_Flux_Shape', this%surf_shape, is_flag = .true.)) &
      call messages_input_error('PES_Flux_Shape')
    if(this%surf_shape == M_SPHERICAL .and. mdim /= 3) then
      message(1) = 'Spherical grid works only in 3d.'
      call messages_fatal(1, namespace=namespace)
    end if
    call messages_print_var_option(stdout, 'PES_Flux_Shape', this%surf_shape)

    !%Variable PES_Flux_Offset
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Shifts the surface for <tt>PES_Flux_Shape = cub</tt>. The syntax is:
    !%
    !% <tt>%PES_Flux_Offset
    !% <br>&nbsp;&nbsp;xshift | yshift | zshift
    !% <br>%
    !% </tt>
    !%End
    offset = M_ZERO
    if(parse_block(namespace, 'PES_Flux_Offset', blk) == 0) then
      do imdim = 1, mdim
        call parse_block_float(blk, 0, imdim - 1, offset(imdim))
      end do
      call parse_block_end(blk)
    end if

    !%Variable PES_Flux_Lmax
    !%Type integer
    !%Default 1
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Maximum order of the spherical harmonic to be integrated on an equidistant spherical 
    !% grid (to be changed to Gauss-Legendre quadrature).
    !%End
    if(this%surf_shape == M_SPHERICAL) then
      call parse_variable(namespace, 'PES_Flux_Lmax', 1, this%lmax)
      if(this%lmax < 1) call messages_input_error('PES_Flux_Lmax', 'must be > 0')
      call messages_print_var_value(stdout, 'PES_Flux_Lmax', this%lmax)
    end if

    this%avoid_ab = .false.
    if(this%surf_shape == M_CUBIC .or. this%surf_shape == M_PLANES) then
      if(hm%bc%abtype /= NOT_ABSORBING) then
        !%Variable PES_Flux_AvoidAB
        !%Type logical
        !%Default yes
        !%Section Time-Dependent::PhotoElectronSpectrum
        !%Description
        !% For <tt>PES_Flux_Shape = cub</tt>, checks whether surface points are inside the
        !% absorbing zone and discards them if set to yes (default).
        !%End
        call parse_variable(namespace, 'PES_Flux_AvoidAB', .true., this%avoid_ab)
        call messages_print_var_value(stdout, 'PES_Flux_AvoidAB', this%avoid_ab)
      end if

      !%Variable PES_Flux_Lsize
      !%Type block
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% For <tt>PES_Flux_Shape = cub</tt> sets the dimensions along each direction. The syntax is:
      !%
      !% <tt>%PES_Flux_Lsize
      !% <br>&nbsp;&nbsp;xsize | ysize | zsize
      !% <br>%
      !% </tt>
      !%
      !% where <tt>xsize</tt>, <tt>ysize</tt>, and <tt>zsize</tt> are with respect to the origin. The surface can
      !% be shifted with <tt>PES_Flux_Offset</tt>.
      !% If <tt>PES_Flux_Shape = pln</tt>, specifies the position of two planes perpendicular to
      !% the non-periodic dimension symmetrically placed at <tt>PES_Flux_Lsize</tt> distance from
      !% the origin.
      !%End
      if(parse_block(namespace, 'PES_Flux_Lsize', blk) == 0) then
        do imdim = 1, mdim
          call parse_block_float(blk, 0, imdim - 1, border(imdim))
        end do
        ! Snap the face to the closest grid point
        border(1:mdim) = floor(border(1:mdim)/mesh%spacing(1:mdim))*mesh%spacing(1:mdim)
        call parse_block_end(blk)

      else if (parse_is_defined(namespace, 'PES_Flux_Lsize')) then 
        border(mdim)  = mesh%sb%lsize(mdim) * M_HALF
        if (simul_box_is_periodic(mesh%sb)) then        
          ! the cube sides along the periodic directions are out of the simulation box
          border(1:pdim)= mesh%sb%lsize(1:pdim) * M_TWO 
          call parse_variable(namespace, 'PES_Flux_Lsize', border(mdim), border(mdim))
          ! Snap the plane to the closest grid point
          border(mdim) = floor(border(mdim)/mesh%spacing(mdim))*mesh%spacing(mdim) 
        else 
          call parse_variable(namespace, 'PES_Flux_Lsize', border(mdim), border(mdim))
          border(1:mdim-1) = border(mdim)
          border(1:mdim) = floor(border(1:mdim)/mesh%spacing(1:mdim))*mesh%spacing(1:mdim)            
        end if
      else
        select case(mesh%sb%box_shape)
        case(PARALLELEPIPED)
          border(1:mdim) = mesh%sb%lsize(1:mdim) * M_HALF
        case(SPHERE)
          border(1:mdim) = mesh%sb%rsize/sqrt(M_TWO) * M_HALF
        case default
          call messages_write('PES_Flux_Lsize not specified. No default values available for this box shape.')
          call messages_new_line()
          call messages_write('Specify the location of the parallelepiped with block PES_Flux_Lsize.')
          call messages_fatal(namespace=namespace)
        end select
        call messages_write('PES_Flux_Lsize not specified. Using default values.')
        call messages_info()
      end if

      call messages_print_var_value(stdout, 'PES_Flux_Lsize', border(1:mdim))

      !%Variable PES_Flux_Face_Dens
      !%Type block
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Define the number of points density per unit of area (in au) on the 
      !% face of the 'cub' surface.
      !%End
      if(parse_is_defined(namespace, 'PES_Flux_Face_Dens')) then
        this%surf_interp = .true.
        call parse_variable(namespace, 'PES_Flux_Face_Dens', M_ONE, fc_ptdens)
        call messages_print_var_value(stdout, 'PES_Flux_Face_Dens', fc_ptdens)
      end if

    else
      
      this%surf_interp = .true.
      
      !%Variable PES_Flux_Radius
      !%Type float
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% The radius of the sphere, if <tt>PES_Flux_Shape == sph</tt>.
      !%End
      if(parse_is_defined(namespace, 'PES_Flux_Radius')) then
        call parse_variable(namespace, 'PES_Flux_Radius', M_ZERO, this%radius)
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
          call messages_fatal(2, namespace=namespace)
        end select
        message(1) = 'PES_Flux_Radius not specified. Using default values.'
        call messages_info(1)
        call messages_print_var_value(stdout, 'PES_Flux_Radius', this%radius)
      end if

    end if

    !%Variable PES_Flux_StepsThetaR
    !%Type integer
    !%Default: 2 <tt>PES_Flux_Lmax</tt> + 1
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in <math>\theta</math> (<math>0 \le \theta \le \pi</math>) for the spherical surface.
    !%End
    if(this%surf_shape == M_SPHERICAL) then
      call parse_variable(namespace, 'PES_Flux_StepsThetaR', 2*this%lmax + 1, nstepsthetar)
      if(nstepsthetar < 0) call messages_input_error('PES_Flux_StepsThetaR')
      call messages_print_var_value(stdout, "PES_Flux_StepsThetaR", nstepsthetar)
    end if

    !%Variable PES_Flux_StepsPhiR
    !%Type integer
    !%Default: 2 <tt>PES_Flux_Lmax</tt> + 1
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Number of steps in <math>\phi</math> (<math>0 \le \phi \le 2 \pi</math>) for the spherical surface.
    !%End
    if(this%surf_shape == M_SPHERICAL) then
      call parse_variable(namespace, 'PES_Flux_StepsPhiR', 2*this%lmax + 1, nstepsphir)
      if(nstepsphir < 0) call messages_input_error('PES_Flux_StepsPhiR')
      call messages_print_var_value(stdout, "PES_Flux_StepsPhiR", nstepsphir)
    end if
    
    
    if(this%surf_interp)  call mesh_interpolation_init(this%interp, mesh)

    ! -----------------------------------------------------------------
    ! Get the surface points
    ! -----------------------------------------------------------------
    if(this%surf_shape == M_CUBIC .or. this%surf_shape == M_PLANES) then
      call pes_flux_getcube(this, mesh, hm, border, offset, fc_ptdens)
    else
      ! equispaced grid in theta & phi (Gauss-Legendre would optimize to nstepsthetar = this%lmax & nstepsphir = 2*this%lmax + 1):
      ! nstepsthetar = M_TWO * this%lmax + 1
      ! nstepsphir   = M_TWO * this%lmax + 1
      call pes_flux_getsphere(this, mesh, nstepsthetar, nstepsphir, offset)
    end if

    ! distribute the surface points on nodes,
    ! since mesh domains may have different numbers of surface points.
    call pes_flux_distribute(1, this%nsrfcpnts, this%nsrfcpnts_start, this%nsrfcpnts_end, mesh%mpi_grp%comm)
    if(debug%info) then
#if defined(HAVE_MPI)
      call MPI_Barrier(mpi_world%comm, mpi_err)
      write(*,*) &
        'Debug: surface points on node ', mpi_world%rank, ' : ', this%nsrfcpnts_start, this%nsrfcpnts_end
      call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
    end if

    if(this%surf_shape == M_PLANES) then
      this%nsrfcpnts_start = 1
      this%nsrfcpnts_end   = this%nsrfcpnts
    end if
      
    if(debug%info .and. mpi_grp_is_root(mpi_world)) then
      do isp = 1, this%nsrfcpnts
        write(223,*) isp, this%rcoords(:, isp), this%srfcnrml(:,isp)
      end do
      if(this%nsrfcpnts > 0) flush(223)
    end if

    ! get the values of the spherical harmonics for the surface points for M_SPHERICAL
    if(this%surf_shape == M_SPHERICAL) then
      SAFE_ALLOCATE(this%ylm_r(0:this%lmax, -this%lmax:this%lmax, this%nsrfcpnts_start:this%nsrfcpnts_end))
      this%ylm_r = M_z0

      if(this%nsrfcpnts_start > 0) then
        do isp = this%nsrfcpnts_start, this%nsrfcpnts_end
          do ll = 0, this%lmax
            do mm = -ll, ll
              call ylmr(this%rcoords(1, isp), this%rcoords(2, isp), this%rcoords(3, isp), ll, mm, this%ylm_r(ll, mm, isp))
              this%ylm_r(ll, mm, isp) = conjg(this%ylm_r(ll, mm, isp))
            end do
          end do
        end do
      end if

    end if

    ! Generate the reciprocal space mesh grid
    call pes_flux_reciprocal_mesh_gen(this, namespace, mesh%sb, st, mesh%mpi_grp%comm)

    ! -----------------------------------------------------------------
    ! Options for time integration 
    ! -----------------------------------------------------------------
    this%max_iter  = max_iter
    this%save_iter = save_iter
    this%itstep    = 0

    if(this%surf_shape == M_CUBIC .or. this%surf_shape == M_PLANES) then
      this%tdsteps = save_iter
    else
      this%tdsteps = 1
#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) this%tdsteps = mesh%mpi_grp%size
#endif
    end if

    !FIXME: CURRENT HACK 
    if(this%surf_shape == M_CUBIC .or. this%surf_shape == M_PLANES) then
      call pes_flux_tabulate_cub_pw(this, mesh, st)
      this%tdsteps = 1
#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) this%tdsteps = mesh%mpi_grp%size
#endif

    end if
   
    ! -----------------------------------------------------------------
    ! Other stuff
    ! -----------------------------------------------------------------
    !%Variable PES_Flux_UseMemory
    !%Type logical
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Use memory to tabulate Volkov plane-wave components on the surface.
    !% This option speeds up calculations by precomputing plane-wave phases
    !% on the suface. 
    !% By default true when <tt>PES_Flux_Shape = cub</tt>.
    !%End
    call parse_variable(namespace, 'PES_Flux_UseMemory', .true., this%usememory)
    call messages_print_var_value(stdout, "PES_Flux_UseMemory", this%usememory)            

    SAFE_ALLOCATE(this%wf(stst:stend, 1:sdim, kptst:kptend, 0:this%nsrfcpnts, 1:this%tdsteps))
    this%wf = M_z0

    SAFE_ALLOCATE(this%gwf(stst:stend, 1:sdim, kptst:kptend, 0:this%nsrfcpnts, 1:this%tdsteps, 1:mdim))
    this%gwf = M_z0

    SAFE_ALLOCATE(this%veca(1:mdim, 1:this%tdsteps))
    this%veca = M_ZERO

    if(this%surf_shape == M_SPHERICAL) then
      SAFE_ALLOCATE(this%spctramp_sph(stst:stend, 1:sdim, kptst:kptend, 1:this%nk, 1:this%nstepsomegak))
      this%spctramp_sph = M_z0

      SAFE_ALLOCATE(this%conjgphase_prev_sph(1:this%nk, 1:this%nstepsomegak))
      this%conjgphase_prev_sph = M_z1

    else
      SAFE_ALLOCATE(this%spctramp_cub(stst:stend, 1:sdim, kptst:kptend, 1:this%nkpnts))
      this%spctramp_cub = M_z0

      SAFE_ALLOCATE(this%conjgphase_prev_cub(1:this%nkpnts, kptst:kptend))
      this%conjgphase_prev_cub = M_z1

      if(this%usememory) then
        SAFE_ALLOCATE(k_dot_aux(1:this%nkpnts))
        SAFE_ALLOCATE(this%conjgplanewf_cub(1:this%nkpnts, this%nsrfcpnts_start:this%nsrfcpnts_end, kptst:kptend))
        this%conjgplanewf_cub = M_z0

        do ik = kptst, kptend
          do isp = this%nsrfcpnts_start, this%nsrfcpnts_end

            k_dot_aux = M_ZERO
            do imdim = 1, mdim
              k_dot_aux(:) = k_dot_aux(:) + this%kcoords_cub(imdim, :, ik) * this%rcoords(imdim, isp)
            end do

            this%conjgplanewf_cub(:, isp, ik) = exp(-M_zI * k_dot_aux(:)) / (M_TWO * M_PI)**(mdim/M_TWO)

          end do
        end do

        SAFE_DEALLOCATE_A(k_dot_aux)

      end if
    end if

    call messages_write('Info: Total number of surface points = ')
    call messages_write(this%nsrfcpnts)
    call messages_info() 

    call messages_write('Info: Total number of momentum points =  ')
    call messages_write(this%nkpnts)
    call messages_info()

    POP_SUB(pes_flux_init)
  end subroutine pes_flux_init

  ! ---------------------------------------------------------
  subroutine pes_flux_end(this)
    type(pes_flux_t), intent(inout) :: this

    PUSH_SUB(pes_flux_end)

    if(this%surf_shape == M_SPHERICAL) then
      SAFE_DEALLOCATE_P(this%kcoords_sph)
      SAFE_DEALLOCATE_P(this%ylm_k)
      SAFE_DEALLOCATE_P(this%j_l)
      SAFE_DEALLOCATE_P(this%ylm_r)
      SAFE_DEALLOCATE_P(this%conjgphase_prev_sph)
      SAFE_DEALLOCATE_P(this%spctramp_sph)
      SAFE_DEALLOCATE_P(this%krad)
    else
      SAFE_DEALLOCATE_P(this%kcoords_cub)
      SAFE_DEALLOCATE_P(this%conjgphase_prev_cub)
      SAFE_DEALLOCATE_P(this%spctramp_cub)
      if(this%usememory) then
        SAFE_DEALLOCATE_P(this%conjgplanewf_cub)
      end if

      SAFE_DEALLOCATE_P(this%srfcpnt)
      SAFE_DEALLOCATE_P(this%rankmin)
      
      SAFE_DEALLOCATE_P(this%face_idx_range)
      SAFE_DEALLOCATE_P(this%expg)
      SAFE_DEALLOCATE_P(this%sinc)
      SAFE_DEALLOCATE_P(this%LLr)
      SAFE_DEALLOCATE_P(this%NN)
    end if

    SAFE_DEALLOCATE_P(this%srfcnrml)
    SAFE_DEALLOCATE_P(this%rcoords)

    SAFE_DEALLOCATE_P(this%wf)
    SAFE_DEALLOCATE_P(this%gwf)
    SAFE_DEALLOCATE_P(this%veca)

    POP_SUB(pes_flux_end)
  end subroutine pes_flux_end

  ! ---------------------------------------------------------
  subroutine pes_flux_reciprocal_mesh_gen(this, namespace, sb, st, comm, post)
    type(pes_flux_t),    intent(inout) :: this
    type(namespace_t),   intent(in)    :: namespace
    type(simul_box_t),   intent(in)    :: sb
    type(states_elec_t), intent(in)    :: st
    integer,             intent(in)    :: comm
    logical, optional,   intent(in)    :: post !< only fill the data needed for postprocessing  

    integer           :: mdim, pdim
    integer           :: kptst, kptend  
    integer           :: isp, ikp, ikpt, ibz1, ibz2
    integer           :: il, ll, mm, idim
    integer           :: ikk, ith, iph, iomk,ie
    FLOAT             :: kmax, kmin, kact, thetak, phik
    type(block_t)     :: blk
      
    FLOAT             :: Emin, Emax, DE , kvec(1:3) 
    integer           :: NBZ(1:2), nkp_out, nkmin, nkmax
    
    integer             :: ig
    FLOAT, allocatable  :: gpoints(:,:), gpoints_reduced(:,:)
    FLOAT               :: dk(1:3), kpoint(1:3)
    logical             :: use_enegy_grid
      
    PUSH_SUB(pes_flux_reciprocal_mesh_gen)  

    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    mdim   = sb%dim
    pdim   = sb%periodic_dim

    this%dim  = mdim
    this%pdim = pdim

    !%Variable PES_Flux_Momenutum_Grid
    !%Type integer
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Decides how the grid in momentum space is generated.
    !%Option polar 1
    !% The grid is in polar coordinates with the zenith axis is along z. 
    !% The grid parameters are defined by PES_Flux_Kmax, PES_Flux_DeltaK, 
    !% PES_Flux_StepsThetaK, PES_Flux_StepsPhiK.
    !% This is the default choice for PES_Flux_Shape = sph or cub.
    !%Option cartesian 2
    !% The grid is in cartesian coordinates with parameters defined by
    !% PES_Flux_ARPES_grid, PES_Flux_EnergyGrid.
    !% This is the default choice for PES_Flux_Shape = sph or cub.
    !%End

    ! default values
    if (this%surf_shape == M_SPHERICAL .or. this%surf_shape == M_CUBIC) then
      this%kgrid = M_POLAR
    else
      this%kgrid = M_CARTESIAN
    end if
    
    call parse_variable(namespace, 'PES_Flux_Momenutum_Grid', this%kgrid, this%kgrid)
    if(.not.varinfo_valid_option('PES_Flux_Momenutum_Grid', this%kgrid, is_flag = .true.)) &
      call messages_input_error('PES_Flux_Momenutum_Grid')
    if(this%surf_shape == M_SPHERICAL .and. mdim /= 3) then
      message(1) = 'Spherical grid works only in 3d.'
      call messages_fatal(1, namespace=namespace)
    end if
    call messages_print_var_option(stdout, 'PES_Flux_Momenutum_Grid', this%kgrid)

    !%Variable PES_Flux_EnergyGrid
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The block <tt>PES_Flux_EnergyGrid</tt> specifies the energy grid 
    !% in momentum space. 
    !% <tt><br>%PES_Flux_EnergyGrid
    !% <br>&nbsp;&nbsp; Emin | Emax | DeltaE
    !% <br>%</tt>
    !%End
    
    Emin = 0
    Emax = 10 
    De   = CNST(0.1)
    use_enegy_grid = .false.
    
    if(parse_block(namespace, 'PES_Flux_EnergyGrid', blk) == 0) then

      call parse_block_float(blk, 0, 0, Emin)
      call parse_block_float(blk, 0, 1, Emax)
      call parse_block_float(blk, 0, 2, DE)
      
      Emin = units_to_atomic(units_inp%energy, Emin)
      Emax = units_to_atomic(units_inp%energy, Emax)
      DE   = units_to_atomic(units_inp%energy, DE)
      
      call parse_block_end(blk)
      use_enegy_grid = .true.
      
    
      call messages_write("Energy grid (Emin, Emax, DE) [") 
      call messages_write(trim(units_abbrev(units_out%energy)))
      call messages_write("]:  (")
      call messages_write(Emin,fmt = '(f8.3)')
      call messages_write(", ")
      call messages_write(Emax, fmt = '(f8.3)')
      call messages_write(", ")
      call messages_write(DE, fmt = '(e9.2)')
      call messages_write(")")
      call messages_info()
    
    end if 

    
!     if (this%surf_shape == M_SPHERICAL .or. this%surf_shape == M_CUBIC) then
    if (this%kgrid == M_POLAR) then  
      ! -----------------------------------------------------------------
      ! Setting up k-mesh
      ! 1D = 2 points, 2D = polar coordinates, 3D = spherical coordinates
      ! -----------------------------------------------------------------
      if (.not. use_enegy_grid) then
        
        !%Variable PES_Flux_Kmax
        !%Type float
        !%Default 1.0
        !%Section Time-Dependent::PhotoElectronSpectrum
        !%Description
        !% The maximum value of |k|.
        !%End
        call parse_variable(namespace, 'PES_Flux_Kmax', M_ONE, kmax)
        call messages_print_var_value(stdout, "PES_Flux_Kmax", kmax)
        if(kmax <= M_ZERO) call messages_input_error('PES_Flux_Kmax')

        !%Variable PES_Flux_Kmin
        !%Type float
        !%Default 1.0
        !%Section Time-Dependent::PhotoElectronSpectrum
        !%Description
        !% The minimum value of |k|.
        !%End
        call parse_variable(namespace, 'PES_Flux_Kmin', M_ONE, kmin)
        call messages_print_var_value(stdout, "PES_Flux_Kmin", kmin)
        if(kmax <= M_ZERO) call messages_input_error('PES_Flux_Kmin')


        !%Variable PES_Flux_DeltaK
        !%Type float
        !%Default 0.002
        !%Section Time-Dependent::PhotoElectronSpectrum
        !%Description
        !% Spacing of the k-mesh in |k| (equidistant).
        !%End
        call parse_variable(namespace, 'PES_Flux_DeltaK', CNST(0.002), this%dk)
        if(this%dk <= M_ZERO) call messages_input_error('PES_Flux_DeltaK')
        call messages_print_var_value(stdout, "PES_Flux_DeltaK", this%dk)

      endif

      !%Variable PES_Flux_StepsThetaK
      !%Type integer
      !%Default 45
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Number of steps in <math>\theta</math> (<math>0 \le \theta \le \pi</math>) for the spherical grid in k.
      !%End
      call parse_variable(namespace, 'PES_Flux_StepsThetaK', 45, this%nstepsthetak)
      if(this%nstepsthetak < 0) call messages_input_error('PES_Flux_StepsThetaK')

      !%Variable PES_Flux_StepsPhiK
      !%Type integer
      !%Default 90
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Number of steps in <math>\phi</math> (<math>0 \le \phi \le 2 \pi</math>) for the spherical grid in k.
      !%End
      call parse_variable(namespace, 'PES_Flux_StepsPhiK', 90, this%nstepsphik)
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
        this%nstepsomegak = this%nstepsphik * (this%nstepsthetak - 1) + 2

      end select

      if(mdim == 3) call messages_print_var_value(stdout, "PES_Flux_StepsThetaK", this%nstepsthetak)
      call messages_print_var_value(stdout, "PES_Flux_StepsPhiK", this%nstepsphik)

      if(use_enegy_grid) then
        this%nk     = nint(abs(Emax-Emin)/DE)
      else  
        this%nk     = nint(abs(kmax-kmin)/this%dk)
      end if
      this%nkpnts = this%nstepsomegak * this%nk

      this%ll(1)      = this%nk  
      this%ll(2)      = this%nstepsphik
      this%ll(3)      = this%nstepsthetak - 1
      this%ll(mdim+1:3) = 1
      
      SAFE_ALLOCATE(this%krad(1:this%nk))
      

    else 
      ! Cartesian

      !%Variable PES_Flux_Gpoint_Upsample
      !%Type integer
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Increase resolution in momentum by adding Gpoints in between adjacent
      !% Kpoints. For each additional Gpoint G an entire Kpoint grid centered at 
      !% G is added to the momentum grid. 
      !%End
      call parse_variable(namespace, 'PES_Flux_Gpoint_Upsample', 1, this%ngpt)
      call messages_print_var_value(stdout, "PES_Flux_Gpoint_Upsample", this%ngpt)
      
      SAFE_ALLOCATE(gpoints(1:mdim,1:this%ngpt))
      SAFE_ALLOCATE(gpoints_reduced(1:mdim,1:this%ngpt)) 
      
      gpoints(:,:) = M_ZERO
      gpoints_reduced(:,:) = M_ZERO
      
      dk(1:mdim) = M_ONE/(sb%kpoints%nik_axis(1:mdim) * this%ngpt)
      do ig =2, this%ngpt
        gpoints_reduced(1:pdim, ig) = gpoints_reduced(1:pdim, ig-1) + dk(1:pdim)
        call kpoints_to_absolute(sb%klattice, gpoints_reduced(1:pdim, ig), gpoints(1:pdim, ig), pdim)
        print *, ig," gpoints_reduced(1:pdim, ig) =", gpoints_reduced(1:pdim, ig)
        print *, ig," gpoints(1:pdim, ig) =", gpoints(1:pdim, ig)
      end do 
      
      
      !%Variable PES_Flux_ARPES_grid
      !%Type logical
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Use a curvilinear momentum space grid that compensates the transformation 
      !% used to obtain ARPES. With this choice ARPES data is laid out on a Cartesian
      !% regular grid.
      !% By default true when <tt>PES_Flux_Shape = pln</tt>.
      !%End
      call parse_variable(namespace, 'PES_Flux_ARPES_grid', .true., this%arpes_grid)
      call messages_print_var_value(stdout, "PES_Flux_ARPES_grid", this%arpes_grid)       
      
                
      
      if (this%arpes_grid) then
  
        nkmax = floor(Emax/DE)
        nkmin = floor(Emin/DE)

      else 
      
        kmax = sqrt(M_TWO*Emax)
        kmin = sqrt(M_TWO*Emin)
        this%dk = sqrt(M_TWO*DE)

        nkmax = floor(kmax/this%dk)
        nkmin = floor(kmin/this%dk)
        
      end if
    
      this%nk = abs(nkmax - nkmin) + 1
    
      !%Variable PES_Flux_BZones
      !%Type block
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% The block <tt>PES_Flux_BZones</tt> specifies the number of Brillouin
      !% zones along each periodic dimensions. This value will be used 
      !% to generate the mesh in reciprocal space.
      !%
      !% <tt><br>%PES_Flux_EnergyGrid
      !% <br>&nbsp;&nbsp; NBZx | NBZy  
      !% <br>%</tt>
      !%
      !% The option can also be specified with the same value for each direction
      !% as <tt>PES_Flux_BZones = NBZ</tt>. 
      !% By default <tt>PES_Flux_BZones = 2</tt>; this means that we will have,
      !% on top of the 1st zone, two additional zones at positive and negative 
      !% reciprocal lattice translation for each periodic dimension. 
      !%End
    
      NBZ(:) = 1
      if (.not. kpoints_have_zero_weight_path(sb%kpoints)) NBZ(1:pdim) = 2
      if(parse_block(namespace, 'PES_Flux_BZones', blk) == 0) then

        call parse_block_integer(blk, 0, 0, NBZ(1))
        call parse_block_integer(blk, 0, 1, NBZ(2))

        call parse_block_end(blk)

      else
        call parse_variable(namespace, 'PES_Flux_BZones', maxval(NBZ(:)), NBZ(1))
        NBZ(:) = NBZ(1)

      end if

      ! If we are using a path in reciprocal space 
      ! we do not need to replicate the BZ in directions 
      ! perpendicular to the path
      if (kpoints_have_zero_weight_path(sb%kpoints) .and. any(NBZ(:)> 1) ) then
        call messages_write("Using a path in reciprocal space with PES_Flux_BZones > 1.")
        call messages_new_line()
        call messages_write("This may cause unphysical results if the path crosses the 1st BZ boundary.")
        call messages_warning(namespace=namespace)
!         call get_kpath_perp_direction(sb%kpoints, idim)
!         if (idim > 0 ) NBZ(idim) = 1
      end if 

      ! This information is needed for postprocessing the data
      this%ll(:)   = 1
      this%ll(1:pdim) = (NBZ(1:pdim)-1)*2+1
      this%ll(mdim)   = this%nk 
      

      call messages_write("Number of Brillouin zones = ")
      do idim = 1 , pdim
        call messages_write( this%ll(idim) )
        if (.not. idim == pdim) call messages_write(" x ")        
      end do 
      call messages_info()
            
      ! Total number of points
      this%nkpnts = product(this%ll(1:mdim))*this%ngpt
      
      

    end if    

  
    this%parallel_in_momentum = .false.

    ! Create the grid
    select case (this%kgrid)

    case (M_POLAR)
    
      if(use_enegy_grid) then
        do ie = 1, this%nk  
          this%krad(ie) = sqrt(M_TWO*(ie * DE + Emin))
        end do
      else  
        do ikk = 1, this%nk  
          this%krad(ikk) = ikk * this%dk + kmin
        end do
      end if
      

      if (this%surf_shape == M_SPHERICAL) then

        if(optional_default(post, .false.)) then 
          POP_SUB(pes_flux_reciprocal_mesh_gen)
          return 
        end if
    
        ! we split the k-mesh in radial & angular part
        call pes_flux_distribute(1, this%nk, this%nk_start, this%nk_end, comm)
        if((this%nk_end - this%nk_start + 1) < this%nk) this%parallel_in_momentum = .true.
        call pes_flux_distribute(1, this%nstepsomegak, this%nstepsomegak_start, this%nstepsomegak_end, comm)

        if(debug%info) then
#if defined(HAVE_MPI)
          call MPI_Barrier(mpi_world%comm, mpi_err)
          write(*,*) &
            'Debug: momentum points on node ', mpi_world%rank, ' : ', this%nk_start, this%nk_end
          call MPI_Barrier(mpi_world%comm, mpi_err)
          write(*,*) &
            'Debug: momentum directions on node ', mpi_world%rank, ' : ', this%nstepsomegak_start, this%nstepsomegak_end
          call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
        end if
        SAFE_ALLOCATE(this%j_l(0:this%lmax, this%nk_start:this%nk_end))
        this%j_l = M_ZERO

        SAFE_ALLOCATE(this%kcoords_sph(1:3, this%nk_start:this%nk_end, 1:this%nstepsomegak))
        this%kcoords_sph = M_ZERO

        SAFE_ALLOCATE(this%ylm_k(0:this%lmax, -this%lmax:this%lmax, this%nstepsomegak_start:this%nstepsomegak_end))
        this%ylm_k = M_z0

        ! spherical harmonics & kcoords_sph
        iomk = 0
        do ith = 0, this%nstepsthetak
          thetak = ith * M_PI / this%nstepsthetak
          do iph = 0, this%nstepsphik - 1
            phik = iph * M_TWO * M_PI / this%nstepsphik
            iomk = iomk + 1
            if(iomk >= this%nstepsomegak_start .and. iomk <= this%nstepsomegak_end) then
              do ll = 0, this%lmax
                do mm = -ll, ll
                  call ylmr(cos(phik) * sin(thetak), sin(phik) * sin(thetak), cos(thetak), ll, mm, this%ylm_k(ll, mm, iomk))
                end do
              end do
            end if
            if(this%nk_start > 0) then
              this%kcoords_sph(1, this%nk_start:this%nk_end, iomk) = cos(phik) * sin(thetak)
              this%kcoords_sph(2, this%nk_start:this%nk_end, iomk) = sin(phik) * sin(thetak)
              this%kcoords_sph(3, this%nk_start:this%nk_end, iomk) = cos(thetak)
            end if
            if(ith == 0 .or. ith == this%nstepsthetak) exit
          end do
        end do

        if(this%nk_start > 0) then
          ! Bessel functions & kcoords_sph
          do ikk = this%nk_start, this%nk_end
            kact = this%krad(ikk)
            do ll = 0, this%lmax
              this%j_l(ll, ikk) = loct_sph_bessel(ll, kact * this%radius) * &
                                  M_TWO * M_PI / (M_TWO * M_PI)**M_THREE/M_TWO
            end do
            this%kcoords_sph(:, ikk, :) = kact * this%kcoords_sph(:, ikk, :)
          end do
        end if

      else 
        !planar or cubic surface
        
        ! no distribution
        this%nkpnts_start = 1 
        this%nkpnts_end   = this%nkpnts


        ! store in the additional kpoint dim the gauge independent grid for final 
        ! momentum representation (i.e. the one at the Gamma point)
        SAFE_ALLOCATE(this%kcoords_cub(1:mdim, 1:this%nkpnts, kptst:kptend+1)) 
        
        this%kcoords_cub = M_ZERO

        select case(mdim)
        case(1)
          ikp = 0
          do ikk = -this%nk, this%nk
            if(ikk == 0) cycle
            ikp = ikp + 1
            kact = this%krad(ikk)
            this%kcoords_cub(1, ikp, kptst:kptend+1) = kact
          end do

        case default
          thetak = M_PI / M_TWO
          do ikpt = kptst, kptend+1
            if (ikpt == kptend+1) then
              kpoint(1:sb%dim) = M_ZERO
            else
              kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, ikpt)
            end if
            ikp = 0            
            do ikk = 1, this%nk
              do ith = 0, this%nstepsthetak
                if(mdim == 3) thetak = ith * M_PI / this%nstepsthetak 
                do iph = 0, this%nstepsphik - 1
                  ikp = ikp + 1
                  phik = iph * M_TWO * M_PI / this%nstepsphik
                  kact = this%krad(ikk)
                                this%kcoords_cub(1, ikp, ikpt) = kact * cos(phik) * sin(thetak) + kpoint(1)
                                this%kcoords_cub(2, ikp, ikpt) = kact * sin(phik) * sin(thetak) + kpoint(2)
                  if(mdim == 3) this%kcoords_cub(3, ikp, ikpt) = kact * cos(thetak)
                  if(mdim == 3 .and. (ith == 0 .or. ith == this%nstepsthetak)) exit
                end do
              end do
            end do
          end do
        end select
      end if

    case (M_CARTESIAN)

  !         call pes_flux_distribute(1, this%nkpnts, this%nkpnts_start, this%nkpnts_end, mpi_world%comm)

      this%nkpnts_start = 1 
      this%nkpnts_end   = this%nkpnts
  

      SAFE_ALLOCATE(this%kcoords_cub(1:mdim, 1:this%nkpnts, kptst:kptend))


      nkp_out = 0 
      do ikpt = kptst, kptend
        ikp = 0
        do ikk = nkmin, nkmax
      
          ! loop over periodic directions
          select case (pdim)
            case (1)
            do ibz1 = -(NBZ(1)-1), (NBZ(1)-1) 
              do ig = 1, this%ngpt
                kvec(1) = ibz1 * sb%klattice(1, 1) + gpoints(1,ig)               
                call fill_non_periodic_dimension(this)
              end do      
            end do

            case (2)
            
            do ibz2 = -(NBZ(2)-1), (NBZ(2)-1) 
              do ibz1 = -(NBZ(1)-1), (NBZ(1)-1) 
                do ig = 1, this%ngpt
                  kvec(1:2) = ibz1 * sb%klattice(1:2, 1) + ibz2 * sb%klattice(1:2, 2) &
                            + gpoints(1:2,ig) 
                  call fill_non_periodic_dimension(this)
                end do
              end do
            end do

          end select

      
        end do
      end do
  
      if (debug%info .and. mpi_grp_is_root(mpi_world)) then
        ! this does not work for parallel in kpoint 
        ! you need to gather kcoords_pln
        write(229,*) "#   ikpt (kpoint index),   ikp (momentum index),   this%kcoords_pln(1:mdim, ikp, ikpt)"
        do ikpt = kptst, kptend
          do ikp = 1, this%nkpnts
            write(229,*) ikpt, ikp, this%kcoords_cub(1:mdim, ikp, ikpt)
          end do
        end do
        flush(229)
      end if

      call messages_write("Number of points with E<p//^2/2 = ")
      call messages_write(nkp_out)
      call messages_write(" [of ")
      call messages_write(this%nkpnts*kpoints_number(sb%kpoints))
      call messages_write("]")
      call messages_info()
      
      SAFE_DEALLOCATE_A(gpoints)
      SAFE_DEALLOCATE_A(gpoints_reduced)

    end select
    
    POP_SUB(pes_flux_reciprocal_mesh_gen)
    
  contains 
    
    ! Fill the non-periodic direction
    subroutine fill_non_periodic_dimension(this)
      type(pes_flux_t),   intent(inout) :: this
        
      integer :: sign
      FLOAT   :: kpar(1:pdim), kpoint(1:3), val

      ikp = ikp + 1
      
      sign = 1         
      if (ikk /= 0) sign= ikk/abs(ikk)        
      
      if (this%arpes_grid) then
        kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, ikpt)
        kpar(1:pdim) = kvec(1:pdim) - kpoint(1:pdim)
        val = abs(ikk)*DE * M_TWO - sum(kpar(1:pdim)**2)
        if (val >= 0) then
          kvec(mdim) =  sign * sqrt(val)
        else  ! if E < p//^2/2
          !FIXME: Should handle the exception doing something smarter than this
          kvec(mdim) = sqrt(val) ! set to NaN
!           kvec(mdim) = M_HUGE ! set to infinity
          nkp_out = nkp_out + 1
        end if
      else         
        kvec(mdim) = ikk * this%dk 
      end if
      
      this%kcoords_cub(1:mdim, ikp, ikpt) =  kvec(1:mdim)
      
    end subroutine fill_non_periodic_dimension
                   
          
    
  end subroutine pes_flux_reciprocal_mesh_gen

  ! ---------------------------------------------------------
  subroutine pes_flux_calc(this, mesh, st, gr, hm, iter, dt, stopping)
    type(pes_flux_t),    intent(inout)    :: this
    type(mesh_t),        intent(in)       :: mesh
    type(states_elec_t), intent(inout)    :: st
    type(grid_t),        intent(in)       :: gr
    type(hamiltonian_elec_t), intent(in)  :: hm
    integer,             intent(in)       :: iter
    FLOAT,               intent(in)       :: dt
    logical,             intent(in)       :: stopping

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim, imdim
    integer            :: isp
    integer            :: il, ip_local
    CMPLX, allocatable :: gpsi(:,:), psi(:)
    CMPLX, allocatable :: interp_values(:)
    logical            :: contains_isp

    PUSH_SUB(pes_flux_calc)

    if(iter > 0) then

      if (debug%info) then
        call messages_write("Debug: Calculating pes_flux")
        call messages_info()
      end if
      
      stst   = st%st_start
      stend  = st%st_end
      kptst  = st%d%kpt%start
      kptend = st%d%kpt%end
      sdim   = st%d%dim
      mdim   = mesh%sb%dim

      SAFE_ALLOCATE(psi(1:mesh%np_part))
      SAFE_ALLOCATE(gpsi(1:mesh%np_part, 1:mdim))

      if(this%surf_interp) then
        SAFE_ALLOCATE(interp_values(1:this%nsrfcpnts))
      end if

      this%itstep = this%itstep + 1

      ! get and save current laser field
      do il = 1, hm%ep%no_lasers
        call laser_field(hm%ep%lasers(il), this%veca(1:mdim, this%itstep), iter*dt)
      end do
      this%veca(:, this%itstep) = - this%veca(:, this%itstep)

      ! save wavefunctions & gradients
      do ik = kptst, kptend
        do ist = stst, stend
          do isdim = 1, sdim
            call states_elec_get_state(st, mesh, isdim, ist, ik, psi)
            call zderivatives_grad(gr%der, psi, gpsi, .true.)

            if(this%surf_interp) then
              call mesh_interpolation_evaluate(this%interp, this%nsrfcpnts, psi(1:mesh%np_part), &
                this%rcoords(1:mdim, 1:this%nsrfcpnts), interp_values(1:this%nsrfcpnts))
              this%wf(ist, isdim, ik, 1:this%nsrfcpnts, this%itstep) = st%occ(ist, ik) * interp_values(1:this%nsrfcpnts)
              do imdim = 1, mdim
                call mesh_interpolation_evaluate(this%interp, this%nsrfcpnts, gpsi(1:mesh%np_part, imdim), &
                  this%rcoords(1:mdim, 1:this%nsrfcpnts), interp_values(1:this%nsrfcpnts))
                this%gwf(ist, isdim, ik, 1:this%nsrfcpnts, this%itstep, imdim) = st%occ(ist, ik) * interp_values(1:this%nsrfcpnts)
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
                  this%wf(ist, isdim, ik, isp, this%itstep) = st%occ(ist, ik) * psi(ip_local)
                  this%gwf(ist, isdim, ik, isp, this%itstep, 1:mdim) = &
                    st%occ(ist, ik) * gpsi(ip_local, 1:mdim) ! * mesh%spacing(1:mdim) ?
                end if
              end do
              if(mesh%parallel_in_domains) then
                call comm_allreduce(mesh%mpi_grp%comm, this%wf(ist, isdim, ik, :, this%itstep))
                do imdim = 1, mdim
                  call comm_allreduce(mesh%mpi_grp%comm, this%gwf(ist, isdim, ik, :, this%itstep, imdim))
                end do
              end if
            end if
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(gpsi)

      if(this%itstep == this%tdsteps .or. mod(iter, this%save_iter) == 0 .or. iter == this%max_iter .or. stopping) then
        if(this%surf_shape == M_CUBIC .or. this%surf_shape == M_PLANES) then
          call pes_flux_integrate_cub_pw(this, mesh, st, dt)
        else
          call pes_flux_integrate_sph(this, mesh, st, dt)
        end if

        ! clean up fields
         this%wf  = M_z0
        this%gwf  = M_z0
        this%veca = M_ZERO
        this%itstep = 0
      end if

    end if

    POP_SUB(pes_flux_calc)
  end subroutine pes_flux_calc



  ! ---------------------------------------------------------
  subroutine pes_flux_integrate_cub(this, mesh, st, dt)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(inout) :: st
    FLOAT,               intent(in)    :: dt

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim, imdim
    integer            :: isp, ikp, itstep
    integer            :: idir
    CMPLX, allocatable :: Jk_cub(:,:,:,:), spctramp_cub(:,:,:,:)
    CMPLX, allocatable :: conjgplanewf_cub(:,:)
    CMPLX, allocatable :: conjgphase_cub(:,:,:)
    FLOAT, allocatable :: k_dot_aux(:)
    integer            :: ikp_start, ikp_end, isp_start, isp_end
    FLOAT              :: vec, kpoint(1:3)

    type(profile_t), save :: prof_init
      
    PUSH_SUB(pes_flux_integrate_cub)

    ! this routine is parallelized over surface points since the number of 
    ! states is in most cases less than the number of surface points

    call profiling_in(prof_init, 'PES_FLUX_INTEGRATE_CUB') 

    if (debug%info) then
      call messages_write("Debug: calculating pes_flux surface integral")
      call messages_info()
    end if


    stst      = st%st_start
    stend     = st%st_end
    kptst     = st%d%kpt%start
    kptend    = st%d%kpt%end
    sdim      = st%d%dim
    mdim      = mesh%sb%dim

    ikp_start = this%nkpnts_start
    ikp_end   = this%nkpnts_end
    isp_start = this%nsrfcpnts_start
    isp_end   = this%nsrfcpnts_end


    SAFE_ALLOCATE(k_dot_aux(1:this%nkpnts))

    SAFE_ALLOCATE(Jk_cub(stst:stend, 1:sdim, kptst:kptend, 1:this%nkpnts))
    SAFE_ALLOCATE(spctramp_cub(stst:stend, 1:sdim, kptst:kptend, 1:this%nkpnts))
    spctramp_cub = M_z0

    if(.not. this%usememory) then
      SAFE_ALLOCATE(conjgplanewf_cub(1:this%nkpnts, kptst:kptend))
    end if

    SAFE_ALLOCATE(conjgphase_cub(1:this%nkpnts, 0:this%tdsteps, kptst:kptend))
    conjgphase_cub = M_z0

    ! calculate Volkov phase using the previous time step
    conjgphase_cub(:, 0,:) = this%conjgphase_prev_cub(:,:)

    do ik = kptst, kptend
      
      kpoint(:) = M_ZERO
      if(simul_box_is_periodic(mesh%sb)) then
        kpoint(1:mdim) = kpoints_get_point(mesh%sb%kpoints, ik)
      end if

      ! integrate over time
      do itstep = 1, this%itstep

        do ikp = ikp_start, ikp_end
          vec = sum((this%kcoords_cub(1:mdim, ikp, ik) - kpoint(1:mdim) - this%veca(1:mdim, itstep) / P_c)**2)
          conjgphase_cub(ikp, itstep, ik) = conjgphase_cub(ikp, itstep - 1, ik) & 
                                            * exp(M_zI * vec * dt / M_TWO)
        end do
        if (this%parallel_in_momentum) call comm_allreduce(mesh%mpi_grp%comm, conjgphase_cub(:,itstep,ik))
      end do

    end do

    this%conjgphase_prev_cub(:,:) = conjgphase_cub(:, this%itstep,:)


    ! integrate over time & surface (on node)
    do isp = isp_start, isp_end
      do idir = 1, mdim
        ! calculate flux only along the surface normal
        if(abs(this%srfcnrml(idir, isp)) <= M_EPSILON) cycle

        Jk_cub = M_z0

        do ik = kptst, kptend
  
          kpoint(:) = M_ZERO
          if(simul_box_is_periodic(mesh%sb)) then
            kpoint(1:mdim) = kpoints_get_point(mesh%sb%kpoints, ik)
          end if
          

          if(.not. this%usememory) then
            k_dot_aux(:) = M_ZERO
            do imdim = 1, mdim
              k_dot_aux(:) = k_dot_aux(:) + this%kcoords_cub(imdim, :, ik) * this%rcoords(imdim, isp)
            end do
            conjgplanewf_cub(:,ik) = exp(-M_zI * k_dot_aux(:)) / (M_TWO * M_PI)**(mdim/M_TWO)
          end if

          do ist = stst, stend
            do isdim = 1, sdim

              ! integrate over time
              do itstep = 1, this%itstep
                Jk_cub(ist, isdim, ik, 1:this%nkpnts) = &
                  Jk_cub(ist, isdim, ik, 1:this%nkpnts) + conjgphase_cub(1:this%nkpnts, itstep, ik) * &
                  (this%wf(ist, isdim, ik, isp, itstep) * &
                   (M_TWO * (this%veca(idir, itstep) / P_c + kpoint(idir)) - this%kcoords_cub(idir, 1:this%nkpnts, ik)) + &
                    this%gwf(ist, isdim, ik, isp, itstep, idir) * M_zI)
              end do

              ! Add the phase contribute at the surface point isp
              if(this%usememory) then
                Jk_cub(ist, isdim, ik, 1:this%nkpnts) = &
                  Jk_cub(ist, isdim, ik, 1:this%nkpnts) * this%conjgplanewf_cub(1:this%nkpnts, isp, ik)
              else
                Jk_cub(ist, isdim, ik, 1:this%nkpnts) = &
                  Jk_cub(ist, isdim, ik, 1:this%nkpnts) * conjgplanewf_cub(1:this%nkpnts, ik)
              end if

            end do ! spin-dimension loop
          end do ! states loop
          
        end do ! kpoint loop 
        spctramp_cub(:,:,:,:) = spctramp_cub(:,:,:,:) + Jk_cub(:,:,:,:) * this%srfcnrml(idir, isp) / M_TWO
        
      end do ! dimension loop
    end do ! surface point loop

    if(this%parallel_in_momentum) then
      do ist = stst, stend
        do isdim = 1, sdim
          do ik = kptst, kptend
            call comm_allreduce(mesh%mpi_grp%comm, spctramp_cub(ist, isdim, ik, :))
          end do
        end do
      end do
    end if

    this%spctramp_cub = this%spctramp_cub + spctramp_cub

    SAFE_DEALLOCATE_A(k_dot_aux)
    SAFE_DEALLOCATE_A(Jk_cub)
    SAFE_DEALLOCATE_A(spctramp_cub)
    SAFE_DEALLOCATE_A(conjgplanewf_cub)
    SAFE_DEALLOCATE_A(conjgphase_cub)

    call profiling_out(prof_init)

    POP_SUB(pes_flux_integrate_cub)
  end subroutine pes_flux_integrate_cub


  ! ---------------------------------------------------------
  subroutine pes_flux_tabulate_cub_pw(this, mesh, st)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: st
      
    integer   :: mdim,nfaces, ikp_start, ikp_end, fdim, ng, nfp, kptst, kptend
    integer   :: ifc, ig, igu, igv,  idim, idir, n_dir, ifp, max_fcpts, ik, ikp
    integer   :: isp_start, isp_end, isp, iguv(1:2)
    FLOAT     :: gvec(this%nsrfcpnts,1:2),  DRS(1:2), shift 
    FLOAT, allocatable :: gg(:,:)

    PUSH_SUB(pes_flux_tabulate_cub_pw)
  
    mdim      = mesh%sb%dim 
    fdim      = mdim - 1

    ikp_start = this%nkpnts_start
    ikp_end   = this%nkpnts_end
    kptst     = st%d%kpt%start
    kptend    = st%d%kpt%end
  
    nfaces = mdim*2
    if(this%surf_shape == M_PLANES) nfaces = 2 ! We only have two planes 
  
  
    ! maximum number of face points
    max_fcpts = maxval(this%face_idx_range(1:nfaces, 2)-this%face_idx_range(1:nfaces, 1)) + 1

    SAFE_ALLOCATE(this%expg(1:max_fcpts,1:max_fcpts,1:nint((nfaces+0.5)/2)))
    SAFE_ALLOCATE(this%sinc(ikp_start:ikp_end, 1:max_fcpts, kptst:kptend, 1:nint((nfaces+0.5)/2)))
  
    this%expg = M_z0
    this%sinc = M_ZERO
    
    do ifc = 1, nint((nfaces+0.5)/2)
      isp_start = this%face_idx_range(ifc, 1)
      isp_end   = this%face_idx_range(ifc, 2)
      ng = isp_end - isp_start + 1 
      nfp = ng ! number of points on face
! print *, "ifc=", ifc, "isp_start/end=", isp_start, isp_end
      
      n_dir = 0 
      do idir = 1, mdim
        if(abs(this%srfcnrml(idir, this%face_idx_range(ifc, 1))) >= M_EPSILON) n_dir = idir
      end do 
            
      SAFE_ALLOCATE(gg(1:2,1:ng))
      
      gg(:,:) = M_ZERO
      
      ig = 1 
      do igu = 1, this%NN(n_dir,1)
        do igv = 1, this%NN(n_dir,2)
          iguv=(/igu,igv/)
          do idim = 1, fdim
            if(this%LLr(n_dir,idim) > M_EPSILON) then
              gg(idim, ig)= (iguv(idim) - this%NN(n_dir,idim)/2 - &
                            (M_ONE + mod(this%NN(n_dir,idim),2))*M_HALF ) * & !even/odd grid centering
                             M_TWO*M_PI/this%LLr(n_dir,idim)          
            end if
          end do
          ig = ig + 1
        end do
      end do
      
! print *, ig, ng
      
      if(debug%info .and. mpi_grp_is_root(mpi_world)) then
        do ig = 1,ng
          write(225,*) ifc, ig, gg(1,ig), gg(2,ig)
        end do
        if(ig > 0) flush(225)
      end if
      
      
! print *, isp_start, isp_end, this%nsrfcpnts
      
      this%expg(:,:,:) = M_z1
      do ig = 1,ng
        do ifp = 1, nfp
          idim = 1
          do idir = 1, mdim 
            if (idir == n_dir ) cycle
            if(this%LLr(n_dir,idim)> M_EPSILON) &
              this%expg(ifp,ig, ifc) = this%expg(ifp,ig, ifc) * exp(M_zI*gg(idim,ig) * &
                                       this%rcoords(idir, isp_start+ifp-1)) * &
                                       (M_PI/this%LLr(n_dir,idim))**(fdim/M_z2)
            idim= idim+1
          end do
        end do
      end do  
        
      this%sinc(:,:,:,:) = M_ONE
      do ik = kptst, kptend
        do ig = 1,ng
          do ikp = ikp_start, ikp_end
            idim = 1
            do idir = 1, mdim 
              if (idir == n_dir ) cycle
              if(abs(gg(idim,ig) - this%kcoords_cub(idir, ikp, ik)) > M_EPSILON) then
                this%sinc(ikp,ig,ik,ifc) = this%sinc(ikp,ig, ik, ifc)* &
                              sin(M_HALF*this%LLr(n_dir,idim) *(gg(idim,ig) - this%kcoords_cub(idir,ikp, ik))) / &
                              (gg(idim,ig) - this%kcoords_cub(idir, ikp, ik))
              else
                this%sinc(ikp,ig,ik,ifc) = M_ONE
              end if   
              if(this%LLr(n_dir,idim)> M_EPSILON) &
                this%sinc(ikp,ig,ik,ifc) = this%sinc(ikp,ig,ik,ifc) *  sqrt(M_TWO/this%LLr(n_dir,idim))            
              idim= idim+1
            end do          
          end do
        end do
      end do      
      
        
      SAFE_DEALLOCATE_A(gg)  
            
    end do
    
    
    
    POP_SUB(pes_flux_tabulate_cub_pw)
  end subroutine pes_flux_tabulate_cub_pw




  ! ---------------------------------------------------------
  subroutine pes_flux_integrate_cub_pw(this, mesh, st, dt)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(inout) :: st
    FLOAT,               intent(in)    :: dt

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim, imdim
    integer            :: isp, ikp, itstep
    integer            :: idir, n_dir, nfaces
    CMPLX, allocatable :: Jk_cub(:,:,:,:), spctramp_cub(:,:,:,:)
    CMPLX, allocatable :: conjgplanewf_cub(:,:)
    CMPLX, allocatable :: conjgphase_cub(:,:,:)
    FLOAT, allocatable :: k_dot_aux(:)
    integer            :: ikp_start, ikp_end, isp_start, isp_end
    FLOAT              :: vec, kpoint(1:3), k_dot_r, expkr 
    integer            :: ifc
    CMPLX, allocatable :: wfpw_node(:,:,:,:), gwfpw_node(:,:,:,:),wfpw(:), gwfpw(:)
    CMPLX, allocatable :: phase(:,:),vphase(:,:)
    
    integer            :: tdstep_on_node
    integer            :: ng, ig, nfp
    
    type(profile_t), save :: prof_init
      
    PUSH_SUB(pes_flux_integrate_cub_pw)

    ! this routine is parallelized over surface points since the number of 
    ! states is in most cases less than the number of surface points

    call profiling_in(prof_init, 'PES_FLUX_INTEGRATE_CUB_PW') 

    if (debug%info) then
      call messages_write("Debug: calculating pes_flux surface integral (pw expression)")
      call messages_info()
    end if


    tdstep_on_node = 1
#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) tdstep_on_node = mesh%mpi_grp%rank + 1
#endif

    stst      = st%st_start
    stend     = st%st_end
    kptst     = st%d%kpt%start
    kptend    = st%d%kpt%end
    sdim      = st%d%dim
    mdim      = mesh%sb%dim

    ikp_start = this%nkpnts_start
    ikp_end   = this%nkpnts_end
    isp_start = this%nsrfcpnts_start
    isp_end   = this%nsrfcpnts_end



    SAFE_ALLOCATE(Jk_cub(stst:stend, 1:sdim, kptst:kptend, ikp_start:ikp_end))
    SAFE_ALLOCATE(spctramp_cub(stst:stend, 1:sdim, kptst:kptend, ikp_start:ikp_end))
    spctramp_cub = M_z0

    
    SAFE_ALLOCATE(vphase(ikp_start:ikp_end, kptst:kptend))
    SAFE_ALLOCATE(phase(ikp_start:ikp_end, kptst:kptend))
    vphase(:,:) = M_z0
    phase(:,:) = M_z0

    
    

!     print *, this%face_idx_range(:, 1)
    
    nfaces = mdim*2
    if(this%surf_shape == M_PLANES) nfaces = 2 ! We only have two planes 
                  
          
    do ifc = 1, nfaces
      
      isp_start = this%face_idx_range(ifc, 1)
      isp_end   = this%face_idx_range(ifc, 2)
      
      ng = isp_end - isp_start + 1 ! faces can have a different number of points
      nfp = ng

      SAFE_ALLOCATE( wfpw_node(stst:stend, 1:sdim, kptst:kptend, 1:ng))
      SAFE_ALLOCATE(gwfpw_node(stst:stend, 1:sdim, kptst:kptend, 1:ng))

      SAFE_ALLOCATE( wfpw(1:ng))
      SAFE_ALLOCATE(gwfpw(1:ng))

      wfpw = M_z0
      gwfpw = M_z0

!       print *,ifc, isp_start, isp_end
      ! get the direction normal to the surface 
      n_dir = 0 
      do idir = 1, mdim
        if(abs(this%srfcnrml(idir, this%face_idx_range(ifc, 1))) >= M_EPSILON) n_dir = idir
      end do 

      itstep = tdstep_on_node
      do ik = kptst, kptend
        do ist = stst, stend
          do isdim = 1, sdim
        
            do ig = 1, ng
              
!                 print *, int(ifc/2), int((ifc+0.5)/2),nint((ifc+0.5)/2), nint((2+0.5)/2), nint((3+0.5)/2)
              gwfpw(ig) = &
                sum(this%gwf(ist, isdim, ik, isp_start:isp_end, itstep, n_dir) &
                  * this%expg(1:nfp,ig,nint((ifc+0.5)/2)) &
                  * this%srfcnrml(n_dir, isp_start:isp_end)) ! surface area element ds
              
              wfpw(ig) = &
                sum(this%wf(ist, isdim, ik, isp_start:isp_end, itstep)        &
                  * this%expg(1:nfp,ig,nint((ifc+0.5)/2)) & 
                  * this%srfcnrml(n_dir, isp_start:isp_end))
        
            end do 
            
            
            if(itstep == tdstep_on_node) then
              gwfpw_node(ist, isdim, ik, :) = gwfpw(:)
               wfpw_node(ist, isdim, ik, :) = wfpw(:)
            end if              
            

          end do 
        end do 
      end do 


      ! get the previous Volkov phase
      vphase(ikp_start:ikp_end,:) = this%conjgphase_prev_cub(ikp_start:ikp_end,:)
      
      Jk_cub(:, :, :, :) = M_z0

      do itstep = 1, this%itstep

        do ik = kptst, kptend
          
          kpoint(1:mdim) = kpoints_get_point(mesh%sb%kpoints, ik)
          do ikp = ikp_start, ikp_end
            vec = sum((this%kcoords_cub(1:mdim, ikp, ik) - kpoint(1:mdim) - this%veca(1:mdim, itstep) / P_c)**2)
            vphase(ikp, ik) = vphase(ikp, ik) * exp(M_zI * vec * dt / M_TWO)

            vec = this%kcoords_cub(n_dir, ikp, ik) * this%rcoords(n_dir, this%face_idx_range(ifc, 1))
            phase(ikp, ik)  = vphase(ikp, ik) * exp(M_zI * vec )/sqrt(M_TWO * M_PI)
          end do
          
          if(itstep /= tdstep_on_node) cycle

          do ist = stst, stend
            do isdim = 1, sdim
              
              ! communicate the current surface integrals
              if(itstep == tdstep_on_node) then
                gwfpw(:) = gwfpw_node(ist, isdim, ik, :)
                 wfpw(:) =  wfpw_node(ist, isdim, ik, :)
              end if


              do ig = 1, ng
                
                Jk_cub(ist, isdim, ik, ikp_start:ikp_end) = Jk_cub(ist, isdim, ik,ikp_start:ikp_end) +  &
                  phase(ikp_start:ikp_end, ik) * ( wfpw(ig) * &
                     (M_TWO * this%veca(n_dir, itstep) / P_c   -  this%kcoords_cub(n_dir, ikp_start:ikp_end, ik)) + &
                      M_zI * gwfpw(ig))* this%sinc(ikp_start:ikp_end, ig, ik, nint((ifc+0.5)/2))
              end do
                            

          end do ! isdim
        end do ! ist     
      end do ! is

    end do !istep
    
    
    spctramp_cub(:,:,:,:) = spctramp_cub(:,:,:,:) + Jk_cub(:,:,:,:) * M_HALF !&
!                             this%srfcnrml(n_dir, this%face_idx_range(ifc, 1)) / M_TWO
      
        
        
        
    SAFE_DEALLOCATE_A(gwfpw)
    SAFE_DEALLOCATE_A(wfpw)
    SAFE_DEALLOCATE_A(gwfpw_node)
    SAFE_DEALLOCATE_A(wfpw_node)
        
    end do ! face loop 


    this%conjgphase_prev_cub(ikp_start:ikp_end,:) = vphase(ikp_start:ikp_end,:)


    if(this%parallel_in_momentum) then
      do ist = stst, stend
        do isdim = 1, sdim
          do ik = kptst, kptend
            call comm_allreduce(mesh%mpi_grp%comm, spctramp_cub(ist, isdim, ik, :))
          end do
        end do
      end do
    end if

    if(mesh%parallel_in_domains) then
!       call comm_allreduce(mesh%mpi_grp%comm, this%conjgphase_prev_cub)
      call comm_allreduce(mesh%mpi_grp%comm, spctramp_cub)
    end if


    this%spctramp_cub = this%spctramp_cub + spctramp_cub

    SAFE_DEALLOCATE_A(Jk_cub)
    SAFE_DEALLOCATE_A(spctramp_cub)
    SAFE_DEALLOCATE_A(phase)
    SAFE_DEALLOCATE_A(vphase)

    call profiling_out(prof_init)

    POP_SUB(pes_flux_integrate_cub_pw)
  end subroutine pes_flux_integrate_cub_pw


  ! ---------------------------------------------------------
  subroutine pes_flux_integrate_sph(this, mesh, st, dt)
    type(pes_flux_t),    intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(inout) :: st
    FLOAT,               intent(in)    :: dt

    integer            :: stst, stend, kptst, kptend, sdim
    integer            :: ist, ik, isdim, imdim
    integer            :: itstep, lmax
    integer            :: ll, mm 
    CMPLX, allocatable :: s1_node(:,:,:,:,:,:), s1_act(:,:,:)
    CMPLX, allocatable :: s2_node(:,:,:,:,:), s2_act(:,:)
    CMPLX, allocatable :: integ11_t(:,:), integ12_t(:,:,:)
    CMPLX, allocatable :: integ21_t(:), integ22_t(:,:)
    CMPLX, allocatable :: spctramp_sph(:,:,:,:,:)
    integer            :: ikk, ikk_start, ikk_end
    integer            :: isp_start, isp_end
    integer            :: iomk, iomk_start, iomk_end
    CMPLX, allocatable :: phase_act(:,:)
    FLOAT              :: vec
    integer            :: tdstep_on_node

    PUSH_SUB(pes_flux_integrate_sph)

    ! this routine is parallelized over time steps

    tdstep_on_node = 1
#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) tdstep_on_node = mesh%mpi_grp%rank + 1
#endif

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim

    lmax       = this%lmax
    ikk_start  = this%nk_start
    ikk_end    = this%nk_end
    isp_start  = this%nsrfcpnts_start
    isp_end    = this%nsrfcpnts_end
    iomk_start = this%nstepsomegak_start
    iomk_end   = this%nstepsomegak_end

    ! surface integral S_lm (for time step on node)
    SAFE_ALLOCATE(s1_node(stst:stend, 1:sdim, kptst:kptend, 0:lmax, -lmax:lmax, 1:3))
    SAFE_ALLOCATE(s2_node(stst:stend, 1:sdim, kptst:kptend, 0:lmax, -lmax:lmax))

    SAFE_ALLOCATE(s1_act(0:lmax, -lmax:lmax, 1:3))
    SAFE_ALLOCATE(s2_act(0:lmax, -lmax:lmax))

    do itstep = 1, this%itstep

      do ik = kptst, kptend
        do ist = stst, stend
          do isdim = 1, sdim

            s2_act = M_z0

            do ll = 0, lmax 
              do mm = -ll, ll
                do imdim = 1, 3
                  ! surface integrals
                  s1_act(ll, mm, imdim) = &
                    sum(this%ylm_r(ll, mm, isp_start:isp_end) * this%srfcnrml(imdim, isp_start:isp_end) * &
                    this%wf(ist, isdim, ik, isp_start:isp_end, itstep))

                  s2_act(ll, mm) = s2_act(ll, mm) + &
                    sum(this%ylm_r(ll, mm, isp_start:isp_end) * this%srfcnrml(imdim, isp_start:isp_end) * &
                    this%gwf(ist, isdim, ik, isp_start:isp_end, itstep, imdim))
                end do
              end do
            end do

            if(mesh%parallel_in_domains) then
              call comm_allreduce(mesh%mpi_grp%comm, s1_act)
              call comm_allreduce(mesh%mpi_grp%comm, s2_act)
            end if

            if(itstep == tdstep_on_node) then
              s1_node(ist, isdim, ik, :, :, :) = s1_act(:,:,:)
              s2_node(ist, isdim, ik, :, :)    = s2_act(:,:)
            end if

          end do
        end do
      end do
    end do

    ! spectral amplitude for k-points on node
    SAFE_ALLOCATE(integ11_t(0:this%nstepsomegak, 1:3))
    SAFE_ALLOCATE(integ21_t(0:this%nstepsomegak))
    SAFE_ALLOCATE(integ12_t(ikk_start:ikk_end, 1:this%nstepsomegak, 1:3))
    SAFE_ALLOCATE(integ22_t(ikk_start:ikk_end, 1:this%nstepsomegak))

    SAFE_ALLOCATE(spctramp_sph(stst:stend, 1:sdim, kptst:kptend, 0:this%nk, 1:this%nstepsomegak))
    spctramp_sph = M_z0

    ! get the previous Volkov phase
    SAFE_ALLOCATE(phase_act(ikk_start:ikk_end, 1:this%nstepsomegak))
    if(ikk_start > 0) then
      phase_act(ikk_start:ikk_end, :) = this%conjgphase_prev_sph(ikk_start:ikk_end, :)
    else
      phase_act(ikk_start:ikk_end, :) = M_z0
    end if

    do itstep = 1, this%itstep
      ! calculate the current Volkov phase
      do ikk = ikk_start, ikk_end
        do iomk = 1, this%nstepsomegak
          vec = sum((this%kcoords_sph(1:3, ikk, iomk) - this%veca(1:3, itstep) / P_c)**M_TWO)
          phase_act(ikk, iomk) = phase_act(ikk, iomk) * exp(M_zI * vec * dt / M_TWO)
        end do
      end do

      do ik = kptst, kptend
        do ist = stst, stend
          do isdim = 1, sdim

            ! communicate the current surface integrals
            if(itstep == tdstep_on_node) then
              s1_act(:,:,:) = s1_node(ist, isdim, ik, :, :, :)
              s2_act(:,:) = s2_node(ist, isdim, ik, :, :)
            end if
#if defined(HAVE_MPI)
            if(mesh%parallel_in_domains) then
              call MPI_Bcast(s1_act, (lmax+1)*(2*lmax+1)*3, MPI_CMPLX, itstep - 1, mesh%mpi_grp%comm, mpi_err)
              call MPI_Bcast(s2_act, (lmax+1)*(2*lmax+1), MPI_CMPLX, itstep - 1, mesh%mpi_grp%comm, mpi_err)
              call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)
            end if
#endif

            integ12_t = M_z0
            integ22_t = M_z0

            do ll = 0, lmax

              integ11_t = M_z0
              integ21_t = M_z0

              do mm = -ll, ll
                ! multiply with spherical harmonics & sum over all mm
                do imdim = 1, 3
                  integ11_t(iomk_start:iomk_end, imdim) = integ11_t(iomk_start:iomk_end, imdim) + &
                    s1_act(ll, mm, imdim) * this%ylm_k(ll, mm, iomk_start:iomk_end)
                end do
                integ21_t(iomk_start:iomk_end) = integ21_t(iomk_start:iomk_end) + &
                  s2_act(ll, mm) * this%ylm_k(ll, mm, iomk_start:iomk_end)
              end do

              if(mesh%parallel_in_domains) then
                call comm_allreduce(mesh%mpi_grp%comm, integ11_t)
                call comm_allreduce(mesh%mpi_grp%comm, integ21_t)
              end if

              ! multiply with Bessel function & sum over all ll
              do ikk = ikk_start, ikk_end
                integ12_t(ikk, 1:this%nstepsomegak, 1:3) = integ12_t(ikk, 1:this%nstepsomegak, 1:3) + &
                  integ11_t(1:this%nstepsomegak, 1:3) * this%j_l(ll, ikk) * ( - M_zI)**ll
                integ22_t(ikk, 1:this%nstepsomegak) = integ22_t(ikk, 1:this%nstepsomegak) + &
                  integ21_t(1:this%nstepsomegak) * this%j_l(ll, ikk) * ( - M_zI)**ll
              end do
            end do
            ! sum over dimensions
            do imdim = 1, 3
              spctramp_sph(ist, isdim, ik, ikk_start:ikk_end, :) = &
                spctramp_sph(ist, isdim, ik, ikk_start:ikk_end, :) + &
                phase_act(ikk_start:ikk_end,:) * (integ12_t(ikk_start:ikk_end,:, imdim) * &
                 (M_TWO * this%veca(imdim, itstep)  / P_c - this%kcoords_sph(imdim, ikk_start:ikk_end,:)))
            end do
            spctramp_sph(ist, isdim, ik, ikk_start:ikk_end, :) = &
              spctramp_sph(ist, isdim, ik, ikk_start:ikk_end,:) + &
              phase_act(ikk_start:ikk_end,:) * integ22_t(ikk_start:ikk_end,:) * M_zI
          end do
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(s1_node)
    SAFE_DEALLOCATE_A(s2_node)
    SAFE_DEALLOCATE_A(s1_act)
    SAFE_DEALLOCATE_A(s2_act)

    SAFE_DEALLOCATE_A(integ11_t)
    SAFE_DEALLOCATE_A(integ12_t)
    SAFE_DEALLOCATE_A(integ21_t)
    SAFE_DEALLOCATE_A(integ22_t)

    ! save the Volkov phase and the spectral amplitude
    this%conjgphase_prev_sph = M_z0

    if(ikk_start > 0) then
      this%conjgphase_prev_sph(ikk_start:ikk_end, :) = phase_act(ikk_start:ikk_end, :)
    end if
    SAFE_DEALLOCATE_A(phase_act)

    if(mesh%parallel_in_domains) then
      call comm_allreduce(mesh%mpi_grp%comm, this%conjgphase_prev_sph)
      call comm_allreduce(mesh%mpi_grp%comm, spctramp_sph)
    end if

    this%spctramp_sph(:,:,:,1:this%nk,:) = this%spctramp_sph(:,:,:,1:this%nk,:) + spctramp_sph(:,:,:,1:this%nk,:)
    SAFE_DEALLOCATE_A(spctramp_sph)

    POP_SUB(pes_flux_integrate_sph)
  end subroutine pes_flux_integrate_sph

  ! ---------------------------------------------------------
  subroutine pes_flux_getcube(this, mesh, hm, border, offset, fc_ptdens)
    type(mesh_t),     intent(in)    :: mesh
    type(pes_flux_t), intent(inout) :: this
    type(hamiltonian_elec_t), intent(in) :: hm
    FLOAT,            intent(in)    :: border(1:MAX_DIM)
    FLOAT,            intent(in)    :: offset(1:MAX_DIM)
    FLOAT,            intent(in)    :: fc_ptdens

    integer, allocatable  :: which_surface(:)
    FLOAT                 :: xx(MAX_DIM), dd, area, dS(MAX_DIM,1:2), factor
    integer               :: mdim, imdim, idir, isp, pm,nface, idim, ndir, iu,iv, iuv(1:2)
    integer               :: ip_global, npface
    integer               :: rankmin, nsurfaces
    logical               :: in_ab
    integer               :: ip_local, nsrfcpnts, NN(MAX_DIM,1:2), idx(MAX_DIM,1:2) 
    integer               :: isp_end, isp_start, ifc, n_dir, nfaces, mindim
    FLOAT                 :: RSmax(1:2),RSmin(1:2),RS(1:2), dRS(1:2)


    PUSH_SUB(pes_flux_getcube)

    ! this routine is parallelized over the mesh in any case

    mdim = mesh%sb%dim
    
    

    SAFE_ALLOCATE(this%face_idx_range(1:mdim*2,1:2))
    SAFE_ALLOCATE(this%LLr(mdim,1:2))
    SAFE_ALLOCATE(this%NN(mdim,1:2))

    this%face_idx_range(:,:) = 0
    this%LLr(:,:) = M_ZERO
    this%NN(:,:) = 1
    NN(:,:) = 1

    if (this%surf_interp) then
      ! Create a surface with points not on the mesh
      
      idx(:,:) = 0 
      
      mindim  = 1 
      factor = M_TWO
      if(this%surf_shape == M_PLANES) then
        mindim = mdim  ! We only have two planes along the non periodic dimension 
        factor = M_ONE
      end if
      
      
      this%nsrfcpnts = 0 
      do ndir = mdim, mindim, -1
        area = M_ONE
        do idir=1, mdim
          if (idir == ndir) cycle
          area = area * border(idir)*factor
        end do
        
        npface = int(fc_ptdens * area) ! number of points on the face
        
        idim = 1 
        do idir=1, mdim
          if (idir == ndir) cycle
          NN(ndir,idim) = int(npface / (border(idir)*factor)) 
          dS(ndir,idim) = border(idir)*factor/NN(ndir,idim)
          this%LLr(ndir,idim) = NN(ndir,idim)*dS(ndir,idim)
          idx(ndir,idim) = idir
          idim = idim + 1
        end do
        this%nsrfcpnts = this%nsrfcpnts + 2*product(NN(ndir,1:mdim-1))
      end do
      
      this%NN(1:mdim,:)  = NN(1:mdim,:)
      
!       print *,  this%nsrfcpnts, NN(1,:),NN(2,:)
!
!       print *,  this%LLr(1,:), this%LLr(2,:)
      
      SAFE_ALLOCATE(this%srfcnrml(1:mdim, 0:this%nsrfcpnts))
      SAFE_ALLOCATE(this%rcoords(1:mdim, 0:this%nsrfcpnts))
    
      this%srfcnrml(:,:) = M_ZERO
      this%rcoords(:,:)  = M_ZERO
    

      isp = 1  
      ifc = 1
      do ndir = mdim, mindim, -1

        !Up face
        this%face_idx_range(ifc,1) = isp
        do iu = 1, NN(ndir,1)
          do iv = 1, NN(ndir,2)
            this%rcoords(ndir, isp)  =  border(ndir)
            this%srfcnrml(ndir, isp) =  product(dS(ndir,1:mdim-1))             
            iuv =(/iu,iv/)
            do idim = 1, mdim-1
              this%rcoords(idx(ndir,idim), isp) = (-NN(ndir,idim)*M_HALF -M_HALF + iuv(idim)) * dS(ndir,idim)
            end do 
            isp = isp + 1
          end do
        end do
        this%face_idx_range(ifc,2) = isp-1

        ifc = ifc + 1

        !Down face
        this%face_idx_range(ifc,1) = isp
        do iu = 1, NN(ndir,1)
          do iv = 1, NN(ndir,2)
            this%rcoords(ndir, isp)  = -border(ndir)
            this%srfcnrml(ndir, isp) = -product(dS(ndir,1:mdim-1)) 
            iuv =(/iu,iv/)              
            do idim = 1, mdim-1
              this%rcoords(idx(ndir,idim), isp) = (-NN(ndir,idim)*M_HALF -M_HALF + iuv(idim)) * dS(ndir,idim)
            end do              
            isp = isp + 1
          end do
        end do
        this%face_idx_range(ifc,2) = isp-1

        ifc = ifc + 1 
      end do
            
      
    else
      ! Surface points are on the mesh
      
      nfaces = mdim*2
      if(this%surf_shape == M_PLANES) nfaces = 2 ! We only have two planes 

      in_ab = .false.

      SAFE_ALLOCATE(which_surface(1:mesh%np_global))
      which_surface = 0

      ! get the surface points
      this%nsrfcpnts = 0
      do ip_local = 1, mesh%np
        if(mesh%parallel_in_domains) then
          ip_global = mesh%vp%local(mesh%vp%xlocal + ip_local - 1)
        else
          ip_global = ip_local
        end if
      
        nsurfaces = 0

        xx(1:MAX_DIM) = mesh%x(ip_local, 1:MAX_DIM) - offset(1:MAX_DIM)

        ! eventually check whether we are in absorbing zone
        if(this%avoid_ab) then
          select case(hm%bc%abtype)
          case(MASK_ABSORBING)
            in_ab = (hm%bc%mf(ip_local) /= M_ONE)
          case(IMAGINARY_ABSORBING)
            in_ab = (hm%bc%mf(ip_local) /= M_ZERO)
          end select
        end if

        ! check whether the point is inside the cube
        if(all(abs(xx(1:mdim)) <= border(1:mdim)) .and. .not. in_ab) then
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

      if(mesh%parallel_in_domains) then
        call comm_allreduce(mesh%mpi_grp%comm, this%nsrfcpnts)
        call comm_allreduce(mesh%mpi_grp%comm, which_surface)
      end if

      SAFE_ALLOCATE(this%srfcpnt(1:this%nsrfcpnts))
      SAFE_ALLOCATE(this%srfcnrml(1:mdim, 0:this%nsrfcpnts))
      SAFE_ALLOCATE(this%rcoords(1:mdim, 0:this%nsrfcpnts))
      SAFE_ALLOCATE(this%rankmin(1:this%nsrfcpnts))

      this%srfcnrml = M_ZERO
      this%rcoords  = M_ZERO

      this%face_idx_range = 0
    
      isp = 0
      nface = 0
      do idir = mdim, 1, -1
        do pm = -1, 1, 2 
          nface = nface + 1
          this%face_idx_range(nface,1) = isp + 1
        
          do ip_global = 1, mesh%np_global
            if(abs(which_surface(ip_global)) == idir .and. sign(1, which_surface(ip_global)) == pm) then
              isp = isp + 1
              ! coordinate of surface point
              xx(1:MAX_DIM) = mesh_x_global(mesh, ip_global)
              this%rcoords(1:mdim, isp) = xx(1:mdim)
              ! local ip & node which has the surface point
              this%srfcpnt(isp) = mesh_nearest_point(mesh, this%rcoords(1:mdim, isp), dd, rankmin)
              this%rankmin(isp) = rankmin
              ! surface normal
              this%srfcnrml(idir, isp) = sign(1, which_surface(ip_global))
              ! add the surface element (of directions orthogonal to the normal vector)
              do imdim = 1, mdim
                if(imdim == idir) cycle
                this%srfcnrml(idir, isp) = this%srfcnrml(idir, isp) * mesh%spacing(imdim)
              end do
            end if
          end do

          this%face_idx_range(nface,2) = isp
        
!           print *, idir, pm, "this%face_idx_range=",  this%face_idx_range(nface,:)
        end do      
      end do
      
      do ifc = 1, nint((nfaces+0.5)/2)
        isp_start = this%face_idx_range(ifc, 1)
        isp_end   = this%face_idx_range(ifc, 2)
      
        !retrieve face normal direction
        n_dir = 0 
        do idir = 1, mdim
          if(abs(this%srfcnrml(idir, this%face_idx_range(ifc, 1))) >= M_EPSILON) n_dir = idir
        end do 
      
! print *, n_dir,  isp_start, isp_start     
        !retrieve the spacing on the face
        dRS(:) = M_ZERO
        idim = 1
        do idir = 1, mdim 
          if (idir == n_dir ) cycle
          dRS(idim)= mesh%spacing(idir)
          idim = idim+1
        end do

      
        !retrieve the dimensions of the face
        RSmin = M_ZERO
        RSmax = M_ZERO
        do isp = isp_start, isp_end
          idim = 1
          do idir = 1, mdim 
            if (idir == n_dir ) cycle
            RS(idim)=this%rcoords(idir,isp)
            if (RS(idim) < RSmin(idim)) RSmin(idim) = RS(idim)
            if (RS(idim) > RSmax(idim)) RSmax(idim) = RS(idim)
            idim = idim+1
          end do        
        end do
        
   
        do idir = 1, mdim -1
          this%LLr(n_dir,idir) = RSmax(idir) - RSmin(idir)
          if(dRS(idir) > M_ZERO) this%NN(n_dir,idir) = int(this%LLr(n_dir,idir)/dRS(idir))+1
        end do        

!         LL(:) = RSmax(:) - RSmin(:)
!         do idir = 1, fdim
!           NN(idir) = int(LL(idir)/dS(idir))+1
!         end do


! print *, ifc, this%LLr(n_dir,:), dRS(:), this%NN(n_dir,:), isp_end - isp_start+1
      end do
      
      
    end if
    
!     print *, this%LLr(:,:), this%NN(:,:), dS(1,1)
!     print *, ifc, this%LLr, dS, this%NN, isp_end - isp_start+1

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

    if(nstepsthetar <= 1) then
      nstepsphir = 1
      nstepsthetar = 1
    end if
    dthetar = M_PI / nstepsthetar
    this%nsrfcpnts  = nstepsphir * (nstepsthetar - 1) + 2

    SAFE_ALLOCATE(this%srfcnrml(1:mdim, 0:this%nsrfcpnts))
    SAFE_ALLOCATE(this%rcoords(1:mdim, 0:this%nsrfcpnts))

    this%srfcnrml = M_ZERO
    this%rcoords  = M_ZERO

    ! initializing spherical grid
    isp = 0
    do ith = 0, nstepsthetar
      thetar  = ith * dthetar

      if(ith == 0 .or. ith == nstepsthetar) then
        weight = (M_ONE - cos(dthetar/M_TWO)) * M_TWO * M_PI
      else
        weight = abs(cos(thetar - dthetar/M_TWO) - cos(thetar + dthetar/M_TWO)) &
          * M_TWO * M_PI / nstepsphir
      end if

      do iph = 0, nstepsphir - 1     ! 2*pi is the same than zero
        isp = isp + 1
        phir = iph * M_TWO * M_PI / nstepsphir
                      this%srfcnrml(1, isp) = cos(phir) * sin(thetar)
        if(mdim >= 2) this%srfcnrml(2, isp) = sin(phir) * sin(thetar)
        if(mdim == 3) this%srfcnrml(3, isp) = cos(thetar)
        this%rcoords(1:mdim, isp) = this%radius * this%srfcnrml(1:mdim, isp)
        ! here we also include the surface elements
        this%srfcnrml(1:mdim, isp) = this%radius**M_TWO * weight * this%srfcnrml(1:mdim, isp)
        if(ith == 0 .or. ith == nstepsthetar) exit
      end do
    end do

    POP_SUB(pes_flux_getsphere)
  end subroutine pes_flux_getsphere


  ! ---------------------------------------------------------
  subroutine pes_flux_distribute(istart_global, iend_global, istart, iend, comm)
    integer,          intent(in)    :: istart_global
    integer,          intent(in)    :: iend_global
    integer,          intent(inout) :: istart
    integer,          intent(inout) :: iend
    integer,          intent(in)    :: comm

#if defined(HAVE_MPI)
    integer, allocatable :: dimrank(:)
    integer              :: mpisize, mpirank, irest, irank
    integer              :: inumber
#endif

    PUSH_SUB(pes_flux_distribute)

#if defined(HAVE_MPI)
    call mpi_comm_size(comm, mpisize, mpi_err)
    call mpi_comm_rank(comm, mpirank, mpi_err)

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

    if(istart > iend) then
      iend = 0
      istart = 0
    end if

    SAFE_DEALLOCATE_A(dimrank)

#else
    istart = istart_global
    iend   = iend_global
#endif

    POP_SUB(pes_flux_distribute)
  end subroutine pes_flux_distribute


#include "pes_flux_out_inc.F90"

end module pes_flux_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
