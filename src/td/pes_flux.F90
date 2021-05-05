!! Copyright (C) 2015 P. Wopperer and U. De Giovannini
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
  use box_parallelepiped_oct_m
  use box_sphere_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use kpoints_oct_m
  use io_oct_m
  use lalg_adv_oct_m
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
  use space_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
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
    pes_flux_out_energy,       &
    pes_flux_out_vmap

  type pes_flux_t
    private
    !< NOTE: On unfortunate choice of nomenclature. In this module we use the variable k to indicate 
    !< the momentum of photoelectrons. This has to be distinguished from the electrons crystal momentum 
    !< in periodic systems and its discretized values (AKA kpoints). 
    
    integer            :: nkpnts                         !< total number of k-points
    integer            :: nkpnts_start, nkpnts_end       !< start/end of index for k-points on the current node
    integer            :: nk
    integer            :: nk_start, nk_end

    ! spherical momentum grid
    integer            :: nstepsthetak, nstepsphik       !< parameters for k-mesh
    FLOAT              :: thetak_rng(1:2)                !< the range of thetak in [0,pi]
    FLOAT              :: phik_rng(1:2)                  !< the range of phik in [0,2pi]
    integer            :: nstepsomegak
    integer            :: nstepsomegak_start, nstepsomegak_end
    FLOAT              :: ktransf(1:3,1:3)               !< transformation matrix for k-mesh 

    FLOAT, allocatable :: klinear(:,:)                   !< A set of liner grids to help define the k-mesh
                                                       !< polar:     klinear(nk,1)    defines the radial part of the k-mesh
                                                       !< cartesian: klinear(nk,1:3)  defines the the k-mesh along each axes

    FLOAT              :: dk                             !< parameters for k-mesh
    FLOAT, allocatable :: kcoords_cub(:,:,:)             !< coordinates of k-points
    FLOAT, allocatable :: kcoords_sph(:,:,:)
    integer            :: kgrid                          !< how is the grid in k: polar/cartesian

    ! Surface relates quantities
    integer, public    :: surf_shape                     !< shape of the surface (= cube/sphere/planes)
    integer            :: nsrfcpnts                      !< total number of surface points
    integer            :: nsrfcpnts_start, nsrfcpnts_end !< for cubic surface: number of surface points on node
    FLOAT, allocatable :: srfcnrml(:,:)                  !< vectors normal to the surface (includes surface element)
    FLOAT, allocatable :: rcoords(:,:)                   !< coordinates of the surface points
    integer, allocatable :: srfcpnt(:)                     !< for cubic surface: returns local index of the surface points
    integer, allocatable :: rankmin(:)                     !< for cubic surface: returns node which has the surface point
    integer            :: lmax                           !< for spherical surface
    CMPLX, allocatable :: ylm_r(:,:,:)                   !< for spherical surface
    CMPLX, allocatable :: ylm_k(:,:,:)                   !< for spherical surface
    FLOAT, allocatable :: j_l(:,:)                       !< for spherical surface
    FLOAT              :: radius
    integer, allocatable :: face_idx_range(:,:)            !< face_idx_start(nface,1:2) 1 (2) start (end) idx of face nface in rcoords(:,:)

    
    CMPLX,   allocatable :: bvk_phase(:,:)                 !< for cubic surface: Born-von Karman phase
    CMPLX,   allocatable :: expkr(:,:,:,:)                 !< for cubic surface: Exponential tabulated on the cube face
    CMPLX,   allocatable :: expkr_perp(:,:)                !< for cubic surface: Exponential tabulated on the direction perpendicular to the face
    FLOAT,   allocatable :: LLr(:,:)                       !< for cubic surface: coordinates of the face edges
    integer, allocatable :: NN(:,:)                        !< for cubic surface: number of points on each face mapping coord
    integer, allocatable :: Lkpuvz_inv(:,:,:)              !< map a point on the momentum mesh into its components u,v,z (u,v parametric on face
                                                       !< and z perpendicular)

    ! Surface and time integration 
    integer          :: tdsteps                        !< = sys%outp%restart_write_interval (PES_PLANE/PES_CUBIC)
                                                       !< = mesh%mpi_grp%size (PES_SPHERICAL)
    integer          :: max_iter                       !< td%max_iter
    integer          :: save_iter                      !< sys%outp%restart_write_interval
    integer          :: itstep                         !< 1 <= itstep <= tdsteps

    CMPLX, allocatable :: wf(:,:,:,:,:)                  !< wavefunction
    CMPLX, allocatable :: gwf(:,:,:,:,:,:)               !< gradient of wavefunction
    FLOAT, allocatable :: veca(:,:)                      !< vector potential
    CMPLX, allocatable :: conjgphase_prev(:,:)           !< Volkov phase for all k-points from previous time step
    CMPLX, allocatable :: spctramp_cub(:,:,:,:)          !< spectral amplitude
    CMPLX, allocatable :: spctramp_sph(:,:,:,:,:)

    integer, public  :: ll(3)                          !< the dimensions of a cubic mesh containing the momentum-space
                                                       !< mesh. Used when working with semi-periodic systems 

    type(mesh_interpolation_t) :: interp
      
    logical          :: parallel_in_momentum           !< whether we are parallelizing over the k-mesh  
    logical          :: arpes_grid
    logical          :: surf_interp                    !< interpolate points on the surface
    logical          :: use_symmetry              
    logical          :: runtime_output              
    logical          :: anisotrpy_correction                   

    integer          :: par_strategy                   !< parallelization strategy 
    integer          :: dim                            !< simulation box dimensions
    integer          :: pdim                           !< periodic dimensions

  end type pes_flux_t

  integer, parameter ::   &
    PES_CUBIC      = 1,     &
    PES_SPHERICAL  = 2,     &
    PES_PLANE      = 3

  integer, parameter ::   &
    PES_POLAR      = 1,     &
    PES_CARTESIAN  = 2


contains


  ! ---------------------------------------------------------
  subroutine pes_flux_init(this, namespace, space, mesh, st, hm, save_iter, max_iter)
    type(pes_flux_t),         intent(inout) :: this
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(mesh_t),             intent(in)    :: mesh
    type(states_elec_t),      intent(in)    :: st
    type(hamiltonian_elec_t), intent(in)    :: hm
    integer,                  intent(in)    :: save_iter
    integer,                  intent(in)    :: max_iter

    type(block_t)      :: blk
    FLOAT              :: border(MAX_DIM)       ! distance of surface from border
    FLOAT              :: offset(MAX_DIM)       ! offset for border
    integer            :: stst, stend, kptst, kptend, sdim, mdim, pdim
    integer            :: imdim
    integer            :: isp
    integer            :: il

    integer            :: nstepsphir, nstepsthetar
    integer            :: ll, mm
    integer            :: default_shape
    FLOAT              :: fc_ptdens        !< density of points on face 
    
    integer            :: par_strategy
    logical            :: use_symmetry

    PUSH_SUB(pes_flux_init)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    sdim   = st%d%dim
    mdim   = space%dim
    pdim   = space%periodic_dim

    this%surf_interp = .false.


    do il = 1, hm%ext_lasers%no_lasers
      if(laser_kind(hm%ext_lasers%lasers(il)) /= E_FIELD_VECTOR_POTENTIAL) then
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
    !% Constructs a plane perpendicular to the non-periodic dimension 
    !% at <tt>PES_Flux_Lsize</tt>.
    !%End

    default_shape = PES_SPHERICAL
    if (space%is_periodic()) then
      default_shape = PES_PLANE
    else if (mdim <= 2) then
      default_shape = PES_CUBIC
    else
      select type (box => mesh%sb%box)
      type is (box_parallelepiped_t)
        default_shape = PES_CUBIC
      end select
    end if
    
    call parse_variable(namespace, 'PES_Flux_Shape', default_shape, this%surf_shape)
    if(.not.varinfo_valid_option('PES_Flux_Shape', this%surf_shape, is_flag = .true.)) &
      call messages_input_error(namespace,'PES_Flux_Shape')
    if(this%surf_shape == PES_SPHERICAL .and. mdim /= 3) then
      message(1) = 'Spherical grid works only in 3d.'
      call messages_fatal(1, namespace=namespace)
    end if
    call messages_print_var_option(stdout, 'PES_Flux_Shape', this%surf_shape)


    !%Variable PES_Flux_AnisotropyCorrection
    !%Type logical
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Apply anisotropy correction. 
    !%  
    !%End
    if(this%surf_shape == PES_CUBIC) then
      call parse_variable(namespace, 'PES_Flux_AnisotropyCorrection', .false., this%anisotrpy_correction)
      call messages_print_var_value(stdout, 'PES_Flux_AnisotropyCorrection', this%anisotrpy_correction)
    end if

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
    !%Default 80
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Maximum order of the spherical harmonic to be integrated on an equidistant spherical 
    !% grid (to be changed to Gauss-Legendre quadrature).
    !%End
    if(this%surf_shape == PES_SPHERICAL) then
      call parse_variable(namespace, 'PES_Flux_Lmax', 80, this%lmax)
      if(this%lmax < 1) call messages_input_error(namespace,'PES_Flux_Lmax', 'must be > 0')
      call messages_print_var_value(stdout, 'PES_Flux_Lmax', this%lmax)
    end if


    if(this%surf_shape == PES_CUBIC .or. this%surf_shape == PES_PLANE) then
      

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
        border(1:mdim) = int(border(1:mdim) / mesh%spacing(1:mdim)) * mesh%spacing(1:mdim)
        call parse_block_end(blk)

      else if (parse_is_defined(namespace, 'PES_Flux_Lsize')) then 
        border(mdim)  = mesh%sb%lsize(mdim) * M_HALF
        if (space%is_periodic()) then        
          ! the cube sides along the periodic directions are out of the simulation box
          border(1:pdim)= mesh%sb%lsize(1:pdim) * M_TWO 
          call parse_variable(namespace, 'PES_Flux_Lsize', border(mdim), border(mdim))
          ! Snap the plane to the closest grid point
          border(mdim) = floor(border(mdim) / mesh%spacing(mdim)) * mesh%spacing(mdim) 
        else 
          call parse_variable(namespace, 'PES_Flux_Lsize', border(mdim), border(mdim))
          border(1:mdim - 1) = border(mdim)
          border(1:mdim)     = floor(border(1:mdim) / mesh%spacing(1:mdim)) * mesh%spacing(1:mdim)            
        end if
      else
        select type (box => mesh%sb%box)
        type is (box_sphere_t)
          border(1:mdim) = box%radius / sqrt(M_TWO) * M_HALF
        type is (box_parallelepiped_t)
          border(1:mdim) = box%half_length(1:mdim) * M_HALF
        class default
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
        if(this%radius <= M_ZERO) call messages_input_error(namespace, 'PES_Flux_Radius')
        call messages_print_var_value(stdout, 'PES_Flux_Radius', this%radius)
      else
        select type (box => mesh%sb%box)
        type is (box_sphere_t)
          this%radius = box%radius
        type is (box_parallelepiped_t)
          this%radius = minval(box%half_length(1:mdim))
        class default
          message(1) = 'PES_Flux_Radius not specified. No default values available for this box shape.'
          message(2) = 'Specify the radius of the sphere with variable PES_Flux_Radius.'
          call messages_fatal(2, namespace=namespace)
        end select
        message(1) = 'PES_Flux_Radius not specified. Using default values.'
        call messages_info(1)
        call messages_print_var_value(stdout, 'PES_Flux_Radius', this%radius)
      end if

      !%Variable PES_Flux_StepsThetaR
      !%Type integer
      !%Default: 2 <tt>PES_Flux_Lmax</tt> + 1
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Number of steps in <math>\theta</math> (<math>0 \le \theta \le \pi</math>) for the spherical surface.
      !%End
      call parse_variable(namespace, 'PES_Flux_StepsThetaR', 2 * this%lmax + 1, nstepsthetar)
      if(nstepsthetar < 0) call messages_input_error(namespace, 'PES_Flux_StepsThetaR')
      call messages_print_var_value(stdout, "PES_Flux_StepsThetaR", nstepsthetar)

      !%Variable PES_Flux_StepsPhiR
      !%Type integer
      !%Default: 2 <tt>PES_Flux_Lmax</tt> + 1
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Number of steps in <math>\phi</math> (<math>0 \le \phi \le 2 \pi</math>) for the spherical surface.
      !%End
      call parse_variable(namespace, 'PES_Flux_StepsPhiR', 2 * this%lmax + 1, nstepsphir)
      if(nstepsphir < 0) call messages_input_error(namespace, 'PES_Flux_StepsPhiR')
      call messages_print_var_value(stdout, "PES_Flux_StepsPhiR", nstepsphir)

    end if

    
    if(this%surf_interp)  call mesh_interpolation_init(this%interp, mesh)

    ! -----------------------------------------------------------------
    ! Get the surface points
    ! -----------------------------------------------------------------
    if(this%surf_shape == PES_CUBIC .or. this%surf_shape == PES_PLANE) then
      call pes_flux_getcube(this, mesh, hm, border, offset, fc_ptdens)
    else
      ! equispaced grid in theta & phi (Gauss-Legendre would optimize to nstepsthetar = this%lmax & nstepsphir = 2*this%lmax + 1):
      ! nstepsthetar = M_TWO * this%lmax + 1
      ! nstepsphir   = M_TWO * this%lmax + 1
      call pes_flux_getsphere(this, mesh, nstepsthetar, nstepsphir)
    end if
    

    !%Variable PES_Flux_Parallelization
    !%Type flag
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The parallelization strategy to be used to calculate the PES spectrum 
    !% using the resources available in the domain parallelization pool. 
    !% This option is not available without domain parallelization.
    !% Parallelization over k-points and states is always enabled when available.
    !%Option pf_none bit(1)
    !% No parallelization.
    !%Option pf_time bit(2)
    !% Parallelize time integration. This requires to store some quantities over a 
    !% number of time steps equal to the number of cores available. 
    !%Option pf_momentum bit(3)
    !% Parallelize over the final momentum grid. This strategy has a much lower
    !% memory footprint than the one above (time) but seems to provide a smaller
    !% speedup.
    !%Option pf_surface bit(4)
    !% Parallelize over surface points.
    !%
    !%
    !% Option pf_time and pf_surface can be combined: pf_time + pf_surface.
    !%   
    !%End

    this%par_strategy      = OPTION__PES_FLUX_PARALLELIZATION__PF_NONE
    if(mesh%parallel_in_domains) then

      if(this%surf_shape == PES_SPHERICAL) then
        this%par_strategy = OPTION__PES_FLUX_PARALLELIZATION__PF_TIME  &
                          + OPTION__PES_FLUX_PARALLELIZATION__PF_SURFACE
      else 
        this%par_strategy = OPTION__PES_FLUX_PARALLELIZATION__PF_SURFACE    
        if (space%dim == 1) this%par_strategy = OPTION__PES_FLUX_PARALLELIZATION__PF_TIME  
      end if
      par_strategy = this%par_strategy
      call parse_variable(namespace, 'PES_Flux_Parallelization', par_strategy, this%par_strategy)
      if(.not.varinfo_valid_option('PES_Flux_Parallelization', this%par_strategy, is_flag = .true.)) &
        call messages_input_error(namespace,'PES_Flux_Parallelization')

    end if
    
    !Sanity check  
    if (bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_SURFACE) /= 0 .and. &
        bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_MOMENTUM) /= 0) then
        call messages_input_error(namespace,'PES_Flux_Parallelization', "Cannot combine pf_surface and pf_momentum")
    end if  

    write(message(1),'(a)') 'Input: [PES_Flux_Parallelization = '
    if (bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_NONE) /= 0) then
      write(message(1),'(a,x,a)') trim(message(1)), 'pf_none '
    end if
    if (bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_TIME) /= 0) then
      write(message(1),'(a,x,a)') trim(message(1)), 'pf_time '
    end if
    if (bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_MOMENTUM) /= 0) then
      write(message(1),'(a,x,a)') trim(message(1)), 'pf_momentum'
    end if
    if (bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_SURFACE) /= 0) then
      write(message(1),'(a,x,a)') trim(message(1)), 'pf_surface'
    end if
    write(message(1),'(2a)') trim(message(1)), ']'
    call messages_info(1)  

    

    ! distribute the surface points on nodes,
    ! since mesh domains may have different numbers of surface points.
    if (bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_SURFACE) /= 0) then
      
      call pes_flux_distribute(1, this%nsrfcpnts, this%nsrfcpnts_start, this%nsrfcpnts_end, mesh%mpi_grp%comm)
      
      if (this%surf_shape /= PES_SPHERICAL) call pes_flux_distribute_facepnts_cub(this, mesh)
      
!   Keep this because is useful for debug but not enough to bother having it polished out
!       if(debug%info) then
! #if defined(HAVE_MPI)
!         call MPI_Barrier(mpi_world%comm, mpi_err)
!         write(*,*) &
!           'Debug: surface points on node ', mpi_world%rank, ' : ', this%nsrfcpnts_start, this%nsrfcpnts_end
!         call MPI_Barrier(mpi_world%comm, mpi_err)
! #endif
!       end if

    else
      this%nsrfcpnts_start = 1
      this%nsrfcpnts_end   = this%nsrfcpnts

    end if

!   Keep this because is useful for debug but not enough to bother having it polished out        
!     if(debug%info .and. mpi_grp_is_root(mpi_world)) then
!       do isp = 1, this%nsrfcpnts
!         write(223,*) isp, this%rcoords(:, isp), this%srfcnrml(:,isp)
!       end do
!       if(this%nsrfcpnts > 0) flush(223)
!     end if

    ! Generate the momentum space mesh grid
    call pes_flux_reciprocal_mesh_gen(this, namespace, space, mesh%sb, st, hm%kpoints, mesh%mpi_grp%comm)


    !%Variable PES_Flux_UseSymmetries
    !%Type logical
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Use surface and momentum grid symmetries to speed up calculation and 
    !% lower memory footprint. 
    !% By default available only when the surface shape matches the grid symmetry i.e.:
    !% PES_Flux_Shape = m_cub or m_pln and PES_Flux_Momenutum_Grid = m_cartesian
    !% or 
    !% PES_Flux_Shape = m_sph and PES_Flux_Momenutum_Grid = m_polar
    !%End
    use_symmetry = .false.
    if ((this%surf_shape == PES_CUBIC .or. this%surf_shape == PES_PLANE) &
         .and. this%kgrid == PES_CARTESIAN .and. mdim == 3) use_symmetry = .true.
    if (this%surf_shape == PES_SPHERICAL .and. this%kgrid == PES_POLAR) use_symmetry = .true.
    call parse_variable(namespace, 'PES_Flux_UseSymmetries', use_symmetry, this%use_symmetry)
    call messages_print_var_value(stdout, 'PES_Flux_UseSymmetries', this%use_symmetry)
    
    

    ! get the values of the spherical harmonics for the surface points for PES_SPHERICAL
    if(this%surf_shape == PES_SPHERICAL) then
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

    else 
      
      call pes_flux_integrate_cub_tabulate(this, space, mesh, st, hm%kpoints)

    end if


    !%Variable PES_Flux_RuntimeOutput
    !%Type logical
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Write output in ascii format at runtime. 
    !%  
    !%End
    call parse_variable(namespace, 'PES_Flux_RuntimeOutput', .false., this%runtime_output)
    call messages_print_var_value(stdout, 'PES_Flux_RuntimeOutput', this%runtime_output)



    ! -----------------------------------------------------------------
    ! Options for time integration 
    ! -----------------------------------------------------------------
    this%max_iter  = max_iter
    this%save_iter = save_iter
    this%itstep    = 0
    this%tdsteps = 1
    
    if (mesh%parallel_in_domains .and. bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_TIME) /= 0) then
      this%tdsteps = mesh%mpi_grp%size
    end if
   
    ! -----------------------------------------------------------------
    ! Allocations 
    ! -----------------------------------------------------------------

    SAFE_ALLOCATE(this%wf(stst:stend, 1:sdim, kptst:kptend, 0:this%nsrfcpnts, 1:this%tdsteps))
    this%wf = M_z0

    SAFE_ALLOCATE(this%gwf(stst:stend, 1:sdim, kptst:kptend, 0:this%nsrfcpnts, 1:this%tdsteps, 1:mdim))
    this%gwf = M_z0

    SAFE_ALLOCATE(this%veca(1:mdim, 1:this%tdsteps))
    this%veca = M_ZERO

    if(this%surf_shape == PES_SPHERICAL) then
      SAFE_ALLOCATE(this%spctramp_sph(stst:stend, 1:sdim, kptst:kptend, 1:this%nk, 1:this%nstepsomegak))
      this%spctramp_sph = M_z0

      SAFE_ALLOCATE(this%conjgphase_prev(1:this%nk, 1:this%nstepsomegak))

    else
      SAFE_ALLOCATE(this%spctramp_cub(stst:stend, 1:sdim, kptst:kptend, 1:this%nkpnts))
      this%spctramp_cub = M_z0

      SAFE_ALLOCATE(this%conjgphase_prev(1:this%nkpnts, kptst:kptend))

    end if

    this%conjgphase_prev = M_z1

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

    if(this%surf_shape == PES_SPHERICAL) then
      SAFE_DEALLOCATE_A(this%kcoords_sph)
      SAFE_DEALLOCATE_A(this%ylm_k)
      SAFE_DEALLOCATE_A(this%j_l)
      SAFE_DEALLOCATE_A(this%ylm_r)
      SAFE_DEALLOCATE_A(this%conjgphase_prev)
      SAFE_DEALLOCATE_A(this%spctramp_sph)
    else
      SAFE_DEALLOCATE_A(this%kcoords_cub)
      SAFE_DEALLOCATE_A(this%conjgphase_prev)
      SAFE_DEALLOCATE_A(this%spctramp_cub)
      
      if(.not. this%surf_interp) then
        SAFE_DEALLOCATE_A(this%srfcpnt)
      end if
      SAFE_DEALLOCATE_A(this%rankmin)
      
      SAFE_DEALLOCATE_A(this%face_idx_range)
      SAFE_DEALLOCATE_A(this%LLr)
      SAFE_DEALLOCATE_A(this%NN)      
      
      SAFE_DEALLOCATE_A(this%expkr)
      SAFE_DEALLOCATE_A(this%expkr_perp)
      SAFE_DEALLOCATE_A(this%bvk_phase)      

      SAFE_DEALLOCATE_A(this%Lkpuvz_inv)      
      
    end if

    SAFE_DEALLOCATE_A(this%klinear)

    SAFE_DEALLOCATE_A(this%srfcnrml)
    SAFE_DEALLOCATE_A(this%rcoords)

    SAFE_DEALLOCATE_A(this%wf)
    SAFE_DEALLOCATE_A(this%gwf)
    SAFE_DEALLOCATE_A(this%veca)

    POP_SUB(pes_flux_end)
  end subroutine pes_flux_end

  ! ---------------------------------------------------------
  subroutine pes_flux_reciprocal_mesh_gen(this, namespace, space, sb, st, kpoints, comm, post)
    type(pes_flux_t),    intent(inout) :: this
    type(namespace_t),   intent(in)    :: namespace
    type(space_t),       intent(in)    :: space
    type(simul_box_t),   intent(in)    :: sb
    type(states_elec_t), intent(in)    :: st
    type(kpoints_t),     intent(in)    :: kpoints
    integer,             intent(in)    :: comm
    logical, optional,   intent(in)    :: post !< only fill the data needed for postprocessing  

    integer           :: mdim, pdim
    integer           :: kptst, kptend  
    integer           :: ikp, ikpt
    integer           :: ll, mm, idim, idir
    integer           :: ikk, ith, iph, iomk,ie, ik1, ik2, ik3, kgrid_block_dim
    FLOAT             :: kmax(1:MAX_DIM), kmin(1:MAX_DIM), kact, thetak, phik
    type(block_t)     :: blk
      
    FLOAT             :: Emin, Emax, DE , kvec(1:3) 
    integer           :: nkp_out, nkmin, nkmax
    
    integer             :: kgrid
    FLOAT, allocatable  :: gpoints(:,:), gpoints_reduced(:,:)
    FLOAT               :: dk(1:3), kpoint(1:3), Dthetak, Dphik
    logical             :: use_enegy_grid, arpes_grid
      
    PUSH_SUB(pes_flux_reciprocal_mesh_gen)  

    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    mdim   = space%dim
    pdim   = space%periodic_dim

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
    select case (this%surf_shape)
    case (PES_SPHERICAL)
      kgrid = PES_POLAR
    case (PES_PLANE)
      kgrid = PES_CARTESIAN
    case (PES_CUBIC)
      kgrid = PES_CARTESIAN
      if (mdim == 1)  kgrid = PES_POLAR

    case default
      ASSERT(.false.)
              
    end select
        
    call parse_variable(namespace, 'PES_Flux_Momenutum_Grid', kgrid, this%kgrid)
    if(.not.varinfo_valid_option('PES_Flux_Momenutum_Grid', this%kgrid, is_flag = .true.)) &
      call messages_input_error(namespace,'PES_Flux_Momenutum_Grid')
    call messages_print_var_option(stdout, 'PES_Flux_Momenutum_Grid', this%kgrid)


    
    !Check availability of the calculation requested
    if (this%surf_shape == PES_SPHERICAL) then
      if (this%kgrid == PES_CARTESIAN) then
        call messages_not_implemented('Cartesian momentum grid with a spherical surface')
      end if
      if (space%is_periodic()) then
        call messages_not_implemented('Spherical surface flux for periodic systems')
      end if
      if (mdim == 1) then
        call messages_not_implemented('Spherical surface flux for one-dimensional systems')
      end if
    end if
    
    if (this%surf_shape == PES_CUBIC) then
      if (space%is_periodic()) then
        call messages_not_implemented('Use of cubic surface for periodic systems (use pln)')                
      end if
    end if
    
    

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
      call messages_write(units_from_atomic(units_out%energy, Emin),fmt = '(f8.3)')
      call messages_write(", ")
      call messages_write(units_from_atomic(units_out%energy, Emax), fmt = '(f8.3)')
      call messages_write(", ")
      call messages_write(units_from_atomic(units_out%energy, DE), fmt = '(e9.2)')
      call messages_write(")")
      call messages_info()
    
      kmax(1:mdim) = sqrt(M_TWO*Emax)
      kmin(1:mdim) = sqrt(M_TWO*Emin)
      this%dk = sqrt(M_TWO*DE)

    end if 

    kgrid_block_dim = 1
    !ugly hack (since this input variable is properly handled below) but effective 
    call parse_variable(namespace, 'PES_Flux_ARPES_grid', .false., this%arpes_grid)
    if (.not. use_enegy_grid .or. this%arpes_grid) then
      
      !%Variable PES_Flux_Kmax
      !%Type float
      !%Default 1.0
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% The maximum value of the photoelectron momentum.
      !% For cartesian momentum grids one can specify a value different 
      !% for cartesian direction using a block input.
      !%End      
      if(parse_block(namespace, 'PES_Flux_Kmax', blk) == 0) then
        if (this%kgrid == PES_CARTESIAN) then
          do idim = 1, mdim 
            call parse_block_float(blk, 0, idim - 1, kmax(idim))
          end do
          kgrid_block_dim = mdim
        else
          message(1) = 'Wrong block format for PES_Flux_Kmax and non-cartesian grid'
          call messages_fatal(1, namespace=namespace)
        end if
      else 
        call parse_variable(namespace, 'PES_Flux_Kmax', M_ONE, kmax(1))
        kmax(:)=kmax(1)
      end if

      !%Variable PES_Flux_Kmin
      !%Type float
      !%Default 0.0
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% The minimum value of the photoelectron momentum.
      !% For cartesian momentum grids one can specify a value different 
      !% for cartesian direction using a block input.
      !%End
      if(parse_block(namespace, 'PES_Flux_Kmin', blk) == 0) then
        if (this%kgrid == PES_CARTESIAN) then
          do idim = 1, mdim 
            call parse_block_float(blk, 0, idim - 1, kmin(idim))
          end do
          kgrid_block_dim = mdim
        else
          message(1) = 'Wrong block format for PES_Flux_Kmin and non-cartesian grid'
          call messages_fatal(1, namespace=namespace)
        end if
      else 
        call parse_variable(namespace, 'PES_Flux_Kmin', M_ZERO, kmin(1))
        kmin(:)=kmin(1)
      end if


      !%Variable PES_Flux_DeltaK
      !%Type float
      !%Default 0.02
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Spacing of the the photoelectron momentum grid.
      !%End
      call parse_variable(namespace, 'PES_Flux_DeltaK', CNST(0.02), this%dk)
      if(this%dk <= M_ZERO) call messages_input_error(namespace,'PES_Flux_DeltaK')

    end if

    do idim = 1, kgrid_block_dim   
      if (kgrid_block_dim == 1) then
        call messages_write("Momentum linear grid (Pmin, Pmax, DP) [") 
      else 
        call messages_write("Momentum linear grid (Pmin, Pmax, DP) "//index2axis(idim)//"-axis [") 
      end if
      call messages_write(trim(units_abbrev(units_out%mass*units_out%velocity)))
      call messages_write("]:  (")
      call messages_write(kmin(idim),fmt = '(f8.3)')
      call messages_write(", ")
      call messages_write(kmax(idim), fmt = '(f8.3)')
      call messages_write(", ")
      call messages_write(this%dk, fmt = '(e9.2)')
      call messages_write(")")
      call messages_info()
    end do


    
!     if (this%surf_shape == PES_SPHERICAL .or. this%surf_shape == PES_CUBIC) then
    if (this%kgrid == PES_POLAR) then  
      ! -----------------------------------------------------------------
      ! Setting up k-mesh
      ! 1D = 2 points, 2D = polar coordinates, 3D = spherical coordinates
      ! -----------------------------------------------------------------
      
      
      !%Variable PES_Flux_ThetaK
      !%Type block
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Define the grid points on theta (<math>0 \le \theta \le \pi</math>) when
      !% using a spherical grid in momentum. 
      !% The block defines the maximum and minimum values of theta and the number of 
      !% of points for the discretization.
      !%
      !% <tt>%PES_Flux_ThetaK
      !% <br>&nbsp;&nbsp; theta_min | theta_max  | npoints
      !% <br>%
      !% </tt>
      !%
      !% By default theta_min=0, theta_max = pi, npoints = 45.
      !%End
      this%nstepsthetak = 45
      this%thetak_rng(1:2) = (/M_ZERO, M_PI/)
      if(parse_block(namespace, 'PES_Flux_ThetaK', blk) == 0) then
        call parse_block_float(blk, 0, 0, this%thetak_rng(1))
        call parse_block_float(blk, 0, 1, this%thetak_rng(2))
        call parse_block_integer(blk, 0, 2, this%nstepsthetak)
        call parse_block_end(blk)
        do idim = 1,2
          if (this%thetak_rng(idim) < M_ZERO .or. this%thetak_rng(idim) > M_PI) &
             call messages_input_error(namespace,'PES_Flux_ThetaK')              
        end do 
        if(this%nstepsthetak < 0) call messages_input_error(namespace,'PES_Flux_ThetaK')
      else 

        !%Variable PES_Flux_StepsThetaK
        !%Type integer
        !%Default 45
        !%Section Time-Dependent::PhotoElectronSpectrum
        !%Description
        !% Number of steps in <math>\theta</math> (<math>0 \le \theta \le \pi</math>) for the spherical grid in k.
        !%End
        call parse_variable(namespace, 'PES_Flux_StepsThetaK', 45, this%nstepsthetak)
        if(this%nstepsthetak < 0) call messages_input_error(namespace,'PES_Flux_StepsThetaK')

        ! should do this at some point 
  !       call messages_obsolete_variable(namespace, 'PES_Flux_StepsThetaK')
      end if
      
    

      !%Variable PES_Flux_PhiK
      !%Type block
      !%Section Time-Dependent::PhotoElectronSpectrum
      !%Description
      !% Define the grid points on theta (<math>0 \le \theta \le 2\pi</math>) when
      !% using a spherical grid in momentum. 
      !% The block defines the maximum and minimum values of theta and the number of 
      !% of points for the discretization.
      !%
      !% <tt>%PES_Flux_PhiK
      !% <br>&nbsp;&nbsp; theta_min | theta_max  | npoints
      !% <br>%
      !% </tt>
      !%
      !% By default theta_min=0, theta_max = pi, npoints = 90.
      !%End
      this%nstepsphik = 90
      this%phik_rng(1:2) = (/M_ZERO, 2 * M_PI/)
      if (mdim == 1) this%phik_rng(1:2) = (/M_PI, M_ZERO/)
      if(parse_block(namespace, 'PES_Flux_PhiK', blk) == 0) then
        call parse_block_float(blk, 0, 0, this%phik_rng(1))
        call parse_block_float(blk, 0, 1, this%phik_rng(2))
        call parse_block_integer(blk, 0, 2, this%nstepsphik)
        call parse_block_end(blk)
        do idim = 1,2
          if (this%phik_rng(idim) < M_ZERO .or. this%phik_rng(idim) > M_TWO * M_PI) &
             call messages_input_error(namespace,'PES_Flux_PhiK')              
        end do 
        if(this%nstepsphik < 0) call messages_input_error(namespace,'PES_Flux_PhiK')

      else

        !%Variable PES_Flux_StepsPhiK
        !%Type integer
        !%Default 90
        !%Section Time-Dependent::PhotoElectronSpectrum
        !%Description
        !% Number of steps in <math>\phi</math> (<math>0 \le \phi \le 2 \pi</math>) for the spherical grid in k.
        !%End
        call parse_variable(namespace, 'PES_Flux_StepsPhiK', 90, this%nstepsphik)
        if(this%nstepsphik < 0) call messages_input_error(namespace,'PES_Flux_StepsPhiK')
        if(this%nstepsphik == 0) this%nstepsphik = 1
      end if

      
      Dthetak  = M_ZERO
      if (mdim ==3)  Dthetak = abs(this%thetak_rng(2) - this%thetak_rng(1)) / (this%nstepsthetak)
      Dphik = abs(this%phik_rng(2) - this%phik_rng(1)) / (this%nstepsphik)


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
        ! count the omegak samples
        iomk = 0
        do ith = 0, this%nstepsthetak
          thetak = ith * Dthetak + this%thetak_rng(1)
          do iph = 0, this%nstepsphik - 1
            phik = iph * Dphik + this%phik_rng(1)
            iomk = iomk + 1
            if(thetak < M_EPSILON .or. abs(thetak-M_PI) < M_EPSILON) exit
          end do
        end do
        this%nstepsomegak  = iomk

      case default
        ASSERT(.false.)
        
      end select

      write(message(1),'(a)') "Polar momentum grid:"
      call messages_info(1)
      if(mdim == 3)  then
        write(message(1),'(a,f12.6,a,f12.6,a, i6)') & 
            "  Theta = (", this%thetak_rng(1), ",",this%thetak_rng(2), &
             "); n = ", this%nstepsthetak
        call messages_info(1)
      end if
      write(message(1),'(a,f12.6,a,f12.6,a, i6)') & 
          "  Phi   = (", this%phik_rng(1), ",",this%phik_rng(2), &
           "); n = ", this%nstepsphik
      call messages_info(1)



      if(use_enegy_grid) then
        this%nk     = nint(abs(Emax - Emin) / DE)
      else  
        this%nk     = nint(abs(kmax(1) - kmin(1)) / this%dk)
      end if
      this%nkpnts = this%nstepsomegak * this%nk

      this%ll(1)      = this%nk  
      this%ll(2)      = this%nstepsphik
      this%ll(3)      = this%nstepsthetak !- 1
      this%ll(mdim+1:3) = 1
      
      SAFE_ALLOCATE(this%klinear(1:this%nk, 1))
      

    else 
      ! Cartesian
      
      dk(1:mdim) = M_ONE/kpoints%nik_axis(1:mdim)

      this%arpes_grid = .false.
      if (space%is_periodic()) then
        !%Variable PES_Flux_ARPES_grid
        !%Type logical
        !%Section Time-Dependent::PhotoElectronSpectrum
        !%Description
        !% Use a curvilinear momentum space grid that compensates the transformation 
        !% used to obtain ARPES. With this choice ARPES data is laid out on a Cartesian
        !% regular grid.
        !% By default true when <tt>PES_Flux_Shape = pln</tt> and a <tt>KPointsPath</tt>
        !% is specified.
        !%End
        arpes_grid = kpoints%have_zero_weight_path()
        call parse_variable(namespace, 'PES_Flux_ARPES_grid', arpes_grid, this%arpes_grid)
        call messages_print_var_value(stdout, "PES_Flux_ARPES_grid", this%arpes_grid)       
      end if      
                
      


      this%ll(:) = 1

      if (kpoints%have_zero_weight_path()) then

        if (this%arpes_grid) then
          nkmax = floor(Emax / DE)
          nkmin = floor(Emin / DE)

        else 
          nkmax = floor(kmax(mdim) / this%dk)
          nkmin = floor(kmin(mdim) / this%dk)
    
        end if

        this%ll(mdim) = abs(nkmax - nkmin) + 1
      
        this%nk = this%ll(mdim) 
      
        
      else

        if (.not. this%arpes_grid) then
          this%ll(1:mdim) = floor(abs(kmax(1:mdim) - kmin(1:mdim))/this%dk) + 1
          this%nk = maxval(this%ll(1:mdim)) 
          
        else
          
          nkmax = floor(Emax / DE)
          nkmin = floor(Emin / DE)
          this%nk = abs(nkmax - nkmin) + 1
          
          this%ll(1:pdim) = floor(abs(kmax(1:pdim) - kmin(1:pdim))/this%dk) + 1
          this%ll(mdim) = this%nk
          
        end if

        SAFE_ALLOCATE(this%klinear(1:maxval(this%ll(1:mdim)), 1:mdim))
        this%klinear = M_ZERO

      end if      
      
      ! Total number of points
      this%nkpnts = product(this%ll(1:mdim))
      

    end if    




  
    this%parallel_in_momentum = .false.

    ! Create the grid
    select case (this%kgrid)

    case (PES_POLAR)
    
      if(use_enegy_grid) then
        do ie = 1, this%nk  
          this%klinear(ie, 1) = sqrt(M_TWO * (ie * DE + Emin))
        end do
      else  
        do ikk = 1, this%nk  
          this%klinear(ikk, 1) = ikk * this%dk + kmin(1)
        end do
      end if
      
      
      if (this%surf_shape == PES_SPHERICAL) then

        if(optional_default(post, .false.)) then 
          POP_SUB(pes_flux_reciprocal_mesh_gen)
          return 
        end if
    
        ! we split the k-mesh in radial & angular part
        call pes_flux_distribute(1, this%nk, this%nk_start, this%nk_end, comm)
        if((this%nk_end - this%nk_start + 1) < this%nk) this%parallel_in_momentum = .true.
        call pes_flux_distribute(1, this%nstepsomegak, this%nstepsomegak_start, this%nstepsomegak_end, comm)

!   Keep this because is useful for debug but not enough to bother having it polished out        
!         if(debug%info) then
! #if defined(HAVE_MPI)
!           call MPI_Barrier(mpi_world%comm, mpi_err)
!           write(*,*) &
!             'Debug: momentum points on node ', mpi_world%rank, ' : ', this%nk_start, this%nk_end
!           call MPI_Barrier(mpi_world%comm, mpi_err)
!           write(*,*) &
!             'Debug: momentum directions on node ', mpi_world%rank, ' : ', this%nstepsomegak_start, this%nstepsomegak_end
!           call MPI_Barrier(mpi_world%comm, mpi_err)
! #endif
!         end if
        SAFE_ALLOCATE(this%j_l(0:this%lmax, this%nk_start:this%nk_end))
        this%j_l = M_ZERO

        SAFE_ALLOCATE(this%kcoords_sph(1:3, this%nk_start:this%nk_end, 1:this%nstepsomegak))
        this%kcoords_sph = M_ZERO

        SAFE_ALLOCATE(this%ylm_k(0:this%lmax, - this%lmax:this%lmax, this%nstepsomegak_start:this%nstepsomegak_end))
        this%ylm_k = M_z0

        ! spherical harmonics & kcoords_sph
        iomk = 0
        do ith = 0, this%nstepsthetak
          thetak = ith * Dthetak + this%thetak_rng(1)
          do iph = 0, this%nstepsphik - 1
            phik = iph * Dphik + this%phik_rng(1)
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
            if(thetak < M_EPSILON .or. abs(thetak-M_PI) < M_EPSILON) exit
          end do
        end do

        if(this%nk_start > 0) then
          ! Bessel functions & kcoords_sph
          do ikk = this%nk_start, this%nk_end
            kact = this%klinear(ikk,1)
            do ll = 0, this%lmax
              this%j_l(ll, ikk) = loct_sph_bessel(ll, kact * this%radius) * &
                                  M_TWO * M_PI / (M_TWO * M_PI)**M_THREE / M_TWO
            end do
            this%kcoords_sph(:, ikk, :) = kact * this%kcoords_sph(:, ikk, :)
          end do
        end if

      else 
        !planar or cubic surface
        

        ! no distribution
        this%nkpnts_start = 1 
        this%nkpnts_end   = this%nkpnts


        SAFE_ALLOCATE(this%kcoords_cub(1:mdim, this%nkpnts_start:this%nkpnts_end, 1)) 
        
        this%kcoords_cub = M_ZERO

        if (mdim == 1) then

          ikp = 0
          do ikk = -this%nk, this%nk
            if(ikk == 0) cycle
            ikp = ikp + 1
            kact = sign(1,ikk) * this%klinear(abs(ikk), 1)
            this%kcoords_cub(1, ikp, 1) = kact
          end do

        else !mdim = 2,3 
          thetak = M_PI / M_TWO

          ikp = 0            
          do ikk = 1, this%nk
            do ith = 0, this%nstepsthetak
              if(mdim == 3) thetak = ith * Dthetak + this%thetak_rng(1)
              do iph = 0, this%nstepsphik - 1
                ikp = ikp + 1
                phik = iph * Dphik + this%phik_rng(1)
                kact = this%klinear(ikk,1)
                              this%kcoords_cub(1, ikp,1) = kact * cos(phik) * sin(thetak) 
                              this%kcoords_cub(2, ikp,1) = kact * sin(phik) * sin(thetak) 
                if(mdim == 3) this%kcoords_cub(3, ikp,1) = kact * cos(thetak)
                if(mdim == 3 .and. (thetak < M_EPSILON .or. abs(thetak-M_PI) < M_EPSILON)) exit
              end do
            end do
          end do
            
        end if
      end if

    case (PES_CARTESIAN)

      this%nkpnts_start = 1 
      this%nkpnts_end   = this%nkpnts
      




      if (kpoints%have_zero_weight_path()) then
        ! its a special case because we are generating a different (1D) grid for 
        ! each kpoint and then we combine it in postprocessing

        SAFE_ALLOCATE(this%kcoords_cub(1:mdim, this%nkpnts_start:this%nkpnts_end, kptst:kptend))
                
        nkp_out = 0 
        do ikpt = kptst, kptend
          ikp = 0
          do ikk = nkmin, nkmax
    
            kvec(1:mdim) = - kpoints%get_point(ikpt)             
            call fill_non_periodic_dimension(this)
        
          end do
        end do

      else

        SAFE_ALLOCATE(this%kcoords_cub(1:mdim, this%nkpnts_start:this%nkpnts_end, 1))

        SAFE_ALLOCATE(this%Lkpuvz_inv(this%ll(1), this%ll(2), this%ll(3)))
      
        do idir = 1, mdim
          do ikk = 1, this%ll(idir)  
            this%klinear(ikk, idir) = ikk * this%dk + kmin(idir)
          end do
        end do

        
        if (.not. this%arpes_grid) then
          ! Normal velocity map

          
          ikpt = 1
          ikp = 0
          kvec(:) = M_ZERO
          do ik3 = 1, this%ll(3)
            if(mdim>2) kvec(3) = this%klinear(ik3, 3)
            do ik2 = 1, this%ll(2)
              if(mdim>1) kvec(2) = this%klinear(ik2, 2)
              do ik1 = 1, this%ll(1)
                ikp = ikp + 1 
                kvec(1) = this%klinear(ik1, 1)
                this%kcoords_cub(1:mdim, ikp, ikpt) =  kvec(1:mdim)
                this%Lkpuvz_inv(ik1, ik2, ik3) = ikp

              end do
            end do
          end do

        else ! we want an ARPES-friendly grid layout 
          
          nkp_out = 0 
          ikpt = 1
          ikp = 0
          kvec(:) = M_ZERO
          do ikk = nkmin, nkmax ! this is going to be turned into energy 
    
            ! loop over periodic directions
            select case (pdim)
              case (1)
              
              do ik1 = 1, this%ll(1)
                kvec(1) = this%klinear(ik1,1)
                kvec(1:pdim) = matmul(sb%latt%klattice_primitive(1:pdim, 1:pdim), kvec(1:pdim))
                call fill_non_periodic_dimension(this)               
              end do

              case (2)
          
              do ik2 = 1, this%ll(2)
                do ik1 = 1, this%ll(1)
                  kvec(1:2) = (/this%klinear(ik1,1), this%klinear(ik2,2)/)
                  kvec(1:pdim) = matmul(sb%latt%klattice_primitive(1:pdim, 1:pdim), kvec(1:pdim))
                  call fill_non_periodic_dimension(this)
                  
                  this%Lkpuvz_inv(ik1, ik2, ikk - nkmin + 1) = ikp
                  
                  
                end do
              end do

            case default
              ASSERT(.false.)
              
            end select
              
          end do
        
        
        end if
        
      end if

!   Keep this because is useful for debug but not enough to bother having it polished out    
!       if (debug%info .and. mpi_grp_is_root(mpi_world)) then
!         ! this does not work for parallel in kpoint
!         ! you need to gather kcoords_pln
!         if (kpoints_have_zero_weight_path(sb%kpoints)) then
!           write(229,*) "#   ikpt (kpoint index),   ikp (momentum index),   this%kcoords_cub(1:mdim, ikp, ikpt)"
!           do ikpt = kptst, kptend
!             do ikp = 1, this%nkpnts
!               write(229,*) ikpt, ikp, this%kcoords_cub(1:mdim, ikp, ikpt)
!             end do
!           end do
!         else
!           write(229,*) "#   ikp (momentum index),   this%kcoords_cub(1:mdim, ikp, ikpt)"
!           do ikp = 1, this%nkpnts
!             write(229,*) ikp, this%kcoords_cub(1:mdim, ikp, 1)
!           end do
!         end if
!         flush(229)
!
!
!         if(mdim == 3) then
!           write(230,*) "#   ik1, ik2, ik3,  this%Lkpuvz_inv(ik1,ik2,ik3) "
!           do ik1 = 1, this%ll(1)
!             do ik2 = 1, this%ll(2)
!               do ik3 = 1, this%ll(3)
!                 write(230,*) ik1, ik2, ik3,  this%Lkpuvz_inv(ik1,ik2,ik3)
!               end do
!             end do
!           end do
!
!           flush(230)
!         end if
!
!       end if
      
      if (this%arpes_grid) then
        call messages_write("Number of points with E<p//^2/2 = ")
        call messages_write(nkp_out)
        call messages_write(" [of ")
        call messages_write(this%nkpnts*kpoints_number(kpoints))
        call messages_write("]")
        call messages_info()
      end if
      
      SAFE_DEALLOCATE_A(gpoints)
      SAFE_DEALLOCATE_A(gpoints_reduced)

    case default
      ASSERT(.false.)

    end select
    
    
    !%Variable PES_Flux_GridTransformMatrix
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Define an optional transformation matrix for the momentum grid.
    !%
    !% <tt>%PES_Flux_GridTransformMatrix
    !% <br>&nbsp;&nbsp; M_11 | M_12  | M_13
    !% <br>&nbsp;&nbsp; M_21 | M_22  | M_23
    !% <br>&nbsp;&nbsp; M_31 | M_32  | M_33
    !% <br>%
    !% </tt>
    !%End
    this%ktransf(:,:) = M_ZERO
    do idim = 1,mdim
      this%ktransf(idim,idim) = M_ONE
    end do
    
    if(parse_block(namespace, 'PES_Flux_GridTransformMatrix', blk) == 0) then
      do idim = 1,mdim
        do idir = 1, mdim
          call parse_block_float(blk, idir - 1, idim - 1, this%ktransf(idim, idir))
        end do
      end do
      call parse_block_end(blk)
      
      write(message(1),'(a)') 'Momentum grid transformation matrix :'
      do idir = 1, space%dim
        write(message(1 + idir),'(9f12.6)') ( this%ktransf(idim, idir), idim = 1, mdim) 
      end do
      call messages_info(1 + mdim)

      
      !Apply the transformation
      if (this%surf_shape == PES_SPHERICAL) then

        iomk = 0
        do ith = 0, this%nstepsthetak
          do iph = 0, this%nstepsphik - 1
            iomk = iomk + 1
            do ikk = this%nk_start, this%nk_end
              kvec(1:mdim) = this%kcoords_sph(1:mdim, ikk, iomk) 
              this%kcoords_sph(1:mdim, ikk, iomk) = matmul(this%ktransf(1:mdim, 1:mdim), kvec(1:mdim)) 
            end do
            if(ith == 0 .or. ith == this%nstepsthetak) exit
          end do
        end do
    
      else !planar or cubic surface
      
        do ikpt = kptst, kptend + 1
          if (ikpt == kptend + 1) then
            kpoint(1:space%dim) = M_ZERO
          else
            kpoint(1:space%dim) = kpoints%get_point(ikpt)
          end if

          do ikp = 1, this%nkpnts
          
            kvec(1:mdim) = this%kcoords_cub(1:mdim, ikp, ikpt) - kpoint(1:mdim) 
            this%kcoords_cub(1:mdim, ikp, ikpt) =  matmul(this%ktransf(1:mdim, 1:mdim), kvec(1:mdim)) &
                                                   + kpoint(1:mdim)
          end do
        end do
      
      end if

      
    end if

    
    
    
    POP_SUB(pes_flux_reciprocal_mesh_gen)
    
  contains 
    
    ! Fill the non-periodic direction
    subroutine fill_non_periodic_dimension(this)
      type(pes_flux_t),   intent(inout) :: this
        
      integer :: sign
      FLOAT   :: kpar(1:pdim), val

      ikp = ikp + 1
      
      sign = 1         
      if (ikk /= 0) sign= ikk / abs(ikk)        
      
      kpar(1:pdim) = kvec(1:pdim) 
      val = abs(ikk) * DE * M_TWO - sum(kpar(1:pdim)**2)
      if (val >= 0) then
        kvec(mdim) =  sign * sqrt(val)
      else  ! if E < p//^2/2
        !FIXME: Should handle the exception doing something smarter than this
        kvec(mdim) = sqrt(val) ! set to NaN
        nkp_out = nkp_out + 1
      end if
      
      this%kcoords_cub(1:mdim, ikp, ikpt) =  kvec(1:mdim)
      
      
    end subroutine fill_non_periodic_dimension
                   
          
    
  end subroutine pes_flux_reciprocal_mesh_gen

  ! ---------------------------------------------------------
  subroutine pes_flux_calc(this, space, mesh, st, gr, hm, iter, dt, stopping)
    type(pes_flux_t),         intent(inout)    :: this
    type(space_t),            intent(in)       :: space
    type(mesh_t),             intent(in)       :: mesh
    type(states_elec_t),      intent(inout)    :: st
    type(grid_t),             intent(in)       :: gr
    type(hamiltonian_elec_t), intent(in)       :: hm
    integer,                  intent(in)       :: iter
    FLOAT,                    intent(in)       :: dt
    logical,                  intent(in)       :: stopping

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim, imdim
    integer            :: isp,ip
    integer            :: il, ip_local
    CMPLX, allocatable :: gpsi(:,:), psi(:)
    CMPLX, allocatable :: interp_values(:)
    logical            :: contains_isp
    FLOAT              :: kpoint(1:3)

    PUSH_SUB(pes_flux_calc)

    if(iter > 0) then

      stst   = st%st_start
      stend  = st%st_end
      kptst  = st%d%kpt%start
      kptend = st%d%kpt%end
      sdim   = st%d%dim
      mdim   = space%dim

      SAFE_ALLOCATE(psi(1:mesh%np_part))
      SAFE_ALLOCATE(gpsi(1:mesh%np_part, 1:mdim))

      if(this%surf_interp) then
        SAFE_ALLOCATE(interp_values(1:this%nsrfcpnts))
      end if

      this%itstep = this%itstep + 1

      ! get and save current laser field
      do il = 1, hm%ext_lasers%no_lasers
        call laser_field(hm%ext_lasers%lasers(il), this%veca(1:mdim, this%itstep), iter*dt)
      end do
      this%veca(:, this%itstep) = - this%veca(:, this%itstep)

!  Ideally one could directly access uniform_vector_potential
!       this%veca(:, this%itstep) = hm%hm_base%uniform_vector_potential(:)

      ! save wavefunctions & gradients
      do ik = kptst, kptend
        do ist = stst, stend
          do isdim = 1, sdim
            call states_elec_get_state(st, mesh, isdim, ist, ik, psi)
            
            if(this%surf_shape == PES_PLANE) then
              ! Apply the phase containing kpoint only
              kpoint(1:mdim) = hm%kpoints%get_point(st%d%get_kpoint_index(ik))

              !$omp parallel do schedule(static)
              do ip = 1, mesh%np_part
                psi(ip) = exp(- M_zI * sum(mesh%x(ip, 1:mdim) * kpoint(1:mdim))) * psi(ip) 
              end do
              !$omp end parallel do
            end if
            
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
                    st%occ(ist, ik) * gpsi(ip_local, 1:mdim) 
                end if
              end do
              if(mesh%parallel_in_domains) then
                call mesh%allreduce(this%wf(ist, isdim, ik, 1:this%nsrfcpnts, this%itstep))
                do imdim = 1, mdim
                  call mesh%allreduce(this%gwf(ist, isdim, ik, 1:this%nsrfcpnts, this%itstep, imdim))
                end do
              end if
            end if
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(gpsi)

      if(this%itstep == this%tdsteps .or. mod(iter, this%save_iter) == 0 .or. iter == this%max_iter .or. stopping) then
        if(this%surf_shape == PES_CUBIC .or. this%surf_shape == PES_PLANE) then
          call pes_flux_integrate_cub(this, space, mesh, st, hm%kpoints, dt)
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
  subroutine pes_flux_integrate_cub_tabulate(this, space, mesh, st, kpoints)
    type(pes_flux_t),    intent(inout) :: this
    type(space_t),       intent(in)    :: space
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: st
    type(kpoints_t),     intent(in)    :: kpoints

    integer            :: kptst, kptend,  mdim
    integer            :: ik, j1, j2, jvec(1:2)
    integer            :: isp, ikp, ikp_start, ikp_end
    integer            :: ik_map
    
    integer            :: fdim, idir,pdim, nfaces, ifc, n_dir
    CMPLX              :: tmp
    FLOAT              :: Jac(1:2,1:2), jdet, kpoint(1:3), vec(1:3), lvec(1:3)
    
    PUSH_SUB(pes_flux_integrate_cub_tabulate) 

    if (kpoints%have_zero_weight_path()) then
      kptst     = st%d%kpt%start
      kptend    = st%d%kpt%end
    else
      kptst     = 1
      kptend    = 1
    end if

    mdim = space%dim
    pdim = space%periodic_dim

    ikp_start = this%nkpnts_start
    ikp_end   = this%nkpnts_end
    
    
    if(this%surf_shape == PES_PLANE) then
      fdim = mdim - 1
      ! This is not general but should work in the specific case where it is relevant
      !i.e. when the system is semiperiodic in <=2 dimensions
      Jac(1:fdim, 1:fdim) = mesh%sb%latt%rlattice_primitive(1:fdim, 1:fdim) !The Jacobian on the surface
      jdet = lalg_determinant(fdim, Jac, preserve_mat = .false.)

!   Keep this because is useful for debug but not enough to bother having it polished out        
!       if(debug%info .and. mpi_grp_is_root(mpi_world)) then
!         print *, "jdet =", jdet
!         print *, "mesh%sb%rlattice_primitive(1:fdim, 1:fdim) ", mesh%sb%rlattice_primitive(1:fdim, 1:fdim)
!       end if
      
      this%srfcnrml(:, 1:this%nsrfcpnts) = this%srfcnrml(:, 1:this%nsrfcpnts) * jdet
    end if
    
    
    
    if (.not. this%use_symmetry .or. kpoints%have_zero_weight_path()) then
      
      SAFE_ALLOCATE(this%expkr(1:this%nsrfcpnts,ikp_start:ikp_end,kptst:kptend,1))

      do ik = kptst, kptend
        do ikp = ikp_start, ikp_end
          do isp = 1, this%nsrfcpnts 
            this%expkr(isp,ikp,ik,1) = exp(- M_zI * sum(this%rcoords(1:mdim,isp) &
                                                * this%kcoords_cub(1:mdim, ikp, ik)) ) &
                                                * (M_TWO * M_PI)**(- mdim / M_TWO)
      
          end do
        end do
      end do


      
    else !do something smart to exploit symmetries

      nfaces = mdim*2
      if(this%surf_shape == PES_PLANE) nfaces = 1 
      
      
      SAFE_ALLOCATE(this%expkr(1:this%nsrfcpnts,maxval(this%ll(1:mdim)),kptst:kptend,1:mdim))
      this%expkr(:,:,:,:) = M_z1

      do ik = kptst, kptend !should only be ik=1
        do idir = 1, mdim
          do ikp = 1, this%ll(idir)
            do isp = 1, this%nsrfcpnts 
              this%expkr(isp,ikp,ik,idir) = exp(- M_zI * this%rcoords(idir,isp) &
                                                  * this%klinear(ikp, idir) ) & 
                                                  * (M_TWO * M_PI)**(-M_ONE / M_TWO)
    
            end do
          end do
        end do
      end do

      SAFE_ALLOCATE(this%expkr_perp(maxval(this%ll(1:mdim)), nfaces))
      this%expkr_perp(:,:) = M_z1

      do ifc = 1, nfaces
        if (this%face_idx_range(ifc, 1) < 0) cycle ! this face have no local surface point
        isp = this%face_idx_range(ifc, 1)
        do idir = 1, mdim
          if(abs(this%srfcnrml(idir, isp)) >= M_EPSILON) n_dir = idir
        end do 

        do ikp = 1, this%ll(n_dir)
          this%expkr_perp(ikp, ifc) = exp(- M_zI * this%rcoords(n_dir,isp) &
                                            * this%klinear(ikp, n_dir) ) & 
                                            * (M_TWO * M_PI)**(- M_ONE / M_TWO)


        end do
      end do
      
    end if
    
    
    if (space%is_periodic()) then
      !Tabulate the Born-von Karman phase 
      SAFE_ALLOCATE(this%bvk_phase(ikp_start:ikp_end,st%d%kpt%start:st%d%kpt%end))

      this%bvk_phase(:,:) = M_z0
      vec(:) = M_ZERO

      do ik = st%d%kpt%start, st%d%kpt%end
        kpoint(1:mdim) = kpoints%get_point(ik)
        if (kpoints%have_zero_weight_path()) then
          ik_map = ik
        else
          ik_map = 1
        end if
        do ikp = ikp_start, ikp_end
          vec(1:pdim) = this%kcoords_cub(1:pdim, ikp, ik_map) + kpoint(1:pdim)
          do j1 = 0, kpoints%nik_axis(1) - 1
            do j2 = 0, kpoints%nik_axis(2) - 1
              jvec(1:2)=(/j1, j2/)
              lvec(1:pdim)=matmul(mesh%sb%latt%rlattice(1:pdim, 1:2), jvec(1:2))
              tmp = sum(lvec(1:pdim) * vec(1:pdim))
              this%bvk_phase(ikp, ik) =  this%bvk_phase(ikp,ik) &
                                     +  exp(M_zI * tmp)
            
            end do
          end do

        end do
      end do
      this%bvk_phase(:,:) = this%bvk_phase(:,:) * M_ONE / product(kpoints%nik_axis(1:pdim))
    

!   Keep this because is useful for debug but not enough to bother having it polished out        
!       if(debug%info .and. mpi_grp_is_root(mpi_world)) then
!         write(225,*) "ik, ikp, this%bvk_phase(ikp,ik)"
!         do ik = st%d%kpt%start, st%d%kpt%end
!           do ikp = ikp_start, ikp_end
!             write(225,*) ik, ikp, this%bvk_phase(ikp,ik)
!           end do
!         end do
!         flush(225)
!       end if
    
    end if
            
    POP_SUB(pes_flux_integrate_cub_tabulate)  
          
  end subroutine pes_flux_integrate_cub_tabulate


  ! ---------------------------------------------------------
  pure function get_ikp(this, ikpu, ikpv, ikpz, n_dir) result(ikp)
    type(pes_flux_t),        intent(in)  :: this
    integer,                 intent(in)  :: ikpu
    integer,                 intent(in)  :: ikpv
    integer,                 intent(in)  :: ikpz
    integer,                 intent(in)  :: n_dir
    integer                              :: ikp
    
    select case (n_dir)
    case (1)
      ikp = this%Lkpuvz_inv(ikpz, ikpu, ikpv)
    case (2)
      ikp = this%Lkpuvz_inv(ikpu, ikpz, ikpv)
    case (3)
      ikp = this%Lkpuvz_inv(ikpu, ikpv, ikpz)
    
    case default  
    ! should die here but cannot use assert in a pure function
      
    end select

  end function get_ikp

  ! ---------------------------------------------------------
  subroutine pes_flux_distribute_facepnts_cub(this, mesh)
    type(pes_flux_t), intent(inout) :: this
    type(mesh_t),        intent(in) :: mesh  

    
    integer :: mdim, nfaces, ifc, ifp_start, ifp_end
      
    PUSH_SUB(pes_flux_distribute_facepnts_cub)

    mdim      = mesh%sb%dim

    nfaces = mdim*2
    if(this%surf_shape == PES_PLANE) nfaces = 1 
    
    do ifc = 1, nfaces
    
      ifp_start = this%face_idx_range(ifc, 1)
      ifp_end   = this%face_idx_range(ifc, 2)
        
        
      if(this%nsrfcpnts_start <= ifp_end) then ! the local domain could be in this face
        
        if(this%nsrfcpnts_start >= ifp_start) then 
          if (this%nsrfcpnts_start <= ifp_end) then
            this%face_idx_range(ifc, 1) =  this%nsrfcpnts_start
          else
            this%face_idx_range(ifc, 1:2) = -1  ! the local domain is not in this face
          end if
        end if
        
        if(this%nsrfcpnts_end   <= ifp_end  ) then 
          if (this%nsrfcpnts_end >= ifp_start) then
            this%face_idx_range(ifc, 2) =  this%nsrfcpnts_end
          else 
             this%face_idx_range(ifc, 1:2) = -1  ! the local domain is not in this face
          end if
        end if
        
      else 

        this%face_idx_range(ifc, 1:2) = -1  ! the local domain is not in this face

      end if
      
!   Keep this because is useful for debug but not enough to bother having it polished out        
!       if(debug%info) then
!         print *, mpi_world%rank, ifc, ifp_start, ifp_end, this%face_idx_range(ifc, 1:2), this%nsrfcpnts_start, this%nsrfcpnts_end
!       end if

    end do
    
    
    POP_SUB(pes_flux_distribute_facepnts_cub)      
  end subroutine pes_flux_distribute_facepnts_cub


  ! ---------------------------------------------------------
  subroutine pes_flux_integrate_cub(this, space, mesh, st, kpoints, dt)
    type(pes_flux_t),    intent(inout) :: this
    type(space_t),       intent(in)    :: space
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(inout) :: st 
    type(kpoints_t),     intent(in)    :: kpoints
    FLOAT,               intent(in)    :: dt

    integer            :: stst, stend, kptst, kptend, sdim, mdim
    integer            :: ist, ik, isdim, imdim, ik_map
    integer            :: ikp, itstep
    integer            :: idir, n_dir, nfaces
    CMPLX, allocatable :: Jk_cub(:,:,:,:), spctramp_cub(:,:,:,:)
    integer            :: ikp_start, ikp_end, isp_start, isp_end
    FLOAT              :: vec, kpoint(1:3)
    integer            :: ifc
    CMPLX, allocatable :: wfpw(:), gwfpw(:)
    CMPLX, allocatable :: phase(:,:),vphase(:,:)
    
    integer            :: tdstep_on_node
    integer            :: nfp
    
    !Symmetry helpers
    integer            :: ikpu, ikpv, ikpz, dir_on_face(1:2)
    CMPLX              :: face_int_gwf, face_int_wf
    
    type(profile_t), save :: prof_init
      
    PUSH_SUB(pes_flux_integrate_cub)

    ! this routine is parallelized over time steps and surface points 
    
    call profiling_in(prof_init, 'pes_flux_integrate_cub') 

    if (debug%info) then
      call messages_write("Debug: calculating pes_flux cub surface integral (accelerated direct expression)")
      call messages_info()
    end if


    tdstep_on_node = 1
    if (bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_TIME) /= 0) then
#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) tdstep_on_node = mesh%mpi_grp%rank + 1
#endif
    end if

    stst      = st%st_start
    stend     = st%st_end
    kptst     = st%d%kpt%start
    kptend    = st%d%kpt%end
    sdim      = st%d%dim
    mdim      = mesh%sb%dim

    ikp_start = this%nkpnts_start
    ikp_end   = this%nkpnts_end



    SAFE_ALLOCATE(Jk_cub(stst:stend, 1:sdim, kptst:kptend, ikp_start:ikp_end))
    SAFE_ALLOCATE(spctramp_cub(stst:stend, 1:sdim, kptst:kptend, ikp_start:ikp_end))
    spctramp_cub = M_z0

    
    SAFE_ALLOCATE(vphase(ikp_start:ikp_end, kptst:kptend))
    SAFE_ALLOCATE(phase(ikp_start:ikp_end, kptst:kptend))
    vphase(:,:) = M_z0
    phase(:,:) = M_z0

    
    nfaces = mdim * 2
    if(this%surf_shape == PES_PLANE) nfaces = 1 

    SAFE_ALLOCATE( wfpw(ikp_start:ikp_end))
    SAFE_ALLOCATE(gwfpw(ikp_start:ikp_end))
          

    do ifc = 1, nfaces
      
      if (this%face_idx_range(ifc, 1)<0) cycle ! this face have no local surface point
      
      
      isp_start = this%face_idx_range(ifc, 1)
      isp_end   = this%face_idx_range(ifc, 2)

  
      nfp = isp_end - isp_start + 1 ! faces can have a different number of points

      wfpw = M_z0
      gwfpw = M_z0

      ! get the directions normal to the surface and parallel to it
      imdim = 1
      do idir = 1, mdim
        if(abs(this%srfcnrml(idir, isp_start)) >= M_EPSILON) then
          n_dir = idir
        else 
          dir_on_face(imdim)=idir
          imdim = imdim + 1
        end if
      end do 


      ! get the previous Volkov phase
      vphase(ikp_start:ikp_end,:) = this%conjgphase_prev(ikp_start:ikp_end,:)
      
      Jk_cub(:, :, :, :) = M_z0

      do itstep = 1, this%itstep

        do ik = kptst, kptend

          if (kpoints%have_zero_weight_path()) then
            ik_map = ik 
          else
            ik_map = 1 
          end if
          
          kpoint(1:mdim) = kpoints%get_point(ik)
          do ikp = ikp_start, ikp_end
            vec = sum((this%kcoords_cub(1:mdim, ikp, ik_map) - this%veca(1:mdim, itstep) / P_c)**2)
            vphase(ikp, ik) = vphase(ikp, ik) * exp(M_zI * vec * dt / M_TWO)


            if (space%is_periodic()) then
              phase(ikp, ik)  = vphase(ikp, ik) *  this%bvk_phase(ikp,ik)
            else 
              phase(ikp, ik)  = vphase(ikp, ik) 
            end if

          end do
          
          if(itstep /= tdstep_on_node) cycle

          do ist = stst, stend
            do isdim = 1, sdim
              
              
              ! calculate the surface integrals
              if (.not. this%use_symmetry .or. kpoints%have_zero_weight_path()) then

                do ikp = ikp_start , ikp_end
            
                  gwfpw(ikp) = &
                    sum(this%gwf(ist, isdim, ik, isp_start:isp_end, itstep, n_dir) &
                      * this%expkr(isp_start:isp_end, ikp, ik_map, 1)              &
                      * this%srfcnrml(n_dir, isp_start:isp_end))


                  wfpw(ikp) = &
                    sum(this%wf(ist, isdim, ik, isp_start:isp_end, itstep)        &
                      * this%expkr(isp_start:isp_end, ikp,ik_map, 1)              &
                      * this%srfcnrml(n_dir, isp_start:isp_end))
                end do 

              else 
                !$omp parallel do private (ikpv,ikpz,face_int_gwf,face_int_wf) shared(gwfpw, wfpw)
                do ikpu = 1, this%ll(dir_on_face(1))
                  do ikpv = 1, this%ll(dir_on_face(2))
              
                  
                    face_int_gwf = sum(this%gwf(ist, isdim, ik, isp_start:isp_end, itstep, n_dir) &
                                 * this%expkr(isp_start:isp_end, ikpu, ik_map, dir_on_face(1)) &
                                 * this%expkr(isp_start:isp_end, ikpv, ik_map, dir_on_face(2)) &
                                 * this%srfcnrml(n_dir, isp_start:isp_end))
                                 
                    face_int_wf  = sum(this%wf(ist, isdim, ik, isp_start:isp_end, itstep) &
                                 * this%expkr(isp_start:isp_end, ikpu, ik_map, dir_on_face(1)) &
                                 * this%expkr(isp_start:isp_end, ikpv, ik_map, dir_on_face(2)) &
                                 * this%srfcnrml(n_dir, isp_start:isp_end))
                                 
                    do ikpz = 1, this%ll(n_dir)

                       gwfpw(get_ikp(this, ikpu, ikpv, ikpz, n_dir)) = face_int_gwf &
                                                                     * this%expkr_perp(ikpz, ifc)  
                       wfpw( get_ikp(this, ikpu, ikpv, ikpz, n_dir)) = face_int_wf &
                                                                     * this%expkr_perp(ikpz, ifc)  

                    end do                 
                                        
                  end do
                end do
 
              
              end if
                            

                ! Sum it up  
                Jk_cub(ist, isdim, ik, ikp_start:ikp_end) = Jk_cub(ist, isdim, ik, ikp_start:ikp_end) +  &
                  phase(ikp_start:ikp_end, ik) * ( wfpw(ikp_start:ikp_end) * &
                     (M_TWO * this%veca(n_dir, itstep) / P_c  - this%kcoords_cub(n_dir, ikp_start:ikp_end, ik_map)) + &
                      M_zI * gwfpw(ikp_start:ikp_end) )


          end do ! isdim
        end do ! ist     
      end do ! is

    end do !istep
    
    
    spctramp_cub(:,:,:,:) = spctramp_cub(:,:,:,:) + Jk_cub(:,:,:,:) * M_HALF 
        
        
    end do ! face loop 


    this%conjgphase_prev(ikp_start:ikp_end, :) = vphase(ikp_start:ikp_end, :)


    if(mesh%parallel_in_domains .and.(     bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_TIME)    /= 0 &
                                      .or. bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_SURFACE) /= 0 )) then
      call mesh%allreduce(spctramp_cub)
    end if



    this%spctramp_cub = this%spctramp_cub + spctramp_cub * dt

    SAFE_DEALLOCATE_A(gwfpw)
    SAFE_DEALLOCATE_A(wfpw)

    SAFE_DEALLOCATE_A(Jk_cub)
    SAFE_DEALLOCATE_A(spctramp_cub)
    SAFE_DEALLOCATE_A(phase)
    SAFE_DEALLOCATE_A(vphase)

    call profiling_out(prof_init)

    POP_SUB(pes_flux_integrate_cub)
  end subroutine pes_flux_integrate_cub



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

    ! this routine is parallelized over time steps and surface points

    if (debug%info) then
      call messages_write("Debug: calculating pes_flux sph surface integral")
      call messages_info()
    end if

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

            if(mesh%parallel_in_domains .and. bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_SURFACE) /= 0) then
              call mesh%allreduce(s1_act)
              call mesh%allreduce(s2_act)
            end if

            if(itstep == tdstep_on_node) then
              s1_node(ist, isdim, ik, :, :, :) = s1_act(:, :, :)
              s2_node(ist, isdim, ik, :, :)    = s2_act(:, :)
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
      phase_act(ikk_start:ikk_end, :) = this%conjgphase_prev(ikk_start:ikk_end, :)
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
              call MPI_Bcast(s1_act, (lmax + 1) * (2 * lmax + 1) * 3, MPI_CMPLX, itstep - 1, mesh%mpi_grp%comm, mpi_err)
              call MPI_Bcast(s2_act, (lmax + 1) * (2 * lmax + 1), MPI_CMPLX, itstep - 1, mesh%mpi_grp%comm, mpi_err)
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
                call mesh%allreduce(integ11_t)
                call mesh%allreduce(integ21_t)
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
    this%conjgphase_prev = M_z0

    if(ikk_start > 0) then
      this%conjgphase_prev(ikk_start:ikk_end, :) = phase_act(ikk_start:ikk_end, :)
    end if
    SAFE_DEALLOCATE_A(phase_act)

    if(mesh%parallel_in_domains .and. bitand(this%par_strategy, OPTION__PES_FLUX_PARALLELIZATION__PF_TIME) /= 0) then
      call mesh%allreduce(this%conjgphase_prev)
      call mesh%allreduce(spctramp_sph)
    end if

    this%spctramp_sph(:, :, :, 1:this%nk, :) = this%spctramp_sph(:, :, :, 1:this%nk, :) & 
                                             + spctramp_sph(:, :, :, 1:this%nk, :) * dt
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
    FLOAT                 :: xx(1:mesh%sb%dim), dd, area, dS(1:MAX_DIM, 1:2), factor
    integer               :: mdim, imdim, idir, isp, pm,nface, idim, ndir, iu,iv, iuv(1:2)
    integer               :: ip_global, npface
    integer               :: rankmin, nsurfaces
    logical               :: in_ab
    integer               :: ip_local, NN(1:MAX_DIM, 1:2), idx(1:MAX_DIM, 1:2) 
    integer               :: isp_end, isp_start, ifc, n_dir, nfaces, mindim
    FLOAT                 :: RSmax(1:2),RSmin(1:2),RS(1:2), dRS(1:2), weight


    PUSH_SUB(pes_flux_getcube)

    ! this routine works on absolute coordinates
    

    mdim = mesh%sb%dim
    
    

    SAFE_ALLOCATE(this%face_idx_range(1:mdim * 2, 1:2))
    SAFE_ALLOCATE(this%LLr(mdim, 1:2))
    SAFE_ALLOCATE(this%NN(mdim, 1:2))

    this%face_idx_range(:, :) = 0
    this%LLr(:, :) = M_ZERO
    this%NN(:, :) = 1
    NN(:, :) = 1

    if (this%surf_interp) then
      ! Create a surface with points not on the mesh
      
      idx(:,:) = 0 
      
      mindim  = 1 
      factor = M_TWO
      if(this%surf_shape == PES_PLANE) then
        mindim = mdim  ! We only have two planes along the non periodic dimension 
        ! this is due to the fact that for semiperiodic systems we multiply border by two to prenvet the creation of surfaces at the edges
        factor = M_ONE 
      end if
      
      
      this%nsrfcpnts = 0 
      
      do ndir = mdim, mindim, -1
        area = M_ONE
        do idir=1, mdim
          if (idir == ndir) cycle
          area = area * border(idir) * factor
        end do
        
        npface = int(fc_ptdens * area) ! number of points on the face

        idim = 1 
        do idir=1, mdim
          if (idir == ndir) cycle
          NN(ndir, idim) = int( sqrt(npface * (border(idir) * factor)**2 / area) )
          dS(ndir, idim) = border(idir) * factor / NN(ndir, idim)
          this%LLr(ndir, idim) = NN(ndir, idim) * dS(ndir, idim)
          idx(ndir, idim) = idir
          idim = idim + 1
        end do
        this%nsrfcpnts = this%nsrfcpnts + 2 * product(NN(ndir, 1:mdim - 1))
      end do
      
      this%NN(1:mdim, :)  = NN(1:mdim, :)
      
      
      ASSERT(this%nsrfcpnts > 0) !at this point should be satisfied otherwise the point density is way too small

      
      SAFE_ALLOCATE(this%srfcnrml(1:mdim, 0:this%nsrfcpnts))
      SAFE_ALLOCATE(this%rcoords(1:mdim, 0:this%nsrfcpnts))
    
      this%srfcnrml(:, :) = M_ZERO
      this%rcoords(:, :)  = M_ZERO
    

      isp = 1  
      ifc = 1
      do ndir = mdim, mindim, -1

        !Up face
        this%face_idx_range(ifc,1) = isp
        do iu = 1, NN(ndir,1)
          do iv = 1, NN(ndir,2)
            this%rcoords(ndir, isp)  =  border(ndir)
            this%srfcnrml(ndir, isp) =  product(dS(ndir,1:mdim-1))             
            iuv =(/iu, iv/)
            do idim = 1, mdim-1
              this%rcoords(idx(ndir, idim), isp) = (-NN(ndir, idim) * M_HALF - M_HALF + iuv(idim)) * dS(ndir, idim)
            end do 
            isp = isp + 1
          end do
        end do
        this%face_idx_range(ifc, 2) = isp-1

        ifc = ifc + 1

        !Down face
        this%face_idx_range(ifc, 1) = isp
        do iu = 1, NN(ndir, 1)
          do iv = 1, NN(ndir, 2)
            this%rcoords(ndir, isp)  = -border(ndir)
            this%srfcnrml(ndir, isp) = -product(dS(ndir, 1:mdim - 1)) 
            iuv =(/iu, iv/)              
            do idim = 1, mdim-1
              this%rcoords(idx(ndir, idim), isp) = (-NN(ndir, idim) * M_HALF -M_HALF + iuv(idim)) * dS(ndir, idim)
            end do              
            isp = isp + 1
          end do
        end do
        this%face_idx_range(ifc, 2) = isp-1

        ifc = ifc + 1 
      end do
      
      do isp = 1, this%nsrfcpnts
        xx(1:mdim) = matmul(mesh%sb%latt%rlattice_primitive(1:mdim, 1:mdim),this%rcoords(1:mdim, isp))
        this%rcoords(1:mdim, isp) = xx(1:mdim)
      end do
      
    else
      ! Surface points are on the mesh
      
      nfaces = mdim * 2
      if(this%surf_shape == PES_PLANE) nfaces = 2 

      in_ab = .false.

      SAFE_ALLOCATE(which_surface(1:mesh%np_global))
      which_surface = 0

      ! get the surface points
      this%nsrfcpnts = 0
      do ip_local = 1, mesh%np
        ip_global = mesh_local2global(mesh, ip_local)
      
        nsurfaces = 0

        xx(1:mdim) = mesh%x(ip_local, 1:mdim) - offset(1:mdim)

        ! eventually check whether we are in absorbing zone
        select case(hm%bc%abtype)
        case(MASK_ABSORBING)
          in_ab = (hm%bc%mf(ip_local) /= M_ONE)
        case(IMAGINARY_ABSORBING)
          in_ab = (hm%bc%mf(ip_local) /= M_ZERO)
        case default
          ASSERT(.false.)          
        end select

        ! check whether the point is inside the cube
        if(all(abs(xx(1:mdim)) <= border(1:mdim)) .and. .not. in_ab) then
          ! check whether the point is close to any border
          do imdim = 1, mdim
            if(this%surf_shape == PES_PLANE) then
              dd = border(imdim) - xx(imdim)
            else
              dd = border(imdim) - abs(xx(imdim))
            end if
            if(dd < mesh%spacing(imdim) / M_TWO) then
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
        call mesh%allreduce(this%nsrfcpnts)
        call mesh%allreduce(which_surface)
      end if

      SAFE_ALLOCATE(this%srfcpnt(1:this%nsrfcpnts))
      SAFE_ALLOCATE(this%srfcnrml(1:mdim, 0:this%nsrfcpnts))
      SAFE_ALLOCATE(this%rcoords(1:mdim, 0:this%nsrfcpnts))
      SAFE_ALLOCATE(this%rankmin(1:this%nsrfcpnts))

      this%srfcnrml = M_ZERO
      this%rcoords  = M_ZERO

      this%face_idx_range(:,:) = 0
    
      isp = 0
      nface = 0
      do idir = mdim, 1, -1 ! Start counting from z to simplify implementation for semiperiodic systems (pln)
        do pm = 1, -1, -2 
          nface = nface + 1
          this%face_idx_range(nface, 1) = isp + 1
        
          do ip_global = 1, mesh%np_global
            if(abs(which_surface(ip_global)) == idir .and. sign(1, which_surface(ip_global)) == pm) then
              isp = isp + 1
              ! coordinate of surface point
              xx(1:mdim) = mesh_x_global(mesh, ip_global)
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
        
        end do      
      end do
      
      
      !Get dimensions and spacing on each (pair of) face 
      do ifc = 1, nint((nfaces + 0.5) / 2)
        isp_start = this%face_idx_range(ifc, 1)
        isp_end   = this%face_idx_range(ifc, 2)
      
        !retrieve face normal direction
        n_dir = 0 
        do idir = 1, mdim
          if(abs(this%srfcnrml(idir, isp_start)) >= M_EPSILON) n_dir = idir
        end do 
      
        !retrieve the spacing on the face
        dRS(:) = M_ZERO
        idim = 1
        do idir = 1, mdim 
          if (idir == n_dir ) cycle
          dRS(idim)= mesh%spacing(idir)
          idim = idim + 1
        end do

      
        !retrieve the dimensions of the face
        RSmin = M_ZERO
        RSmax = M_ZERO
        do isp = isp_start, isp_end
          !measure in reduced coordinates 
          xx(1:mdim) = matmul(this%rcoords(1:mdim, isp), mesh%sb%latt%klattice_primitive(1:mdim, 1:mdim))
          idim = 1
          do idir = 1, mdim 
            if (idir == n_dir ) cycle
            RS(idim)=xx(idir)
            if (RS(idim) < RSmin(idim)) RSmin(idim) = RS(idim)
            if (RS(idim) > RSmax(idim)) RSmax(idim) = RS(idim)
            idim = idim + 1
          end do        
        end do
        
   
        do idir = 1, mdim - 1
          this%LLr(n_dir,idir) = RSmax(idir) - RSmin(idir) + dRS(idir)
          if(dRS(idir) > M_ZERO) this%NN(n_dir,idir) = int(this%LLr(n_dir,idir) / dRS(idir))
        end do        

      end do
      
      
    end if
    
    if(this%anisotrpy_correction) then
      !Compensate the fact that the angular distribution of points is not uniform
      !and have peaks in correspondence to the edges and corners of the parallelepiped
      do isp = 1, this%nsrfcpnts
        weight = sum(this%rcoords(1:mdim, isp) * this%srfcnrml(1:mdim, isp))
        weight = weight/sum(this%rcoords(1:mdim, isp)**2) / sum(this%srfcnrml(1:mdim, isp)**2)
        this%srfcnrml(1:mdim, isp) = this%srfcnrml(1:mdim, isp) * abs(weight)                
      end do        
      
    end if


    SAFE_DEALLOCATE_A(which_surface)

    POP_SUB(pes_flux_getcube)
  end subroutine pes_flux_getcube

  ! ---------------------------------------------------------
  subroutine pes_flux_getsphere(this, mesh, nstepsthetar, nstepsphir)
    type(pes_flux_t), intent(inout) :: this
    type(mesh_t),     intent(in)    :: mesh
    integer,          intent(inout) :: nstepsthetar, nstepsphir
    
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
        weight = (M_ONE - cos(dthetar / M_TWO)) * M_TWO * M_PI
      else
        weight = abs(cos(thetar - dthetar / M_TWO) - cos(thetar + dthetar / M_TWO)) &
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
