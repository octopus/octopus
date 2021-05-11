!! Copyright (C) 2006-2011 U. De Giovannini, M. Marques 
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

module pes_mask_oct_m
  use batch_oct_m
  use box_sphere_oct_m
  use box_parallelepiped_oct_m
  use boundary_op_oct_m
  use comm_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use density_oct_m
  use fft_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_binary_oct_m
  use io_function_oct_m
  use io_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lasers_oct_m
  use loct_oct_m
  use math_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
#if defined(HAVE_NETCDF)
  use netcdf
#endif  
  use output_oct_m
  use parser_oct_m
  use profiling_oct_m
  use qshep_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use space_oct_m
  use sort_oct_m
  use states_elec_dim_oct_m
  use states_elec_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use vtk_oct_m
  use wfs_elec_oct_m
  
  implicit none

  private
  
  public ::                             &
    pes_mask_t,                         &
    pes_mask_init,                      &
    pes_mask_end,                       &
    pes_mask_calc,                      &
    pes_mask_mesh_to_cube,              &
    pes_mask_cube_to_mesh,              &
    pes_mask_generate_mask_function,    &
    pes_mask_x_to_k,                    &
    pes_mask_k_to_x,                    &
    pes_mask_read_info,                 &
    pes_mask_output_full_mapm,            &
    pes_mask_output_ar_spherical_cut_m,   &
    pes_mask_output_ar_plane_m,           &
    pes_mask_output_ar_polar_m,           &
    pes_mask_output_full_mapm_cut,        &
    pes_mask_output_power_totalm,         &
    pes_mask_load,                      &
    pes_mask_dump,                      &
    pes_mask_output,                    &
    pes_mask_pmesh,                     &
    pes_mask_map_from_states
    
  
  type pes_mask_t
    private
    CMPLX, allocatable, public :: k(:,:,:,:,:,:) !< The states in momentum space
                                                 !< mask%k(ll(1),ll(2),ll(3),st%d%dim, st%nst, st%d%nik)
                                               
    ! mesh- and cube-related stuff      
    integer          :: np                     !< number of mesh points associated with the mesh
    !< (either mesh%np or mesh%np_global)
    integer          :: ll(3)                  !< the size of the parallelepiped mesh
    integer          :: fs_n_global(1:3)       !< the dimensions of the cube in fourier space
    integer          :: fs_n(1:3)              !< the dimensions of the local portion of the cube in fourier space
    integer          :: fs_istart(1:3)         !< where does the local portion of the cube start in fourier space
    
    FLOAT            :: spacing(3)       !< the spacing
    
    type(mesh_t), pointer, public  :: mesh => NULL()   !< a pointer to the mesh
    type(cube_t)                   :: cube             !< the cubic mesh
    
    FLOAT, allocatable, public :: vec_pot(:,:) !< external time-dependent potential i.e. the lasers
    
    FLOAT, allocatable, public :: Mk(:,:,:)    !< the momentum space filter
    type(cube_function_t)      :: cM           !< the mask cube function
    FLOAT, allocatable, public :: mask_R(:)    !< the mask inner (component 1) and outer (component 2) radius
    integer                    :: shape        !< which mask function?
    FLOAT, allocatable         :: ufn(:)       !< user-defined mask function
    logical                    :: user_def
    
    FLOAT, allocatable, public :: Lk(:,:)      !< associate a k value to a cube index Lk(i,{1,2,3})={kx,ky,kz}(i)
    
    FLOAT            :: enlarge(3)             !< Fourier space enlargement
    FLOAT            :: enlarge_2p(3)          !< Two-point space enlargement
    
    FLOAT :: start_time              !< the time we switch on the photoelectron detector   
    FLOAT :: energyMax 
    FLOAT :: energyStep 
    
    integer :: sw_evolve             !< choose the time propagator for the continuum wfs
    logical :: back_action           !< whether to enable back action from B to A
    logical :: add_psia              !< add the contribution of Psi_A in the buffer region to the output
    logical :: filter_k              !< whether to filter the wavefunctions in momentum space
    
    integer :: mode                  !< calculation mode
    integer :: pw_map_how            !< how to perform projection on plane waves
    
    type(fft_t)    :: fft            !< FFT plan
        
    type(mesh_cube_parallel_map_t) :: mesh_cube_map  !< The parallel map
    
    
  end type pes_mask_t

  
  integer, public, parameter ::        &
    PES_MASK_MODE_MASK         =   1,  &  
    PES_MASK_MODE_BACKACTION   =   2,  &  
    PES_MASK_MODE_PASSIVE      =   3
  
  integer, parameter ::       &
    PW_MAP_FFT    =  3,  &    !< FFT - normally from fftw3
    PW_MAP_NFFT        =  5,  &    !< non-equispaced fft (NFFT)
    PW_MAP_PFFT        =  6,  &    !< use PFFT
    PW_MAP_PNFFT       =  7        !< use PNFFT
  
  integer, parameter ::       &
    M_SIN2            =  1,   &  
    M_STEP            =  2,   & 
    M_ERF             =  3    
  
  integer, parameter ::       &
    IN                =  1,   &  
    OUT               =  2

  integer, public, parameter ::   &
    INTEGRATE_NONE    = -1,       &
    INTEGRATE_PHI     =  1,       &  
    INTEGRATE_THETA   =  2,       &  
    INTEGRATE_R       =  3,       &  
    INTEGRATE_KX      =  4,       &  
    INTEGRATE_KY      =  5,       &  
    INTEGRATE_KZ      =  6

  
contains 


  ! ---------------------------------------------------------
  subroutine pes_mask_init(mask, namespace, space, mesh, sb, st, hm, max_iter,dt)
    type(pes_mask_t),         intent(out) :: mask
    type(namespace_t),        intent(in)  :: namespace
    type(space_t),            intent(in)  :: space
    type(mesh_t), target,     intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    type(states_elec_t),      intent(in)  :: st
    type(hamiltonian_elec_t), intent(in)  :: hm
    integer,                  intent(in)  :: max_iter
    FLOAT,                    intent(in)  :: dt
    
    type(block_t) :: blk
    
    integer :: il, it, ll(3)
    FLOAT :: field(3)
    FLOAT :: DeltaE, MaxE, pCutOff, tmp
    integer :: defaultMask,k1,k2,st1,st2
    integer :: cols_pesmask_block, idim, ip

    FLOAT :: xx(space%dim), r
    FLOAT :: ufn_re, ufn_im
    character(len=1024) :: user_def_expr
   
    PUSH_SUB(pes_mask_init)
        
    mask%mesh => mesh
    
    if (space%is_periodic()) &
      call messages_experimental("PES_mask with periodic dimensions")
    
    
    write(message(1),'(a,i1,a)') 'Info: Calculating PES using mask technique.'
    call messages_info(1)
    
    
    select type (box => sb%box)
    type is (box_sphere_t)
    class default
      if (.not. space%is_periodic()) then
        message(1) = 'PhotoElectronSpectrum = pes_mask usually requires BoxShape = sphere.'
        message(2) = 'Unless you know what you are doing modify this parameter and rerun.'
        call messages_warning(2, namespace=namespace)
      end if
    end select

    if(hm%bc%abtype /= NOT_ABSORBING) then
      message(1) = 'PhotoElectronSpectrum = pes_mask already contains absorbing boundaries.'
      message(2) = 'Set AbsorbingBoundaries = no and rerun.'
      call messages_fatal(2, namespace=namespace)
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Calculation mode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !%Variable PESMaskMode
    !%Type integer
    !%Default mask_mode
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% PES calculation mode.
    !%Option mask_mode 1
    !% Mask method. 
    !%Option fullmask_mode 2
    !% Full mask method. This includes a back action of the momentum-space states on the 
    !% interaction region. This enables electrons to come back from the continuum. 
    !%Option passive_mode 3 
    !% Passive analysis of the wf. Simply analyze the plane-wave components of the 
    !% wavefunctions on the region <i>r</i> > <i>R1</i>. This mode employs a step masking function by default.
    !%End
    call parse_variable(namespace, 'PESMaskMode', PES_MASK_MODE_MASK, mask%mode)
    if(.not.varinfo_valid_option('PESMaskMode', mask%mode)) call messages_input_error(namespace, 'PESMaskMode')
    call messages_print_var_option(stdout, "PESMaskMode", mask%mode)
    
    select case(mask%mode)
    case(PES_MASK_MODE_PASSIVE)
      defaultMask = M_STEP   
      mask%back_action = .false.   
      
    case(PES_MASK_MODE_BACKACTION)
      defaultMask = M_SIN2   
      mask%back_action = .true.
      mask%mode = PES_MASK_MODE_MASK
    case default
      defaultMask = M_SIN2
      mask%back_action = .false.   
    end select
    
    
    !%Variable PESMaskStartTime 
    !%Type float
    !%Default -1.0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The time photoelectrons start to be recorded. In pump-probe simulations, this allows
    !% getting rid of an unwanted ionization signal coming from the pump.
    !% NOTE: This will enforce the mask boundary conditions for all times. 
    !%End
    call parse_variable(namespace, 'PESMaskStartTime', -M_ONE, mask%start_time, unit = units_inp%time)

    !%Variable PESMaskPlaneWaveProjection
    !%Type integer
    !%Default fft_map
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% With the mask method, wavefunctions in the continuum are treated as plane waves.
    !% This variable sets how to calculate the plane-wave projection in the buffer 
    !% region. We perform discrete Fourier transforms (DFT) in order to approximate  
    !% a continuous Fourier transform. The major drawback of this approach is the built-in
    !% periodic boundary condition of DFT. Choosing an appropriate plane-wave projection 
    !% for a given simulation in addition to <tt>PESMaskEnlargeFactor</tt> and 
    !% <tt>PESMask2PEnlargeFactor</tt>will help to converge the results.   
    !%
    !% NOTE: depending on the value of <tt>PESMaskMode</tt> <tt>PESMaskPlaneWaveProjection</tt>,
    !% may affect not only performance but also the time evolution of the density. 
    !%Option fft_out 2 
    !% FFT filtered in order to keep only outgoing waves. 1D only. 
    !%Option fft_map 3 
    !% FFT transform.
    !%Option nfft_map 5
    !% Non-equispaced FFT map. 
    !%Option pfft_map 6
    !% Use PFFT library. 
    !%Option pnfft_map 7
    !% Use PNFFT library. 
    !%End
    call parse_variable(namespace, 'PESMaskPlaneWaveProjection', PW_MAP_FFT, mask%pw_map_how)
    
    if(.not.varinfo_valid_option('PESMaskPlaneWaveProjection', mask%pw_map_how)) then
      call messages_input_error(namespace, 'PESMaskPlaneWaveProjection')
    end if
    
    call messages_print_var_option(stdout, "PESMaskPlaneWaveProjection", mask%pw_map_how)

    if (mask%pw_map_how ==  PW_MAP_PFFT .and. (.not. mask%mesh%parallel_in_domains)) then
      message(1)= "Trying to use PESMaskPlaneWaveProjection = pfft_map with no domain parallelization."
      message(2)= "Projection method changed to more efficient fft_map."
      call messages_warning(2, namespace=namespace)
      mask%pw_map_how = PW_MAP_FFT
    end if

    if (mask%pw_map_how ==  PW_MAP_PNFFT .and. (.not. mask%mesh%parallel_in_domains)) then
      message(1)= "Trying to use PESMaskPlaneWaveProjection = pnfft_map with no domain parallelization."
      message(2)= "Projection method changed to more efficient nfft_map."
      call messages_warning(2, namespace=namespace)
      mask%pw_map_how = PW_MAP_NFFT
    end if
    
#if !defined(HAVE_NFFT) 
    if (mask%pw_map_how ==  PW_MAP_NFFT) then
      message(1) = "PESMaskPlaneWaveProjection = nfft_map requires NFFT but that library was not linked."
      call messages_fatal(1, namespace=namespace)
    end if
#endif
    
#if !defined(HAVE_PFFT) 
    if (mask%pw_map_how ==  PW_MAP_PFFT) then
      message(1) = "PESMaskPlaneWaveProjection = pfft_map requires PFFT but that library was not linked."
      call messages_fatal(1, namespace=namespace)
    end if
#endif

#if !defined(HAVE_PNFFT) 
    if (mask%pw_map_how ==  PW_MAP_PNFFT) then
      message(1) = "PESMaskPlaneWaveProjection = pnfft_map requires PNFFT but that library was not linked."
      call messages_fatal(1, namespace=namespace)
    end if
#endif
    
    
    !%Variable PESMaskEnlargeFactor
    !%Type float
    !%Default 1
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Mask box enlargement level. Enlarges the mask bounding box by a <tt>PESMaskEnlargeFactor</tt>.
    !% This helps to avoid wavefunction wrapping at the boundaries.
    !%End

    mask%enlarge = M_ONE
    call parse_variable(namespace, 'PESMaskEnlargeFactor', M_ONE, mask%enlarge(1))
    
    if ( mask%enlarge(1) /= M_ONE ) then

      mask%enlarge(space%periodic_dim + 1:space%dim) = mask%enlarge(1)
      mask%enlarge(1:space%periodic_dim) = M_ONE
      
      if (space%is_periodic()) then
        call messages_print_var_value(stdout, "PESMaskEnlargeFactor", mask%enlarge(1:space%dim))
      else
        call messages_print_var_value(stdout, "PESMaskEnlargeFactor", mask%enlarge(1))
      end if
      
    end if
    if( mask%enlarge(1) < M_ONE ) then
      message(1) = "PESMaskEnlargeFactor must be bigger than one."
      call messages_fatal(1, namespace=namespace)
    end if
 
    call messages_obsolete_variable(namespace, 'PESMaskEnlargeLev', 'PESMaskEnlargeFactor')
    
    !%Variable PESMask2PEnlargeFactor
    !%Type float
    !%Default 1.0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Mask two points enlargement factor. Enlarges the mask box by adding two
    !% points at the edges of the box in each direction (x,y,z) at a distance  
    !% L=Lb*<tt>PESMask2PEnlargeFactor</tt> where <i>Lb</i> is the box size.
    !% This allows to run simulations with an additional void space at a price of
    !% adding few points. The Fourier space associated with the new box is restricted 
    !% by the same factor.
    !% 
    !% Note: needs <tt> PESMaskPlaneWaveProjection = nfft_map or pnfft_map </tt>.
    !%End
    
    mask%enlarge_2p = M_ONE
    call parse_variable(namespace, 'PESMask2PEnlargeFactor', M_ONE, mask%enlarge_2p(1))

    
    if ( mask%enlarge_2p(1) /= M_ONE ) then

      mask%enlarge_2p(space%periodic_dim + 1:space%dim) = mask%enlarge_2p(1)
      mask%enlarge_2p(1:space%periodic_dim) = M_ONE

      if (space%is_periodic()) then
        call messages_print_var_value(stdout, "PESMask2PEnlargeFactor", mask%enlarge_2p(1:space%dim))
      else
        call messages_print_var_value(stdout, "PESMask2PEnlargeFactor", mask%enlarge_2p(1))
      end if
                          
      if (mask%pw_map_how /=  PW_MAP_NFFT .and. mask%pw_map_how /=  PW_MAP_PNFFT) then
        message(1) = "PESMask2PEnlargeFactor requires PESMaskPlaneWaveProjection = nfft_map"
        message(2) = "or pnfft_map in order to run properly." 
        call messages_warning(2, namespace=namespace)
      end if        
    end if
    
    if( mask%enlarge_2p(1) < M_ONE ) then
      message(1) = "PESMask2PEnlargeFactor must be bigger than one."
      call messages_fatal(1, namespace=namespace)
    end if
    
    call messages_obsolete_variable(namespace, 'PESMaskNFFTEnlargeLev', 'PESMask2PEnlargeFactor')
    
    
    mask%ll = 1
    mask%spacing = -M_ONE
    
    mask%spacing(1:3) = mesh%spacing(1:3)
    mask%ll(1:3) = mesh%idx%ll(1:3)    
    
    !Enlarge the cube region
    mask%ll(1:space%dim) = int(mask%ll(1:space%dim) * mask%enlarge(1:space%dim))
    
    select case(mask%pw_map_how)
      
    case(PW_MAP_PFFT)
      ASSERT(mask%mesh%parallel_in_domains)
      call cube_init(mask%cube, mask%ll, namespace, space, &
        fft_type = FFT_COMPLEX, fft_library = FFTLIB_PFFT, nn_out = ll, &
        mpi_grp = mask%mesh%mpi_grp, need_partition=.true., spacing = mesh%spacing)
      !        print *,mpi_world%rank, "mask%mesh%mpi_grp%comm", mask%mesh%mpi_grp%comm, mask%mesh%mpi_grp%size
      !         print *,mpi_world%rank, "mask%cube%mpi_grp%comm", mask%cube%mpi_grp%comm, mask%cube%mpi_grp%size
      

!       mask%ll(1) = mask%cube%fs_n(3)
!       mask%ll(2) = mask%cube%fs_n(1)
!       mask%ll(3) = mask%cube%fs_n(2)
!     Note: even if tempting, setting       
!     mask%ll(1:3) = mask%cube%fs_n(1:3) 
!     results in the wrong index mapping! (1->3, 2->1, 3->2) 
!     Well.. very much against intuition it turns out to be 
      mask%ll(1:3) = mask%cube%rs_n(1:3) 

      mask%fft = mask%cube%fft
      mask%np = mesh%np_part ! the mask is local
      if ( mask%mesh%parallel_in_domains .and. mask%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_init(mask%mesh_cube_map, mask%mesh, mask%cube)
      end if
      
    case(PW_MAP_FFT)
      call cube_init(mask%cube, mask%ll, namespace, space, fft_type = FFT_COMPLEX, fft_library = FFTLIB_FFTW, &
        nn_out = ll, spacing = mesh%spacing)
      mask%ll = ll 
      mask%fft = mask%cube%fft
      mask%np = mesh%np_part_global 

            
    case(PW_MAP_NFFT)
      
      !NFFT initialization
      ! we just add 2 points for the enlarged region
      if (mask%enlarge_2p(1) /= 1) mask%ll(1:space%dim) = mask%ll(1:space%dim) + 2 

      call cube_init(mask%cube, mask%ll, namespace, space,  fft_type = FFT_COMPLEX, fft_library = FFTLIB_NFFT, &
        nn_out = ll, spacing = mesh%spacing, tp_enlarge = mask%enlarge_2p)

      mask%ll = ll 
      mask%fft = mask%cube%fft
      mask%np = mesh%np_part_global       
 
      
    case(PW_MAP_PNFFT)  
    
      if (mask%enlarge_2p(1) /= 1) mask%ll(1:space%dim) = mask%ll(1:space%dim) + 2 

      call cube_init(mask%cube, mask%ll, namespace, space, fft_type = FFT_COMPLEX, fft_library = FFTLIB_PNFFT, &
        nn_out = ll, spacing = mesh%spacing, tp_enlarge = mask%enlarge_2p, &
        mpi_grp = mask%mesh%mpi_grp, need_partition=.true.)

      mask%ll(1:3) = mask%cube%fs_n(1:3) 

      mask%fft = mask%cube%fft
      mask%np = mesh%np_part ! the mask is local
      if ( mask%mesh%parallel_in_domains .and. mask%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_init(mask%mesh_cube_map, mask%mesh, mask%cube)
      end if

    case default 
      !Program should die before coming here
      write(message(1),'(a)') "PESMaskPlaneWaveProjection unrecognized option." 
      call messages_fatal(1, namespace=namespace)
      
    end select
    
    !Indices  
    
    if (mask%pw_map_how == PW_MAP_PFFT) then
      mask%fs_istart = mask%cube%rs_istart 
      mask%fs_n = mask%cube%rs_n 
      mask%fs_n_global = mask%cube%rs_n_global 
    else
      mask%fs_istart = mask%cube%fs_istart 
      mask%fs_n = mask%cube%fs_n 
      mask%fs_n_global = mask%cube%fs_n_global 
    end if 

    if(debug%info) then
      print *,mpi_world%rank, "mask%ll                  ", mask%ll(:)
      print *,mpi_world%rank, "mask%cube%fs_n_global(:) ", mask%cube%fs_n_global(:)      
      print *,mpi_world%rank, "mask%cube%fs_n(:)        ", mask%cube%fs_n(:)
      print *,mpi_world%rank, "mask%cube%fs_istart(:)   ", mask%cube%fs_istart(:)
      print *,mpi_world%rank, "mask%cube%rs_n_global(:) ", mask%cube%rs_n_global(:)      
      print *,mpi_world%rank, "mask%cube%rs_n(:)        ", mask%cube%rs_n(:)
      print *,mpi_world%rank, "mask%cube%rs_istart(:)   ", mask%cube%rs_istart(:)
    end if 

!     print *, mpi_world%rank, " mask%ll", mask%ll(1:3), "states -",st%st_start,st%st_end
    !Allocations
    call zcube_function_alloc_RS(mask%cube, mask%cM, force_alloc = .true.)
    
    SAFE_ALLOCATE(mask%Lk(1:maxval(mask%fs_n_global(:)),1:3))
    mask%Lk(:,:) = M_ZERO
    
    st1 = st%st_start
    st2 = st%st_end
    k1 = st%d%kpt%start
    k2 = st%d%kpt%end   
    SAFE_ALLOCATE(mask%k(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3),1:st%d%dim,st1:st2,k1:k2))
    mask%k = M_z0
    
    
    ! generate the map between mesh and cube
    call  pes_mask_generate_Lk(mask, sb) ! generate the physical momentum vector
    
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Mask Function options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SAFE_ALLOCATE(mask%mask_R(1:2))
    
    !%Variable PESMaskShape
    !%Type integer
    !%Default m_sin2
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The mask function shape.
    !%Option m_sin2 1
    !% sin2 mask.
    !%Option m_step 2 
    !%step function.  
    !%Option m_erf 3 
    !%Error function. Not Implemented.
    !%End
    call parse_variable(namespace, 'PESMaskShape', defaultMask, mask%shape)
    if(.not.varinfo_valid_option('PESMaskShape', mask%shape)) call messages_input_error(namespace, 'PESMaskShape')
    call messages_print_var_option(stdout, "PESMaskShape", mask%shape)
    
    !%Variable PESMaskSize
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Set the size of the mask function. 
    !% Here you can set the inner (R1) and outer (R2) radius by setting
    !% the block as follows:
    !%
    !% <tt>%PESMaskSize
    !% <br>&nbsp;&nbsp; R1 | R2 | "user-defined"
    !% <br>%</tt>
    !%
    !% The optional 3rd column is a user-defined expression for the mask 
    !% function. For example, <i>r</i> creates a spherical mask (which is the 
    !% default for <tt>BoxShape = sphere</tt>). Note, values R2 larger than 
    !% the box size may lead in this case to unexpected reflection 
    !% behaviours.
    !%End
    cols_pesmask_block = 0
    if (parse_block(namespace, 'PESMaskSize', blk) == 0) then
      cols_pesmask_block = parse_block_cols(blk, 0)
    end if

    mask%user_def = .false.

    select case(cols_pesmask_block)
    case(0)
      select type (box => sb%box)
      type is (box_sphere_t)
        mask%mask_R(1) = box%radius/M_TWO
        mask%mask_R(2) = box%radius
        message(1) = "Input: PESMaskSize R(1) and R(2) not specified. Using default values for spherical mask."
      type is (box_parallelepiped_t)
        mask%mask_R(1) = box%half_length(1)/M_TWO
        mask%mask_R(2) = box%half_length(1)
        message(1) = "Input: PESMaskSize R(1) and R(2) not specified. Using default values for cubic mask."
      end select
      call messages_info(1)
    case(1)
      call parse_block_float(blk, 0, 0, mask%mask_R(1), units_inp%length)
      select type (box => sb%box)
      type is (box_sphere_t)
        mask%mask_R(2) = box%radius
        message(1) = "Input: PESMaskSize R(2) not specified. Using default value for spherical mask."
      type is (box_parallelepiped_t)
        mask%mask_R(2) = box%half_length(1)
        message(1) = "Input: PESMaskSize R(2) not specified. Using default value for cubic mask."
      end select
      call messages_info(1)
    case(2)
      call parse_block_float(blk, 0, 0, mask%mask_R(1), units_inp%length)
      call parse_block_float(blk, 0, 1, mask%mask_R(2), units_inp%length)

      select type (box => sb%box)
      type is (box_sphere_t)
        if (mask%mask_R(2) > box%radius)  mask%mask_R(2) = box%radius
        message(1) = "Info: using spherical mask."
      type is (box_parallelepiped_t)
        if (mask%mask_R(2) > box%half_length(1)) mask%mask_R(2) = box%half_length(1)
        message(1) = "Info: using cubic mask."
      end select
      call messages_info(1)

    case(3)
      mask%user_def = .true.
      SAFE_ALLOCATE(mask%ufn(1:mask%np))
      mask%ufn = M_ZERO
      call parse_block_float(blk, 0, 0, mask%mask_R(1), units_inp%length)
      call parse_block_float(blk, 0, 1, mask%mask_R(2), units_inp%length)
      call parse_block_string(blk, 0, 2, user_def_expr)
      do ip = 1, mask%np
        xx = M_ZERO
        xx(1:space%dim) = mesh%x(ip, 1:space%dim)
        r = units_from_atomic(units_inp%length, sqrt(sum(xx(1:space%dim)**2)))
        do idim = 1, space%dim
          xx(idim) = units_from_atomic(units_inp%length, xx(idim))
        end do
        call parse_expression(ufn_re, ufn_im, space%dim, xx, r, M_ZERO, user_def_expr)
        mask%ufn(ip) = ufn_re
      end do
      message(1) = "Input: using user-defined mask function from expression:"
      write(message(2),'(a,a)') '   R = ', trim(user_def_expr) 
      call messages_info(2)
    end select

    call parse_block_end(blk)

    write(message(1),'(a,es10.3,3a)') & 
      "  R1 = ", units_from_atomic(units_inp%length, mask%mask_R(1) ),&
      ' [', trim(units_abbrev(units_inp%length)), ']'
    write(message(2),'(a,es10.3,3a)') & 
      "  R2 = ", units_from_atomic(units_inp%length, mask%mask_R(2) ),&
      ' [', trim(units_abbrev(units_inp%length)), ']'
    call messages_info(2)
    
    
    call pes_mask_generate_mask(mask, namespace, mesh)
    
    !%Variable PESMaskFilterCutOff 
    !%Type float
    !%Default -1
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% In calculation with <tt>PESMaskMode = fullmask_mode</tt> and NFFT, spurious frequencies 
    !% may lead to numerical instability of the algorithm. This option gives the possibility 
    !% to filter out the unwanted components by setting an energy cut-off. 
    !% If <tt>PESMaskFilterCutOff = -1</tt> no filter is applied.
    !%End
    call parse_variable(namespace, 'PESMaskFilterCutOff', -M_ONE, pCutOff, unit = units_inp%energy)
    
    mask%filter_k = .false.
    
    if(pCutOff > M_ZERO) then       
      call messages_print_var_value(stdout, "PESMaskFilterCutOff", pCutOff, unit = units_out%energy)
      mask%filter_k = .true.
      
      SAFE_ALLOCATE(mask%Mk(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3)))
      
      call pes_mask_generate_filter(mask,pCutOff)
    end if
    
!!  Output

    !%Variable PESMaskIncludePsiA
    !%Type logical
    !%Default false
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Add the contribution of <math>\Psi_A</math> in the mask region to the photo-electron spectrum.
    !% Literally adds the Fourier components of:
    !% <math>\Theta(r-R1) \Psi_A(r)</math>
    !% with <math>\Theta</math> being the Heaviside step function. 
    !% With this option PES will contain all the contributions starting from the inner 
    !% radius <math>R1</math>. Use this option to improve convergence with respect to the box size 
    !% and total simulation time. 
    !% Note: Carefully choose <math>R1</math> in order to avoid contributions from returning electrons. 
    !%End
    call parse_variable(namespace, 'PESMaskIncludePsiA', .false., mask%add_psia)
    if(mask%add_psia) then
      message(1)= "Input: Include contribution from Psi_A."
      call messages_info(1)
    end if
    
    
    !%Variable PESMaskSpectEnergyMax 
    !%Type float
    !%Default maxval(mask%Lk)<math>^2/2</math>
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The maximum energy for the PES spectrum.
    !%End
    MaxE = M_EPSILON
    do idim = 1, space%dim
      tmp = maxval(mask%Lk(1:mask%ll(idim),1:space%dim))**M_TWO/M_TWO
      if (tmp > MaxE) MaxE = tmp
    end do
    call parse_variable(namespace, 'PESMaskSpectEnergyMax', MaxE, mask%energyMax, unit = units_inp%energy)
    call messages_print_var_value(stdout, "PESMaskSpectEnergyMax", mask%energyMax, unit = units_out%energy)

    !%Variable PESMaskSpectEnergyStep 
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The PES spectrum energy step.
    !%End
    DeltaE = minval(mask%Lk(2,1:space%dim) - mask%Lk(1,1:space%dim))**M_TWO/M_TWO
    call parse_variable(namespace, 'PESMaskSpectEnergyStep', DeltaE, mask%energyStep, unit = units_inp%energy)
    call messages_print_var_value(stdout, "PESMaskSpectEnergyStep", mask%energyStep, unit = units_out%energy)
    

!!  Set external fields 

    SAFE_ALLOCATE(mask%vec_pot(0:max_iter,1:3))
    mask%vec_pot=M_ZERO

    do il = 1, hm%ext_lasers%no_lasers
      select case(laser_kind(hm%ext_lasers%lasers(il)))
      case(E_FIELD_VECTOR_POTENTIAL)
        do it = 1, max_iter
          field=M_ZERO
          call laser_field(hm%ext_lasers%lasers(il), field, it*dt)
          ! We must sum with a -1 sign to account for the 
          ! electron charge.
          mask%vec_pot(it,:)= mask%vec_pot(it,:) - field(:)           
        end do
        
      case default 
        write(message(1),'(a)') 'PESMask should work only with TDExternalFields = vector_potential.'
        write(message(2),'(a)') 'Unless PESMaskMode = passive_mode the results are likely to be wrong. '
        call messages_warning(2, namespace=namespace)
        
      end select
    end do

    

    ! Compensate the sign for the forward <-> backward FT inversion 
    ! used with NFFT. See pes_mask_X_to_K comment for more info
    if(mask%pw_map_how == PW_MAP_PNFFT .or. &
       mask%pw_map_how == PW_MAP_NFFT) mask%Lk = - mask%Lk
        
        
    POP_SUB(pes_mask_init)
  end subroutine pes_mask_init


  ! ---------------------------------------------------------
  subroutine pes_mask_end(mask)
    type(pes_mask_t), intent(inout) :: mask
    
    PUSH_SUB(pes_mask_end)
    
    SAFE_DEALLOCATE_A(mask%k)
    
    SAFE_DEALLOCATE_A(mask%vec_pot)
    SAFE_DEALLOCATE_A(mask%mask_R)
    SAFE_DEALLOCATE_A(mask%Lk)
    
    
    if (mask%filter_k) then
      SAFE_DEALLOCATE_A(mask%Mk)
    end if
    
    if (mask%mesh%parallel_in_domains .and. mask%cube%parallel_in_domains) then
      call mesh_cube_parallel_map_end(mask%mesh_cube_map)
    end if
    
    call zcube_function_free_RS(mask%cube, mask%cM)
    call cube_end(mask%cube)   

    if (mask%user_def) then 
      SAFE_DEALLOCATE_A(mask%ufn)
    end if
    
    POP_SUB(pes_mask_end)
  end subroutine pes_mask_end
  

  ! --------------------------------------------------------
  subroutine pes_mask_generate_Lk(mask, sb)
    type(pes_mask_t),  intent(inout) :: mask
    type(simul_box_t),   intent(in)  :: sb

    integer :: ii, dim

    PUSH_SUB(pes_mask_generate_Lk)

    dim = mask%mesh%sb%dim

    do ii = 1, maxval(mask%ll(:))
      mask%Lk(ii,1:dim)= matmul(sb%latt%klattice_primitive(1:dim,1:dim), mask%cube%Lfs(ii,1:dim))
    end do

    POP_SUB(pes_mask_generate_Lk)
  end subroutine pes_mask_generate_Lk

  ! ======================================
  !>  Generate the momentum-space filter 
  ! ======================================
  subroutine pes_mask_generate_filter(mask,cutOff)
    type(pes_mask_t), intent(inout) :: mask
    FLOAT,            intent(in)    ::cutOff
    
    
    integer :: kx,ky,kz, power
    FLOAT   :: KK(3), EE, Emax

    PUSH_SUB(pes_mask_generate_filter)
    
    mask%Mk = M_ZERO

    Emax = maxval(mask%Lk(:,:))**2 / M_TWO

    power = 8 

    do kx = 1, mask%ll(1)
      KK(1) = mask%Lk(kx,1) 
      do ky = 1, mask%ll(2)
        KK(2) = mask%Lk(ky,2)
        do kz = 1, mask%ll(3)
          KK(3) = mask%Lk(kz,3)
          
          EE = sum(KK(1:mask%mesh%sb%dim)**2) / M_TWO

          if ( EE > cutOff .and. EE < Emax ) then
            mask%Mk(kx,ky,kz) = M_ONE * sin((Emax-EE) * M_PI / (M_TWO * (cutOff) ))**power
          else 
            if( EE <= cutOff )  mask%Mk(kx,ky,kz) = M_ONE 
          end if
          
        end do
      end do
    end do


    POP_SUB(pes_mask_generate_filter)
  end subroutine pes_mask_generate_filter


  ! --------------------------------------------------------
  !>  Generate the mask function on the cubic mesh containing 
  !!  the simulation box
  ! ---------------------------------------------------------
  subroutine pes_mask_generate_mask(mask, namespace, mesh)
    type(pes_mask_t),  intent(inout) :: mask
    type(namespace_t), intent(in)    :: namespace
    type(mesh_t),      intent(in)    :: mesh

    PUSH_SUB(pes_mask_generate_mask)
    
    call pes_mask_generate_mask_function(mask, namespace, mesh, mask%shape, mask%mask_R)
    
    POP_SUB(pes_mask_generate_mask)
    
  end subroutine pes_mask_generate_mask

  ! --------------------------------------------------------
  !>  Generate the mask function on the cubic mesh containing 
  !!  the simulation box
  ! ---------------------------------------------------------
  subroutine pes_mask_generate_mask_function(mask, namespace, mesh, shape, R, mask_sq)
    type(pes_mask_t),  intent(inout) :: mask
    type(namespace_t), intent(in)    :: namespace
    type(mesh_t),      intent(in)    :: mesh
    integer,           intent(in)    :: shape
    FLOAT,             intent(in)    :: R(2)
    FLOAT, optional,   intent(out)   :: mask_sq(:,:,:)

    integer :: ip, dir
    FLOAT   :: width
    FLOAT   :: xx(1:mesh%sb%dim), rr, dd, ddv(1:mesh%sb%dim), tmp(1:mesh%sb%dim)
    CMPLX, allocatable :: mask_fn(:)
    logical :: local_

    PUSH_SUB(pes_mask_generate_mask_function)


    ! generate the mask function on the mesh 
    SAFE_ALLOCATE(mask_fn(1:mask%np))

    mask_fn = M_ZERO
    width = R(2) - R(1)
    xx = M_ZERO

    !We want the mask cube function to be divided on the nodes?
    local_= mask%cube%parallel_in_domains

    select case(shape)
    case(M_SIN2)
      do ip = 1, mask%np
        if(local_) then
          call mesh_r(mesh, ip, rr, coords=xx)
        else
          xx = mesh_x_global(mesh, ip) 
          rr = sqrt(dot_product(xx, xx))
        end if

        if(mask%user_def) then
          dd = mask%ufn(ip) - R(1)
          if(dd > M_ZERO) then
            if(mask%ufn(ip) < R(2) ) then
              mask_fn(ip) = M_ONE * sin(dd * M_PI / (M_TWO * (width) ))**2
            else
              mask_fn(ip) = M_ONE
            end if
          end if

        else ! mask%user_def == .false.

          select type (box => mesh%sb%box)
          type is (box_sphere_t)
            dd = rr -  R(1) 
            if(dd > M_ZERO ) then 
              if (dd  <  width) then
                mask_fn(ip) = M_ONE * sin(dd * M_PI / (M_TWO * (width) ))**2
              else 
                mask_fn(ip) = M_ONE 
              end if
            end if

          type is (box_parallelepiped_t)
          
            ! We are filling from the center opposite to the spherical case
            tmp = M_ONE
            mask_fn(ip) = M_ONE
            ddv = abs(xx) -  R(1) 
            do dir=1, mesh%sb%dim 
              if(ddv(dir) > M_ZERO ) then 
                if (ddv(dir)  <  width) then
                  tmp(dir) = M_ONE - sin(ddv(dir) * M_PI / (M_TWO * (width) ))**2
                else 
                  tmp(dir) = M_ZERO
                end if
              end if        
            mask_fn(ip) = mask_fn(ip)*tmp(dir)
            end do
            mask_fn(ip) = M_ONE - mask_fn(ip)
          
          end select
        end if
      end do

    case(M_STEP)
      do ip = 1, mask%np
        if(local_) then
          call mesh_r(mesh, ip, rr, coords=xx)
        else
          xx = mesh_x_global(mesh, ip) 
          rr = sqrt(dot_product(xx, xx))
        end if
        dd = rr - R(1) 
        if(dd > M_ZERO ) then 
          if (dd  <  width) then
            mask_fn(ip) = M_ONE 
          else 
            mask_fn(ip) = M_ZERO
          end if
        end if
      end do
      
    case(M_ERF)
      
    case default
      message(1)="PhotoElectronSpectrum = pes_mask. Unrecognized mask type."
      call messages_fatal(1, namespace=namespace)
    end select



    mask_fn(:) = M_ONE - mask_fn(:)


    call pes_mask_mesh_to_cube(mask, mask_fn, mask%cM, local = local_)

    if(present(mask_sq)) mask_sq = TOFLOAT(mask%cM%zRS)



    SAFE_DEALLOCATE_A(mask_fn)
    
    POP_SUB(pes_mask_generate_mask_function)
    
  end subroutine pes_mask_generate_mask_function


  ! --------------------------------------------------------
  subroutine pes_mask_apply_mask(mask, st)
    type(states_elec_t), intent(inout) :: st
    type(pes_mask_t),    intent(in)    :: mask
    
    integer :: ik, ist, idim
    CMPLX, allocatable :: mmask(:), psi(:)
    
    PUSH_SUB(pes_mask_apply_mask)
    SAFE_ALLOCATE(mmask(1:mask%mesh%np))
    SAFE_ALLOCATE(psi(1:mask%mesh%np))
    
    call pes_mask_cube_to_mesh(mask, mask%cM, mmask)
    
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          call states_elec_get_state(st, mask%mesh, idim, ist, ik, psi)
          psi(1:mask%mesh%np) = psi(1:mask%mesh%np) * mmask(1:mask%mesh%np)
          call states_elec_set_state(st, mask%mesh, idim, ist, ik, psi)
        end do
      end do
    end do
    
    SAFE_DEALLOCATE_A(mmask)
    
    POP_SUB(pes_mask_apply_mask)
  end subroutine pes_mask_apply_mask






  ! ---------------------------------------------------------
  !>  Propagate in time a wavefunction in momentum space with 
  !!  the Volkov Hamiltonian 
  !!\f[
  !!     Hv=(p - A/c)^2/2
  !!\f]
  !! the time evolution becomes
  !!\f[     
  !!     \psi(p,t+dt)=<p|Uv(dt)|\psi(t)> = exp(-i(p-A/c)^2/2 dt) \psi(p,t)
  !!\f]
  !!  with 
  !!
  ! ---------------------------------------------------------
  subroutine pes_mask_Volkov_time_evolution_wf(mask, space, mesh, kpoints, dt, iter, wf, ikpoint)
    type(pes_mask_t), intent(in)    :: mask
    type(space_t),    intent(in)    :: space
    type(mesh_t),     intent(in)    :: mesh
    type(kpoints_t),  intent(in)    :: kpoints
    FLOAT,            intent(in)    :: dt
    integer,          intent(in)    :: iter
    CMPLX,            intent(inout) :: wf(:,:,:)
    integer,          intent(in)    :: ikpoint
    
    integer ::  ix, iy, iz
    FLOAT :: vec
    FLOAT :: KK(1:3), kpoint(1:3)

    PUSH_SUB(pes_mask_Volkov_time_evolution_wf)

    kpoint = M_ZERO
    if (space%is_periodic()) then
      kpoint(1:mesh%sb%dim) = kpoints%get_point(ikpoint)
    end if
  
    do ix = 1, mask%ll(1)
      KK(1) = mask%Lk(ix + mask%fs_istart(1) - 1, 1)
      do iy = 1, mask%ll(2)
        KK(2) = mask%Lk(iy + mask%fs_istart(2) - 1, 2)
        do iz = 1, mask%ll(3)
          KK(3) = mask%Lk(iz + mask%fs_istart(3) - 1, 3)
          ! The k-points have the same sign as the vector potential consistently 
          ! with what is done to generate the phase (hm%phase) in hamiltonian_elec_update()
          vec = sum(( KK(1:mesh%sb%dim) &
                - kpoint(1:mesh%sb%dim) &
                - mask%vec_pot(iter,1:mesh%sb%dim)/P_C)**2) / M_TWO
          wf(ix, iy, iz) = wf(ix, iy, iz) * exp(-M_zI * dt * vec)
        end do
      end do
    end do


!     do ix = 1, mask%ll(1)
!       KK(1) = mask%Lk(ix + mask%fs_istart(1) - 1)
!       do iy = 1, mask%ll(2)
!         KK(2) = mask%Lk(iy + mask%fs_istart(2) - 1)
!         do iz = 1, mask%ll(3)
!           KK(3) = mask%Lk(iz + mask%fs_istart(3) - 1)
!           vec = sum(( KK(1:mesh%sb%dim) - mask%vec_pot(iter,1:mesh%sb%dim)/P_C)**2) / M_TWO
!           wf(ix, iy, iz) = wf(ix, iy, iz) * exp(-M_zI * dt * vec)
!         end do
!       end do
!     end do
    
    POP_SUB(pes_mask_Volkov_time_evolution_wf)
  end subroutine pes_mask_Volkov_time_evolution_wf




  !---------------------------------------------------------
  !> Project the wavefunction on plane waves
  !---------------------------------------------------------
  subroutine pes_mask_X_to_K(mask, wfin, wfout)
    type(pes_mask_t), intent(in)    :: mask
    CMPLX,            intent(inout) :: wfin(:,:,:)
    CMPLX,            intent(out)   :: wfout(:,:,:)
    
    type(profile_t), save :: prof
    type(cube_function_t) :: cf_tmp
    FLOAT                 :: norm 
    
    call profiling_in(prof, "PESMASK_X_to_K")
    
    PUSH_SUB(pes_mask_X_to_K)
    
    wfout=M_z0
    
    select case(mask%pw_map_how)
      
    case(PW_MAP_FFT)
      call zfft_forward(mask%cube%fft, wfin, wfout)
            
    case(PW_MAP_NFFT)
    ! By definition NFFT forward transform generates a function defined on an 
    ! unstructured grid from its Fourier components (defined on a cubic equispaced grid).
    ! This is not what we want since for us it is the real-space grid the one that can be 
    ! unstructured. 
    ! In order to preserve the possibility to have an unstructured rs-grid we use the  
    ! backward transform and flip the sign of all the momenta (Lk = -Lk).
    
      call zfft_backward(mask%cube%fft, wfin, wfout, norm)
      wfout = wfout * norm
      
!       call zfft_forward(mask%cube%fft, wfin, wfout)
      
    case(PW_MAP_PFFT)
      call zcube_function_alloc_RS(mask%cube, cf_tmp)
      call cube_function_alloc_fs(mask%cube, cf_tmp)
      cf_tmp%zRs = wfin
      call zfft_forward(mask%cube%fft, cf_tmp%zRs, cf_tmp%fs)
      wfout = cf_tmp%fs
      call zcube_function_free_RS(mask%cube, cf_tmp)
      call cube_function_free_fs(mask%cube, cf_tmp)

    case(PW_MAP_PNFFT)
      call zfft_backward(mask%cube%fft, wfin, wfout, norm)
      wfout = wfout * norm

!       call zfft_forward(mask%cube%fft, wfin, wfout)

!       call zfft_backward(mask%cube%fft, M_zI*conjg(wfin), wfout)
!       wfout =  M_zI*conjg(wfout)
      
    case default
      
    end select
    
    
    POP_SUB(pes_mask_X_to_K)
    call profiling_out(prof)
    
  end subroutine pes_mask_X_to_K
  
  ! ------------------------------------------------
  subroutine pes_mask_K_to_X(mask, wfin, wfout)
    type(pes_mask_t), intent(in)    :: mask
    CMPLX,            intent(inout) :: wfin(:,:,:)
    CMPLX,            intent(out)   :: wfout(:,:,:)
    
    type(profile_t), save :: prof
    type(cube_function_t) :: cf_tmp
    FLOAT :: norm
    
    call profiling_in(prof, "PESMASK_K_toX")
    
    PUSH_SUB(pes_mask_K_to_X)
    
    wfout=M_z0
    
    select case(mask%pw_map_how)
            
    case(PW_MAP_FFT)
      call zfft_backward(mask%cube%fft, wfin,wfout)
            
    case(PW_MAP_NFFT)
      call zfft_forward(mask%cube%fft, wfin, wfout, norm)
      wfout = wfout / norm
!       call zfft_backward(mask%cube%fft, wfin,wfout)
      
    case(PW_MAP_PFFT)
      
      call zcube_function_alloc_RS(mask%cube, cf_tmp)
      call cube_function_alloc_fs(mask%cube, cf_tmp)
      cf_tmp%fs  = wfin
      call zfft_backward(mask%cube%fft, cf_tmp%fs, cf_tmp%zRs)
      wfout = cf_tmp%zRs
      call zcube_function_free_RS(mask%cube, cf_tmp)
      call cube_function_free_fs(mask%cube, cf_tmp)

    case(PW_MAP_PNFFT)
      call zfft_forward(mask%cube%fft, wfin, wfout, norm)
      wfout = wfout / norm
      
!       call zfft_backward(mask%cube%fft, wfin,wfout)

!       call zfft_forward(mask%cube%fft, M_zI*conjg(wfin),wfout)
!       wfout =  M_zI*conjg(wfout)
      

    case default

    end select


    POP_SUB(pes_mask_K_to_X)
    
    call profiling_out(prof)

  end subroutine pes_mask_K_to_X

  !---------------------------------------------------------
  subroutine pes_mask_mesh_to_cube(mask, mf, cf, local)
    type(pes_mask_t),      intent(in)    :: mask
    CMPLX,                 intent(in)    :: mf(:)
    type(cube_function_t), intent(inout) :: cf
    logical, optional,     intent(in)    :: local
    
    logical :: local_
    
    PUSH_SUB(pes_mask_mesh_to_cube)
    
    local_ = optional_default(local, .true.)
    
    if (mask%cube%parallel_in_domains) then
      call zmesh_to_cube_parallel(mask%mesh, mf, mask%cube, cf, mask%mesh_cube_map)
    else
      if(mask%mesh%parallel_in_domains) then
        call zmesh_to_cube(mask%mesh, mf, mask%cube, cf, local = local_)
      else 
        call zmesh_to_cube(mask%mesh, mf, mask%cube, cf)
      end if
    end if
    
    POP_SUB(pes_mask_mesh_to_cube)
  end subroutine pes_mask_mesh_to_cube


  !---------------------------------------------------------
  subroutine pes_mask_cube_to_mesh(mask, cf, mf)
    type(pes_mask_t),      intent(in) :: mask
    CMPLX,                 intent(out):: mf(:)
    type(cube_function_t), intent(in) :: cf
    
    PUSH_SUB(pes_mask_cube_to_mesh)
    
    if (mask%cube%parallel_in_domains) then
      call zcube_to_mesh_parallel(mask%cube, cf, mask%mesh, mf, mask%mesh_cube_map)
    else
      if(mask%mesh%parallel_in_domains) then
        call zcube_to_mesh(mask%cube, cf, mask%mesh, mf, local = .true.)
      else 
        call zcube_to_mesh(mask%cube, cf, mask%mesh, mf)
      end if
    end if
    
    POP_SUB(pes_mask_cube_to_mesh)
  end subroutine pes_mask_cube_to_mesh


  !---------------------------------------------------------
  !
  !            Performs all the dirty work 
  !
  !---------------------------------------------------------
  subroutine pes_mask_calc(mask, namespace, space, mesh, st, kpoints, dt, iter)
    type(pes_mask_t),    intent(inout) :: mask
    type(namespace_t),   intent(in)    :: namespace
    type(space_t),       intent(in)    :: space
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(inout) :: st
    type(kpoints_t),     intent(in)    :: kpoints
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter

    integer :: idim, ist, ik
    type(cube_function_t):: cf1, cf2
    CMPLX, allocatable :: mf(:), psi(:)

    FLOAT :: time

    type(profile_t), save :: prof
    
    call profiling_in(prof, "PESMASK_calc")
    
    PUSH_SUB(pes_mask_calc)
    
    time = iter *dt

    if (time > mask%start_time) then ! record photoelectrons only after mask%start_time
      
      call zcube_function_alloc_RS(mask%cube, cf1, force_alloc = .true.) 
      call  cube_function_alloc_FS(mask%cube, cf1, force_alloc = .true.) 
      call zcube_function_alloc_RS(mask%cube, cf2, force_alloc = .true.)
      call  cube_function_alloc_FS(mask%cube, cf2, force_alloc = .true.)
      
      select case(mask%mode) 
      case(PES_MASK_MODE_MASK)
        if(mask%back_action .eqv. .true.) then
          SAFE_ALLOCATE(mf(1:mask%mesh%np_part))
        end if
      end select
      
      SAFE_ALLOCATE(psi(1:mask%mesh%np_part))
      
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim

            call states_elec_get_state(st, mask%mesh, idim, ist, ik, psi)
            
            cf1%zRs(:,:,:) = M_z0
            cf2%zRS(:,:,:) = M_z0
            cf1%Fs(:,:,:)  = M_z0
            cf2%Fs(:,:,:)  = M_z0
            
            call pes_mask_mesh_to_cube(mask, psi, cf1)
            
            select case(mask%mode)
              !----------------------------------------- 
              ! Mask Method
              !----------------------------------------
            case(PES_MASK_MODE_MASK)
              
              cf1%zRs = (M_ONE - mask%cM%zRs) * cf1%zRs                               ! cf1 =(1-M)*\Psi_A(x,t2)
              call pes_mask_X_to_K(mask, cf1%zRs, cf2%Fs)                             ! cf2 = \tilde{\Psi}_A(k,t2)


              if ( mask%filter_k ) then ! apply a filter to the Fourier transform to remove unwanted energies
                ASSERT(allocated(mask%Mk))
                cf2%Fs = cf2%Fs * mask%Mk 
              end if
              
              
              cf1%Fs(:,:,:) = mask%k(:,:,:, idim, ist, ik)                            ! cf1 = \Psi_B(k,t1)
              mask%k(:,:,:, idim, ist, ik) =  cf2%Fs(:,:,:)                           ! mask%k = \tilde{\Psi}_A(k,t2)
              call pes_mask_Volkov_time_evolution_wf(mask, space, mesh, kpoints, dt,iter-1,cf1%Fs, &   ! cf1 = \tilde{\Psi}_B(k,t2)
                                                     st%d%get_kpoint_index(ik))
                                                     
              mask%k(:,:,:, idim, ist, ik) =  mask%k(:,:,:, idim, ist, ik)&
                + cf1%Fs(:,:,:)      ! mask%k = \tilde{\Psi}_A(k,t2) + \tilde{\Psi}_B(k,t2)
              
              if(mask%back_action .eqv. .true.) then
                
                ! Apply Back-action to wavefunction in A
                call pes_mask_K_to_X(mask, cf1%Fs, cf2%zRs)                       ! cf2 = \Psi_B(x,t2)
                call pes_mask_cube_to_mesh(mask, cf2, mf)
                psi(1:mask%mesh%np) = psi(1:mask%mesh%np) + mf(1:mask%mesh%np)
                call states_elec_set_state(st, mask%mesh, idim, ist, ik, psi)
                
                ! Apply correction to wavefunction in B
                cf2%zRs= (mask%cM%zRs) * cf2%zRs                                     ! cf2 = M*\Psi_B(x,t1)
                call pes_mask_X_to_K(mask, cf2%zRs, cf1%Fs)
                
                mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) - cf1%Fs
                
              end if


              !----------------------------------------- 
              ! Passive Mask method
              !----------------------------------------
            case(PES_MASK_MODE_PASSIVE)
              
              
              cf1%zRs = (M_ONE-mask%cM%zRs) * cf1%zRs
              call pes_mask_X_to_K(mask, cf1%zRs, cf2%Fs) 
              
              mask%k(:,:,:, idim, ist, ik) = cf2%Fs(:,:,:)
              
              
            case default
              !Program should die before coming here
              write(message(1),'(a)') "PhotoElectroSpectrum = pes_mask. Unrecognized calculation mode." 
              call messages_fatal(1, namespace=namespace)
              
            end select
            
            
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(psi)

      call zcube_function_free_RS(mask%cube, cf1)
      call  cube_function_free_FS(mask%cube, cf1)
      call zcube_function_free_RS(mask%cube, cf2)
      call  cube_function_free_FS(mask%cube, cf2)

      select case(mask%mode) 
      case(PES_MASK_MODE_MASK)
        if(mask%back_action .eqv. .true.) then
          SAFE_DEALLOCATE_A(mf)
        end if
      end select
      
    end if ! time > mask%start_time
    

    if(mask%mode  ==  PES_MASK_MODE_MASK ) call pes_mask_apply_mask(mask, st)  !apply the mask to all the KS orbitals



    
    POP_SUB(pes_mask_calc)
    
    call profiling_out(prof)
  end subroutine pes_mask_calc
  
#include "pes_mask_out_inc.F90"
  
end module pes_mask_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
