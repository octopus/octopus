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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module pes_mask_m
  use comm_m
  use cube_function_m
  use cube_m
  use datasets_m
  use density_m
  use fft_m
  use fourier_space_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_function_m
  use io_m
  use lasers_m
  use math_m
  use mesh_m
  use mesh_cube_parallel_map_m
  use messages_m
  use mpi_m
#if defined(HAVE_NFFT) 
  use nfft_m
#endif
#if defined(HAVE_NETCDF)
  use netcdf
#endif
  use output_m
  use parser_m
  use poisson_m
  use profiling_m
  use simul_box_m
  use states_io_m
  use states_m
  use string_m
  use system_m
  use tdpsf_m
  use unit_m
  use unit_system_m
  use varinfo_m
  
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
    pes_mask_k_to_x
  
  type PES_mask_t
    CMPLX, pointer :: k(:,:,:,:,:,:) => NULL() !< The states in momentum space
    
    ! mesh- and cube-related stuff      
    integer          :: np                     !< number of mesh points associated with the mesh
    !< (either mesh%np or mesh%np_global)
    integer          :: ll(3)                  !< the size of the square mesh
    integer          :: fs_n_global(1:3)       !< the dimensions of the cube in fourier space
    integer          :: fs_n(1:3)              !< the dimensions of the local portion of the cube in fourier space
    integer          :: fs_istart(1:3)         !< where does the local portion of the cube start in fourier space
    
    FLOAT            :: spacing(3)       !< the spacing
    
    type(mesh_t), pointer  :: mesh             !< a pointer to the mesh
    type(cube_t)     :: cube                   !< the cubic mesh
    
    FLOAT, pointer :: ext_pot(:,:) => NULL()   !< external time-dependent potential i.e. the lasers
    
    FLOAT, pointer :: Mk(:,:,:) => NULL()      !< the momentum space filter
    type(cube_function_t) :: cM                !< the mask cube function
    FLOAT, pointer :: mask_R(:) => NULL()      !< the mask inner (component 1) and outer (component 2) radius
    integer        :: shape                    !< which mask function?
    
    FLOAT, pointer :: Lk(:) => NULL()          !< associate a k value to an cube index
    !< we implicitly assume k to be the same for all directions
    
    integer          :: resample_lev           !< resampling level
    integer          :: enlarge                !< Fourier space enlargement
    integer          :: enlarge_nfft           !< NFFT space enlargement
    !     integer          :: llr(MAX_DIM)           !< the size of the rescaled cubic mesh
    
    FLOAT :: start_time              !< the time we switch on the photoelectron detector   
    FLOAT :: energyMax 
    FLOAT :: energyStep 
    
    integer :: sw_evolve             !< choose the time propagator for the continuum wfs
    logical :: back_action           !< whether to enable back action from B to A
    logical :: add_psia              !< add the contribution of Psi_A in the buffer region to the output
    logical :: interpolate_out       !< whether to apply interpolation on the output files
    logical :: filter_k              !< whether to filter the wavefunctions in momentum space
    
    integer :: mode                  !< calculation mode
    integer :: pw_map_how            !< how to perform projection on plane waves
    
    type(fft_t)    :: fft            !< FFT plan
    
    type(tdpsf_t) :: psf             !< Phase-space filter struct reference
    
    type(mesh_cube_parallel_map_t) :: mesh_cube_map  !< The parallel map
    
    
  end type PES_mask_t

  integer, parameter ::       &
    FREE             =  1,    &    !< The scattering waves evolve in time as free plane waves
    VOLKOV           =  2,    &    !< The scattering waves evolve with exp(i(p-A(t)/c)^2*dt/2)
    CORRECTED1D      =  3,    &
    EMBEDDING1D      =  4,    &     
    VOLKOV_CORRECTED =  5
  
  integer, public, parameter ::        &
    PES_MASK_MODE_MASK         =   1,  &  
    PES_MASK_MODE_BACKACTION   =   2,  &  
    PES_MASK_MODE_PASSIVE      =   3,  &
    PES_MASK_MODE_PSF          =   4
  
  integer, parameter ::       &
    PW_MAP_INTEGRAL    =  1,  &    !< projection on outgoing waves by direct integration
    PW_MAP_FFT         =  2,  &    !< FFT on outgoing waves (1D only)
    PW_MAP_BARE_FFT    =  3,  &    !< FFT - normally from fftw3
    PW_MAP_TDPSF       =  4,  &    !< time-dependent phase-space filter
    PW_MAP_NFFT        =  5,  &    !< non-equispaced fft (NFFT)
    PW_MAP_PFFT        =  6        !< use PFFT
  
  integer, parameter ::       &
    M_SIN2            =  1,   &  
    M_STEP            =  2,   & 
    M_ERF             =  3    
  
  integer, parameter ::       &
    IN                =  1,   &  
    OUT               =  2
  
contains 


  ! ---------------------------------------------------------
  subroutine PES_mask_init(mask, mesh, sb, st, hm, max_iter,dt)
    type(PES_mask_t), target, intent(out) :: mask
    type(mesh_t),target,      intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    type(states_t),           intent(in)  :: st
    type(hamiltonian_t),      intent(in)  :: hm
    integer,                  intent(in)  :: max_iter
    FLOAT,                    intent(in)  :: dt
    
    type(block_t) :: blk
    
    integer :: il, it, ii, ll(3)
    FLOAT :: field(3)
    FLOAT :: DeltaE, MaxE, pCutOff
    FLOAT :: width 
    integer :: defaultMask,k1,k2,st1,st2, optimize_parity(3)
    logical :: optimize(3)
    FLOAT, allocatable  ::  XX(:,:)  
    
    PUSH_SUB(PES_mask_init)
    
    
    mask%mesh => mesh
    
    call messages_experimental('Photo-electron spectrum')  
    
    write(message(1),'(a,i1,a)') 'Info: Calculating PES using mask technique.'
    call messages_info(1)
    
    
    
    if(sb%box_shape /= SPHERE) then
      message(1) = 'PhotoElectronSpectrum = pes_mask requires BoxShape = sphere'
      message(2) = 'Modify the parameter and rerun.'
      call messages_fatal(2)
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
    !%Option psf_mode 4
    !% Phase-space filter. Implementation not complete.
    !%End
    call parse_integer(datasets_check('PESMaskMode'),PES_MASK_MODE_MASK,mask%mode)
    if(.not.varinfo_valid_option('PESMaskMode', mask%mode)) call input_error('PESMaskMode')
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
    

    !%Variable PESMaskPropagator 
    !%Type integer
    !%Default volkov
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Photoelectron waves time-propagation operator in momentum space.
    !%Option volkov 2
    !% Plane wave evolves with exp(i(p-A(t)/c)^2*dt/2).
    !%Option free 1
    !% Free plane-wave propagation.   
    !%End
    call parse_integer(datasets_check('PESMaskPropagator'),VOLKOV,mask%sw_evolve)
    if(.not.varinfo_valid_option('PESMaskPropagator',mask%sw_evolve)) call input_error('PESMaskPropagator')
    call messages_print_var_option(stdout, "PESMaskPropagator",mask%sw_evolve)
    
    !%Variable PESMaskStartTime 
    !%Type float
    !%Default -1.0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The time photoelectrons start to be recorded. In pump-probe simulations this allow to 
    !% get rid of unwanted ionization signal coming from the pump.
    !% NOTE: this will enforce the mask boundary conditions for all the times. 
    !%End
    call parse_float(datasets_check('PESMaskStartTime'),&
      units_to_atomic(units_inp%time, - M_ONE), mask%start_time)



    ! ========================================
    ! = Optimization and numerical stability =
    ! ========================================
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
    !% for a given simulation in addition to <tt>PESMaskEnlargeLev</tt> and 
    !% <tt>PESMaskNFFTEnlargeLev</tt>will help to converge the results.   
    !%
    !% NOTE: depending on the value of <tt>PESMaskMode</tt> <tt>PESMaskPlaneWaveProjection</tt>,
    !% may affect not only performance but also the time evolution of the density. 
    !%Option integral 1
    !% Direct integration_map.
    !%Option fft_out 2 
    !% FFT filtered in order to keep only outgoing waves. 1D only.  
    !%Option fft_map 3 
    !% FFT transform.
    !%Option tdpsf_map 4
    !% Time-dependent phase-space filter map.
    !%Option nfft_map 5
    !% Non-equispaced FFT map. 
    !%Option pfft_map 6
    !% Use PFFT libraries. 
    !%End
    call parse_integer(datasets_check('PESMaskPlaneWaveProjection'),PW_MAP_BARE_FFT,mask%pw_map_how)
    
    if(.not.varinfo_valid_option('PESMaskPlaneWaveProjection', mask%pw_map_how)) call input_error('PESMaskPlaneWaveProjection')
    call messages_print_var_option(stdout, "PESMaskPlaneWaveProjection", mask%pw_map_how)
    
#if !defined(HAVE_NFFT) 
    if (mask%pw_map_how ==  PW_MAP_NFFT) then
      message(1) = "PESMaskPlaneWaveProjection = nfft_map requires libnfft3. Recompile and try again." 
      call messages_fatal(1) 
    endif
#endif
    
#if !defined(HAVE_PFFT) 
    if (mask%pw_map_how ==  PW_MAP_PFFT) then
      message(1) = "PESMaskPlaneWaveProjection = pfft_map requires PFFT. Recompile and try again." 
      call messages_fatal(1) 
    endif
#endif
    
    if (mask%pw_map_how ==  PW_MAP_PFFT .and. (.not. mask%mesh%parallel_in_domains)) then
      message(1)= "Trying to use PESMaskPlaneWaveProjection = pfft_map with no domain parallelization."
      message(2)= "Projection method changed to more efficient fft_map."
      call messages_warning(2)
      mask%pw_map_how = PW_MAP_BARE_FFT
    end if
    
    !%Variable PESMaskEnlargeLev
    !%Type integer
    !%Default 0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Mask box enlargement level. Enlarges the mask bounding box by a factor 2**<tt>PESMaskEnlargeLev</tt>.
    !% This will avoid wavefunction wrapping at the boundaries.
    !%End

    call parse_integer(datasets_check('PESMaskEnlargeLev'),0,mask%enlarge)
    
    if ( mask%enlarge /= 0 ) then
      call messages_print_var_value(stdout, "PESMaskEnlargeLev", mask%enlarge)
    end if
    
    !%Variable PESMaskNFFTEnlargeLev
    !%Type integer
    !%Default 0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Mask box enlargement level. Enlarges the mask box by a factor 2**<tt>PESMaskEnlargeLev</tt>
    !% using NFFT. This way we add two points in each direction at a distance
    !% L=Lb*2**<tt>PESMaskEnlargeLev</tt> where <i>Lb</i> is the box size.
    !% 
    !% Note: the corresponding Fourier space is restricted by the same factor.
    !%End
    
    call parse_integer(datasets_check('PESMaskNFFTEnlargeLev'),0,mask%enlarge_nfft)
    
    if ( mask%enlarge_nfft /= 0 ) then
      call messages_print_var_value(stdout, "PESMaskNFFTEnlargeLev", mask%enlarge_nfft)
    end if
    
    
    mask%ll = 1
    mask%spacing = -M_ONE
    
    if(sb%mr_flag) then ! multiresolution 
      mask%spacing(1:sb%dim) = mesh%spacing(1:sb%dim)*2**(sb%hr_area%num_radii)       
      mask%ll(1:sb%dim) = int(M_TWO*sb%rsize/mask%spacing(1:sb%dim)) + 1
    else 
      mask%spacing(1:3) = mesh%spacing(1:3)
      mask%ll(1:3) = mesh%idx%ll(1:3)    
    end if
    
    !Enlarge the bounding box region
    mask%ll(1:sb%dim)= mask%ll(1:sb%dim)*M_TWO**mask%enlarge 
    
    
    select case(mask%pw_map_how)
    case(PW_MAP_FFT)
      call cube_init(mask%cube, mask%ll, mesh%sb, fft_type = FFT_COMPLEX, fft_library = FFTLIB_FFTW, nn_out = ll)
      mask%ll = ll !FFT optimization may change this values
      mask%fft = mask%cube%fft
      mask%np = mesh%np_part_global 
      
    case(PW_MAP_PFFT)
      ASSERT(mask%mesh%parallel_in_domains)
      call cube_init(mask%cube, mask%ll, mesh%sb, fft_type = FFT_COMPLEX, fft_library = FFTLIB_PFFT, nn_out = ll, &
        mpi_grp = mask%mesh%mpi_grp)
      !        print *,mpi_world%rank, "mask%mesh%mpi_grp%comm", mask%mesh%mpi_grp%comm, mask%mesh%mpi_grp%size
      !         print *,mpi_world%rank, "mask%cube%mpi_grp%comm", mask%cube%mpi_grp%comm, mask%cube%mpi_grp%size
      mask%ll(1) = mask%cube%fs_n(3)
      mask%ll(2) = mask%cube%fs_n(1)
      mask%ll(3) = mask%cube%fs_n(2)
      mask%fft = mask%cube%fft
      mask%np = mesh%np_part ! the mask is local
      if ( mask%mesh%parallel_in_domains .and. mask%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_init(mask%mesh_cube_map, mask%mesh, mask%cube)
      end if
      
    case(PW_MAP_BARE_FFT)
      call cube_init(mask%cube, mask%ll, mesh%sb, fft_type = FFT_COMPLEX, fft_library = FFTLIB_FFTW, nn_out = ll)
      mask%ll = ll 
      mask%fft = mask%cube%fft
      mask%np = mesh%np_part_global 
      
    case(PW_MAP_TDPSF)
      call cube_init(mask%cube, mask%ll, mesh%sb, fft_type = FFT_COMPLEX, fft_library = FFTLIB_FFTW, nn_out = ll)
      mask%ll = ll 
      mask%fft = mask%cube%fft
      mask%np = mesh%np_part_global 
      
    case(PW_MAP_NFFT)
      !FIXME: this part is a bit messy and should be integrated into cube_init
      
      !NFFT initialization
      
      ! we just add 2 points for the enlarged region
      if (mask%enlarge_nfft .ne. 1) mask%ll(1:sb%dim) = mask%ll(1:sb%dim) + 2 
      
#ifdef HAVE_NFFT    
      !Set NFFT defaults to values that are optimal for PES (at least for the cases I have tested)
      !These values are overridden by the NFFT options in the input file 
      mask%fft%nfft%set_defaults = .true.
      mask%fft%nfft%guru = .true.
      mask%fft%nfft%mm = 2 
      mask%fft%nfft%sigma = CNST(1.1)
      mask%fft%nfft%precompute = NFFT_PRE_PSI
#endif
      
      ! These options should not affect NFFT scheme  
      optimize(1:3) = .false.
      optimize(sb%periodic_dim+1:sb%dim) = .true.
      optimize_parity(1:sb%periodic_dim) = 0
      optimize_parity(sb%periodic_dim+1:sb%dim) = 1
      
      call fft_init(mask%fft, mask%ll, sb%dim, FFT_COMPLEX, FFTLIB_NFFT, optimize, optimize_parity )
      
      SAFE_ALLOCATE(XX(1:mask%ll(1),3))
      
      !Generate the NFFT-enlarged node grid
      if (mask%enlarge_nfft .gt. 0) then
        do ii=2, mask%ll(1)-1 
          XX(ii,1)= (ii - int(mask%ll(1)/2) -1)*mask%spacing(1)
        end do
        XX(1,1)= (-int(mask%ll(1)/2))*mask%spacing(1)*M_TWO**mask%enlarge_nfft 
        XX(mask%ll(1),1)= (int(mask%ll(1)/2))*mask%spacing(1)*M_TWO**mask%enlarge_nfft 
        
      else
        do ii=1, mask%ll(1) 
          XX(ii,1)= (ii - int(mask%ll(1)/2) -1)*mask%spacing(1)
        end do
      end if
      
      XX(:,2) = XX(:,1)
      XX(:,3) = XX(:,1)
      
      !Set the node points and precompute the NFFT plan
      call fft_init_stage1(mask%fft, XX)
      
      SAFE_DEALLOCATE_A(XX)
      
      call cube_init(mask%cube, mask%ll, mesh%sb)  
      SAFE_ALLOCATE(mask%cube%fft)
      mask%cube%fft = mask%fft
      call fft_get_dims(mask%cube%fft, mask%cube%rs_n_global, mask%cube%fs_n_global, mask%cube%rs_n, mask%cube%fs_n, &
        mask%cube%rs_istart, mask%cube%fs_istart)
      
      mask%np = mesh%np_part_global 
      
    case default 
      !Program should die before coming here
      write(message(1),'(a)') "PESMaskPlaneWaveProjection unrecognized option." 
      call messages_fatal(1)
      
    end select
    
    !Indices  
    
    mask%fs_istart = mask%cube%fs_istart 
    mask%fs_n = mask%cube%fs_n 
    mask%fs_n_global = mask%cube%fs_n_global 
    
    !   print *, mpi_world%rank, " mask%ll", mask%ll(1:3), "states -",st%st_start,st%st_end
    !Allocations
    call cube_function_null(mask%cM)    
    call zcube_function_alloc_RS(mask%cube, mask%cM, force_alloc = .true.)
    
    SAFE_ALLOCATE(mask%Lk(1:mask%fs_n_global(1)))
    
    st1 = st%st_start
    st2 = st%st_end
    k1 = st%d%kpt%start
    k2 = st%d%kpt%end   
    SAFE_ALLOCATE(mask%k(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3),1:st%d%dim,st1:st2,k1:k2))
    mask%k = M_z0
    
    
    ! generate the map between mesh and cube
    call  PES_mask_generate_Lk(mask) ! generate the physical momentum vector
    
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  Mask Function options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SAFE_ALLOCATE(mask%mask_R(2))
    
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
    call parse_integer(datasets_check('PESMaskShape'),defaultMask,mask%shape)
    if(.not.varinfo_valid_option('PESMaskShape', mask%shape)) call input_error('PESMaskShape')
    call messages_print_var_option(stdout, "PESMaskShape", mask%shape)
    
    !%Variable PESMaskSize
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Set the size of the mask function. 
    !% Here you can set the inner (R1) and outer (R2) radius by setting
    !% the block as follows:
    !%
    !%<tt>%PESMaskSize 
    !% <br>&nbsp;&nbsp; R1 | R2 
    !% <br>%</tt>
    !%
    !%End
    if (parse_block(datasets_check('PESMaskSize'), blk) < 0)then
      mask%mask_R(1)=mesh%sb%rsize/M_TWO
      mask%mask_R(2)=mesh%sb%rsize
      message(1)="PESMaskSize not specified. Using default values."
      call messages_info(1)
    else 
      call parse_block_float(blk, 0, 0, mask%mask_R(1), units_inp%length)
      call parse_block_float(blk, 0, 1, mask%mask_R(2), units_inp%length)
    end if
    
    if(mask%mask_R(2) .gt. mesh%sb%rsize)  mask%mask_R(2) = mesh%sb%rsize 
    
    
    
    write(message(1),'(a,es10.3,3a)') & 
      "Input: Mask  R1 = ", units_from_atomic(units_inp%length, mask%mask_R(1) ),&
      ' [', trim(units_abbrev(units_inp%length)), ']'
    write(message(2),'(a,es10.3,3a)') & 
      "             R2 = ", units_from_atomic(units_inp%length, mask%mask_R(2) ),&
      ' [', trim(units_abbrev(units_inp%length)), ']'
    call messages_info(2)
    
    
    if (mask%mode .eq. PES_MASK_MODE_PSF ) then 
      width = mask%mask_R(2)-mask%mask_R(1)
      call tdpsf_init(mask%psf,mask%fft, mesh, dt,width)
    else
      call PES_mask_generate_mask(mask,mesh)
    end if
    



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
    call parse_float(datasets_check('PESMaskFilterCutOff'),&
      units_to_atomic(units_inp%energy, - M_ONE), pCutOff)
    
    nullify(mask%Mk)
    mask%filter_k = .false.
    
    if(pCutOff > M_ZERO) then       
      call messages_print_var_value(stdout, "PESMaskFilterCutOff",pCutOff)
      mask%filter_k = .true.
      
      SAFE_ALLOCATE(mask%Mk(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3)))
      
      call PES_mask_generate_filter(mask,pCutOff)
    endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  Output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !%Variable PESMaskIncludePsiA
    !%Type logical
    !%Default false
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Add the contribution of \Psi_A in the mask region to the photo-electron spectrum.
    !% Literally adds the Fourier components of: 
    !% \Theta(r-R1)*\Psi_A(r)
    !% with \Theta being the Heaviside step function. 
    !% With this option PES will contain all the contributions starting from the inner 
    !% radius R1. Use this option to improve convergence with respect to the box size 
    !% and total simulation time. 
    !% Note: carefully choose R1 in order to avoid contributions from returning electrons. 
    !%End
    call parse_logical(datasets_check('PESMaskIncludePsiA'),.false.,mask%add_psia)
    if(mask%add_psia) then
      message(1)= "Input: Include contribution from Psi_A."
      call messages_info(1)
    end if
    
    
    !%Variable PESMaskSpectEnergyMax 
    !%Type float
    !%Default maxval(mask%Lk)**2/2
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The maximum energy for the PES spectrum.
    !%End
    MaxE = maxval(mask%Lk)**2/2
    call parse_float(datasets_check('PESMaskSpectEnergyMax'),&
      units_to_atomic(units_inp%energy,MaxE),mask%energyMax)
    call messages_print_var_value(stdout, "PESMaskSpectEnergyMax",mask%energyMax)



    !%Variable PESMaskSpectEnergyStep 
    !%Type float
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The PES spectrum energy step.
    !%End
    DeltaE = (mask%Lk(2)-mask%Lk(1))**2/M_TWO
    call parse_float(datasets_check('PESMaskSpectEnergyStep'),&
      units_to_atomic(units_inp%energy,DeltaE),mask%energyStep)
    call messages_print_var_value(stdout, "PESMaskSpectEnergyStep",mask%energyStep)
    
    !%Variable PESMaskOutputInterpolate 
    !%Type logical
    !%Default false
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% Use interpolation to evaluate the quantities in polar coordinates.
    !% NOTE: In 3D this is practically prohibitive in the present implementation.
    !% We suggest to use the postprocessing tool <tt>oct-photoelectron_spectrum</tt> in this case. 
    !%End
    call parse_logical(datasets_check('PESMaskOutputInterpolate'),.false.,mask%interpolate_out)
    if(mask%interpolate_out) then
      message(1)= "Input: output interpolation ENABLED."
      call messages_info(1)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  Set external fields 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SAFE_ALLOCATE(mask%ext_pot(0:max_iter,1:3))
    mask%ext_pot=M_ZERO

    if(mask%sw_evolve .eq. VOLKOV) then
      ! Precalculate the potential vector for all the simulation time
      do il = 1, hm%ep%no_lasers
        select case(laser_kind(hm%ep%lasers(il)))
        case(E_FIELD_MAGNETIC, E_FIELD_ELECTRIC)
          write(message(1),'(a)') 'PESMask works only with vector_potential unless in passive_mode.'
          call messages_warning(1)
        case(E_FIELD_VECTOR_POTENTIAL)
          do it = 1, max_iter
            field=M_ZERO
            call laser_field(hm%ep%lasers(il), field, it*dt)
            mask%ext_pot(it,:)= mask%ext_pot(it,:)-field(:) !Sum up all the fields (for some reason needs negative sign)
          end do
        end select
      end do
    else 
    end if



    POP_SUB(PES_mask_init)
  end subroutine PES_mask_init


  ! ---------------------------------------------------------
  subroutine PES_mask_end(mask)
    type(PES_mask_t), intent(inout) :: mask
    
    PUSH_SUB(PES_mask_end)
    
    SAFE_DEALLOCATE_P(mask%k)
    
    SAFE_DEALLOCATE_P(mask%ext_pot)
    SAFE_DEALLOCATE_P(mask%mask_R)
    SAFE_DEALLOCATE_P(mask%Lk)
    
    if(mask%mode == PES_MASK_MODE_PSF) then 
      call tdpsf_end(mask%psf)
    end if

    if ( mask%filter_k ) then
      SAFE_DEALLOCATE_P(mask%Mk)
    end if
    
    if (mask%mesh%parallel_in_domains .and. mask%cube%parallel_in_domains) then
      call mesh_cube_parallel_map_end(mask%mesh_cube_map)
    end if
    
    call zcube_function_free_RS(mask%cube, mask%cM)
    call cube_end(mask%cube)   
    
    POP_SUB(PES_mask_end)
  end subroutine PES_mask_end
  

  ! --------------------------------------------------------
  subroutine PES_mask_generate_Lk(mask)
    type(pes_mask_t), intent(inout) :: mask

    integer :: ii,nn
    FLOAT   :: temp

    PUSH_SUB(PES_mask_generate_Lk)

    temp = M_TWO * M_PI / (mask%fs_n_global(1) * mask%spacing(1))
    nn = mask%fs_n_global(1)

    do ii = 1, mask%fs_n_global(1)

      if (mask%pw_map_how .eq.  PW_MAP_NFFT) then
        !The Fourier space is shrunk by the factor M_TWO**mask%enlarge_nfft
        mask%Lk(ii) = (ii - nn/2 - 1)*temp/(M_TWO**mask%enlarge_nfft)
      else
        mask%Lk(ii) = pad_feq(ii,nn, .true.) * temp
      end if

    end do

    POP_SUB(PES_mask_generate_Lk)
  end subroutine PES_mask_generate_Lk

  ! ======================================
  !>  Generate the momentum-space filter 
  ! ======================================
  subroutine PES_mask_generate_filter(mask,cutOff)
    type(PES_mask_t), intent(inout) :: mask
    FLOAT,            intent(in)    ::cutOff
    
    
    integer :: kx,ky,kz, power
    FLOAT   :: KK(3), EE, Emax

    PUSH_SUB(PES_mask_generate_filter)
    
    mask%Mk = M_ZERO

    Emax = maxval(mask%Lk)**2 / M_TWO

    power = 8 

    do kx = 1, mask%ll(1)
      KK(1) = mask%Lk(kx) 
      do ky = 1, mask%ll(2)
        KK(2) = mask%Lk(ky)
        do kz = 1, mask%ll(3)
          KK(3) = mask%Lk(kz)
          
          EE = sum(KK(1:mask%mesh%sb%dim)**2) / M_TWO

          if ( EE > cutOff .and. EE < Emax ) then
            mask%Mk(kx,ky,kz) = M_ONE * sin((Emax-EE) * M_PI / (M_TWO * (cutOff) ))**power
          else 
            if( EE <= cutOff )  mask%Mk(kx,ky,kz) = M_ONE 
          end if
          
        end do
      end do
    end do


    POP_SUB(PES_mask_generate_filter)
  end subroutine PES_mask_generate_filter


  ! --------------------------------------------------------
  !>  Generate the mask function on the cubic mesh containing 
  !!  the simulation box
  ! ---------------------------------------------------------
  subroutine PES_mask_generate_mask(mask,mesh)
    type(PES_mask_t), intent(inout) :: mask
    type(mesh_t),     intent(in)    :: mesh



    PUSH_SUB(PES_mask_generate_mask)
    
    call PES_mask_generate_mask_function(mask,mesh, mask%shape, mask%mask_R)
    
    POP_SUB(PES_mask_generate_mask)
    
  end subroutine PES_mask_generate_mask

  ! --------------------------------------------------------
  !>  Generate the mask function on the cubic mesh containing 
  !!  the simulation box
  ! ---------------------------------------------------------
  subroutine PES_mask_generate_mask_function(mask,mesh, shape, R, mask_sq)
    type(PES_mask_t), intent(inout) :: mask
    type(mesh_t),     intent(in)    :: mesh
    integer,          intent(in)    :: shape
    FLOAT,            intent(in)    :: R(2)
    FLOAT, optional,  intent(out)   :: mask_sq(:,:,:)

    integer :: ip
    FLOAT   :: width
    FLOAT   :: xx(1:MAX_DIM), rr, dd
    CMPLX,allocatable :: mask_fn(:)
    logical :: local_

    PUSH_SUB(PES_mask_generate_mask_function)


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
          xx=mesh_x_global(mesh, ip) 
          rr = sqrt(dot_product(xx(1:mesh%sb%dim), xx(1:mesh%sb%dim)))
        end if
        dd = rr -  R(1) 
        if(dd .gt. M_ZERO ) then 
          if (dd .lt. width) then
            mask_fn(ip) = M_ONE * sin(dd * M_PI / (M_TWO * (width) ))**2
          else 
            mask_fn(ip) = M_ONE 
          end if
        end if
      end do

    case(M_STEP)
      do ip = 1, mask%np
        if(local_) then
          call mesh_r(mesh, ip, rr, coords=xx)
        else
          xx=mesh_x_global(mesh, ip) 
          rr = sqrt(dot_product(xx(1:mesh%sb%dim), xx(1:mesh%sb%dim)))
        end if
        dd = rr - R(1) 
        if(dd .gt. M_ZERO ) then 
          if (dd .lt. width) then
            mask_fn(ip) = M_ONE 
          else 
            mask_fn(ip) = M_ZERO
          end if
        end if
      end do
      
    case(M_ERF)
      
    case default
      message(1)="PhotoElectronSpectrum = pes_mask. Unrecognized mask type."
      call messages_fatal(1) 
    end select



    mask_fn(:) = M_ONE - mask_fn(:)

    call PES_mask_mesh_to_cube(mask, mask_fn, mask%cM, local = local_)

    if(present(mask_sq)) mask_sq = real(mask%cM%zRS)



    SAFE_DEALLOCATE_A(mask_fn)
    
    POP_SUB(PES_mask_generate_mask_function)
    
  end subroutine PES_mask_generate_mask_function


  ! --------------------------------------------------------
  subroutine PES_mask_apply_mask(mask, st)
    type(states_t),   intent(inout) :: st
    type(PES_mask_t), intent(in)    :: mask
    
    integer :: ik, ist, idim
    CMPLX, allocatable :: mmask(:)
    
    PUSH_SUB(PES_mask_apply_mask)
    SAFE_ALLOCATE(mmask(1:mask%mesh%np_part))
    
    call PES_mask_cube_to_mesh(mask, mask%cM, mmask)
    
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          st%zpsi(1:mask%mesh%np_part, idim, ist, ik) = st%zpsi(1:mask%mesh%np_part, idim, ist, ik)*mmask(1:mask%mesh%np_part)
        end do
      end do
    end do
    
    SAFE_DEALLOCATE_A(mmask)
    
    POP_SUB(PES_mask_apply_mask)
  end subroutine PES_mask_apply_mask






  ! ---------------------------------------------------------
  !>  Propagate in time a wavefunction in momentum space with 
  !!  the Volkov Hamiltonian 
  !!\f[     
  !!     wf(p,t+dt)=exp(-i*Hv*dt)*wf(p,t)
  !!\f]
  !!  with 
  !!\f[
  !!     Hv=(p^2-A)^2/2
  !!\f]
  !! \note velocity gauge is implied 
  ! ---------------------------------------------------------
  subroutine PES_mask_Volkov_time_evolution_wf(mask, mesh, dt, iter, wf)
    type(PES_mask_t), intent(in)    :: mask
    type(mesh_t),     intent(in)    :: mesh
    FLOAT,            intent(in)    :: dt
    integer,          intent(in)    :: iter
    CMPLX,            intent(inout) :: wf(:,:,:)
    
    integer ::  ix, iy, iz
    FLOAT :: vec
    FLOAT :: KK(1:3)

    PUSH_SUB(PES_mask_Volkov_time_evolution_wf)


    ! propagate wavefunction in momentum space in presence of a td field (in the velocity gauge)
    if(mask%cube%fft%library == FFTLIB_PFFT) then !PFFT FS indices are transposed

      do ix = 1, mask%fs_n(1)
        KK(3) = mask%Lk(ix + mask%fs_istart(1) - 1)
        do iy = 1, mask%fs_n(2)
          KK(1) = mask%Lk(iy + mask%fs_istart(2) - 1)
          do iz = 1, mask%fs_n(3)
            KK(2) = mask%Lk(iz + mask%fs_istart(3) - 1)
            
            vec = sum(( KK(1:mesh%sb%dim) - mask%ext_pot(iter,1:mesh%sb%dim)/P_C)**2) / M_TWO
            wf(iz, ix, iy) = wf(iz, ix, iy) * exp(-M_zI * dt * vec)
            
            
          end do
        end do
      end do
    else
      
      do ix = 1, mask%ll(1)
        KK(1) = mask%Lk(ix + mask%fs_istart(1) - 1)
        do iy = 1, mask%ll(2)
          KK(2) = mask%Lk(iy + mask%fs_istart(2) - 1)
          do iz = 1, mask%ll(3)
            KK(3) = mask%Lk(iz + mask%fs_istart(3) - 1)
            vec = sum(( KK(1:mesh%sb%dim) - mask%ext_pot(iter,1:mesh%sb%dim)/P_C)**2) / M_TWO
            wf(ix, iy, iz) = wf(ix, iy, iz) * exp(-M_zI * dt * vec)
          end do
        end do
      end do
      
    end if
    

    POP_SUB(PES_mask_Volkov_time_evolution_wf)
  end subroutine PES_mask_Volkov_time_evolution_wf



  !---------------------------------------------------------
  subroutine fft_X_to_K(mask, wfin, wfout)
    type(PES_mask_t), intent(in)  :: mask
    CMPLX,            intent(in)  :: wfin(:,:,:)
    CMPLX,            intent(out) :: wfout(:,:,:)
    
    integer :: ix, iy, iz, ixx(3)
    integer :: ll(3)
    CMPLX, allocatable :: wftmp(:,:,:),wfPlus(:,:,:),wfMinus(:,:,:) 
    
    PUSH_SUB(fft_X_to_K)
    
    
    ll(1:3) = mask%ll(1:3)
    SAFE_ALLOCATE(wftmp(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
    SAFE_ALLOCATE(wfPlus(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
    SAFE_ALLOCATE(wfMinus(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))


    wfout = M_z0
    wfPlus= M_z0
    wfMinus= M_z0
    wftmp = M_z0


    do ix = 1, mask%ll(1)
      do iy = 1, mask%ll(2)
        do iz = 1, mask%ll(3)
          
          if(ix <= ll(1)/2+1)   then
            wfMinus(ix, iy, iz)  =  wfin(ix, iy, iz)
            wfPlus(ix,iy,iz) = M_z0
          else
            wfPlus(ix, iy, iz)   =  wfin(ix, iy, iz)
            wfMinus(ix,iy,iz)  = M_z0
          end if
          
        end do
      end do
    end do


    wftmp =  M_z0
    call zfft_forward(mask%fft, wfPlus,wftmp)
    wfPlus=wftmp
    
    wftmp =  M_z0
    call zfft_forward(mask%fft, wfMinus,wftmp)
    wfMinus=wftmp
    
    do ix = 1, mask%ll(1)
      ixx(1) = pad_feq(ix, mask%ll(1), .true.)
      do iy = 1, mask%ll(2)
        ixx(2) = pad_feq(iy, mask%ll(2), .true.)
        do iz = 1, mask%ll(3)
          ixx(3) = pad_feq(iz, mask%ll(3), .true.)
          
          if(ixx(1) > 0 ) then
            wfMinus(ix, iy, iz) = M_z0
          else
            wfPlus(ix, iy, iz) = M_z0
          end if
          
        end do
      end do
    end do
    
    wfout= wfPlus + wfMinus
    
    
    SAFE_DEALLOCATE_A(wftmp)
    SAFE_DEALLOCATE_A(wfPlus)
    SAFE_DEALLOCATE_A(wfMinus)


    POP_SUB(fft_X_to_K)
  end subroutine fft_X_to_K

  ! ------------------------------------------------
  subroutine fft_K_to_X(mask,wfin,wfout,inout)
    type(PES_mask_t), intent(in)  :: mask
    CMPLX,            intent(in)  :: wfin(:,:,:)
    CMPLX,            intent(out) :: wfout(:,:,:)
    integer,          intent(in)  :: inout
    
    integer :: ix, iy, iz, ixx(3)
    integer :: ll(3)
    CMPLX, allocatable :: wftmp(:,:,:),wfPlus(:,:,:),wfMinus(:,:,:) 

    PUSH_SUB(fft_K_to_X)
    

    ll(1:3) = mask%ll(1:3)
    SAFE_ALLOCATE(wftmp(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
    SAFE_ALLOCATE(wfPlus(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
    SAFE_ALLOCATE(wfMinus(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))


    wfout = M_z0
    wfPlus= M_z0
    wfMinus= M_z0
    wftmp = M_z0



    do ix = 1, mask%ll(1)
      ixx(1) = pad_feq(ix, mask%ll(1), .true.)
      do iy = 1, mask%ll(2)
        ixx(2) = pad_feq(iy, mask%ll(2), .true.)
        do iz = 1, mask%ll(3)
          ixx(3) = pad_feq(iz, mask%ll(3), .true.)
          
          select case(inout)
          case(IN)
            if(ixx(1) < 0 ) then
              wfPlus(ix,iy,iz) = wfin(ix,iy,iz)
              wfMinus(ix, iy, iz) = M_z0
            else
              wfPlus(ix, iy, iz) = M_z0
              wfMinus(ix,iy,iz) = wfin(ix,iy,iz)
            end if
          case(OUT)
            if(ixx(1) > 0 ) then
              wfPlus(ix,iy,iz) = wfin(ix,iy,iz)
              wfMinus(ix, iy, iz) = M_z0
            else
              wfPlus(ix, iy, iz) = M_z0
              wfMinus(ix,iy,iz) = wfin(ix,iy,iz)
            end if

          end select

        end do
      end do
    end do

    
    wftmp =  M_z0
    call zfft_backward(mask%fft, wfPlus,wftmp)
    wfPlus=wftmp

    wftmp =  M_z0
    call zfft_backward(mask%fft, wfMinus,wftmp)
    wfMinus=wftmp



    do ix = 1, mask%ll(1)
      do iy = 1, mask%ll(2)
        do iz = 1, mask%ll(3)
          select case(inout)
          case(IN)
            if(ix <= ll(1)/2+1)   then
              wfPlus(ix,iy,iz) = M_z0
            else
              wfMinus(ix,iy,iz)  = M_z0
            end if
          case (OUT)
            if(ix <= ll(1)/2+1)   then
              wfPlus(ix,iy,iz) = M_z0
            else
              wfMinus(ix,iy,iz)  = M_z0
            end if
          end select
        end do
      end do
    end do


    wfout= wfPlus + wfMinus


    SAFE_DEALLOCATE_A(wftmp)
    SAFE_DEALLOCATE_A(wfPlus)
    SAFE_DEALLOCATE_A(wfMinus)


    POP_SUB(fft_K_to_X)
  end subroutine fft_K_to_X

  !---------------------------------------------------------
  subroutine integral_X_to_K(mask,mesh,wfin,wfout)
    type(PES_mask_t), intent(in)  :: mask
    type(mesh_t),     intent(in)  :: mesh
    CMPLX,            intent(in)  :: wfin(:,:,:)
    CMPLX,            intent(out) :: wfout(:,:,:)
    
    integer :: ix, iy, iz, ixx(3)
    FLOAT   :: temp(3)
    FLOAT   :: k_dot_r
    integer :: kx,ky,kz,ikk(3)
    
    PUSH_SUB(integral_X_to_K)
    
    wfout = M_z0
    
    temp(:) = M_TWO * M_PI / (mask%ll(:) * mask%spacing(:))
    do kx = 1, mask%ll(1)
      ikk(1) = pad_feq(kx, mask%ll(1), .true.)
      do ky = 1, mask%ll(2)
        ikk(2) = pad_feq(ky, mask%ll(2), .true.)
        do kz = 1, mask%ll(3)
          ikk(3) = pad_feq(kz, mask%ll(3), .true.)
          
          do ix = 1, mask%ll(1)
            ixx(1) =  ix-mask%ll(1)/2+1
            do iy = 1, mask%ll(2)
              ixx(2) = iy-mask%ll(2)/2+1
              do iz = 1, mask%ll(3)
                ixx(3) = iz-mask%ll(3)/2+1
                
                
                k_dot_r=sum(ikk(1:mesh%sb%dim)*temp(1:mesh%sb%dim)*ixx(1:mesh%sb%dim)*mask%spacing(1:mesh%sb%dim))  
                

                !               if(k_dot_r .gt. 0) then 
                wfout(kx,ky,kz)=wfout(kx,ky,kz)+wfin(ix,iy,iz)*exp(-M_zI*k_dot_r)
                !               end if 
                
              end do
            end do
          end do
          
        end do
      end do
    end do
    
    !  wfout=wfout*(mask%spacing(1)*mask%spacing(2)*mask%spacing(3))/((M_TWO*M_PI)**(mesh%sb%dim/2))
    
    
    POP_SUB(integral_X_to_K)
  end subroutine integral_X_to_K

  ! ------------------------------------------------
  subroutine integral_K_to_X(mask,mesh,wfin,wfout)
    type(PES_mask_t), intent(in)  :: mask
    type(mesh_t),     intent(in)  :: mesh
    CMPLX,            intent(in)  :: wfin(:,:,:)
    CMPLX,            intent(out) :: wfout(:,:,:)
    
    integer :: ix, iy, iz, ixx(3)
    FLOAT   :: temp(3)
    FLOAT   :: k_dot_r
    integer :: kx,ky,kz,ikk(3)
    
    PUSH_SUB(integral_K_to_X)
    
    wfout = M_z0
    
    
    temp(:) = M_TWO * M_PI / (mask%ll(:) * mask%spacing(:))
    
    do ix = 1, mask%ll(1)
      ixx(1) = ix-mask%ll(1)/2+1
      do iy = 1, mask%ll(2)
        ixx(2) = iy-mask%ll(2)/2+1
        do iz = 1, mask%ll(3)
          ixx(3) = iz-mask%ll(3)/2+1
          
          do kx = 1, mask%ll(1)
            ikk(1) = pad_feq(kx, mask%ll(1), .true.)
            do ky = 1, mask%ll(2)
              ikk(2) = pad_feq(ky, mask%ll(2), .true.)
              do kz = 1, mask%ll(3)
                ikk(3) = pad_feq(kz, mask%ll(3), .true.)
                
                k_dot_r=sum(ikk(1:mesh%sb%dim)*temp(1:mesh%sb%dim)*ixx(1:mesh%sb%dim)*mask%spacing(1:mesh%sb%dim))  
                
                !               if(k_dot_r .gt. 0) then 
                wfout(ix,iy,iz)=wfout(ix,iy,iz)+wfin(kx,ky,kz)*exp(M_zI*k_dot_r)
                !               end if 
                
              end do
            end do
          end do
          
        end do
      end do
    end do
    
    !  wfout=wfout*(temp(1)*temp(2)*temp(3))/((M_TWO*M_PI)**(mesh%sb%dim/2))
    
    POP_SUB(integral_K_to_X)
  end subroutine integral_K_to_X



  !---------------------------------------------------------
  !> Project the wavefunction on plane waves
  !---------------------------------------------------------
  subroutine PES_mask_X_to_K(mask,mesh,wfin,wfout)
    type(PES_mask_t), intent(in)  :: mask
    type(mesh_t),     intent(in)  :: mesh
    CMPLX,              intent(in):: wfin(:,:,:)
    CMPLX,             intent(out):: wfout(:,:,:)
    
    type(profile_t), save :: prof
    type(cube_function_t) :: cf_tmp
    
    call profiling_in(prof, "PESMASK_X_to_K")
    
    PUSH_SUB(PES_mask_X_to_K)
    
    wfout=M_z0
    
    select case(mask%pw_map_how)
    case(PW_MAP_INTEGRAL)
      call integral_X_to_K(mask,mesh,wfin, wfout)
      
    case(PW_MAP_FFT)
      call fft_X_to_K(mask, wfin, wfout)
      
    case(PW_MAP_BARE_FFT)
      call zfft_forward(mask%cube%fft, wfin, wfout)
      
    case(PW_MAP_TDPSF)
      call tdpsf_X_to_K(mask%psf, wfin, wfout)
      
    case(PW_MAP_NFFT)
      call zfft_forward(mask%cube%fft, wfin, wfout)
      
    case(PW_MAP_PFFT)
      call cube_function_null(cf_tmp)    
      call zcube_function_alloc_RS(mask%cube, cf_tmp)
      call cube_function_alloc_fs(mask%cube, cf_tmp)
      cf_tmp%zRs = wfin
      call zfft_forward(mask%cube%fft, cf_tmp%zRs, cf_tmp%fs)
      wfout = cf_tmp%fs
      call zcube_function_free_RS(mask%cube, cf_tmp)
      call cube_function_free_fs(mask%cube, cf_tmp)
      
    case default
      
    end select
    
    
    POP_SUB(PES_mask_X_to_K)
    call profiling_out(prof)
    
  end subroutine PES_mask_X_to_K
  
  ! ------------------------------------------------
  subroutine PES_mask_K_to_X(mask,mesh,wfin,wfout)
    type(PES_mask_t), intent(in)  :: mask
    type(mesh_t),     intent(in)  :: mesh
    CMPLX,            intent(in)  :: wfin(:,:,:)
    CMPLX,            intent(out) :: wfout(:,:,:)
    
    type(profile_t), save :: prof
    type(cube_function_t) :: cf_tmp
    
    call profiling_in(prof, "PESMASK_K_toX")
    
    PUSH_SUB(PES_mask_K_to_X)
    
    wfout=M_z0
    
    select case(mask%pw_map_how)
    case(PW_MAP_INTEGRAL)
      call integral_K_to_X(mask,mesh,wfin,wfout)
      
    case(PW_MAP_FFT)
      call fft_K_to_X(mask,wfin,wfout,OUT)
      
    case(PW_MAP_BARE_FFT)
      call zfft_backward(mask%cube%fft, wfin,wfout)
      
    case(PW_MAP_TDPSF)
      call tdpsf_K_to_X(mask%psf, wfin,wfout)
      
    case(PW_MAP_NFFT)
      call zfft_backward(mask%cube%fft, wfin,wfout)
      
    case(PW_MAP_PFFT)
      
      call cube_function_null(cf_tmp)    
      call zcube_function_alloc_RS(mask%cube, cf_tmp)
      call cube_function_alloc_fs(mask%cube, cf_tmp)
      cf_tmp%fs  = wfin
      call zfft_backward(mask%cube%fft, cf_tmp%fs, cf_tmp%zRs)
      wfout = cf_tmp%zRs
      call zcube_function_free_RS(mask%cube, cf_tmp)
      call cube_function_free_fs(mask%cube, cf_tmp)
      

    case default

    end select


    POP_SUB(PES_mask_K_to_X)
    
    call profiling_out(prof)

  end subroutine PES_mask_K_to_X

  !---------------------------------------------------------
  subroutine PES_mask_mesh_to_cube(mask, mf, cf, local)
    type(PES_mask_t),      intent(in) :: mask
    CMPLX,                 intent(in) :: mf(:)
    type(cube_function_t), intent(out):: cf
    logical, optional,     intent(in) :: local
    
    logical :: local_
    
    PUSH_SUB(PES_mask_mesh_to_cube)
    
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
    
    POP_SUB(PES_mask_mesh_to_cube)
  end subroutine PES_mask_mesh_to_cube


  !---------------------------------------------------------
  subroutine PES_mask_cube_to_mesh(mask, cf, mf)
    type(PES_mask_t),      intent(in) :: mask
    CMPLX,                 intent(out):: mf(:)
    type(cube_function_t), intent(in) :: cf
    
    PUSH_SUB(PES_mask_cube_to_mesh)
    
    if (mask%cube%parallel_in_domains) then
      call zcube_to_mesh_parallel(mask%cube, cf, mask%mesh, mf, mask%mesh_cube_map)
    else
      if(mask%mesh%parallel_in_domains) then
        call zcube_to_mesh(mask%cube, cf, mask%mesh, mf, local = .true.)
      else 
        call zcube_to_mesh(mask%cube, cf, mask%mesh, mf)
      end if
    end if
    
    POP_SUB(PES_mask_cube_to_mesh)
  end subroutine PES_mask_cube_to_mesh


  !---------------------------------------------------------
  !
  !            Performs all the dirty work 
  !
  !---------------------------------------------------------
  subroutine PES_mask_calc(mask, mesh, st, dt, iter)
    type(PES_mask_t),    intent(inout) :: mask
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter

    integer :: idim, ist, ik
    type(cube_function_t):: cf1,cf2,cf3,cf4
    CMPLX, allocatable :: mf(:)

    FLOAT :: time

    type(profile_t), save :: prof
    
    call profiling_in(prof, "PESMASK_calc")
    
    PUSH_SUB(PES_mask_calc)
    
    time = iter *dt

    if (time > mask%start_time) then ! record photoelectrons only after mask%start_time
      
      call cube_function_null(cf1)    
      call zcube_function_alloc_RS(mask%cube, cf1, force_alloc = .true.) 
      call  cube_function_alloc_FS(mask%cube, cf1, force_alloc = .true.) 
      call cube_function_null(cf2)    
      call zcube_function_alloc_RS(mask%cube, cf2, force_alloc = .true.)
      call  cube_function_alloc_FS(mask%cube, cf2, force_alloc = .true.)
      
      select case(mask%mode) 
      case(PES_MASK_MODE_MASK)
        if(mask%back_action .eqv. .true.) then
          SAFE_ALLOCATE(mf(1:mask%mesh%np_part))
        end if
      case(PES_MASK_MODE_PSF)
        call cube_function_null(cf3)    
        call zcube_function_alloc_RS(mask%cube, cf3, force_alloc = .true.) 
        call  cube_function_alloc_FS(mask%cube, cf3, force_alloc = .true.) 
        call cube_function_null(cf4)    
        call zcube_function_alloc_RS(mask%cube, cf4, force_alloc = .true.) 
        call  cube_function_alloc_FS(mask%cube, cf4, force_alloc = .true.) 
      end select
      

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            
            cf1%zRs(:,:,:) = M_z0
            cf2%zRS(:,:,:) = M_z0
            cf1%Fs(:,:,:)  = M_z0
            cf2%Fs(:,:,:)  = M_z0
            
            call PES_mask_mesh_to_cube(mask, st%zpsi(:, idim, ist, ik), cf1)
            
            select case(mask%mode)
              !----------------------------------------- 
              ! Mask Method
              !----------------------------------------
            case(PES_MASK_MODE_MASK)
              
              cf1%zRs = (M_ONE - mask%cM%zRs) * cf1%zRs                               ! cf1 =(1-M)*U(t2,t1)*\Psi_A(x,t1)
              call PES_mask_X_to_K(mask,mesh,cf1%zRs,cf2%Fs)                          ! cf2 = \tilde{\Psi}_A(k,t2)


              if ( mask%filter_k ) then ! apply a filter to the Fourier transform to remove unwanted energies
                ASSERT(associated(mask%Mk))
                cf2%Fs = cf2%Fs * mask%Mk 
              end if
              
              
              cf1%Fs(:,:,:) = mask%k(:,:,:, idim, ist, ik)                            ! cf1 = \Psi_B(k,t1)
              mask%k(:,:,:, idim, ist, ik) =  cf2%Fs(:,:,:)                           ! mask%k = \tilde{\Psi}_A(k,t2)
              call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter-1,cf1%Fs)     ! cf1 = \tilde{\Psi}_B(k,t2)
              
              mask%k(:,:,:, idim, ist, ik) =  mask%k(:,:,:, idim, ist, ik)&
                + cf1%Fs(:,:,:)      ! mask%k = \tilde{\Psi}_A(k,t2) + \tilde{\Psi}_B(k,t2)
              
              if(mask%back_action .eqv. .true.) then
                
                ! Apply Back-action to wavewunction in A
                call PES_mask_K_to_X(mask,mesh,cf1%Fs,cf2%zRs)                       ! cf2 = \Psi_B(x,t2)
                call PES_mask_cube_to_mesh(mask, cf2, mf)  
                st%zpsi(:, idim, ist, ik) = st%zpsi(:, idim, ist, ik) + mf
                
                
                ! Apply correction to wavefunciton in B
                cf2%zRs= (mask%cM%zRs) * cf2%zRs                                     ! cf2 = M*\Psi_B(x,t1)
                call PES_mask_X_to_K(mask,mesh,cf2%zRs,cf1%Fs)
                
                mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) - cf1%Fs
                
              end if


              !----------------------------------------- 
              ! Passive Mask method
              !----------------------------------------
            case(PES_MASK_MODE_PASSIVE)
              
              
              cf1%zRs = (M_ONE-mask%cM%zRs) * cf1%zRs
              call PES_mask_X_to_K(mask,mesh,cf1%zRs,cf2%Fs) 
              
              mask%k(:,:,:, idim, ist, ik) = cf2%Fs(:,:,:)
              
              
              !----------------------------------------- 
              ! Phase Space Filter
              !----------------------------------------
              
            case(PES_MASK_MODE_PSF)
              
              call tdpsf_X_to_K(mask%psf,cf1%zRs,cf3%Fs)
              call tdpsf_K_to_X(mask%psf,cf3%Fs,cf2%zRs)
              
              cf2%zRs = cf1%zRs - cf2%zRs
              
              if(mask%back_action) then
                cf1%Fs = mask%k(:,:,:, idim, ist, ik) - cf4%Fs
                call PES_mask_K_to_X(mask,mesh,cf1%Fs, mask%k(:,:,:, idim, ist, ik))
                cf2%Fs = cf2%Fs + mask%k(:,:,:, idim, ist, ik)
              end if
              
              !substitute the KS wf with the filtered one
              call PES_mask_cube_to_mesh(mask, cf2, st%zpsi(:, idim, ist, ik))
              !the out-going part of the wf
              cf4%zRs = cf1%zRs - cf2%zRs
              call PES_mask_X_to_K(mask,mesh,cf4%zRs,cf3%Fs) 
              
              call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter,mask%k(:,:,:, idim, ist, ik) )
              cf4%Fs = mask%k(:,:,:, idim, ist, ik)            
              
              mask%k(:,:,:, idim, ist, ik) =  cf4%Fs + cf3%Fs
              
              
            case default
              !Program should die before coming here
              write(message(1),'(a)') "PhotoElectroSpectrum = pes_mask. Unrecognized calculation mode." 
              call messages_fatal(1)
              
            end select
            
            
          end do
        end do
      end do

      call zcube_function_free_RS(mask%cube, cf1)
      call  cube_function_free_FS(mask%cube, cf1)
      call zcube_function_free_RS(mask%cube, cf2)
      call  cube_function_free_FS(mask%cube, cf2)

      select case(mask%mode) 
      case(PES_MASK_MODE_MASK)
        if(mask%back_action .eqv. .true.) then
          SAFE_DEALLOCATE_A(mf)
        end if
      case(PES_MASK_MODE_PSF)
        call zcube_function_free_RS(mask%cube, cf3)
        call  cube_function_free_FS(mask%cube, cf3)
        call zcube_function_free_RS(mask%cube, cf4)
        call  cube_function_free_FS(mask%cube, cf4)
      end select
      
    end if ! time > mask%start_time
    

    if(mask%mode .eq. PES_MASK_MODE_MASK ) call PES_mask_apply_mask(mask, st)  !apply the mask to all the KS orbitals



    
    POP_SUB(PES_mask_calc)
    
    call profiling_out(prof)
  end subroutine PES_mask_calc
end module pes_mask_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
