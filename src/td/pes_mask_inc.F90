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

! ---------------------------------------------------------
subroutine PES_mask_init(mask, mesh, sb, st, hm, max_iter,dt)
  type(PES_mask_t),         intent(out) :: mask
  type(mesh_t),target,      intent(in)  :: mesh
  type(simul_box_t),        intent(in)  :: sb
  type(states_t),           intent(in)  :: st
  type(hamiltonian_t),      intent(in)  :: hm
  integer,                  intent(in)  :: max_iter
  FLOAT,                    intent(in)  :: dt

  type(block_t) :: blk

  integer :: il,it,ii
  FLOAT :: field(MAX_DIM)
  FLOAT :: DeltaE,MaxE,MaxDR, pCutOff
  FLOAT :: width 
  integer :: dir, defaultMask,k1,k2,st1,st2, optimize_parity(3)
  logical :: optimize(3)
  FLOAT, allocatable  :: X1(:),X2(:),X3(:)  

  PUSH_SUB(PES_mask_init)

 
  mask%mesh => mesh

  call messages_experimental('Photo-electron spectrum')  

  write(message(1),'(a,i1,a)') 'Info: Calculating PES using mask technique.'
  call messages_info(1)



  if(mesh%parallel_in_domains .and. st%parallel_in_states) then
    write(message(1),'(a)') "PES_mask: simultaneous parallelization on mesh and states not supported"
    write(message(2),'(a)') "Modify ParallelizationStrategy and rerun." 
    call messages_fatal(2) 
  end if    


  if(sb%box_shape /= SPHERE) then
     message(1) = 'PhotoElectronSpectrum = pes_mask requires BoxShape = sphere'
     message(2) = 'Modify the parameter and rerun.'
     call messages_warning(1)
   end if 


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Calculation mode
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !%Variable PESMaskMode
  !%Type integer
  !%Default mask_mode
  !%Section Time-Dependent::PES
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
  call parse_integer(datasets_check('PESMaskMode'),MODE_MASK,mask%mode)
  if(.not.varinfo_valid_option('PESMaskMode', mask%mode)) call input_error('PESMaskMode')
  call messages_print_var_option(stdout, "PESMaskMode", mask%mode)

  select case(mask%mode)
    case(MODE_PASSIVE)
    defaultMask = M_STEP   
    mask%back_action = .false.   

    case(MODE_BACKACTION)
    defaultMask = M_SIN2   
    mask%back_action = .true.
    mask%mode = MODE_MASK
    case default
    defaultMask = M_SIN2
    mask%back_action = .false.   
  end select


  !%Variable PESMaskPropagator 
  !%Type integer
  !%Default volkov
  !%Section Time-Dependent::PES
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



! ========================================
! = Optimization and numerical stability =
! ========================================
  !%Variable PESMaskPlaneWaveProjection
  !%Type integer
  !%Default fft_map
  !%Section Time-Dependent::PES
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

  
  !%Variable PESMaskEnlargeLev
  !%Type integer
  !%Default 0
  !%Section Time-Dependent::PES
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
  !%Section Time-Dependent::PES
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

  if(sb%mr_flag) then
    mask%spacing(1:sb%dim) = mesh%spacing(1:sb%dim)*2**(sb%hr_area%num_radii)       
    mask%ll(1:sb%dim) = int(M_TWO*sb%rsize/mask%spacing(1:sb%dim)) + 1
  else 
    mask%spacing = mesh%spacing
    mask%ll = mesh%idx%ll    
  end if 

  !Enlarge the bounding box region
  mask%ll(1:sb%dim)= mask%ll(1:sb%dim)*M_TWO**mask%enlarge 
   
  if (mask%pw_map_how .ne.  PW_MAP_NFFT) then

    ! allocate FFTs in case they are not allocated yet
    optimize(1:sb%periodic_dim) = .false.
    optimize(sb%periodic_dim+1:sb%dim) = .true.
    optimize_parity(1:sb%periodic_dim) = 0
    optimize_parity(sb%periodic_dim+1:sb%dim) = 1
    call fft_init(mask%fft, mask%ll, sb%dim, FFT_COMPLEX, FFTLIB_FFTW, optimize, optimize_parity)

  else
#if defined(HAVE_NFFT) 


    !NFFT initialization


    ! we just add 2 points for the enlarged region
    if (mask%enlarge_nfft .ne. 1) mask%ll(1:sb%dim) = mask%ll(1:sb%dim) + 2 
    
    call nfft_init(mask%ll,sb%dim,mask%ll(1) ,nfft_complex, mask%nfft,optimize = .true.)
    
    SAFE_ALLOCATE(X1(1:mask%ll(1)))
    
    !generate the grid
    if (mask%enlarge_nfft .gt. 0) then
      do ii=2, mask%ll(1)-1 
        X1(ii)= (ii - int(mask%ll(1)/2) -1)*mask%spacing(1)
      end do
      X1(1)= (-int(mask%ll(1)/2))*mask%spacing(1)*M_TWO**mask%enlarge_nfft 
      X1(mask%ll(1))= (int(mask%ll(1)/2))*mask%spacing(1)*M_TWO**mask%enlarge_nfft 
      
    else
      do ii=1, mask%ll(1) 
        X1(ii)= (ii - int(mask%ll(1)/2) -1)*mask%spacing(1)
      end do
    end if
    
    !!Set the node points and precompute the NFFT plan
    call nfft_precompute(mask%nfft, X1,X1,X1)

    SAFE_DEALLOCATE_A(X1)
#endif
  end if

  !!ALLOCATIONS

  mask%np = mesh%np_part_global !we do not divide the cube objects in this implementation 

 

  SAFE_ALLOCATE(mask%Lxyz_inv(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3)))
  call cube_init(mask%cube, mask%ll, mesh%sb)  

  SAFE_ALLOCATE(mask%M(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3)))
	call cube_function_null(mask%cM)    
	call dcube_function_alloc_RS(mask%cube, mask%cM)

  SAFE_ALLOCATE(mask%Lk(1:mask%ll(1)))

  st1 = st%st_start
  st2 = st%st_end
  k1 = st%d%kpt%start
  k2 = st%d%kpt%end   
  SAFE_ALLOCATE(mask%k(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3),1:st%d%dim,st1:st2,k1:k2))
  mask%k = M_z0


  ! generate the map between mesh and square mesh
  call  PES_mask_generate_Lxyz_inv(mask)
  call  PES_mask_generate_Lk(mask)



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Mask Function options
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SAFE_ALLOCATE(mask%mask_R(2))

  !%Variable PESMaskShape
  !%Type integer
  !%Default m_sin2
  !%Section Time-Dependent::PES
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
  !%Section Time-Dependent::PES
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


  if (mask%mode .eq. MODE_PSF ) then 
    width = mask%mask_R(2)-mask%mask_R(1)
    call tdpsf_init(mask%psf,mask%fft, mesh, max_iter,dt,width)
  else
    call PES_mask_generate_mask(mask,mesh)
  end if

  !%Variable PESMaskFilterCutOff 
  !%Type float
  !%Default -1
  !%Section Time-Dependent::PES
  !%Description
  !% In calculation with <tt>PESMaskMode = fullmask_mode<\tt> and NFFT, spurious frequencies 
  !% may lead to numerical instability of the algorithm. This option gives the possibility 
  !% to filter out the unwanted components by setting an energy cut-off. 
  !% If <tt>PESMaskFilterCutOff = -1<\tt> no filter is applied.
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
  !%Section Time-Dependent::PES
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
  !%Section Time-Dependent::PES
  !%Description
  !% The maximum energy for the PES spectrum.
  !%End
  MaxE = maxval(mask%Lk)**2/2
  call parse_float(datasets_check('PESMaskSpectEnergyMax'),&
       units_to_atomic(units_inp%energy,MaxE),mask%energyMax)
  call messages_print_var_value(stdout, "PESMaskSpectEnergyMax",mask%energyMax)



  !%Variable PESMaskSpectEnergyStep 
  !%Type float
  !%Section Time-Dependent::PES
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
  !%Section Time-Dependent::PES
  !%Description
  !% Use interpolation to evaluate the quantities in polar coordinates.
  !% NOTE: In 3D this is practically prohibitive in the present implemetation.
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

  SAFE_ALLOCATE(mask%ext_pot(0:max_iter,1:MAX_DIM))
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
          mask%ext_pot(it,:)= mask%ext_pot(it,:)-field(:) !Sum up all the fields
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
    SAFE_DEALLOCATE_P(mask%M)
    SAFE_DEALLOCATE_P(mask%mask_R)
    SAFE_DEALLOCATE_P(mask%Lxyz_inv)
    SAFE_DEALLOCATE_P(mask%Lk)
    
    if(mask%pw_map_how .ne. PW_MAP_NFFT)then
      call fft_end(mask%fft)
    else
#if defined(HAVE_NFFT) 
      call nfft_end(mask%nfft)
#endif
    end if 

    if(mask%mode == MODE_PSF) then 
      call tdpsf_end(mask%psf)
    end if
   
    if ( mask%filter_k ) then
      SAFE_DEALLOCATE_P(mask%Mk)
    end if

   call cube_end(mask%cube)   
   call dcube_function_free_RS(mask%cube, mask%cM)

  POP_SUB(PES_mask_end)
end subroutine PES_mask_end

! ---------------------------------------------------------
subroutine PES_mask_generate_Lxyz_inv(mask)
  type(PES_mask_t), intent(inout) :: mask
  
  integer :: ix,iy,iz,ip,rankmin,dir
  FLOAT :: dmin! , ixx(MAX_DIM)
  integer :: ixx(MAX_DIM)

  PUSH_SUB(PES_mask_generate_Lxyz_inv)


  mask%Lxyz_inv = -1


  do ix = 1, mask%ll(1)
    ixx(1)= ix - int(mask%ll(1)/2) -1  
    do iy= 1, mask%ll(2) 
      ixx(2)= iy - int(mask%ll(2)/2) -1  
      do iz = 1, mask%ll(3)
        ixx(3)= iz - int(mask%ll(3)/2) -1 
        
         
        !!Support multiresolution
        if( all(ixx(1:mask%mesh%sb%dim) >  mask%mesh%idx%nr(1,1:mask%mesh%sb%dim) &
                              + mask%mesh%idx%enlarge(1:mask%mesh%sb%dim) ) .and. &
            all(ixx(1:mask%mesh%sb%dim) <  mask%mesh%idx%nr(2,1:mask%mesh%sb%dim) &
                               - mask%mesh%idx%enlarge(1:mask%mesh%sb%dim) ) ) then

          ip = mask%mesh%idx%Lxyz_inv(ixx(1),ixx(2),ixx(3)) 
          if (ip > 0 .and. ip < mask%mesh%np_global) then
            mask%Lxyz_inv(ix,iy,iz) = mask%mesh%idx%Lxyz_inv(ixx(1),ixx(2),ixx(3))
          end if 
        end if
      
      end do
    end do
  end do


  POP_SUB(PES_mask_generate_Lxyz_inv)
end subroutine PES_mask_generate_Lxyz_inv

! --------------------------------------------------------
subroutine PES_mask_generate_Lk(mask)
  type(pes_mask_t), intent(inout) :: mask
  
  integer :: ii,nn
  FLOAT   :: temp

  PUSH_SUB(PES_mask_generate_Lk)

  temp = M_TWO * M_PI / (mask%ll(1) * mask%spacing(1))
  nn = mask%ll(1)

  do ii = 1, mask%ll(1)

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
! = Generate the momentum-space filter =
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
!  Generate the mask function on the cubic mesh containing 
!  the simulation box
! ---------------------------------------------------------
subroutine PES_mask_generate_mask(mask,mesh)
  type(PES_mask_t), intent(inout) :: mask
  type(mesh_t),     intent(in)    :: mesh



  PUSH_SUB(PES_mask_generate_mask)

  call PES_mask_generate_mask_function(mask,mesh, mask%shape, mask%mask_R, mask%M)

  POP_SUB(PES_mask_generate_mask)

end subroutine PES_mask_generate_mask

! --------------------------------------------------------
!  Generate the mask function on the cubic mesh containing 
!  the simulation box
! ---------------------------------------------------------
subroutine PES_mask_generate_mask_function(mask,mesh, shape, R, mask_sq, mask_m)
  type(PES_mask_t),     intent(inout)    :: mask
  type(mesh_t),     intent(in)    :: mesh
  integer,          intent(in)    :: shape
  FLOAT,            intent(in)    :: R(2)
  FLOAT,            intent(out)   :: mask_sq(:,:,:)
  FLOAT, optional,  intent(out)   :: mask_m(:)

  integer :: ip, ix3(MAX_DIM)
  integer :: ip_local
  FLOAT   :: dd1,dd2,width
  FLOAT   :: xx(1:MAX_DIM), rr, dd, radius
  integer :: ix,iy,iz, ii
  FLOAT,allocatable :: mask_fn(:)
  logical :: local_

  PUSH_SUB(PES_mask_generate_mask_function)


  ! generate the mask function on the mesh 
  SAFE_ALLOCATE(mask_fn(1:mask%np))

  mask_fn = M_ZERO
  width = R(2) - R(1)
  xx = M_ZERO
 
  !We want the mask cube function to be divided on the nodes?
  local_=(mask%np .eq. mask%mesh%np) .or. (mask%np .eq. mask%mesh%np_part)

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
 
  call dmesh_to_cube(mask%mesh, mask_fn, mask%cube, mask%cM)

  mask_sq = mask%cM%dRS


  if(present(mask_m)) then 
    mask_m = mask_fn
  end if 

  SAFE_DEALLOCATE_A(mask_fn)

  POP_SUB(PES_mask_generate_mask_function)

end subroutine PES_mask_generate_mask_function


! --------------------------------------------------------
subroutine PES_mask_apply_mask(mask,st,mesh)
  type(states_t),   intent(inout) :: st
  type(PES_mask_t), intent(in)    :: mask
  type(mesh_t),     intent(in)    :: mesh

  integer :: ik, ist, idim
  FLOAT, allocatable :: mmask(:)

  PUSH_SUB(PES_mask_apply_mask)
  SAFE_ALLOCATE(mmask(1:mask%mesh%np_part))

  call dcube_to_mesh(mask%cube, mask%cM, mask%mesh, mmask, local = .true.)

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
!  Propagate in time a wavefunction in momentum space with 
!  the Volkov Hamiltonian 
!     
!     wf(p,t+dt)=exp(-i*Hv*dt)*wf(p,t)
!
!  with 
!
!     Hv=(p^2-A)^2/2
!
! NOTE: velocity gauge is implied 
! ---------------------------------------------------------
subroutine PES_mask_Volkov_time_evolution_wf(mask, mesh, dt, iter, wf)
  type(PES_mask_t), intent(in)    :: mask
  type(mesh_t),     intent(in)    :: mesh
  FLOAT,            intent(in)    :: dt
  integer,          intent(in)    :: iter
  CMPLX,            intent(inout) :: wf(:,:,:)

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  FLOAT :: temp(MAX_DIM), vec
  FLOAT :: dd,KK(MAX_DIM)
  integer :: il,ll(MAX_DIM)

  PUSH_SUB(PES_mask_Volkov_time_evolution_wf)


  ! propagate wavefunction in momentum space in presence of a td field (in the velocity gauge)
  do ix = 1, mask%ll(1)
    KK(1) = mask%Lk(ix)
    do iy = 1, mask%ll(2)
      KK(2) = mask%Lk(iy)
      do iz = 1, mask%ll(3)
        KK(3) = mask%Lk(iz)

        vec = sum(( KK(1:mesh%sb%dim) - mask%ext_pot(iter,1:mesh%sb%dim)/P_C)**2) / M_TWO
        wf(ix, iy, iz) = wf(ix, iy, iz) * exp(-M_zI * dt * vec)
        
      end do
    end do
  end do


  POP_SUB(PES_mask_Volkov_time_evolution_wf)
end subroutine PES_mask_Volkov_time_evolution_wf


! ---------------------------------------------------------
subroutine PES_mask_backaction_wf_apply(mask, mesh, wfB,state)
  type(PES_mask_t), intent(in)    :: mask
  type(mesh_t),     intent(in)    :: mesh
  CMPLX,            intent(inout) :: state(:)
  CMPLX,            intent(in)    :: wfB(:,:,:)

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  type(cube_function_t) :: cf
  CMPLX, allocatable :: mf(:)
  FLOAT :: temp(MAX_DIM), vec
  FLOAT :: dd
  integer :: il,ll(MAX_DIM)

  integer:: ip_local

  type(profile_t), save :: prof
  call profiling_in(prof, "PESMASK_back_action")


  PUSH_SUB(PES_mask_backaction_wf_apply)


  ll(1:MAX_DIM) = mask%ll(1:MAX_DIM)
  call cube_function_null(cf)    
  call zcube_function_alloc_RS(mask%cube, cf) 
  SAFE_ALLOCATE(mf(1:mask%mesh%np_part))

  
  cf%zRS = M_z0
  call PES_mask_K_to_X(mask,mesh, wfB , cf%zRs)
  
  call zcube_to_mesh(mask%cube, cf, mask%mesh, mf, local = .true.)  
  state = state + mf    



 
  SAFE_DEALLOCATE_A(mf)
  call zcube_function_free_RS(mask%cube,cf)
 
  POP_SUB(PES_mask_backaction_wf_apply)
 
  call profiling_out(prof)
  
end subroutine PES_mask_backaction_wf_apply







!---------------------------------------------------------
! Local frontend to similar functions defined in cube_function_inc.F90
!---------------------------------------------------------

subroutine PES_mask_mesh_to_cube(mask, mf, cf, local)
  type(pes_mask_t),      intent(in)    :: mask
  CMPLX,  target,        intent(in)    :: mf(:) 
  type(cube_function_t), intent(out) :: cf
  logical, optional,     intent(in)    :: local 

  logical :: local_

  PUSH_SUB(PES_mask_mesh_to_cube)

  local_ = optional_default(local, .true.) .and. mask%mesh%parallel_in_domains
  
  call zmesh_to_cube(mask%mesh, mf, mask%cube, cf, local_)

  POP_SUB(PES_mask_mesh_to_cube)

end subroutine PES_mask_mesh_to_cube


subroutine PES_mask_cube_to_mesh(mask, cf, mf, local)
  type(pes_mask_t),      intent(in)    :: mask
  CMPLX,  target,        intent(out)   :: mf(:) 
  type(cube_function_t), intent(in)    :: cf
  logical, optional,     intent(in)    :: local 

  logical :: local_

  PUSH_SUB(PES_mask_cube_to_mesh)

  local_ = optional_default(local, .true.) .and. mask%mesh%parallel_in_domains

  call zcube_to_mesh(mask%cube, cf, mask%mesh, mf, local_)


  POP_SUB(PES_mask_cube_to_mesh)

end subroutine PES_mask_cube_to_mesh


!---------------------------------------------------------
subroutine fft_X_to_K(mask, mesh, wfin, wfout)
  type(PES_mask_t), intent(in)  :: mask
  type(mesh_t),     intent(in)  :: mesh
  CMPLX,            intent(in)  :: wfin(:,:,:)
  CMPLX,            intent(out) :: wfout(:,:,:)

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  FLOAT   :: temp(MAX_DIM), vec
  FLOAT   :: dd
  integer :: il,ll(MAX_DIM)
  CMPLX, allocatable :: wftmp(:,:,:),wfPlus(:,:,:),wfMinus(:,:,:) 
  integer :: ip_local


  PUSH_SUB(fft_X_to_K)


  ll(1:MAX_DIM) = mask%ll(1:MAX_DIM)
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
subroutine fft_K_to_X(mask,mesh,wfin,wfout,inout)
  type(PES_mask_t), intent(in)  :: mask
  type(mesh_t),     intent(in)  :: mesh
  CMPLX,            intent(in)  :: wfin(:,:,:)
  CMPLX,            intent(out) :: wfout(:,:,:)
  integer,          intent(in)  :: inout


  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  FLOAT   :: temp(MAX_DIM), vec
  FLOAT   :: dd
  integer :: il,ll(MAX_DIM)
  CMPLX, allocatable :: wftmp(:,:,:),wfPlus(:,:,:),wfMinus(:,:,:) 
  integer :: ip_local


  PUSH_SUB(fft_K_to_X)


  ll(1:MAX_DIM) = mask%ll(1:MAX_DIM)
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

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  FLOAT   :: temp(MAX_DIM), vec
  FLOAT   :: k_dot_r
  integer :: il,ll(MAX_DIM)
  integer :: ip_local,kx,ky,kz,ikk(MAX_DIM)


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

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  FLOAT   :: temp(MAX_DIM), vec
  FLOAT   :: dd,k_dot_r
  integer :: il,ll(MAX_DIM)
  integer :: ip_local,kx,ky,kz,ikk(MAX_DIM)


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
! Project the wavefunction on plane waves
!---------------------------------------------------------
subroutine PES_mask_X_to_K(mask,mesh,wfin,wfout)
  type(PES_mask_t), intent(in)  :: mask
  type(mesh_t),     intent(in)  :: mesh
  type(cube_function_t), intent(in):: wfin
  type(cube_function_t), intent(out):: wfout

  FLOAT :: Norm1,Norm2
  type(profile_t), save :: prof
 
  integer ::i


  call profiling_in(prof, "PESMASK_X_to_K")


  PUSH_SUB(PES_mask_X_to_K)

  wfout%zRs=M_z0


  select case(mask%pw_map_how)
    case(PW_MAP_INTEGRAL)
      call integral_X_to_K(mask,mesh,wfin%zRs,wfout%zRs)
      
    case(PW_MAP_FFT)
      call fft_X_to_K(mask,mesh,wfin%zRs,wfout%zRs)
      
    case(PW_MAP_BARE_FFT)
      call zfft_forward(mask%fft, wfin%zRs,wfout%zRs)
      
    case(PW_MAP_TDPSF)
      call tdpsf_X_to_K(mask%psf, wfin%zRs,wfout%zRs)
      
#if defined(HAVE_NFFT) 
    case(PW_MAP_NFFT)
      call znfft_forward(mask%nfft,wfin%zRs,wfout%zRs)
#endif


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

  CMPLX :: DK(MAX_DIM)

  type(profile_t), save :: prof

  call profiling_in(prof, "PESMASK_K_toX")

  PUSH_SUB(PES_mask_K_to_X)
 
  wfout=M_z0

  select case(mask%pw_map_how)
    case(PW_MAP_INTEGRAL)
      call integral_K_to_X(mask,mesh,wfin,wfout)

    case(PW_MAP_FFT)
      call fft_K_to_X(mask,mesh,wfin,wfout,OUT)
      
    case(PW_MAP_BARE_FFT)
      call zfft_backward(mask%fft, wfin,wfout)
!!      DK(:) = M_TWO * M_PI / (mask%ll(:) * mask%spacing(:))

!      wfout = wfout *sqrt(mask%ll(1)*mask%ll(2)*mask%ll(3)/M_PI**(mesh%sb%dim))

!!      wfout = wfout *DK(1)*DK(2)*DK(3)*sqrt(mask%ll(1)*mask%ll(2)*mask%ll(3)/M_PI**(mesh%sb%dim))

!      wfout = wfout *DK(1)*DK(2)*DK(3)/sqrt(M_PI**(mesh%sb%dim))

    case(PW_MAP_TDPSF)
      call tdpsf_K_to_X(mask%psf, wfin,wfout)

#if defined(HAVE_NFFT) 
    case(PW_MAP_NFFT)
      call znfft_backward(mask%nfft,wfin,wfout)
      wfout=wfout/(mask%ll(1)*mask%ll(2)*mask%ll(3)* M_TWO**(mask%enlarge_nfft*mask%mesh%sb%dim))
#endif
      
    case default

  end select

 
  POP_SUB(PES_mask_K_to_X)

  call profiling_out(prof)

end subroutine PES_mask_K_to_X


!---------------------------------------------------------
!
!            Performs all the dirty work 
!
!---------------------------------------------------------
subroutine PES_mask_calc(mask, mesh, st, dt, mask_fn,hm,geo,iter)
  type(PES_mask_t),    intent(inout) :: mask
  type(mesh_t),        intent(in)    :: mesh
  type(states_t),      intent(inout) :: st
  FLOAT,               intent(in)    :: dt
  FLOAT,               intent(in)    :: mask_fn(:) !namely hm%ab_pot
  integer,             intent(in)    :: iter
  type(hamiltonian_t), intent(in)    :: hm
  type(geometry_t),    intent(in)    :: geo

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  type(cube_function_t):: cf1,cf2,cf3,cf4
  FLOAT :: temp(MAX_DIM), vec
  FLOAT :: dd
  integer :: il,i

  FLOAT :: dmin
  integer :: rankmin,ip_local,size

#if defined(HAVE_MPI)
  integer :: status(MPI_STATUS_SIZE)
  integer :: iproc, dataSize
#endif

  type(profile_t), save :: prof

  call profiling_in(prof, "PESMASK_calc")


  PUSH_SUB(PES_mask_calc)



  call cube_function_null(cf1)    
  call zcube_function_alloc_RS(mask%cube, cf1) 
  call cube_function_null(cf2)    
  call zcube_function_alloc_RS(mask%cube, cf2) 
  select case(mask%mode) 
  case(MODE_MASK)
    if(mask%back_action .eqv. .true.) then
      call cube_function_null(cf3)    
      call zcube_function_alloc_RS(mask%cube, cf3) 
    end if
  case(MODE_PSF)
    call cube_function_null(cf3)    
    call zcube_function_alloc_RS(mask%cube, cf3) 
    call cube_function_null(cf4)    
    call zcube_function_alloc_RS(mask%cube, cf4) 
  end select
  

  size = (mask%ll(1))*(mask%ll(2))*(mask%ll(3)) 


  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        cf1%zRs = M_z0
        cf2%zRS = M_z0
       
        call zmesh_to_cube(mask%mesh, st%zpsi(:, idim, ist, ik), mask%cube, cf1, local=.true.)

        select case(mask%mode)
     !----------------------------------------- 
     ! Mask Method
     !----------------------------------------
          case(MODE_MASK)
           
            
            cf1%zRs = (M_ONE-mask%M)*cf1%zRs                                        ! cf1 =(1-M)*U(t2,t1)*\Psi_A(x,t1)
            call PES_mask_X_to_K(mask,mesh,cf1,cf2)                                 ! cf2 = \tilde{\Psi}_A(k,t2)

            if ( mask%filter_k ) then ! apply a filter to the Fourier transform to remove unwanted energies
              
              ASSERT(associated(mask%Mk))

              cf2%zRs= cf2%zRs * mask%Mk(:,:,:) 
            end if
            
            
            if(mask%back_action .eqv. .true.) then
              cf3%zRs = mask%k(:,:,:,idim,ist,ik)                                   ! cf3 = \tilde{\Psi}_B(k,t1) 
              call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter-1,cf3%zRs)  ! cf3 = \tilde{\Psi}_B(k,t2))
              call PES_mask_backaction_wf_apply(mask, mesh,cf3%zRs,st%zpsi(:, idim, ist, ik))
            end if
            
            
            
            cf1%zRs = mask%k(:,:,:, idim, ist, ik)                                  ! cf1 = \Psi_B(k,t1)
            mask%k(:,:,:, idim, ist, ik) =  cf2%zRs                                 ! mask%k = \tilde{\Psi}_A(k,t2)          
            call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter-1,cf1%zRs)    ! cf1 = \tilde{\Psi}_B(k,t2)
            
            
            
            if(mask%back_action .eqv. .true.) then
              ! Apply correction to \Psi_B that enforces it to have 0 components in  deep region A
              !-----                                  ----!
              call PES_mask_K_to_X(mask,mesh,cf1%zRs,cf2%zRs)                       ! cf2 = \Psi_B(x,t2)
              cf2%zRs= (mask%M)*cf2%zRs                                             ! cf2 = M*\Psi_B(x,t1)    
              call PES_mask_X_to_K(mask,mesh,cf2,cf3)                           
              cf1%zRs=cf1%zRs-cf3%zRs
              !-----                                  ----!
            end if
            
            
            cf2%zRs=cf1%zRs                                                         ! cf2 = \tilde{\Psi}_B(k,t2)
            
            
            ! and add to our spectrum
            mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) + cf2%zRs 


     !----------------------------------------- 
     ! Passive Mask method
     !----------------------------------------
          case(MODE_PASSIVE)
            

            cf1%zRs = (M_ONE-mask%M)*cf1%zRs
            call PES_mask_X_to_K(mask,mesh,cf1,cf2) 
            
            mask%k(:,:,:, idim, ist, ik) = cf2%zRs(:,:,:)


     !----------------------------------------- 
     ! Phase Space Filter
     !----------------------------------------
     
          case(MODE_PSF)
! !          if (MOD(iter*dt,mask%psf%Tstep) .eq. M_ZERO) then 
! !            write (*,*) "APPLY TDPSF!!"
! !!            call tdpsf_filter_out(mask%psf,wf1,wf2)
! !            call PES_mask_X_to_K(mask,mesh,wf1,wf3)
! !            call PES_mask_K_to_X(mask,mesh,wf3,wf2)
! 
!             call tdpsf_X_to_K(mask%psf,wf1,wf3)
!             call tdpsf_K_to_X(mask%psf,wf3,wf2)
! 
!             wf2 = wf1 - wf2
!             
!             
!             if(mask%back_action) then
!               wf1 = mask%k(:,:,:, idim, ist, ik) - wf4
!               call PES_mask_K_to_X(mask,mesh,wf1,mask%k(:,:,:, idim, ist, ik))
!               wf2 = wf2 + mask%k(:,:,:, idim, ist, ik)
!             end if
! 
!             !substitute the KS wf with the filtered one
!             call PES_mask_square_to_mesh(mask,mesh,st%zpsi(:, idim, ist, ik),wf2,MaskHow = 3,Const = M_ONE)
! 
!             !the out-going part of the wf
!             wf4 = wf1 - wf2
!             call PES_mask_X_to_K(mask,mesh,wf4,wf3) 
! 
! !            do ix=1,mask%ll(1)
! !              write (*,*) ix, wf1(ix,1,1),wf2(ix,1,1),wf3(ix,1,1)
! !            end do             
! !!            call PES_mask_X_to_K(mask,mesh,wf3,wf1)
! !            call PES_mask_K_to_X(mask,mesh,mask%k(:,:,:, idim, ist, ik),wf3)
! !            call PES_mask_X_to_K(mask,mesh,wf3,mask%k(:,:,:, idim, ist, ik))
! 
!             call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter,mask%k(:,:,:, idim, ist, ik) )
!             wf4 = mask%k(:,:,:, idim, ist, ik)            
! 
! !            mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) + wf1(:,:,:)
!             mask%k(:,:,:, idim, ist, ik) =  wf4 + wf3
!             
! !          end if 

            call tdpsf_X_to_K(mask%psf,cf1%zRs,cf3%zRs)
            call tdpsf_K_to_X(mask%psf,cf3%zRs,cf2%zRs)

            cf2%zRs = cf1%zRs - cf2%zRs

            if(mask%back_action) then
              cf1%zRs = mask%k(:,:,:, idim, ist, ik) - cf4%zRs
              call PES_mask_K_to_X(mask,mesh,cf1%zRs,mask%k(:,:,:, idim, ist, ik))
              cf2%zRs = cf2%zRs + mask%k(:,:,:, idim, ist, ik)
            end if

            !substitute the KS wf with the filtered one
            call zcube_to_mesh(mask%cube, cf2, mask%mesh, st%zpsi(:, idim, ist, ik), local = .true.)
            !the out-going part of the wf
            cf4%zRs = cf1%zRs - cf2%zRs
            call PES_mask_X_to_K(mask,mesh,cf4,cf3) 

            call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter,mask%k(:,:,:, idim, ist, ik) )
            cf4%zRs = mask%k(:,:,:, idim, ist, ik)            

            mask%k(:,:,:, idim, ist, ik) =  cf4%zRs + cf3%zRs


          case default
           !Program should die before coming here
            write(message(1),'(a)') "PhotoElectroSpectrum = pes_mask. Unrecognized calculation mode." 
            call messages_fatal(1)
 
       end select


      end do
    end do
  end do


  if(mask%mode .eq. MODE_MASK ) call PES_mask_apply_mask(mask,st,mesh)  !apply the mask to all the KS orbitals


  call zcube_function_free_RS(mask%cube, cf1)
  call zcube_function_free_RS(mask%cube, cf2)
  select case(mask%mode) 
  case(MODE_MASK)
    if(mask%back_action) then
      call zcube_function_free_RS(mask%cube, cf3)
    end if
  case(MODE_PSF)
    call zcube_function_free_RS(mask%cube, cf3)
    call zcube_function_free_RS(mask%cube, cf4)
  end select





  POP_SUB(PES_mask_calc)

  call profiling_out(prof)
end subroutine PES_mask_calc


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
