!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
  type(PES_mask_t),    intent(out) :: mask
  type(mesh_t),target, intent(in)  :: mesh
  type(simul_box_t),   intent(in)  :: sb
  type(states_t),      intent(in)  :: st
  type(hamiltonian_t), intent(in)  :: hm
  integer,             intent(in)  :: max_iter
  FLOAT,               intent(in)  :: dt

  type(block_t) :: blk

  integer :: ll(MAX_DIM), il,it,ii
  FLOAT :: field(MAX_DIM)
  FLOAT :: DeltaE,MaxE,MaxDR
  FLOAT :: width 
  integer :: dir
  FLOAT, allocatable  :: X1(:),X2(:),X3(:)  

  PUSH_SUB(PES_mask_init)

 
  mask%mesh => mesh

  call messages_experimental('Photo-electron spectrum')  

  write(message(1),'(a,i1,a)') 'Info: Calculating PES using mask technique.'
  call messages_info(1)




  if (st%parallel_in_states) then
    call messages_experimental("PES_mask parallelization on states")

    if(mesh%parallel_in_domains) then
      write(message(1),'(a)') "PES_mask: simultaneous parallelization on mesh and states not supported"
      write(message(2),'(a)') "Modify ParallelizationStrategy and rerun." 
      call messages_fatal(2) 
      end if    
  endif



  !%Variable PESMaskOutWaveProjection
  !%Type integer
  !%Default fft_bare
  !%Section Time-Dependent::PES
  !%Description
  !% How to calculate the plane waves projection.
  !%Option integral 1
  !% Direct integration_map.
  !%Option fft_map 2 
  !%1D only.  
  !%Option fft_bare 3 
  !%Bare FFT map. This will contain also ingoing terms.
  !%Option tdpsf_map 4
  !%Time-dependent phase-space filter map.
  !%Option nfft_map 5
  !%Non-equispaced FFT map.
  !%End
  call parse_integer(datasets_check('PESMaskOutWaveProjection'),PW_MAP_BARE_FFT,mask%pw_map_how)

  if(.not.varinfo_valid_option('PESMaskOutWaveProjection', mask%pw_map_how)) call input_error('PESMaskOutWaveProjection')
  call messages_print_var_option(stdout, "PESMaskOutWaveProjection", mask%pw_map_how)

#if !defined(HAVE_NFFT) 
  if (mask%pw_map_how ==  PW_MAP_NFFT) then
    message(1) = "PESMaskOutWaveProjection = nfft_map requires libnfft3. Recompile and try again." 
    call messages_fatal(1) 
  endif 
#endif

  
  !%Variable PESMaskEnlargeLev
  !%Type integer
  !%Default 0
  !%Section Time-Dependent::PES
  !%Description
  !% Mask box enlargement level. Enlarges the mask bounding box by a factor 2**<tt>PESMaskEnlargeLev</tt>
  !% in order to avoid wrapping at the boundaries.
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
  !% Note: the corresponding Fourier space is shrunk by the same factor.
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
    call fft_init(mask%ll,sb%dim,fft_complex,mask%fft, optimize = .not.simul_box_is_periodic(sb))

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


  ll(1:MAX_DIM) = mask%ll(1:MAX_DIM)
  
  SAFE_ALLOCATE(mask%k(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3),1:st%d%dim,1:st%nst,1:st%d%nik))
  SAFE_ALLOCATE(mask%Lxyz_inv(1:ll(1),1:ll(2),1:ll(3)))
  
  SAFE_ALLOCATE(mask%M(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3)))
  SAFE_ALLOCATE(mask%Lk(1:mask%ll(1)))

  mask%k = M_z0


  ! generate the map between mesh and square mesh
  call  PES_mask_generate_Lxyz_inv(mask)
  !
  call  PES_mask_generate_Lk(mask)


  !%Variable PESMaskMode
  !%Type integer
  !%Default mask_mode
  !%Section Time-Dependent::PES
  !%Description
  !% PES calculation mode.
  !%Option mask_mode 2
  !% Mask function method. 
  !%Option passive_mode 4 
  !% Passive analysis of the wf.
  !%Option psf_mode 8
  !% Phase-space filter.
  !%End
  call parse_integer(datasets_check('PESMaskMode'),MODE_MASK,mask%mode)

  if(.not.varinfo_valid_option('PESMaskMode', mask%mode)) call input_error('PESMaskMode')
  call messages_print_var_option(stdout, "PESMaskMode", mask%mode)

  !%Variable PESMaskEnableBackAction 
  !%Type logical
  !%Default false
  !%Section Time-Dependent::PES
  !%Description
  !% Enable photoelectron scattering-waves correction to the evolution of 
  !% the bound states.  
  !%End
  call parse_logical(datasets_check('PESMaskEnableBackAction'),.false.,mask%back_action)

  if(mask%back_action .eqv. .true.) then
    message(1)= "Input: PES_mask back action correction ENABLED."
    call messages_info(1)

  end if




  SAFE_ALLOCATE(mask%mask_fn(1:mesh%np_global))
  SAFE_ALLOCATE(mask%mask_R(2))

  mask%mask_fn=M_ZERO
  ! these initializations have no effect, they are overwritten by parse_block_float below. DAS
  mask%mask_R(1)=mesh%sb%rsize/M_TWO
  mask%mask_R(2)=mesh%sb%rsize

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
    !call input_error('PESMaskSize')
  else 
    call parse_block_float(blk, 0, 0, mask%mask_R(1), units_inp%length)
    call parse_block_float(blk, 0, 1, mask%mask_R(2), units_inp%length)
  end if
  
  if(mask%mask_R(2) .gt. mesh%sb%rsize)  mask%mask_R(2) = mesh%sb%rsize 


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
  !%error function.
  !%End
  call parse_integer(datasets_check('PESMaskShape'),M_SIN2,mask%shape)

  if(.not.varinfo_valid_option('PESMaskShape', mask%shape)) call input_error('PESMaskShape')
  call messages_print_var_option(stdout, "PESMaskShape", mask%shape)

  write(message(1),'(a,es10.3,3a)') & 
           "Input: Mask  R1 = ", units_from_atomic(units_inp%length, mask%mask_R(1) )," [a.u.]"
  write(message(2),'(a,es10.3,3a)') & 
           "             R2 = ", units_from_atomic(units_inp%length, mask%mask_R(2) )," [a.u.]"
  call messages_info(2)


  if (mask%mode .eq. MODE_PSF ) then 
    width = mask%mask_R(2)-mask%mask_R(1)
    call tdpsf_init(mask%psf,mask%fft, mesh, max_iter,dt,width)
  else
    call PES_mask_generate_mask(mask,mesh)
  end if



  !%Variable PESMaskPropagator 
  !%Type integer
  !%Default volkov
  !%Section Time-Dependent::PES
  !%Description
  !% Photoelectron waves propagator in momentum space.
  !%Option volkov 2
  !% Plane wave evolves with exp(i(p-A(t)/c)^2*dt/2).
  !%Option free 1
  !% Free plane-wave propagation.   
  !%End
  call parse_integer(datasets_check('PESMaskPropagator'),VOLKOV,mask%sw_evolve)

  if(.not.varinfo_valid_option('PESMaskPropagator',mask%sw_evolve)) call input_error('PESMaskPropagator')
  call messages_print_var_option(stdout, "PESMaskPropagator",mask%sw_evolve)




  SAFE_ALLOCATE(mask%ext_pot(0:max_iter,1:MAX_DIM))
  mask%ext_pot=M_ZERO

  if(mask%sw_evolve .eq. VOLKOV) then
  ! Precalculate the potential vector for all the simulation time
    do it = 1, max_iter
      do il = 1, hm%ep%no_lasers
        field=M_ZERO
        select case(laser_kind(hm%ep%lasers(il)))
        case(E_FIELD_MAGNETIC, E_FIELD_ELECTRIC)
          write(message(1),'(a)') 'PESMask works only with vector_potential unless in passive_mode.'
          call messages_warning(1)
        case(E_FIELD_VECTOR_POTENTIAL)
          call laser_field(hm%ep%lasers(il), field, it*dt)
       end select
       mask%ext_pot(it,:)= mask%ext_pot(it,:)-field(:) !Sum up all the fields
     end do
    end do
  else 
  end if


  !%Variable PESMaskIncludePsiA
  !%Type logical
  !%Default false
  !%Section Time-Dependent::PES
  !%Description
  !% Add the contribution of Psi_A to the photo-electron spectrum.
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



  POP_SUB(PES_mask_init)
end subroutine PES_mask_init


! ---------------------------------------------------------
subroutine PES_mask_end(mask)
  type(PES_mask_t), intent(inout) :: mask

  PUSH_SUB(PES_mask_end)

    SAFE_DEALLOCATE_P(mask%k)

    SAFE_DEALLOCATE_P(mask%ext_pot)
    SAFE_DEALLOCATE_P(mask%M)
    SAFE_DEALLOCATE_P(mask%mask_fn)
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


! --------------------------------------------------------
!  Generate the mask function on the cubic mesh containing 
!  the simulation box
! ---------------------------------------------------------
subroutine PES_mask_generate_mask(mask,mesh)
  type(PES_mask_t), intent(inout) :: mask
  type(mesh_t),     intent(in)    :: mesh

  integer :: ip, ix3(MAX_DIM)
  integer :: ip_local
  FLOAT   :: dd1,dd2,width
  FLOAT   :: xx(MAX_DIM), rr, dd, radius
  integer :: ix,iy,iz


  PUSH_SUB(PES_mask_generate_mask)

  ! generate the mask function on the mesh 
  mask%mask_fn = M_ZERO
  width = mask%mask_R(2)-mask%mask_R(1)

  select case(mask%shape)
    case(M_SIN2)
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, coords=xx)
        dd = rr -  mask%mask_R(1) 
        if(dd .gt. M_ZERO ) then 
          if (dd .lt. width) then
            mask%mask_fn(ip) = M_ONE * sin(dd * M_PI / (M_TWO * (width) ))**2
          else 
            mask%mask_fn(ip) = M_ONE 
          end if
        end if
      end do
      
    case(M_STEP)
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, coords=xx)
        dd = rr -  mask%mask_R(1) 
        if(dd .gt. M_ZERO ) then 
          if (dd .lt. width) then
            mask%mask_fn(ip) = M_ONE 
          else 
            mask%mask_fn(ip) = M_ZERO
          end if
        end if
      end do
      
    case(M_ERF)

    case default
    !Program should quit before coming here
  end select



  !the mask is zero in the points of the cube not contained in the
  !simulation box
  mask%M = M_z0


  do ix=1,mask%ll(1)
    do iy=1,mask%ll(2)
      do iz=1,mask%ll(3)
        
        ip= mask%Lxyz_inv(ix,iy,iz)
        if (ip > 0) then
          mask%M(ix,iy,iz) = (1-mask%mask_fn(ip))
        end if
      end do
    end do
  end do

  POP_SUB(PES_mask_generate_mask)

end subroutine PES_mask_generate_mask


! --------------------------------------------------------
subroutine PES_mask_apply_mask(mask,st,mesh)
  type(states_t),   intent(inout) :: st
  type(PES_mask_t), intent(in)    :: mask
  type(mesh_t),     intent(in)    :: mesh

  integer :: ik, ist, idim

  PUSH_SUB(PES_mask_apply_mask)

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim
        st%zpsi(1:mesh%np, idim, ist, ik) = st%zpsi(1:mesh%np, idim, ist, ik)*(M_ONE - mask%mask_fn(1:mesh%np))
      end do
    end do
  end do

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
  CMPLX, allocatable ::  wf2(:,:,:)
  FLOAT :: temp(MAX_DIM), vec
  FLOAT :: dd
  integer :: il,ll(MAX_DIM)

  integer:: ip_local

  type(profile_t), save :: prof
  call profiling_in(prof, "PESMASK_back_action")


  PUSH_SUB(PES_mask_backaction_wf_apply)


  ll(1:MAX_DIM) = mask%ll(1:MAX_DIM)
  SAFE_ALLOCATE(wf2(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
 

  wf2 = M_z0
  call PES_mask_K_to_X(mask,mesh, wfB , wf2)
  
  
  do ix=1,mask%ll(1)
    do iy=1,mask%ll(2)
      do iz=1,mask%ll(3)
        
        ip= mask%Lxyz_inv(ix,iy,iz)
        if (ip > 0) then
          state(ip) = state(ip) + wf2(ix, iy, iz)
        end if
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(wf2)
 

  POP_SUB(PES_mask_backaction_wf_apply)
 
  call profiling_out(prof)
  
end subroutine PES_mask_backaction_wf_apply





!---------------------------------------------------------------
! Maps a wavefunction on the mesh to a function on the smallest
! bounding box containing the mesh. 
! Optionally can also apply a operation on the wavefunction such
! as multiplication for a mask function M of the mesh and a constant 
! factor.
!
!          wf_m(mesh)[*M(mesh)*C] -> wf_sq(square)  
!
! It is also possible to apply the complement mask (1-M) via the 
! optional parameter CompM.
!--------------------------------------------------------------
subroutine PES_mask_mesh_to_square(mask,mesh, wf_m,wf_sq, MaskHow,Const)
  type(PES_mask_t), intent(in)  :: mask
  type(mesh_t),     intent(in)  :: mesh
  CMPLX,            intent(in)  :: wf_m(:)
  CMPLX,            intent(out) :: wf_sq(:,:,:) 
  integer,          intent(in)  :: MaskHow
  FLOAT,            intent(in)  :: Const

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM)
  FLOAT :: dd
  integer :: il,ip_local,ixx,iyy,izz
  INTEGER :: idx(MAX_DIM)

  type(profile_t), save :: prof
  call profiling_in(prof, "PESMASK_mesh_to_square")



  ! no push_sub, called too frequently
  
  wf_sq= M_z0

  do ix=1, mask%ll(1)
    ixx = ix
    do iy=1,mask%ll(2)
      iyy = iy
      do iz=1,mask%ll(3)
        izz = iz
        
        ip= mask%Lxyz_inv(ix,iy,iz)
        if (ip > 0) then
          select case(MaskHow)
            case(1)
              wf_sq(ixx,iyy,izz) = Const*mask%M(ixx,iyy,izz) * wf_m(ip)
            case(2)
              wf_sq(ixx,iyy,izz) = Const*(1-mask%M(ixx,iyy,izz)) * wf_m(ip)
            case(3)
              wf_sq(ixx,iyy,izz) = Const * wf_m(ip)
            case default
              wf_sq(ixx,iyy,izz) = Const * wf_m(ip)
              
          end select
            
        end if
      end do
    end do
  end do

  call profiling_out(prof)

end subroutine PES_mask_mesh_to_square

!---------------------------------------------------------------
! Inverse function of PES_mask_mesh_to_square
!
!          wf_sq(square) -> wf_sq(mesh)*[*M(mesh)*C]  
!
!--------------------------------------------------------------
subroutine PES_mask_square_to_mesh(mask,mesh, wf_m,wf_sq, MaskHow,Const)
  type(PES_mask_t), intent(in)  :: mask
  type(mesh_t),     intent(in)  :: mesh
  CMPLX,            intent(out) :: wf_m(:)
  CMPLX,            intent(in)  :: wf_sq(:,:,:)
  integer,          intent(in)  :: MaskHow
  FLOAT,            intent(in)  :: Const

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM)
  FLOAT :: dd
  integer :: il,ip_local,idx(MAX_DIM),ixx,iyy,izz
  FLOAT :: CConst

  type(profile_t), save :: prof
  call profiling_in(prof, "PESMASK_square_to_mesh")


  ! no push_sub, called too frequently
  
  wf_m = M_z0
  
  do ix=1, mask%ll(1)
    ixx = ix
    do iy=1,mask%ll(2)
      iyy = +iy
      do iz=1,mask%ll(3)
        izz = iz
        
        ip= mask%Lxyz_inv(ix,iy,iz)
        if (ip > 0) then
          select case(MaskHow)
            case(1)
              wf_m(ip) = Const*mask%M(ixx,iyy,izz) * wf_sq(ixx,iyy,izz)
            case(2)
              wf_m(ip) = Const*(1-mask%M(ixx,iyy,izz)) * wf_sq(ixx,iyy,izz)
            case(3)
              wf_m(ip) = Const * wf_sq(ixx,iyy,izz)
            case default
              wf_m(ip) = Const * wf_sq(ixx,iyy,izz)
              
          end select

        end if
      end do
    end do
  end do
  
  call profiling_out(prof)
  
end subroutine PES_mask_square_to_mesh

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
  CMPLX,            intent(in)  :: wfin(:,:,:)
  CMPLX,            intent(out) :: wfout(:,:,:)

  FLOAT :: Norm1,Norm2

  type(profile_t), save :: prof
  call profiling_in(prof, "PESMASK_X_to_K")

  PUSH_SUB(PES_mask_X_to_K)

  wfout=M_z0

  select case(mask%pw_map_how)
    case(PW_MAP_INTEGRAL)
      call integral_X_to_K(mask,mesh,wfin,wfout)
      
    case(PW_MAP_FFT)
      call fft_X_to_K(mask,mesh,wfin,wfout)
      
    case(PW_MAP_BARE_FFT)
      call zfft_forward(mask%fft, wfin,wfout)
      
    case(PW_MAP_TDPSF)
      call tdpsf_X_to_K(mask%psf, wfin,wfout)
      
#if defined(HAVE_NFFT) 
    case(PW_MAP_NFFT)
      call znfft_forward(mask%nfft,wfin,wfout)
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
  CMPLX, allocatable :: wf1(:,:,:), wf2(:,:,:),wf3(:,:,:),wf4(:,:,:)
  FLOAT :: temp(MAX_DIM), vec,aa(1:mesh%np_part,1:mesh%sb%dim)
  FLOAT :: dd
  integer :: il

  FLOAT :: dmin
  integer :: rankmin,ip_local,size

#if defined(HAVE_MPI)
  integer :: status(MPI_STATUS_SIZE)
  integer :: iproc, dataSize
#endif

  type(profile_t), save :: prof

  call profiling_in(prof, "PESMASK_calc")


  PUSH_SUB(PES_mask_calc)


  !Allocate memory for the different schemes
  SAFE_ALLOCATE(wf1(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
  SAFE_ALLOCATE(wf2(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
  select case(mask%mode) 
    case(MODE_MASK)
      if(mask%back_action .eqv. .true.) then
        SAFE_ALLOCATE(wf3(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
      end if
    case(MODE_PSF)
      SAFE_ALLOCATE(wf3(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
      SAFE_ALLOCATE(wf4(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
  end select

  size = (mask%ll(1))*(mask%ll(2))*(mask%ll(3)) 


  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        wf1 = M_z0
        wf2 = M_z0
        
        call PES_mask_mesh_to_square(mask,mesh,st%zpsi(:, idim, ist, ik),wf1,MaskHow = 3,Const = M_ONE)
 
#if defined(HAVE_MPI)

        if(mesh%parallel_in_domains)then

          !send all the wavefunctions to root node
          if(mesh%mpi_grp%rank .gt. 0 ) then
            call MPI_Send(wf1,size,MPI_CMPLX,0,666, mesh%mpi_grp%comm, mpi_err)
          else
            do iproc = 1, mesh%mpi_grp%size-1
              call MPI_Recv(wf2,size,MPI_CMPLX,iproc,666, mesh%mpi_grp%comm, status, mpi_err)
              ! add contribution for other nodes
              wf1 = wf1 +wf2 
            end do
            
          end if
        end if
#endif

        select case(mask%mode)
     !----------------------------------------- 
     ! Mask Method
     !----------------------------------------
          case(MODE_MASK)
           
            wf1 = (M_ONE-mask%M)*wf1                                            ! wf1 =(1-M)*U(t2,t1)*\Psi_A(x,t1)
            call PES_mask_X_to_K(mask,mesh,wf1,wf2)                             ! wf2 = \tilde{\Psi}_A(k,t2)
            
            if(mask%back_action .eqv. .true.) then
              wf3 = mask%k(:,:,:,idim,ist,ik)                                   ! wf3 = \tilde{\Psi}_B(k,t1) 
              call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter-1,wf3)  ! wf3 = \tilde{\Psi}_B(k,t2))
              call PES_mask_backaction_wf_apply(mask, mesh,wf3,st%zpsi(:, idim, ist, ik))
            end if
            
            wf1 = mask%k(:,:,:, idim, ist, ik)                                  ! wf1 = \Psi_B(k,t1)
            mask%k(:,:,:, idim, ist, ik) =  wf2                                 ! mask%k = \tilde{\Psi}_A(k,t2)          
            call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter-1,wf1)    ! wf1 = \tilde{\Psi}_B(k,t2)
            
            if(mask%back_action .eqv. .true.) then
              ! Apply correction to \Psi_B that enforces it to have 0 components in  deep region A
              !-----                                  ----!
              call PES_mask_K_to_X(mask,mesh,wf1,wf2)                           ! wf2 = \Psi_B(x,t2)
              wf2= (mask%M)*wf2                                                 ! wf2 = M*\Psi_B(x,t1)    
              call PES_mask_X_to_K(mask,mesh,wf2,wf3)                           
              wf1=wf1-wf3
              !-----                                  ----!
            end if

            wf2=wf1                                                             ! wf2 = \tilde{\Psi}_B(k,t2)

            ! and add to our spectrum
            mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) + wf2 


     !----------------------------------------- 
     ! Passive Mask method
     !----------------------------------------
          case(MODE_EXACT)
            
            wf1 = (M_ONE-mask%M)*wf1
            call PES_mask_X_to_K(mask,mesh,wf1,wf2) 
            
            mask%k(:,:,:, idim, ist, ik) = wf2(:,:,:)


     !----------------------------------------- 
     ! Phase Space Filter
     !----------------------------------------
     
          case(MODE_PSF)
!          if (MOD(iter*dt,mask%psf%Tstep) .eq. M_ZERO) then 
!            write (*,*) "APPLY TDPSF!!"
!!            call tdpsf_filter_out(mask%psf,wf1,wf2)
!            call PES_mask_X_to_K(mask,mesh,wf1,wf3)
!            call PES_mask_K_to_X(mask,mesh,wf3,wf2)

            call tdpsf_X_to_K(mask%psf,wf1,wf3)
            call tdpsf_K_to_X(mask%psf,wf3,wf2)

            wf2 = wf1 - wf2
            
            
            if(mask%back_action) then
              wf1 = mask%k(:,:,:, idim, ist, ik) - wf4
              call PES_mask_K_to_X(mask,mesh,wf1,mask%k(:,:,:, idim, ist, ik))
              wf2 = wf2 + mask%k(:,:,:, idim, ist, ik)
            end if

            !substitute the KS wf with the filtered one
            call PES_mask_square_to_mesh(mask,mesh,st%zpsi(:, idim, ist, ik),wf2,MaskHow = 3,Const = M_ONE)

            !the out-going part of the wf
            wf4 = wf1 - wf2
            call PES_mask_X_to_K(mask,mesh,wf4,wf3) 

!            do ix=1,mask%ll(1)
!              write (*,*) ix, wf1(ix,1,1),wf2(ix,1,1),wf3(ix,1,1)
!            end do             
!!            call PES_mask_X_to_K(mask,mesh,wf3,wf1)
!            call PES_mask_K_to_X(mask,mesh,mask%k(:,:,:, idim, ist, ik),wf3)
!            call PES_mask_X_to_K(mask,mesh,wf3,mask%k(:,:,:, idim, ist, ik))

            call PES_mask_Volkov_time_evolution_wf(mask, mesh,dt,iter,mask%k(:,:,:, idim, ist, ik) )
            wf4 = mask%k(:,:,:, idim, ist, ik)            

!            mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) + wf1(:,:,:)
            mask%k(:,:,:, idim, ist, ik) =  wf4 + wf3
            
!          end if 


          case default
           !Program should die before coming here
           ! Then a fatal error should be generated. -DAS
        end select


      end do
    end do
  end do
 

  if(mask%mode .eq. MODE_MASK ) call PES_mask_apply_mask(mask,st,mesh)  !apply the mask to all the states


#ifdef HAVE_MPI
  ! wait for all processors to finish!
  ! is this really necessary? I doubt it. --DAS
  if(st%mpi_grp%size .gt. 1 ) then
    call MPI_Barrier(st%mpi_grp%comm, mpi_err)
  end if
#endif



  SAFE_DEALLOCATE_A(wf1)
  SAFE_DEALLOCATE_A(wf2)
  select case(mask%mode) 
  case(MODE_MASK)
    if(mask%back_action) then
      SAFE_DEALLOCATE_A(wf3)
    end if
  case(MODE_PSF)
    SAFE_DEALLOCATE_A(wf3)
    SAFE_DEALLOCATE_A(wf4)
  end select

  POP_SUB(PES_mask_calc)

  call profiling_out(prof)
end subroutine PES_mask_calc


!---------------------------------------------------------------------------
! Collect the states from all the nodes when the code run parallel on states
! --------------------------------------------------------------------------
subroutine PES_mask_collect(mask, st,mesh)
  type(PES_mask_t), intent(inout) :: mask
  type(states_t),   intent(in)    :: st
  type(mesh_t),     intent(in)    :: mesh

  CMPLX, allocatable :: wf(:,:,:)
  FLOAT, allocatable :: wfr(:,:,:)
  integer ::  idim, ist, ik

#ifdef HAVE_MPI
  integer :: iproc, size
  integer :: status(MPI_STATUS_SIZE)
#endif

 PUSH_SUB(PES_mask_collect)

#ifdef HAVE_MPI

  SAFE_ALLOCATE(wf(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
  SAFE_ALLOCATE(wfr(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))

  size = (mask%ll(1))*(mask%ll(2))*(mask%ll(3))

 

  do ik = 1, st%d%nik
    do ist = 1, st%nst
      
      do idim = 1, st%d%dim
        
        !send all the wavefunctions to root node
        if(st%mpi_grp%rank .gt. 0 ) then
          wf= mask%k(:,:,:, idim, ist, ik) 
          call MPI_Send(wf, size, MPI_CMPLX,0, 1, st%mpi_grp%comm, mpi_err)
        else
          !root node collects all the data
          do iproc = 1, st%mpi_grp%size-1 
            call MPI_Recv(wf, size ,MPI_CMPLX, iproc, 1, st%mpi_grp%comm, status, mpi_err)
            mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) + wf
          end do
        end if


      end do
    end do
  end do
  
  !wait for all threads to finish
  ! again, I doubt this is necessary. --DAS
  if(st%mpi_grp%size .gt. 1 ) then
    call MPI_Barrier(st%mpi_grp%comm, mpi_err)
  end if
  
  SAFE_DEALLOCATE_A(wf)
  SAFE_DEALLOCATE_A(wfr)

#endif

  POP_SUB(PES_mask_collect)
end subroutine PES_mask_collect

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
