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
!! $Id: pes.F90 7095 2010-11-25 14:48:16Z umberto $

#include "global.h"

module tdpsf_m
  use datasets_m
  use fft_m
  use global_m
  use io_m
  use io_binary_m
  use mesh_m
  use index_m
  use messages_m
  use parser_m
  use profiling_m
  use simul_box_m
  use states_m
  use unit_m
  use unit_system_m
  use mpi_m
  use hamiltonian_m
  use geometry_m
!  use h_sys_output_m
  use output_m
  use grid_m
  use states_io_m
  use io_function_m
  use density_m
  use batch_m
  use varinfo_m
  use loct_math_m

  implicit none

!  private
  public ::                             &
    tdpsf_init,                         &                
    tdpsf_filter_out

  integer, parameter ::   &
    PLUS           =  1,  &   
    MINUS          =  2 



  type tdpsf_t
    
    type(fft_t)  :: fft
    type(mesh_t), pointer :: mesh      

    FLOAT :: kb 
    FLOAT :: Tstep 
    FLOAT :: kmin 
    FLOAT :: kmax 
    FLOAT :: width
    FLOAT :: sigma
    FLOAT :: delta

    FLOAT, pointer :: XFltr(:,:,:)     !< Filter on real space       
    FLOAT, pointer :: KFltr(:,:,:)     !< Filter on K space
    
    INTEGER        ::ll(MAX_DIM)
      

  end type tdpsf_t


contains

  ! ---------------------------------------------------------
  subroutine tdpsf_init(psf,fft, mesh, max_iter,dt,width)
    type(tdpsf_t),       intent(out)     :: psf
    type(fft_t),         intent(in)      :: fft
    type(mesh_t),target, intent(in)      :: mesh
    integer,             intent(in)      :: max_iter
    FLOAT,               intent(in)      :: dt
    FLOAT,               intent(in)      :: width

    integer :: ll(MAX_DIM),lmax,i,dim
    FLOAT   :: sigma,kmax,kmin,dlmax
    FLOAT   :: sigma_max,sigma_min
    

    PUSH_SUB(tdpsf_init)

    ll(1:MAX_DIM) = mesh%idx%ll(1:MAX_DIM)
    dim = mesh%sb%dim 


    lmax = M_ZERO
    dlmax= M_ZERO
    kmax = M_ZERO
    do i = 1, MAX_DIM
      if (ll(i) .gt. lmax) lmax=ll(i)
      if(mesh%spacing(i) .gt. dlmax) dlmax=mesh%spacing(i)
      if(M_ONE/(mesh%spacing(i)) .gt. kmax) kmax= M_PI/(mesh%spacing(i))
    end do   
 
    SAFE_ALLOCATE(psf%XFltr(1:2,1:MAX_DIM,1:lmax))
    SAFE_ALLOCATE(psf%KFltr(1:2,1:MAX_DIM,1:lmax))
 
    psf%XFltr=M_ZERO
    psf%KFltr=M_ZERO
   
    psf%ll=ll

    psf%mesh => mesh  ! set a pointer to the mesh 
    
    psf%fft = fft
     
    !waves cannot cross buffer region for in a time shorter than Tstep
    psf%Tstep= width/(3*M_PI/dlmax)
    psf%Tstep=2.0E-3    
    
    psf%Tstep= int(psf%Tstep/dt)*dt ! cast Tstep as an integer multiple of dt

    psf%width = width

    write(message(1),*) 'TDPSF:  time step = ',units_from_atomic(units_inp%time,psf%Tstep)
    call messages_info(1)
    write(message(1),*) 'TDPSF:  Buffer region = ',units_from_atomic(units_inp%length,psf%width)
    call messages_info(1)
    write (*,*) "dlmax",dlmax,"lmax",lmax,"kmax",kmax


    !%Variable TDPSFSigma 
    !%Type float
    !%Section Time-Dependent::TDPSF
    !%Description
    !% Standard deviation of the phase space filter.
    !%End
    call parse_float(datasets_check('TDPSFSigma'),sqrt(M_TWO),psf%sigma)
    call messages_print_var_value(stdout, "TDPSFSigma",psf%sigma)


    !%Variable TDPSFDelta 
    !%Type float
    !%Section Time-Dependent::TDPSF
    !%Description
    !% Filter error treshold.
    !%End
    call parse_float(datasets_check('TDPSFDelta'),M_ONE*1.0E-4,psf%delta)
    call messages_print_var_value(stdout, "TDPSFDelta",psf%delta)

    !%Variable TDPSFKmin
    !%Type float
    !%Section Time-Dependent::TDPSF
    !%Description
    !% k-space filter width.
    !%End
    call parse_float(datasets_check('TDPSFKmin'),M_PI/width,psf%kmin)
    call messages_print_var_value(stdout, "TDPSFKmin",psf%kmin)

    kmin = psf%kmin

    sigma_min = sqrt(log(1/psf%delta)+log(width**2*lmax**(dim-1)*(M_TWO*psf%sigma)**(3*dim)/(M_PI**(3.0*dim/2.0))) )/kmin 
    
    sigma_max = width /sqrt(log(1/psf%delta+dim*log(2*psf%sigma/sqrt(M_PI))))
    
    write (*,*) "sigma_min= ",sigma_min,"sigma_max= ",sigma_max
 


    sigma = psf%sigma
    kmax = kmax/4
    call tdpsf_generate_filters(psf,width,kmin,sigma)


    POP_SUB(tdpsf_init)
  end subroutine tdpsf_init


  ! ---------------------------------------------------------
  subroutine tdpsf_end(psf)
    type(tdpsf_t), intent(inout) :: psf

    PUSH_SUB(tdpsf_end)
    SAFE_DEALLOCATE_P(psf%XFltr)
    SAFE_DEALLOCATE_P(psf%KFltr)

    POP_SUB(tdpsf_end)
  end subroutine tdpsf_end


! ---------------------------------------------------------
  subroutine tdpsf_generate_filters(psf,ww,kmin,sigma)
    type(tdpsf_t), intent(inout) :: psf
    FLOAT,         intent(in)    :: ww
    FLOAT,         intent(in)    :: kmin
    FLOAT,         intent(in)    :: sigma


    INTEGER :: ii,idir,dim
    FLOAT   :: xx,kk
    FLOAT   :: normX,normK
    FLOAT   :: LL 

    PUSH_SUB(tdpsf_generate_filters)

    dim = psf%mesh%sb%dim
    normX=M_TWO**dim *  M_PI**(-dim/M_TWO) * (M_ONE/sigma)**dim

    write (*,*) sigma 

    ! The filter in real space
    do idir=1, dim
      LL =  (psf%ll(idir)/M_TWO)*psf%mesh%spacing(idir)-ww
      do ii=1,psf%ll(idir)
        xx = (ii - psf%ll(idir)/M_TWO +M_ONE)*psf%mesh%spacing(idir)
         
        psf%XFltr(PLUS,idir,ii)  = loct_erf((LL+M_TWO*ww/3.0-xx)/sigma) - loct_erf((LL+ww/3.0-xx)/sigma) 
        psf%XFltr(MINUS,idir,ii) = loct_erf((LL+M_TWO*ww/3.0+xx)/sigma) - loct_erf((LL+ww/3.0+xx)/sigma) 

        write(*,*) ii,xx,psf%XFltr(PLUS,idir,ii)/M_TWO,psf%XFltr(MINUS,idir,ii)/M_TWO
      end do
    end do  

    psf%XFltr= psf%XFltr/M_TWO

   ! The filer in reciprocal space
    normK=M_TWO**((dim-M_ONE)/M_TWO) *  M_PI**((dim-M_ONE)/M_TWO) * (sigma)**(dim+M_ONE)

     do idir=1, dim
      do ii=1,psf%ll(idir)
        kk = pad_feq(ii, psf%mesh%idx%ll(dim), .true.)*M_PI/(psf%ll(idir)*psf%mesh%spacing(idir))

        psf%KFltr(PLUS,idir,ii)  = M_ONE + loct_erf((kk-kmin)/sigma)
        psf%KFltr(MINUS,idir,ii) = M_ONE + loct_erf((-kk-kmin)/sigma)

        write(*,*) ii,kk,psf%KFltr(PLUS,idir,ii)/M_TWO,psf%KFltr(MINUS,idir,ii)/M_TWO
      end do
    end do

    psf%KFltr= psf%KFltr/M_TWO


    POP_SUB(tdpsf_generate_filters)
  end subroutine tdpsf_generate_filters

  ! ---------------------------------------------------------
  subroutine tdpsf_project(psf,wfin,wfout,axis,direction)
    type(tdpsf_t), intent(inout) :: psf
    CMPLX,         intent(in)    :: wfin(:,:,:) 
    CMPLX,         intent(out)   :: wfout(:,:,:)
    INTEGER,       intent(in)    :: axis
    INTEGER,       intent(in)    :: direction 

    INTEGER :: idir, ix,iy,iz, dim,ivec(MAX_DIM)
    CMPLX,allocatable :: wftmp(:,:,:)


    PUSH_SUB(tdpsf_project)

    SAFE_ALLOCATE(wftmp(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))

    wfout = M_z0
    wftmp = M_z0

    do ix=1,psf%ll(1)
      do iy=1,psf%ll(2)
        do iz=1,psf%ll(3)
          ivec(1)=ix
          ivec(2)=iy
          ivec(3)=iz
!          write(*,*) "vec",ivec,"ax",axis,"dir",direction,"filter",psf%XFltr(direction,axis,ivec(axis)),wfin(ix,iy,iz)
          wfout(ix,iy,iz)=wfin(ix,iy,iz)*psf%XFltr(direction,axis,ivec(axis))
        end do
      end do
    end do

    wftmp =  M_z0
    call zfft_forward(psf%fft, wfout,wftmp)
    wfout=wftmp

    do ix=1,psf%ll(1)
      do iy=1,psf%ll(2)
        do iz=1,psf%ll(3)
          ivec(1)=ix  
          ivec(2)=iy
          ivec(3)=iz
          wfout(ix,iy,iz)=wfout(ix,iy,iz)*psf%KFltr(direction,axis,ivec(axis))
        end do
      end do
    end do

    wftmp =  M_z0
    call zfft_backward(psf%fft, wfout,wftmp)
    wfout = M_z0   
  
    do ix=1,psf%ll(1)
      do iy=1,psf%ll(2)
        do iz=1,psf%ll(3)
          ivec(1)=ix
          ivec(2)=iy
          ivec(3)=iz
          wfout(ix,iy,iz)=wfin(ix,iy,iz) - wftmp(ix,iy,iz)*psf%XFltr(direction,axis,ivec(axis))
        end do
      end do
    end do


    SAFE_DEALLOCATE_A(wftmp)

    POP_SUB(tdpsf_project)
  end subroutine tdpsf_project


  ! ---------------------------------------------------------
  subroutine tdpsf_filter_out(psf,wfin,wfout)
    type(tdpsf_t), intent(inout) :: psf
    CMPLX,         intent(in)    :: wfin(:,:,:)  
    CMPLX,         intent(out)   :: wfout(:,:,:)

    INTEGER :: idir, ii, dim
    CMPLX,allocatable :: wftmp(:,:,:)


    PUSH_SUB(tdpsf_filter_out)
      
    SAFE_ALLOCATE(wftmp(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))

    wfout = M_z0
    wftmp = M_z0

    !-----  X axis ------
    ! Plus
    call tdpsf_project(psf,wfin,wftmp,1,PLUS)
    ! Minus
    call tdpsf_project(psf,wftmp,wfout,1,MINUS)

    if(psf%mesh%sb%dim .gt. 1) then 
      !-----  Y axis ------
      ! Plus
      call tdpsf_project(psf,wfout,wftmp,2,PLUS)
      ! Minus
      call tdpsf_project(psf,wftmp,wfout,2,MINUS)

      if(psf%mesh%sb%dim .gt. 2) then 
        !-----  Z axis ------
        ! Plus
        call tdpsf_project(psf,wfout,wftmp,3,PLUS)
        ! Minus
        call tdpsf_project(psf,wftmp,wfout,3,MINUS)
      end if 
    end if 



    SAFE_DEALLOCATE_A(wftmp)


    POP_SUB(tdpsf_filter_out)
  end subroutine tdpsf_filter_out


  ! ---------------------------------------------------------
  subroutine tdpsf_project_X_to_K(psf,wfin,wfout,axis,direction)
    type(tdpsf_t), intent(in)    :: psf
    CMPLX,         intent(in)    :: wfin(:,:,:) 
    CMPLX,         intent(out)   :: wfout(:,:,:)
    INTEGER,       intent(in)    :: axis
    INTEGER,       intent(in)    :: direction 

    INTEGER :: idir, ix,iy,iz, dim,ivec(MAX_DIM)
    CMPLX,allocatable :: wftmp(:,:,:)


    PUSH_SUB(tdpsf_project_X_to_K)

    SAFE_ALLOCATE(wftmp(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))

    wfout = M_z0
    wftmp = M_z0

    do ix=1,psf%ll(1)
      do iy=1,psf%ll(2)
        do iz=1,psf%ll(3)
          ivec(1)=ix
          ivec(2)=iy
          ivec(3)=iz
!          write(*,*) "vec",ivec,"ax",axis,"dir",direction,"filter",psf%XFltr(direction,axis,ivec(axis)),wfin(ix,iy,iz)
          wfout(ix,iy,iz)=wfin(ix,iy,iz)*psf%XFltr(direction,axis,ivec(axis))
        end do
      end do
    end do

    wftmp =  M_z0
    call zfft_forward(psf%fft, wfout,wftmp)
    wfout=wftmp

    do ix=1,psf%ll(1)
      do iy=1,psf%ll(2)
        do iz=1,psf%ll(3)
          ivec(1)=ix  
          ivec(2)=iy
          ivec(3)=iz
          wfout(ix,iy,iz)=wfout(ix,iy,iz)*psf%KFltr(direction,axis,ivec(axis))
        end do
      end do
    end do


    SAFE_DEALLOCATE_A(wftmp)

    POP_SUB(tdpsf_project_X_to_K)
  end subroutine tdpsf_project_X_to_K


  ! ---------------------------------------------------------
  subroutine tdpsf_X_to_K(psf,wfin,wfout)
    type(tdpsf_t), intent(in) :: psf
    CMPLX,         intent(in)    :: wfin(:,:,:)  
    CMPLX,         intent(out)   :: wfout(:,:,:)

    INTEGER :: idir, ii, dim
    CMPLX,allocatable :: wf1(:,:,:),wf2(:,:,:),wffltrd(:,:,:)


    PUSH_SUB(tdpsf_X_to_K)
      
    SAFE_ALLOCATE(wf1(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))
    SAFE_ALLOCATE(wf2(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))
    SAFE_ALLOCATE(wffltrd(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))

    wfout = M_z0
    wf1 = M_z0
    wf2 = M_z0

    !-----  X axis ------
    ! Plus
    call tdpsf_project_X_to_K(psf,wfin,wf1,1,PLUS)
    wfout = wfout + wf1
    call tdpsf_project_K_to_X(psf,wf1 ,wf2,1,PLUS)
    wffltrd = wfin - wf2      
  
    ! Minus
    call tdpsf_project_X_to_K(psf,wffltrd,wf1,1,MINUS)
    wfout = wfout + wf1

    if(psf%mesh%sb%dim .gt. 1) then 
      !-----  Y axis ------
      call tdpsf_project_K_to_X(psf,wf1,wf2,1,MINUS)
      wffltrd = wffltrd - wf2           
      ! Plus
      call tdpsf_project_X_to_K(psf,wffltrd,wf1,2,PLUS)
      wfout = wfout + wf1
      call tdpsf_project_K_to_X(psf,wf1,wf2,2,PLUS)
      wffltrd = wffltrd - wf2
      ! Minus
      call tdpsf_project_X_to_K(psf,wffltrd,wf1,2,MINUS)
      wfout = wfout + wf1

      if(psf%mesh%sb%dim .gt. 2) then 
        !-----  Z axis ------
        call tdpsf_project_K_to_X(psf,wf1 ,wf2,2,MINUS)
        wffltrd = wffltrd - wf2
        ! Plus
        call tdpsf_project_X_to_K(psf,wffltrd,wf1,3,PLUS)
        wfout = wfout + wf1
        call tdpsf_project_K_to_X(psf,wf1,wf2,3,PLUS)
        wffltrd = wffltrd - wf2
        ! Minus
        call tdpsf_project_X_to_K(psf,wffltrd,wf1,3,MINUS)
        wfout = wfout + wf1
      end if 
    end if 



    SAFE_DEALLOCATE_A(wf1)
    SAFE_DEALLOCATE_A(wf2)
    SAFE_DEALLOCATE_A(wffltrd)


    POP_SUB(tdpsf_X_to_K)
  end subroutine tdpsf_X_to_K



  ! ---------------------------------------------------------
  subroutine tdpsf_project_K_to_X(psf,wfin,wfout,axis,direction)
    type(tdpsf_t), intent(in) :: psf
    CMPLX,         intent(in)    :: wfin(:,:,:) 
    CMPLX,         intent(out)   :: wfout(:,:,:)
    INTEGER,       intent(in)    :: axis
    INTEGER,       intent(in)    :: direction 

    INTEGER :: idir, ix,iy,iz, dim,ivec(MAX_DIM)


    PUSH_SUB(tdpsf_project_K_to_X)


    wfout =  M_z0
    call zfft_backward(psf%fft, wfin,wfout)
  
    do ix=1,psf%ll(1)
      do iy=1,psf%ll(2)
        do iz=1,psf%ll(3)
          ivec(1)=ix
          ivec(2)=iy
          ivec(3)=iz
          wfout(ix,iy,iz)= wfout(ix,iy,iz)*psf%XFltr(direction,axis,ivec(axis))
        end do
      end do
    end do



    POP_SUB(tdpsf_project_K_to_X)
  end subroutine tdpsf_project_K_to_X


  ! ---------------------------------------------------------
  subroutine tdpsf_K_to_X(psf,wfin,wfout)
    type(tdpsf_t), intent(in) :: psf
    CMPLX,         intent(in)    :: wfin(:,:,:)  
    CMPLX,         intent(out)   :: wfout(:,:,:)

    INTEGER :: idir, ii, dim
    CMPLX,allocatable :: wf1(:,:,:),wf2(:,:,:),wffltrd(:,:,:)


    PUSH_SUB(tdpsf_K_to_X)
      
    SAFE_ALLOCATE(wf1(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))
    SAFE_ALLOCATE(wf2(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))
    SAFE_ALLOCATE(wffltrd(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))

    wfout = M_z0
    wf1 = M_z0
    wf2 = M_z0

    !-----  X axis ------
    ! Plus
    call tdpsf_project_K_to_X(psf,wfin,wf1,1,PLUS)
    wfout = wfout + wf1
    call tdpsf_project_X_to_K(psf,wf1 ,wf2,1,PLUS)
    wffltrd = wfin - wf2      
  
    ! Minus
    call tdpsf_project_K_to_X(psf,wffltrd,wf1,1,MINUS)
    wfout = wfout + wf1

    if(psf%mesh%sb%dim .gt. 1) then 
      !-----  Y axis ------
      call tdpsf_project_X_to_K(psf,wf1,wf2,1,MINUS)
      wffltrd = wffltrd - wf2           
      ! Plus
      call tdpsf_project_K_to_X(psf,wffltrd,wf1,2,PLUS)
      wfout = wfout + wf1
      call tdpsf_project_X_to_K(psf,wf1,wf2,2,PLUS)
      wffltrd = wffltrd - wf2
      ! Minus
      call tdpsf_project_K_to_X(psf,wffltrd,wf1,2,MINUS)
      wfout = wfout + wf1

      if(psf%mesh%sb%dim .gt. 2) then 
        !-----  Z axis ------
        call tdpsf_project_X_to_K(psf,wf1 ,wf2,2,MINUS)
        wffltrd = wffltrd - wf2
        ! Plus
        call tdpsf_project_K_to_X(psf,wffltrd,wf1,3,PLUS)
        wfout = wfout + wf1
        call tdpsf_project_X_to_K(psf,wf1,wf2,3,PLUS)
        wffltrd = wffltrd - wf2
        ! Minus
        call tdpsf_project_K_to_X(psf,wffltrd,wf1,3,MINUS)
        wfout = wfout + wf1
      end if 
    end if 



    SAFE_DEALLOCATE_A(wf1)
    SAFE_DEALLOCATE_A(wf2)
    SAFE_DEALLOCATE_A(wffltrd)


    POP_SUB(tdpsf_K_to_X)
  end subroutine tdpsf_K_to_X


  ! ---------------------------------------------------------
  subroutine tdpsf_project_K(psf,wfin,wfout,axis,direction)
    type(tdpsf_t), intent(inout) :: psf
    CMPLX,         intent(in)    :: wfin(:,:,:) 
    CMPLX,         intent(out)   :: wfout(:,:,:)
    INTEGER,       intent(in)    :: axis
    INTEGER,       intent(in)    :: direction 

    INTEGER :: idir, ix,iy,iz, dim,ivec(MAX_DIM)
    CMPLX,allocatable :: wftmp(:,:,:)


    PUSH_SUB(tdpsf_project_K)

    SAFE_ALLOCATE(wftmp(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))

    wfout = M_z0
!    wftmp = M_z0

!    do ix=1,psf%ll(1)
!     do iy=1,psf%ll(2)
!       do iz=1,psf%ll(3)
!         ivec(1)=ix
!         ivec(2)=iy
!         ivec(3)=iz
!         wfout(ix,iy,iz)=wfin(ix,iy,iz)*psf%KFltr(direction,axis,ivec(axis))
!       end do
!     end do
!   end do

!   wftmp =  M_z0
!   call zfft_backward(psf%fft, wfout,wftmp)
!   wfout=wftmp

!   do ix=1,psf%ll(1)
!     do iy=1,psf%ll(2)
!       do iz=1,psf%ll(3)
!         ivec(1)=ix  
!         ivec(2)=iy
!         ivec(3)=iz
!         wfout(ix,iy,iz)=wfout(ix,iy,iz)*psf%XFltr(direction,axis,ivec(axis))
!       end do
!     end do
!   end do

!   wftmp =  M_z0
!   call zfft_forward(psf%fft, wfout,wftmp)
!   wfout = M_z0   
  
    do ix=1,psf%ll(1)
      do iy=1,psf%ll(2)
        do iz=1,psf%ll(3)
          ivec(1)=ix
          ivec(2)=iy
          ivec(3)=iz
!          wfout(ix,iy,iz)=wfin(ix,iy,iz) - wftmp(ix,iy,iz)*psf%KFltr(direction,axis,ivec(axis))
          wfout(ix,iy,iz)=wfin(ix,iy,iz) - wfin(ix,iy,iz)*psf%KFltr(direction,axis,ivec(axis))
        end do
      end do
    end do


    SAFE_DEALLOCATE_A(wftmp)

    POP_SUB(tdpsf_project_K)
  end subroutine tdpsf_project_K




  ! ---------------------------------------------------------
  subroutine tdpsf_filter_out_K(psf,wfin,wfout)
    type(tdpsf_t), intent(inout) :: psf
    CMPLX,         intent(in)    :: wfin(:,:,:)  
    CMPLX,         intent(out)   :: wfout(:,:,:)

    INTEGER :: idir, ii, dim
    CMPLX,allocatable :: wftmp(:,:,:)


    PUSH_SUB(tdpsf_filter_out_K)
      
    SAFE_ALLOCATE(wftmp(1:psf%mesh%idx%ll(1), 1:psf%mesh%idx%ll(2), 1:psf%mesh%idx%ll(3)))

    wfout = M_z0
    wftmp = M_z0

    !-----  X axis ------
    ! Plus
    call tdpsf_project_K(psf,wfin,wftmp,1,PLUS)
    ! Minus
    call tdpsf_project_K(psf,wftmp,wfout,1,MINUS)

    if(psf%mesh%sb%dim .gt. 1) then 
      !-----  Y axis ------
      ! Plus
      call tdpsf_project_K(psf,wfout,wftmp,2,PLUS)
      ! Minus
      call tdpsf_project_K(psf,wftmp,wfout,2,MINUS)

      if(psf%mesh%sb%dim .gt. 2) then 
        !-----  Z axis ------
        ! Plus
        call tdpsf_project_K(psf,wfout,wftmp,3,PLUS)
        ! Minus
        call tdpsf_project_K(psf,wftmp,wfout,3,MINUS)
      end if 
    end if 



    SAFE_DEALLOCATE_A(wftmp)


    POP_SUB(tdpsf_filter_out_K)
  end subroutine tdpsf_filter_out_K





end module tdpsf_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
