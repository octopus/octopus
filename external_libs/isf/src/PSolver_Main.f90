!> @file
!!    Main routine to perform Poisson solver calculation
!! @author
!!    Creation date: February 2007
!!    Luigi Genovese
!!    Copyright (C) 2002-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculate the Hartree potential by solving the Poisson equation 
!! @f$\nabla^2 V(x,y,z)=-4 \pi \rho(x,y,z)@f$
!! from a given @f$\rho@f$, 
!! for different boundary conditions an for different data distributions.
!! Following the boundary conditions, it applies the Poisson Kernel previously calculated.
!!    
!! @warning
!!    The dimensions of the arrays must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
!!    Moreover, for the cases with the exchange and correlation the density must be initialised
!!    to 10^-20 and not to zero.
!!
!! @todo
!!    Wire boundary condition is missing
subroutine H_potential(datacode,kernel,rhopot,pot_ion,eh,offset,sumpion,&
      quiet,stress_tensor) !optional argument
   use yaml_output
   use time_profiling, only: f_timing
   use dynamic_memory
   implicit none

   !> kernel of the Poisson equation. It is provided in distributed case, with
   !! dimensions that are related to the output of the PS_dim4allocation routine
   !! it MUST be created by following the same geocode as the Poisson Solver.
   type(coulomb_operator), intent(in) :: kernel

   !> @copydoc poisson_solver::doc::datacode
   !! To be used only in the periodic case, ignored for other boundary conditions.
   character(len=1), intent(in) :: datacode

   !> Logical value which states whether to sum pot_ion to the final result or not
   !!   .true.  rhopot will be the Hartree potential + pot_ion
   !!           pot_ion will be untouched
   !!   .false. rhopot will be only the Hartree potential
   !!           pot_ion will be ignored
   logical, intent(in) :: sumpion

   !> Total integral on the supercell of the final potential on output
   real(dp), intent(in) :: offset

   !> Hartree Energy (Hartree)
   real(gp), intent(out) :: eh

   !> On input, it represents the density values on the grid points
   !! On output, it is the Hartree potential (and maybe also pot_ion)
   real(dp), dimension(*), intent(inout) :: rhopot

   !> Additional external potential that is added to the output, 
   !! when the XC parameter ixc/=0 and sumpion=.true.
   !! When sumpion=.true., it is always provided in the distributed form,
   !! clearly without the overlapping terms which are needed only for the XC part
   real(wp), dimension(*), intent(inout) :: pot_ion

   !> Optional argument to avoid output writings
   character(len=3), intent(in), optional :: quiet

   !> Stress tensor: Add the stress tensor part from the Hartree potential
   real(dp), dimension(6), intent(out), optional :: stress_tensor

   !local variables
   character(len=*), parameter :: subname='H_potential'
   logical :: wrtmsg,cudasolver
   integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
   integer :: i_stat,ierr,ind,ind2,ind3,indp,ind2p,ind3p,i
   integer :: i1,i2,i3,j2,istart,iend,i3start,jend,jproc,i3xcsh
   integer :: nxc,istden,istglo
   real(dp) :: scal,ehartreeLOC,pot
   real(dp), dimension(6) :: strten
   real(dp), dimension(:,:,:), allocatable :: zf
   real(dp), dimension(:), allocatable :: zf1
   integer, dimension(:,:), allocatable :: gather_arr
   integer, dimension(3) :: n
   integer :: size1,size2,switch_alg

   call f_routine(id='H_potential')
   
   cudasolver=.false.
   
   !do not write anything on screen if quiet is set to yes
   if (present(quiet)) then
      if(quiet == 'yes' .or. quiet == 'YES') then
         wrtmsg=.false.
      else if(trim(quiet) == 'no' .or. trim(quiet) == 'NO') then
         wrtmsg=.true.
      else
         call yaml_warning('ERROR: Unrecognised value for "quiet" option: ' // trim(quiet))
         !write(*,*)'ERROR: Unrecognised value for "quiet" option:',quiet
         stop
      end if
   else
      wrtmsg=.true.
   end if
   wrtmsg=wrtmsg .and. kernel%mpi_env%iproc==0 .and. kernel%mpi_env%igroup==0
   ! rewrite
   if (wrtmsg) call yaml_mapping_open('Poisson Solver')
   
   !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','ON')
   call f_timing(TCAT_PSOLV_COMPUT,'ON')
   !calculate the dimensions wrt the geocode
   if (kernel%geocode == 'P') then
      if (wrtmsg) &
           call yaml_map('BC','Periodic')
      call P_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
           md1,md2,md3,nd1,nd2,nd3,kernel%mpi_env%nproc,.false.)
   else if (kernel%geocode == 'S') then
      if (wrtmsg) &
           call yaml_map('BC','Surface')
      call S_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
           md1,md2,md3,nd1,nd2,nd3,&
           kernel%mpi_env%nproc,kernel%igpu,.false.)
   else if (kernel%geocode == 'F') then
      if (wrtmsg) &
           call yaml_map('BC','Free')
      call F_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
           md1,md2,md3,nd1,nd2,nd3,&
           kernel%mpi_env%nproc,kernel%igpu,.false.)
   else if (kernel%geocode == 'W') then
      if (wrtmsg) &
           call yaml_map('BC','Wires')
      call W_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
           md1,md2,md3,nd1,nd2,nd3,kernel%mpi_env%nproc,kernel%igpu,.false.)
   else
      stop 'PSolver: geometry code not admitted'
   end if
   
   cudasolver= (kernel%igpu==1 .and. .not. present(stress_tensor))
   
   if (wrtmsg) then
      call yaml_map('Box',kernel%ndims,fmt='(i5)')
      call yaml_map('MPI tasks',kernel%mpi_env%nproc,fmt='(i5)')
      if (cudasolver) call yaml_map('GPU acceleration',.true.)
      call yaml_mapping_close()
!      call yaml_newline()
   end if
   
   if(kernel%geocode == 'P') then
      !no powers of hgrid because they are incorporated in the plane wave treatment
      scal=1.0_dp/(real(n1,dp)*real(n2*n3,dp)) !to reduce chances of integer overflow
   else if (kernel%geocode == 'S') then
      !only one power of hgrid 
      !factor of -4*pi for the definition of the Poisson equation
      scal=-16.0_dp*atan(1.0_dp)*real(kernel%hgrids(2),dp)/real(n1*n2,dp)/real(n3,dp)
   else if (kernel%geocode == 'F' .or. kernel%geocode == 'H') then
      !hgrid=max(hx,hy,hz)
      scal=product(kernel%hgrids)/real(n1*n2,dp)/real(n3,dp)
   else if (kernel%geocode == 'W') then
      !only one power of hgrid 
      !factor of -1/(2pi) already included in the kernel definition
      scal=-2.0_dp*kernel%hgrids(1)*kernel%hgrids(2)/real(n1*n2,dp)/real(n3,dp)
   end if
   !here the case ncplx/= 1 should be added
   
   !array allocations
   zf = f_malloc((/ md1, md3, 2*md2/kernel%mpi_env%nproc /),id='zf')
   !initalise to zero the zf array
   call to_zero(md1*md3*(md2/kernel%mpi_env%nproc),zf(1,1,1))
   
   istart=kernel%mpi_env%iproc*(md2/kernel%mpi_env%nproc)
   iend=min((kernel%mpi_env%iproc+1)*md2/kernel%mpi_env%nproc,m2)
   if (istart <= m2-1) then
      nxc=iend-istart
   else
      nxc=0
   end if
   
   if (datacode=='G') then
      !starting address of rhopot in the case of global i/o
      i3start=istart+1
   else if (datacode == 'D') then
      !distributed i/o
      i3start=1
   else
      stop 'PSolver: datacode not admitted'
   end if
   
   !this routine builds the values for each process of the potential (zf), multiplying by scal 
   
   !fill the array with the values of the charge density
   !no more overlap between planes
   !still the complex case should be defined
   
   do i3 = 1, nxc
      !$omp parallel do default(shared) private(i2, i1, i)
      do i2=1,m3
         do i1=1,m1
            i=i1+(i2-1)*m1+(i3+i3start-2)*m1*m3
            zf(i1,i2,i3)=rhopot(i)
         end do
      end do
      !$omp end parallel do
   end do
   
   if (.not. cudasolver) then !CPU case

      !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','OF')
      call f_timing(TCAT_PSOLV_COMPUT,'OF')
      call G_PoissonSolver(kernel%mpi_env%iproc,kernel%mpi_env%nproc,&
           kernel%part_mpi%mpi_comm,kernel%inplane_mpi%iproc,kernel%inplane_mpi%mpi_comm,kernel%geocode,1,&
           n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,kernel%kernel,&
           zf(1,1,1),&
           scal,kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3),offset,strten)
      call f_timing(TCAT_PSOLV_COMPUT,'ON')
      !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','ON')

      !check for the presence of the stress tensor
      if (present(stress_tensor)) then
         call vcopy(6,strten(1),1,stress_tensor(1),1)
      end if
   
   else !GPU case
   
      n(1)=n1!kernel%ndims(1)*(2-kernel%geo(1))
      n(2)=n3!kernel%ndims(2)*(2-kernel%geo(2))
      n(3)=n2!kernel%ndims(3)*(2-kernel%geo(3))
   
      size1=md1*md2*md3! nproc always 1 kernel%ndims(1)*kernel%ndims(2)*kernel%ndims(3)
   
      if (kernel%keepGPUmemory == 0) then
        size2=2*n1*n2*n3
        call cudamalloc(size2,kernel%work1_GPU,i_stat)
        if (i_stat /= 0) print *,'error cudamalloc',i_stat
        call cudamalloc(size2,kernel%work2_GPU,i_stat)
        if (i_stat /= 0) print *,'error cudamalloc',i_stat
      endif
   
    if (kernel%mpi_env%nproc > 1) then
      zf1 = f_malloc(md1*md3*md2,id='zf1')
   
      call mpi_gather(zf,size1/kernel%mpi_env%nproc,mpidtypd,zf1,size1/kernel%mpi_env%nproc, &
           mpidtypd,0,kernel%mpi_env%mpi_comm,ierr)
   
      if (kernel%mpi_env%iproc == 0) then
       !fill the GPU memory
   
       call reset_gpu_data(size1,zf1,kernel%work1_GPU)
   
       switch_alg=0
   
       if (kernel%initCufftPlan == 1) then
         call cuda_3d_psolver_general(n,kernel%plan,kernel%work1_GPU,kernel%work2_GPU, &
           kernel%k_GPU,switch_alg,kernel%geo,scal)
       else
         call cuda_3d_psolver_plangeneral(n,kernel%work1_GPU,kernel%work2_GPU, &
           kernel%k_GPU,kernel%geo,scal)
       endif
   
       !take data from GPU
       call get_gpu_data(size1,zf1,kernel%work1_GPU)
       endif
   
       call MPI_Scatter(zf1,size1/kernel%mpi_env%nproc,mpidtypd,zf,size1/kernel%mpi_env%nproc, &
            mpidtypd,0,kernel%mpi_env%mpi_comm,ierr)
   
       call f_free(zf1)
   
    else
   
      !fill the GPU memory
      call reset_gpu_data(size1,zf,kernel%work1_GPU)
   
      switch_alg=0
   
      if (kernel%initCufftPlan == 1) then
         call cuda_3d_psolver_general(n,kernel%plan,kernel%work1_GPU,kernel%work2_GPU, &
           kernel%k_GPU,switch_alg,kernel%geo,scal)
      else
         call cuda_3d_psolver_plangeneral(n,kernel%work1_GPU,kernel%work2_GPU, &
           kernel%k_GPU,kernel%geo,scal)
      endif
   
   
      !take data from GPU
      call get_gpu_data(size1,zf,kernel%work1_GPU)
  
   
    endif
   
    if (kernel%keepGPUmemory == 0) then
      call cudafree(kernel%work1_GPU)
      call cudafree(kernel%work2_GPU)
    endif
   
    if (kernel%keepGPUmemory == 0) then
      call cudafree(kernel%work1_GPU)
      call cudafree(kernel%work2_GPU)
    endif
   
   endif
   
   !the value of the shift depends on the distributed i/o or not
   if (datacode=='G') then
      i3xcsh=istart !beware on the fact that this is not what represents its name!!!
      !is_step=n01*n02*n03
   else if (datacode=='D') then
      i3xcsh=0 !shift not needed anymore
   end if
   
   !if (iproc == 0) print *,'n03,nxc,kernel%geocode,datacode',n03,nxc,kernel%geocode,datacode
   
   ehartreeLOC=0.0_dp
   !recollect the final data
   !this part can be eventually removed once the zf disappears
   if (sumpion) then
      do j2=1,nxc
         i2=j2+i3xcsh 
         ind3=(i2-1)*kernel%ndims(1)*kernel%ndims(2)
         ind3p=(j2-1)*kernel%ndims(1)*kernel%ndims(2)
         !$omp parallel do default(shared) private(i3, ind2, ind2p, i1, ind, indp, pot) &
         !$omp reduction(+:ehartreeLOC)
         do i3=1,m3
            ind2=(i3-1)*kernel%ndims(1)+ind3
            ind2p=(i3-1)*kernel%ndims(1)+ind3p
            do i1=1,m1
               ind=i1+ind2
               indp=i1+ind2p
               pot=zf(i1,i3,j2)
               ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
               rhopot(ind)=real(pot,wp)+real(pot_ion(indp),wp)
            end do
         end do
         !$omp end parallel do
      end do
   else
      do j2=1,nxc
         i2=j2+i3xcsh 
         ind3=(i2-1)*kernel%ndims(1)*kernel%ndims(2)
         !$omp parallel do default(shared) private(i3, ind2, i1, ind, pot) &
         !$omp reduction(+:ehartreeLOC)
         do i3=1,m3
            ind2=(i3-1)*kernel%ndims(1)+ind3
            do i1=1,m1
               ind=i1+ind2
               pot=zf(i1,i3,j2)
               ehartreeLOC=ehartreeLOC+rhopot(ind)*pot
               rhopot(ind)=real(pot,wp)
            end do
         end do
         !$omp end parallel do
      end do
   end if
   
   ehartreeLOC=ehartreeLOC*0.5_dp*product(kernel%hgrids)!hx*hy*hz
   
   call f_free(zf)
   
   !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','OF')
   !call f_timing(TCAT_PSOLV_COMPUT,'OF')

   !gathering the data to obtain the distribution array
   !evaluating the total ehartree
   eh=real(ehartreeLOC,gp)
   if (kernel%mpi_env%nproc > 1) then
      !call timing(kernel%mpi_env%iproc,'PSolv_commun  ','ON')
      !call f_timing(TCAT_PSOLV_COMPUT,'ON')

      eh=ehartreeLOC
      call mpiallred(eh,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
      !reduce also the value of the stress tensor
   
      if (present(stress_tensor)) then
         call mpiallred(stress_tensor(1),6,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
      end if
   
      !call timing(kernel%mpi_env%iproc,'PSolv_commun  ','OF')
      !call f_timing(TCAT_PSOLV_COMPUT,'OF')

      if (datacode == 'G') then
         !building the array of the data to be sent from each process
         !and the array of the displacement
   
         !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','ON')
         !call f_timing(TCAT_PSOLV_COMPUT,'ON')
         gather_arr = f_malloc((/ 0.to.kernel%mpi_env%nproc-1, 1.to.2 /),id='gather_arr')
         do jproc=0,kernel%mpi_env%nproc-1
            istart=min(jproc*(md2/kernel%mpi_env%nproc),m2-1)
            jend=max(min(md2/kernel%mpi_env%nproc,m2-md2/kernel%mpi_env%nproc*jproc),0)
            gather_arr(jproc,1)=m1*m3*jend
            gather_arr(jproc,2)=m1*m3*istart
         end do
         !gather all the results in the same rhopot array
         istart=min(kernel%mpi_env%iproc*(md2/kernel%mpi_env%nproc),m2-1)
   
         !call f_timing(TCAT_PSOLV_COMPUT,'OF')
         !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','OF')
         !call timing(kernel%mpi_env%iproc,'PSolv_commun  ','ON')
         !call f_timing(TCAT_PSOLV_COMPUT,'ON')

         istden=1+kernel%ndims(1)*kernel%ndims(2)*istart
         istglo=1
!!$        call MPI_ALLGATHERV(rhopot(istden),gather_arr(kernel%mpi_env%iproc,1),mpidtypw,&
!!$             rhopot(istglo),gather_arr(0,1),gather_arr(0,2),mpidtypw,&
!!$             kernel%mpi_env%mpi_comm,ierr)
         call mpiallgatherv(rhopot(istglo), gather_arr(:,1), gather_arr(:,2), &
              & kernel%mpi_env%iproc, kernel%mpi_env%mpi_comm,ierr)
!!$         call MPI_ALLGATHERV(MPI_IN_PLACE,gather_arr(kernel%mpi_env%iproc,1),mpidtypw,&
!!$              rhopot(istglo),gather_arr(0,1),gather_arr(0,2),mpidtypw,&
!!$              kernel%mpi_env%mpi_comm,ierr)
         !call f_timing(TCAT_PSOLV_COMPUT,'OF')
         !call timing(kernel%mpi_env%iproc,'PSolv_commun  ','OF')
         !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','ON')
         !call f_timing(TCAT_PSOLV_COMPUT,'ON')
         call f_free(gather_arr)
         !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','OF')
     
      end if
   end if

   call f_timing(TCAT_PSOLV_COMPUT,'OF')
   call f_release_routine()

END SUBROUTINE H_potential


!> Calculate the dimensions needed for the allocation of the arrays 
!! related to the Poisson Solver
!!
!! @warning
!!    The XC enlarging due to GGA part is not present for surfaces and 
!!    periodic boundary condition. This is related to the fact that the calculation of the
!!    gradient and the White-Bird correction are not yet implemented for non-isolated systems
!! @author Luigi Genovese
!! @date February 2007
subroutine PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,use_gradient,use_wb_corr,&
      n3d,n3p,n3pi,i3xcsh,i3s)
   implicit none
   character(len=1), intent(in) :: geocode  !< @copydoc poisson_solver::doc::geocode
   character(len=1), intent(in) :: datacode !< @copydoc poisson_solver::doc::datacode
   integer, intent(in) :: iproc        !< Process Id
   integer, intent(in) :: nproc        !< Number of processes
   integer, intent(in) :: n01,n02,n03  !< Dimensions of the real space grid to be hit with the Poisson Solver
   logical, intent(in) :: use_gradient !< .true. if functional is using the gradient.
   logical, intent(in) :: use_wb_corr  !< .true. if functional is using WB corrections.
   !> Third dimension of the density. For distributed data, it takes into account 
   !! the enlarging needed for calculating the XC functionals.
   !! For global data it is simply equal to n03. 
   !! When there are too many processes and there is no room for the density n3d=0
   integer, intent(out) :: n3d
   !> Third dimension for the potential. The same as n3d, but without 
   !! taking into account the enlargment for the XC part. For non-GGA XC, n3p=n3d.
   integer, intent(out) :: n3p
   !> Dimension of the pot_ion array, always with distributed data. 
   !! For distributed data n3pi=n3p
   integer, intent(out) :: n3pi
   !> Shift of the density that must be performed to enter in the 
   !! non-overlapping region. Useful for recovering the values of the potential
   !! when using GGA XC functionals. If the density starts from rhopot(1,1,1),
   !! the potential starts from rhopot(1,1,i3xcsh+1). 
   !! For non-GGA XCs and for global distribution data i3xcsh=0
   integer, intent(out) :: i3xcsh
   !> Starting point of the density effectively treated by each processor 
   !! in the third direction.
   !! It takes into account also the XC enlarging. The array rhopot will correspond
   !! To the planes of third coordinate from i3s to i3s+n3d-1. 
   !! The potential to the planes from i3s+i3xcsh to i3s+i3xcsh+n3p-1
   !! The array pot_ion to the planes from i3s+i3xcsh to i3s+i3xcsh+n3pi-1
   !! For global disposition i3s is equal to distributed case with i3xcsh=0.
   integer, intent(out) :: i3s
   !local variables
   !n(c) integer, parameter :: nordgr=4
   integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
   integer :: istart,iend,nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr
   integer :: n3pr1,n3pr2
   
   
   !calculate the dimensions wrt the geocode
   if (geocode == 'P') then
      call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
       if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   else if (geocode == 'S') then
      call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   else if (geocode == 'F' .or. geocode == 'H') then
      call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2/2) then
           md2=(md2/nproc+1)*nproc 
        endif
      endif
   else if (geocode == 'W') then
      call W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   else
      write(*,*) geocode
      stop 'PS_dim4allocation: geometry code not admitted'
   end if
   
   !formal start and end of the slice
   istart=iproc*(md2/nproc)
   iend=min((iproc+1)*md2/nproc,m2)
   
   if (datacode == 'D') then
      call xc_dimensions(geocode,use_gradient,use_wb_corr,istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh)
   
      nwb=nxcl+nxc+nxcr-2
      nxt=nwbr+nwb+nwbl
   
      n3p=nxc
      n3d=nxt
      n3pi=n3p
   else if (datacode == 'G') then
      n3d=n03
      n3p=n03
      i3xcsh=0
      i3s=min(istart,m2-1)+1
      n3pi=max(iend-istart,0)
   else
      print *,datacode
      stop 'PS_dim4allocation: data code not admitted'
   end if

!!  print *,'P4,iproc',iproc,'nxc,ncxl,ncxr,nwbl,nwbr',nxc,nxcl,nxcr,nwbl,nwbr,&
!!       'ixc,n3d,n3p,i3xcsh,i3s',ixc,n3d,n3p,i3xcsh,i3s

END SUBROUTINE PS_dim4allocation


!> Calculate the dimensions to be used for the XC part, taking into account also
!! the White-bird correction which should be made for some GGA functionals
!!
!! SYNOPSIS
!!    @param use_gradient .true. if functional is using the gradient.
!!    @param use_wb_corr  .true. if functional is using WB corrections.
!!    @param m2        dimension to be parallelised
!!    @param nxc       size of the parallelised XC potential
!!    @param ncxl,ncxr left and right buffers for calculating the WB correction after call drivexc
!!    @param nwbl,nwbr left and right buffers for calculating the gradient to pass to drivexc    
!!    @param i3s       starting addres of the distributed dimension
!!    @param i3xcsh    shift to be applied to i3s for having the striting address of the potential
!!
!! @warning It is imperative that iend <=m2
subroutine xc_dimensions(geocode,use_gradient,use_wb_corr,&
     & istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh)
  implicit none

  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: use_gradient, use_wb_corr
  integer, intent(in) :: istart,iend,m2
  integer, intent(out) :: nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh
  !local variables
  integer, parameter :: nordgr=4 !the order of the finite-difference gradient (fixed)

  if (istart <= m2-1) then
     nxc=iend-istart
     if (use_gradient .and. geocode == 'F') then
        if (.not. use_wb_corr) then
           !now the dimension of the part required for the gradient
           nwbl=min(istart,nordgr)
           nwbr=min(m2-iend,nordgr) !always m2 < iend
           nxcl=1
           nxcr=1
        else
           !now the dimension of the part required for the gradient
           if(istart<=nordgr) then
              nxcl=istart+1
              nwbl=0
           else
              nxcl=nordgr+1
              nwbl=min(nordgr,istart-nordgr)
           end if
           if(iend>=m2-nordgr+1) then
              nxcr=m2-iend+1
              nwbr=0
           else
              nxcr=nordgr+1
              nwbr=min(nordgr,m2-nordgr-iend)
           end if
        end if
     else if (geocode /= 'F' .and. use_gradient .and. nxc /= m2) then
        if (.not. use_wb_corr) then
           !now the dimension of the part required for the gradient
           nwbl=nordgr
           nwbr=nordgr
           nxcl=1
           nxcr=1
        else
           nxcl=nordgr+1
           nwbl=nordgr
           nxcr=nordgr+1
           nwbr=nordgr
        end if
     !this case is also considered below
     !else if (geocode /= 'F' .and. use_gradient .and. nxc == m2) then
     else 
        nwbl=0
        nwbr=0
        nxcl=1
        nxcr=1
     end if
     i3xcsh=nxcl+nwbl-1
     i3s=istart+1-i3xcsh
  else
     nwbl=0
     nwbr=0
     nxcl=1
     nxcr=1
     nxc=0
     i3xcsh=0
     i3s=m2
  end if
END SUBROUTINE xc_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the periodic system
!!
!!    @param n01,n02,n03 original real dimensions (input)
!!
!!    @param m1,m2,m3 original real dimension, with m2 and m3 exchanged
!!
!!    @param n1,n2,n3 the first FFT dimensions (even for the moment - the medium point being n/2+1)
!!
!!    @param md1,md2,md3 the n1,n2,n3 dimensions. They contain the real unpadded space.
!!           !!           md2 is further enlarged to be a multiple of nproc
!!
!!    @param nd1,nd2,nd3 fourier dimensions for which the kernel is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the other geometries of the Poisson Solver.
!!    The dimensions 2 and 3 are exchanged.
!! @author Luigi Genovese
!! @date October 2006
subroutine P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,enlarge_md2)
  !use module_defs, only: md2plus
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3
 
 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 

 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3

 call fourier_dim(l1,n1)
 if (n1 /= m1) then
    print *,'the FFT in the x direction is not allowed'
    print *,'n01 dimension',n01
    stop
 end if

 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    print *,'the FFT in the z direction is not allowed'
    print *,'n03 dimension',n03
    stop
 end if
 
 call fourier_dim(l3,n3)
 if (n3 /= m3) then
    print *,'the FFT in the y direction is not allowed'
    print *,'n02 dimension',n02
    stop
 end if

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3

 !enlarge the md2 dimension to be compatible with MPI_ALLTOALL communication
 do while(nproc*(md2/nproc) < n2)
    !151 if (nproc*(md2/nproc) < n2) then
    md2=md2+1
 end do
!    goto 151
 !endif
 
 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1 
 nd2=n2/2+1
 nd3=n3/2+1
 
 !enlarge the md2 dimension to be compatible with MPI_ALLTOALL communication
 do while(modulo(nd3,nproc) /= 0)
!250 if (modulo(nd3,nproc) /= 0) then
    nd3=nd3+1
!    goto 250
! endif
 end do

END SUBROUTINE P_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the surface system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                which has been properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!! @author Luigi Genovese
!! @date October 2006
subroutine S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2)
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc,gpu
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
 else
  l3=2*m3
 endif

 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3

 call fourier_dim(l1,n1)
 if (n1 /= m1) then
    print *,'the FFT in the x direction is not allowed'
    print *,'n01 dimension',n01
    stop
 end if
 
 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    print *,'the FFT in the z direction is not allowed'
    print *,'n03 dimension',n03
    stop
 end if

 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2)
    !151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    !goto 151
    !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1
 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !goto 250
    !endif
 end do

END SUBROUTINE S_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the Wires BC system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                which has been properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!! @author Luigi Genovese
!! @date October 2006
subroutine W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2)
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc,gpu
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
 else
  l3=2*m3
 endif

 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
       exit
    end if
    l1=l1+1
 end do

 
 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    print *,'the FFT in the z direction is not allowed'
    print *,'n03 dimension',n03
    stop
 end if

 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do

 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2)
    !151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    !goto 151
    !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1
 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !goto 250
    !endif
 end do

END SUBROUTINE W_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! zero-padded convolution
!!
!!    @param n01,n02,n03 original real dimensions (input)
!!
!!    @param m1,m2,m3 original real dimension with the dimension 2 and 3 exchanged
!!
!!    @param n1,n2 the first FFT even dimensions greater that 2*m1, 2*m2
!!    @param n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    @param md1,md2,md3 half of n1,n2,n3 dimension. They contain the real unpadded space,
!!                which has been properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    @param nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    The dimension m2 and m3 correspond to n03 and n02 respectively
!!    this is needed since the convolution routine manage arrays of dimension
!!    (md1,md3,md2/nproc)
!! @author Luigi Genovese
!! @date February 2006
subroutine F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2)
  !use module_defs, only: md2plus
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc,gpu
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space, inverted for convenience
 m1=n01
 m2=n03
 m3=n02
 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=2*m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
 else
  l3=2*m3
 endif
 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3
 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
       exit
    end if
    l1=l1+1
 end do
 do
    call fourier_dim(l2,n2)
    if (modulo(n2,2) == 0) then
       exit
    end if
    l2=l2+1
 end do
 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2/2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2/2)
   !151 if (nproc*(md2/nproc).lt.n2/2) then
    md2=md2+1
   !goto 151
   !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1

 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !    goto 250
    ! endif
 end do

END SUBROUTINE F_FFT_dimensions
