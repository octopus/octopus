!> @file 
!! @brief Check convolutions on GPU
!! @author
!!    Copyright (C) 2008-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! @date
!!    Septembre 2008


!> Program test for the convolution in GPU
program conv_check_fft
  use module_base
  use module_types
  use Poisson_Solver
  implicit none
  integer  :: n1,n2,n3 
  real(gp) :: hx,r2,sigma2,x,y,z,maxdiff,arg
  real(wp), dimension(:,:,:,:), allocatable :: psi_in,psi_out 
  !local variables
  character(len=*), parameter :: subname='conv_check_fft'
  character(len=50) :: chain
  integer :: i,i_stat,j,i1,i2,i3,ntimes,i_all
  integer :: l,ierror
  real(wp) :: tt,scale
  real(gp) :: CPUtime,GPUtime,ekin,ehartree
  integer, dimension(:), allocatable :: modarr
  real(kind=8), dimension(:), allocatable :: psi 
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda,v_cuda,v_cuda_str,v_cuda_str1
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda_str 
  real(kind=4), dimension(:,:,:,:), allocatable :: psi_cuda_l,v_cuda_l  
  real(kind=8) :: ekinGPUd,ehartreeGPUd
  real(kind=8) :: psi_GPU,work_GPU,work2_GPU,k_GPU
  real(kind=8) :: context,queue
  integer, parameter :: lowfilK=-14,lupfilK=14 
  real(kind=8), dimension(lowfilK:lupfilK) :: fil
  integer(kind=8) :: tsc0, tsc1
  !objects for the 3D Poisson solver
  real(kind=8), dimension(:), allocatable :: rhopot,rhopot2,rhopot1
  integer :: plan,plan1,plan1_,plan2,plan3,plan3_,size1,size2,sizek,switch_alg
  real(kind=8) :: scal
  integer :: count_0
  integer, dimension(3) :: ndim
  integer, dimension(3) :: n
  integer, dimension(3) :: geo
  real(kind=8), dimension(3) :: hgriddim
  type(coulomb_operator) :: pkernel,pkernel2

  write(chain,'(a6)') 'fort.1'
  open(unit=1,file=trim(chain),status='old',iostat=ierror)

  if (ierror==0) read(unit=1,fmt=*,iostat=ierror) n1,n2,n3,ntimes
  if (ierror /= 0) then
     write(*,*) "In a file 'fort.1', put a line with:"
     write(*,*) "n1 n2 n3 ntimes"
     write(*,*) "where:"
     write(*,*) "- n1 n2 n3 are the dimensions of the real space"
     write(*,*) "- ntimes is the number of repeated runs"
     stop
  end if

 !n1=150
 !n2=150
 !n3=150

 !ntimes=1
 

  hx=0.2d0
  !n(c) hy=0.1e0_gp
  !n(c) hz=0.1e0_gp

  hgriddim(1)=hx
  hgriddim(2)=hx
  hgriddim(3)=hx
  
   !allocate arrays
   psi_in = f_malloc((/ 2, n1, n2*n3, 1 /),id='psi_in')
   psi_out = f_malloc((/ 2, n2*n3, n1, 1 /),id='psi_out')

   psi_cuda = f_malloc((/ 2, n1, n2*n3, 1 /),id='psi_cuda')
   psi_cuda_str = f_malloc((/ 2, n2*n3, n1, 1 /),id='psi_cuda_str')
   v_cuda = f_malloc((/ 2, n2*n3, n1, 1 /),id='v_cuda')
   v_cuda_str = f_malloc((/ 2, n1, n2*n3, 2 /),id='v_cuda_str')
   v_cuda_str1 = f_malloc((/ 2, n2*n3, n1, 2 /),id='v_cuda_str1')
   psi_cuda_l = f_malloc((/ 2, n1, n2*n3, 1 /),id='psi_cuda_l')
   v_cuda_l = f_malloc((/ 2, n2*n3, n1, 1 /),id='v_cuda_l')
   rhopot = f_malloc(n1*n2*n3,id='rhopot')
   rhopot1 = f_malloc(n1*n2*n3,id='rhopot1')
   rhopot2 = f_malloc(n1*n2*n3,id='rhopot2')

   !initialise array
   sigma2=0.25d0*((n1*hx)**2)
   do i=1,n2*n3
      do i1=1,n1
        do i2=1,2
          x=hx*real(i1-n1/2-1,kind=8)
          !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
          r2=x**2
          arg=0.5d0*r2/sigma2
          tt=dexp(-arg)
          call random_number(tt)
          psi_in(i2,i1,i,1)=tt
        end do
      end do
   end do

   !initialize rhopots
   call vcopy(n1*n2*n3,psi_in(1,1,1,1),2,rhopot(1),1)
   call vcopy(n1*n2*n3,psi_in(1,1,1,1),2,rhopot2(1),1)

   !the input and output arrays must be reverted in this implementation
   do i=1,n2*n3
      do i1=1,n1
        do i2=1,2
          v_cuda_str1(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_str1(i2,i,i1,2)=0.0
          v_cuda(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_l(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
        end do
      end do
   end do

if (ntimes /=0) then

  write(*,'(a,i6,i6)')'CPU FFT, dimensions:',n1,n2*n3

   call nanosec(tsc0);
   i3=1
   call fft_1d_ctoc(-1,n2*n3,n1,v_cuda_str1,i3)
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

   write(*,'(a,i6,i6)')'GPU FFT, dimensions:',n1,n2*n3

   size1=n1*n2*n3
   size2=2*n1*n2*n3
   call cuda_1d_plan(n1,n2*n3,plan)
   call cudamalloc(size2,work_GPU,ierror)
   call cudamalloc(size2,psi_GPU,ierror)
   call reset_gpu_data(size2,psi_in,work_GPU)   

   call nanosec(tsc0);
   do i=1,ntimes
    call cuda_1d_forward(plan,work_GPU,psi_GPU)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size2,psi_cuda,psi_GPU)
   call cudafree(work_GPU)
   call cudafree(psi_GPU)
   call cufftDestroy(plan);

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_2D_cplx_results_t( n1, n2*n3,  psi_cuda, v_cuda_str1(1,1,1,i3), maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   do i=1,n2*n3
      do i1=1,n1
        do i2=1,2
          v_cuda_str(i2,i1,i,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_str(i2,i1,i,2)=0.0
          v_cuda(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_l(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
        end do
      end do
   end do


  write(*,'(a,i6,i6,i6)')'CPU 3D FFT, dimensions:',n1,n2,n3

   call nanosec(tsc0);
   i3=1
   call FFT(n1,n2,n3,n1,n2,n3,v_cuda_str,-1,i3)
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)


   write(*,'(a,i6,i6,i6)')'GPU 3D FFT, dimensions:',n1,n2,n3

   size1=n1*n2*n3
   size2=2*n1*n2*n3
   call cuda_3d_plan(n1,n2,n3,plan)
   call cudamalloc(size2,work_GPU,ierror)
   call cudamalloc(size2,psi_GPU,ierror)
   call reset_gpu_data(size2,psi_in,work_GPU)

   call nanosec(tsc0);
   do i=1,ntimes
     call cuda_3d_forward(plan,work_GPU,psi_GPU)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size2,psi_cuda,psi_GPU)

   call cudafree(work_GPU)
   call cudafree(psi_GPU)
   call cufftDestroy(plan);

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

   call compare_3D_cplx_results(n1, n2, n3, v_cuda_str(1,1,1,i3), psi_cuda, maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,5 * (log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)
 
  write(*,'(a,i6,i6,i6)')'CPU 3D Reverse FFT, dimensions:',n1,n2,n3

   call nanosec(tsc0);
   i3=1
   call FFT(n1,n2,n3,n1,n2,n3,v_cuda_str,1,i3)
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)

   write(*,'(a,i6,i6,i6)')'GPU 3D Reverse FFT, dimensions:',n1,n2,n3

   size1=n1*n2*n3
   size2=2*n1*n2*n3
   call cuda_3d_plan(n1,n2,n3,plan)
   call cudamalloc(size2,work_GPU,ierror)
   call cudamalloc(size2,psi_GPU,ierror)
   call reset_gpu_data(size2,psi_cuda,work_GPU)

   call nanosec(tsc0);
   do i=1,ntimes
     call cuda_3d_inverse(n1,n2,n3,plan,work_GPU,psi_GPU)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size2,psi_cuda,psi_GPU)

   call cudafree(work_GPU)
   call cudafree(psi_GPU)
   call cufftDestroy(plan);

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

   call compare_3D_cplx_results(n1, n2, n3, psi_in, psi_cuda, maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,5 * (log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

end if
ntimes=1   
   !Poisson Solver - periodic boundary
   ndim(1)=n1
   ndim(2)=n2
   ndim(3)=n3

   !calculate the kernel in parallel for each processor
   pkernel=pkernel_init(.true.,0,1,0,'P',ndim,hgriddim,16)
   call pkernel_set(pkernel,verbose >1)
   !pkernel%igpu=0
   !call createKernel(0,1,'P',ndim,hgriddim,16,pkernel,(verbose > 1))

   call nanosec(tsc0);
   call H_potential('D',pkernel, &
        rhopot,rhopot,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   write(*,'(a,i6,i6,i6)')'CPU 3D Poisson Solver (Periodic), dimensions:',n1,n2,n3
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)


   !here the GPU part
   pkernel2=pkernel_init(.true.,0,1,1,'P',ndim,hgriddim,16)
   call pkernel_set(pkernel2,verbose >1)

   !pkernel2%igpu=1
   !call createKernel(0,1,'P',ndim,hgriddim,16,pkernel2,(verbose > 1))

   call nanosec(tsc0);
   call H_potential('D',pkernel2, &
        rhopot2,rhopot2,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   write(*,'(a,i6,i6,i6)')'GPU 3D Poisson Solver (Periodic), dimensions:',n1,n2,n3
   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_3D_results(n1, n2, n3, rhopot(1), rhopot2(1), maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,2*5 * (log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   call pkernel_free(pkernel)
   call pkernel_free(pkernel2)

   !Poisson Solver - Free boundary

   !initialize rhopots
   call vcopy(n1*n2*n3/8,psi_in(1,1,1,1),2,rhopot(1),1)
   call vcopy(n1*n2*n3/8,psi_in(1,1,1,1),2,rhopot2(1),1)

   !calculate the kernel in parallel for each processor
   ndim(1)=n1/2
   ndim(2)=n2/2
   ndim(3)=n3/2

   pkernel=pkernel_init(.true.,0,1,0,'F',ndim,hgriddim,16)
   call pkernel_set(pkernel,verbose >1)

   !pkernel%igpu=0
   !call createKernel(0,1,'F',ndim,hgriddim,16,pkernel,(verbose > 1))
 
   call nanosec(tsc0);
    call H_potential('D',pkernel, &
        rhopot,rhopot,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);
print *,'ehartree',ehartree
   write(*,'(a,i6,i6,i6)')'CPU 3D Poisson Solver (Free), dimensions:',n1,n2,n3
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)


   !here the GPU part
   pkernel2=pkernel_init(.true.,0,1,1,'F',ndim,hgriddim,16)
   call pkernel_set(pkernel2,verbose >1)

   !pkernel2%igpu=1
   !call createKernel(0,1,'F',ndim,hgriddim,16,pkernel2,(verbose > 1))

   call nanosec(tsc0);
    call H_potential('D',pkernel2, &
        rhopot2,rhopot2,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   write(*,'(a,i6,i6,i6)')'GPU 3D Poisson Solver (Free), dimensions:',n1,n2,n3
   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
print *,'ehartree',ehartree

   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_3D_results(n1/2, n2/2, n3/2, rhopot(1), rhopot2(1), maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,2*5 * (log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   call pkernel_free(pkernel)
   call pkernel_free(pkernel2)

   !Poisson Solver - Surface boundary

   !initialize rhopots
   call vcopy(n1*n2*n3/2,psi_in(1,1,1,1),2,rhopot(1),1)
   call vcopy(n1*n2*n3/2,psi_in(1,1,1,1),2,rhopot2(1),1)

   !calculate the kernel in parallel for each processor
   ndim(1)=n1
   ndim(2)=n2/2
   ndim(3)=n3

   pkernel=pkernel_init(.true.,0,1,0,'S',ndim,hgriddim,16)
   call pkernel_set(pkernel,verbose >1)

   !pkernel%igpu=0
   !call createKernel(0,1,'S',ndim,hgriddim,16,pkernel,(verbose > 1))

   call nanosec(tsc0);
    call H_potential('D',pkernel, &
        rhopot,rhopot,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   write(*,'(a,i6,i6,i6)')'CPU 3D Poisson Solver (Surface), dimensions:',n1,n2,n3
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)

   !here the GPU part
   pkernel2=pkernel_init(.true.,0,1,1,'S',ndim,hgriddim,16)
   call pkernel_set(pkernel2,verbose >1)

!!$   !pkernel2%igpu=1
!!$   !call createKernel(0,1,'S',ndim,hgriddim,16,pkernel2,(verbose > 1))

   call nanosec(tsc0);
    call H_potential('D',pkernel2, &
        rhopot2,rhopot2,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   write(*,'(a,i6,i6,i6)')'GPU 3D Poisson Solver (Surface), dimensions:',n1,n2,n3
   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_3D_results(n1, n2/2, n3, rhopot(1), rhopot2(1), maxdiff, 3.0d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,2*5 * (log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   call pkernel_free(pkernel)
   call pkernel_free(pkernel2)


   !Poisson Solver - Wire boundary

   !initialize rhopots
   call vcopy(n1*n2*n3/4,psi_in(1,1,1,1),2,rhopot(1),1)
   call vcopy(n1*n2*n3/4,psi_in(1,1,1,1),2,rhopot2(1),1)

   !calculate the kernel in parallel for each processor
   ndim(1)=n1/2
   ndim(2)=n2/2
   ndim(3)=n3
   pkernel=pkernel_init(.true.,0,1,0,'W',ndim,hgriddim,16)
   call pkernel_set(pkernel,verbose >1)

!!$  ! pkernel%igpu=0
!!$  ! call createKernel(0,1,'W',ndim,hgriddim,16,pkernel,(verbose > 1))

   call nanosec(tsc0);
    call H_potential('D',pkernel, &
        rhopot,rhopot,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   write(*,'(a,i6,i6,i6)')'CPU 3D Poisson Solver (Wire), dimensions:',n1,n2,n3
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)

   !here the GPU part
   pkernel2=pkernel_init(.true.,0,1,1,'W',ndim,hgriddim,16)
   call pkernel_set(pkernel2,verbose >1)

!!$  ! pkernel2%igpu=1
!!$  ! call createKernel(0,1,'W',ndim,hgriddim,16,pkernel2,(verbose > 1))

   call nanosec(tsc0);
    call H_potential('D',pkernel2, &
        rhopot2,rhopot2,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);
   
  write(*,'(a,i6,i6,i6)')'GPU 3D Poisson Solver (Wire), dimensions:',n1,n2,n3
   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_3D_results(n1/2, n2/2, n3, rhopot(1), rhopot2(1), maxdiff, 3.0d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,2*5 * (log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   call pkernel_free(pkernel)
   call pkernel_free(pkernel2)

  call f_free(rhopot)
  call f_free(rhopot2)

contains

  subroutine print_time(time,nbelem,nop,ntimes)
    implicit none
    real(gp),intent(in)::time
    integer,intent(in)::nbelem,ntimes
    real(gp),intent(in)::nop

    write(*,'(a,f9.4,1pe12.5)')'Finished. Time(ms), GFlops',&
      time*1.d3/real(ntimes,kind=8),&
      real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(time*1.d9)

  end subroutine print_time

  subroutine compare_time(REFtime,TESTtime,nbelem,nop,ntimes,maxdiff,threshold)
    implicit none
    real(gp),intent(in)::REFtime,TESTtime,maxdiff,threshold
    integer,intent(in)::nbelem,ntimes
    real(gp),intent(in)::nop

    write(*,'(a,i10,f9.5,1pe12.5,2(0pf12.4,0pf12.4))',advance='no')&
      'nbelem,REF/TEST ratio,Time,Gflops: REF,TEST',&
       nbelem,REFtime/TESTtime,maxdiff,&
       REFtime*1.d3/real(ntimes,kind=8),&
       real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(REFtime*1.d9),&
       TESTtime*1.d3/real(ntimes,kind=8),&
       real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(TESTtime*1.d9)
    if (maxdiff <= threshold) then
      write(*,'(a)')''
    else
      write(*,'(a)')'<<<< WARNING' 
    end if
  end subroutine compare_time

  subroutine compare_3D_results(dim1, dim2, dim3, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2, dim3
    real(gp),intent(in):: psi_ref(dim1,dim2,dim3), psi(dim1,dim2,dim3)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2,i3

    maxdiff=0.d0
    do i3=1,dim3
      do i2=1,dim2
        do i1=1,dim1
          comp=abs(psi_ref(i1,i2,i3)-psi(i1,i2,i3))
          if(comp > printdiff) then
            write(*,*)i3,i2,i1,psi_ref(i1,i2,i3),psi(i1,i2,i3)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if
        end do
      end do
    end do
  end subroutine compare_3D_results

  subroutine compare_3D_cplx_results(dim1, dim2, dim3, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2, dim3
    real(gp),intent(in):: psi_ref(2,dim1,dim2,dim3), psi(2,dim1,dim2,dim3)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2,i3

    maxdiff=0.d0
    do i3=1,dim3
      do i2=1,dim2
        do i1=1,dim1
          comp=abs(psi_ref(1,i1,i2,i3)-psi(1,i1,i2,i3))
          if(comp > printdiff) then
            write(*,*)i3,i2,i1,1,psi_ref(1,i1,i2,i3),psi(1,i1,i2,i3)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if
          comp=abs(psi_ref(2,i1,i2,i3)-psi(2,i1,i2,i3))
          if(comp > printdiff) then
            write(*,*)i3,i2,i1,2,psi_ref(2,i1,i2,i3),psi(2,i1,i2,i3)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if

        end do
      end do
    end do
  end subroutine compare_3D_cplx_results


  subroutine compare_2D_results(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(dim1,dim2), psi(dim1,dim2)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2

    maxdiff=0.d0
    do i2=1,dim2
      do i1=1,dim1
        comp=abs(psi_ref(i1,i2)-psi(i1,i2))
        if(comp > printdiff) then
          write(*,*)i2,i1,psi_ref(i1,i2),psi(i1,i2)
        endif
        if (comp > maxdiff) then
          maxdiff=comp
        end if
      end do
    end do
  end subroutine compare_2D_results

  subroutine compare_2D_cplx_results(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(2,dim1,dim2), psi(2,dim1,dim2)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2,c

    maxdiff=0.d0
    do i2=1,dim2
      do i1=1,dim1
        do c=1,2
          comp=abs(psi_ref(c,i1,i2)-psi(c,i1,i2))
          if(comp > printdiff) then
            write(*,*)c,i2,i1,psi_ref(c,i1,i2),psi(c,i1,i2)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if
        end do
      end do
    end do
  end subroutine compare_2D_cplx_results


  subroutine compare_2D_results_t(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(dim1,dim2), psi(dim2,dim1)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2

    maxdiff=0.d0
    do i2=1,dim2
      do i1=1,dim1
        comp=abs(psi_ref(i1,i2)-psi(i2,i1))
        if(comp > printdiff) then
          write(*,*)i2,i1,psi_ref(i1,i2),psi(i2,i1)
        endif
        if (comp > maxdiff) then
          maxdiff=comp
        end if
      end do
    end do
  end subroutine compare_2D_results_t

  subroutine compare_2D_cplx_results_t(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(2,dim1,dim2), psi(2,dim2,dim1)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2,c

    maxdiff=0.d0
    do i2=1,dim2
      do i1=1,dim1
        do c=1,2
          comp=abs(psi_ref(c,i1,i2)-psi(c,i2,i1))
          if(comp > printdiff) then
            write(*,*)c,i2,i1,psi_ref(c,i1,i2),psi(c,i2,i1)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if
        enddo
      end do
    end do
  end subroutine compare_2D_cplx_results_t


  subroutine compare_1D_results(dim1, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1
    real(gp),intent(in):: psi_ref(dim1), psi(dim1)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1

    maxdiff=0.d0
    do i1=1,dim1
      comp=abs(psi_ref(i1)-psi(i1))
      if(comp > printdiff) then
        write(*,*)i1,psi_ref(i1),psi(i1)
      endif
      if (comp > maxdiff) then
        maxdiff=comp
      end if
    end do
  end subroutine compare_1D_results


  subroutine transpose_kernel_forGPU(geocode,n0,pkernel,pkernel2,offset)
    use module_base
    use Poisson_Solver
    implicit none
    integer, intent(in) :: offset
    integer, dimension(3), intent(in) :: n0
    character(len=*), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    real(dp), dimension(*), intent(in) :: pkernel
    real(dp), dimension(:), pointer :: pkernel2
    !local variables
    character(len=*), parameter :: subname='transpose_kernel_forGPU'
    integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3,nproc,i1,i2,i3,ind,indt
    integer :: j1,j2,j3,i_stat

    !only makes sense in serial for the moment
    nproc=1
    !calculate the dimensions wrt the geocode
    if (geocode == 'P') then
       call P_FFT_dimensions(n0(1),n0(2),n0(3),m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
    else if (geocode == 'S') then
       call S_FFT_dimensions(n0(1),n0(2),n0(3),m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,1,.false.)
    else if (geocode == 'F') then
       call F_FFT_dimensions(n0(1),n0(2),n0(3),m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,1,.false.)
     else if (geocode == 'W') then
       call W_FFT_dimensions(n0(1),n0(2),n0(3),m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,1,.false.)
    else
       stop 'ERROR(transpose_kernel_forGPU): geometry code not admitted'
    end if

    pkernel2 = f_malloc_ptr((n1/2+1)*n2*n3,id='pkernel2')

    do i3=1,n3
       j3=i3+(i3/(n3/2+2))*(n3+2-2*i3)!injective dimension
       do i2=1,n2
          j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)!injective dimension
          do i1=1,n1
             j1=i1+(i1/(n1/2+2))*(n1+2-2*i1)!injective dimension
             !injective index
             ind=j1+(j2-1)*nd1+(j3-1)*nd1*nd2
             !unfolded index
             !indt=j1+(i3-1)*nd1+(i2-1)*nd1*n3             
             !indt=i2+(i1-1)*n2+(i3-1)*n1*n2             
             indt=i2+(j1-1)*n2+(i3-1)*nd1*n2
             pkernel2(indt)=pkernel(ind)
          end do
       end do
    end do
    !offset to zero
    if (offset.eq.1) pkernel2(1)=0.0_dp

  end subroutine transpose_kernel_forGPU


  subroutine transpose_kernel_forGPU_cufft3D(geocode,n0,pkernel,pkernel2,offset)
    use module_base
    use Poisson_Solver
    implicit none
    integer, intent(in) :: offset
    integer, dimension(3), intent(in) :: n0
    character(len=*), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    real(dp), dimension(*), intent(in) :: pkernel
    real(dp), dimension(:), pointer :: pkernel2
    !local variables
    character(len=*), parameter :: subname='transpose_kernel_forGPU'
    integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3,nproc,i1,i2,i3,ind,indt
    integer :: j1,j2,j3,i_stat

    !only makes sense in serial for the moment
    nproc=1
    !calculate the dimensions wrt the geocode
    if (geocode == 'P') then
       call P_FFT_dimensions(n0(1),n0(2),n0(3),m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
    else if (geocode == 'S') then
       call S_FFT_dimensions(n0(1),n0(2),n0(3),m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,1,.false.)
    else if (geocode == 'F') then
       call F_FFT_dimensions(n0(1),n0(2),n0(3),m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,1,.false.)
     else if (geocode == 'W') then
       call W_FFT_dimensions(n0(1),n0(2),n0(3),m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,1,.false.)
    else
       stop 'ERROR(transpose_kernel_forGPU): geometry code not admitted'
    end if

    pkernel2 = f_malloc_ptr(n1*n2*n3,id='pkernel2')

    do i3=1,n3
       j3=i3+(i3/(n3/2+2))*(n3+2-2*i3)!injective dimension
       do i2=1,n2
          j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)!injective dimension
          do i1=1,n1
             j1=i1+(i1/(n1/2+2))*(n1+2-2*i1)!injective dimension
             !injective index
             ind=j1+(j2-1)*nd1+(j3-1)*nd1*nd2
             !unfolded index
             indt=j1+(i3-1)*nd1+(i2-1)*nd1*n3
             !indt=i2+(i1-1)*n2+(i3-1)*n1*n2             
             !indt=i2+(j1-1)*n2+(i3-1)*nd1*n2
             pkernel2(indt)=pkernel(ind)
          end do
       end do
    end do
    !offset to zero
    if (offset.eq.1) pkernel2(1)=0.0_dp

  end subroutine transpose_kernel_forGPU_cufft3D

end program conv_check_fft
