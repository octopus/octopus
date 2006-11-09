!  Copyright (C) Luigi Genovese, Thierry Deutsch, CEA Grenoble, 2006
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

!here we have all the subroutines needed for the parallel convolution

subroutine PARtest_kernel(n01,n02,n03,nfft1,nfft2,nfft3,&
    hgrid,karray,rhopot,iproc,nproc)
 implicit none
 !Arguments
 integer :: n01,n02,n03,nfft1,nfft2,nfft3,iproc,nproc
 real*8 :: hgrid
 real*8, dimension(nfft1,nfft2,nfft3/nproc) :: karray
 real*8, dimension(n01,n02,n03) :: rhopot
 !Local variables
 real*8 :: a_gauss,a2
 real*8 :: rhotot,shft1,shft2,shft3,ehart
 real*8 :: pi,x1,x2,x3,r,r2,factor,derf,max_diff,diff,tt
 integer :: i1,i2,i3,ii1,ii2,ii3

 a_gauss=4.d0*hgrid
 a2 = a_gauss**2

 write(*,*) 'test_kernel, dim kernel',nfft1/2+1,nfft2/2+1,nfft3/2+1

 !Shift the center of the Gaussian
 !away from central grid point to break symmetries
 !shft1=1.3d0*hgrid
 !shft2=0.5d0*hgrid
 !shft3=0.1d0*hgrid
 shft1=0.d0
 shft2=0.d0
 shft3=0.d0

 !Initialization
 pi = 4.d0*atan(1.d0)
 !Normalization
 factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
 !Gaussian function
 rhotot=0.d0
 do i3=1,n03
    x3 = hgrid*(i3-n03/2)-shft3
    do i2=1,n02
       x2 = hgrid*(i2-n02/2)-shft2
       do i1=1,n01
          x1 = hgrid*(i1-n01/2)-shft1
          r2 = x1*x1+x2*x2+x3*x3
          rhopot(i1,i2,i3) = factor*exp(-r2/a2)
          rhotot=rhotot+rhopot(i1,i2,i3)
       end do
    end do
 end do
 rhotot=rhotot*hgrid**3

 !! Plot values along x axis
 !   open(unit=11,file='rho.dat')
 !   do i3=1,n03
 !      x1 =                            - shft1
 !      x2 =                            - shft2
 !      x3 = hgrid*(i3-n03/2) - shft3
 !      r=sqrt(x1**2+x2**2+x3**2)
 !      write(unit=11,fmt="(e10.3,e12.5,2(e21.14),e9.2,2(e12.5))") &
 !           r, rhopot(n01/2,n02/2,i3)
 !   end do
 !   close(unit=11)

  !Calculate potential using Poisson Solver
     write(*,*) 'testing poisson solver'
     call ParPSolver_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,&
          hgrid,karray,rhopot,ehart,iproc,nproc)

  !! Plot values along x axis
  !   open(unit=11,file='pot.dat')
  !   do i3=1,n03
  !      x1 =                            - shft1
  !      x2 =                            - shft2
  !      x3 = hgrid*(i3-n03/2) - shft3
  !      r=sqrt(x1**2+x2**2+x3**2)
  !      if (r == 0.d0) then
  !         !limit_{x -> 0} erf(x/x) = 2/sqrt(pi)
  !         factor = 2.d0/(sqrt(pi)*a_gauss)
  !         tt=0.d0
  !      else
  !         factor = derf(r/a_gauss)/r
  !         tt=abs(1.d0/r)
  !      end if
  !      write(unit=11,fmt="(e10.3,3(e21.14),e9.2,2(e12.5))") &
  !           r, rhopot(n01/2,n02/2,i3),factor,tt
  !   end do
  !   close(unit=11)



  ! Global error
  max_diff = 0.d0
  do i3=1,n03
     x3 = hgrid*(i3-n03/2) - shft3
     do i2=1,n02
        x2 = hgrid*(i2-n02/2) - shft2
        do i1=1,n01
           x1 = hgrid*(i1-n01/2) - shft1
           r=sqrt(x1**2+x2**2+x3**2)
           if (r == 0.d0) then
              !limit_{x -> 0} erf(x/x) = 2/sqrt(pi)
              factor = 2.d0/(sqrt(pi)*a_gauss)
           else
              factor = derf(r/a_gauss)/r
           end if
           diff=abs(rhopot(i1,i2,i3)-factor)
           if (diff.gt.max_diff) then
              max_diff=diff
              ii1=i1
              ii2=i2
              ii3=i3
           endif
        end do
     end do
  end do

  write(*,*) 'Testing Poisson Solver for a_gauss=',a_gauss
  write(*,'(1x,a,f7.2,1x,e10.3,1x,e10.3)') &
       'hgridh,Deltarho,max_diff',hgrid,rhotot-1.d0,max_diff
  write(*,*) 'Max diff at : ',ii1,ii2,ii3
  write(*,*) 'Poisson Solver test finished'


end subroutine PARtest_kernel
!!***

!!****h* BigDFT/calculate_pardimensions
!! NAME
!!   calculate_pardimensions
!!
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    zero-padded convolution
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension with the dimension 2 and 3 exchanged
!!
!!    n1,n2 the first FFT even dimensions greater that 2*m1, 2*m2
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2,md3 half of n1,n2,n3 dimension. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    The dimension m2 and m3 correspond to n03 and n02 respectively
!!    this is needed since the convolution routine manage arrays of dimension
!!    (md1,md3,md2/nproc)
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
subroutine calculate_pardimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3,l1A,l3A

 !dimensions of the density in the real space, inverted for convenience

 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)

 !     n1=2*m1; n2=2*m2; n3=2*m3

 l1=2*m1
 l2=2*m2
 l3=m3 !beware of the half dimension
 do
    !this is for the FFT of the kernel
    !one can erase it when the kernel is parallelized
    call fourier_dim(l1,n1)
    !call fourier_dim(n1/2,l1A)
    if (modulo(n1,2) == 0&! .and. 2*l1A == n1
         ) then
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
    !call fourier_dim(n3/2,l3A)
    if (modulo(n3,2) == 0 &!.and. 2*l3A == n3 .and. modulo(l3A,2) == 0
         ) then
       exit
    end if
    l3=l3+1
 end do
 n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2/2
 md3=n3/2
151 if (nproc*(md2/nproc).lt.n2/2) then
    md2=md2+1
    goto 151
 endif


 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 nd1=n1/2+1;  nd2=n2/2+1
 nd3=n3/2+1
250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine calculate_pardimensions
!!***


!!****h* BigDFT/ParPSolver_Kernel
!! NAME
!!   ParPSolver_Kernel
!!
!! FUNCTION
!!    Solver of Poisson equation applying a kernel, parallel computation
!!
!! SYNOPSIS
!!    Poisson solver applying a kernel and
!!    using Fourier transform for the convolution.
!!    rhopot : input  -> the density
!!             output -> the Hartree potential + pot_ion
!!    All the processes manage the same global rhopot array
!!    The potential pot_ion is ADDED in the array rhopot.
!!    Calculate also the Hartree potential
!!
!!    Replaces the charge density contained in rhopot
!!    by the Hartree stored as well in rhopot.
!!    If xc_on is true, it also adds the XC potential and
!!    ionic potential pot_ion
!!
!!    kernelLOC: the kernel in fourier space, calculated from ParBuild_Kernel routine
!!               it is a local vector (each process have its own part)
!!
!!    We double the size of the mesh except in one dimension
!!    in order to use the property of the density to be real.
!! WARNING
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
subroutine ParPSolver_Kernel(n01,n02,n03,nd1,nd2,nd3, &
    hgrid,kernelLOC,rhopot,ehartree,iproc,nproc)
 implicit none
 !Arguments
 integer, intent(in)  :: n01,n02,n03,nd1,nd2,nd3,iproc,nproc
 real*8, intent(in) :: hgrid
 real*8, intent(in), dimension(nd1,nd2,nd3/nproc) :: kernelLOC
 real*8, intent(inout), dimension(n01,n02,n03) :: rhopot
 real*8, intent(out) :: ehartree
 !Local variables
 integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3

 call calculate_pardimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)

 call  pconvxc_off(m1,m2,m3,n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,iproc,nproc,&
      rhopot,kernelLOC,hgrid,ehartree)
 if (iproc.eq.0) print *,"ehartree is",ehartree

end subroutine ParPSolver_Kernel
!!***

!!****h* BigDFT/pconvxc_off
!! NAME
!!   pconvxc_off
!!
!! FUNCTION
!!    Calculate the parallel convolution with the kernel
!!    without the exchange-correlation part
!!
!! SYNOPSIS
!!    Poisson solver applying a kernel and
!!    using Fourier transform for the convolution.
!!    rhopot : input  -> the density
!!             output -> the Hartree potential + pot_ion
!!    All the processes manage the same global rhopot array
!!    The potential pot_ion is ADDED in the array rhopot.
!!    Calculate also the Hartree potential
!!
!!    Replaces the charge density contained in rhopot
!!    by the Hartree stored as well in rhopot.
!!
!!    kernelLOC: the kernel in fourier space, calculated from ParBuild_Kernel routine
!!               it is a local vector (each process have its own part)
!!
!!    We double the size of the mesh except in one dimension
!!    in order to use the property of the density to be real.
!! WARNING
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
subroutine pconvxc_off(m1,m2,m3,n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,iproc,nproc,&
    rhopot,kernelloc,hgrid,ehartree)
 implicit none
 include 'mpif.h'
 integer, intent(in) :: m1,m2,m3,n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,iproc,nproc
 real(kind=8), dimension(nd1,nd2,nd3/nproc), intent(in) :: kernelloc
 real(kind=8), dimension(m1,m3,m2), intent(inout) :: rhopot
 real(kind=8), intent(in) :: hgrid
 real(kind=8), intent(out) :: ehartree
 !Local variables
 integer :: ierr,istart,iend,jend,jproc,i_allocated,i_stat
 real(kind=8) :: ehartreeLOC,scal
 real(kind=8), dimension(:,:,:), allocatable :: zf
 integer, dimension(:,:), allocatable :: gather_arr
 integer count1,count2,count_rate,count_max
 real(kind=8) t1,t0,tel

 !factor to be used to keep unitarity
 scal=hgrid**3/real(n1*n2*n3,kind=8)

 i_allocated=0
 allocate(zf(md1,md3,md2/nproc),stat=i_stat)
 i_allocated=i_allocated+i_stat
 allocate(gather_arr(0:nproc-1,2),stat=i_stat)
 i_allocated=i_allocated+i_stat
 if (i_allocated /= 0) then
    print *,"pconvxc_off:Problem of memory allocation"
    stop
 end if

 !Here we insert the process-related values of the density, starting from the total density
 call enterdensity(rhopot(1,1,1),m1,m2,m3,md1,md2,md3,iproc,nproc,zf(1,1,1))

 !this routine builds the values for each process of the potential (zf), multiplying by the factor
 call convolxc_off(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,kernelloc,zf&
      ,scal,hgrid,ehartreeLOC)

 !evaluating the total ehartree
 if (nproc.gt.1) then
      call cpu_time(t0)
      call system_clock(count1,count_rate,count_max)
    call MPI_ALLREDUCE(ehartreeLOC,ehartree,1,MPI_double_precision,  &
         MPI_SUM,MPI_COMM_WORLD,ierr)
      call cpu_time(t1)
      call system_clock(count2,count_rate,count_max)
      tel=dble(count2-count1)/dble(count_rate)
      write(78,*) 'PSolver: ALLREDUCE TIME',iproc,t1-t0,tel
      write(78,*) '----------------------------------------------'
 else
    ehartree=ehartreeLOC
 end if

 !building the array of the data to be sent from each process
 !and the array of the displacement
 do jproc=0,nproc-1
    istart=min(jproc*(md2/nproc),m2-1)
    jend=max(min(md2/nproc,m2-md2/nproc*jproc),0)
    gather_arr(jproc,1)=m1*m3*jend
    gather_arr(jproc,2)=istart*m1*m3
 end do

 !assign the distributed density to the rhopot array
 istart=min(iproc*(md2/nproc),m2-1)
 jend=max(min(md2/nproc,m2-md2/nproc*iproc),0)
 iend=istart+jend
 rhopot(:,:,istart+1:iend)=&
      zf(1:m1,1:m3,1:jend)

 !gather all the results in the same rhopot array
      call cpu_time(t0)
      call system_clock(count1,count_rate,count_max)
 call MPI_ALLGATHERV(rhopot(1,1,istart+1),gather_arr(iproc,1),MPI_double_precision,rhopot,gather_arr(:,1),&
      gather_arr(:,2),MPI_double_precision,MPI_COMM_WORLD,ierr)
      call cpu_time(t1)
      call system_clock(count2,count_rate,count_max)
      tel=dble(count2-count1)/dble(count_rate)
      write(78,*) 'PSolver: ALLGATHERV TIME',iproc,t1-t0,tel
      write(78,*) '----------------------------------------------'

 deallocate(zf,gather_arr)

end subroutine pconvxc_off
!!***


!!****h* BigDFT/enterdensity
!! NAME
!!   enterdensity
!!
!! FUNCTION
!!
!!   Define a real space process-dependent vector zf with the global dimensions that are half of the FFT grid
!!   in order to perform convolution. The dimension md2 is a multiple of nproc
!!   Can be used also to define the local part of pot_ion
!!
!! AUTHOR
!!    L. Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
subroutine enterdensity(rhopot,m1,m2,m3,md1,md2,md3,iproc,nproc,zf)
 implicit none
 !Arguments
 integer, intent(in) :: m1,m2,m3,md1,md2,md3,iproc,nproc
 real*8, dimension(0:md1-1,0:md3-1,0:md2/nproc-1), intent(out) :: zf
 real*8, dimension(0:m1-1,0:m3-1,0:m2-1), intent(in) :: rhopot
 !Local Variables
 integer :: j1,j2,j3,jp2
 real*8 :: flag

 !Body
 do jp2=0,md2/nproc-1
    j2=iproc*(md2/nproc)+jp2
    if (j2.le.m2-1) then
       do j3=0,m3-1
          do j1=0,m1-1
             zf(j1,j3,jp2)=rhopot(j1,j3,j2)
          end do
          do j1=m1,md1-1
             zf(j1,j3,jp2)=0.d0
          end do
       end do
       do j3=m3,md3-1
          do j1=0,md1-1
             zf(j1,j3,jp2)=0.d0
          end do
       end do
    else
       do j3=0,md3-1
          do j1=0,md1-1
             zf(j1,j3,jp2)=0.d0
          end do
       end do
    end if
 end do

end subroutine enterdensity
!!***


!!****h* BigDFT/ParBuild_Kernel
!! NAME
!!   ParBuild_Kernel
!!
!! FUNCTION
!!    Build the kernel of a gaussian function
!!    for interpolating scaling functions.
!!    Do the parallel HalFFT of the symmetrized function and stores into
!!    memory only 1/8 of the grid divided by the number of processes nproc
!!
!! SYNOPSIS
!!    Build the kernel (karray) of a gaussian function
!!    for interpolating scaling functions
!!    $$ K(j) = \sum_k \omega_k \int \int \phi(x) g_k(x'-x) \delta(x'- j) dx dx' $$
!!
!!    n01,n02,n03        Mesh dimensions of the density
!!    nfft1,nfft2,nfft3  Dimensions of the FFT grid (HalFFT in the third direction)
!!    n1k,n2k,n3k        Dimensions of the kernel FFT
!!    hgrid              Mesh step
!!    itype_scf          Order of the scaling function (8,14,16)
!!
!! AUTHORS
!!    T. Deutsch, L. Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
subroutine ParBuild_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,hgrid,itype_scf,iproc,nproc,karrayoutLOC)

 implicit none

 !Arguments
 integer, intent(in) :: n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,itype_scf,iproc,nproc
 real(kind=8), intent(in) :: hgrid
 real(kind=8), dimension(n1k,n2k,n3k/nproc), intent(out) :: karrayoutLOC

 !Local variables
 !Do not touch !!!!
 integer, parameter :: n_gauss = 89
 !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
 integer, parameter :: n_points = 2**6

 !Better p_gauss for calculation
 !(the support of the exponential should be inside [-n_range/2,n_range/2])
 real(kind=8), parameter :: p0_ref = 1.d0
 real(kind=8), dimension(n_gauss) :: p_gauss,w_gauss

 real(kind=8), dimension(:), allocatable :: kernel_scf,kern_1_scf
 real(kind=8), dimension(:), allocatable :: x_scf ,y_scf
 real(kind=8), dimension(:,:,:,:), allocatable :: karrayfour
 real(kind=8), dimension(:,:,:), allocatable :: karray

 real(kind=8) :: ur_gauss,dr_gauss,acc_gauss,pgauss,kern,a_range,kern_tot
 real(kind=8) :: pi,factor,factor2,urange,dx,absci,p0gauss,weight,p0_cell
 real(kind=8) :: a1,a2,a3,amax
 integer :: n_scf,nker1,nker2,nker3
 integer :: i_gauss,n_range,n_cell,istart,iend,istart1,istart2,iend1,iend2
 integer :: i,j,n_iter,i_iter,ind,i1,i2,i3,i_kern,i_stat,i_allocated
 integer :: i01,i02,i03,n1h,n2h,n3h

 !Number of integration points : 2*itype_scf*n_points
 n_scf=2*itype_scf*n_points
 !Set karray
 !dimension test

 !here we must set the dimensions for the fft part, starting from the nfft
 !remember that actually nfft2 is associated to n03 and viceversa

 !dimensions that define the center of symmetry
 n1h=nfft1/2
 n2h=nfft2/2
 n3h=nfft3/2

 !Auxiliary dimensions only for building the FFT part
 nker1=nfft1
 nker2=nfft2
 nker3=nfft3/2+1

 !adjusting the last two dimensions to be multiples of nproc
 do
    if(modulo(nker2,nproc) == 0) exit
    nker2=nker2+1
 end do
 do
    if(modulo(nker3,nproc) == 0) exit
    nker3=nker3+1
 end do

 !this will be the array of the kernel in the real space
 allocate(karray(nker1,nfft3,nker2/nproc))

 !defining proper extremes for the calculation of the
 !local part of the kernel

 istart=iproc*nker2/nproc+1
 iend=min((iproc+1)*nker2/nproc,n2h+n03)

 istart1=istart
 if(iproc .eq. 0) istart1=n2h-n03+2

 iend2=iend

 iend1=n2h
 istart2=n2h+1
 if(istart .gt. n2h) then
    iend1=istart1-1
    istart2=istart
 end if
 if(iend .le. n2h) then
    istart2=iend2+1
    iend1=iend
 end if

!!!!!START KERNEL CONSTRUCTION
 if(iproc .eq. 0) then
    write(unit=*,fmt="(1x,a,i0,a)") &
         "Build the kernel in parallel using a sum of ",n_gauss," gaussians"
    write(unit=*,fmt="(1x,a,i0,a)") &
         "Use interpolating scaling functions of ",itype_scf," order"
 end if

 !Allocations
 i_allocated = 0
 allocate(x_scf(0:n_scf),stat=i_stat)
 i_allocated = i_allocated + i_stat
 allocate(y_scf(0:n_scf),stat=i_stat)
 i_allocated = i_allocated + i_stat
 if (i_allocated /= 0) then
    print *,"Build_Kernel: Problem of memory allocation"
    stop
 end if

 !Build the scaling function
 call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
 !Step grid for the integration
 dx = real(n_range,kind=8)/real(n_scf,kind=8)
 !Extend the range (no more calculations because fill in by 0.d0)
 n_cell = max(n01,n02,n03)
 n_range = max(n_cell,n_range)

 !Allocations
 allocate(kernel_scf(-n_range:n_range),stat=i_stat)
 i_allocated = i_allocated + i_stat
 allocate(kern_1_scf(-n_range:n_range),stat=i_stat)
 i_allocated = i_allocated + i_stat
 if (i_allocated /= 0) then
    print *,"Build_Kernel: Problem of memory allocation"
    stop
 end if

 !Lengthes of the box (use FFT dimension)
 a1 = hgrid * real(n01,kind=8)
 a2 = hgrid * real(n02,kind=8)
 a3 = hgrid * real(n03,kind=8)

 x_scf(:) = hgrid * x_scf(:)
 y_scf(:) = 1.d0/hgrid * y_scf(:)
 dx = hgrid * dx
 !To have a correct integration
 p0_cell = p0_ref/(hgrid*hgrid)

 !Initialization of the gaussian (Beylkin)
 call gequad(n_gauss,p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
 !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3)
 !(biggest length in the cube)
 !We divide the p_gauss by a_range**2 and a_gauss by a_range
 a_range = sqrt(a1*a1+a2*a2+a3*a3)
 factor = 1.d0/a_range
 !factor2 = factor*factor
 factor2 = 1.d0/(a1*a1+a2*a2+a3*a3)
 do i_gauss=1,n_gauss
    p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
 end do
 do i_gauss=1,n_gauss
    w_gauss(i_gauss) = factor*w_gauss(i_gauss)
 end do

 karray(:,:,:) = 0.d0
 !Use in this order (better for accuracy).
 loop_gauss: do i_gauss=n_gauss,1,-1
    !Gaussian
    pgauss = p_gauss(i_gauss)

    !We calculate the number of iterations to go from pgauss to p0_ref
    n_iter = nint((log(pgauss) - log(p0_cell))/log(4.d0))
    if (n_iter <= 0)then
       n_iter = 0
       p0gauss = pgauss
    else
       p0gauss = pgauss/4.d0**n_iter
    end if

    !Stupid integration
    !Do the integration with the exponential centered in i_kern
    kernel_scf(:) = 0.d0
    do i_kern=0,n_range
       kern = 0.d0
       do i=0,n_scf
          absci = x_scf(i) - real(i_kern,kind=8)*hgrid
          absci = absci*absci
          kern = kern + y_scf(i)*exp(-p0gauss*absci)*dx
       end do
       kernel_scf(i_kern) = kern
       kernel_scf(-i_kern) = kern
       if (abs(kern) < 1.d-18) then
          !Too small not useful to calculate
          exit
       end if
    end do

    !Start the iteration to go from p0gauss to pgauss
    call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,kern_1_scf)

    !Add to the kernel (only the local part)

    do i3=istart1,iend1
       i03 =  n2h - i3 + 1
       do i2=1,n02
          i02 = i2-1
          do i1=1,n01
             i01 = i1-1
             karray(i1+n1h,i2+n3h,i3-istart+1) = karray(i1+n1h,i2+n3h,i3-istart+1) + w_gauss(i_gauss)* &
                  kernel_scf(i01)*kernel_scf(i02)*kernel_scf(i03)
          end do
       end do
    end do
    do i3=istart2,iend2
       i03 = i3 - n2h -1
       do i2=1,n02
          i02 = i2-1
          do i1=1,n01
             i01 = i1-1
             karray(i1+n1h,i2+n3h,i3-istart+1) = karray(i1+n1h,i2+n3h,i3-istart+1) + w_gauss(i_gauss)* &
                  kernel_scf(i01)*kernel_scf(i02)*kernel_scf(i03)
          end do
       end do
    end do


 end do loop_gauss

 !Build the kernel in the real space as an even function, thus having a real FFT

 do i3=istart1,iend2
    do i2=1,n02
       do i1=2,n01
          karray(n1h+2-i1,i2+n3h,i3-istart+1) = karray(i1+n1h,i2+n3h,i3-istart+1)
       end do
    end do
    do i2=2,n02
       do i1=1,nker1
          karray(i1,n3h+2-i2,i3-istart+1) = karray(i1,i2+n3h,i3-istart+1)
       end do
    end do
 end do


 !De-allocations
 deallocate(kernel_scf)
 deallocate(kern_1_scf)
 deallocate(x_scf)
 deallocate(y_scf)

!!!!END KERNEL CONSTRUCTION



 allocate(karrayfour(2,nker1,nker2,nker3/nproc),stat=i_stat)
 if (i_stat /= 0) then
    print *,"Build_Kernel: Problem of memory allocation (kernel FFT)"
    stop
 end if

 if(iproc .eq. 0) print *,"Do a 3D PHalFFT for the kernel"

 call kernelfft(nfft1,nfft2,nfft3,nker1,nker2,nker3,nproc,iproc,karray,karrayfour)

!!$  amax=0.d0
!!$  do i1=1,nker1*nker2*nker3/nproc
!!$     amax=max(amax,abs(karrayfour(2,i1)))
!!$  end do
!!$  print *,"reality in the process",iproc,"up to",amax

 !Reconstruct the real kernel FFT
 do i3=1,n3k/nproc
    do i2=1,n2k
       do i1=1,n1k
          karrayoutLOC(i1,i2,i3)=karrayfour(1,i1,i2,i3)
       end do
    end do
 end do

 !De-allocations
 deallocate(karray,karrayfour)

end subroutine ParBuild_Kernel
!!***
