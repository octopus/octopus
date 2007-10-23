!  Copyright (C) Luigi Genovese, Thierry Deutsch, CEA Grenoble, 2006
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .




!!****h* BigDFT/build_kernel
!! NAME
!!   build_kernel
!!
!! FUNCTION
!!    Build the kernel of a gaussian function 
!!    for interpolating scaling functions.
!!
!! SYNOPSIS
!!    Build the kernel (karrayout) of a gaussian function 
!!    for interpolating scaling functions
!!    $$ K(j) = \int \int \phi(x) g(x'-x) \delta(x'- j) dx dx' $$
!!
!!    n01,n02,n03     Mesh dimensions of the density
!!    n1k,n2k,n3k     Dimensions of the kernel
!!    hgrid           Mesh step
!!    itype_scf       Order of the scaling function (8,14,16)
!!
!! AUTHORS
!!    T. Deutsch, L. Genovese
!! COPYRIGHT
!!    Copyright (C) 2005 CEA
!! CREATION DATE
!!    13/07/2005
!!
!! MODIFICATION HISTORY
!!    13/09/2005 Use hgrid instead of acell
!!    09/12/2005 Real kernel, stocked only half components
!!    13/12/2005 Add routines to simplify the port into Stefan's code
!!
!! SOURCE
!!
subroutine build_kernel(n01,n02,n03,nfft1,nfft2,nfft3,hgrid,itype_scf,karrayout)
  
  implicit none
  
  !Arguments
  integer, intent(in) :: n01,n02,n03,nfft1,nfft2,nfft3,itype_scf
  real(kind=8), intent(in) :: hgrid
  real(kind=8), dimension(nfft1/2+1,nfft2/2+1,nfft3/2+1), intent(out) :: karrayout
  
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
  real(kind=8), dimension(:,:,:), allocatable :: karrayhalf
  
  real(kind=8) :: ur_gauss,dr_gauss,acc_gauss,pgauss,kern,a_range
  real(kind=8) :: factor,factor2,dx,absci,p0gauss,p0_cell
  real(kind=8) :: a1,a2,a3
  integer :: nd1,nd2,nd3,n1k,n2k,n3k,n_scf
  integer :: i_gauss,n_range,n_cell
  integer :: i,n_iter,i1,i2,i3,i_kern,i_stat,i_allocated
  integer :: i01,i02,i03,inkee,n1h,n2h,n3h,nd1h
  
  !Number of integration points : 2*itype_scf*n_points
  n_scf=2*itype_scf*n_points
  !Dimensions of Kernel
  n1k=nfft1/2+1 ; n2k=nfft2/2+1 ; n3k=nfft3/2+1
  n1h=nfft1/2  ; n2h=nfft2/2 ; n3h=nfft3/2
  nd1 = nfft1 + modulo(nfft1+1,2)
  nd2 = nfft2 + modulo(nfft2+1,2)
  nd3 = nfft3 + modulo(nfft3+1,2)
  
  !Half size for the half FFT
  nd1h=(nd1+1)/2
  
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
  
  !Initialisation of the gaussian (Beylkin)
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
  
  karrayout(:,:,:) = 0.d0
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

     !Add to the kernel.
     do i3=1,n03
        i03 = i3-1
        do i2=1,n02
           i02 = i2-1
           do i1=1,n01
              i01 = i1-1
              karrayout(i1,i2,i3) = karrayout(i1,i2,i3) + w_gauss(i_gauss)* & 
                   kernel_scf(i01)*kernel_scf(i02)*kernel_scf(i03)
           end do
        end do
     end do

  end do loop_gauss
  !De-allocations
  deallocate(kernel_scf)
  deallocate(kern_1_scf)
  deallocate(x_scf)
  deallocate(y_scf)

  !Set karray
  allocate(karrayhalf(2,nd1h*nd2*nd3,2),stat=i_stat)
  if (i_stat /= 0) then
     print *,"Build_Kernel: Problem of memory allocation (karrayhalf)"
     stop
  end if
  !Set karray : use mirror symmetries
  inkee=1
  call karrayhalf_in(n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,&
       karrayout,karrayhalf)
  call fft(n1h,nfft2,nfft3,nd1h,nd2,nd3,karrayhalf,1,inkee)
  !Reconstruct the real kernel
  call kernel_recon(n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,&
       karrayhalf(1,1,inkee),karrayout)
  !De-allocations
  deallocate(karrayhalf)
end subroutine Build_Kernel
!!***


!!****h* BigDFT/calculate_dimensions
!! NAME
!!   calculate_dimensions
!!
!! FUNCTION
!!   Give the dimensions of the FFT
!!
!! SOURCE
!!
subroutine calculate_dimensions(n01,n02,n03,nfft1,nfft2,nfft3)
  implicit none
  !Arguments
  integer, intent(in) :: n01,n02,n03
  integer, intent(out) :: nfft1,nfft2,nfft3
  !Local variables
  integer :: i1,i2,i3,l1
  !Test 2*n01, 2*n02, 2*n03
  !write(*,*) 'in dimensions_fft',n01,n02,n03
  i1=2*n01
  i2=2*n02
  i3=2*n03
  do
     call fourier_dim(i1,nfft1)
     call fourier_dim(nfft1/2,l1)
     if (modulo(nfft1,2) == 0 .and. modulo(l1,2) == 0 .and. 2*l1 == nfft1) then
        exit
     end if
     i1=i1+1
  end do
  do
     call fourier_dim(i2,nfft2)
     if (modulo(nfft2,2) == 0) then
        exit
     end if
     i2=i2+1
  end do
  do
     call fourier_dim(i3,nfft3)
     if (modulo(nfft3,2) == 0) then
        exit
     end if
     i3=i3+1
  end do
  !nd1 = nfft1 + modulo(nfft1+1,2)
  !nd2 = nfft2 + modulo(nfft2+1,2)
  !nd3 = nfft3 + modulo(nfft3+1,2)
  !write(*,*) 'out dimensions_fft',nfft1,nfft2,nfft3
end subroutine calculate_dimensions
!!***


!!****h* BigDFT/karrayhalf_in
!! NAME
!!   karrayhalf_in
!!
!! FUNCTION
!!    Put in the array for4446666444 FFT
!!
!! SOURCE
!!
subroutine karrayhalf_in(n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,&
     kernel,karrayhalf)
  implicit none
  !Arguments
  integer, intent(in) :: n01,n02,n03,n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3
  real(kind=8), dimension(n1k,n2k,n3k), intent(in) :: kernel
  real(kind=8), dimension(2,(nd1+1)/2,nd2,nd3), intent(out) :: karrayhalf
  !Local variables
  real(kind=8), dimension(:), allocatable :: karray
  integer :: i1,i2,i3,i_stat,nd1h,n1h,n2h,n3h
  !Body
  n1h=nfft1/2
  n2h=nfft2/2
  n3h=nfft3/2
  allocate(karray(nfft1),stat=i_stat)
  if (i_stat /= 0) then
     print *,"Problem of memory allocation"
     stop
  end if
  nd1h=(nd1+1)/2
  karrayhalf(:,:,:,:) = 0.d0
  do i3=1,n03
     do i2=1,n02
        karray(:) = 0.d0
        do i1=1,n01
           karray(i1+n1h) = kernel(i1,i2,i3)
        end do
        do i1=2,n01
           karray(n1h-i1+1+nd1-nfft1) = kernel(i1,i2,i3)
        end do
        do i1=1,n1h
           karrayhalf(1,i1,i2+n2h,i3+n3h) = karray(2*i1-1)
           karrayhalf(2,i1,i2+n2h,i3+n3h) = karray(2*i1)
        end do
     end do
     do i2=2,n02
        do i1=1,nd1h
           karrayhalf(:,i1,n2h-i2+1+nd2-nfft2,i3+n3h) = &
                karrayhalf(:,i1,i2+n2h,i3+n3h)
        end do
     end do
  end do
  do i3=2,n03
     do i2=1,nd2
        do i1=1,nd1h
           karrayhalf(:,i1,i2,n3h-i3+1+nd3-nfft3) = karrayhalf(:,i1,i2,i3+n3h)
        end do
     end do
  end do
  !De-allocation
  deallocate(karray)
end subroutine karrayhalf_in
!!***


!!****h* BigDFT/kernel_recon
!! NAME
!!   kernel_recon
!!
!! FUNCTION
!!    Reconstruction of the kernel from the FFT array zarray
!!    We keep only the half kernel in each direction (x,y,z).
!!
!! SOURCE
!!
subroutine kernel_recon(n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3,zarray,karray)
  implicit none
  !Arguments
  integer, intent(in) :: n1k,n2k,n3k,nfft1,nfft2,nfft3,nd1,nd2,nd3
  real(kind=8), dimension(2,(nd1+1)/2*nd2*nd3), intent(in) :: zarray
  real(kind=8), dimension(n1k,n2k,n3k), intent(out) :: karray
  !Local variables
  real(kind=8), dimension(:), allocatable :: cos_array,sin_array
  integer :: i1,i2,i3,ind1,ind2,nd1h,n1h,n2h,n3h
  real(kind=8) :: rfe,ife,rfo,ifo,cp,sp,rk,ik,a,b,c,d,pi2
  !Body
  n1h=nfft1/2
  n2h=nfft2/2
  n3h=nfft3/2
  nd1h=(nd1+1)/2
  pi2=8.d0*datan(1.d0)
  pi2=pi2/real(nfft1,kind=8)
  allocate(cos_array(1:nd1h))
  allocate(sin_array(1:nd1h))
  do i1=1,nd1h
     cos_array(i1)= dcos(pi2*(i1-1))
     sin_array(i1)=-dsin(pi2*(i1-1))
  end do
  do i3=1,n3h+1
     do i2=1,n2h+1
        do i1=1,nd1h
           call norm_ind(nd1h,nd2,nd3,i1,i2,i3,ind1)
           call symm_ind(nd1h,nd2,nd3,i1,i2,i3,ind2)
           a=zarray(1,ind1)
           b=zarray(2,ind1)
           c=zarray(1,ind2)
           d=zarray(2,ind2)
           rfe=0.5d0*(a+c)
           ife=0.5d0*(b-d)
           rfo=0.5d0*(a-c)
           ifo=0.5d0*(b+d) 
           cp=cos_array(i1)
           sp=sin_array(i1)
           rk=rfe+cp*ifo-sp*rfo
           ik=ife-cp*rfo-sp*ifo
           !For big dimension 1.d-9 otherwise 1.d-10
           !Remove the test
           !if(abs(ik) >= 1.d-10) then
           !   print *,"non real kernel FFT",i1,i2,i3,ik  
           !   stop
           !end if
           !Build the intermediate FFT convolution (full)
           !call norm_ind(nd1,nd2,nd3,i1,i2,i3,indA)
           karray(i1,i2,i3)=rk 
        end do
     end do
  end do
  !De-allocations
  deallocate(cos_array)
  deallocate(sin_array)
end subroutine kernel_recon
!!***
