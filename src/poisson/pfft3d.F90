! FFT PART RELATED TO THE CONVOLUTION--------------------------------------------

#include "global.h"

integer function ncache_optimal()
  ncache_optimal = 8*1024
end function ncache_optimal


!!!HERE POT MUST BE THE KERNEL (BEWARE THE HALF DIMENSION)

!!****h* BigDFT/convolxc_off
!! NAME
!!   convolxc_off
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Applies the local FFT space Kernel to the density in Real space.
!!     Does NOT calculate the LDA exchange-correlation terms
!!
!! SYNOPSIS
!!     zf:          Density (input/output)
!!                  ZF(i1,i3,i2)
!!                  i1=1,md1 , i2=1,md2/nproc , i3=1,md3 
!!     pot:         Kernel, only the distributed part (REAL)
!!                  POT(i1,i2,i3)
!!                  i1=1,nd1 , i2=1,nd2 , i3=1,nd3/nproc 
!!     nproc:       number of processors used as returned by MPI_COMM_SIZE
!!     iproc:       [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3:    logical dimension of the transform. As transform lengths 
!!	            most products of the prime factors 2,3,5 are allowed.
!!                  The detailed table with allowed transform lengths can 
!!                  be found in subroutine CTRIG
!!     md1,md2,md3: Dimension of ZF
!!     nd1,nd2,nd3: Dimension of POT
!!     scal:        factor of renormalization of the FFT in order to acheve unitarity 
!!                  and the correct dimension
!!     hgrid:       grid spacing, used for integrating eharthree
!!     comm:        MPI communicator to use
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine convolxc_off(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3, &
  nproc,iproc,pot,zf,scal,hgrid,comm)
#if defined(MPI_MOD)
  use mpi
#endif

#if defined(MPI_H)
include "mpif.h"
#endif

  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc
  real(kind=8), intent(in) :: scal,hgrid
  real(kind=8), dimension(nd1,nd2,nd3/nproc), intent(in) :: pot
  real(kind=8), dimension(md1,md3,md2/nproc), intent(inout) :: zf
  integer, intent(in) :: comm

#if defined(HAVE_MPI)
  !Local variables
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  integer :: j2,j3,i1,i3,i,j,inzee,ierr
  real(kind=8) :: twopion !!! ,ehartreetmp
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: btrig1,btrig2,btrig3, &
       ftrig1,ftrig2,ftrig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3

  integer :: ncache_optimal

  !Body
  ! check input
  if (mod(n1,2).ne.0) stop 'Parallel convolution:ERROR:n1'
  if (mod(n2,2).ne.0) stop 'Parallel convolution:ERROR:n2'
  if (mod(n3,2).ne.0) stop 'Parallel convolution:ERROR:n3'
  if (nd1.lt.n1/2+1) stop 'Parallel convolution:ERROR:nd1'
  if (nd2.lt.n2/2+1) stop 'Parallel convolution:ERROR:nd2'
  if (nd3.lt.n3/2+1) stop 'Parallel convolution:ERROR:nd3'
  if (md1.lt.n1/2) stop 'Parallel convolution:ERROR:md1'
  if (md2.lt.n2/2) stop 'Parallel convolution:ERROR:md2'
  if (md3.lt.n3/2) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,nproc).ne.0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,nproc).ne.0) stop 'Parallel convolution:ERROR:md2'
  
  !defining work arrays dimensions
  ncache=ncache_optimal()
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4

  ! if (timing_flag == 1 .and. iproc ==0) print *,'parallel ncache=',ncache

  lzt=n2/2
  if (mod(n2/2,2).eq.0) lzt=lzt+1
  if (mod(n2/2,4).eq.0) lzt=lzt+1
  
  !Allocations
  allocate(btrig1(2,8192),ftrig1(2,8192),after1(7),now1(7),before1(7), &
       btrig2(2,8192),ftrig2(2,8192),after2(7),now2(7),before2(7), &
       btrig3(2,8192),ftrig3(2,8192),after3(7),now3(7),before3(7), &
       zw(2,ncache/4,2),zt(2,lzt,n1),zmpi2(2,n1,md2/nproc,nd3),cosinarr(2,n3/2))
  if (nproc.gt.1) allocate(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc))

  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,btrig3,after3,before3,now3,1,ic3)
  call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
  call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
  do  j=1,n1
     ftrig1(1,j)= btrig1(1,j)
     ftrig1(2,j)=-btrig1(2,j)
  enddo
  do  j=1,n2
     ftrig2(1,j)= btrig2(1,j)
     ftrig2(2,j)=-btrig2(2,j)
  enddo
  do  j=1,n3
     ftrig3(1,j)= btrig3(1,j)
     ftrig3(2,j)=-btrig3(2,j)
  enddo

  !Calculating array of phases for HalFFT decoding
  twopion=8.d0*datan(1.d0)/real(n3,kind=8)
  do i3=1,n3/2
     cosinarr(1,i3)=dcos(twopion*(i3-1))
     cosinarr(2,i3)=-dsin(twopion*(i3-1))
  end do

  !initializing integral
  !!! ehartree=0.d0

  ! transform along z axis
  lot=ncache/(2*n3)
  if (lot.lt.1) then  
     write(6,*) & 
          'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
          'least one 1-d FFT of this size even though this will' // & 
          'reduce the performance for shorter transform lengths'
     stop
  endif
  
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2/2) then
        do i1=1,(n1/2),lot
           ma=i1
           mb=min(i1+(lot-1),(n1/2))
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           call halfill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

           !performing FFT
           !input: I1,I3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result, 
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1/2,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     endif
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
     !communication scheduling
     call MPI_Alltoall(zmpi2(:, 1, 1, 1),n1*(md2/nproc)*(nd3/nproc), &
          MPI_DOUBLE_PRECISION, &
          zmpi1(:, 1, 1, 1, 1),n1*(md2/nproc)*(nd3/nproc), &
          MPI_DOUBLE_PRECISION,comm,ierr)
  endif
  !output: I1,J2,j3,Jp2,(jp3)

  !now each process perform complete convolution of its planes
  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
	Jp2stb=1
	J2stb=1
	Jp2stf=1
	J2stf=1
        
        ! transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif
        
        do j=1,n2/2,lot
           ma=j
           mb=min(j+(lot-1),n2/2)
           nfft=mb-ma+1

           !reverse index ordering, leaving the planes to be transformed at the end
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)
           
           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo

           !storing the last step into zt array
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                btrig1,after1(i),now1(i),before1(i),1)           
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering 
           !input: I2,i1,j3,(jp3)
           call switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)
           
           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !output: i1,i2,j3,(jp3)
           
           !Multiply with kernel in fourier space
           call multkernel(nd1,nd2,n1,n2,lot,nfft,j,pot(1,1,j3),zw(1,1,inzee))
           
           !TRANSFORM BACK IN REAL SPACE
           
           !transform along y axis
           !input: i1,i2,j3,(jp3)
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           enddo

           !reverse ordering
           !input: i1,I2,j3,(jp3)
           call unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           !output: I2,i1,j3,(jp3)
        end do
        
        !transform along x axis
        !input: I2,i1,j3,(jp3)
        lot=ncache/(4*n1)
        do j=1,n2/2,lot
           ma=j
           mb=min(j+(lot-1),n2/2)
           nfft=mb-ma+1

           !performing FFT
           i=1
           call fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                ftrig1,after1(i),now1(i),before1(i),-1)
           
           inzee=1
           do i=2,ic1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig1,after1(i),now1(i),before1(i),-1)
              inzee=3-inzee
           enddo
           !output: I2,I1,j3,(jp3)

           !reverse ordering
           !input: J2,Jp2,I1,j3,(jp3)
           if (nproc.eq.1) then
              call unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
           else
              call unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
           endif
           ! output: I1,J2,j3,Jp2,(jp3)
        end do
     endif
  end do
  
  
  !Interprocessor data transposition
  !input: I1,J2,j3,Jp2,(jp3)
  if (nproc.gt.1) then
     !communication scheduling
     call MPI_ALLTOALL(zmpi1(:, 1, 1, 1, 1),n1*(md2/nproc)*(nd3/nproc), &
          MPI_DOUBLE_PRECISION, &
          zmpi2(:, 1, 1, 1),n1*(md2/nproc)*(nd3/nproc), &
          MPI_DOUBLE_PRECISION,comm,ierr)
     !output: I1,J2,j3,jp3,(Jp2)
  endif

  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(2*n3)
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2/2) then
        do i1=1,(n1/2),lot
           ma=i1
           mb=min(i1+(lot-1),(n1/2))
           nfft=mb-ma+1

           !reverse ordering and repack the FFT data in order to be backward HalFFT transformed
           !input: I1,J2,i3,(Jp2)
           call unscramble_pack(i1,j2,lot,nfft,n1/2,n3,md2,nproc,nd3,zmpi2,zw(1,1,1),cosinarr)
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)           
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !calculates the Hartree energy locally and rebuild the output array
           call unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2)&
                ,scal) !!! ,ehartreetmp)

           !integrate local pieces together
           !!! ehartree=ehartree+0.5d0*ehartreetmp*hgrid**3
        end do
     endif
  end do
  
  !De-allocations  
  deallocate(btrig1,ftrig1,after1,now1,before1, &
       btrig2,ftrig2,after2,now2,before2, &
       btrig3,ftrig3,after3,now3,before3, &
       zmpi2,zw,zt,cosinarr)
  if (nproc.gt.1) deallocate(zmpi1)

#endif  
end subroutine convolxc_off


!!****h* BigDFT/multkernel
!! NAME
!!   multkernel
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Multiply with the kernel taking into account its symmetry
!!     Conceived to be used into convolution loops
!!
!! SYNOPSIS
!!     pot:      Kernel, symmetric and real, half the length
!!     zw:       Work array (input/output)
!!     n1,n2:    logical dimension of the FFT transform, reference for zw
!!     nd1,nd2:  Dimensions of POT
!!     jS, nfft: starting point of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine multkernel(nd1,nd2,n1,n2,lot,nfft,jS,pot,zw)
  !Argments
  integer, intent(in) :: nd1,nd2,n1,n2,lot,nfft,jS
  real(kind=8), dimension(nd1,nd2), intent(in) :: pot
  real(kind=8), dimension(2,lot,n2), intent(inout) :: zw
  !Local variables
  integer :: j,j1,i2,j2

  !Body
  
  !case i2=1
  do j=1,nfft
     j1=n1/2+1-abs(n1/2+2-jS-j)!this stands for j1=min(jS-1+j,n1+3-jS-j)
     zw(1,j,1)=zw(1,j,1)*pot(j1,1)
     zw(2,j,1)=zw(2,j,1)*pot(j1,1)
  end do

  !generic case
  do i2=2,n2/2
     do j=1,nfft
        j1=n1/2+1-abs(n1/2+2-jS-j)
        j2=n2+2-i2
        zw(1,j,i2)=zw(1,j,i2)*pot(j1,i2)
        zw(2,j,i2)=zw(2,j,i2)*pot(j1,i2)
        zw(1,j,j2)=zw(1,j,j2)*pot(j1,i2)
        zw(2,j,j2)=zw(2,j,j2)*pot(j1,i2)
     end do
  end do
  
  !case i2=n2/2+1
  do j=1,nfft
     j1=n1/2+1-abs(n1/2+2-jS-j)
     j2=n2/2+1
     zw(1,j,j2)=zw(1,j,j2)*pot(j1,j2)
     zw(2,j,j2)=zw(2,j,j2)*pot(j1,j2)
  end do
  
end subroutine multkernel


 	subroutine switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
	implicit real(8) (a-h,o-z)
        integer :: lot, n1, n2, lzt, j, nfft, i
	dimension zw(2,lot,n2),zt(2,lzt,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n2 is multiple of 2

! Low frequencies 
	do 100,j=1,nfft
        do 100,i=n2/2+1,n2
	zw(1,j,i)=zt(1,i-n2/2,j)
	zw(2,j,i)=zt(2,i-n2/2,j)
100	continue

! High frequencies 
        do 90,i=1,n2/2
	do 90,j=1,nfft
	zw(1,j,i)=0.d0
	zw(2,j,i)=0.d0
90	continue

	return
	end subroutine switch_upcorn

        
 	subroutine mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
	implicit real(8) (a-h,o-z)
        integer :: n1, nd2, nd3, nproc, lot, mfft, jp2, jp2stb, j2
        integer :: j2stb, nfft, i1, j3, md2
        
        dimension zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n1 is multiple of 2

	mfft=0
	do 300,Jp2=Jp2stb,nproc
	do 200,J2=J2stb,md2/nproc
	mfft=mfft+1
	if (mfft.gt.nfft) then
	Jp2stb=Jp2
	J2stb=J2
	return
	endif
        do 90,I1=1,n1/2
	zw(1,mfft,I1)=0.d0
	zw(2,mfft,I1)=0.d0
90	continue
        do 100,I1=n1/2+1,n1
	zw(1,mfft,I1)=zmpi1(1,I1-n1/2,J2,j3,Jp2)
	zw(2,mfft,I1)=zmpi1(2,I1-n1/2,J2,j3,Jp2)
100	continue
200	continue
	J2stb=1
300	continue
	end subroutine mpiswitch_upcorn

        subroutine halfill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
	implicit real(8) (a-h,o-z)
        integer :: lot, n3, md1, md3, i3, i1, nfft
	dimension zw(2,lot,n3/2),zf(md1,md3)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n3 is multiple of 4
!in principle we can relax this condition
	
        do 90,i3=1,n3/4
	do i1=1,nfft
	zw(1,i1,i3)=0.d0
	zw(2,i1,i3)=0.d0
	enddo
90      continue
        do 100,i3=n3/4+1,n3/2
	do i1=1,nfft
	zw(1,i1,i3)=zf(i1,2*i3-1-n3/2)
	zw(2,i1,i3)=zf(i1,2*i3-n3/2)
	enddo
100     continue


	return
	end subroutine halfill_upcorn


!!****h* BigDFT/scramble_unpack
!! NAME
!!   scramble_unpack
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Assign the correct planes to the work array zmpi2
!!     in order to prepare for interprocessor data transposition.
!!     In the meanwhile, it unpacks the data of the HalFFT in order to prepare for
!!     multiplication with the kernel
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (output)
!!     zw:             Work array (input)
!!     cosinarr:      Array of the phases needed for unpacking
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(2,n3/2), intent(in) :: cosinarr
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(out) :: zmpi2
  !Local variables
  integer :: i3,i,ind1,ind2
  real(kind=8) ::  a,b,c,d,cp,sp,feR,feI,foR,foI,fR,fI
  
  !Body

  !case i3=1 and i3=n3/2+1
  do i=0,nfft-1
     a=zw(1,i+1,1)
     b=zw(2,i+1,1)
     zmpi2(1,i1+i,j2,1)=a+b
     zmpi2(2,i1+i,j2,1)=0.d0
     zmpi2(1,i1+i,j2,n3/2+1)=a-b
     zmpi2(2,i1+i,j2,n3/2+1)=0.d0
  end do
  !case 2<=i3<=n3/2
  do i3=2,n3/2
     ind1=i3
     ind2=n3/2-i3+2
     cp=cosinarr(1,i3)
     sp=cosinarr(2,i3)
     do i=0,nfft-1
        a=zw(1,i+1,ind1)
        b=zw(2,i+1,ind1)
        c=zw(1,i+1,ind2)
        d=zw(2,i+1,ind2)
        feR=.5d0*(a+c)
        feI=.5d0*(b-d)
        foR=.5d0*(a-c)
        foI=.5d0*(b+d) 
        fR=feR+cp*foI-sp*foR
        fI=feI-cp*foR-sp*foI
        zmpi2(1,i1+i,j2,ind1)=fR
        zmpi2(2,i1+i,j2,ind1)=fI
     end do
  end do

end subroutine scramble_unpack

 
!!****h* BigDFT/unscramble_pack
!! NAME
!!   unscramble_pack
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Insert the correct planes of the work array zmpi2
!!     in order to prepare for backward FFT transform
!!     In the meanwhile, it packs the data in order to be transformed with the HalFFT 
!!     procedure
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (input)
!!     zw:             Work array (output)
!!     cosinarr:       Array of the phases needed for packing
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine unscramble_pack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3/2), intent(out) :: zw
  real(kind=8), dimension(2,n3/2), intent(in) :: cosinarr
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(in) :: zmpi2
  !Local variables
  integer :: i3,i,indA,indB
  real(kind=8) ::  a,b,c,d,cp,sp,re,ie,ro,io,rh,ih

  !Body

  do i3=1,n3/2
     indA=i3
     indB=n3/2+2-i3
     cp=cosinarr(1,i3)
     sp=cosinarr(2,i3)
     do i=0,nfft-1
        a=zmpi2(1,i1+i,j2,indA)
        b=zmpi2(2,i1+i,j2,indA)
        c= zmpi2(1,i1+i,j2,indB)
        d=-zmpi2(2,i1+i,j2,indB)
        re=(a+c)
        ie=(b+d)
        ro=(a-c)*cp-(b-d)*sp
        io=(a-c)*sp+(b-d)*cp
        rh=re-io 
        ih=ie+ro
        zw(1,i+1,indA)=rh
        zw(2,i+1,indA)=ih
     end do
  end do

end subroutine unscramble_pack




 	subroutine unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
	implicit real(8) (a-h,o-z)
        integer :: lot, n2, lzt, n1, j, nfft, i
	dimension zw(2,lot,n2),zt(2,lzt,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n2 is multiple of 2

! Low frequencies
	do 100,j=1,nfft
        do 100,i=1,n2/2
	zt(1,i,j)=zw(1,j,i)
	zt(2,i,j)=zw(2,j,i)
100	continue
	return
	end subroutine unswitch_downcorn


 	subroutine unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
	implicit real(8) (a-h,o-z)
        integer :: n1, md2, nproc, nd3, lot, mfft, i1, j3
        integer :: jp2, j2, nfft, jp2stf, j2stf
        dimension zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n1 is multiple of 2

	mfft=0
	do 300,Jp2=Jp2stf,nproc
	do 200,J2=J2stf,md2/nproc
	mfft=mfft+1
	if (mfft.gt.nfft) then
	Jp2stf=Jp2
	J2stf=J2
	return
	endif
        do 100,I1=1,n1/2
	zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
	zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
100	continue
200	continue
	J2stf=1
300	continue
	end subroutine unmpiswitch_downcorn


!!****h* BigDFT/unfill_downcorn
!! NAME
!!   unfill_downcorn
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Restore data into output array, calculating in the meanwhile
!!     Hartree energy of the potential 
!!
!! SYNOPSIS
!!     zf:          Original distributed density as well as
!!                  Distributed solution of the poisson equation (inout)
!!     zw:          FFT work array
!!     n3:          (twice the) dimension of the last FFTtransform.
!!     md1,md3:     Dimensions of the undistributed part of the real grid
!!     nfft:        number of planes
!!     scal:        Needed to achieve unitarity and correct dimensions
!!
!! WARNING
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf&
     ,scal) !!! ,ehartreetmp)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(md1,md3),intent(inout) :: zf
  real(kind=8), intent(in) :: scal
  !!! real(kind=8), intent(out) :: ehartreetmp
  !Local variables
  integer :: i3,i1
  real(kind=8) :: pot1

  !Body
  !!! ehartreetmp=0.d0
  do i3=1,n3/4
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        !!!ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3-1)
	zf(i1,2*i3-1)= pot1 
	pot1 = scal*zw(2,i1,i3)
        !!!ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3)
	zf(i1,2*i3)= pot1 
     enddo
  end do
  
end subroutine unfill_downcorn



! FFT PART RELATED TO THE KERNEL -----------------------------------------------------------------
!!****h* BigDFT/kernelfft
!! NAME
!!   kernelfft
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Calculates the FFT of the distributed kernel
!!
!! SYNOPSIS
!!     zf:          Real kernel (input)
!!                  zf(i1,i2,i3)
!!                  i1=1,nd1 , i2=1,nd2/nproc , i3=1,nd3 
!!     zr:          Distributed Kernel FFT 
!!                  zr(2,i1,i2,i3)
!!                  i1=1,nd1 , i2=1,nd2 , i3=1,nd3/nproc
!!     nproc:       number of processors used as returned by MPI_COMM_SIZE
!!     iproc:       [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3:    logical dimension of the transform. As transform lengths 
!!	            most products of the prime factors 2,3,5 are allowed.
!!                  The detailed table with allowed transform lengths can 
!!                  be found in subroutine CTRIG
!!     nd1,nd2,nd3: Dimensions of zr
!!     comm:        MPI communicator to use
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine kernelfft(n1,n2,n3,nd1,nd2,nd3,nproc,iproc,zf,zr,comm)
  use mpi_m

  implicit none

  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,nproc,iproc
  real(kind=8), dimension(nd1,n3,nd2/nproc), intent(in) :: zf
  real(kind=8), dimension(2,nd1,nd2,nd3/nproc), intent(out) :: zr
  integer, intent(in) :: comm

#if defined(HAVE_MPI)
  !Local variables
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2st,J2st
  integer :: j2,j3,i1,i3,i,j,inzee,ierr
  real(kind=8) :: twopion
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: trig1,trig2,trig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3

  integer ncache_optimal

  !Body

  !check input
  if (nd1.lt.n1) stop 'ERROR:nd1'
  if (nd2.lt.n2) stop 'ERROR:nd2'
  if (nd3.lt.n3/2+1) stop 'ERROR:nd3'
  if (mod(nd3,nproc).ne.0) stop 'ERROR:nd3'
  if (mod(nd2,nproc).ne.0) stop 'ERROR:nd2'
  
  !defining work arrays dimensions
  ncache = ncache_optimal()
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4
  lzt=n2
  if (mod(n2,2).eq.0) lzt=lzt+1
  if (mod(n2,4).eq.0) lzt=lzt+1
  
  !Allocations
  allocate(trig1(2,8192),after1(7),now1(7),before1(7), &
       trig2(2,8192),after2(7),now2(7),before2(7), &
       trig3(2,8192),after3(7),now3(7),before3(7), &
       zw(2,ncache/4,2),zt(2,lzt,n1), &
       zmpi2(2,n1,nd2/nproc,nd3),cosinarr(2,n3/2))
  if (nproc.gt.1) allocate(zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc))
  
  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,trig3,after3,before3,now3,1,ic3)
  call ctrig(n1,trig1,after1,before1,now1,1,ic1)
  call ctrig(n2,trig2,after2,before2,now2,1,ic2)
  
  !Calculating array of phases for HalFFT decoding
  twopion=8.d0*datan(1.d0)/real(n3,kind=8)
  do i3=1,n3/2
     cosinarr(1,i3)=dcos(twopion*(i3-1))
     cosinarr(2,i3)=-dsin(twopion*(i3-1))
  end do
  
  !transform along z axis

  lot=ncache/(2*n3)
  if (lot.lt.1) stop 'kernelfft:enlarge ncache for z'
  
  do j2=1,nd2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     !        if (iproc*(nd2/nproc)+j2.le.n2) then
     do i1=1,n1,lot
        ma=i1
        mb=min(i1+(lot-1),n1)
        nfft=mb-ma+1

        !inserting real data into complex array of half lenght
        !input: I1,I3,J2,(Jp2)
        call inserthalf(nd1,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))
        
        !performing FFT
        inzee=1
        do i=1,ic3
           call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                trig3,after3(i),now3(i),before3(i),1)
           inzee=3-inzee
        enddo
        !output: I1,i3,J2,(Jp2)

        !unpacking FFT in order to restore correct result, 
        !while exchanging components
        !input: I1,i3,J2,(Jp2)
        call scramble_unpack(i1,j2,lot,nfft,n1,n3,nd2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
        !output: I1,J2,i3,(Jp2)
     end do
!        endif
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
     !communication scheduling
     call MPI_ALLTOALL(zmpi2(:, 1, 1, 1),2*n1*(nd2/nproc)*(nd3/nproc), &
          MPI_DOUBLE_PRECISION, &
          zmpi1(:, 1, 1, 1, 1),2*n1*(nd2/nproc)*(nd3/nproc), &
          MPI_double_precision,comm,ierr)
     ! output: I1,J2,j3,Jp2,(jp3)
  endif


  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
        Jp2st=1
        J2st=1
        
        !transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) stop 'kernelfft:enlarge ncache for x'
        
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse ordering
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)          
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo
           !storing the last step into zt
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                trig1,after1(i),now1(i),before1(i),1)
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) stop 'kernelfft:enlarge ncache for y'

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering
           !input: I2,i1,j3,(jp3)
           call switch(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)

           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2-1
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !storing the last step into output array
           i=ic2
           call fftstp(lot,nfft,n2,nd1,nd2,zw(1,1,inzee),zr(1,j,1,j3), &
                trig2,after2(i),now2(i),before2(i),1)
           
        end do
        !output: i1,i2,j3,(jp3)
     endif
  end do

  !De-allocations
  deallocate(trig1,after1,now1,before1, &
       trig2,after2,now2,before2, &
       trig3,after3,now3,before3, &
       zmpi2,zw,zt,cosinarr)
  if (nproc.gt.1) deallocate(zmpi1)
#endif
end subroutine kernelfft


         subroutine switch(nfft,n2,lot,n1,lzt,zt,zw)
        implicit real(kind=8)        (a-h,o-z)
        integer :: lot, n2, lzt, n1, j, nfft, i
        dimension zw(2,lot,n2),zt(2,lzt,n1)

        do 200,j=1,nfft
        do 100,i=1,n2
        zw(1,j,i)=zt(1,i,j)
        zw(2,j,i)=zt(2,i,j)
100        continue
200        continue
        return
        end subroutine switch

         subroutine mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw)
        implicit real(kind=8)        (a-h,o-z)
        integer :: n1, nd2, nproc, nd3, lot, mfft, jp2, jp2st
        integer :: j2st, nfft, i1, j3, j2
        dimension zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc),zw(2,lot,n1)

        mfft=0
        do 300,Jp2=Jp2st,nproc
        do 200,J2=J2st,nd2/nproc
        mfft=mfft+1
        if (mfft.gt.nfft) then
        Jp2st=Jp2
        J2st=J2
        return
        endif
        do 100,I1=1,n1
        zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
        zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
100        continue
200        continue
        J2st=1
300        continue
        end subroutine mpiswitch

        subroutine inserthalf(nd1,lot,nfft,n3,zf,zw)
        implicit real(kind=8)        (a-h,o-z)
        integer :: lot, n3, nd1, i3, i1, nfft
        dimension zw(2,lot,n3/2),zf(nd1,n3)

        do 100,i3=1,n3/2
        do 100,i1=1,nfft
        zw(1,i1,i3)=zf(i1,2*i3-1)
        zw(2,i1,i3)=zf(i1,2*i3)
100     continue

        return
        end subroutine inserthalf
