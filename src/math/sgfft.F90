#include <global.h> 
! This routines are part of the ISF poisson solver, eventually they
! will be integrated with the other FFT. Do not use them for other
! purposes.

module sgfft_m
  use global_m
  use messages_m
  use mpi_m
#ifdef HAVE_OMP
  use omp
#endif
  use profiling_m

  implicit none
  
  private

  public ::           &
    fft,              &
    fourier_dim,      &
    kernelfft,        &
    convolxc_off
  
contains

  !!****h* BigDFT/fourier_dim
  !! NAME
  !!   fourier_dim
  !!
  !! FUNCTION
  !!   Give a number n_next > n compatible for the FFT
  !!
  !! SOURCE
  !!
  subroutine fourier_dim(n,n_next)
    integer, intent(in) :: n
    integer, intent(out) :: n_next

    !Local variables
    integer, parameter :: ndata1024 = 149, ndata = 149
    !Multiple of 2,3,5
    integer, dimension(ndata), parameter :: idata = (/   &
      3,    4,   5,     6,    8,    9,   12,   15,   16,   18, &
      20,   24,   25,   27,   30,   32,   36,   40,   45,   48, &
      54,   60,   64,   72,   75,   80,   81,   90,   96,  100, &
      108,  120,  125,  128,  135,  144,  150,  160,  162,  180, &
      192,  200,  216,  225,  240,  243,  256,  270,  288,  300, &
      320,  324,  360,  375,  384,  400,  405,  432,  450,  480, &
      486,  500,  512,  540,  576,  600,  625,  640,  648,  675, &
      720,  729,  750,  768,  800,  810,  864,  900,  960,  972, &
      1000, 1024, 1080, 1125, 1152, 1200, 1215, 1280, 1296, 1350,&
      1440, 1458, 1500, 1536, 1600, 1620, 1728, 1800, 1875, 1920,&
      1944, 2000, 2025, 2048, 2160, 2250, 2304, 2400, 2430, 2500,&
      2560, 2592, 2700, 2880, 3000, 3072, 3125, 3200, 3240, 3375,&
      3456, 3600, 3750, 3840, 3888, 4000, 4050, 4096, 4320, 4500,&
      4608, 4800, 5000, 5120, 5184, 5400, 5625, 5760, 6000, 6144,&
      6400, 6480, 6750, 6912, 7200, 7500, 7680, 8000, 8192 /)
    integer :: i

    loop_data: do i=1,ndata1024
      if (n <= idata(i)) then
        n_next = idata(i)
        return
      end if
    end do loop_data
    write(unit=*,fmt=*) "fourier_dim: ",n," is bigger than ",idata(ndata1024)
    stop
  end subroutine fourier_dim
  !!***


  !  Copyright (C) Stefan Goedecker, CEA Grenoble, 2002
  !  This file is distributed under the terms of the
  !  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
  
  
  ! --------------------------------------------------------------
  !   3-dimensional complex-complex FFT routine: 
  !   When compared to the best vendor implementations on RISC architectures 
  !   it gives close to optimal performance (perhaps loosing 20 percent in speed)
  !   and it is significanly faster than many not so good vendor implementations 
  !   as well as other portable FFT`s. 
  !   On all vector machines tested so far (Cray, NEC, Fujitsu) is 
  !   was significantly faster than the vendor routines
  ! The theoretical background is described in :
  ! 1) S. Goedecker: Rotating a three-dimensional array in optimal
  ! positions for vector processing: Case study for a three-dimensional Fast
  ! Fourier Transform, Comp. Phys. Commun. \underline{76}, 294 (1993)
  ! Citing of this reference is greatly appreciated if the routines are used 
  ! for scientific work.
  
  
  ! Presumably good compiler flags:
  ! IBM, serial power 2: xlf -qarch=pwr2 -O2 -qmaxmem=-1
  ! with OpenMP: IBM: xlf_r -qfree -O4 -qarch=pwr3 -qtune=pwr3 -qsmp=omp -qmaxmem=-1 ; 
  !                   a.out
  ! DEC: f90 -O3 -arch ev67 -pipeline
  ! with OpenMP: DEC: f90 -O3 -arch ev67 -pipeline -omp -lelan ; 
  !                   prun -N1 -c4 a.out
  
  
  !-----------------------------------------------------------
  
  ! FFT PART -----------------------------------------------------------------
  !        CALCULATES THE DISCRETE FOURIER TRANSFORM F(I1,I2,I3)=
  !        S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) R(j1,j2,j3)
  !       with optimal performance on vector computer, workstations and 
  !       multiprocessor shared memory computers using OpenMP compiler directives
  !        INPUT:
  !            n1,n2,n3:physical dimension of the transform. It must be a 
  !                     product of the prime factors 2,3,5, but greater than 3. 
  !                    If two ni`s are equal it is recommended to place them 
  !                    behind each other.
  !            nd1,nd2,nd3:memory dimension of Z. ndi must always be greater or 
  !                        equal than ni. On a vector machine, it is recomended 
  !                       to chose ndi=ni if ni is odd and ndi=ni+1 if ni is 
  !                       even to obtain optimal execution speed. On RISC 
  !                       machines ndi=ni is usually fine for odd ni, for even 
  !                       ni one should try ndi=ni+1, ni+2, ni+4 to find the 
  !                       optimal performance. 
  !           inzee=1: first part of Z is data (input) array, 
  !                    second part work array
  !           inzee=2: first part of Z is work array, second part data array
  !                Z(1,i1,i2,i3,inzee)=real(R(i1,i2,i3))
  !                Z(2,i1,i2,i3,inzee)=imag(R(i1,i2,i3))
  !        OUTPUT:
  !           inzee=1: first part of Z is data (output) array, 
  !                    second part work array
  !           inzee=2: first part of Z is work array, second part data array
  !                real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
  !                imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
  !           inzee on output is in general different from inzee on input
  !        The input data are always overwritten independently of the 
  !       value of inzee.
  ! PERFORMANCE AND THE NCACHE
  !       The most important feature for performance is the right choice of 
  !       the parameter ncache. On a vector machine ncache has to be put to 0.
  !       On a RISC machine with cache, it is very important to find the optimal 
  !       value of NCACHE. NCACHE determines the size of the work array zw, that
  !       has to fit into cache. It has therefore to be chosen to equal roughly 
  !        half the size of the physical cache in units of real*8 numbers.
  !       If the machine has 2 cache levels it can not be predicted which 
  !       cache level will be the most relevant one for choosing ncache. 
  !       The optimal value of ncache can easily be determined by numerical 
  !       experimentation. A too large value of ncache leads to a dramatic 
  !       and sudden decrease of performance, a too small value to a to a 
  !       slow and less dramatic decrease of performance. If NCACHE is set 
  !       to a value so small, that not even a single one dimensional transform 
  !       can be done in the workarray zw, the program stops with an error 
  !       message.
  !  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
  !  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1995, 1999
  !  Copyright (C) Stefan Goedecker, CEA Grenoble, 2002
  !  This file is distributed under the terms of the
  !  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

  subroutine fft(n1, n2, n3, nd1, nd2, nd3, z, isign, inzee)
    integer, intent(in)    :: n1, n2, n3, nd1, nd2, nd3
    real(8), intent(inout) :: z(1:2, 1:nd1*nd2*nd3, 1:2)
    integer, intent(in)    :: isign
    integer, intent(inout) :: inzee

    real(8), allocatable :: zw(:, :, :)
    real(8), allocatable :: trig(:, :)
    integer, allocatable :: after(:), now(:), before(:)

    integer :: ncache, nfft, mm, i, ic
    integer :: npr, iam, inzet, m, lot, nn, lotomp, j, jompa, jompb
    integer :: n, ma, mb, jj, inzeep
    type(profile_t), save :: prof

    call profiling_in(prof, "SGFFT")

    if (max(n1,n2,n3).gt.8192) stop 'fft:8192 limit reached'

    ! some reasonable values of ncache: 
    !   IBM/RS6000/590: 16*1024 ; IBM/RS6000/390: 3*1024 ; 
    !   IBM/PwPC: 1*1024 ; SGI/MIPS/R8000: 16*1024 ; DEC/Alpha/EV5 and EV6 6*1024
    !   But if you care about performance find the optimal value of ncache yourself!
    !       On all vector machines: ncache=0

    ncache = 8192
    if (ncache /= 0 .and. ncache <= max(n1,n2,n3)*4) ncache=max(n1,n2,n3/2)*4

    ! check whether input values are reasonable
    if (inzee.le.0 .or. inzee.ge.3) stop 'fft:wrong inzee'
    if (isign.ne.1 .and. isign.ne.-1) stop 'fft:wrong isign'
    if (n1.gt.nd1) stop 'fft:n1>nd1'
    if (n2.gt.nd2) stop 'fft:n2>nd2'
    if (n3.gt.nd3) stop 'fft:n3>nd3'


    ! vector computer with memory banks:
    if (ncache.eq.0) then
      allocate(trig(2,8192),after(20),now(20),before(20))

      call ctrig(n3,trig,after,before,now,isign,ic)
      nfft=nd1*n2
      mm=nd1*nd2

      do i=1,ic-1
        call fftstp(mm,nfft,nd3,mm,nd3,z(1,1,inzee),z(1,1,3-inzee), &
          trig,after(i),now(i),before(i),isign)
        inzee=3-inzee
      end do

      i=ic

      call fftrot(mm,nfft,nd3,mm,nd3,z(1,1,inzee),z(1,1,3-inzee), &
        trig,after(i),now(i),before(i),isign)

      inzee=3-inzee

      if (n2.ne.n3) call ctrig(n2,trig,after,before,now,isign,ic)
      nfft=nd3*n1
      mm=nd3*nd1

      do i=1,ic-1
        call fftstp(mm,nfft,nd2,mm,nd2,z(1,1,inzee),z(1,1,3-inzee), &
          trig,after(i),now(i),before(i),isign)
        inzee=3-inzee
      end do

      i=ic
      call fftrot(mm,nfft,nd2,mm,nd2,z(1,1,inzee),z(1,1,3-inzee), &
        trig,after(i),now(i),before(i),isign)
      inzee=3-inzee

      if (n1.ne.n2) call ctrig(n1,trig,after,before,now,isign,ic)
      nfft=nd2*n3
      mm=nd2*nd3

      do i=1,ic-1
        call fftstp(mm,nfft,nd1,mm,nd1,z(1,1,inzee),z(1,1,3-inzee), &
          trig,after(i),now(i),before(i),isign)
        inzee=3-inzee
      end do

      i=ic

      call fftrot(mm,nfft,nd1,mm,nd1,z(1,1,inzee),z(1,1,3-inzee), &
        trig,after(i),now(i),before(i),isign)

      inzee=3-inzee

      ! RISC machine with cache:
    else
      ! INtel IFC does not understand default(private)
      
      !$omp parallel & 
      !$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
      !$omp shared(n1,n2,n3,nd1,nd2,nd3,z,isign,inzee,ncache) 
      npr = 1
      ! npr = omp_get_num_threads()
      iam = 0
      ! iam = omp_get_thread_num()

      allocate(zw(2,ncache/4,2),trig(2,1024),after(20),now(20),before(20))

      inzet=inzee
      ! TRANSFORM ALONG Z AXIS

      mm = nd1*nd2
      m = nd3
      lot = max(1,ncache/(4*n3))
      nn = lot
      n=n3
      if (2*n*lot*2.gt.ncache) stop 'fft:ncache1'

      call ctrig(n3,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd1*n2)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd1*n2)
        nfft=mb-ma+1
        j=ma
        jj=j*nd3-nd3+1
        call fftrot(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
          trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd1*n2)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd1*n2)
        do j=jompa,jompb,lot
          ma=j
          mb=min(j+(lot-1),jompb)
          nfft=mb-ma+1
          jj=j*nd3-nd3+1

          i=1
          inzeep=2
          call fftstp(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
            trig,after(i),now(i),before(i),isign)
          inzeep=1

          do i=2,ic-1
            call fftstp(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
              trig,after(i),now(i),before(i),isign)
            inzeep=3-inzeep
          end do
          i=ic
          call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
            trig,after(i),now(i),before(i),isign)
        end do
      endif

      inzet=3-inzet

      !$omp barrier

      ! TRANSFORM ALONG Y AXIS
      mm=nd3*nd1
      m=nd2
      lot=max(1,ncache/(4*n2))
      nn=lot
      n=n2
      if (2*n*lot*2.gt.ncache) stop 'fft:ncache2'

      if (n2.ne.n3) call ctrig(n2,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd3*n1)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd3*n1)
        nfft=mb-ma+1
        j=ma
        jj=j*nd2-nd2+1
        call fftrot(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
          trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd3*n1)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd3*n1)
        do j=jompa,jompb,lot
          ma=j
          mb=min(j+(lot-1),jompb)
          nfft=mb-ma+1
          jj=j*nd2-nd2+1

          i=1
          inzeep=2
          call fftstp(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
            trig,after(i),now(i),before(i),isign)
          inzeep=1

          do i=2,ic-1
            call fftstp(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
              trig,after(i),now(i),before(i),isign)
            inzeep=3-inzeep
          end do

          i=ic
          call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
            trig,after(i),now(i),before(i),isign)
        end do
      endif
      inzet=3-inzet

      !$omp barrier

      ! TRANSFORM ALONG X AXIS
      mm=nd2*nd3
      m=nd1
      lot=max(1,ncache/(4*n1))
      nn=lot
      n=n1
      if (2*n*lot*2.gt.ncache) stop 'fft:ncache3'

      if (n1.ne.n2) call ctrig(n1,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd2*n3)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd2*n3)
        nfft=mb-ma+1
        j=ma
        jj=j*nd1-nd1+1
        call fftrot(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
          trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd2*n3)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd2*n3)
        do j=jompa,jompb,lot
          ma=j
          mb=min(j+(lot-1),jompb)
          nfft=mb-ma+1
          jj=j*nd1-nd1+1

          i=1
          inzeep=2
          call fftstp(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
            trig,after(i),now(i),before(i),isign)
          inzeep=1

          do i=2,ic-1
            call fftstp(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
              trig,after(i),now(i),before(i),isign)
            inzeep=3-inzeep
          end do
          i=ic
          call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
            trig,after(i),now(i),before(i),isign)
        end do
      endif
      inzet=3-inzet

      deallocate(zw,trig,after,now,before)
      if (iam.eq.0) inzee=inzet
      !$omp end parallel  

    endif
    
    call profiling_out(prof)

  end subroutine fft

  ! ---------------------------------------------------------------------------------

  !  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
  !  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
  !  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
  !  This file is distributed under the terms of the
  !  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
  
  !     Different factorizations affect the performance
  !     Factoring 64 as 4*4*4 might for example be faster on some machines than 8*8.
  subroutine ctrig(n,trig,after,before,now,isign,ic)
    integer :: n, isign, ic
    integer :: now(7), after(7), before(7)
    real(8) :: trig(2,8192)

    integer :: i, j, itt, nh
    real(8) :: angle, trigc, trigs, twopi
    INTEGER, DIMENSION(7,149) :: idata 

    ! The factor 6 is only allowed in the first place!
    data ((idata(i,j),i=1,7),j=1,74) /                        &
      3,   3, 1, 1, 1, 1, 1,       4,   4, 1, 1, 1, 1, 1,   &
      5,   5, 1, 1, 1, 1, 1,       6,   6, 1, 1, 1, 1, 1,   &
      8,   8, 1, 1, 1, 1, 1,       9,   3, 3, 1, 1, 1, 1,   &
      12,   4, 3, 1, 1, 1, 1,      15,   5, 3, 1, 1, 1, 1,   &
      16,   4, 4, 1, 1, 1, 1,      18,   6, 3, 1, 1, 1, 1,   &
      20,   5, 4, 1, 1, 1, 1,      24,   8, 3, 1, 1, 1, 1,   &
      25,   5, 5, 1, 1, 1, 1,      27,   3, 3, 3, 1, 1, 1,   &
      30,   6, 5, 1, 1, 1, 1,      32,   8, 4, 1, 1, 1, 1,   &
      36,   4, 3, 3, 1, 1, 1,      40,   8, 5, 1, 1, 1, 1,   &
      45,   5, 3, 3, 1, 1, 1,      48,   4, 4, 3, 1, 1, 1,   &
      54,   6, 3, 3, 1, 1, 1,      60,   5, 4, 3, 1, 1, 1,   &
      64,   8, 8, 1, 1, 1, 1,      72,   8, 3, 3, 1, 1, 1,   &
      75,   5, 5, 3, 1, 1, 1,      80,   5, 4, 4, 1, 1, 1,   &
      81,   3, 3, 3, 3, 1, 1,      90,   6, 5, 3, 1, 1, 1,   &
      96,   8, 4, 3, 1, 1, 1,     100,   5, 5, 4, 1, 1, 1,   &
      108,   4, 3, 3, 3, 1, 1,     120,   8, 5, 3, 1, 1, 1,   &
      125,   5, 5, 5, 1, 1, 1,     128,   8, 4, 4, 1, 1, 1,   &
      135,   5, 3, 3, 3, 1, 1,     144,   6, 8, 3, 1, 1, 1,   &
      150,   6, 5, 5, 1, 1, 1,     160,   8, 5, 4, 1, 1, 1,   &
      162,   6, 3, 3, 3, 1, 1,     180,   5, 4, 3, 3, 1, 1,   &
      192,   6, 8, 4, 1, 1, 1,     200,   8, 5, 5, 1, 1, 1,   &
      216,   8, 3, 3, 3, 1, 1,     225,   5, 5, 3, 3, 1, 1,   &
      240,   6, 8, 5, 1, 1, 1,     243,   3, 3, 3, 3, 3, 1,   &
      256,   8, 8, 4, 1, 1, 1,     270,   6, 5, 3, 3, 1, 1,   &
      288,   8, 4, 3, 3, 1, 1,     300,   5, 5, 4, 3, 1, 1,   &
      320,   5, 4, 4, 4, 1, 1,     324,   4, 3, 3, 3, 3, 1,   &
      360,   8, 5, 3, 3, 1, 1,     375,   5, 5, 5, 3, 1, 1,   &
      384,   8, 4, 4, 3, 1, 1,     400,   5, 5, 4, 4, 1, 1,   &
      405,   5, 3, 3, 3, 3, 1,     432,   4, 4, 3, 3, 3, 1,   &
      450,   6, 5, 5, 3, 1, 1,     480,   8, 5, 4, 3, 1, 1,   &
      486,   6, 3, 3, 3, 3, 1,     500,   5, 5, 5, 4, 1, 1,   &
      512,   8, 8, 8, 1, 1, 1,     540,   5, 4, 3, 3, 3, 1,   &
      576,   4, 4, 4, 3, 3, 1,     600,   8, 5, 5, 3, 1, 1,   &
      625,   5, 5, 5, 5, 1, 1,     640,   8, 5, 4, 4, 1, 1,   &
      648,   8, 3, 3, 3, 3, 1,     675,   5, 5, 3, 3, 3, 1,   &
      720,   5, 4, 4, 3, 3, 1,     729,   3, 3, 3, 3, 3, 3,   &
      750,   6, 5, 5, 5, 1, 1,     768,   4, 4, 4, 4, 3, 1/

    data ((idata(i,j),i=1,7),j=75,149) /                        &
      800,   8, 5, 5, 4, 1, 1,     810,   6, 5, 3, 3, 3, 1,   &
      864,   8, 4, 3, 3, 3, 1,     900,   5, 5, 4, 3, 3, 1,   &
      960,   5, 4, 4, 4, 3, 1,     972,   4, 3, 3, 3, 3, 3,   &
      1000,   8, 5, 5, 5, 1, 1,    1024,   4, 4, 4, 4, 4, 1,   &
      1080,   6, 5, 4, 3, 3, 1,    1125,   5, 5, 5, 3, 3, 1,   &
      1152,   6, 4, 4, 4, 3, 1,    1200,   6, 8, 5, 5, 1, 1,   &
      1215,   5, 3, 3, 3, 3, 3,    1280,   8, 8, 5, 4, 1, 1,   &
      1296,   6, 8, 3, 3, 3, 1,    1350,   6, 5, 5, 3, 3, 1,   &
      1440,   6, 5, 4, 4, 3, 1,    1458,   6, 3, 3, 3, 3, 3,   &
      1500,   5, 5, 5, 4, 3, 1,    1536,   6, 8, 8, 4, 1, 1,   &
      1600,   8, 8, 5, 5, 1, 1,    1620,   5, 4, 3, 3, 3, 3,   &
      1728,   6, 8, 4, 3, 3, 1,    1800,   6, 5, 5, 4, 3, 1,   &
      1875,   5, 5, 5, 5, 3, 1,    1920,   6, 5, 4, 4, 4, 1,   &
      1944,   6, 4, 3, 3, 3, 3,    2000,   5, 5, 5, 4, 4, 1,   &
      2025,   5, 5, 3, 3, 3, 3,    2048,   8, 4, 4, 4, 4, 1,   &
      2160,   6, 8, 5, 3, 3, 1,    2250,   6, 5, 5, 5, 3, 1,   &
      2304,   6, 8, 4, 4, 3, 1,    2400,   6, 5, 5, 4, 4, 1,   &
      2430,   6, 5, 3, 3, 3, 3,    2500,   5, 5, 5, 5, 4, 1,   &
      2560,   8, 5, 4, 4, 4, 1,    2592,   6, 4, 4, 3, 3, 3,   &
      2700,   5, 5, 4, 3, 3, 3,    2880,   6, 8, 5, 4, 3, 1,   &
      3000,   6, 5, 5, 5, 4, 1,    3072,   6, 8, 4, 4, 4, 1,   &
      3125,   5, 5, 5, 5, 5, 1,    3200,   8, 5, 5, 4, 4, 1,   &
      3240,   6, 5, 4, 3, 3, 3,    3375,   5, 5, 5, 3, 3, 3,   &
      3456,   6, 4, 4, 4, 3, 3,    3600,   6, 8, 5, 5, 3, 1,   &
      3750,   6, 5, 5, 5, 5, 1,    3840,   6, 8, 5, 4, 4, 1,   &
      3888,   6, 8, 3, 3, 3, 3,    4000,   8, 5, 5, 5, 4, 1,   &
      4050,   6, 5, 5, 3, 3, 3,    4096,   8, 8, 4, 4, 4, 1,   &
      4320,   6, 5, 4, 4, 3, 3,    4500,   5, 5, 5, 4, 3, 3,   &
      4608,   6, 8, 8, 4, 3, 1,    4800,   6, 8, 5, 5, 4, 1,   &
      5000,   8, 5, 5, 5, 5, 1,    5120,   8, 8, 5, 4, 4, 1,   &
      5184,   6, 8, 4, 3, 3, 3,    5400,   6, 5, 5, 4, 3, 3,   &
      5625,   5, 5, 5, 5, 3, 3,    5760,   6, 8, 8, 5, 3, 1,   &
      6000,   6, 8, 5, 5, 5, 1,    6144,   6, 8, 8, 4, 4, 1,   &
      6400,   8, 8, 5, 5, 4, 1,    6480,   6, 8, 5, 3, 3, 3,   &
      6750,   6, 5, 5, 5, 3, 3,    6912,   6, 8, 4, 4, 3, 3,   &
      7200,   6, 5, 5, 4, 4, 3,    7500,   5, 5, 5, 5, 4, 3,   &
      7680,   6, 8, 8, 5, 4, 1,    8000,   8, 8, 5, 5, 5, 1,   &
      8192,   8, 8, 8, 4, 4, 1 /

    do i = 1, 149
      if (n.eq.idata(1,i)) then
        ic=0
        do j=1,6
          itt=idata(1+j,i)
          if (itt.gt.1) then
            ic=ic+1
            now(j)=idata(1+j,i)
          else
            goto 1000
          endif
        end do
        goto 1000
      end if
    end do

    print*,'VALUE OF',n,'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
37  format(15(i5))
    write(6,37) (idata(1,i),i=1,149)
    stop
1000 continue

    after(1)=1
    before(ic)=1
    do i = 2, ic
      after(i)=after(i-1)*now(i-1)
      before(ic-i+1)=before(ic-i+2)*now(ic-i+2)
    end do

    twopi=6.283185307179586d0
    angle=isign*twopi/n

    if (mod(n,2).eq.0) then
      nh=n/2
      trig(1,1)=1.0_8
      trig(2,1)=0.0_8
      trig(1,nh+1)=-1.0_8
      trig(2,nh+1)=0.0_8
      do i = 1, nh - 1
        trigc=cos(i*angle)
        trigs=sin(i*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
      end do
    else
      nh=(n-1)/2
      trig(1,1)=1.0_8
      trig(2,1)=0.0_8
      do i = 1, nh
        trigc=cos(i*angle)
        trigs=sin(i*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
      end do
    end if

  end subroutine ctrig


  ! ------------------------------------------------------------------------
  !  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
  !  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1995, 1999
  !  This file is distributed under the terms of the
  !  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
  
  subroutine fftstp(mm,nfft,m,nn,n,zin,zout,trig,after,now,before,isign)
    integer :: mm, nfft, m, nn, n, isign
    integer :: after, before, atn, atb, now
    real(8) :: trig(2,8192), zin(2, mm, m), zout(2, nn, n)

    real(8) :: rt2i, dp, cp, cm, ci5, cr5, ci6, cr6, am, ap, ci8, cr8
    real(8) :: r, r1, r2, r3, r4, r5, r6, r7, r8, r25, r34
    real(8) :: s, s1, s2, s3, s4, s5, s6, s7, s8, s25, s34
    real(8) :: bb, bm, dm, bp
    real(8) :: cr2, ci2, cr3, ci3, cr4, ci4, cr7, ci7
    real(8) :: cos2, sin2, cos4, sin4
    real(8) :: ur1, ur2, ur3, ui1, ui2, ui3
    real(8) :: vr1, vr2, vr3, vi1, vi2, vi3
    integer :: ia, ib, j, ias, itrig, itt
    integer :: nin1, nin2, nin3, nin4, nin5, nin6, nin7, nin8
    integer :: nout1, nout2, nout3, nout4, nout5, nout6, nout7, nout8

    atn=after*now
    atb=after*before

    !         sqrt(.5d0)
    rt2i=0.7071067811865475d0
    if (now.eq.2) then
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib = 1, before
        nin1=nin1+after
        nin2=nin1+atb
        nout1=nout1+atn
        nout2=nout1+after
        do j = 1, nfft
          r1=zin(1,j,nin1)
          s1=zin(2,j,nin1)
          r2=zin(1,j,nin2)
          s2=zin(2,j,nin2)
          zout(1,j,nout1)= r2 + r1
          zout(2,j,nout1)= s2 + s1
          zout(1,j,nout2)= r1 - r2
          zout(2,j,nout2)= s1 - s2
        end do
      end do
      do ia = 2, after
        ias=ia-1
        if (2*ias.eq.after) then
          if (isign.eq.1) then
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nout1=nout1+atn
              nout2=nout1+after
              do j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                zout(1,j,nout1)= r1 - r2
                zout(2,j,nout1)= s2 + s1
                zout(1,j,nout2)= r2 + r1
                zout(2,j,nout2)= s1 - s2
              end do
            end do
          else
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nout1=nout1+atn
              nout2=nout1+after
              do j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                zout(1,j,nout1)= r2 + r1
                zout(2,j,nout1)= s1 - s2
                zout(1,j,nout2)= r1 - r2
                zout(2,j,nout2)= s2 + s1
              end do
            end do
          endif
        else if (4*ias.eq.after) then
          if (isign.eq.1) then
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nout1=nout1+atn
              nout2=nout1+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                zout(1,j,nout1)= r2 + r1
                zout(2,j,nout1)= s2 + s1
                zout(1,j,nout2)= r1 - r2
                zout(2,j,nout2)= s1 - s2
              end do
            end do
          else
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nout1=nout1+atn
              nout2=nout1+after
              do j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(s - r)*rt2i
                zout(1,j,nout1)= r2 + r1
                zout(2,j,nout1)= s2 + s1
                zout(1,j,nout2)= r1 - r2
                zout(2,j,nout2)= s1 - s2
              end do
            end do
          endif
        else if (4*ias.eq.3*after) then
          if (isign.eq.1) then
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nout1=nout1+atn
              nout2=nout1+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(r - s)*rt2i
                zout(1,j,nout1)= r1 - r2
                zout(2,j,nout1)= s2 + s1
                zout(1,j,nout2)= r2 + r1
                zout(2,j,nout2)= s1 - s2
              end do
            end do
          else
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nout1=nout1+atn
              nout2=nout1+after
              do j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(s - r)*rt2i
                s2=(r + s)*rt2i
                zout(1,j,nout1)= r2 + r1
                zout(2,j,nout1)= s1 - s2
                zout(1,j,nout2)= r1 - r2
                zout(2,j,nout2)= s2 + s1
              end do
            end do
          endif
        else
          itrig=ias*before+1
          cr2=trig(1,itrig)
          ci2=trig(2,itrig)
          nin1=ia-after
          nout1=ia-atn
          do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nout1=nout1+atn
            nout2=nout1+after
            do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r=zin(1,j,nin2)
              s=zin(2,j,nin2)
              r2=r*cr2 - s*ci2
              s2=r*ci2 + s*cr2
              zout(1,j,nout1)= r2 + r1
              zout(2,j,nout1)= s2 + s1
              zout(1,j,nout2)= r1 - r2
              zout(2,j,nout2)= s1 - s2
            end do
          end do
        endif
      end do
    else if (now.eq.4) then
      if (isign.eq.1) then 
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib = 1, before
          nin1=nin1+after
          nin2=nin1+atb
          nin3=nin2+atb
          nin4=nin3+atb
          nout1=nout1+atn
          nout2=nout1+after
          nout3=nout2+after
          nout4=nout3+after
          do j = 1, nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            r3=zin(1,j,nin3)
            s3=zin(2,j,nin3)
            r4=zin(1,j,nin4)
            s4=zin(2,j,nin4)
            r=r1 + r3
            s=r2 + r4
            zout(1,j,nout1) = r + s
            zout(1,j,nout3) = r - s
            r=r1 - r3
            s=s2 - s4
            zout(1,j,nout2) = r - s 
            zout(1,j,nout4) = r + s
            r=s1 + s3
            s=s2 + s4
            zout(2,j,nout1) = r + s 
            zout(2,j,nout3) = r - s
            r=s1 - s3
            s=r2 - r4
            zout(2,j,nout2) = r + s 
            zout(2,j,nout4) = r - s
          end do
        end do
        do ia = 2, after
          ias=ia-1
          if (2*ias.eq.after) then
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r-s)*rt2i
                s2=(r+s)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=(r + s)*rt2i
                s4=(r - s)*rt2i
                r=r1 - r3
                s=r2 - r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 + r3
                s=s2 - s4
                zout(1,j,nout2) = r - s 
                zout(1,j,nout4) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s 
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 + r4
                zout(2,j,nout2) = r + s 
                zout(2,j,nout4) = r - s
              end do
            end do
          else
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            itrig=itrig+itt
            cr4=trig(1,itrig)
            ci4=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                r=zin(1,j,nin3)
                s=zin(2,j,nin3)
                r3=r*cr3 - s*ci3
                s3=r*ci3 + s*cr3
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=r*cr4 - s*ci4
                s4=r*ci4 + s*cr4
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r - s 
                zout(1,j,nout4) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s 
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r + s 
                zout(2,j,nout4) = r - s
              end do
            end do
          endif
        end do
      else
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib = 1, before
          nin1=nin1+after
          nin2=nin1+atb
          nin3=nin2+atb
          nin4=nin3+atb
          nout1=nout1+atn
          nout2=nout1+after
          nout3=nout2+after
          nout4=nout3+after
          do j = 1, nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            r3=zin(1,j,nin3)
            s3=zin(2,j,nin3)
            r4=zin(1,j,nin4)
            s4=zin(2,j,nin4)
            r=r1 + r3
            s=r2 + r4
            zout(1,j,nout1) = r + s
            zout(1,j,nout3) = r - s
            r=r1 - r3
            s=s2 - s4
            zout(1,j,nout2) = r + s
            zout(1,j,nout4) = r - s
            r=s1 + s3
            s=s2 + s4
            zout(2,j,nout1) = r + s
            zout(2,j,nout3) = r - s
            r=s1 - s3
            s=r2 - r4
            zout(2,j,nout2) = r - s
            zout(2,j,nout4) = r + s
          end do
        end do
        do ia = 2, after
          ias=ia-1
          if (2*ias.eq.after) then
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              do j = 1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(s - r)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=(s - r)*rt2i
                s4=(r + s)*rt2i
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 + s4
                zout(1,j,nout2) = r + s
                zout(1,j,nout4) = r - s
                r=s1 - s3
                s=s2 - s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 + s3
                s=r2 - r4
                zout(2,j,nout2) = r - s
                zout(2,j,nout4) = r + s
              end do
            end do
          else
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            itrig=itrig+itt
            cr4=trig(1,itrig)
            ci4=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                r=zin(1,j,nin3)
                s=zin(2,j,nin3)
                r3=r*cr3 - s*ci3
                s3=r*ci3 + s*cr3
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=r*cr4 - s*ci4
                s4=r*ci4 + s*cr4
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r + s
                zout(1,j,nout4) = r - s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r - s
                zout(2,j,nout4) = r + s
              end do
            end do
          endif
        end do
      endif
    else if (now.eq.8) then
      if (isign.eq.-1) then 
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib = 1, before
          nin1=nin1+after
          nin2=nin1+atb
          nin3=nin2+atb
          nin4=nin3+atb
          nin5=nin4+atb
          nin6=nin5+atb
          nin7=nin6+atb
          nin8=nin7+atb
          nout1=nout1+atn
          nout2=nout1+after
          nout3=nout2+after
          nout4=nout3+after
          nout5=nout4+after
          nout6=nout5+after
          nout7=nout6+after
          nout8=nout7+after
          do j = 1, nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            r3=zin(1,j,nin3)
            s3=zin(2,j,nin3)
            r4=zin(1,j,nin4)
            s4=zin(2,j,nin4)
            r5=zin(1,j,nin5)
            s5=zin(2,j,nin5)
            r6=zin(1,j,nin6)
            s6=zin(2,j,nin6)
            r7=zin(1,j,nin7)
            s7=zin(2,j,nin7)
            r8=zin(1,j,nin8)
            s8=zin(2,j,nin8)
            r=r1 + r5
            s=r3 + r7
            ap=r + s
            am=r - s
            r=r2 + r6
            s=r4 + r8
            bp=r + s
            bm=r - s
            r=s1 + s5
            s=s3 + s7
            cp=r + s
            cm=r - s
            r=s2 + s6
            s=s4 + s8
            dp=r + s
            dm=r - s
            zout(1,j,nout1) = ap + bp
            zout(2,j,nout1) = cp + dp
            zout(1,j,nout5) = ap - bp
            zout(2,j,nout5) = cp - dp
            zout(1,j,nout3) = am + dm
            zout(2,j,nout3) = cm - bm
            zout(1,j,nout7) = am - dm
            zout(2,j,nout7) = cm + bm
            r=r1 - r5
            s=s3 - s7
            ap=r + s
            am=r - s
            r=s1 - s5
            s=r3 - r7
            bp=r + s
            bm=r - s
            r=s4 - s8
            s=r2 - r6
            cp=r + s
            cm=r - s
            r=s2 - s6
            s=r4 - r8
            dp=r + s
            dm=r - s
            r = ( cp + dm)*rt2i
            s = ( dm - cp)*rt2i
            cp= ( cm + dp)*rt2i
            dp = ( cm - dp)*rt2i
            zout(1,j,nout2) = ap + r
            zout(2,j,nout2) = bm + s
            zout(1,j,nout6) = ap - r
            zout(2,j,nout6) = bm - s
            zout(1,j,nout4) = am + cp
            zout(2,j,nout4) = bp + dp
            zout(1,j,nout8) = am - cp
            zout(2,j,nout8) = bp - dp
          end do
        end do
        do ia = 2, after
          ias=ia-1
          itt=ias*before
          itrig=itt+1
          cr2=trig(1,itrig)
          ci2=trig(2,itrig)
          itrig=itrig+itt
          cr3=trig(1,itrig)
          ci3=trig(2,itrig)
          itrig=itrig+itt
          cr4=trig(1,itrig)
          ci4=trig(2,itrig)
          itrig=itrig+itt
          cr5=trig(1,itrig)
          ci5=trig(2,itrig)
          itrig=itrig+itt
          cr6=trig(1,itrig)
          ci6=trig(2,itrig)
          itrig=itrig+itt
          cr7=trig(1,itrig)
          ci7=trig(2,itrig)
          itrig=itrig+itt
          cr8=trig(1,itrig)
          ci8=trig(2,itrig)
          nin1=ia-after
          nout1=ia-atn
          do ib = 1, before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nin5=nin4+atb
            nin6=nin5+atb
            nin7=nin6+atb
            nin8=nin7+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            nout5=nout4+after
            nout6=nout5+after
            nout7=nout6+after
            nout8=nout7+after
            do j = 1, nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r=zin(1,j,nin2)
              s=zin(2,j,nin2)
              r2=r*cr2 - s*ci2
              s2=r*ci2 + s*cr2
              r=zin(1,j,nin3)
              s=zin(2,j,nin3)
              r3=r*cr3 - s*ci3
              s3=r*ci3 + s*cr3
              r=zin(1,j,nin4)
              s=zin(2,j,nin4)
              r4=r*cr4 - s*ci4
              s4=r*ci4 + s*cr4
              r=zin(1,j,nin5)
              s=zin(2,j,nin5)
              r5=r*cr5 - s*ci5
              s5=r*ci5 + s*cr5
              r=zin(1,j,nin6)
              s=zin(2,j,nin6)
              r6=r*cr6 - s*ci6
              s6=r*ci6 + s*cr6
              r=zin(1,j,nin7)
              s=zin(2,j,nin7)
              r7=r*cr7 - s*ci7
              s7=r*ci7 + s*cr7
              r=zin(1,j,nin8)
              s=zin(2,j,nin8)
              r8=r*cr8 - s*ci8
              s8=r*ci8 + s*cr8
              r=r1 + r5
              s=r3 + r7
              ap=r + s
              am=r - s
              r=r2 + r6
              s=r4 + r8
              bp=r + s
              bm=r - s
              r=s1 + s5
              s=s3 + s7
              cp=r + s
              cm=r - s
              r=s2 + s6
              s=s4 + s8
              dp=r + s
              dm=r - s
              zout(1,j,nout1) = ap + bp
              zout(2,j,nout1) = cp + dp
              zout(1,j,nout5) = ap - bp
              zout(2,j,nout5) = cp - dp
              zout(1,j,nout3) = am + dm
              zout(2,j,nout3) = cm - bm
              zout(1,j,nout7) = am - dm
              zout(2,j,nout7) = cm + bm
              r=r1 - r5
              s=s3 - s7
              ap=r + s
              am=r - s
              r=s1 - s5
              s=r3 - r7
              bp=r + s
              bm=r - s
              r=s4 - s8
              s=r2 - r6
              cp=r + s
              cm=r - s
              r=s2 - s6
              s=r4 - r8
              dp=r + s
              dm=r - s
              r = ( cp + dm)*rt2i
              s = ( dm - cp)*rt2i
              cp= ( cm + dp)*rt2i
              dp = ( cm - dp)*rt2i
              zout(1,j,nout2) = ap + r
              zout(2,j,nout2) = bm + s
              zout(1,j,nout6) = ap - r
              zout(2,j,nout6) = bm - s
              zout(1,j,nout4) = am + cp
              zout(2,j,nout4) = bp + dp
              zout(1,j,nout8) = am - cp
              zout(2,j,nout8) = bp - dp
            end do
          end do
        end do
      else
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib = 1, before
          nin1=nin1+after
          nin2=nin1+atb
          nin3=nin2+atb
          nin4=nin3+atb
          nin5=nin4+atb
          nin6=nin5+atb
          nin7=nin6+atb
          nin8=nin7+atb
          nout1=nout1+atn
          nout2=nout1+after
          nout3=nout2+after
          nout4=nout3+after
          nout5=nout4+after
          nout6=nout5+after
          nout7=nout6+after
          nout8=nout7+after
          do j = 1, nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            r3=zin(1,j,nin3)
            s3=zin(2,j,nin3)
            r4=zin(1,j,nin4)
            s4=zin(2,j,nin4)
            r5=zin(1,j,nin5)
            s5=zin(2,j,nin5)
            r6=zin(1,j,nin6)
            s6=zin(2,j,nin6)
            r7=zin(1,j,nin7)
            s7=zin(2,j,nin7)
            r8=zin(1,j,nin8)
            s8=zin(2,j,nin8)
            r=r1 + r5
            s=r3 + r7
            ap=r + s
            am=r - s
            r=r2 + r6
            s=r4 + r8
            bp=r + s
            bm=r - s
            r=s1 + s5
            s=s3 + s7
            cp=r + s
            cm=r - s
            r=s2 + s6
            s=s4 + s8
            dp=r + s
            dm=r - s
            zout(1,j,nout1) = ap + bp
            zout(2,j,nout1) = cp + dp
            zout(1,j,nout5) = ap - bp
            zout(2,j,nout5) = cp - dp
            zout(1,j,nout3) = am - dm
            zout(2,j,nout3) = cm + bm
            zout(1,j,nout7) = am + dm
            zout(2,j,nout7) = cm - bm
            r= r1 - r5
            s=-s3 + s7
            ap=r + s
            am=r - s
            r=s1 - s5
            s=r7 - r3
            bp=r + s
            bm=r - s
            r=-s4 + s8
            s= r2 - r6
            cp=r + s
            cm=r - s
            r=-s2 + s6
            s= r4 - r8
            dp=r + s
            dm=r - s
            r = ( cp + dm)*rt2i
            s = ( cp - dm)*rt2i
            cp= ( cm + dp)*rt2i
            dp= ( dp - cm)*rt2i
            zout(1,j,nout2) = ap + r
            zout(2,j,nout2) = bm + s
            zout(1,j,nout6) = ap - r
            zout(2,j,nout6) = bm - s
            zout(1,j,nout4) = am + cp
            zout(2,j,nout4) = bp + dp
            zout(1,j,nout8) = am - cp
            zout(2,j,nout8) = bp - dp
          end do
        end do

        do ia = 2, after
          ias=ia-1
          itt=ias*before
          itrig=itt+1
          cr2=trig(1,itrig)
          ci2=trig(2,itrig)
          itrig=itrig+itt
          cr3=trig(1,itrig)
          ci3=trig(2,itrig)
          itrig=itrig+itt
          cr4=trig(1,itrig)
          ci4=trig(2,itrig)
          itrig=itrig+itt
          cr5=trig(1,itrig)
          ci5=trig(2,itrig)
          itrig=itrig+itt
          cr6=trig(1,itrig)
          ci6=trig(2,itrig)
          itrig=itrig+itt
          cr7=trig(1,itrig)
          ci7=trig(2,itrig)
          itrig=itrig+itt
          cr8=trig(1,itrig)
          ci8=trig(2,itrig)
          nin1=ia-after
          nout1=ia-atn
          do ib = 1, before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nin5=nin4+atb
            nin6=nin5+atb
            nin7=nin6+atb
            nin8=nin7+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            nout5=nout4+after
            nout6=nout5+after
            nout7=nout6+after
            nout8=nout7+after
            do j = 1, nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r=zin(1,j,nin2)
              s=zin(2,j,nin2)
              r2=r*cr2 - s*ci2
              s2=r*ci2 + s*cr2
              r=zin(1,j,nin3)
              s=zin(2,j,nin3)
              r3=r*cr3 - s*ci3
              s3=r*ci3 + s*cr3
              r=zin(1,j,nin4)
              s=zin(2,j,nin4)
              r4=r*cr4 - s*ci4
              s4=r*ci4 + s*cr4
              r=zin(1,j,nin5)
              s=zin(2,j,nin5)
              r5=r*cr5 - s*ci5
              s5=r*ci5 + s*cr5
              r=zin(1,j,nin6)
              s=zin(2,j,nin6)
              r6=r*cr6 - s*ci6
              s6=r*ci6 + s*cr6
              r=zin(1,j,nin7)
              s=zin(2,j,nin7)
              r7=r*cr7 - s*ci7
              s7=r*ci7 + s*cr7
              r=zin(1,j,nin8)
              s=zin(2,j,nin8)
              r8=r*cr8 - s*ci8
              s8=r*ci8 + s*cr8
              r=r1 + r5
              s=r3 + r7
              ap=r + s
              am=r - s
              r=r2 + r6
              s=r4 + r8
              bp=r + s
              bm=r - s
              r=s1 + s5
              s=s3 + s7
              cp=r + s
              cm=r - s
              r=s2 + s6
              s=s4 + s8
              dp=r + s
              dm=r - s
              zout(1,j,nout1) = ap + bp
              zout(2,j,nout1) = cp + dp
              zout(1,j,nout5) = ap - bp
              zout(2,j,nout5) = cp - dp
              zout(1,j,nout3) = am - dm
              zout(2,j,nout3) = cm + bm
              zout(1,j,nout7) = am + dm
              zout(2,j,nout7) = cm - bm
              r= r1 - r5
              s=-s3 + s7
              ap=r + s
              am=r - s
              r=s1 - s5
              s=r7 - r3
              bp=r + s
              bm=r - s
              r=-s4 + s8
              s= r2 - r6
              cp=r + s
              cm=r - s
              r=-s2 + s6
              s= r4 - r8
              dp=r + s
              dm=r - s
              r = ( cp + dm)*rt2i
              s = ( cp - dm)*rt2i
              cp= ( cm + dp)*rt2i
              dp= ( dp - cm)*rt2i
              zout(1,j,nout2) = ap + r
              zout(2,j,nout2) = bm + s
              zout(1,j,nout6) = ap - r
              zout(2,j,nout6) = bm - s
              zout(1,j,nout4) = am + cp
              zout(2,j,nout4) = bp + dp
              zout(1,j,nout8) = am - cp
              zout(2,j,nout8) = bp - dp
            end do
          end do
        end do

      end if
    else if (now.eq.3) then 
      !         .5d0*sqrt(3.d0)
      bb=isign*0.8660254037844387d0
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib = 1, before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do j = 1, nfft
          r1=zin(1,j,nin1)
          s1=zin(2,j,nin1)
          r2=zin(1,j,nin2)
          s2=zin(2,j,nin2)
          r3=zin(1,j,nin3)
          s3=zin(2,j,nin3)
          r=r2 + r3
          s=s2 + s3
          zout(1,j,nout1) = r + r1
          zout(2,j,nout1) = s + s1
          r1=r1 - .5d0*r
          s1=s1 - .5d0*s
          r2=bb*(r2-r3)
          s2=bb*(s2-s3)
          zout(1,j,nout2) = r1 - s2 
          zout(2,j,nout2) = s1 + r2
          zout(1,j,nout3) = r1 + s2 
          zout(2,j,nout3) = s1 - r2
        end do
      end do
      do ia = 2, after
        ias=ia-1
        if (4*ias.eq.3*after) then
          if (isign.eq.1) then
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r3 + r2
                s=s2 - s3
                zout(1,j,nout1) = r1 - r
                zout(2,j,nout1) = s + s1
                r1=r1 + .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2-r3)        
                s2=bb*(s2+s3)
                zout(1,j,nout2) = r1 - s2 
                zout(2,j,nout2) = s1 - r2
                zout(1,j,nout3) = r1 + s2 
                zout(2,j,nout3) = s1 + r2
              end do
            end do
          else
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r2 - r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s1 - s
                r1=r1 - .5d0*r
                s1=s1 + .5d0*s        
                r2=bb*(r2+r3)        
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 + s2 
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 - s2 
                zout(2,j,nout3) = s1 - r2
              end do
            end do
          endif
        else if (8*ias.eq.3*after) then
          if (isign.eq.1) then
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3) 
                r=r2 - r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2+r3)        
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2 
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2 
                zout(2,j,nout3) = s1 - r2
              end do
            end do
          else
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(s - r)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=r2 + r3
                s=s2 - s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2-r3)        
                s2=bb*(s2+s3)
                zout(1,j,nout2) = r1 - s2 
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2 
                zout(2,j,nout3) = s1 - r2
              end do
            end do
          endif
        else
          itt=ias*before
          itrig=itt+1
          cr2=trig(1,itrig)
          ci2=trig(2,itrig)
          itrig=itrig+itt
          cr3=trig(1,itrig)
          ci3=trig(2,itrig)
          nin1=ia-after
          nout1=ia-atn
          do ib = 1, before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            do j = 1, nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r=zin(1,j,nin2)
              s=zin(2,j,nin2)
              r2=r*cr2 - s*ci2
              s2=r*ci2 + s*cr2
              r=zin(1,j,nin3)
              s=zin(2,j,nin3)
              r3=r*cr3 - s*ci3
              s3=r*ci3 + s*cr3
              r=r2 + r3
              s=s2 + s3
              zout(1,j,nout1) = r + r1
              zout(2,j,nout1) = s + s1
              r1=r1 - .5d0*r
              s1=s1 - .5d0*s
              r2=bb*(r2-r3)
              s2=bb*(s2-s3)
              zout(1,j,nout2) = r1 - s2 
              zout(2,j,nout2) = s1 + r2
              zout(1,j,nout3) = r1 + s2 
              zout(2,j,nout3) = s1 - r2
            end do
          end do
        end if
      end do
    else if (now.eq.5) then
      !         cos(2.d0*pi/5.d0)
      cos2=0.3090169943749474d0
      !         cos(4.d0*pi/5.d0)
      cos4=-0.8090169943749474d0
      !        sin(2.d0*pi/5.d0)
      sin2=isign*0.9510565162951536d0
      !         sin(4.d0*pi/5.d0)
      sin4=isign*0.5877852522924731d0
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib = 1, before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        do j = 1, nfft
          r1=zin(1,j,nin1)
          s1=zin(2,j,nin1)
          r2=zin(1,j,nin2)
          s2=zin(2,j,nin2)
          r3=zin(1,j,nin3)
          s3=zin(2,j,nin3)
          r4=zin(1,j,nin4)
          s4=zin(2,j,nin4)
          r5=zin(1,j,nin5)
          s5=zin(2,j,nin5)
          r25 = r2 + r5
          r34 = r3 + r4
          s25 = s2 - s5
          s34 = s3 - s4
          zout(1,j,nout1) = r1 + r25 + r34
          r = r1 + cos2*r25 + cos4*r34
          s = sin2*s25 + sin4*s34
          zout(1,j,nout2) = r - s
          zout(1,j,nout5) = r + s
          r = r1 + cos4*r25 + cos2*r34
          s = sin4*s25 - sin2*s34
          zout(1,j,nout3) = r - s
          zout(1,j,nout4) = r + s
          r25 = r2 - r5
          r34 = r3 - r4
          s25 = s2 + s5
          s34 = s3 + s4
          zout(2,j,nout1) = s1 + s25 + s34
          r = s1 + cos2*s25 + cos4*s34
          s = sin2*r25 + sin4*r34
          zout(2,j,nout2) = r + s
          zout(2,j,nout5) = r - s
          r = s1 + cos4*s25 + cos2*s34
          s = sin4*r25 - sin2*r34
          zout(2,j,nout3) = r + s
          zout(2,j,nout4) = r - s
        end do
      end do
      do ia = 2, after
        ias=ia-1
        if (8*ias.eq.5*after) then
          if (isign.eq.1) then
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nin5=nin4+atb        
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              nout5=nout4+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3) 
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=(r + s)*rt2i
                s4=(r - s)*rt2i
                r5=zin(1,j,nin5)
                s5=zin(2,j,nin5)
                r25 = r2 - r5
                r34 = r3 + r4
                s25 = s2 + s5
                s34 = s3 - s4
                zout(1,j,nout1) = r1 + r25 - r34
                r = r1 + cos2*r25 - cos4*r34 
                s = sin2*s25 + sin4*s34
                zout(1,j,nout2) = r - s
                zout(1,j,nout5) = r + s
                r = r1 + cos4*r25 - cos2*r34 
                s = sin4*s25 - sin2*s34
                zout(1,j,nout3) = r - s
                zout(1,j,nout4) = r + s
                r25 = r2 + r5
                r34 = r4 - r3
                s25 = s2 - s5
                s34 = s3 + s4
                zout(2,j,nout1) = s1 + s25 + s34
                r = s1 + cos2*s25 + cos4*s34
                s = sin2*r25 + sin4*r34
                zout(2,j,nout2) = r + s
                zout(2,j,nout5) = r - s
                r = s1 + cos4*s25 + cos2*s34
                s = sin4*r25 - sin2*r34
                zout(2,j,nout3) = r + s
                zout(2,j,nout4) = r - s
              end do
            end do
          else
            nin1=ia-after
            nout1=ia-atn
            do ib = 1, before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nin5=nin4+atb        
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              nout5=nout4+after
              do j = 1, nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(s - r)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=(s - r)*rt2i
                s4=(r + s)*rt2i
                r5=zin(1,j,nin5)
                s5=zin(2,j,nin5)
                r25 = r2 - r5
                r34 = r3 + r4
                s25 = s2 + s5
                s34 = s4 - s3
                zout(1,j,nout1) = r1 + r25 + r34
                r = r1 + cos2*r25 + cos4*r34
                s = sin2*s25 + sin4*s34
                zout(1,j,nout2) = r - s
                zout(1,j,nout5) = r + s
                r = r1 + cos4*r25 + cos2*r34
                s = sin4*s25 - sin2*s34
                zout(1,j,nout3) = r - s
                zout(1,j,nout4) = r + s
                r25 = r2 + r5
                r34 = r3 - r4
                s25 = s2 - s5
                s34 = s3 + s4
                zout(2,j,nout1) = s1 + s25 - s34
                r = s1 + cos2*s25 - cos4*s34
                s = sin2*r25 + sin4*r34
                zout(2,j,nout2) = r + s
                zout(2,j,nout5) = r - s
                r = s1 + cos4*s25 - cos2*s34
                s = sin4*r25 - sin2*r34
                zout(2,j,nout3) = r + s
                zout(2,j,nout4) = r - s
              end do
            end do
          endif
        else
          ias=ia-1
          itt=ias*before
          itrig=itt+1
          cr2=trig(1,itrig)
          ci2=trig(2,itrig)
          itrig=itrig+itt
          cr3=trig(1,itrig)
          ci3=trig(2,itrig)
          itrig=itrig+itt
          cr4=trig(1,itrig)
          ci4=trig(2,itrig)
          itrig=itrig+itt
          cr5=trig(1,itrig)
          ci5=trig(2,itrig)
          nin1=ia-after
          nout1=ia-atn
          do ib = 1, before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nin5=nin4+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            nout5=nout4+after
            do j = 1, nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r=zin(1,j,nin2)
              s=zin(2,j,nin2)
              r2=r*cr2 - s*ci2
              s2=r*ci2 + s*cr2
              r=zin(1,j,nin3)
              s=zin(2,j,nin3)
              r3=r*cr3 - s*ci3
              s3=r*ci3 + s*cr3
              r=zin(1,j,nin4)
              s=zin(2,j,nin4)
              r4=r*cr4 - s*ci4
              s4=r*ci4 + s*cr4
              r=zin(1,j,nin5)
              s=zin(2,j,nin5)
              r5=r*cr5 - s*ci5
              s5=r*ci5 + s*cr5
              r25 = r2 + r5
              r34 = r3 + r4
              s25 = s2 - s5
              s34 = s3 - s4
              zout(1,j,nout1) = r1 + r25 + r34
              r = r1 + cos2*r25 + cos4*r34
              s = sin2*s25 + sin4*s34
              zout(1,j,nout2) = r - s
              zout(1,j,nout5) = r + s
              r = r1 + cos4*r25 + cos2*r34
              s = sin4*s25 - sin2*s34
              zout(1,j,nout3) = r - s
              zout(1,j,nout4) = r + s
              r25 = r2 - r5
              r34 = r3 - r4
              s25 = s2 + s5
              s34 = s3 + s4
              zout(2,j,nout1) = s1 + s25 + s34
              r = s1 + cos2*s25 + cos4*s34
              s = sin2*r25 + sin4*r34
              zout(2,j,nout2) = r + s
              zout(2,j,nout5) = r - s
              r = s1 + cos4*s25 + cos2*s34
              s = sin4*r25 - sin2*r34
              zout(2,j,nout3) = r + s
              zout(2,j,nout4) = r - s
            end do
          end do
        endif
      end do
    else if (now.eq.6) then
      !         .5d0*sqrt(3.d0)
      bb=isign*0.8660254037844387d0

      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib = 1, before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        do j = 1, nfft
          r2=zin(1,j,nin3)
          s2=zin(2,j,nin3)
          r3=zin(1,j,nin5)
          s3=zin(2,j,nin5)
          r=r2 + r3
          s=s2 + s3
          r1=zin(1,j,nin1)
          s1=zin(2,j,nin1)
          ur1 = r + r1
          ui1 = s + s1
          r1=r1 - .5d0*r
          s1=s1 - .5d0*s
          r=r2-r3
          s=s2-s3
          ur2 = r1 - s*bb
          ui2 = s1 + r*bb
          ur3 = r1 + s*bb
          ui3 = s1 - r*bb

          r2=zin(1,j,nin6)
          s2=zin(2,j,nin6)
          r3=zin(1,j,nin2)
          s3=zin(2,j,nin2)
          r=r2 + r3
          s=s2 + s3
          r1=zin(1,j,nin4)
          s1=zin(2,j,nin4)
          vr1 = r + r1
          vi1 = s + s1
          r1=r1 - .5d0*r
          s1=s1 - .5d0*s
          r=r2-r3
          s=s2-s3
          vr2 = r1 - s*bb
          vi2 = s1 + r*bb
          vr3 = r1 + s*bb
          vi3 = s1 - r*bb

          zout(1,j,nout1)=ur1+vr1
          zout(2,j,nout1)=ui1+vi1
          zout(1,j,nout5)=ur2+vr2
          zout(2,j,nout5)=ui2+vi2
          zout(1,j,nout3)=ur3+vr3
          zout(2,j,nout3)=ui3+vi3
          zout(1,j,nout4)=ur1-vr1
          zout(2,j,nout4)=ui1-vi1
          zout(1,j,nout2)=ur2-vr2
          zout(2,j,nout2)=ui2-vi2
          zout(1,j,nout6)=ur3-vr3
          zout(2,j,nout6)=ui3-vi3
        end do
      end do
    else 
      stop 'error fftstp'
    endif

  end subroutine fftstp



  !  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
  !  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1995, 1999
  !  This file is distributed under the terms of the
  !  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
  subroutine fftrot(mm,nfft,m,nn,n,zin,zout,trig,after,now,before,isign)
    integer :: mm, nfft, m, nn, n, isign
    integer :: after, before, atn, atb, now
    real(8) :: trig(2,8192), zin(2, mm, m), zout(2, n, nn)

    real(8) :: rt2i, dp, cp, cm, ci5, cr5, ci6, cr6, am, ap, ci8, cr8
    real(8) :: r, r1, r2, r3, r4, r5, r6, r7, r8, r25, r34
    real(8) :: s, s1, s2, s3, s4, s5, s6, s7, s8, s25, s34
    real(8) :: bb, bm, dm, bp
    real(8) :: cr2, ci2, cr3, ci3, cr4, ci4, cr7, ci7
    real(8) :: cos2, sin2, cos4, sin4
    real(8) :: ur1, ur2, ur3, ui1, ui2, ui3
    real(8) :: vr1, vr2, vr3, vi1, vi2, vi3
    integer :: ia, ib, j, ias, itrig, itt
    integer :: nin1, nin2, nin3, nin4, nin5, nin6, nin7, nin8
    integer :: nout1, nout2, nout3, nout4, nout5, nout6, nout7, nout8
    
    atn=after*now
    atb=after*before
    
    !         sqrt(.5d0)
    rt2i=0.7071067811865475d0
    if (now.eq.2) then
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do 2001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nout1=nout1+atn
        nout2=nout1+after
        do 2001,j=1,nfft
          r1=zin(1,j,nin1)
          s1=zin(2,j,nin1)
          r2=zin(1,j,nin2)
          s2=zin(2,j,nin2)
          zout(1,nout1,j)= r2 + r1
          zout(2,nout1,j)= s2 + s1
          zout(1,nout2,j)= r1 - r2
          zout(2,nout2,j)= s1 - s2
2001      continue
        do 2000,ia=2,after
        ias=ia-1
        if (2*ias.eq.after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(2,j,nin2)
                        s2=zin(1,j,nin2)
                        zout(1,nout1,j)= r1 - r2
                        zout(2,nout1,j)= s2 + s1
                        zout(1,nout2,j)= r2 + r1
                        zout(2,nout2,j)= s1 - s2
2010                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(2,j,nin2)
                        s2=zin(1,j,nin2)
                        zout(1,nout1,j)= r2 + r1
                        zout(2,nout1,j)= s1 - s2
                        zout(1,nout2,j)= r1 - r2
                        zout(2,nout2,j)= s2 + s1
2020                        continue
                endif
        else if (4*ias.eq.after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2030,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2030,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        zout(1,nout1,j)= r2 + r1
                        zout(2,nout1,j)= s2 + s1
                        zout(1,nout2,j)= r1 - r2
                        zout(2,nout2,j)= s1 - s2
2030                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2040,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2040,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        zout(1,nout1,j)= r2 + r1
                        zout(2,nout1,j)= s2 + s1
                        zout(1,nout2,j)= r1 - r2
                        zout(2,nout2,j)= s1 - s2
2040                        continue
                endif
        else if (4*ias.eq.3*after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2050,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2050,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(r - s)*rt2i
                        zout(1,nout1,j)= r1 - r2
                        zout(2,nout1,j)= s2 + s1
                        zout(1,nout2,j)= r2 + r1
                        zout(2,nout2,j)= s1 - s2
2050                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2060,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2060,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(s - r)*rt2i
                        s2=(r + s)*rt2i
                        zout(1,nout1,j)= r2 + r1
                        zout(2,nout1,j)= s1 - s2
                        zout(1,nout2,j)= r1 - r2
                        zout(2,nout2,j)= s2 + s1
2060                        continue
                endif
        else
                itrig=ias*before+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
                do 2090,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nout1=nout1+atn
                nout2=nout1+after
                do 2090,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                zout(1,nout1,j)= r2 + r1
                zout(2,nout1,j)= s2 + s1
                zout(1,nout2,j)= r1 - r2
                zout(2,nout2,j)= s1 - s2
2090                continue
        endif
2000        continue
        else if (now.eq.4) then
        if (isign.eq.1) then 
                ia=1
                nin1=ia-after
                nout1=ia-atn
                do 4001,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                do 4001,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,nout1,j) = r + s
                zout(1,nout3,j) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,nout2,j) = r - s 
                zout(1,nout4,j) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,nout1,j) = r + s 
                zout(2,nout3,j) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,nout2,j) = r + s 
                zout(2,nout4,j) = r - s
4001                continue
                do 4000,ia=2,after
                ias=ia-1
                if (2*ias.eq.after) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 4010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r-s)*rt2i
                        s2=(r+s)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r=r1 - r3
                        s=r2 - r4
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 + r3
                        s=s2 - s4
                        zout(1,nout2,j) = r - s 
                        zout(1,nout4,j) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s 
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 + r4
                        zout(2,nout2,j) = r + s 
                        zout(2,nout4,j) = r - s
4010                        continue
                else
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        do 4020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,nout2,j) = r - s 
                        zout(1,nout4,j) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s 
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,nout2,j) = r + s 
                        zout(2,nout4,j) = r - s
4020                        continue
                endif
4000                continue
        else
                ia=1
                nin1=ia-after
                nout1=ia-atn
                do 4101,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                do 4101,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,nout1,j) = r + s
                zout(1,nout3,j) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,nout2,j) = r + s
                zout(1,nout4,j) = r - s
                r=s1 + s3
                s=s2 + s4
                zout(2,nout1,j) = r + s
                zout(2,nout3,j) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,nout2,j) = r - s
                zout(2,nout4,j) = r + s
4101                continue
                do 4100,ia=2,after
                ias=ia-1
                if (2*ias.eq.after) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 4110,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4110,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=(r + s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 + s4
                        zout(1,nout2,j) = r + s
                        zout(1,nout4,j) = r - s
                        r=s1 - s3
                        s=s2 - s4
                        zout(2,nout1,j) = r + s
                        zout(2,nout3,j) = r - s
                        r=s1 + s3
                        s=r2 - r4
                        zout(2,nout2,j) = r - s
                        zout(2,nout4,j) = r + s
4110                        continue
                else
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        do 4120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4120,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,nout2,j) = r + s
                        zout(1,nout4,j) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,nout2,j) = r - s
                        zout(2,nout4,j) = r + s
4120                        continue
                endif
4100                continue
        endif
        else if (now.eq.8) then
        if (isign.eq.-1) then 
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        do 8120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        do 8120,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp
                        zout(1,nout3,j) = am + dm
                        zout(2,nout3,j) = cm - bm
                        zout(1,nout7,j) = am - dm
                        zout(2,nout7,j) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( dm - cp)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp= ( cm - dp)*rt2i
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp
8120                        continue
                do 8000,ia=2,after
                ias=ia-1
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        itrig=itrig+itt
                        cr5=trig(1,itrig)
                        ci5=trig(2,itrig)
                        itrig=itrig+itt
                        cr6=trig(1,itrig)
                        ci6=trig(2,itrig)
                        itrig=itrig+itt
                        cr7=trig(1,itrig)
                        ci7=trig(2,itrig)
                        itrig=itrig+itt
                        cr8=trig(1,itrig)
                        ci8=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        do 8020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        do 8020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=zin(1,j,nin5)
                        s=zin(2,j,nin5)
                        r5=r*cr5 - s*ci5
                        s5=r*ci5 + s*cr5
                        r=zin(1,j,nin6)
                        s=zin(2,j,nin6)
                        r6=r*cr6 - s*ci6
                        s6=r*ci6 + s*cr6
                        r=zin(1,j,nin7)
                        s=zin(2,j,nin7)
                        r7=r*cr7 - s*ci7
                        s7=r*ci7 + s*cr7
                        r=zin(1,j,nin8)
                        s=zin(2,j,nin8)
                        r8=r*cr8 - s*ci8
                        s8=r*ci8 + s*cr8
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp
                        zout(1,nout3,j) = am + dm
                        zout(2,nout3,j) = cm - bm
                        zout(1,nout7,j) = am - dm
                        zout(2,nout7,j) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( dm - cp)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp= ( cm - dp)*rt2i
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp

8020                        continue
8000                continue

        else
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        do 8121,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        do 8121,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp
                        zout(1,nout3,j) = am - dm
                        zout(2,nout3,j) = cm + bm
                        zout(1,nout7,j) = am + dm
                        zout(2,nout7,j) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp= ( dp - cm)*rt2i
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp
8121                        continue

                do 8001,ia=2,after
                ias=ia-1
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        itrig=itrig+itt
                        cr5=trig(1,itrig)
                        ci5=trig(2,itrig)
                        itrig=itrig+itt
                        cr6=trig(1,itrig)
                        ci6=trig(2,itrig)
                        itrig=itrig+itt
                        cr7=trig(1,itrig)
                        ci7=trig(2,itrig)
                        itrig=itrig+itt
                        cr8=trig(1,itrig)
                        ci8=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        do 8021,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        do 8021,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=zin(1,j,nin5)
                        s=zin(2,j,nin5)
                        r5=r*cr5 - s*ci5
                        s5=r*ci5 + s*cr5
                        r=zin(1,j,nin6)
                        s=zin(2,j,nin6)
                        r6=r*cr6 - s*ci6
                        s6=r*ci6 + s*cr6
                        r=zin(1,j,nin7)
                        s=zin(2,j,nin7)
                        r7=r*cr7 - s*ci7
                        s7=r*ci7 + s*cr7
                        r=zin(1,j,nin8)
                        s=zin(2,j,nin8)
                        r8=r*cr8 - s*ci8
                        s8=r*ci8 + s*cr8
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dp=r + s
                        dm=r - s
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp
                        zout(1,nout3,j) = am - dm
                        zout(2,nout3,j) = cm + bm
                        zout(1,nout7,j) = am + dm
                        zout(2,nout7,j) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dp)*rt2i
                        dp= ( dp - cm)*rt2i
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp
8021                        continue
8001                continue

        endif
        else if (now.eq.3) then 
!         .5d0*sqrt(3.d0)
        bb=isign*0.8660254037844387d0
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 3001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do 3001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r=r2 + r3
        s=s2 + s3
        zout(1,nout1,j) = r + r1
        zout(2,nout1,j) = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,nout2,j) = r1 - s2 
        zout(2,nout2,j) = s1 + r2
        zout(1,nout3,j) = r1 + s2 
        zout(2,nout3,j) = s1 - r2
3001        continue
        do 3000,ia=2,after
        ias=ia-1
        if (4*ias.eq.3*after) then
        if (isign.eq.1) then
                nin1=ia-after
                nout1=ia-atn
                do 3010,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3010,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r2 + r3
                s=s2 - s3
                zout(1,nout1,j) = r1 - r
                zout(2,nout1,j) = s + s1
                r1=r1 + .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2-r3)        
                s2=bb*(s2+s3)
                zout(1,nout2,j) = r1 - s2 
                zout(2,nout2,j) = s1 - r2
                zout(1,nout3,j) = r1 + s2 
                zout(2,nout3,j) = s1 + r2
3010                continue
        else
                nin1=ia-after
                nout1=ia-atn
                do 3020,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3020,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r2 - r3
                s=s2 + s3
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s1 - s
                r1=r1 - .5d0*r
                s1=s1 + .5d0*s        
                r2=bb*(r2+r3)        
                s2=bb*(s2-s3)
                zout(1,nout2,j) = r1 + s2 
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 - s2 
                zout(2,nout3,j) = s1 - r2
3020                continue
        endif
        else if (8*ias.eq.3*after) then
        if (isign.eq.1) then
                nin1=ia-after
                nout1=ia-atn
                do 3030,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3030,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3) 
                r=r2 - r3
                s=s2 + s3
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s + s1
                r1=r1 - .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2+r3)        
                s2=bb*(s2-s3)
                zout(1,nout2,j) = r1 - s2 
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 + s2 
                zout(2,nout3,j) = s1 - r2
3030                continue
        else
                nin1=ia-after
                nout1=ia-atn
                do 3040,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3040,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(s - r)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=r2 + r3
                s=s2 - s3
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s + s1
                r1=r1 - .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2-r3)        
                s2=bb*(s2+s3)
                zout(1,nout2,j) = r1 - s2 
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 + s2 
                zout(2,nout3,j) = s1 - r2
3040                continue
        endif
        else
        itt=ias*before
        itrig=itt+1
        cr2=trig(1,itrig)
        ci2=trig(2,itrig)
        itrig=itrig+itt
        cr3=trig(1,itrig)
        ci3=trig(2,itrig)
        nin1=ia-after
        nout1=ia-atn
        do 3090,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do 3090,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r=zin(1,j,nin2)
        s=zin(2,j,nin2)
        r2=r*cr2 - s*ci2
        s2=r*ci2 + s*cr2
        r=zin(1,j,nin3)
        s=zin(2,j,nin3)
        r3=r*cr3 - s*ci3
        s3=r*ci3 + s*cr3
        r=r2 + r3
        s=s2 + s3
        zout(1,nout1,j) = r + r1
        zout(2,nout1,j) = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,nout2,j) = r1 - s2 
        zout(2,nout2,j) = s1 + r2
        zout(1,nout3,j) = r1 + s2 
        zout(2,nout3,j) = s1 - r2
3090        continue
        endif
3000        continue
        else if (now.eq.5) then
!         cos(2.d0*pi/5.d0)
        cos2=0.3090169943749474d0
!         cos(4.d0*pi/5.d0)
        cos4=-0.8090169943749474d0
!        sin(2.d0*pi/5.d0)
        sin2=isign*0.9510565162951536d0
!         sin(4.d0*pi/5.d0)
        sin4=isign*0.5877852522924731d0
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 5001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        do 5001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r4=zin(1,j,nin4)
        s4=zin(2,j,nin4)
        r5=zin(1,j,nin5)
        s5=zin(2,j,nin5)
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout(1,nout1,j) = r1 + r25 + r34
        r = r1 + cos2*r25 + cos4*r34
        s = sin2*s25 + sin4*s34
        zout(1,nout2,j) = r - s
        zout(1,nout5,j) = r + s
        r = r1 + cos4*r25 + cos2*r34
        s = sin4*s25 - sin2*s34
        zout(1,nout3,j) = r - s
        zout(1,nout4,j) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout(2,nout1,j) = s1 + s25 + s34
        r = s1 + cos2*s25 + cos4*s34
        s = sin2*r25 + sin4*r34
        zout(2,nout2,j) = r + s
        zout(2,nout5,j) = r - s
        r = s1 + cos4*s25 + cos2*s34
        s = sin4*r25 - sin2*r34
        zout(2,nout3,j) = r + s
        zout(2,nout4,j) = r - s
5001        continue
        do 5000,ia=2,after
        ias=ia-1
        if (8*ias.eq.5*after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 5010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb        
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        do 5010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3) 
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r25 = r2 - r5
                        r34 = r3 + r4
                        s25 = s2 + s5
                        s34 = s3 - s4
                        zout(1,nout1,j) = r1 + r25 - r34
                        r = r1 + cos2*r25 - cos4*r34
                        s = sin2*s25 + sin4*s34
                        zout(1,nout2,j) = r - s
                        zout(1,nout5,j) = r + s
                        r = r1 + cos4*r25 - cos2*r34
                        s = sin4*s25 - sin2*s34
                        zout(1,nout3,j) = r - s
                        zout(1,nout4,j) = r + s
                        r25 = r2 + r5
                        r34 = r4 - r3
                        s25 = s2 - s5
                        s34 = s3 + s4
                        zout(2,nout1,j) = s1 + s25 + s34
                        r = s1 + cos2*s25 + cos4*s34
                        s = sin2*r25 + sin4*r34
                        zout(2,nout2,j) = r + s
                        zout(2,nout5,j) = r - s
                        r = s1 + cos4*s25 + cos2*s34
                        s = sin4*r25 - sin2*r34
                        zout(2,nout3,j) = r + s
                        zout(2,nout4,j) = r - s
5010                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 5020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb        
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        do 5020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=(r + s)*rt2i
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r25 = r2 - r5
                        r34 = r3 + r4
                        s25 = s2 + s5
                        s34 = s4 - s3
                        zout(1,nout1,j) = r1 + r25 + r34
                        r = r1 + cos2*r25 + cos4*r34
                        s = sin2*s25 + sin4*s34
                        zout(1,nout2,j) = r - s
                        zout(1,nout5,j) = r + s
                        r = r1 + cos4*r25 + cos2*r34
                        s = sin4*s25 - sin2*s34
                        zout(1,nout3,j) = r - s
                        zout(1,nout4,j) = r + s
                        r25 = r2 + r5
                        r34 = r3 - r4
                        s25 = s2 - s5
                        s34 = s3 + s4
                        zout(2,nout1,j) = s1 + s25 - s34
                        r = s1 + cos2*s25 - cos4*s34
                        s = sin2*r25 + sin4*r34
                        zout(2,nout2,j) = r + s
                        zout(2,nout5,j) = r - s
                        r = s1 + cos4*s25 - cos2*s34
                        s = sin4*r25 - sin2*r34
                        zout(2,nout3,j) = r + s
                        zout(2,nout4,j) = r - s
5020                        continue
                endif
        else
                ias=ia-1
                itt=ias*before
                itrig=itt+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                itrig=itrig+itt
                cr3=trig(1,itrig)
                ci3=trig(2,itrig)
                itrig=itrig+itt
                cr4=trig(1,itrig)
                ci4=trig(2,itrig)
                itrig=itrig+itt
                cr5=trig(1,itrig)
                ci5=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
                do 5100,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nin5=nin4+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                nout5=nout4+after
                do 5100,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                r=zin(1,j,nin3)
                s=zin(2,j,nin3)
                r3=r*cr3 - s*ci3
                s3=r*ci3 + s*cr3
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=r*cr4 - s*ci4
                s4=r*ci4 + s*cr4
                r=zin(1,j,nin5)
                s=zin(2,j,nin5)
                r5=r*cr5 - s*ci5
                s5=r*ci5 + s*cr5
                r25 = r2 + r5
                r34 = r3 + r4
                s25 = s2 - s5
                s34 = s3 - s4
                zout(1,nout1,j) = r1 + r25 + r34
                r = r1 + cos2*r25 + cos4*r34
                s = sin2*s25 + sin4*s34
                zout(1,nout2,j) = r - s
                zout(1,nout5,j) = r + s
                r = r1 + cos4*r25 + cos2*r34
                s = sin4*s25 - sin2*s34
                zout(1,nout3,j) = r - s
                zout(1,nout4,j) = r + s
                r25 = r2 - r5
                r34 = r3 - r4
                s25 = s2 + s5
                s34 = s3 + s4
                zout(2,nout1,j) = s1 + s25 + s34
                r = s1 + cos2*s25 + cos4*s34
                s = sin2*r25 + sin4*r34
                zout(2,nout2,j) = r + s
                zout(2,nout5,j) = r - s
                r = s1 + cos4*s25 + cos2*s34
                s = sin4*r25 - sin2*r34
                zout(2,nout3,j) = r + s
                zout(2,nout4,j) = r - s
5100                continue
        endif
5000        continue
       else if (now.eq.6) then
!         .5d0*sqrt(3.d0)
        bb=isign*0.8660254037844387d0

        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 6120,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        do 6120,j=1,nfft
        r2=zin(1,j,nin3)
        s2=zin(2,j,nin3)
        r3=zin(1,j,nin5)
        s3=zin(2,j,nin5)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        ur1 = r + r1
        ui1 = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r=r2-r3
        s=s2-s3
        ur2 = r1 - s*bb
        ui2 = s1 + r*bb
        ur3 = r1 + s*bb
        ui3 = s1 - r*bb

        r2=zin(1,j,nin6)
        s2=zin(2,j,nin6)
        r3=zin(1,j,nin2)
        s3=zin(2,j,nin2)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin4)
        s1=zin(2,j,nin4)
        vr1 = r + r1
        vi1 = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r=r2-r3
        s=s2-s3
        vr2 = r1 - s*bb
        vi2 = s1 + r*bb
        vr3 = r1 + s*bb
        vi3 = s1 - r*bb

        zout(1,nout1,j)=ur1+vr1
        zout(2,nout1,j)=ui1+vi1
        zout(1,nout5,j)=ur2+vr2
        zout(2,nout5,j)=ui2+vi2
        zout(1,nout3,j)=ur3+vr3
        zout(2,nout3,j)=ui3+vi3
        zout(1,nout4,j)=ur1-vr1
        zout(2,nout4,j)=ui1-vi1
        zout(1,nout2,j)=ur2-vr2
        zout(2,nout2,j)=ui2-vi2
        zout(1,nout6,j)=ur3-vr3
        zout(2,nout6,j)=ui3-vi3
6120        continue

       else
        stop 'error fftrot'
       endif

  end subroutine


  ! FFT PART RELATED TO THE CONVOLUTION--------------------------------------------
  
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
  !!	              most products of the prime factors 2,3,5 are allowed.
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
    integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc
    real(8), intent(in) :: scal,hgrid
    real(8), dimension(nd1,nd2,nd3/nproc), intent(in) :: pot
    real(8), dimension(md1,md3,md2/nproc), intent(inout) :: zf
    integer, intent(in) :: comm
    
#if defined(HAVE_MPI)
    integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
    integer :: j2,j3,i1,i3,i,j,inzee,ierr
    real(8) :: twopion
    !work arrays for transpositions
    real(8), allocatable :: zt(:,:,:)
    !work arrays for MPI
    real(8), allocatable :: zmpi1(:,:,:,:,:)
    real(8), allocatable :: zmpi2(:,:,:,:)
    !cache work array
    real(8), allocatable :: zw(:,:,:)
    !FFT work arrays
    real(8), allocatable :: cosinarr(:,:)
    real(8) :: btrig1(2,8192)
    real(8) :: ftrig1(2,8192)
    integer :: after1(7)
    integer :: now1(7)
    integer :: before1(7)
    real(8) :: btrig2(2,8192)
    real(8) :: ftrig2(2,8192)
    integer :: after2(7)
    integer :: now2(7)
    integer :: before2(7)
    real(8) :: btrig3(2,8192)
    real(8) :: ftrig3(2,8192)
    integer :: after3(7)
    integer :: now3(7)
    integer :: before3(7)
    type(profile_t), save :: prof, prof_comm

    call profiling_in(prof, "SG_PCONV")

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
    ncache = ncache_optimal()
    if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4

    ! if (timing_flag == 1 .and. iproc ==0) print *,'parallel ncache=',ncache

    lzt = n2/2
    if (mod(n2/2,2).eq.0) lzt=lzt+1
    if (mod(n2/2,4).eq.0) lzt=lzt+1

    SAFE_ALLOCATE(zw(1:2, 1:ncache/4, 1:2))
    SAFE_ALLOCATE(zt(1:2, 1:lzt, 1:n1))
    SAFE_ALLOCATE(zmpi2(1:2, 1:n1, 1:md2/nproc, 1:nd3))
    SAFE_ALLOCATE(cosinarr(1:2, 1:n3/2))

    if (nproc.gt.1) then
      SAFE_ALLOCATE(zmpi1(1:2, 1:n1, 1:md2/nproc, 1:nd3/nproc, 1:nproc))
    end if

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
    twopion=8.d0*datan(1.0_8)/real(n3,8)
    do i3=1,n3/2
      cosinarr(1,i3)=dcos(twopion*(i3-1))
      cosinarr(2,i3)=-dsin(twopion*(i3-1))
    end do

    !initializing integral
!!! ehartree=0.0_8

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
      call profiling_in(prof_comm, "SG_ALLTOALL")
      !communication scheduling
      call MPI_Alltoall(zmpi2(1, 1, 1, 1), n1*(md2/nproc)*(nd3/nproc), &
        MPI_DOUBLE_PRECISION, &
        zmpi1(1, 1, 1, 1, 1), n1*(md2/nproc)*(nd3/nproc), &
        MPI_DOUBLE_PRECISION, comm, ierr)
      call profiling_out(prof_comm)
    end if
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
      call profiling_in(prof_comm, "SG_ALLTOALL")
      !communication scheduling
      call MPI_ALLTOALL(zmpi1(1, 1, 1, 1, 1), n1*(md2/nproc)*(nd3/nproc), &
        MPI_DOUBLE_PRECISION, &
        zmpi2(1, 1, 1, 1), n1*(md2/nproc)*(nd3/nproc), &
        MPI_DOUBLE_PRECISION, comm, ierr)
      !output: I1,J2,j3,jp3,(Jp2)
      call profiling_out(prof_comm)
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
          call unfill_downcorn(md1, md3, lot, nfft, n3, zw(1,1,inzee), zf(i1,1,j2), scal)

        end do
      endif
    end do

    SAFE_DEALLOCATE_A(zmpi2)
    SAFE_DEALLOCATE_A(zw)
    SAFE_DEALLOCATE_A(zt)
    SAFE_DEALLOCATE_A(cosinarr)
    SAFE_DEALLOCATE_A(zmpi1)

    call profiling_out(prof)

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
    integer, intent(in) :: nd1,nd2,n1,n2,lot,nfft,jS
    real(8), dimension(nd1,nd2), intent(in) :: pot
    real(8), dimension(2,lot,n2), intent(inout) :: zw

    !Local variables
    integer :: j,j1,i2,j2

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
    integer :: lot, n1, n2, lzt, j, nfft, i
    real(8) :: zw(2,lot,n2), zt(2,lzt,n1)

    ! WARNING: Assuming that high frequencies are in the corners 
    !          and that n2 is multiple of 2
    
    ! Low frequencies 
    do j = 1, nfft
      do i = n2/2 + 1, n2
	zw(1,j,i)=zt(1,i-n2/2,j)
        zw(2,j,i)=zt(2,i-n2/2,j)
      end do
    end do
 
    ! High frequencies 
    do i = 1, n2/2
      do j = 1, nfft
        zw(1,j,i)=0.0_8
        zw(2,j,i)=0.0_8
      end do
    end do
     
  end subroutine switch_upcorn

  ! ----------------------------------------------------------------------------

  subroutine mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
    integer :: n1, nd2, nd3, nproc, lot, mfft, jp2, jp2stb, j2
    integer :: j2stb, nfft, i1, j3, md2
    real(8) :: zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
    
    ! WARNING: Assuming that high frequencies are in the corners 
    !          and that n1 is multiple of 2
    
    mfft=0
    do Jp2 = Jp2stb, nproc
      do J2 = J2stb, md2/nproc
	mfft = mfft + 1
        if (mfft.gt.nfft) then
          Jp2stb=Jp2
          J2stb=J2
          return
        end if
        do I1=1,n1/2
          zw(1,mfft,I1)=0.0_8
          zw(2,mfft,I1)=0.0_8
        end do
        do I1=n1/2+1,n1
          zw(1,mfft,I1)=zmpi1(1,I1-n1/2,J2,j3,Jp2)
          zw(2,mfft,I1)=zmpi1(2,I1-n1/2,J2,j3,Jp2)
        end do
      end do
      J2stb=1
    end do

  end subroutine mpiswitch_upcorn

  ! -------------------------------------------------------------------

  subroutine halfill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
    integer :: lot, n3, md1, md3, i3, i1, nfft
    real(8) :: zw(2,lot,n3/2),zf(md1,md3)
    
    ! WARNING: Assuming that high frequencies are in the corners 
    !          and that n3 is multiple of 4
    !in principle we can relax this condition
    
    do i3 = 1, n3/4
      do i1=1,nfft
	zw(1,i1,i3)=0.0_8
        zw(2,i1,i3)=0.0_8
      end do
    end do

    do i3 = n3/4+1,n3/2
      do i1 = 1,nfft
        zw(1,i1,i3)=zf(i1,2*i3-1-n3/2)
	zw(2,i1,i3)=zf(i1,2*i3-n3/2)
      end do
    end do

  end subroutine halfill_upcorn

  ! -------------------------------------------------------------------
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
    integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
    real(8), dimension(2,lot,n3/2), intent(in) :: zw
    real(8), dimension(2,n3/2), intent(in) :: cosinarr
    real(8), dimension(2,n1,md2/nproc,nd3), intent(out) :: zmpi2

    integer :: i3,i,ind1,ind2
    real(8) ::  a,b,c,d,cp,sp,feR,feI,foR,foI,fR,fI

    !case i3=1 and i3=n3/2+1
    do i=0,nfft-1
      a=zw(1,i+1,1)
      b=zw(2,i+1,1)
      zmpi2(1,i1+i,j2,1)=a+b
      zmpi2(2,i1+i,j2,1)=0.0_8
      zmpi2(1,i1+i,j2,n3/2+1)=a-b
      zmpi2(2,i1+i,j2,n3/2+1)=0.0_8
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
    integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
    real(8), dimension(2,lot,n3/2), intent(out) :: zw
    real(8), dimension(2,n3/2), intent(in) :: cosinarr
    real(8), dimension(2,n1,md2/nproc,nd3), intent(in) :: zmpi2

    integer :: i3,i,indA,indB
    real(8) ::  a,b,c,d,cp,sp,re,ie,ro,io,rh,ih

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

  ! --------------------------------------------------------------------

  subroutine unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
    integer :: lot, n2, lzt, n1, j, nfft, i
    real(8) :: zw(2,lot,n2),zt(2,lzt,n1)
    ! WARNING: Assuming that high frequencies are in the corners 
    !          and that n2 is multiple of 2
    
    ! Low frequencies
    do j = 1, nfft
      do i = 1, n2/2
	zt(1,i,j)=zw(1,j,i)
        zt(2,i,j)=zw(2,j,i)
      end do
    end do
    
  end subroutine unswitch_downcorn


  subroutine unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
    integer :: n1, md2, nproc, nd3, lot, mfft, i1, j3
    integer :: jp2, j2, nfft, jp2stf, j2stf
    real(8) :: zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
    ! WARNING: Assuming that high frequencies are in the corners 
    !          and that n1 is multiple of 2

    mfft=0
    do Jp2=Jp2stf,nproc
      do J2 = J2stf, md2/nproc
	mfft=mfft+1
        if (mfft.gt.nfft) then
          Jp2stf=Jp2
          J2stf=J2
          return
        endif
        do I1 = 1, n1/2
          zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
          zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
      end do
      J2stf=1
    end do
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
  subroutine unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf, scal)
    integer, intent(in) :: md1,md3,lot,nfft,n3
    real(8), dimension(2,lot,n3/2), intent(in) :: zw
    real(8), dimension(md1,md3),intent(inout) :: zf
    real(8), intent(in) :: scal

    integer :: i3,i1
    real(8) :: pot1

    do i3=1,n3/4
      do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        zf(i1, 2*i3 - 1) = pot1 
        pot1 = scal*zw(2, i1, i3)
        zf(i1, 2*i3) = pot1 
      end do
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
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,nproc,iproc
  real(8), dimension(nd1,n3,nd2/nproc), intent(in) :: zf
  real(8), dimension(2,nd1,nd2,nd3/nproc), intent(out) :: zr
  integer, intent(in) :: comm

#if defined(HAVE_MPI)
  !Local variables
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2st,J2st
  integer :: j2,j3,i1,i3,i,j,inzee,ierr
  real(8) :: twopion
  !work arrays for transpositions
  real(8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(8), dimension(:,:), allocatable :: trig1,trig2,trig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3

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
  twopion=8.d0*datan(1.0_8)/real(n3,8)
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

  ! ---------------------------------------------------------------

  subroutine switch(nfft,n2,lot,n1,lzt,zt,zw)
    integer :: lot, n2, lzt, n1, j, nfft, i
    real(8) :: zw(1:2, 1:lot, 1:n2), zt(1:2, 1:lzt, 1:n1)
    
    do j = 1, nfft
      do i = 1, n2
        zw(1,j,i) = zt(1,i,j)
        zw(2,j,i) = zt(2,i,j)
      end do
    end do
  end subroutine switch

  ! ---------------------------------------------------------------

  subroutine mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw)
    integer :: n1, nd2, nproc, nd3, lot, mfft, jp2, jp2st
    integer :: j2st, nfft, i1, j3, j2
    real(8) :: zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
    
    mfft=0
    do Jp2 = Jp2st, nproc
      do J2 = J2st, nd2/nproc
        mfft = mfft + 1
        if (mfft .gt. nfft) then
          Jp2st = Jp2
          J2st = J2
          return
        endif
        do I1 = 1, n1
          zw(1,mfft,I1) = zmpi1(1, I1, J2, j3, Jp2)
          zw(2,mfft,I1) = zmpi1(2, I1, J2, j3, Jp2)
        end do
      end do
      J2st=1
    end do
  end subroutine mpiswitch

  ! ---------------------------------------------------------------

  subroutine inserthalf(nd1,lot,nfft,n3,zf,zw)
    integer :: lot, n3, nd1, i3, i1, nfft
    real(8) :: zw(1:2, 1:lot, 1:n3/2), zf(1:nd1, 1:n3)
    
    do i3 = 1, n3/2
      do i1 = 1, nfft
        zw(1, i1, i3) = zf(i1, 2*i3 - 1)
        zw(2, i1, i3) = zf(i1, 2*i3)
      end do
    end do
        
  end subroutine inserthalf

end module sgfft_m
