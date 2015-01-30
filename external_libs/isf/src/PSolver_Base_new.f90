!> @file
!!  New routines for Poisson solver
!! @author
!! Copyright (C) 2002-2011 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!> Parallel version of Poisson Solver
!! General version, for each boundary condition
subroutine G_PoissonSolver(iproc,nproc,planes_comm,iproc_inplane,inplane_comm,geocode,ncplx,&
     n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,pot,zf,scal,hx,hy,hz,offset,strten)
  use Poisson_Solver, only: dp, gp, TCAT_PSOLV_COMMUN,TCAT_PSOLV_COMPUT
  use wrapper_mpi
  !use memory_profiling
  use time_profiling, only: f_timing
  use dynamic_memory
  implicit none
  !to be preprocessed
  include 'perfdata.inc'
  !Arguments
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(inout) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,ncplx
  integer, intent(in) :: planes_comm,inplane_comm,iproc_inplane
  real(gp), intent(in) :: scal,hx,hy,hz,offset
  real(dp), dimension(nd1,nd2,nd3/nproc), intent(in) :: pot
  real(dp), dimension(ncplx,md1,md3,md2/nproc), intent(inout) :: zf
  real(dp), dimension(6), intent(out) :: strten !< non-symmetric components of Ha stress tensor
  !Local variables
  character(len=*), parameter :: subname='G_Poisson_Solver'
  logical :: perx,pery,perz,halffty,cplx
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer :: ncache,lzt,lot,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  integer :: j2,j3,i1,i3,i,j,inzee,ierr,n1dim,n2dim,n3dim,ntrig, i_stat
  !integer :: i_all
  real(kind=8) :: twopion
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  real(kind=8), dimension(:,:,:), allocatable :: zt_t
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: btrig1,btrig2,btrig3
  real(kind=8), dimension(:,:), allocatable :: ftrig1,ftrig2,ftrig3,cosinarr
  integer, dimension(7) :: after1,now1,before1
  integer, dimension(7) :: after2,now2,before2,after3,now3,before3
  real(gp), dimension(6) :: strten_omp
  !integer :: ncount0,ncount1,ncount_max,ncount_rate

  integer :: maxIter
  integer :: n3pr1,n3pr2,j1start,n1p,n2dimp

  call f_routine(id='G_PoissonSolver')

  !Initialize stress tensor no matter of the BC
  !call to_zero(6,strten(1))

  strten=0.d0

  !conditions for periodicity in the three directions
  !perx=(geocode /= 'F' .and. geocode /= 'W' .and. geocode /= 'H')
  perx=(geocode == 'P' .or. geocode == 'S')
  pery=(geocode == 'P')
  perz=(geocode /= 'F' .and. geocode /= 'H')

  cplx= (ncplx == 2)

  !also for complex input this should be eliminated
  halffty=.not. pery .and. .not. cplx

  !defining work arrays dimensions
  ncache=ncache_optimal

  !n3/2 if the dimension is real and isolated
  if (halffty) then
     n3dim=n3/2
  else
     n3dim=n3
  end if

  if (perz) then
     n2dim=n2
  else
     n2dim=n2/2
  end if

  if (perx) then
     n1dim=n1
  else
     n1dim=n1/2
  end if
  call f_timing(TCAT_PSOLV_COMPUT,'ON')
  !call timing(iproc,'PSolv_comput  ','ON')
  ! check input
  if (mod(n1,2) /= 0 .and. .not. perx) stop 'Parallel convolution:ERROR:n1' !this can be avoided
  if (mod(n2,2) /= 0 .and. .not. perz) stop 'Parallel convolution:ERROR:n2' !this can be avoided
  if (mod(n3,2) /= 0 .and. .not. pery) stop 'Parallel convolution:ERROR:n3' !this can be avoided
  if (nd1 < n1/2+1) stop 'Parallel convolution:ERROR:nd1' 
  if (nd2 < n2/2+1) stop 'Parallel convolution:ERROR:nd2' 
  if (nd3 < n3/2+1) stop 'Parallel convolution:ERROR:nd3' 
  if (md1 < n1dim) stop 'Parallel convolution:ERROR:md1'
  if (md2 < n2dim) stop 'Parallel convolution:ERROR:md2'
  if (md3 < n3dim) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,nproc) /= 0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,nproc) /= 0) stop 'Parallel convolution:ERROR:md2'
  
  if (ncache <= max(n1,n2,n3dim)*4) ncache=max(n1,n2,n3dim)*4

  if (timing_flag == 1 .and. iproc ==0) print *,'parallel ncache=',ncache

  
  ntrig=max(n3dim,n1,n2)

  n1p=n1

  if (nproc>2*(n3/2+1)-1 .and. .false.) then
    n3pr1=nproc/(n3/2+1)
    n3pr2=n3/2+1
    if (mod(n1,n3pr1) .ne. 0) n1p=((n1/n3pr1)+1)*n3pr1
  else
    n3pr1=1
    n3pr2=nproc
  endif

  if (n3pr1>1) then
   lzt=(md2/nproc)*n3pr1*n3pr2
   n2dimp=lzt
  else
   lzt=n2dim
   if (mod(n2dim,2) == 0) lzt=lzt+1
   if (mod(n2dim,4) == 0) lzt=lzt+1 !maybe this is useless
   n2dimp=n2dim
  endif

  !if (iproc==0) print*,'psolver param',md1,md2,md3,nd1,nd2,nd3,n1dim,n2dim,n3dim,n1,n2,n3

  if (lzt < n2dim) then
    if (iproc==0) print*,'PSolver error: lzt < n2dim',lzt,n2dim
    call mpi_finalize(i_stat)
    stop
  endif
  !if (iproc==0) print*,'lzt n2dim',lzt,n2dim


  !Allocations
  btrig1 = f_malloc((/ 2, ntrig /),id='btrig1')
  ftrig1 = f_malloc((/ 2, ntrig /),id='ftrig1')
!  after1 = f_malloc(7,id='after1')
!  now1 = f_malloc(7,id='now1')
!  before1 = f_malloc(7,id='before1')
  btrig2 = f_malloc((/ 2, ntrig /),id='btrig2')
  ftrig2 = f_malloc((/ 2, ntrig /),id='ftrig2')
!  after2 = f_malloc(7,id='after2')
!  now2 = f_malloc(7,id='now2')
!  before2 = f_malloc(7,id='before2')
  btrig3 = f_malloc((/ 2, ntrig /),id='btrig3')
  ftrig3 = f_malloc((/ 2, ntrig /),id='ftrig3')
!  after3 = f_malloc(7,id='after3')
!  now3 = f_malloc(7,id='now3')
!  before3 = f_malloc(7,id='before3')
  !allocate(zw(2,ncache/4,2+ndebug),stat=i_stat)
  !call memocc(i_stat,zw,'zw',subname)
  !allocate(zt(2,lzt,n1+ndebug),stat=i_stat)
  !call memocc(i_stat,zt,'zt',subname)
  zmpi2 = f_malloc((/ 2, n1, md2/nproc, nd3 /),id='zmpi2')
  !also for complex input this should be eliminated
  if (halffty) then
     cosinarr = f_malloc((/ 2, n3/2 /),id='cosinarr')
  end if
  if (nproc > 1) then
     zmpi1 = f_malloc((/ 2, n1, md2/nproc, nd3/nproc, n3pr2 /),id='zmpi1')
  end if

  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)

  !$omp parallel sections default(shared)
  !$omp section
    call ctrig_sg(n3dim,ntrig,btrig3,after3,before3,now3,1,ic3)
    do j = 1, n3dim
      ftrig3(1, j) = btrig3(1, j)
      ftrig3(2, j) = -btrig3(2, j)
    enddo
  !$omp section
    call ctrig_sg(n1,ntrig,btrig1,after1,before1,now1,1,ic1)
    do j = 1, n1
      ftrig1(1, j) = btrig1(1, j)
      ftrig1(2, j) = -btrig1(2, j)
    enddo
  !$omp section
    call ctrig_sg(n2,ntrig,btrig2,after2,before2,now2,1,ic2)
    do j = 1, n2
      ftrig2(1, j) = btrig2(1, j)
      ftrig2(2, j) = -btrig2(2, j)
    enddo
  !$omp section
    if (halffty) then
      !Calculating array of phases for HalFFT decoding
      twopion=8.d0*datan(1.d0)/real(n3,kind=8)
      do i3=1,n3/2
        cosinarr(1,i3)= dcos(twopion*real(i3-1,kind=8))
        cosinarr(2,i3)=-dsin(twopion*real(i3-1,kind=8))
      end do
    end if
  !$omp end parallel sections

  ! transform along z axis
  lot=ncache/(4*n3dim)
  if (lot < 1) then  
     write(6,*) & 
          'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
          'least one 1-d FFT of this size even though this will' // & 
          'reduce the performance for shorter transform lengths'
     stop
  endif


  !put to zero the zw array
  !this should not be done at each time since the
  !array is refilled always the same way
  !zw=0.0_dp
  !call razero(4*(ncache/4),zw)
  !different loop if halfft or not (output part)

  maxIter = min(md2 /nproc, n2dim - iproc *(md2 /nproc))

  if (n3pr1 > 1) zt_t = f_malloc((/ 2, lzt/n3pr1, n1p /),id='zt_t')

!$omp parallel default(shared)&
  !$omp private(nfft,inzee,Jp2stb,J2stb,Jp2stf,J2stf,i3,strten_omp, zw, zt) &
  !$omp private(j2,i1,i,j3,j) &
  !$omp firstprivate(lot, maxIter)
!  !$omp firstprivate(before3, now3, after3)
  
  allocate( zw(2, ncache/4, 2), stat=i_stat )
  allocate( zt(2,lzt/n3pr1, n1p), stat=i_stat )
  !$omp do schedule(static)
  do j2 = 1, maxIter
     !this condition ensures that we manage only the interesting part for the FFT
     !if (iproc*(md2/nproc)+j2 <= n2dim) then

        do i1=1,n1dim,lot
           nfft = min(i1 + (lot -1), n1dim) - i1 +1

            if (halffty) then
              !inserting real data into complex array of half lenght
              call halfill_upcorn(md1,md3,lot,nfft,n3,zf(1,i1,1,j2),zw(1,1,1))
            else if (cplx) then
              !zf should have four indices
              call C_fill_upcorn(md1,md3,lot,nfft,n3,zf(1,i1,1,j2),zw(1,1,1))
            else
              call P_fill_upcorn(md1,md3,lot,nfft,n3,zf(1,i1,1,j2),zw(1,1,1))
            end if
           !performing FFT

           !input: I1,I3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp_sg(lot,nfft,n3dim,lot,n3dim,zw(1,1,inzee), &
                zw(1,1,3-inzee),ntrig,btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo

           !output: I1,i3,J2,(Jp2)
           !exchanging components
           !input: I1,i3,J2,(Jp2)
           if (halffty) then
              call scramble_unpack(i1,j2,lot,nfft,n1dim,n3,md2,nproc,nd3,&
                   zw(1,1,inzee),zmpi2,cosinarr)
           else
              call scramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2)
           end if
           !output: I1,J2,i3,(Jp2)
        end do
     !end if
  end do
  !$omp end do
  ! DO NOT USE NOWAIT, removes the implicit barrier

  !$omp master
    nd3=nd3/n3pr1
    !Interprocessor data transposition
    !input: I1,J2,j3,jp3,(Jp2)
    if (nproc > 1 .and. iproc < n3pr1*n3pr2) then
       call f_timing(TCAT_PSOLV_COMPUT,'OF')
       !call timing(iproc,'PSolv_comput  ','OF')
       !call timing(iproc,'PSolv_commun  ','ON')
       call f_timing(TCAT_PSOLV_COMMUN,'ON')
       !communication scheduling
       call MPI_ALLTOALL(zmpi2,2*n1dim*(md2/nproc)*(nd3/n3pr2), &
            MPI_double_precision, &
            zmpi1,2*n1dim*(md2/nproc)*(nd3/n3pr2), &
            MPI_double_precision,planes_comm,ierr)
       call f_timing(TCAT_PSOLV_COMMUN,'OF')
        !call timing(iproc,'PSolv_commun  ','OF')
       !call timing(iproc,'PSolv_comput  ','ON')
       call f_timing(TCAT_PSOLV_COMPUT,'ON')
    endif

    !output: I1,J2,j3,Jp2,(jp3)
  !$omp end master
  !$omp barrier

  !now each process perform complete convolution of its planes
 
  if (n3pr1==1) then
    maxIter = min(nd3 /nproc, n3/2+1 - iproc*(nd3/nproc))
  else
    if (iproc < n3pr1*n3pr2) then
        maxIter = nd3/n3pr2
    else
        maxIter = 0
    endif
  endif
 
  strten_omp=0

  !$omp do schedule(static)
  do j3 = 1, maxIter
       !this condition ensures that we manage only the interesting part for the FFT
     !if (iproc*(nd3/nproc)+j3 <= n3/2+1) then
          Jp2stb=1
          J2stb=1
          Jp2stf=1
          J2stf=1
          ! transform along x axis
          lot=ncache/(4*n1)
          if (lot < 1) then  
             write(6,*)&
                  'convolxc_off:ncache has to be enlarged to be able to hold at' //&
                  'least one 1-d FFT of this size even though this will' //&
                  'reduce the performance for shorter transform lengths'
             stop
          endif

          do j=1,n2dimp/n3pr1,lot
           nfft=min(j+(lot-1), n2dimp/n3pr1) -j +1

             !reverse index ordering, leaving the planes to be transformed at the end
             !input: I1,J2,j3,Jp2,(jp3)
             if (nproc > 1) then
                call G_mpiswitch_upcorn2(j3,nfft,Jp2stb,J2stb,lot,&
                     n1,n1dim,md2,nd3,n3pr1,n3pr2,zmpi1,zw(1,1,1))
             else
                call G_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,&
                     n1,n1dim,md2,nd3,nproc,zmpi2,zw(1,1,1))
             endif

             !output: J2,Jp2,I1,j3,(jp3)
             !performing FFT
             !input: I2,I1,j3,(jp3)
             inzee=1
             do i=1,ic1-1
                call fftstp_sg(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee),&
                     ntrig,btrig1,after1(i),now1(i),before1(i),1)
                inzee=3-inzee
             enddo

             !storing the last step into zt array
             i=ic1
             call fftstp_sg(lot,nfft,n1,lzt/n3pr1,n1,zw(1,1,inzee),zt(1,j,1),&
                  ntrig,btrig1,after1(i),now1(i),before1(i),1)           
             !output: I2,i1,j3,(jp3)
          end do

          !transform along y axis
          lot=ncache/(4*n2)
          if (lot < 1) then  
             write(6,*)&
                  'convolxc_off:ncache has to be enlarged to be able to hold at' //&
                  'least one 1-d FFT of this size even though this will' //&
                  'reduce the performance for shorter transform lengths'
             stop
          endif

          !LG: I am not convinced that putting MPI_ALLTOALL inside a loop
          !!   is a good idea: this will break OMP parallelization, moreover
          !!   the granularity of the external loop will impose too many communications
          if (n3pr1 > 1 .and. inplane_comm/=MPI_COMM_NULL) then
            call f_timing(TCAT_PSOLV_COMPUT,'OF')
            call f_timing(TCAT_PSOLV_COMMUN,'ON')

            call MPI_ALLTOALL(zt,2*(n1p/n3pr1)*(lzt/n3pr1),MPI_double_precision,zt_t,2*(n1p/n3pr1)*(lzt/n3pr1), &
                            MPI_double_precision,inplane_comm,ierr)

            call f_timing(TCAT_PSOLV_COMMUN,'OF')
            call f_timing(TCAT_PSOLV_COMPUT,'ON')
          endif
        
          do j=1,n1p/n3pr1,lot
           nfft=min(j+(lot-1),n1p/n3pr1)-j+1
             !reverse ordering 
             !input: I2,i1,j3,(jp3)
             if (n3pr1 >1) then
               call G_switch_upcorn2(nfft,n2,n2dim,lot,n1p,lzt,zt_t(1,1,1),zw(1,1,1),n3pr1,j)
             else
               call G_switch_upcorn(nfft,n2,n2dim,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
             endif
             !output: i1,I2,j3,(jp3)
             !performing FFT
             !input: i1,I2,j3,(jp3)
             inzee=1

             do i=1,ic2
                call fftstp_sg(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee),&
                     ntrig,btrig2,after2(i),now2(i),before2(i),1)
                inzee=3-inzee
             enddo
             !output: i1,i2,j3,(jp3)
             !Multiply with kernel in fourier space
             i3=mod(iproc,n3pr2)*(nd3/n3pr2)+j3

            j1start=0
            if (n3pr1>1) j1start=(n1p/n3pr1)*iproc_inplane

            if (geocode == 'P') then
              call P_multkernel(nd1,nd2,n1,n2,n3,lot,nfft,j+j1start,pot(1,1,j3),zw(1,1,inzee),&
                   i3,hx,hy,hz,offset,scal,strten_omp)
             else
                !write(*,*) 'pot(1,1,j3) = ', pot(1,1,j3)
                call multkernel(nd1,nd2,n1,n2,lot,nfft,j+j1start,pot(1,1,j3),zw(1,1,inzee))
             end if
            
!TRANSFORM BACK IN REAL SPACE
             !transform along y axis
             !input: i1,i2,j3,(jp3)
             do i=1,ic2
                call fftstp_sg(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee),&
                     ntrig,ftrig2,after2(i),now2(i),before2(i),-1)
               !zw(:,:,3-inzee)=zw(:,:,inzee)
                inzee=3-inzee
             end do

             !reverse ordering
             !input: i1,I2,j3,(jp3)
           if (n3pr1 == 1) then
             call G_unswitch_downcorn(nfft,n2,n2dim,lot,n1,lzt, &
               zw(1,1,inzee),zt(1,1,j))
           else
             call G_unswitch_downcorn2(nfft,n2,n2dim,lot,n1p,lzt, &
               zw(1,1,inzee),zt_t(1,1,1),n3pr1,j)
           endif
             !output: I2,i1,j3,(jp3)
        end do
          !transform along x axis
          !input: I2,i1,j3,(jp3)

        !LG: this MPI_ALLTOALL is inside a loop. I think that it will rapidly become unoptimal
        if (n3pr1 > 1 .and. inplane_comm/=MPI_COMM_NULL) then
            call f_timing(TCAT_PSOLV_COMPUT,'OF')
            call f_timing(TCAT_PSOLV_COMMUN,'ON')

            call MPI_ALLTOALL(zt_t,2*(n1p/n3pr1)*(lzt/n3pr1),MPI_double_precision,zt,2*(n1p/n3pr1)*(lzt/n3pr1), &
                            MPI_double_precision,inplane_comm,ierr)

            call f_timing(TCAT_PSOLV_COMMUN,'OF')
            call f_timing(TCAT_PSOLV_COMPUT,'ON')
          endif

          lot=ncache/(4*n1)
          do j=1,n2dimp/n3pr1,lot
           nfft=min(j+(lot-1),n2dimp/n3pr1)-j+1

            !performing FFT
             i=1
             if (n3pr1 > 1) then
               call fftstp_sg(lzt/n3pr1,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1),&
                  ntrig,ftrig1,after1(i),now1(i),before1(i),-1)
             else
               call fftstp_sg(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1),&
                  ntrig,ftrig1,after1(i),now1(i),before1(i),-1)
             endif
             inzee=1
             do i=2,ic1
                call fftstp_sg(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee),&
                     ntrig,ftrig1,after1(i),now1(i),before1(i),-1)
                inzee=3-inzee
             enddo
             
             !output: I2,I1,j3,(jp3)
             !reverse ordering
             !input: J2,Jp2,I1,j3,(jp3)
             if (nproc == 1) then
                call G_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,&
                     n1dim,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
             else
                call G_unmpiswitch_downcorn2(j3,nfft,Jp2stf,J2stf,lot,n1,&
                     n1dim,md2,nd3,n3pr1,n3pr2,zw(1,1,inzee),zmpi1)  
             endif
             ! output: I1,J2,j3,Jp2,(jp3)
          end do
     !endif

!END OF TRANSFORM FOR X AND Z

    end do
  !$omp end do
    !$omp critical
    !do i = 1, 6
    !  strten(j) = strten(j) + strten_omp(j)
    !enddo
    strten = strten + strten_omp
    !$omp end critical
  ! DO NOT USE NOWAIT, removes the implicit barrier

  !$omp master

!TRANSFORM BACK IN Y
    !Interprocessor data transposition
    !input: I1,J2,j3,Jp2,(jp3)
    if (nproc > 1 .and. iproc < n3pr1*n3pr2) then
       call f_timing(TCAT_PSOLV_COMPUT,'OF')
       call f_timing(TCAT_PSOLV_COMMUN,'ON')
       !communication scheduling
       call MPI_ALLTOALL(zmpi1,2*n1dim*(md2/nproc)*(nd3/n3pr2), &
            MPI_double_precision, &
            zmpi2,2*n1dim*(md2/nproc)*(nd3/n3pr2), &
            MPI_double_precision,planes_comm,ierr)
       call f_timing(TCAT_PSOLV_COMMUN,'OF')
       call f_timing(TCAT_PSOLV_COMPUT,'ON')
    endif
    !output: I1,J2,j3,jp3,(Jp2)
    nd3=nd3*n3pr1
  !$omp end master
  !$omp barrier

  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(4*n3dim)

  maxIter = min(md2/nproc, n2dim - iproc *(md2/nproc))

  !$omp do schedule(static)
  do j2 = 1, maxIter
     !this condition ensures that we manage only the interesting part for the FFT
        do i1=1,n1dim,lot
           nfft=min(i1+(lot-1),n1dim)-i1+1

           !reverse ordering
           !input: I1,J2,i3,(Jp2)
           if (halffty) then
              call unscramble_pack(i1,j2,lot,nfft,n1dim,n3,md2,nproc,nd3,zmpi2, &
                zw(1,1,1),cosinarr)
           else
              call unscramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw(1,1,1))
           end if
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)           
           inzee=1
           do i=1,ic3
              call fftstp_sg(lot,nfft,n3dim,lot,n3dim,zw(1,1,inzee), &
                zw(1,1,3-inzee),ntrig,ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !rebuild the output array

            if (halffty) then
              call unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(1,i1,1,j2)&
                   ,scal)!,ehartreetmp)
            else if (cplx) then
              call C_unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(1,i1,1,j2),scal)
            else
              call P_unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(1,i1,1,j2),scal)
            end if

           !integrate local pieces together
           !ehartree=ehartree+0.5d0*ehartreetmp*hx*hy*hz

        end do
  end do
  !$omp end do

  deallocate(zw, stat=i_stat)
  deallocate(zt, stat=i_stat)

  !$omp end parallel

  if (n3pr1 > 1) call f_free(zt_t)

!END OF TRANSFORM IN Y DIRECTION

  !De-allocations  
  call f_free(btrig1)
  call f_free(ftrig1)
!  call f_free(after1)
!  call f_free(now1)
!  call f_free(before1)
  call f_free(btrig2)
  call f_free(ftrig2)
!  call f_free(after2)
!  call f_free(now2)
!  call f_free(before2)
  call f_free(btrig3)
  call f_free(ftrig3)
!  call f_free(after3)
!  call f_free(now3)
!  call f_free(before3)
  call f_free(zmpi2)
  !i_all=-product(shape(zw))*kind(zw)
  !deallocate(zw,stat=i_stat)
  !call memocc(i_stat,i_all,'zw',subname)
  !i_all=-product(shape(zt))*kind(zt)
  !deallocate(zt,stat=i_stat)
  !call memocc(i_stat,i_all,'zt',subname)
  if (halffty) then
     call f_free(cosinarr)
  end if

  if (nproc > 1) then
     call f_free(zmpi1)
  end if
  call f_timing(TCAT_PSOLV_COMPUT,'OF')
  call f_release_routine()
  !call system_clock(ncount1,ncount_rate,ncount_max)
  !write(*,*) 'TIMING:PS ', real(ncount1-ncount0)/real(ncount_rate)
END SUBROUTINE G_PoissonSolver


!> General routine, takes into account the free boundary conditions
subroutine G_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,&
     n1,n1dim,md2,nd3,nproc,zmpi1,zw)
  use Poisson_Solver, only: dp
  implicit none
  !Arguments
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc,n1dim
  integer, intent(inout) :: Jp2stb,J2stb
  real(dp),intent(inout) ::  zmpi1(2,n1dim,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
  !Local variables
  integer :: mfft,Jp2,J2,I1,ish, imfft

  
  !shift
  ish=n1-n1dim
  mfft=0

  loop_Jp2: do Jp2=Jp2stb,nproc
     do J2=J2stb,md2/nproc
        mfft=mfft+1
        if (mfft > nfft) then
           Jp2stb=Jp2
           J2stb=J2
           !return
           exit loop_Jp2
        end if
        do I1=1,ish
           zw(1,mfft,I1)=0.0_dp
           zw(2,mfft,I1)=0.0_dp
        end do
        do I1=1,n1dim
           zw(1,mfft,I1+ish)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1+ish)=zmpi1(2,I1,J2,j3,Jp2)
        end do
     end do
     J2stb=1
  end do loop_Jp2

  do I1 = 1, ish
     do imfft = 1, mfft-1
        zw(1, imfft, I1) = 0.0_dp
        zw(2, imfft, I1) = 0.0_dp
     end do
  end do

END SUBROUTINE G_mpiswitch_upcorn

subroutine G_mpiswitch_upcorn2(j3,nfft,Jp2stb,J2stb,lot,&
     n1,n1dim,md2,nd3,n3pr1,n3pr2,zmpi1,zw)
  !use module_base
  use Poisson_Solver, only: dp
  implicit none
  !Arguments
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,n3pr1,n3pr2,n1dim
  integer, intent(inout) :: Jp2stb,J2stb
  real(dp),intent(inout) ::  zmpi1(2,n1dim,md2/(n3pr1*n3pr2),nd3/n3pr2,n3pr2),zw(2,lot,n1)
  !Local variables
  integer :: mfft,Jp2,J2,I1,ish, imfft

  !shift
  ish=n1-n1dim
  mfft=0

  loop_Jp2: do Jp2=Jp2stb,n3pr2
     do J2=J2stb,md2/(n3pr1*n3pr2)
        mfft=mfft+1
        if (mfft > nfft) then
           Jp2stb=Jp2
           J2stb=J2
           !return
           exit loop_Jp2
        end if
        do I1=1,ish
           zw(1,mfft,I1)=0.0_dp
           zw(2,mfft,I1)=0.0_dp
        end do
        do I1=1,n1dim
           zw(1,mfft,I1+ish)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1+ish)=zmpi1(2,I1,J2,j3,Jp2)
        end do
     end do
     J2stb=1
  end do loop_Jp2

  do I1 = 1, ish
     do imfft = 1, mfft-1
        zw(1, imfft, I1) = 0.0_dp
        zw(2, imfft, I1) = 0.0_dp
     end do
  end do

END SUBROUTINE G_mpiswitch_upcorn2


subroutine G_switch_upcorn(nfft,n2,n2dim,lot,n1,lzt,zt,zw)
  use Poisson_Solver, only: dp
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt,n2dim
  real(dp), dimension(2,lzt,n1), intent(in) :: zt
  real(dp), dimension(2,lot,n2), intent(inout) :: zw
  !local variables
  integer :: i,j,ish

  !shift
  ish=n2-n2dim
  ! Low frequencies 
  do j=1,nfft
     do i=1,n2dim
        zw(1,j,i+ish)=zt(1,i,j)
        zw(2,j,i+ish)=zt(2,i,j)
     end do
  end do

  ! High frequencies 
  do i=1,ish
     do j=1,nfft
        zw(1,j,i)=0.0_dp
        zw(2,j,i)=0.0_dp
     end do
  end do

END SUBROUTINE G_switch_upcorn


subroutine G_switch_upcorn2(nfft,n2,n2dim,lot,n1,lzt,zt,zw,n3pr1,zt_init)
  use Poisson_Solver, only: dp
  !use module_base
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt,n2dim,n3pr1,zt_init
  real(dp), dimension(2,lzt/n3pr1,n1/n3pr1,n3pr1), intent(in) :: zt
  real(dp), dimension(2,lot,n2), intent(inout) :: zw
  !local variables
  integer :: i,j,ish

  !shift
  ish=n2-n2dim
  ! Low frequencies 
  do j=1,nfft
     do i=1,n2dim
        zw(1,j,i+ish)=zt(1,mod(i-1,lzt/n3pr1)+1,j+zt_init-1,(i-1)/(lzt/n3pr1)+1)
        zw(2,j,i+ish)=zt(2,mod(i-1,lzt/n3pr1)+1,j+zt_init-1,(i-1)/(lzt/n3pr1)+1)
     end do
  end do

  ! High frequencies 
  do i=1,ish
     do j=1,nfft
        zw(1,j,i)=0.0_dp
        zw(2,j,i)=0.0_dp
     end do
  end do

END SUBROUTINE G_switch_upcorn2


!> (Based on suitable modifications of S.Goedecker routines)
!! Restore data into output array
!!
!! SYNOPSIS
!!   @param   zf        Original distributed density as well as
!!                      Distributed solution of the poisson equation (inout)
!!   @param  zw         FFT work array
!!   @param  n3         (twice the) dimension of the last FFTtransform.
!!   @param  md1,md3   Dimensions of the undistributed part of the real grid
!!   @param  nfft      number of planes
!!   @param  scal      Needed to achieve unitarity and correct dimensions
!!
!! @warning
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
!!
!! @author
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! Author:S
!!    S. Goedecker, L. Genovese
subroutine P_unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf,scal)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,lot,n3), intent(in) :: zw
  real(kind=8), dimension(md1,md3),intent(inout) :: zf
  real(kind=8), intent(in) :: scal
  !Local variables
  integer :: i3,i1
  real(kind=8) :: pot1

  do i3=1,n3
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        zf(i1,i3)= pot1 
     end do
  end do

END SUBROUTINE P_unfill_downcorn


!> Complex output
subroutine C_unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf,scal)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,lot,n3), intent(in) :: zw
  real(kind=8), dimension(2,md1,md3),intent(inout) :: zf
  real(kind=8), intent(in) :: scal
  !Local variables
  integer :: i3,i1
  real(kind=8) :: pot1,pot2

  !axpy can be used here, if the dimensions were equal
  do i3=1,n3
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        pot2 = scal*zw(2,i1,i3)
        zf(1,i1,i3)= pot1 
        zf(2,i1,i3)= pot2 
     end do
  end do

END SUBROUTINE C_unfill_downcorn


subroutine P_fill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
  implicit none
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(md1,md3), intent(in) :: zf
  real(kind=8), dimension(2,lot,n3), intent(out) :: zw
  !local variables
  integer :: i1,i3

  do i3=1,n3
     do i1=1,nfft
      zw(1,i1,i3)=zf(i1,i3)
      zw(2,i1,i3)=0.d0
     end do
  end do

END SUBROUTINE P_fill_upcorn


!> To be used for complex input
subroutine C_fill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
  implicit none
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,md1,md3), intent(in) :: zf
  real(kind=8), dimension(2,lot,n3), intent(out) :: zw
  !local variables
  integer :: i1,i3

  do i3=1,n3
     do i1=1,nfft
        zw(1,i1,i3)=zf(1,i1,i3)
        zw(2,i1,i3)=zf(2,i1,i3)
     end do
  end do

END SUBROUTINE C_fill_upcorn


!> (Based on suitable modifications of S.Goedecker routines)
!! Assign the correct planes to the work array zmpi2
!! in order to prepare for interprocessor data transposition.
!!
!! SYNOPSIS
!!   @param  zmpi2          Work array for multiprocessor manipulation (output)
!!   @param  zw             Work array (input)
!!   @param  n1,n3          logical dimension of the FFT transform, reference for work arrays
!!   @param  md2,nd3        Dimensions of real grid and of the kernel, respectively
!!   @param  i1,j2,lot,nfft Starting points of the plane and number of remaining lines
!!
!! @author
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! Author:S
!!    S. Goedecker, L. Genovese
subroutine scramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3), intent(in) :: zw
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(inout) :: zmpi2
  !Local variables
  integer :: i3,i

  do i3=1,n3/2+1
     do i=0,nfft-1
        zmpi2(1,i1+i,j2,i3)=zw(1,i+1,i3)
        zmpi2(2,i1+i,j2,i3)=zw(2,i+1,i3)
     end do
  end do

END SUBROUTINE scramble_P


!> (Based on suitable modifications of S.Goedecker routines)
!! Insert the correct planes of the work array zmpi2
!! in order to prepare for backward FFT transform
!!
!! SYNOPSIS
!!   @param  zmpi2:          Work array for multiprocessor manipulation (input)
!!   @param  zw:             Work array (output)
!!   @param  n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!   @param  md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!   @param  i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! @author
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! Author:S
!!    S. Goedecker, L. Genovese
subroutine unscramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3), intent(out) :: zw
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(in) :: zmpi2
  !Local variables
  integer :: i3,i,j3

  i3=1
  do i=0,nfft-1
     zw(1,i+1,i3)=zmpi2(1,i1+i,j2,i3)
     zw(2,i+1,i3)=zmpi2(2,i1+i,j2,i3)
  end do

  do i3=2,n3/2+1
     j3=n3+2-i3
     do i=0,nfft-1
        zw(1,i+1,j3)= zmpi2(1,i1+i,j2,i3)
        zw(2,i+1,j3)=-zmpi2(2,i1+i,j2,i3)
        zw(1,i+1,i3)= zmpi2(1,i1+i,j2,i3)
        zw(2,i+1,i3)= zmpi2(2,i1+i,j2,i3)
     end do
  end do

END SUBROUTINE unscramble_P


!> (Based on suitable modifications of S.Goedecker routines)
!! Multiply with the kernel taking into account its symmetry
!! Conceived to be used into convolution loops
!!
!! SYNOPSIS
!!   @param  pot:      Kernel, symmetric and real, half the length
!!   @param  zw:       Work array (input/output)
!!   @param  n1,n2:    logical dimension of the FFT transform, reference for zw
!!   @param  nd1,nd2:  Dimensions of POT
!!   @param  jS, nfft: starting point of the plane and number of remaining lines
!!   @param  offset  : Offset to be defined for periodic BC (usually 0)
!!   @param  strten  : Components of the Hartree stress tensor order: (11,22,33,23,13,12) !
!!
!! @author
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     Copyright (C) 2009 Luigi Genovese, ESRF Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! Author:S
!!    S. Goedecker, L. Genovese
subroutine P_multkernel(nd1,nd2,n1,n2,n3,lot,nfft,jS,pot,zw,j3,hx,hy,hz,offset,scal,strten)
  use Poisson_Solver, only: dp, gp
  implicit none
  !Argments
  integer, intent(in) :: nd1,nd2,n1,n2,n3,lot,nfft,jS,j3
  real(dp), intent(in) :: hx,hy,hz,offset,scal
  real(kind=8), dimension(nd1,nd2), intent(in) :: pot
  real(kind=8), dimension(2,lot,n2), intent(inout) :: zw
  real(dp), dimension(6), intent(inout) :: strten
  !real(kind=8), dimension(0:n1/2), intent(in) :: fourisf
  !Local variables
  real(gp), parameter :: pi=3.14159265358979323846_gp
  integer :: i1,j1,i2,j2
  real(gp) :: rhog2,g2
  real(dp), dimension(3) :: pxyz

  !acerioni --- triclinic cell ::: can the problem be here?!
  !write (*,*) 'P_multkernel) nfft = ', nfft
  !Body
  !running recip space coordinates
  pxyz(3)=real(j3-1,dp)/(real(n3,dp)*hy)

  !j3=1 case (it might contain the zero component)
  if (j3==1) then
     if (jS==1) then
        !zero fourier component (no strten contribution)
        i2=1
        zw(1,1,1)=offset/(hx*hy*hz) 
        zw(2,1,1)=0.d0              
        !running recip space coordinates
        pxyz(2)=0.0_dp
        do i1=2,nfft
           j1=i1+jS-1
           !running recip space coordinate
           pxyz(1)=real(j1-(j1/(n1/2+2))*n1-1,dp)/(real(n1,dp)*hx) 
           !square of modulus of recip space coordinate
           g2=pxyz(1)**2+pxyz(2)**2+pxyz(3)**2
           !density squared over modulus
           rhog2=(zw(1,i1,i2)**2+zw(2,i1,i2)**2)/g2
           rhog2=rhog2/pi*scal**2
           !stress tensor components (diagonal part)
           strten(1)=strten(1)+(pxyz(1)**2/g2-0.5_dp)*rhog2
           strten(3)=strten(3)+(pxyz(3)**2/g2-0.5_dp)*rhog2
           strten(2)=strten(2)+(pxyz(2)**2/g2-0.5_dp)*rhog2
           strten(5)=strten(5)+(pxyz(1)*pxyz(2)/g2)*rhog2
           strten(6)=strten(6)+(pxyz(1)*pxyz(3)/g2)*rhog2
           strten(4)=strten(4)+(pxyz(2)*pxyz(3)/g2)*rhog2

           j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
           j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)
           zw(1,i1,i2)=zw(1,i1,i2)*pot(j1,j2)
           zw(2,i1,i2)=zw(2,i1,i2)*pot(j1,j2)
        end do
        do i2=2,n2
           !running recip space coordinate
           pxyz(2)=real(i2-(i2/(n2/2+2))*n2-1,dp)/(real(n2,dp)*hz) 

           do i1=1,nfft
              j1=i1+jS-1
              !running recip space coordinate
              pxyz(1)=real(j1-(j1/(n1/2+2))*n1-1,dp)/(real(n1,dp)*hx) 

              !square of modulus of recip space coordinate
              g2=pxyz(1)**2+pxyz(2)**2+pxyz(3)**2
              !density squared over modulus
              rhog2=(zw(1,i1,i2)**2+zw(2,i1,i2)**2)/g2
              rhog2=rhog2/pi*scal**2
              !stress tensor components (diagonal part)
              strten(1)=strten(1)+(pxyz(1)**2/g2-0.5_dp)*rhog2
              strten(3)=strten(3)+(pxyz(2)**2/g2-0.5_dp)*rhog2
              strten(2)=strten(2)+(pxyz(3)**2/g2-0.5_dp)*rhog2
              strten(5)=strten(5)+(pxyz(1)*pxyz(2)/g2)*rhog2
              strten(6)=strten(6)+(pxyz(1)*pxyz(3)/g2)*rhog2
              strten(4)=strten(4)+(pxyz(2)*pxyz(3)/g2)*rhog2

              j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
              j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)
              zw(1,i1,i2)=zw(1,i1,i2)*pot(j1,j2)
              zw(2,i1,i2)=zw(2,i1,i2)*pot(j1,j2)
           end do
        end do
     else
        !generic case
        do i2=1,n2
           !running recip space coordinate
           pxyz(2)=real(i2-(i2/(n2/2+2))*n2-1,dp)/(real(n2,dp)*hz) 

           do i1=1,nfft
              j1=i1+jS-1
              !running recip space coordinate
              pxyz(1)=real(j1-(j1/(n1/2+2))*n1-1,dp)/(real(n1,dp)*hx) 

              !square of modulus of recip space coordinate
              g2=pxyz(1)**2+pxyz(2)**2+pxyz(3)**2
              !density squared over modulus
              rhog2=(zw(1,i1,i2)**2+zw(2,i1,i2)**2)/g2
              rhog2=rhog2/pi*scal**2
              !stress tensor components (diagonal part)
              strten(1)=strten(1)+(pxyz(1)**2/g2-0.5_dp)*rhog2
              strten(3)=strten(3)+(pxyz(2)**2/g2-0.5_dp)*rhog2
              strten(2)=strten(2)+(pxyz(3)**2/g2-0.5_dp)*rhog2
              strten(5)=strten(5)+(pxyz(1)*pxyz(2)/g2)*rhog2
              strten(6)=strten(6)+(pxyz(1)*pxyz(3)/g2)*rhog2
              strten(4)=strten(4)+(pxyz(2)*pxyz(3)/g2)*rhog2

              j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
              j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)
        !acerioni
        !write(14,*) i1,j1,i2,j2
        !acerioni
              zw(1,i1,i2)=zw(1,i1,i2)*pot(j1,j2)
              zw(2,i1,i2)=zw(2,i1,i2)*pot(j1,j2)
           end do
        end do
     end if
           !acerioni
           !write(*,*) 'P_multkernel) nd2 =', nd2
           !write(19,*) i1,i2,zw(1,i1,i2),zw(2,i1,i2),1/pot(j1,j2)
           !write(19,*) i1,i2,1/pot(i1,i2)
  else
     !generic case
     do i2=1,n2
        !running recip space coordinate
        pxyz(2)=real(i2-(i2/(n2/2+2))*n2-1,dp)/(real(n2,dp)*hz) 
        do i1=1,nfft
           j1=i1+jS-1
           !running recip space coordinate
           pxyz(1)=real(j1-(j1/(n1/2+2))*n1-1,dp)/(real(n1,dp)*hx) 

           !square of modulus of recip space coordinate
           g2=pxyz(1)**2+pxyz(2)**2+pxyz(3)**2
           !density squared over modulus
           rhog2=(zw(1,i1,i2)**2+zw(2,i1,i2)**2)/g2
           rhog2=rhog2/pi*scal**2
           !stress tensor components (diagonal part)
           strten(1)=strten(1)+(pxyz(1)**2/g2-0.5_dp)*rhog2
           strten(3)=strten(3)+(pxyz(2)**2/g2-0.5_dp)*rhog2
           strten(2)=strten(2)+(pxyz(3)**2/g2-0.5_dp)*rhog2
           strten(5)=strten(5)+(pxyz(1)*pxyz(2)/g2)*rhog2
           strten(6)=strten(6)+(pxyz(1)*pxyz(3)/g2)*rhog2
           strten(4)=strten(4)+(pxyz(2)*pxyz(3)/g2)*rhog2
           !avoid mirroring for last j3 point
           if (j3 /= n3/2+1) then
              !stress tensor components (diagonal part)
              strten(1)=strten(1)+(pxyz(1)**2/g2-0.5_dp)*rhog2
              strten(3)=strten(3)+(pxyz(2)**2/g2-0.5_dp)*rhog2
              strten(2)=strten(2)+(pxyz(3)**2/g2-0.5_dp)*rhog2
              strten(5)=strten(5)+(pxyz(1)*pxyz(2)/g2)*rhog2
              strten(6)=strten(6)+(pxyz(1)*pxyz(3)/g2)*rhog2
              strten(4)=strten(4)+(pxyz(2)*pxyz(3)/g2)*rhog2
           end if

           j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
           j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)
           zw(1,i1,i2)=zw(1,i1,i2)*pot(j1,j2)
           zw(2,i1,i2)=zw(2,i1,i2)*pot(j1,j2)


           !acerioni /// attempt to use PBCs to mimick monoclinic cell
           !zw(1,i1,i2)=zw(1,i1,i2)*pot(i1,i2)
           !zw(2,i1,i2)=zw(2,i1,i2)*pot(i1,i2)
           !acerioni
        end do
     end do
  end if
END SUBROUTINE P_multkernel


!> (Based on suitable modifications of S.Goedecker routines)
!! Multiply with the kernel taking into account its symmetry
!! Conceived to be used into convolution loops
!!
!! SYNOPSIS
!!   @param  pot:      Kernel, symmetric and real, half the length
!!   @param  zw:       Work array (input/output)
!!   @param  n1,n2:    logical dimension of the FFT transform, reference for zw
!!   @param  nd1,nd2:  Dimensions of POT
!!   @param  jS, nfft: starting point of the plane and number of remaining lines
!!
!! @author
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! Author:S
!!    S. Goedecker, L. Genovese
subroutine multkernel(nd1,nd2,n1,n2,lot,nfft,jS,pot,zw)
  implicit none
  !Argments
  integer, intent(in) :: nd1,nd2,n1,n2,lot,nfft,jS
  real(kind=8), dimension(nd1,nd2), intent(in) :: pot
  real(kind=8), dimension(2,lot,n2), intent(inout) :: zw
  !Local variables
  integer :: j,j1,i2,j2

  !Body
  
  !case i2=1
  do j=1,nfft
     j1=j+jS-1
     !isign=(j1/(n1/2+2))
     !j1=(1-2*isign)*j1+isign*(n1+2) !n1/2+1-abs(n1/2+2-jS-i1)
     j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
     !!     j1=n1/2+1-abs(n1/2+2-jS-j)!this stands for j1=min(jS-1+j,n1+3-jS-j)
     zw(1,j,1)=zw(1,j,1)*pot(j1,1)
     zw(2,j,1)=zw(2,j,1)*pot(j1,1)
  end do


  !generic case
  do i2=2,n2/2
     do j=1,nfft
        j1=j+jS-1
        j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!        j1=n1/2+1-abs(n1/2+2-jS-j)
        j2=n2+2-i2
        zw(1,j,i2)=zw(1,j,i2)*pot(j1,i2)
        zw(2,j,i2)=zw(2,j,i2)*pot(j1,i2)
        zw(1,j,j2)=zw(1,j,j2)*pot(j1,i2)
        zw(2,j,j2)=zw(2,j,j2)*pot(j1,i2)
      !  write(*,*) 'Re[zw(i2=', i2, ')] = ', zw(1,j,i2)
      !  write(*,*) 'Im[zw(i2=', i2, ')] = ', zw(2,j,i2)
     end do
  end do
  
  !case i2=n2/2+1
  do j=1,nfft
     j1=j+jS-1
     j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!     j1=n1/2+1-abs(n1/2+2-jS-j)
     j2=n2/2+1
     zw(1,j,j2)=zw(1,j,j2)*pot(j1,j2)
     zw(2,j,j2)=zw(2,j,j2)*pot(j1,j2)
  end do

END SUBROUTINE multkernel


subroutine G_unswitch_downcorn(nfft,n2,n2dim,lot,n1,lzt,zw,zt)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt,n2dim
  real(kind=8), dimension(2,lot,n2), intent(in) :: zw
  real(kind=8), dimension(2,lzt,n1), intent(inout) :: zt
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2dim
      zt(1,i,j)=zw(1,j,i)
      zt(2,i,j)=zw(2,j,i)
     end do
  end do

END SUBROUTINE G_unswitch_downcorn


subroutine G_unswitch_downcorn2(nfft,n2,n2dim,lot,n1,lzt,zw,zt,n3pr1,zt_init)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt,n2dim,n3pr1,zt_init
  real(kind=8), dimension(2,lot,n2), intent(in) :: zw
  real(kind=8), dimension(2,lzt/n3pr1,n1/n3pr1,n3pr1), intent(inout) :: zt
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2dim
      zt(1,mod(i-1,lzt/n3pr1)+1,j+zt_init-1,(i-1)/(lzt/n3pr1)+1)=zw(1,j,i)
      zt(2,mod(i-1,lzt/n3pr1)+1,j+zt_init-1,(i-1)/(lzt/n3pr1)+1)=zw(2,j,i)
     end do
  end do

END SUBROUTINE G_unswitch_downcorn2


subroutine G_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,&
    n1dim,md2,nd3,nproc,zw,zmpi1)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc,n1dim
  integer, intent(inout) :: Jp2stf,J2stf
  real(kind=8), dimension(2,lot,n1), intent(in) :: zw
  real(kind=8), dimension(2,n1dim,md2/nproc,nd3/nproc,nproc), intent(inout) :: zmpi1
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stf,nproc
     do J2=J2stf,md2/nproc
        mfft=mfft+1
        if (mfft > nfft) then
           Jp2stf=Jp2
           J2stf=J2
           return
        end if
        do I1=1,n1dim
           zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
           zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
     end do
     J2stf=1
  end do
END SUBROUTINE G_unmpiswitch_downcorn


subroutine G_unmpiswitch_downcorn2(j3,nfft,Jp2stf,J2stf,lot,n1,&
    n1dim,md2,nd3,n3pr1,n3pr2,zw,zmpi1)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,n3pr1,n3pr2,n1dim
  integer, intent(inout) :: Jp2stf,J2stf
  real(kind=8), dimension(2,lot,n1), intent(in) :: zw
  real(kind=8), dimension(2,n1dim,md2/(n3pr1*n3pr2),nd3/n3pr2,n3pr2), intent(inout) :: zmpi1
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stf,n3pr2
     do J2=J2stf,md2/(n3pr1*n3pr2)
        mfft=mfft+1
        if (mfft > nfft) then
           Jp2stf=Jp2
           J2stf=J2
           return
        end if
        do I1=1,n1dim
           zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
           zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
     end do
     J2stf=1
  end do
END SUBROUTINE G_unmpiswitch_downcorn2


!>  (Based on suitable modifications of S.Goedecker routines)
!!  Restore data into output array, calculating in the meanwhile
!!  Hartree energy of the potential 
!!
!! SYNOPSIS
!!   @param  zf:          Original distributed density as well as
!!   @param               Distributed solution of the poisson equation (inout)
!!   @param  zw:          FFT work array
!!   @param  n3:          (twice the) dimension of the last FFTtransform.
!!   @param  md1,md3:     Dimensions of the undistributed part of the real grid
!!   @param  nfft:        number of planes
!!   @param  scal:        Needed to achieve unitarity and correct dimensions
!!   @param  ehartreetmp: Hartree energy
!!
!! @warning
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
!!
!! @author
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! Author:S
!!    S. Goedecker, L. Genovese
subroutine unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf&
     ,scal)!,ehartreetmp)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(md1,md3),intent(inout) :: zf
  real(kind=8), intent(in) :: scal
  !real(kind=8), intent(out) :: ehartreetmp
  !Local variables
  integer :: i3,i1
  real(kind=8) :: pot1

  !Execution
  !ehartreetmp=0.d0
  do i3=1,n3/4
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        !ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3-1)
      zf(i1,2*i3-1)= pot1 
      pot1 = scal*zw(2,i1,i3)
        !ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3)
      zf(i1,2*i3)= pot1 
     enddo
  end do
  
END SUBROUTINE unfill_downcorn


subroutine halfill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
  implicit real(kind=8) (a-h,o-z)
!Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8) ::  zw(2,lot,n3/2),zf(md1,md3)
!Local variables
  integer :: i1,i3
! WARNING: Assuming that high frequencies are in the corners 
!          and that n3 is multiple of 4
!in principle we can relax this condition
      
  do i3=1,n3/4
     do i1=1,nfft
        zw(1,i1,i3)=0.d0
        zw(2,i1,i3)=0.d0
     end do
  end do
  do i3=n3/4+1,n3/2
     do i1=1,nfft
        zw(1,i1,i3)=zf(i1,2*i3-1-n3/2)
        zw(2,i1,i3)=zf(i1,2*i3-n3/2)
     end do
  end do
      
END SUBROUTINE halfill_upcorn


!> (Based on suitable modifications of S.Goedecker routines)
!! Assign the correct planes to the work array zmpi2
!! in order to prepare for interprocessor data transposition.
!! In the meanwhile, it unpacks the data of the HalFFT in order to prepare for
!! multiplication with the kernel
!!
!! SYNOPSIS
!!   @param  zmpi2:          Work array for multiprocessor manipulation (output)
!!   @param  zw:             Work array (input)
!!   @param  cosinarr:      Array of the phases needed for unpacking
!!   @param  n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!   @param  md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!   @param  i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! @author
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! Author:S
!!    S. Goedecker, L. Genovese
subroutine scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(2,n3/2), intent(in) :: cosinarr
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(inout) :: zmpi2
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

END SUBROUTINE scramble_unpack

 
!> (Based on suitable modifications of S.Goedecker routines)
!! Insert the correct planes of the work array zmpi2
!! in order to prepare for backward FFT transform
!! In the meanwhile, it packs the data in order to be transformed with the HalFFT 
!! procedure
!!
!! SYNOPSIS
!!   @param  zmpi2:          Work array for multiprocessor manipulation (input)
!!   @param  zw:             Work array (output)
!!   @param  cosinarr:       Array of the phases needed for packing
!!   @param  n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!   @param  md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!   @param  i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
!!
!! @author
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! Author:S
!!    S. Goedecker, L. Genovese
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

END SUBROUTINE unscramble_pack
