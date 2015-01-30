!> @file
!!  Performs a check of the Poisson Solver suite by running with different regimes
!!  and for different choices of the XC functionals
!! @author
!!    Copyright (C) 2002-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program to test the Poisson Solver
program PS_Check
   use wrapper_mpi
   use Poisson_Solver
   use yaml_output
   use dynamic_memory
   use dictionaries
   use time_profiling
   implicit none
   !Length of the box
   character(len=*), parameter :: subname='PS_Check'
   real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
   real(kind=8), parameter :: acell = 10.d0
   character(len=50) :: chain
   character(len=1) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
   character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
   real(kind=8), dimension(:), allocatable :: density,rhopot,potential,pot_ion,extra_ref
   type(coulomb_operator) :: pkernel,pkernelseq
   real(kind=8) :: hx,hy,hz,offset
   real(kind=8) :: ehartree
   real(kind=8) :: tel
   real :: tcpu0,tcpu1
   integer :: ncount0,ncount1,ncount_rate,ncount_max
   integer :: n01,n02,n03,itype_scf!,i_all,i_stat
   integer :: iproc,nproc,namelen,ierr,ispden
   integer :: n_cell
   integer, dimension(3) :: nxyz
   integer, dimension(3) :: ndims
   real(wp), dimension(:,:,:,:), pointer :: rhocore
   real(dp), dimension(3) :: hgrids
   type(mpi_environment) :: bigdft_mpi
   character(len = *), parameter :: package_version = "PSolver 1.7-dev.25"
   external :: gather_timings

   call f_lib_initialize() 

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

   !initialize categories for the Poisson Solver
   call PS_initialize_timing_categories()
  

   call f_malloc_set_status(memory_limit=0.e0,iproc=iproc)
   call f_routine(id='PS_Check')

   bigdft_mpi%mpi_comm=MPI_COMM_WORLD !workaround to be removed

   if (iproc ==0) then
      call yaml_set_stream(record_length=92,tabbing=30)!unit=70,filename='log.yaml')
      call yaml_new_document()

      call yaml_map('Reference Paper','The Journal of Chemical Physics 137, 134108 (2012)')
      call yaml_map('Version Number', "PSolver " // trim(PS_getVersion()))
      call yaml_map('Timestamp of this run',yaml_date_and_time_toa())
      call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
      if (ierr ==0) call yaml_map('Root process Hostname',trim(nodename_local))
   end if

   !initialize memory counting and timings
   !call memocc(0,iproc,'count','start')
   !call timing(nproc,'time.yaml','IN')
   call f_timing_reset(filename='time.yaml',master=iproc==0)


   !Start global timing
   call cpu_time(tcpu0)
   call system_clock(ncount0,ncount_rate,ncount_max)

   !the first proc read the data and then send them to the others
   if (iproc==0) then
      !Use arguments
      call get_command_argument(1,value=chain)
      read(unit=chain,fmt=*) nxyz(1)
      call get_command_argument(2,value=chain)
      read(unit=chain,fmt=*) nxyz(2)
      call get_command_argument(3,value=chain)
      read(unit=chain,fmt=*) nxyz(3)
      call get_command_argument(4,value=chain)
      read(unit=chain,fmt=*) geocode
   end if

   call MPI_BCAST(nxyz,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(geocode,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

   n01=nxyz(1)
   n02=nxyz(2)
   n03=nxyz(3)

   !print *,iproc,n01,n02,n03

   !Step size
   n_cell = max(n01,n02,n03)
   hx=acell/real(n01,kind=8)
   hy=acell/real(n02,kind=8)
   hz=acell/real(n03,kind=8)

   !order of the scaling functions chosen
   itype_scf=16

   if (iproc==0) then
      select case(geocode)
      case('F')
         call yaml_map('Boundary Conditions','Isolated')
      case('S')
         call yaml_map('Boundary Conditions','Surface')
      case('W')
         call yaml_map('Boundary Conditions','Wire')
      case('P')
         call yaml_map('Boundary Conditions','Periodic')
      end select
      call yaml_mapping_open('Multiprocessor run',label='MPIrun')
   end if

   !calculate the kernel in parallel for each processor
   ndims=(/n01,n02,n03/)
   hgrids=(/hx,hy,hz/)

   pkernel=pkernel_init(.true.,iproc,nproc,0,&
        geocode,ndims,hgrids,itype_scf,taskgroup_size=nproc/2)
   call pkernel_set(pkernel,.true.)

   !Allocations, considering also spin density
   !Density
   density=f_malloc(n01*n02*n03*2,id='density')
   !Density then potential
   potential=f_malloc(n01*n02*n03,id='potential')
   !ionic potential
   pot_ion=f_malloc(n01*n02*n03,id='pot_ion')

   extra_ref=f_malloc(n01*n02*n03,id='extra_ref')

   nullify(rhocore)

   ! No XC, test only ispden = 1.
   ispden=1

      if (iproc == 0) then
         call yaml_map('Number of Spins',ispden,advance='no')
         call yaml_comment('nspden:'//trim(yaml_toa(ispden)),hfill='-')
      end if

      !if (iproc == 0) call yaml_comment('nspden:'//yaml_toa(ispden,fmt='(i0)'),hfill='=')
      !write(unit=*,fmt="(1x,a,i0)")  '===================== nspden:  ',ispden
      !then assign the value of the analytic density and the potential
      !allocate the rhopot also for complex routines
      rhopot=f_malloc(n01*n02*n03*2,id='rhopot')
!!$      allocate(rhopot(n01*n02*n03*2+ndebug),stat=i_stat)
!!$      call memocc(i_stat,rhopot,'rhopot',subname)

      call test_functions(geocode,n01,n02,n03,ispden,acell,a_gauss,hx,hy,hz,&
      density,potential,rhopot,pot_ion,offset)
      !calculate the Poisson potential in parallel
      !with the global data distribution (also for xc potential)
      call H_potential('G',pkernel,rhopot,pot_ion,ehartree,offset,.false.) !optional argument

      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0) then
         call yaml_mapping_open('Energies',flow=.true.)
         call yaml_map('Hartree',ehartree,fmt='(1pe20.12)')
         call yaml_mapping_close()
         call yaml_mapping_open('Comparison with a reference run')
      end if
      !write(unit=*,fmt="(1x,a,3(1pe20.12))") 'Energies:',ehartree,eexcu,vexcu
      !stop
      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0) then
         !compare the values of the analytic results (pkernel%mpi_env%nproc == -1 indicates that it is serial)
         call compare (0,-1,pkernel%mpi_env%mpi_comm,n01,n02,n03,1,potential,rhopot,'ANALYTIC')
      end if
      !if the latter test pass, we have a reference for all the other calculations
      !build the reference quantities (based on the numerical result, not the analytic)
      potential(:)=rhopot(1:n01*n02*n03)
      extra_ref=potential

      !now the parallel calculation part
      call f_free(rhopot)
!!$      i_all=-product(shape(rhopot))*kind(rhopot)
!!$      deallocate(rhopot,stat=i_stat)
!!$      call memocc(i_stat,i_all,'rhopot',subname)

      call compare_with_reference(pkernel%mpi_env%nproc,geocode,'G',n01,n02,n03,ispden,offset,ehartree,&
      density,potential,pot_ion,pkernel)

      call compare_with_reference(pkernel%mpi_env%nproc,geocode,'D',n01,n02,n03,ispden,offset,ehartree,&
      density,potential,pot_ion,pkernel)

      !test for the serial solver (always done to have a simpler comparison)
      !if (pkernel%iproc == 0 .and. pkernel%mpi_env%nproc > 1 ) then
!      if (pkernel%iproc == 0) then
!         i_all=-product(shape(pkernel))*kind(pkernel)
!         deallocate(pkernel,stat=i_stat)
!         call memocc(i_stat,i_all,'pkernel',subname)
!
!         !calculate the kernel 
!         call createKernel(0,1,geocode,n01,n02,n03,hx,hy,hz,itype_scf,pkernel,.false.)
!
!         call compare_with_reference(1,geocode,'G',n01,n02,n03,ixc,ispden,offset,ehartree,eexcu,vexcu,&
!         density,potential,pot_ion,xc_pot,pkernel)
!
!         call compare_with_reference(1,geocode,'D',n01,n02,n03,ixc,ispden,offset,ehartree,eexcu,vexcu,&
!         density,potential,pot_ion,xc_pot,pkernel)
!      end if

      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0) call yaml_mapping_close() !comparison


   if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0) call yaml_mapping_close() !MPI

      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0) call yaml_mapping_open('Complex run')
      !compare the calculations in complex
      call compare_cplx_calculations(pkernel%mpi_env%iproc,pkernel%mpi_env%nproc,geocode,'G',n01,n02,n03,ehartree,offset,&
      density,potential,pkernel)

      call compare_cplx_calculations(pkernel%mpi_env%iproc,pkernel%mpi_env%nproc,geocode,'D',n01,n02,n03,ehartree,offset,&
      density,potential,pkernel)
      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0)call yaml_mapping_close()

      call f_timing_checkpoint('Parallel',mpi_comm=MPI_COMM_WORLD,nproc=nproc,gather_routine=gather_timings)
      !call timing(MPI_COMM_WORLD,'Parallel','PR')

   call pkernel_free(pkernel)

   if (pkernel%mpi_env%nproc == 1 .and.pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0 )&
        call yaml_map('Monoprocess run','*MPIrun')

   !do not do the sequential calculation if it has been already done
   if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0 .and. pkernel%mpi_env%nproc > 1 ) then
      call yaml_mapping_open('Monoprocess run')
      rhopot=f_malloc(n01*n02*n03*2,id='rhopot')
!!$      allocate(rhopot(n01*n02*n03*2+ndebug),stat=i_stat)
!!$      call memocc(i_stat,rhopot,'rhopot',subname)

      ispden = 1
       call yaml_map('Number of Spins',ispden)

       call test_functions(geocode,n01,n02,n03,ispden,acell,a_gauss,hx,hy,hz,&
            density,potential,rhopot,pot_ion,offset)
       potential=extra_ref !use the previoulsy defined reference
      !calculate the Poisson potential in parallel
      !with the global data distribution (also for xc potential)
       pkernelseq=pkernel_init(.true.,0,1,0,geocode,ndims,hgrids,itype_scf)
       call pkernel_set(pkernelseq,.true.)

!!$       call createKernel(0,1,geocode,(/n01,n02,n03/),(/hx,hy,hz/),itype_scf,pkernelseq,.true.)
       call yaml_mapping_open('Comparison with a reference run')

       call compare_with_reference(1,geocode,'G',n01,n02,n03,ispden,offset,ehartree,&
         density,potential,pot_ion,pkernelseq)

       call compare_with_reference(1,geocode,'D',n01,n02,n03,ispden,offset,ehartree,&
         density,potential,pot_ion,pkernelseq)

       call pkernel_free(pkernelseq)
       call yaml_mapping_close() !comparison

     call f_free(rhopot)
!!$     i_all=-product(shape(rhopot))*kind(rhopot)
!!$     deallocate(rhopot,stat=i_stat)
!!$     call memocc(i_stat,i_all,'rhopot',subname)

     call yaml_mapping_close()
   endif

   call f_timing_checkpoint('Serial',mpi_comm=MPI_COMM_WORLD,&
        nproc=nproc,gather_routine=gather_timings)
   !call timing(MPI_COMM_WORLD,'Serial','PR')

   !call f_malloc_dump_status()

   call f_free(density,potential,pot_ion,extra_ref)

   call f_timing_stop(mpi_comm=MPI_COMM_WORLD,nproc=nproc,gather_routine=gather_timings)

   !Final timing
   call cpu_time(tcpu1)
   call system_clock(ncount1,ncount_rate,ncount_max)
   tel=real(ncount1-ncount0,kind=gp)/real(ncount_rate,kind=gp)
   if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0) then
      call yaml_mapping_open('Timings for root process')
        call yaml_map('CPU time (s)',tcpu1-tcpu0,fmt='(f12.2)')
        call yaml_map('Elapsed time (s)',tel,fmt='(f12.2)')
      call yaml_mapping_close()
   end if
   !call yaml_stream_attributes()
   !&   write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time/ELAPSED time for root process ', pkernel%iproc,tel,tcpu1-tcpu0
   
   call f_release_routine()
   if (iproc==0) then
      call yaml_release_document()
      call yaml_close_all_streams()
   end if

   call MPI_FINALIZE(ierr)

   call f_lib_finalize()

   !END PROGRAM PS_Check

   contains

   subroutine compare_cplx_calculations(iproc,nproc,geocode,distcode,n01,n02,n03,ehref,offset,&
      density,potential,pkernel)
      use Poisson_Solver
      use dynamic_memory
      implicit none
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
      character(len=1), intent(in) :: distcode
      integer, intent(in) :: iproc,nproc,n01,n02,n03
      real(kind=8), intent(in) :: ehref,offset
      real(kind=8), dimension(n01*n02*n03), intent(in) :: potential
      real(kind=8), dimension(n01*n02*n03*2), intent(in) :: density
      type(coulomb_operator), intent(in) :: pkernel
      !local varaibles
      character(len=*), parameter :: subname='compare_cplx_calculations'
      character(len=20) :: message
      integer :: n3d,n3p,n3pi,i3xcsh,i3s,i3sd,i3,i2,i1,istden,istpot,isp,i
      real(kind=8) :: ehartree
      real(kind=8), dimension(:,:,:,:), allocatable :: rhopot
      
      call f_routine(id=subname)

!      offset=0.d0

      !this is performed always without XC since a complex
      !charge density makes no sense
      write(message,'(1x,a,1x,a)') geocode,distcode

      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup==0) then
         select case(distcode)
         case('D')
            call yaml_mapping_open('Distributed data')
         case('G')
            call yaml_mapping_open('Global data')
         end select
      end if


      call PS_dim4allocation(geocode,distcode,iproc,nproc,n01,n02,n03,.false.,.false.,&
      n3d,n3p,n3pi,i3xcsh,i3s)

      !starting point of the three-dimensional arrays
      if (distcode == 'D') then
         istden=n01*n02*(i3s-1)+1
         istpot=n01*n02*(i3s+i3xcsh-1)+1
         i3sd=i3s
      else if (distcode == 'G') then
         istden=1
         istpot=1
         i3sd=1
      end if

      !input poisson solver, complex distribution
      rhopot=f_malloc((/n01,n02,n3d,2/),id='rhopot')

      !allocate(rhopot(n01,n02,n3d,2+ndebug),stat=i_stat)
      !call memocc(i_stat,rhopot,'rhopot',subname)

      !do isp=1,2
      isp=1
      do i3=1,n3d
         do i2=1,n02
            do i1=1,n01
               i=i1+(i2-1)*n01+(modulo(i3sd+i3-2,n03))*n01*n02!+(isp-1)*n01*n02*n03
               rhopot(i1,i2,i3,isp)=density(i)
            end do
         end do
      end do
      !end do

      !perform the calculation in complex, with distributed and gathered distribution
      call H_potential(distcode,pkernel,rhopot,rhopot,ehartree,offset,.false.,quiet='YES')

      call compare(pkernel%mpi_env%iproc +pkernel%mpi_env%igroup,-1,pkernel%mpi_env%mpi_comm,&
           n01,n02,n3d,1,potential(istpot),rhopot,'CPLXREAL')

      isp=2
      do i3=1,n3d
         do i2=1,n02
            do i1=1,n01
               i=i1+(i2-1)*n01+(modulo(i3sd+i3-2,n03))*n01*n02!+(isp-1)*n01*n02*n03
               rhopot(i1,i2,i3,isp)=density(i)
            end do
         end do
      end do

      !perform the calculation in complex, with distributed and gathered distribution
      call H_potential(distcode,pkernel,rhopot(1,1,1,2),rhopot,ehartree,offset,.false.,quiet='YES')

      call compare(pkernel%mpi_env%iproc +pkernel%mpi_env%igroup,-1,pkernel%mpi_env%mpi_comm,&
           n01,n02,n3d,1,potential(istpot),rhopot(1,1,1,2),'CPLXIMAG')


      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup==0) then
         call yaml_mapping_open('Energy differences')
         call yaml_map('Hartree',ehref-ehartree,fmt="(1pe20.12)")
         call yaml_mapping_close()
         call yaml_mapping_close() !run
      end if

      ! write(unit=*,fmt="(1x,a,1pe20.12)")'Energy diff:',ehref-ehartree

      call f_free(rhopot)
!!$      i_all=-product(shape(rhopot))*kind(rhopot)
!!$      deallocate(rhopot,stat=i_stat)
!!$      call memocc(i_stat,i_all,'rhopot',subname)
      call f_release_routine()
    
   END SUBROUTINE compare_cplx_calculations


   !> Compare with the given reference
   subroutine compare_with_reference(nproc,geocode,distcode,n01,n02,n03,&
      nspden,offset,ehref,&
      density,potential,pot_ion,pkernel)
      use Poisson_Solver
      use dynamic_memory
      implicit none
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
      character(len=1), intent(in) :: distcode
      integer, intent(in) :: nproc,n01,n02,n03,nspden
      real(kind=8), intent(in) :: offset,ehref
      real(kind=8), dimension(n01*n02*n03), intent(in) :: potential
      real(kind=8), dimension(n01*n02*n03*nspden), intent(in) :: density
      real(kind=8), dimension(n01*n02*n03), intent(inout) :: pot_ion
      type(coulomb_operator), intent(in) :: pkernel
      !local variables
      character(len=*), parameter :: subname='compare_with_reference'
      character(len=100) :: message
      integer :: n3d,n3p,n3pi,i3xcsh,i3s,istden,istpot,istpoti,i!,i_stat,i_all
      integer :: i1,i2,i3,isp,i3sd
      real(kind=8) :: ehartree!,vexcu,eexcu
      real(kind=8), dimension(:), allocatable :: test
      real(kind=8), dimension(:,:,:,:), allocatable :: rhopot
      real(dp), dimension(:,:,:,:), pointer :: rhocore

      call f_routine(id=subname)

      nullify(rhocore)

      call PS_dim4allocation(geocode,distcode,pkernel%mpi_env%iproc,pkernel%mpi_env%nproc,n01,n02,n03,.false.,.false.,&
      n3d,n3p,n3pi,i3xcsh,i3s)

      !starting point of the three-dimensional arrays
      if (distcode == 'D') then
         istden=n01*n02*(i3s-1)+1
         istpot=n01*n02*(i3s+i3xcsh-1)+1
         i3sd=i3s
      else if (distcode == 'G') then
         istden=1
         istpot=1
         i3sd=1
      end if
      istpoti=n01*n02*(i3s+i3xcsh-1)+1

      !test arrays for comparison
      test=f_malloc(n01*n02*n03*nspden,id='test')
      !input poisson solver
      rhopot=f_malloc((/n01,n02,n3d,nspden/),id='rhopot')

!!$      !test arrays for comparison
!!$      allocate(test(n01*n02*n03*nspden+ndebug),stat=i_stat)
!!$      call memocc(i_stat,test,'test',subname)
!!$      !XC potential
!!$      allocate(test_xc(n01*n02*n03*nspden+ndebug),stat=i_stat)
!!$      call memocc(i_stat,test_xc,'test_xc',subname)
!!$      !input poisson solver
!!$      allocate(rhopot(n01,n02,n3d,nspden+ndebug),stat=i_stat)
!!$      call memocc(i_stat,rhopot,'rhopot',subname)
    
      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup == 0) &
           write(message,'(1x,a,1x,i0,1x,a,1x,i0)') geocode,0,distcode,nspden

      if (pkernel%mpi_env%iproc +pkernel%mpi_env%igroup==0) then
         select case(distcode)
         case('D')
            call yaml_mapping_open('Distributed data')
         case('G')
            call yaml_mapping_open('Global data')
         end select
      end if


      test(1:n01*n02*n03)=potential(1:n01*n02*n03)!+pot_ion

      do isp=1,nspden
         !add the initialisation of the density for the periodic GGA case
         do i3=1,n3d
            do i2=1,n02
               do i1=1,n01
                  i=i1+(i2-1)*n01+(modulo(i3sd+i3-2,n03))*n01*n02+(isp-1)*n01*n02*n03
                  rhopot(i1,i2,i3,isp)=density(i)
               end do
            end do
         end do
      end do

      call H_potential(distcode,pkernel,rhopot,rhopot,ehartree,offset,.false.,quiet='yes') !optional argument
      !compare the values of the analytic results (no dependence on spin)
      call compare(pkernel%mpi_env%iproc +pkernel%mpi_env%igroup,nproc,pkernel%mpi_env%mpi_comm,&
           n01,n02,n3p,1,potential(istpot),rhopot(1,1,1,1),&
      'ANACOMPLET ')!//message)

      !!$    call PSolver(geocode,distcode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
      !!$         rhopot(1,1,1,1),pkernel,test_xc,ehartree,eexcu,vexcu,offset,.false.,nspden,quiet='yes')
      !!$    
      !!$    !compare the values of the analytic results (no dependence on spin)
      !!$    call compare(iproc,nproc,n01,n02,n3p,1,potential(istpot),rhopot(1,1,i3xcsh+1,1),&
      !!$         'ANACOMPLET '//message)

      if (pkernel%mpi_env%iproc + pkernel%mpi_env%igroup==0) then
            call yaml_mapping_open('Energy differences')
            call yaml_map('Hartree',ehref-ehartree,fmt="(1pe20.12)")
            call yaml_mapping_close()
         end if
!!$         write(unit=*,fmt="(1x,a,3(1pe20.12))") &
!!$      'Energies diff:',ehref-ehartree,excref-eexcu,vxcref-vexcu

      do isp=1,nspden
         do i3=1,n3d
            do i2=1,n02
               do i1=1,n01
                  i=i1+(i2-1)*n01+(modulo(i3sd+i3-2,n03))*n01*n02+(isp-1)*n01*n02*n03
                  rhopot(i1,i2,i3,isp)=density(i)
               end do
            end do
         end do
      end do

      call H_potential(distcode,pkernel,rhopot(1,1,1,1),pot_ion(istpoti),ehartree,offset,.false.,quiet='yes') !optional argument
      !fill the other part, for spin, polarised
      if (nspden == 2) then
         !the starting point is not so simple
         if (n3d > n3p) then
            call dcopy(n01*n02*n3p,rhopot(1,1,1,1),1,rhopot(1,1,n3p+1,1),1)
         else
            call dcopy(n01*n02*n3p,rhopot(1,1,1,1),1,rhopot(1,1,1,nspden),1)
         end if
      end if
!      print *,'n3d,n3p,isp,i3s,i3sd',n3d,n3p,isp,i3s,i3sd
      !spin up and down together with the XC part
      !call axpy(n01*n02*n3p*nspden,1.0_dp,test_xc(1),1,rhopot(1,1,1,1),1)
      !then compare again, but the complete result
      call compare(pkernel%mpi_env%iproc + pkernel%mpi_env%igroup,-1,pkernel%mpi_env%mpi_comm,n01,n02,nspden*n3p,1,test(istpot),&
      rhopot(1,1,1,1),'COMPLETE   ')!)//message)

      !!$    call PSolver(geocode,distcode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
      !!$         rhopot(1,1,1,1),pkernel,pot_ion(istpoti),ehartree,eexcu,vexcu,offset,.true.,nspden,quiet='yes')
      !!$    !then compare again, but the complete result
      !!$    call compare(iproc,nproc,n01,n02,nspden*n3p,1,test(istpot),&
      !!$         rhopot(1,1,i3xcsh+1,1),'COMPLETE   '//message)
      if (pkernel%mpi_env%iproc + pkernel%mpi_env%igroup==0) then
         call yaml_mapping_open('Energy differences')
         call yaml_map('Hartree',ehref-ehartree,fmt="(1pe20.12)")
         call yaml_mapping_close()
         call yaml_mapping_close() !Run
      end if
!!$      if (iproc==0) write(unit=*,fmt="(1x,a,3(1pe20.12))") &
!!$      'Energies diff:',ehref-ehartree,excref-eexcu,vxcref-vexcu
      

      call f_free(test)
      call f_free(rhopot)

!!$      i_all=-product(shape(test))*kind(test)
!!$      deallocate(test,stat=i_stat)
!!$      call memocc(i_stat,i_all,'test',subname)
!!$      i_all=-product(shape(test_xc))*kind(test_xc)
!!$      deallocate(test_xc,stat=i_stat)
!!$      call memocc(i_stat,i_all,'test_xc',subname)
!!$
!!$      i_all=-product(shape(rhopot))*kind(rhopot)
!!$      deallocate(rhopot,stat=i_stat)
!!$      call memocc(i_stat,i_all,'rhopot',subname)

      call f_release_routine()

   END SUBROUTINE compare_with_reference


   !> Compare arrays potential and density
   !! if nproc == -1: serial version i.e. special comparison
   subroutine compare(iproc,nproc,mpi_comm,n01,n02,n03,nspden,potential,density,description)
      use wrapper_mpi
      implicit none
      !include 'mpif.h'
      character(len=*), intent(in) :: description
      integer, intent(in) :: iproc,nproc,n01,n02,n03,nspden,mpi_comm
      real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,density
      !local variables
      integer :: i1,i2,i3,ierr,i1_max,i2_max,i3_max
      real(kind=8) :: factor,diff_par,max_diff
      max_diff = 0.d0
      i1_max = 1
      i2_max = 1
      i3_max = 1
      do i3=1,n03
         do i2=1,n02 
            do i1=1,n01
               factor=abs(real(nspden,kind=8)*potential(i1,i2,i3)-density(i1,i2,i3))
               if (max_diff < factor) then
                  max_diff = factor
                  i1_max = i1
                  i2_max = i2
                  i3_max = i3
               end if
            end do
         end do
      end do

      !!!  print *,'iproc,i3xcsh,i3s,max_diff',iproc,i3xcsh,i3s,max_diff

      if (nproc > 1) then
         !extract the max
         call MPI_ALLREDUCE(max_diff,diff_par,1,MPI_double_precision,  &
         MPI_MAX,mpi_comm,ierr)
      else
         diff_par=max_diff
      end if

      if (iproc == 0) then
         call yaml_mapping_open(trim(description),flow=.false.)
         !call yaml_mapping_open('Result comparison')
         !call yaml_map('Run description',trim(description))
         call yaml_map('Difference in Inf. Norm',diff_par,fmt='(1pe20.12)')
         if (diff_par > 1.e-10) call yaml_warning('Calculation possibly wrong, check if the diff is meaningful')
         if (nproc == -1) then
            call yaml_map('Max. diff coordinates',(/i1_max,i2_max,i3_max/),fmt='(i0)')
            call yaml_map('Result',density(i1_max,i2_max,i3_max),fmt='(1pe20.12)')
            call yaml_map('Original',potential(i1_max,i2_max,i3_max),fmt='(1pe20.12e3)')
         end if
         call yaml_mapping_close()

!!$         if (nproc == -1) then
!!$            if (diff_par > 1.e-10) then
!!$               write(unit=*,fmt="(1x,a,1pe20.12,a)") &
!!$               trim(description) // '    Max diff:',diff_par,'   <<<< WARNING'
!!$               write(unit=*,fmt="(1x,a,1pe20.12)") &
!!$               '      result:',density(i1_max,i2_max,i3_max),&
!!$               '    original:',potential(i1_max,i2_max,i3_max)
!!$               write(*,'(a,3(i0,1x))') '  Max diff at: ',i1_max,i2_max,i3_max
!!$               !!!           i3=i3_max
!!$               !!!           i1=i1_max
!!$               !!!           do i2=1,n02
!!$               !!!              !do i1=1,n01
!!$               !!!                 write(20,*)i1,i2,potential(i1,i2,i3),density(i1,i2,i3)
!!$               !!!              !end do
!!$               !!!           end do
!!$               !!!           stop
!!$            else
!!$               write(unit=*,fmt="(1x,a,1pe20.12)") &
!!$               trim(description) // '    Max diff:',diff_par
!!$            end if
!!$         else
!!$            if (diff_par > 1.e-10) then
!!$               write(unit=*,fmt="(1x,a,1pe20.12,a)") &
!!$               trim(description) // '    Max diff:',diff_par,'   <<<< WARNING'
!!$            else
!!$               write(unit=*,fmt="(1x,a,1pe20.12)") &
!!$               trim(description) //'    Max diff:',diff_par
!!$            end if
!!$         end if
      end if

      max_diff=diff_par

   END SUBROUTINE compare


   !> This subroutine builds some analytic functions that can be used for 
   !! testing the poisson solver.
   !! The default choice is already well-tuned for comparison.
   !! WARNING: not all the test functions can be used for all the boundary conditions of
   !! the poisson solver, in order to have a reliable analytic comparison.
   !! The parameters of the functions must be adjusted in order to have a sufficiently localized
   !! function in the isolated direction and an explicitly periodic function in the periodic ones.
   !! Beware of the high-frequency components that may false the results when hgrid is too high.
   subroutine test_functions(geocode,n01,n02,n03,nspden,acell,a_gauss,hx,hy,hz,&
      density,potential,rhopot,pot_ion,offset)
      implicit none
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
      integer, intent(in) :: n01,n02,n03,nspden
      real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz
      real(kind=8), intent(out) :: offset
      real(kind=8), dimension(n01,n02,n03), intent(out) :: pot_ion,potential
      real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: density,rhopot
      !local variables
      integer :: i1,i2,i3,ifx,ify,ifz,i
      real(kind=8) :: x1,x2,x3,length,denval,pi,a2,derf_tt,factor,r,r2
      real(kind=8) :: fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt

    if (trim(geocode) == 'P' .or. trim(geocode)=='W') then

         !parameters for the test functions
         length=acell
         a=0.5d0/a_gauss**2
         !test functions in the three directions
         ifx=5
         ify=5
         ifz=5
         !parameters of the test functions
         ax=length
         ay=length
         az=length
         bx=2.d0!real(nu,kind=8)
         by=2.d0!real(nu,kind=8)
         bz=2.d0

         !!!     !plot of the functions used
         !!!     do i1=1,n03
         !!!        x = hx*real(i1,kind=8)!valid if hy=hz
         !!!        y = hz*real(i1,kind=8) 
         !!!        call functions(x,ax,bx,fx,fx2,ifx)
         !!!        call functions(y,az,bz,fz,fz2,ifz)
         !!!        write(20,*)i1,fx,fx2,fz,fz2
         !!!     end do

         !Initialization of density and potential
         denval=0.d0 !value for keeping the density positive
         do i3=1,n03
            x3 = hz*real(i3-n03/2-1,kind=8)
            call functions(x3,az,bz,fz,fz2,ifz)
            do i2=1,n02
               x2 = hy*real(i2-n02/2-1,kind=8)
               call functions(x2,ay,by,fy,fy2,ify)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2-1,kind=8)
                  call functions(x1,ax,bx,fx,fx2,ifx)
                  do i=1,nspden
                     density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)
                  end do
                  potential(i1,i2,i3) = -16.d0*datan(1.d0)*fx*fy*fz
                  denval=max(denval,-density(i1,i2,i3,1))
               end do
            end do
         end do

         denval=0.d0

      else if (trim(geocode) == 'S') then
         !parameters for the test functions
         length=acell
         a=0.5d0/a_gauss**2
         !test functions in the three directions
         ifx=5
         ifz=5
         !non-periodic dimension
         ify=6
         !parameters of the test functions
         ax=length
         az=length
         bx=2.d0!real(nu,kind=8)
         bz=2.d0!real(nu,kind=8)
         !non-periodic dimension
         ay=length
         by=a

         !Initialisation of density and potential
         denval=0.d0 !value for keeping the density positive
         do i3=1,n03
            x3 = hz*real(i3-n03/2-1,kind=8)
            call functions(x3,az,bz,fz,fz2,ifz)
            do i2=1,n02
               x2 = hy*real(i2-n02/2-1,kind=8)
               call functions(x2,ay,by,fy,fy2,ify)
               do i1=1,n01
                  x1 = hx*real(i1-n02/2-1,kind=8)
                  call functions(x1,ax,bx,fx,fx2,ifx)
                  do i=1,nspden
                     density(i1,i2,i3,i) = &
                     1.d0/real(nspden,kind=8)/(16.d0*datan(1.d0))*(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)
                  end do
                  potential(i1,i2,i3) = -fx*fy*fz
                  denval=max(denval,-density(i1,i2,i3,1))
               end do
            end do
         end do

         denval=0.d0

      else if (trim(geocode) == 'F') then

         !grid for the free BC case
         !hgrid=max(hx,hy,hz)

         pi = 4.d0*atan(1.d0)
         a2 = a_gauss**2

         !Normalization
         factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
         !gaussian function
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
                  r2 = x1*x1+x2*x2+x3*x3
                  do i=1,nspden
                     density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*max(factor*exp(-r2/a2),1d-24)
                  end do
                  r = sqrt(r2)
                  !Potential from a gaussian
                  if (r == 0.d0) then
                     potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
                  else
                     call derf_local(derf_tt,r/a_gauss)
                     potential(i1,i2,i3) = derf_tt/r
                  end if
               end do
            end do
         end do

         denval=0.d0

      else

         print *,'geometry code not admitted',geocode
         stop

      end if

      ! For ixc/=0 the XC potential is added to the solution, and an analytic comparison is no more
      ! possible. In that case the only possible comparison is between the serial and the parallel case
      ! To ease the comparison between the serial and the parallel case we add a random pot_ion
      ! to the potential.

      if (denval /= 0.d0) then
         rhopot(:,:,:,:) = density(:,:,:,:) + denval +1.d-14
      else
         rhopot(:,:,:,:) = density(:,:,:,:) 
      end if

      offset=0.d0
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
               pot_ion(i1,i2,i3)=tt
               offset=offset+potential(i1,i2,i3)
               !add the case for offset in the surfaces case 
               !(for periodic case it is absorbed in offset)
               if (geocode == 'S' .and. denval /= 0.d0) then
                  x2 = hy*real(i2-1,kind=8)-0.5d0*acell+0.5d0*hy
                  potential(i1,i2,i3)=potential(i1,i2,i3)&
                  -8.d0*datan(1.d0)*denval*real(nspden,kind=8)*(x2**2+0.25d0*acell**2)
                  !this stands for
                  !denval*2pi*Lx*Lz/Ly^2(y^2-Ly^2/4), less accurate in hgrid
               end if

               !!!           if (rhopot(i1,i2,i3,1) <= 0.d0) then
               !!!              print *,i1,i2,i3,rhopot(i1,i2,i3,1),denval
               !!!           end if
            end do
         end do
      end do
      if (denval /= 0.d0) density=rhopot
      offset=offset*hx*hy*hz

      !print *,'offset',offset

   END SUBROUTINE test_functions


   subroutine functions(x,a,b,f,f2,whichone)
      implicit none
      integer, intent(in) :: whichone
      real(kind=8), intent(in) :: x,a,b
      real(kind=8), intent(out) :: f,f2
      !local variables
      real(kind=8) :: r,r2,y,yp,ys,factor,pi,g,h,g1,g2,h1,h2
      real(kind=8) :: length,frequency,nu,sigma,agauss

      pi = 4.d0*datan(1.d0)
      select case(whichone)
      case(1)
         !constant
         f=1.d0
         f2=0.d0
      case(2)
         !gaussian of sigma s.t. a=1/(2*sigma^2)
         r2=a*x**2
         f=dexp(-r2)
         f2=(-2.d0*a+4.d0*a*r2)*dexp(-r2)
      case(3)
         !gaussian "shrinked" with a=length of the system
         length=a
         r=pi*x/length
         y=dtan(r)
         yp=pi/length*1.d0/(dcos(r))**2
         ys=2.d0*pi/length*y*yp
         factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
         f2=factor*dexp(-y**2)
         f=dexp(-y**2)
      case(4)
         !cosine with a=length, b=frequency
         length=a
         frequency=b
         r=frequency*pi*x/length
         f=dcos(r)
         f2=-(frequency*pi/length)**2*dcos(r)
      case(5)
         !exp of a cosine, a=length
         nu=2.d0
         r=pi*nu/a*x
         y=dcos(r)
         yp=dsin(r)
         f=dexp(y)
         factor=(pi*nu/a)**2*(-y+yp**2)
         f2=factor*f
      case(6)
         !gaussian times "shrinked" gaussian, sigma=length/10
         length=a
         r=pi*x/length
         y=dtan(r)
         yp=pi/length*1.d0/(dcos(r))**2
         ys=2.d0*pi/length*y*yp
         factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
         g=dexp(-y**2)
         g1=-2.d0*y*yp*g
         g2=factor*dexp(-y**2)

         sigma=length/10
         agauss=0.5d0/sigma**2
         r2=agauss*x**2
         h=dexp(-r2)
         h1=-2.d0*agauss*x*h
         h2=(-2.d0*agauss+4.d0*agauss*r2)*dexp(-r2)
         f=g*h
         f2=g2*h+g*h2+2.d0*g1*h1
      case(7)
         !sine with a=length, b=frequency
         length=a
         frequency=b
         r=frequency*pi*x/length
         f=dsin(r)
         f2=-(frequency*pi/length)**2*dsin(r)
      end select

      END SUBROUTINE functions


   END PROGRAM PS_Check


!> Error function in double precision
subroutine derf_local(derf_yy,yy)
  use poisson_solver
 implicit none
 real(dp),intent(in) :: yy
 real(dp),intent(out) :: derf_yy
 integer          ::  done,ii,isw
 real(dp), parameter :: &
       ! coefficients for 0.0 <= yy < .477
       &  pp(5)=(/ 113.8641541510502e0_dp, 377.4852376853020e0_dp,  &
       &           3209.377589138469e0_dp, .1857777061846032e0_dp,  &
       &           3.161123743870566e0_dp /)
  real(dp), parameter :: &
       &  qq(4)=(/ 244.0246379344442e0_dp, 1282.616526077372e0_dp,  &
       &           2844.236833439171e0_dp, 23.60129095234412e0_dp/)
  ! coefficients for .477 <= yy <= 4.0
  real(dp), parameter :: &
       &  p1(9)=(/ 8.883149794388376e0_dp, 66.11919063714163e0_dp,  &
       &           298.6351381974001e0_dp, 881.9522212417691e0_dp,  &
       &           1712.047612634071e0_dp, 2051.078377826071e0_dp,  &
       &           1230.339354797997e0_dp, 2.153115354744038e-8_dp, &
       &           .5641884969886701e0_dp /)
  real(dp), parameter :: &
       &  q1(8)=(/ 117.6939508913125e0_dp, 537.1811018620099e0_dp,  &
       &           1621.389574566690e0_dp, 3290.799235733460e0_dp,  &
       &           4362.619090143247e0_dp, 3439.367674143722e0_dp,  &
       &           1230.339354803749e0_dp, 15.74492611070983e0_dp/)
  ! coefficients for 4.0 < y,
  real(dp), parameter :: &
       &  p2(6)=(/ -3.603448999498044e-01_dp, -1.257817261112292e-01_dp,   &
       &           -1.608378514874228e-02_dp, -6.587491615298378e-04_dp,   &
       &           -1.631538713730210e-02_dp, -3.053266349612323e-01_dp/)
  real(dp), parameter :: &
       &  q2(5)=(/ 1.872952849923460e0_dp   , 5.279051029514284e-01_dp,    &
       &           6.051834131244132e-02_dp , 2.335204976268692e-03_dp,    &
       &           2.568520192289822e0_dp /)
  real(dp), parameter :: &
       &  sqrpi=.5641895835477563e0_dp, xbig=13.3e0_dp, xlarge=6.375e0_dp, xmin=1.0e-10_dp
  real(dp) ::  res,xden,xi,xnum,xsq,xx

 xx = yy
 isw = 1
!Here change the sign of xx, and keep track of it thanks to isw
 if (xx<0.0e0_dp) then
  isw = -1
  xx = -xx
 end if

 done=0

!Residual value, if yy < -6.375e0_dp
 res=-1.0e0_dp

!abs(yy) < .477, evaluate approximation for erfc
 if (xx<0.477e0_dp) then
! xmin is a very small number
  if (xx<xmin) then
   res = xx*pp(3)/qq(3)
  else
   xsq = xx*xx
   xnum = pp(4)*xsq+pp(5)
   xden = xsq+qq(4)
   do ii = 1,3
    xnum = xnum*xsq+pp(ii)
    xden = xden*xsq+qq(ii)
   end do
   res = xx*xnum/xden
  end if
  if (isw==-1) res = -res
  done=1
 end if

!.477 < abs(yy) < 4.0 , evaluate approximation for erfc
 if (xx<=4.0e0_dp .and. done==0 ) then
  xsq = xx*xx
  xnum = p1(8)*xx+p1(9)
  xden = xx+q1(8)
  do ii=1,7
   xnum = xnum*xx+p1(ii)
   xden = xden*xx+q1(ii)
  end do
  res = xnum/xden
  res = res* exp(-xsq)
  if (isw.eq.-1) then
     res = res-1.0e0_dp
  else
     res=1.0e0_dp-res
  end if
  done=1
 end if

!y > 13.3e0_dp
 if (isw > 0 .and. xx > xbig .and. done==0 ) then
  res = 1.0e0_dp
  done=1
 end if

!4.0 < yy < 13.3e0_dp  .or. -6.375e0_dp < yy < -4.0
!evaluate minimax approximation for erfc
 if ( ( isw > 0 .or. xx < xlarge ) .and. done==0 ) then
  xsq = xx*xx
  xi = 1.0e0_dp/xsq
  xnum= p2(5)*xi+p2(6)
  xden = xi+q2(5)
  do ii = 1,4
   xnum = xnum*xi+p2(ii)
   xden = xden*xi+q2(ii)
  end do
  res = (sqrpi+xi*xnum/xden)/xx
  res = res* exp(-xsq)
  if (isw.eq.-1) then
     res = res-1.0e0_dp
  else
     res=1.0e0_dp-res
  end if
 end if

!All cases have been investigated
 derf_yy = res

end subroutine derf_local
