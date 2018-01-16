!> @file
!!    Routines to create the kernel for Poisson solver
!! @author
!!    Copyright (C) 2002-2017 BigDFT group  (LG)<br/>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Initialization of the Poisson kernel
function pkernel_init_old(verb,iproc,nproc,igpu,geocode,ndims,hgrids,itype_scf,&
     alg,cavity,mu0_screening,angrad,mpi_env,taskgroup_size) result(kernel)
  use yaml_output
  use yaml_strings, only: f_strcpy
  use f_precisions, only: f_loc
  implicit none
  logical, intent(in) :: verb       !< Verbosity
  integer, intent(in) :: itype_scf  !< Type of interpolating scaling function
  integer, intent(in) :: iproc      !< Proc Id
  integer, intent(in) :: nproc      !< Number of processes
  integer, intent(in) :: igpu
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, dimension(3), intent(in) :: ndims
  real(gp), dimension(3), intent(in) :: hgrids
  !> algorithm of the Solver. Might accept the values "VAC" (default), "PCG", or "PI"
  character(len=*), intent(in), optional :: alg
  !> cavity used, possible values "none" (default), "rigid", "sccs"
  character(len=*), intent(in), optional :: cavity
  real(kind=8), intent(in), optional :: mu0_screening
  real(gp), dimension(3), intent(in), optional :: angrad
  type(mpi_environment), intent(in), optional :: mpi_env
  integer, intent(in), optional :: taskgroup_size
  type(coulomb_operator) :: kernel
  !local variables
  real(dp) :: alphat,betat,gammat,mu0t
  integer :: nthreads,group_size
  !integer :: ierr
  !$ integer :: omp_get_max_threads


  !nullification
  kernel=pkernel_null()

  !geocode and ISF family
  kernel%geocode=geocode
  !dimensions and grid spacings
  kernel%ndims=ndims
  kernel%hgrids=hgrids

  if (present(angrad)) then
!     kernel%angrad=angrad
  else
     alphat = 2.0_dp*datan(1.0_dp)
     betat = 2.0_dp*datan(1.0_dp)
     gammat = 2.0_dp*datan(1.0_dp)
!     kernel%angrad=(/alphat,betat,gammat/)
  end if


  !old approach of input variables, before dictionary
  if (.not. present(mu0_screening)) then
     mu0t=0.0_gp
  else
     mu0t=mu0_screening
  end if



  kernel%mu=mu0t

  if (present(alg)) then
     select case(trim(alg))
     case('VAC')
        kernel%method=PS_VAC_ENUM
     case('PI')
        kernel%method=PS_PI_ENUM
        kernel%nord=16
        !here the parameters can be specified from command line
        kernel%max_iter=50
        kernel%minres=1.0e-8_dp!
        kernel%PI_eta=0.6_dp
     case('PCG')
        kernel%method=PS_PCG_ENUM
        kernel%nord=16
        kernel%max_iter=50
        kernel%minres=1.0e-8_dp!
     case default
        call f_err_throw('Error, kernel algorithm '//trim(alg)//&
             'not valid')
     end select
  else
     kernel%method=PS_VAC_ENUM
  end if

  if (present(cavity)) then
     select case(trim(cavity))
     case('vacuum')
        call f_enum_attr(kernel%method,PS_NONE_ENUM)
     case('rigid')
        call f_enum_attr(kernel%method,PS_RIGID_ENUM)
     case('sccs')
        call f_enum_attr(kernel%method,PS_SCCS_ENUM)
     case default
        call f_err_throw('Error, cavity method '//trim(cavity)//&
             ' not valid')
     end select
  else
     call f_enum_attr(kernel%method,PS_NONE_ENUM)
  end if

  kernel%itype_scf=itype_scf


  !gpu acceleration
  kernel%igpu=igpu

  kernel%initCufftPlan = 1
  kernel%keepGPUmemory = 1
  kernel%keepzf = 1

  if (iproc == 0 .and. verb) then
     if (mu0t==0.0_gp) then
        call yaml_comment('Kernel Initialization',hfill='-')
        call yaml_mapping_open('Poisson Kernel Initialization')
     else
        call yaml_mapping_open('Helmholtz Kernel Initialization')
         call yaml_map('Screening Length (AU)',1/mu0t,fmt='(g25.17)')
     end if
  end if

  group_size=nproc
  if (present(taskgroup_size)) then
     !if the taskgroup size is not a divisor of nproc do not create taskgroups
     if (nproc >1 .and. taskgroup_size > 0 .and. taskgroup_size < nproc .and.&
          mod(nproc,taskgroup_size)==0) then
        group_size=taskgroup_size
     end if
  end if

  !import the mpi_environment if present
  if (present(mpi_env)) then
     call copy_mpi_environment(src=mpi_env,dest=kernel%mpi_env)
     !ernel%mpi_env=mpi_env
  else
     call mpi_environment_set(kernel%mpi_env,iproc,nproc,MPI_COMM_WORLD,group_size)
  end if

  !gpu can be used only for one nproc
  if (nproc > 1) kernel%igpu=0

  !-------------------
  nthreads=0
  if (kernel%mpi_env%iproc == 0 .and. kernel%mpi_env%igroup == 0 .and. verb) then
     !$ nthreads = omp_get_max_threads()
     call yaml_map('MPI tasks',kernel%mpi_env%nproc)
     if (nthreads /=0) call yaml_map('OpenMP threads per MPI task',nthreads)
     if (kernel%igpu==1) call yaml_map('Kernel copied on GPU',.true.)
     if (kernel%method /= 'VAC') call yaml_map('Iterative method for Generalised Equation',str(kernel%method))
     if (kernel%method .hasattr. PS_RIGID_ENUM) call yaml_map('Cavity determination','rigid')
     if (kernel%method .hasattr. PS_SCCS_ENUM) call yaml_map('Cavity determination','sccs')
     call yaml_mapping_close() !kernel
  end if

end function pkernel_init_old


!> Allocate a pointer which corresponds to the zero-padded FFT slice needed for
!! calculating the convolution with the kernel expressed in the interpolating scaling
!! function basis. The kernel pointer is unallocated on input, allocated on output.
subroutine pkernel_set(kernel,eps,dlogeps,oneoeps,oneosqrteps,corr,verbose) !optional arguments
  use yaml_output
  use dynamic_memory
  use time_profiling, only: f_timing
  use dictionaries, only: f_err_throw
  use yaml_strings, only: operator(+)
  use numerics
  implicit none
  !Arguments
  type(coulomb_operator), intent(inout) :: kernel
  !> dielectric function. Needed for non VAC methods, given in full dimensions
  real(dp), dimension(:,:,:), intent(in), optional :: eps
  !> logarithmic derivative of epsilon. Needed for PCG method.
  !! if absent, it will be calculated from the array of epsilon
  real(dp), dimension(:,:,:,:), intent(in), optional :: dlogeps
  !> inverse of epsilon. Needed for PI method.
  !! if absent, it will be calculated from the array of epsilon
  real(dp), dimension(:,:,:), intent(in), optional :: oneoeps
  !> inverse square root of epsilon. Needed for PCG method.
  !! if absent, it will be calculated from the array of epsilon
  real(dp), dimension(:,:,:), intent(in), optional :: oneosqrteps
  !> correction term of the Generalized Laplacian
  !! if absent, it will be calculated from the array of epsilon
  real(dp), dimension(:,:,:), intent(in), optional :: corr
  real(dp) :: alpha
  logical, intent(in), optional :: verbose
  !local variables
  logical :: dump,wrtmsg
  character(len=*), parameter :: subname='createKernel'
  integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,i_stat
  integer :: jproc,nlimd,nlimk,jfd,jhd,jzd,jfk,jhk,jzk,npd,npk
  real(kind=8) :: mu0t
  real(kind=8), dimension(:), allocatable :: pkernel2
  integer :: i1,i2,i3,j1,j2,j3,ind,indt,switch_alg,size2,sizek,kernelnproc,size3
  integer :: n3pr1,n3pr2,istart,jend,i23,i3s,n23,displ,gpuPCGRed
  integer :: myiproc_node, mynproc_node
  integer,dimension(3) :: n
  !call timing(kernel%mpi_env%iproc+kernel%mpi_env%igroup*kernel%mpi_env%nproc,'PSolvKernel   ','ON')
  call f_timing(TCAT_PSOLV_KERNEL,'ON')
  call f_routine(id='pkernel_set')
  !pi=4.0_dp*atan(1.0_dp)
  wrtmsg=.true.
  if (present(verbose)) wrtmsg=verbose

  dump=wrtmsg .and. kernel%mpi_env%iproc+kernel%mpi_env%igroup==0

  mu0t=kernel%mu

  if (dump) then
     if (mu0t==0.0_gp) then
        call yaml_mapping_open('Poisson Kernel Creation')
     else
        call yaml_mapping_open('Helmholtz Kernel Creation')
        call yaml_map('Screening Length (AU)',1/mu0t,fmt='(g25.17)')
     end if
  end if

  kernelnproc=kernel%mpi_env%nproc
  if (kernel%igpu == 1) kernelnproc=1
  kernel%stay_on_gpu=0

  select case(kernel%geocode)
     !if (kernel%geocode == 'P') then
  case('P')

     kernel%geo=[1,1,1]

     if (dump) then
        call yaml_map('Boundary Conditions','Periodic')
     end if
     call P_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,kernel%mpi_env%nproc,.false.)

     if (kernel%igpu > 0) then
       kernel%kernel = f_malloc_ptr((n1/2+1)*n2*n3/kernelnproc,id='kernel%kernel')
     else
       kernel%kernel = f_malloc_ptr(nd1*nd2*nd3/kernelnproc,id='kernel%kernel')
     endif
     !!! PSolver n1-n2 plane mpi partitioning !!!
     call inplane_partitioning(kernel%mpi_env,md2,n2,n3/2+1,kernel%part_mpi,kernel%inplane_mpi,n3pr1,n3pr2)
!!$     if (kernel%mpi_env%nproc>2*(n3/2+1)-1) then
!!$       n3pr1=kernel%mpi_env%nproc/(n3/2+1)
!!$       n3pr2=n3/2+1
!!$       md2plus=.false.
!!$       if ((md2/kernel%mpi_env%nproc)*n3pr1*n3pr2 < n2) then
!!$           md2plus=.true.
!!$       endif
!!$
!!$       if (kernel%mpi_env%iproc==0 .and. n3pr1>1) call yaml_map('PSolver n1-n2 plane mpi partitioning activated:',&
!!$          trim(yaml_toa(n3pr1,fmt='(i5)'))//' x'//trim(yaml_toa(n3pr2,fmt='(i5)'))//&
!!$          ' taskgroups')
!!$       if (kernel%mpi_env%iproc==0 .and. md2plus) &
!!$            call yaml_map('md2 was enlarged for PSolver n1-n2 plane mpi partitioning, md2=',md2)
!!$
!!$       !!$omp master !commented out, no parallel region
!!$       if (n3pr1>1) call mpi_environment_set1(kernel%inplane_mpi,kernel%mpi_env%iproc,kernel%mpi_env%nproc, &
!!$                                             kernel%mpi_env%mpi_comm,n3pr1,n3pr2)
!!$       !!$omp end master
!!$       !!$omp barrier
!!$     else
!!$       n3pr1=1
!!$       n3pr2=kernel%mpi_env%nproc
!!$     endif
!!$
!!$     !!$omp master
!!$     call mpi_environment_set(kernel%part_mpi,kernel%mpi_env%iproc,kernel%mpi_env%nproc,kernel%mpi_env%mpi_comm,n3pr2)
!!$     !!$omp end master
!!$     !!$omp barrier

     ! n3pr1, n3pr2 are sent to Free_Kernel subroutine below
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call Periodic_Kernel(n1,n2,n3,nd1,nd2,nd3,&
          kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3),&
          kernel%itype_scf,kernel%kernel,kernel%mpi_env%iproc,kernelnproc,&
          !mu0t,alphat,betat,gammat,&
          n3pr2,n3pr1)

     nlimd=n2
     nlimk=n3/2+1

     !no powers of hgrid because they are incorporated in the plane wave treatment
     kernel%grid%scal=1.0_dp/(real(n1,dp)*real(n2*n3,dp)) !to reduce chances of integer overflow


  !else if (kernel%geocode == 'S') then
  case('S')

     kernel%geo=[1,0,1]

     if (dump) then
        call yaml_map('Boundary Conditions','Surface')
     end if
     !Build the Kernel
     call S_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,kernel%mpi_env%nproc,&
          kernel%igpu,.false.,non_ortho=.not. kernel%mesh%orthorhombic)


     if (kernel%igpu > 0) then
       kernel%kernel = f_malloc_ptr((n1/2+1)*n2*n3/kernelnproc,id='kernel%kernel')
     else
       kernel%kernel = f_malloc_ptr(nd1*nd2*nd3/kernelnproc,id='kernel%kernel')
     endif

     !!! PSolver n1-n2 plane mpi partitioning !!!
     call inplane_partitioning(kernel%mpi_env,md2,n2,n3/2+1,kernel%part_mpi,kernel%inplane_mpi,n3pr1,n3pr2)

!!$     if (kernel%mpi_env%nproc>2*(n3/2+1)-1) then
!!$       n3pr1=kernel%mpi_env%nproc/(n3/2+1)
!!$       n3pr2=n3/2+1
!!$       md2plus=.false.
!!$       if ((md2/kernel%mpi_env%nproc)*n3pr1*n3pr2 < n2) then
!!$           md2plus=.true.
!!$       endif
!!$
!!$       if (kernel%mpi_env%iproc==0 .and. n3pr1>1) &
!!$            call yaml_map('PSolver n1-n2 plane mpi partitioning activated:',&
!!$            trim(yaml_toa(n3pr1,fmt='(i5)'))//' x'//trim(yaml_toa(n3pr2,fmt='(i5)'))//&
!!$            ' taskgroups')
!!$       if (kernel%mpi_env%iproc==0 .and. md2plus) &
!!$            call yaml_map('md2 was enlarged for PSolver n1-n2 plane mpi partitioning, md2=',md2)
!!$
!!$       !$omp master
!!$       if (n3pr1>1) &
!!$            call mpi_environment_set1(kernel%inplane_mpi,kernel%mpi_env%iproc,kernel%mpi_env%nproc, &
!!$            kernel%mpi_env%mpi_comm,n3pr1,n3pr2)
!!$       !$omp end master
!!$       !$omp barrier
!!$     else
!!$       n3pr1=1
!!$       n3pr2=kernel%mpi_env%nproc
!!$     endif
!!$
!!$     !$omp master
!!$     call mpi_environment_set(kernel%part_mpi,kernel%mpi_env%iproc,&
!!$          kernel%mpi_env%nproc,kernel%mpi_env%mpi_comm,n3pr2)
!!$     !$omp end master
!!$     !$omp barrier

     ! n3pr1, n3pr2 are sent to Free_Kernel subroutine below
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     !the kernel must be built and scattered to all the processes
     call Surfaces_Kernel(kernel%mpi_env%iproc,kernelnproc,&
          kernel%mpi_env%mpi_comm,kernel%inplane_mpi%mpi_comm,&
          n1,n2,n3,m3,nd1,nd2,nd3,&
          kernel%mesh,&
          kernel%itype_scf,kernel%kernel,mu0t)!,alphat)!,betat,gammat)!,n3pr2,n3pr1)

     !last plane calculated for the density and the kernel
     nlimd=n2
     nlimk=n3/2+1

     !only one power of hgrid
     !factor of -4*pi for the definition of the Poisson equation
     kernel%grid%scal=-16.0_dp*atan(1.0_dp)*real(kernel%hgrids(2),dp)/real(n1*n2,dp)/real(n3,dp)

  !else if (kernel%geocode == 'F') then
  case('F')
     kernel%geo=[0,0,0]

     if (dump) then
        call yaml_map('Boundary Conditions','Free')
     end if
!     print *,'debug',kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3)
     !Build the Kernel
     call F_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
          md1,md2,md3,nd1,nd2,nd3,kernel%mpi_env%nproc,kernel%igpu,.false.)

     if (kernel%igpu > 0) then
       kernel%kernel = f_malloc_ptr((n1/2+1)*n2*n3/kernelnproc,id='kernel%kernel')
     else
       !allocate(kernel%kernel(nd1*nd2*nd3/kernelnproc+ndebug),stat=i_stat)
       kernel%kernel = f_malloc_ptr(nd1*nd2*(nd3/kernelnproc),id='kernel%kernel')
     endif

     !!! PSolver n1-n2 plane mpi partitioning !!!
     call inplane_partitioning(kernel%mpi_env,md2,n2/2,n3/2+1,kernel%part_mpi,kernel%inplane_mpi,n3pr1,n3pr2)
!!$     if (kernel%mpi_env%nproc>2*(n3/2+1)-1) then
!!$       n3pr1=kernel%mpi_env%nproc/(n3/2+1)
!!$       n3pr2=n3/2+1
!!$       md2plus=.false.
!!$       if ((md2/kernel%mpi_env%nproc)*n3pr1*n3pr2 < n2/2) then
!!$           md2plus=.true.
!!$       endif
!!$
!!$       if (kernel%mpi_env%iproc==0 .and. n3pr1>1) call yaml_map('PSolver n1-n2 plane mpi partitioning activated:',&
!!$          trim(yaml_toa(n3pr1,fmt='(i5)'))//' x'//trim(yaml_toa(n3pr2,fmt='(i5)'))//&
!!$          ' taskgroups')
!!$       if (kernel%mpi_env%iproc==0 .and. md2plus) &
!!$            call yaml_map('md2 was enlarged for PSolver n1-n2 plane mpi partitioning, md2=',md2)
!!$
!!$       !$omp master
!!$       if (n3pr1>1) call mpi_environment_set1(kernel%inplane_mpi,kernel%mpi_env%iproc,kernel%mpi_env%nproc, &
!!$                                             kernel%mpi_env%mpi_comm,n3pr1,n3pr2)
!!$       !$omp end master
!!$       !$omp barrier
!!$     else
!!$       n3pr1=1
!!$       n3pr2=kernel%mpi_env%nproc
!!$     endif
!!$     !$omp master
!!$     call mpi_environment_set(kernel%part_mpi,kernel%mpi_env%iproc,kernel%mpi_env%nproc,kernel%mpi_env%mpi_comm,n3pr2)
!!$     !$omp end master
!!$     !$omp barrier

     ! n3pr1, n3pr2 are sent to Free_Kernel subroutine below
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !the kernel must be built and scattered to all the processes
     call Free_Kernel(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          n1,n2,n3,nd1,nd2,nd3,kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3),&
          kernel%itype_scf,kernel%mpi_env%iproc,kernelnproc,kernel%kernel,mu0t,n3pr2,n3pr1)

     !last plane calculated for the density and the kernel
     nlimd=n2/2
     nlimk=n3/2+1

     !hgrid=max(hx,hy,hz)
     kernel%grid%scal=product(kernel%hgrids)/real(n1*n2,dp)/real(n3,dp)

  !else if (kernel%geocode == 'W') then
  case('W')

     kernel%geo=[0,0,1]

     if (dump) then
        call yaml_map('Boundary Conditions','Wire')
     end if
     call W_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,kernel%mpi_env%nproc,kernel%igpu,.false.)

     if (kernel%igpu > 0) then
       kernel%kernel = f_malloc_ptr((n1/2+1)*n2*n3/kernelnproc,id='kernel%kernel')
     else
       kernel%kernel = f_malloc_ptr(nd1*nd2*(nd3/kernelnproc),id='kernel%kernel')
     endif

     !!! PSolver n1-n2 plane mpi partitioning !!!
     call inplane_partitioning(kernel%mpi_env,md2,n2,n3/2+1,kernel%part_mpi,kernel%inplane_mpi,n3pr1,n3pr2)

!!$     if (kernel%mpi_env%nproc>2*(n3/2+1)-1) then
!!$       n3pr1=kernel%mpi_env%nproc/(n3/2+1)
!!$       n3pr2=n3/2+1
!!$       md2plus=.false.
!!$       if ((md2/kernel%mpi_env%nproc)*n3pr1*n3pr2 < n2) then
!!$           md2plus=.true.
!!$       endif
!!$
!!$       if (kernel%mpi_env%iproc==0 .and. n3pr1>1) call yaml_map('PSolver n1-n2 plane mpi partitioning activated:',&
!!$          trim(yaml_toa(n3pr1,fmt='(i5)'))//' x'//trim(yaml_toa(n3pr2,fmt='(i5)'))//&
!!$          ' taskgroups')
!!$       if (md2plus) call yaml_map('md2 was enlarged for PSolver n1-n2 plane mpi partitioning, md2=',md2)
!!$
!!$       !$omp master
!!$       if (n3pr1>1) call mpi_environment_set1(kernel%inplane_mpi,kernel%mpi_env%iproc,kernel%mpi_env%nproc, &
!!$                                             kernel%mpi_env%mpi_comm,n3pr1,n3pr2)
!!$       !$omp end master
!!$       !$omp barrier
!!$     else
!!$       n3pr1=1
!!$       n3pr2=kernel%mpi_env%nproc
!!$     endif
!!$
!!$     !$omp master
!!$     call mpi_environment_set(kernel%part_mpi,kernel%mpi_env%iproc,kernel%mpi_env%nproc,kernel%mpi_env%mpi_comm,n3pr2)
!!$     !$omp end master
!!$     !$omp barrier
!!$
!!$     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     call Wires_Kernel(kernel%mpi_env%iproc,kernelnproc,&
          kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          n1,n2,n3,nd1,nd2,nd3,kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3),&
          kernel%itype_scf,kernel%kernel,mu0t)

     nlimd=n2
     nlimk=n3/2+1

     !only one power of hgrid
     !factor of -1/(2pi) already included in the kernel definition
     kernel%grid%scal=-2.0_dp*kernel%hgrids(1)*kernel%hgrids(2)/real(n1*n2,dp)/real(n3,dp)

  case default
     !if (iproc==0)
     !write(*,'(1x,a,3a)')'createKernel, geocode not admitted',kernel%geocode
     call f_err_throw('createKernel, geocode '//trim(kernel%geocode)//&
          'not admitted')
  end select
  !print *,'thereAAA',iproc,nproc,kernel%mpi_env%iproc,kernel%nproc,kernel%mpi_env%mpi_comm
!call MPI_BARRIER(kernel%mpi_env%mpi_comm,ierr)

  if (dump) then
     call yaml_mapping_open('Memory Requirements per MPI task')
       call yaml_map('Density (MB)',8.0_gp*real(md1*md3,gp)*real(md2/kernel%mpi_env%nproc,gp)/(1024.0_gp**2),fmt='(f8.2)')
       call yaml_map('Kernel (MB)',8.0_gp*real(nd1*nd3,gp)*real(nd2/kernel%mpi_env%nproc,gp)/(1024.0_gp**2),fmt='(f8.2)')
       call yaml_map('Full Grid Arrays (MB)',&
            8.0_gp*real(kernel%ndims(1)*kernel%ndims(2),gp)*real(kernel%ndims(3),gp)/(1024.0_gp**2),fmt='(f8.2)')
       !print the load balancing of the different dimensions on screen
     if (kernel%mpi_env%nproc > 1) then
        call yaml_mapping_open('Load Balancing of calculations')
        jhd=10000
        jzd=10000
        npd=0
        load_balancing: do jproc=0,kernel%mpi_env%nproc-1
           !print *,'jproc,jfull=',jproc,jproc*md2/nproc,(jproc+1)*md2/nproc
           if ((jproc+1)*md2/kernel%mpi_env%nproc <= nlimd) then
              jfd=jproc
           else if (jproc*md2/kernel%mpi_env%nproc <= nlimd) then
              jhd=jproc
              npd=nint(real(nlimd-(jproc)*md2/kernel%mpi_env%nproc,kind=8)/real(md2/kernel%mpi_env%nproc,kind=8)*100.d0)
           else
              jzd=jproc
              exit load_balancing
           end if
        end do load_balancing
        call yaml_mapping_open('Density')
         call yaml_map('MPI tasks 0-'//jfd**'(i5)','100%')
         if (jfd < kernel%mpi_env%nproc-1) &
              call yaml_map('MPI task '//jhd**'(i5)',npd**'(i5)'+'%')
         if (jhd < kernel%mpi_env%nproc-1) &
              call yaml_map('MPI tasks'//jhd**'(i5)'+'-'+&
              (kernel%mpi_env%nproc-1)**'(i3)','0%')
        call yaml_mapping_close()
        jhk=10000
        jzk=10000
        npk=0
       ! if (geocode /= 'P') then
           load_balancingk: do jproc=0,kernel%mpi_env%nproc-1
              !print *,'jproc,jfull=',jproc,jproc*nd3/kernel%mpi_env%nproc,(jproc+1)*nd3/kernel%mpi_env%nproc
              if ((jproc+1)*nd3/kernel%mpi_env%nproc <= nlimk) then
                 jfk=jproc
              else if (jproc*nd3/kernel%mpi_env%nproc <= nlimk) then
                 jhk=jproc
                 npk=nint(real(nlimk-(jproc)*nd3/kernel%mpi_env%nproc,kind=8)/real(nd3/kernel%mpi_env%nproc,kind=8)*100.d0)
              else
                 jzk=jproc
                 exit load_balancingk
              end if
           end do load_balancingk
           call yaml_mapping_open('Kernel')
           call yaml_map('MPI tasks 0-'//trim(yaml_toa(jfk,fmt='(i5)')),'100%')
!           print *,'here,npk',npk
           if (jfk < kernel%mpi_env%nproc-1) &
                call yaml_map('MPI task'//trim(yaml_toa(jhk,fmt='(i5)')),trim(yaml_toa(npk,fmt='(i5)'))//'%')
           if (jhk < kernel%mpi_env%nproc-1) &
                call yaml_map('MPI tasks'//trim(yaml_toa(jhk,fmt='(i5)'))+'-'+&
                yaml_toa(kernel%mpi_env%nproc-1,fmt='(i3)'),'0%')
           call yaml_mapping_close()
        call yaml_map('Complete LB per task','1/3 LB_density + 2/3 LB_kernel')
        call yaml_mapping_close()
     end if
     call yaml_mapping_close() !memory

  end if

  if(kernel%keepzf == 1) then
    if(kernel%igpu /= 1) then
    !  kernel%w%zf = f_malloc_ptr([md1, md3, md2/kernel%mpi_env%nproc],id='zf')
    !else
      kernel%w%zf = f_malloc_ptr([md1, md3, 2*md2/kernel%mpi_env%nproc],id='zf')
    end if
  end if

  kernel%gpuPCGRed=0
  if (kernel%igpu >0) then
    if(trim(str(kernel%method))=='PCG') kernel%gpuPCGRed=1
    n(1)=n1!kernel%ndims(1)*(2-kernel%geo(1))
    n(2)=n3!kernel%ndims(2)*(2-kernel%geo(2))
    n(3)=n2!kernel%ndims(3)*(2-kernel%geo(3))
    !perform the estimation of the processors
    call mpinoderanks(kernel%mpi_env%iproc,kernel%mpi_env%nproc,kernel%mpi_env%mpi_comm,&
         myiproc_node,mynproc_node)

    call cuda_estimate_memory_needs(kernel, n, &
         int(myiproc_node,kind=8), int(mynproc_node,kind=8)) !LG: why longs?

    size2=2*n1*n2*n3
    sizek=(n1/2+1)*n2*n3
    size3=n1*n2*n3
  if (kernel%gpuPCGRed==1) then
    if (kernel%keepGPUmemory == 1) then
      call cudamalloc(size3,kernel%w%z_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc z_GPU (GPU out of memory ?) ')
      call cudamalloc(size3,kernel%w%r_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc r_GPU (GPU out of memory ?) ')
      call cudamalloc(size3,kernel%w%oneoeps_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc oneoeps_GPU (GPU out of memory ?) ')
      call cudamalloc(size3,kernel%w%p_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc p_GPU (GPU out of memory ?) ')
      call cudamalloc(size3,kernel%w%q_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc q_GPU (GPU out of memory ?) ')
      call cudamalloc(size3,kernel%w%x_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc x_GPU (GPU out of memory ?) ')
      call cudamalloc(size3,kernel%w%corr_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc corr_GPU (GPU out of memory ?) ')
      call cudamalloc(sizeof(alpha),kernel%w%alpha_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc alpha_GPU (GPU out of memory ?) ')
      call cudamalloc(sizeof(alpha),kernel%w%beta_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc beta_GPU (GPU out of memory ?) ')
      call cudamalloc(sizeof(alpha),kernel%w%beta0_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc beta0_GPU (GPU out of memory ?) ')
      call cudamalloc(sizeof(alpha),kernel%w%kappa_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc kappa_GPU (GPU out of memory ?) ')
    end if
  end if
   if (kernel%mpi_env%iproc == 0) then
     if (kernel%igpu == 1) then
      call cudacreatestream(i_stat)
      if (i_stat /= 0) call f_err_throw('error creating stream ')
      call cudacreatecublashandle()
      if (kernel%keepGPUmemory == 1) then
        call cudamalloc(size2,kernel%w%work1_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc work1_GPU (GPU out of memory ?) ')
        call cudamalloc(size2,kernel%w%work2_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc work2_GPU (GPU out of memory ?) ')
        call cudamalloc(size3,kernel%w%rho_GPU,i_stat)
        if (i_stat /= 0) call f_err_throw('error cudamalloc rho_GPU (GPU out of memory ?) ')
      endif
      call cudamalloc(sizek,kernel%w%k_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc k_GPU (GPU out of memory ?) ')
      call cudamalloc(1,kernel%w%ehart_GPU,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc ehart_GPU (GPU out of memory ?) ')
      call cudamalloc(1,kernel%w%eexctX_GPU,i_stat)
      call cudamemset(kernel%w%eexctX_GPU,0,1,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc eexctX_GPU (GPU out of memory ?) ')
      call cudamalloc(64,kernel%w%reduc_GPU,i_stat)
!      call cudamemset(kernel%w%reduc_GPU,0,64,i_stat)
      if (i_stat /= 0) call f_err_throw('error cudamalloc reduc_GPU (GPU out of memory ?) ')
    endif

    pkernel2 = f_malloc((n1/2+1)*n2*n3,id='pkernel2')

    ! transpose kernel for GPU
    do i3=1,n3
       j3=i3+(i3/(n3/2+2))*(n3+2-2*i3)!injective dimension
       do i2=1,n2
          j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)!injective dimension
          do i1=1,n1
             j1=i1+(i1/(n1/2+2))*(n1+2-2*i1)!injective dimension
             !injective index
             ind=j1+(j2-1)*nd1+(j3-1)*nd1*nd2
             !unfolded index
             indt=i2+(j1-1)*n2+(i3-1)*nd1*n2
             pkernel2(indt)=kernel%kernel(ind)
          end do
       end do
    end do
    !offset to zero
    if (kernel%geocode == 'P') pkernel2(1)=0.0_dp

    if (kernel%igpu == 2) kernel%kernel=pkernel2
   endif

   if (kernel%mpi_env%iproc == 0) then
    if (kernel%igpu == 1) then
      call reset_gpu_data((n1/2+1)*n2*n3,pkernel2,kernel%w%k_GPU)

      if (dump) call yaml_map('Kernel Copied on GPU',.true.)

      if (kernel%initCufftPlan == 1) then
        call cuda_3d_psolver_general_plan(n,kernel%plan,switch_alg,kernel%geo)
      endif
    endif

    call f_free(pkernel2)
 endif

endif

!print *,'there',iproc,nproc,kernel%iproc,kernel%mpi_env%nproc,kernel%mpi_env%mpi_comm
!call MPI_BARRIER(kernel%mpi_comm,ierr)
!print *,'okcomm',kernel%mpi_comm,kernel%iproc
!call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)

  if (dump) call yaml_mapping_close() !kernel

  !here the FFT_metadata routine can be filled
  kernel%grid%m1 =m1
  kernel%grid%m2 =m2
  kernel%grid%m3 =m3
  kernel%grid%n1 =n1
  kernel%grid%n2 =n2
  kernel%grid%n3 =n3
  kernel%grid%md1=md1
  kernel%grid%md2=md2
  kernel%grid%md3=md3
  kernel%grid%nd1=nd1
  kernel%grid%nd2=nd2
  kernel%grid%nd3=nd3
  kernel%grid%istart=kernel%mpi_env%iproc*(md2/kernel%mpi_env%nproc)
  kernel%grid%iend=min((kernel%mpi_env%iproc+1)*md2/kernel%mpi_env%nproc,kernel%grid%m2)
  if (kernel%grid%istart <= kernel%grid%m2-1) then
     kernel%grid%n3p=kernel%grid%iend-kernel%grid%istart
  else
     kernel%grid%n3p=0
  end if

if (kernel%igpu == 0) then
  !add the checks that are done at the beginning of the Poisson Solver
  if (mod(kernel%grid%n1,2) /= 0 .and. kernel%geo(1)==0) &
       call f_err_throw('Parallel convolution:ERROR:n1') !this can be avoided
  if (mod(kernel%grid%n2,2) /= 0 .and. kernel%geo(3)==0) &
       call f_err_throw('Parallel convolution:ERROR:n2') !this can be avoided
  if (mod(kernel%grid%n3,2) /= 0 .and. kernel%geo(2)==0) &
       call f_err_throw('Parallel convolution:ERROR:n3') !this can be avoided
  if (kernel%grid%nd1 < kernel%grid%n1/2+1) call f_err_throw('Parallel convolution:ERROR:nd1')
  if (kernel%grid%nd2 < kernel%grid%n2/2+1) call f_err_throw('Parallel convolution:ERROR:nd2')
  if (kernel%grid%nd3 < kernel%grid%n3/2+1) call f_err_throw('Parallel convolution:ERROR:nd3')
  !these can be relaxed
!!$  if (kernel%grid%md1 < n1dim) call f_err_throw('Parallel convolution:ERROR:md1')
!!$  if (kernel%grid%md2 < n2dim) call f_err_throw('Parallel convolution:ERROR:md2')
!!$  if (kernel%grid%md3 < n3dim) call f_err_throw('Parallel convolution:ERROR:md3')
  if (mod(kernel%grid%nd3,kernel%mpi_env%nproc) /= 0) &
       call f_err_throw('Parallel convolution:ERROR:nd3')
  if (mod(kernel%grid%md2,kernel%mpi_env%nproc) /= 0) &
       call f_err_throw('Parallel convolution:ERROR:md2'+ &
    yaml_toa(kernel%mpi_env%nproc)+yaml_toa(kernel%grid%md2))
end if
  !allocate and set the distributions for the Poisson Solver
  kernel%counts = f_malloc_ptr([0.to.kernel%mpi_env%nproc-1],id='counts')
  kernel%displs = f_malloc_ptr([0.to.kernel%mpi_env%nproc-1],id='displs')
  do jproc=0,kernel%mpi_env%nproc-1
     istart=min(jproc*(md2/kernel%mpi_env%nproc),m2-1)
     jend=max(min(md2/kernel%mpi_env%nproc,m2-md2/kernel%mpi_env%nproc*jproc),0)
     kernel%counts(jproc)=kernel%grid%m1*kernel%grid%m3*jend
     kernel%displs(jproc)=kernel%grid%m1*kernel%grid%m3*istart
  end do

  ! multi-gpu poisson distribution
  if (kernel%igpu>0 .and. kernel%mpi_env%iproc ==0) then
    displ=0
    kernel%rhocounts=f_malloc_ptr([0.to.kernel%mpi_env%nproc-1], id='rhocounts')
    kernel%rhodispls=f_malloc_ptr([0.to.kernel%mpi_env%nproc-1], id='rhodispls')
    do jproc=0,kernel%mpi_env%nproc-1
      kernel%rhodispls(jproc)=displ
      istart=jproc*( kernel%grid%md2/kernel%mpi_env%nproc)
      jend=min((jproc+1)* kernel%grid%md2/kernel%mpi_env%nproc,kernel%grid%m2)
      if (istart <= kernel%grid%m2-1) then
         kernel%rhocounts(jproc)=(jend-istart)*kernel%grid%md3*kernel%grid%md1
      else
         kernel%rhocounts(jproc)=0
      end if
      displ=displ+kernel%rhocounts(jproc)
    end do
  end if

  !allocate cavity if needed
  n1=kernel%ndims(1)
  n23=kernel%ndims(2)*kernel%grid%n3p
  call PS_allocate_cavity_workarrays(n1,n23,kernel%ndims,&
       kernel%method,kernel%w)

  !>>>>set the treatment of the cavity
  select case(trim(str(kernel%method)))
  case('PCG')
     if (present(eps)) then
     if (present(oneosqrteps)) then
        call pkernel_set_epsilon(kernel,eps=eps,oneosqrteps=oneosqrteps)
     else if (present(corr)) then
        call pkernel_set_epsilon(kernel,eps=eps,corr=corr)
     else
        call pkernel_set_epsilon(kernel,eps=eps)
     end if
  else if (present(oneosqrteps) .and. present(corr)) then
     call pkernel_set_epsilon(kernel,oneosqrteps=oneosqrteps,corr=corr)
  else if (present(oneosqrteps) .neqv. present(corr)) then
     call f_err_throw('For PCG method either eps, oneosqrteps and/or corr should be present')
  end if
  case('PI')
     if (present(eps)) then
        if (present(oneoeps)) then
           call pkernel_set_epsilon(kernel,eps=eps,oneoeps=oneoeps)
        else if (present(dlogeps)) then
           call pkernel_set_epsilon(kernel,eps=eps,dlogeps=dlogeps)
        else
           call pkernel_set_epsilon(kernel,eps=eps)
        end if
     else if (present(oneoeps) .and. present(dlogeps)) then
        call pkernel_set_epsilon(kernel,oneoeps=oneoeps,dlogeps=dlogeps)
     else if (present(oneoeps) .neqv. present(dlogeps)) then
        call f_err_throw('For PI method either eps, oneoeps and/or dlogeps should be present')
     end if
  end select

  call f_release_routine()
  call f_timing(TCAT_PSOLV_KERNEL,'OF')

END SUBROUTINE pkernel_set


subroutine cuda_estimate_memory_needs(kernel, n,iproc_node, nproc_node)
  use iso_c_binding
!  use module_base
  implicit none
  !Arguments
  type(coulomb_operator), intent(inout) :: kernel
  integer,dimension(3), intent(in) :: n
  integer(kind=8), intent(in) :: iproc_node,nproc_node
  !Local variables
  integer(kind=C_SIZE_T) :: maxPlanSize, freeGPUSize, totalGPUSize
  integer(kind=8) :: size2,sizek,size3,NX,NY,NZ
  integer(kind=8) :: kernelSize, PCGRedSize, plansSize
  real(dp) alpha

  kernelSize=0
  PCGRedSize=0
  plansSize=0
  maxPlanSize=0
  freeGPUSize=0
  totalGPUSize=0

 !estimate with CUDA the free memory size, and the size of the plans
 call cuda_estimate_memory_needs_cu(kernel%mpi_env%iproc,n,&
    kernel%geo,plansSize,maxPlanSize,freeGPUSize, totalGPUSize )

   size2=2*n(1)*n(2)*n(3)*sizeof(alpha)
   sizek=(n(1)/2+1)*n(2)*n(3)*sizeof(alpha)
   size3=n(1)*n(2)*n(3)*sizeof(alpha)

!only the first MPI process of the group needs the GPU for apply_kernel
 if(kernel%mpi_env%iproc==0) then
   kernelSize =2*size2+size3+sizek
 end if

!all processes can use the GPU for apply_reductions
 if((kernel%gpuPCGRed)==1) then
   !add a 10% margin, because we use a little bit more
   PCGRedSize=int(real(7*size3+4*sizeof(alpha),kind=8)*1.1d0,kind=8)
   !print *,"PCG reductions size : %lu\n", PCGRedSize
 end if


!print *,"free mem",freeGPUSize,", total",totalGPUSize,". Trying Total : ",kernelSize+plansSize+PCGRedSize,&
!" with kernel ",kernelSize," plans ",plansSize, "maxplan",&
!maxPlanSize, "and red ",PCGRedSize, "nprocs/node", nproc_node

 if(freeGPUSize<nproc_node*(kernelSize+maxPlanSize)) then
     if(kernel%mpi_env%iproc==0)then
       call f_err_throw('Not Enough memory on the card to allocate GPU kernels, free Memory :' // &
       trim(yaml_toa(freeGPUSize)) // ", total Memory :"// trim(yaml_toa(totalGPUSize)) //&
       ", minimum needed memory :"// trim(yaml_toa(nproc_node*(kernelSize+maxPlanSize))) )
     end if
 else if(freeGPUSize <nproc_node*(kernelSize+plansSize)) then
     if(kernel%mpi_env%iproc==0) &
       call yaml_warning( "WARNING: not enough free memory for cufftPlans on GPU, performance will be degraded")
     kernel%initCufftPlan=0
     kernel%gpuPCGRed=0
 else if((kernel%gpuPCGRed == 1) .and. (freeGPUSize < nproc_node*(kernelSize+plansSize+PCGRedSize))) then
     if(kernel%mpi_env%iproc==0) &
      call yaml_warning( "WARNING: not enough free memory for GPU PCG reductions, performance will be degraded")
     kernel%gpuPCGRed=0;
 else
     !call yaml_comment("Memory on the GPU is sufficient for" // trim(yaml_toa(nproc_node)) // " processes/node")
     kernel%initCufftPlan=1;
 end if

call mpibarrier()

end subroutine cuda_estimate_memory_needs

!>put in depsdrho array the extra potential
subroutine sccs_extra_potential(kernel,pot,depsdrho,dsurfdrho,eps0)
  use FDder
  use yaml_output
  implicit none
  type(coulomb_operator), intent(in) :: kernel
  !>complete potential, needed to calculate the derivative
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(in) :: pot
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(inout) :: depsdrho
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(in) :: dsurfdrho
  real(dp), intent(in) :: eps0
  !local variables
  integer :: i3,i3s,i2,i1,i23,i,n01,n02,n03,unt
  real(dp) :: d2,x,pi,gammaSau,alphaSau,betaVau
  real(dp), dimension(:,:,:,:), allocatable :: nabla2_pot
  !real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)) :: pot2,depsdrho1,depsdrho2

  gammaSau=kernel%cavity%gammaS*5.291772109217d-9/8.238722514d-3 ! in atomic unit
  alphaSau=kernel%cavity%alphaS*5.291772109217d-9/8.238722514d-3 ! in atomic unit
  betaVau=kernel%cavity%betaV/2.942191219d4 ! in atomic unit
  pi = 4.d0*datan(1.d0)
  n01=kernel%ndims(1)
  n02=kernel%ndims(2)
  n03=kernel%ndims(3)
  !starting point in third direction
  i3s=kernel%grid%istart+1

  nabla2_pot=f_malloc([n01,n02,n03,3],id='nabla_pot')
  !calculate derivative of the potential
  !call nabla_u_square(kernel%geocode,n01,n02,n03,pot,nabla2_pot,kernel%nord,kernel%hgrids)
  call nabla_u(kernel%geocode,n01,n02,n03,pot,nabla2_pot,kernel%nord,kernel%hgrids)

  i23=1
  do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
     do i2=1,n02
        do i1=1,n01
           !this section has to be inserted into a optimized calculation of the derivative
           d2=0.0_dp
           do i=1,3
              d2 = d2+nabla2_pot(i1,i2,i3,i)**2
           end do
!!$           !depsdrho1(i1,i2,i3)=depsdrho(i1,i23)
!!$           d2=nabla2_pot(i1,i2,i3)
           depsdrho(i1,i23)=-0.125d0*depsdrho(i1,i23)*d2/pi!&
                            !+(alphaSau+gammaSau)*dsurfdrho(i1,i23)&
                            !+betaVau*depsdrho(i1,i23)/(1.d0-eps0)
           !depsdrho(i1,i23)=depsdrho(i1,i23)*d2
        end do
        i23=i23+1
     end do
  end do

  call f_free(nabla2_pot)

  if (kernel%mpi_env%iproc==0 .and. kernel%mpi_env%igroup==0) then
       call yaml_map('Extra SCF potential calculated',.true.)
  end if

end subroutine sccs_extra_potential

!>put in pol_charge array the polarization charge
subroutine polarization_charge(kernel,pot,rho)
  use FDder
  implicit none
  type(coulomb_operator), intent(inout) :: kernel
  !>complete potential, needed to calculate the derivative
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(in) :: pot
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(in) :: rho
  !local variables
  integer :: i3,i3s,i2,i1,i23,i,n01,n02,n03
  real(dp) :: d2,pi
  real(dp), dimension(:,:,:,:), allocatable :: nabla_pot
  real(dp), dimension(:,:,:), allocatable :: lapla_pot

  pi=4.0_dp*atan(1.0_dp)
  n01=kernel%ndims(1)
  n02=kernel%ndims(2)
  n03=kernel%ndims(3)
  !starting point in third direction
  i3s=kernel%grid%istart+1

  nabla_pot=f_malloc([n01,n02,n03,3],id='nabla_pot')
  lapla_pot=f_malloc([n01,n02,n03],id='lapla_pot')
  !calculate derivative of the potential

  call nabla_u(kernel%geocode,n01,n02,n03,pot,nabla_pot,kernel%nord,kernel%hgrids)
  call div_u_i(kernel%geocode,n01,n02,n03,nabla_pot,lapla_pot,kernel%nord,kernel%hgrids)
  i23=1
  do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
     do i2=1,n02
        do i1=1,n01
           !this section has to be inserted into a optimized calculation of the
           !derivative
           kernel%w%rho_pol(i1,i23)=(-0.25_dp/pi)*lapla_pot(i1,i2,i3)-rho(i1,i23)
        end do
        i23=i23+1
     end do
  end do

  call f_free(nabla_pot)
  call f_free(lapla_pot)

  if (kernel%mpi_env%iproc==0 .and. kernel%mpi_env%igroup==0) &
       call yaml_map('Polarization charge calculated',.true.)

end subroutine polarization_charge

!>build the needed arrays of the cavity from a given density
!!according to the SCF cavity definition given by Andreussi et al. JCP 136, 064102 (2012)
!! @warning: for the moment the density is supposed to be not distributed as the
!! derivatives are calculated sequentially
subroutine pkernel_build_epsilon(kernel,edens,eps0,depsdrho,dsurfdrho)
  use numerics, only: safe_exp
  use f_utils
  use yaml_output
  use FDder
  implicit none
  !> Poisson Solver kernel
  real(dp), intent(in) :: eps0
  type(coulomb_operator), intent(inout) :: kernel
  !> electronic density in the full box. This is needed because of the calculation of the
  !! gradient
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(inout) :: edens
  !> functional derivative of the sc epsilon with respect to
  !! the electronic density, in distributed memory
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(out) :: depsdrho
  !> functional derivative of the surface integral with respect to
  !! the electronic density, in distributed memory
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(out) :: dsurfdrho
  !local variables
  logical, parameter :: dumpeps=.true.  !.true.
  real(kind=8), parameter :: edensmax = 0.005d0 !0.0050d0
  real(kind=8), parameter :: edensmin = 0.0001d0
  real(kind=8), parameter :: innervalue = 0.9d0
  integer :: n01,n02,n03,i,i1,i2,i3,i23,i3s,unt
  real(dp) :: oneoeps0,oneosqrteps0,pi,coeff,coeff1,fact1,fact2,fact3,r,t,d2,dtx,dd,x,y,z
  real(dp) :: de,d
  !real(dp) :: dde,ddtx,c1,c2
  real(dp), dimension(:,:,:), allocatable :: ddt_edens,epscurr,epsinner,depsdrho1,cc
  real(dp), dimension(:,:,:,:), allocatable :: nabla_edens
  real(dp), parameter :: gammaS = 72.d0 ![dyn/cm]
  real(dp), parameter :: alphaS = -22.0d0 ![dyn/cm]
  real(dp), parameter :: betaV = -0.35d0 ![GPa]
  real(dp) :: gammaSau, alphaSau,betaVau,epsm1,rho,logepspr,surf,zeta
  real(gp) :: IntSur,IntVol,noeleene,Cavene,Repene,Disene

  gammaSau=kernel%cavity%gammaS*5.291772109217d-9/8.238722514d-3 ! in atomic unit
  alphaSau=kernel%cavity%alphaS*5.291772109217d-9/8.238722514d-3 ! in atomic unit
  betaVau=kernel%cavity%betaV/2.942191219d4 ! in atomic unit
  IntSur=0.d0
  IntVol=0.d0

  n01=kernel%ndims(1)
  n02=kernel%ndims(2)
  n03=kernel%ndims(3)
  !starting point in third direction
  i3s=kernel%grid%istart+1

  !allocate the work arrays
  nabla_edens=f_malloc([n01,n02,n03,3],id='nabla_edens')
  ddt_edens=f_malloc(kernel%ndims,id='ddt_edens')
  epsinner=f_malloc(kernel%ndims,id='epsinner')
  depsdrho1=f_malloc(kernel%ndims,id='depsdrho1')
  cc=f_malloc(kernel%ndims,id='cc')
  if (dumpeps) epscurr=f_malloc(kernel%ndims,id='epscurr')

  !build the gradients and the laplacian of the density
  !density gradient in du
  call nabla_u(kernel%geocode,n01,n02,n03,edens,nabla_edens,kernel%nord,kernel%hgrids)
  !density laplacian in d2u
  call div_u_i(kernel%geocode,n01,n02,n03,nabla_edens,ddt_edens,kernel%nord,kernel%hgrids,cc)

  pi = 4.d0*datan(1.d0)
  oneoeps0=1.d0/eps0
  oneosqrteps0=1.d0/dsqrt(eps0)
  fact1=2.d0*pi/(dlog(edensmax)-dlog(edensmin))
  fact2=(dlog(eps0))/(2.d0*pi)
  fact3=(dlog(eps0))/(dlog(edensmax)-dlog(edensmin))

  epsm1=(kernel%cavity%epsilon0-1.0_gp)
  if (kernel%mpi_env%iproc==0 .and. kernel%mpi_env%igroup==0) &
       call yaml_map('Rebuilding the cavity for method',trim(str(kernel%method)))

  !now fill the pkernel arrays according the the chosen method
  !if ( trim(PSol)=='PCG') then
  select case(trim(str(kernel%method)))
  case('PCG')
     !in PCG we only need corr, oneosqrtepsilon
     i23=1
     do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
        !do i3=1,n03
        do i2=1,n02
           do i1=1,n01
              if (kernel%w%epsinnersccs(i1,i23).gt.innervalue) then ! Check for inner sccs cavity value to fix as vacuum
                 !eps(i1,i2,i3)=1.d0
                 if (dumpeps) epscurr(i1,i2,i3)=1.d0
                 kernel%w%oneoeps(i1,i23)=1.d0 !oneosqrteps(i1,i2,i3)
                 kernel%w%corr(i1,i23)=0.d0 !corr(i1,i2,i3)
                 depsdrho(i1,i23)=0.d0
 !                depsdrho1(i1,i2,i3)=0.d0
                 dsurfdrho(i1,i23)=0.d0
              else
                rho=edens(i1,i2,i3)
                d2=0.0_dp
                do i=1,3
                   d2 = d2+nabla_edens(i1,i2,i3,i)**2
                end do
                d=sqrt(d2)
                dd = ddt_edens(i1,i2,i3)
                de=epsprime(rho,kernel%cavity)
                if (dumpeps) epscurr(i1,i2,i3)=eps(rho,kernel%cavity)
                depsdrho(i1,i23)=de
                kernel%w%oneoeps(i1,i23)=oneosqrteps(rho,kernel%cavity)
                kernel%w%corr(i1,i23)=corr_term(rho,d2,dd,kernel%cavity)
!!$                !c1=(cc(i1,i2,i3)/d2-dd)/d
!!$                dsurfdrho(i1,i23)=surf_term(rho,d2,dd,cc(i1,i2,i3),kernel%cavity)/epsm1
!!$                !evaluate surfaces and volume integrals
                IntSur=IntSur + de*d/epsm1
                IntVol=IntVol + (kernel%cavity%epsilon0-eps(rho,kernel%cavity))/epsm1
!!$
!!$
!!$                 if (dabs(edens(i1,i2,i3)).gt.edensmax) then
!!$                    !eps(i1,i2,i3)=1.d0
!!$                    if (dumpeps) epscurr(i1,i2,i3)=1.d0
!!$!                    kernel%oneoeps(i1,i23)=1.d0 !oneosqrteps(i1,i2,i3)
!!$!                    kernel%corr(i1,i23)=0.d0 !corr(i1,i2,i3)
!!$!                    depsdrho(i1,i23)=0.d0
!!$!                    depsdrho1(i1,i2,i3)=0.d0
!!$!                    dsurfdrho(i1,i23)=0.d0
!!$                 else if (dabs(edens(i1,i2,i3)).lt.edensmin) then
!!$                    !eps(i1,i2,i3)=eps0
!!$                    if (dumpeps) epscurr(i1,i2,i3)=eps0
!!$!                    kernel%oneoeps(i1,i23)=oneosqrteps0 !oneosqrteps(i1,i2,i3)
!!$!                    kernel%corr(i1,i23)=0.d0 !corr(i1,i2,i3)
!!$!                    depsdrho(i1,i23)=0.d0
!!$!                    depsdrho1(i1,i2,i3)=0.d0
!!$                    dsurfdrho(i1,i23)=0.d0
!!$                 else
!!$                    r=fact1*(log(edensmax)-log(dabs(edens(i1,i2,i3))))
!!$                    t=fact2*(r-sin(r))
!!$                    !eps(i1,i2,i3)=exp(t)
!!$                    if (dumpeps) epscurr(i1,i2,i3)=safe_exp(t)
!!$
!!$!                    kernel%oneoeps(i1,i23)=safe_exp(-0.5d0*t) !oneosqrteps(i1,i2,i3)
!!$!zeta=2.0_gp*pi*log(kernel%cavity%edensmax/rho)/log(kernel%cavity%edensmax/kernel%cavity%edensmin)
!!$!print *,t,(zeta-sin(zeta))/(2.d0*pi)*log(eps0)
!!$!print *,r,zeta,(zeta-sin(zeta))/(2.d0*pi)
!!$!print *,safe_exp((zeta-sin(zeta))/(2.d0*pi)*log(eps0)),eps0,safe_exp((zeta-sin(zeta))/(2.d0*pi))
!!$!print *,safe_exp(t),eps(rho,kernel%cavity)
!!$! print *,oneosqrteps(rho,kernel%cavity),safe_exp(-0.5d0*t),oneosqrteps(rho,kernel%cavity)-safe_exp(-0.5d0*t) !wrong
!!$                    coeff=fact3*(1.d0-cos(r))
!!$                    dtx=-coeff/dabs(edens(i1,i2,i3))  !first derivative of t wrt rho
!!$                    de=safe_exp(t)*dtx ! derivative of epsilon wrt rho
!!$!                    depsdrho(i1,i23)=de
!!$!                    depsdrho1(i1,i2,i3)=de
!!$                    ddtx=fact3*(1.d0-cos(r)+fact1*sin(r))/((edens(i1,i2,i3))**2) !second derivative of t wrt rho
!!$                    dde=de*dtx+safe_exp(t)*ddtx
!!$                    d2=0.d0
!!$                    do i=1,3
!!$                       !dlogeps(i,i1,i2,i3)=dtx*nabla_edens(i1,i2,i3,isp,i)
!!$                       d2 = d2+nabla_edens(i1,i2,i3,i)**2
!!$                    end do
!!$                    d=dsqrt(d2)
!!$                    dd = ddt_edens(i1,i2,i3)
!!$                    coeff1=(0.5d0*(coeff**2)+fact3*fact1*sin(r)+coeff)/((edens(i1,i2,i3))**2)
!!$!                    kernel%corr(i1,i23)=(0.125d0/pi)*safe_exp(t)*(coeff1*d2+dtx*dd) !corr(i1,i2,i3)
!!$                    c1=(cc(i1,i2,i3)/d2-dd)/d
!!$                    dsurfdrho(i1,i23)=(de*c1)/(eps0-1.d0)
!!$                    !dsurfdrho(i1,i23)=(de*c1+dde*c2)/(eps0-1.d0)
!!$!                    IntSur=IntSur + de*d
!!$!                    IntVol=IntVol + (eps0-safe_exp(t))/(eps0-1.d0)
!!$                 end if

              end if
           end do
           i23=i23+1
        end do
     end do
  case('PI')
     !for PI we need  dlogeps,oneoeps
     !first oneovereps
     i23=1
     do i3=i3s,i3s+kernel%grid%n3p-1!kernel%ndims(3)
        do i2=1,n02
           do i1=1,n01
             epsinner(i1,i2,i3)=kernel%w%epsinnersccs(i1,i23)
             if (kernel%w%epsinnersccs(i1,i23).gt.innervalue) then ! Check for inner sccs cavity value to fix as vacuum
                 !eps(i1,i2,i3)=1.d0
                 if (dumpeps) epscurr(i1,i2,i3)=1.d0
                 kernel%w%oneoeps(i1,i23)=1.d0 !oneoeps(i1,i2,i3)
                 depsdrho(i1,i23)=0.d0
                 dsurfdrho(i1,i23)=0.d0
             else
                rho=edens(i1,i2,i3)
                d2=0.0_dp
                do i=1,3
                   !dlogeps(i,i1,i2,i3)=dtx*nabla_edens(i1,i2,i3,isp,i)
                   d2 = d2+nabla_edens(i1,i2,i3,i)**2
                end do
                d=sqrt(d2)
                dd = ddt_edens(i1,i2,i3)
                de=epsprime(rho,kernel%cavity)
                if (dumpeps) epscurr(i1,i2,i3)=eps(rho,kernel%cavity)
                depsdrho(i1,i23)=de
                kernel%w%oneoeps(i1,i23)=oneoeps(rho,kernel%cavity) !in pi
                !c1=(cc(i1,i2,i3)/d2-dd)/d
                dsurfdrho(i1,i23)=surf_term(rho,d2,dd,cc(i1,i2,i3),kernel%cavity)/epsm1

                !evaluate surfaces and volume integrals
                IntSur=IntSur + de*d/epsm1
                IntVol=IntVol + (kernel%cavity%epsilon0-eps(rho,kernel%cavity))/epsm1


!!$              if (dabs(edens(i1,i2,i3)).gt.edensmax) then
!!$                 !eps(i1,i2,i3)=1.d0
!!$                 if (dumpeps) epscurr(i1,i2,i3)=1.d0
!!$                 kernel%oneoeps(i1,i23)=1.d0 !oneoeps(i1,i2,i3)
!!$                 depsdrho(i1,i23)=0.d0
!!$                 dsurfdrho(i1,i23)=0.d0
!!$              else if (dabs(edens(i1,i2,i3)).lt.edensmin) then
!!$                 !eps(i1,i2,i3)=eps0
!!$                 if (dumpeps) epscurr(i1,i2,i3)=eps0
!!$                 kernel%oneoeps(i1,i23)=oneoeps0 !oneoeps(i1,i2,i3)
!!$                 depsdrho(i1,i23)=0.d0
!!$                 dsurfdrho(i1,i23)=0.d0
!!$              else
!!$                 r=fact1*(log(edensmax)-log(dabs(edens(i1,i2,i3))))
!!$                 t=fact2*(r-sin(r))
!!$                 coeff=fact3*(1.d0-cos(r))
!!$                 dtx=-coeff/dabs(edens(i1,i2,i3))
!!$                 de=safe_exp(t)*dtx ! derivative of epsilon wrt rho
!!$                 depsdrho(i1,i23)=de
!!$                 !eps(i1,i2,i3)=dexp(t)
!!$                 if (dumpeps) epscurr(i1,i2,i3)=safe_exp(t)
!!$                 kernel%oneoeps(i1,i23)=safe_exp(-t) !oneoeps(i1,i2,i3)
!!$                 d2=0.d0
!!$                 do i=1,3
!!$                    !dlogeps(i,i1,i2,i3)=dtx*nabla_edens(i1,i2,i3,isp,i)
!!$                    d2 = d2+nabla_edens(i1,i2,i3,i)**2
!!$                 end do
!!$                 d=dsqrt(d2)
!!$                 dd = ddt_edens(i1,i2,i3)
!!$                 c1=(cc(i1,i2,i3)/d2-dd)/d
!!$                 dsurfdrho(i1,i23)=(de*c1)/(eps0-1.d0)
!!$                 !dsurfdrho(i1,i23)=(de*c1+dde*c2)/(eps0-1.d0)
!!$                 IntSur=IntSur + de*d
!!$                 IntVol=IntVol + (eps0-safe_exp(t))/(eps0-1.d0)
!!$              end if

             end if
           end do
           i23=i23+1
        end do
     end do

     !then dlogeps
     do i3=1,n03
        do i2=1,n02
           do i1=1,n01
             if (epsinner(i1,i2,i3).gt.innervalue) then ! Check for inner sccs cavity value to fix as vacuum
              do i=1,3
               kernel%w%dlogeps(i,i1,i2,i3)=0.d0 !dlogeps(i,i1,i2,i3)
              end do
             else
                logepspr=logepsprime(edens(i1,i2,i3),kernel%cavity)
                do i=1,3
                   kernel%w%dlogeps(i,i1,i2,i3)=nabla_edens(i1,i2,i3,i)*logepspr
                end do

!!$              if (dabs(edens(i1,i2,i3)).gt.edensmax) then
!!$                 do i=1,3
!!$                    kernel%dlogeps(i,i1,i2,i3)=0.d0 !dlogeps(i,i1,i2,i3)
!!$                 end do
!!$              else if (dabs(edens(i1,i2,i3)).lt.edensmin) then
!!$                 do i=1,3
!!$                    kernel%dlogeps(i,i1,i2,i3)=0.d0 !dlogeps(i,i1,i2,i3)
!!$                 end do
!!$              else
!!$                 r=fact1*(log(edensmax)-log(dabs(edens(i1,i2,i3))))
!!$                 coeff=fact3*(1.d0-cos(r))
!!$                 dtx=-coeff/dabs(edens(i1,i2,i3))
!!$                 do i=1,3
!!$                    kernel%dlogeps(i,i1,i2,i3)=dtx*nabla_edens(i1,i2,i3,i) !dlogeps(i,i1,i2,i3)
!!$                 end do
!!$              end if

             end if
           end do
        end do
     end do

  end select

  IntSur=IntSur*kernel%hgrids(1)*kernel%hgrids(2)*kernel%hgrids(3)!/(eps0-1.d0)
  IntVol=IntVol*kernel%hgrids(1)*kernel%hgrids(2)*kernel%hgrids(3)

  Cavene= gammaSau*IntSur*627.509469d0
  Repene= alphaSau*IntSur*627.509469d0
  Disene=  betaVau*IntVol*627.509469d0
  noeleene=Cavene+Repene+Disene

  if (kernel%mpi_env%iproc==0 .and. kernel%mpi_env%igroup==0) then
     call yaml_map('Surface integral',IntSur)
     call yaml_map('Volume integral',IntVol)
     call yaml_map('Cavity energy',Cavene)
     call yaml_map('Repulsion energy',Repene)
     call yaml_map('Dispersion energy',Disene)
     call yaml_map('Total non-electrostatic energy',noeleene)
  end if

  if (dumpeps) then

     unt=f_get_free_unit(21)
     call f_open_file(unt,file='epsilon_sccs.dat')
     i1=1!n03/2
     do i2=1,n02
        do i3=1,n03
           write(unt,'(2(1x,I4),3(1x,e14.7))')i2,i3,epscurr(i1,i2,i3),epscurr(n01/2,i2,i3),edens(n01/2,i2,i3)
        end do
     end do
     call f_close(unt)

     unt=f_get_free_unit(22)
     call f_open_file(unt,file='epsilon_line_sccs_x.dat')
     do i1=1,n01
        x=i1*kernel%hgrids(1)
        write(unt,'(1x,I8,4(1x,e22.15))')i1,x,epscurr(i1,n02/2,n03/2),edens(i1,n02/2,n03/2),depsdrho1(i1,n02/2,n03/2)
     end do
     call f_close(unt)

     unt=f_get_free_unit(23)
     call f_open_file(unt,file='epsilon_line_sccs_y.dat')
     do i2=1,n02
        y=i2*kernel%hgrids(2)
        write(unt,'(1x,I8,3(1x,e22.15))')i2,y,epscurr(n01/2,i2,n03/2),edens(n01/2,i2,n03/2)
     end do
     call f_close(unt)

     unt=f_get_free_unit(24)
     call f_open_file(unt,file='epsilon_line_sccs_z.dat')
     do i3=1,n03
        z=i3*kernel%hgrids(3)
        write(unt,'(1x,I8,3(1x,e22.15))')i3,z,epscurr(n01/2,n02/2,i3),edens(n01/2,n02/2,i3)
     end do
     call f_close(unt)

     call f_free(epscurr)

  end if

  call f_free(ddt_edens)
  call f_free(nabla_edens)
  call f_free(epsinner)
  call f_free(depsdrho1)
  call f_free(cc)

end subroutine pkernel_build_epsilon


!> New version of the pkernel_build epsilon Routine,
!! with explicit workarrays.
!! This version is supposed not to allocate any array
subroutine rebuild_cavity_from_rho(rho_full,nabla_rho,nabla2_rho,delta_rho,cc_rho,depsdrho,dsurfdrho,&
     kernel,IntSur,IntVol)
  use FDder
  implicit none
  type(coulomb_operator), intent(inout) :: kernel
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(in) :: rho_full
  real(gp), intent(out) :: IntSur,IntVol
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(out) :: delta_rho
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(out) :: cc_rho,nabla2_rho
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),3), intent(out) :: nabla_rho
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(out) :: depsdrho,dsurfdrho
  !local variables
  integer :: n01,n02,n03,i3s

  n01=kernel%ndims(1)
  n02=kernel%ndims(2)
  n03=kernel%ndims(3)

  !calculate the derivatives of the density
  !build the gradients and the laplacian of the density
  !density gradient in du
  !call nabla_u(kernel%geocode,n01,n02,n03,rho_full,nabla_rho,kernel%nord,kernel%hgrids)
  call nabla_u_and_square(kernel%geocode,n01,n02,n03,rho_full,nabla_rho,nabla2_rho,&
       kernel%nord,kernel%hgrids)
  !density laplacian in delta_rho
  call div_u_i(kernel%geocode,n01,n02,n03,nabla_rho,delta_rho,kernel%nord,kernel%hgrids,cc_rho)

  !aliasing for the starting point
  i3s=kernel%grid%istart+1
  if (kernel%grid%n3p == 0 ) i3s=1
  call build_cavity_from_rho(rho_full(1,1,i3s),nabla2_rho(1,1,i3s),delta_rho(1,1,i3s),cc_rho(1,1,i3s),&
       kernel,depsdrho,dsurfdrho,IntSur,IntVol)

  if (kernel%method=='PI') then
     !form the inner cavity with a gathering
     call PS_gather(src=kernel%w%epsinnersccs,dest=delta_rho,kernel=kernel)
     call dlepsdrho_sccs(kernel%ndims,rho_full,nabla_rho,delta_rho,kernel%w%dlogeps,kernel%cavity)
  end if

end subroutine rebuild_cavity_from_rho

subroutine inplane_partitioning(mpi_env,mdz,n2wires,n3planes,part_mpi,inplane_mpi,n3pr1,n3pr2)
  use wrapper_mpi
  use yaml_output
  implicit none
  integer, intent(in) :: mdz              !< dimension of the density in the z direction
  integer, intent(in) :: n2wires,n3planes !<number of interesting wires in a plane and number of interesting planes
  type(mpi_environment), intent(in) :: mpi_env !< global env of Psolver
  integer, intent(out) :: n3pr1,n3pr2     !< mpi grid of processors based on mpi_env
  type(mpi_environment), intent(out) :: inplane_mpi,part_mpi !<internal environments for the partitioning
  !local variables

  !condition for activation of the inplane partitioning (inactive for the moment due to breaking of OMP parallelization)
  if (mpi_env%nproc>2*(n3planes)-1 .and. .false.) then
     n3pr1=mpi_env%nproc/(n3planes)
     n3pr2=n3planes
!!$     md2plus=.false.
!!$     if ((mdz/mpi_env%nproc)*n3pr1*n3pr2 < n2wires) then
!!$        md2plus=.true.
!!$     endif

     if (mpi_env%iproc==0 .and. n3pr1>1 .and. mpi_env%igroup==0 ) then
          call yaml_map('PSolver n1-n2 plane mpi partitioning activated:',&
          trim(yaml_toa(n3pr1,fmt='(i5)'))//' x'//trim(yaml_toa(n3pr2,fmt='(i5)'))//&
          ' taskgroups')
        call yaml_map('md2 was enlarged for PSolver n1-n2 plane mpi partitioning, md2=',mdz)
     end if

     if (n3pr1>1) &
          call mpi_environment_set1(inplane_mpi,mpi_env%iproc, &
          mpi_env%mpi_comm,n3pr1,n3pr2)
  else
     n3pr1=1
     n3pr2=mpi_env%nproc
     inplane_mpi=mpi_environment_null()
  endif

  call mpi_environment_set(part_mpi,mpi_env%iproc,&
       mpi_env%nproc,mpi_env%mpi_comm,n3pr2)

end subroutine inplane_partitioning
