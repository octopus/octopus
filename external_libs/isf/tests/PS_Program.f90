!> @file
!!  Program test for Poisson
!!  Laplacian V = 4pi rho
!!  May work either in parallel or in serial case
!!  And for different geometries
!! @author
!!    Copyright (C) 2006-2012 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test program for the Poisson Solver
program PSolver_Program
  use Poisson_Solver
  use wrapper_mpi
  use time_profiling
  use dynamic_memory
  implicit none
  !include 'mpif.h'
  !Order of interpolating scaling function
  !integer, parameter :: itype_scf=8
  character(len=*), parameter :: subname='Poisson_Solver'
  real(kind=8), parameter :: a_gauss = 1.0d-2, a2 = a_gauss**2
  !Length of the box
  real(kind=8), parameter :: acell = 10.0d0
  real(kind=8), parameter :: EulerGamma = 0.5772156649015328d0
  character(len=50) :: chain
  character(len=1) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  character(len=1) :: datacode
  character(len=30) :: mode
  real(kind=8), dimension(:,:,:), allocatable :: density,rhopot,potential,pot_ion
  type(coulomb_operator) :: karray
  real(kind=8) :: hx,hy,hz,max_diff,eh,exc,vxc,hgrid,diff_parser,offset,mu0
  real(kind=8) :: ehartree,eexcu,vexcu,diff_par,diff_ser,e1
  integer :: n01,n02,n03,itype_scf,i_all,i_stat
  integer :: i1,i2,i3,j1,j2,j3,i1_max,i2_max,i3_max,iproc,nproc,ierr,i3sd,ncomp
  integer :: n_cell,ixc,n3d,n3p,n3pi,i3xcsh,i3s
  logical :: alsoserial,onlykernel
  integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
  !triclinic lattice
  real(kind=8) :: alpha,beta,gamma,detg
  real(kind=8), dimension(:,:,:,:), pointer :: rhocore_fake
  external :: gather_timings  
  nullify(rhocore_fake)

  alpha = 2.0_dp*datan(1.0_dp)
  beta  = 2.0_dp*datan(1.0_dp)
  gamma = 2.0_dp*datan(1.0_dp)
  !alpha = 1.0_dp*datan(1.0_dp)
  !beta  = 1.0_dp*datan(1.0_dp)
  !gamma = 1.0_dp*datan(1.0_dp)

  detg = 1.0_dp - dcos(alpha)**2 - dcos(beta)**2 - dcos(gamma)**2 + 2.0_dp*dcos(alpha)*dcos(beta)*dcos(gamma)


  !mode = "charged_thin_wire"
  !mode="cylindrical_capacitor"
  !mode="monopolar"
  mode="zigzag_model_wire"
    

  call f_lib_initialize()

  !Use arguments
  call get_command_argument(1,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Program n01 n02 n03 ixc geocode datacode itype_scf [mu0_screening]'
     stop
  end if
  read(unit=chain,fmt=*) n01

  call get_command_argument(2,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Program n01 n02 n03 ixc geocode datacode itype_scf [mu0_screening]'
     stop
  end if
  read(unit=chain,fmt=*) n02

  call get_command_argument(3,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Program n01 n02 n03 ixc geocode datacode itype_scf [mu0_screening]'
     stop
  end if
  read(unit=chain,fmt=*) n03

  call get_command_argument(4,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Program n01 n02 n03 ixc geocode datacode itype_scf [mu0_screening]'
     stop
  end if
  read(unit=chain,fmt=*) ixc

  call get_command_argument(5,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Program n01 n02 n03 ixc geocode datacode itype_scf [mu0_screening]'
     stop
  end if
  read(unit=chain,fmt=*) geocode

  call get_command_argument(6,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Program n01 n02 n03 ixc geocode datacode itype_scf [mu0_screening]'
     stop
  end if
  read(unit=chain,fmt=*) datacode

  call get_command_argument(7,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Program n01 n02 n03 ixc geocode datacode itype_scf [mu0_screening]'
     stop
  end if
  read(unit=chain,fmt=*) itype_scf

  call get_command_argument(8,chain)
  if(trim(chain)=='') then
     mu0 = 0.0_dp
  else
     read(unit=chain,fmt=*) mu0
  end if


  !write(*,*) 'mu0 =', mu0

  !perform also the comparison with the serial case
  alsoserial=.false.
  onlykernel=.false.
  !code for the Poisson Solver in the parallel case
  !datacode='G'

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  if (geocode == 'P') then

     if (iproc==0) print *,"PSolver, periodic BC: ",n01,n02,n03,'processes',nproc
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
  
  else if (geocode == 'S') then

     if (iproc==0) print *,"PSolver for surfaces: ",n01,n02,n03,'processes',nproc
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  
  else if (geocode == 'F') then

     if (iproc==0) print *,"PSolver, free BC: ",n01,n02,n03,'processes',nproc
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  
  else if (geocode == 'W') then
    
     if (iproc==0) print *,"PSolver, wires BC: ",n01,n02,n03,'processes',nproc
     call W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
     
  else if (geocode == 'H') then
   
     if (iproc==0) print *,"PSolver, Helmholtz Equation Solver: ",n01,n02,n03,'processes',nproc
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  
  end if

  !write(*,*) n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3

  !initialize memory counting
  !call memocc(0,iproc,'count','start')

  !Step size
  n_cell = max(n01,n02,n03)
  hx=acell/real(n01,kind=8)
  hy=acell/real(n02,kind=8)
  hz=acell/real(n03,kind=8)

  write(*,'(a5,1pe20.12)') ' hx = ', hx

  !grid for the free BC case
  hgrid=max(hx,hy,hz)
  !hgrid=hx

  !we must choose properly a test case with a positive density
  !itype_scf=16

  write(*,'(a12,i4)') ' itype_scf = ', itype_scf 

  call f_timing_reset(filename='time.yaml',master=iproc==0)
  !call timing(nproc,'time.prc','IN')

  karray=pkernel_init(.true.,iproc,nproc,0,&
       geocode,(/n01,n02,n03/),(/hx,hy,hz/),itype_scf,mu0,(/alpha,beta,gamma/))
  call pkernel_set(karray,.true.)

  !call createKernel(iproc,nproc,geocode,(/n01,n02,n03/),(/hx,hy,hz/),itype_scf,karray,.true.,mu0,(/alpha,beta,gamma/))
  !print *,'sum',sum(karray%kernel)
  if (.not. onlykernel) then
     !Allocations
     !Density
     density = f_malloc((/ n01, n02, n03 /),id='density')
     !Density then potential
     rhopot = f_malloc((/ n01, n02, n03 /),id='rhopot')
     potential = f_malloc((/ n01, n02, n03 /),id='potential')
     !ionic potential
     pot_ion = f_malloc((/ n01, n02, n03 /),id='pot_ion')

     call test_functions(geocode,ixc,n01,n02,n03,acell,a_gauss,hx,hy,hz,&
          density,potential,rhopot,pot_ion,mu0,alpha,beta,gamma)

     ! i2=n02/2
     ! do i3=1,n03
     !    do i1=1,n01
     !       j1=n01/2+1-abs(n01/2+1-i1)
     !       j2=n02/2+1-abs(n02/2+1-i2)
     !       j3=n03/2+1-abs(n03/2+1-i3)
     !       write(110,*)i1,i3,rhopot(i1,i2,i3),potential(i1,i2,i3)               
     !    end do
     ! end do


     !offset, used only for the periodic solver case
     if (ixc==0) then
        offset=0.0_gp!potential(1,1,1)!-pot_ion(1,1,1)
        do i3=1,n03
           do i2=1,n02
              do i1=1,n01
                 offset=offset+potential(i1,i2,i3)
              end do
           end do
        end do
        offset=offset*hx*hy*hz*sqrt(detg) ! /// to be fixed ///
        write(*,*) 'offset = ',offset
     end if

     !dimension needed for allocations
     call PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,(ixc>10),.false.,n3d,n3p,n3pi,i3xcsh,i3s)

     !dimension for comparison in the global or distributed poisson solver
     if (datacode == 'G') then
        i3sd=1
        ncomp=n03
     else if (datacode == 'D') then
        i3sd=i3s
        ncomp=n3p
     end if

!!  print *,'iproc,i3xcsh,i3s',iproc,i3xcsh,i3s

!!$     print *,'density',density(25,25,25)
!!$     print *,'potential',potential(25,25,25)
!!$     print *,'rhopot',rhopot(25,25,25)
!!$     print *,'pot_ion',pot_ion(25,25,25)


     !apply the Poisson Solver (case with distributed potential)
     eexcu=0.0_gp
     vexcu=0.0_gp
     call H_potential(datacode,karray,density(1,1,i3sd),pot_ion(1,1,i3s+i3xcsh),ehartree,offset,.false.)
!!$     call PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
!!$          density(1,1,i3sd),karray%kernel,pot_ion(1,1,i3s+i3xcsh),ehartree,eexcu,vexcu,offset,.true.,1,alpha,beta,gamma)

     print *,'potential integral',sum(density)

     i3=n03/2
     do i2=1,n02
        do i1=1,n01
           !j1=n01/2+1-abs(n01/2+1-i1)
           !j2=n02/2+1-abs(n02/2+1-i2)
           !j3=n03/2+1-abs(n03/2+1-i3)
           write(111,*)i1*hx,i2*hy,rhopot(i1,i2,i3),potential(i1,i2,i3),density(i1,i2,i3)
           !write(111,*) i1*hx+hy*i2*dcos(alpha)+i3*hz*dcos(beta), &
           !     i2*hy*dsin(alpha)+i3*hz*(-dcos(alpha)*dcos(beta)+dcos(gamma))/dsin(alpha), &
           !     rhopot(i1,i2,i3),potential(i1,i2,i3), &
           !     density(i1,i2,i3)
        end do
     end do

     i2=n02/2
     do i3=1,n03
        !do i2=1,n02
           do i1=1,n01
              !j1=n01/2+1-abs(n01/2+1-i1)
              !j2=n02/2+1-abs(n02/2+1-i2)
              !j3=n03/2+1-abs(n03/2+1-i3)
              write(112,*)i1*hx,i2*hy,i3*hz,rhopot(i1,i2,i3),potential(i1,i2,i3),&
                   density(i1,i2,i3)
           end do
        !end do
     end do
!!$     print *,'density2',density(25,25,25)
!!$     print *,'potential2',potential(25,25,25)
!!$     print *,'rhopot2',rhopot(25,25,25)
!!$     print *,'pot_ion2',pot_ion(25,25,25)

  end if

  
!!$  if (geocode == 'P') then
!!$     open(unit=65,file='karrayP.dump')
!!$  elseif (geocode == 'S') then
!!$     open(unit=65,file='karrayS.dump')
!!$  elseif (geocode == 'F') then
!!$     open(unit=65,file='karrayF.dump')
!!$  elseif (geocode == 'W') then
!!$     open(unit=65,file='karrayW.dump')
!!$  end if
!!$
!!$ 
!!$  do i1 = 1, nd1-1
!!$     do i3 = 1, nd3-1
!!$        write(65,fmt="(3(1pe20.12e3))") i1*1.0_dp, i3*1.0_dp, karray%kernel(i1 + (nd2-1)*nd1 + i3*nd1*nd2)
!!$     end do
!!$  end do
!!$  close(65)
  call pkernel_free(karray)
!!$  i_all=-product(shape(karray))*kind(karray)
!!$  deallocate(karray,stat=i_stat)
!!$  call memocc(i_stat,i_all,'karray',subname)

  call f_timing_stop(mpi_comm=karray%mpi_env%mpi_comm,nproc=karray%mpi_env%nproc,gather_routine=gather_timings)

  if (.not. onlykernel) then

     !comparison (each process compare its own part)
     call compare(n01,n02,ncomp,potential(1,1,i3sd+i3xcsh),density(1,1,i3sd+i3xcsh),&
          i1_max,i2_max,i3_max,max_diff)

!!  print *,'iproc,i3xcsh,i3s,max_diff',iproc,i3xcsh,i3s,max_diff

     !extract the max
     if (nproc > 1) then
        call MPI_ALLREDUCE(max_diff,diff_par,1,MPI_double_precision,  &
             MPI_MAX,MPI_COMM_WORLD,ierr)
     else
        diff_par=max_diff
     end if


     if (iproc == 0) then

        write(*,*) '--------------------'
        write(*,*) 'Parallel calculation '
        write(unit=*,fmt="(1x,a,3(1pe20.12))") "eht, exc, vxc:",ehartree,eexcu,vexcu
        write(*,'(a,3(i0,1x))') '  Max diff at: ',i1_max,i2_max,i3_max
        write(unit=*,fmt="(1x,a,1pe20.12)") '    Max diff:',diff_par,&
             '      result:',density(i1_max,i2_max,i3_max),&
             '    original:',potential(i1_max,i2_max,i3_max)
     end if

  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !Serial case
  if (alsoserial) then
     call f_timing_reset(filename='time_serial.yaml',master=iproc==0)
     !call timing(0,'             ','IN')

     karray=pkernel_init(.true.,0,1,0,&
          geocode,(/n01,n02,n03/),(/hx,hy,hz/),itype_scf,mu0,(/alpha,beta,gamma/))

     call pkernel_set(karray,.true.)

!!$     call createKernel(0,1,geocode,(/n01,n02,n03/),(/hx,hy,hz/),itype_scf,karray,.true.,mu0,&
!!$          (/alpha,beta,gamma/))

     if (.not. onlykernel) then
        !offset, used only for the periodic solver case
        if(ixc==0) offset=potential(1,1,1)!-pot_ion(1,1,1)
        
        
!!$        !this is how it should be called. Temporarily desactivated
!!$        call XC_potential(geocode,'G',karray%mpi_env%iproc,karray%mpi_env%nproc,&
!!$             karray%mpi_env%mpi_comm,n01,n02,n03,ixc,hx,hy,hz,&
!!$             rhopot,exc,vxc,1,rhocore_fake,V_XC,xcstr)
!!$        
!!$        call H_potential('G',karray,rhopot,pot_ion,eh,0.0_dp,.true.)
!!$
!!$        !apply the Poisson Solver (case with distributed potential
!!$        !this is the old call
!!$        call PSolver(geocode,'G',0,1,n01,n02,n03,ixc,hx,hy,hz,&
!!$             rhopot,karray%kernel,pot_ion,eh,exc,vxc,offset,.true.,1,alpha,beta,gamma)
        
     end if
     call pkernel_free(karray)
!!$     i_all=-product(shape(karray))*kind(karray)
!!$     deallocate(karray,stat=i_stat)
!!$     call memocc(i_stat,i_all,'karray',subname)

     !call timing(,'              ','RE')
     call f_timing_stop(mpi_comm=karray%mpi_env%mpi_comm)

     if (.not. onlykernel) then
        !Maximum difference
        call compare(n01,n02,n03,potential,rhopot,i1_max,i2_max,i3_max,diff_ser)

        !! print *,'iproc,diff_ser',iproc,diff_ser

        if (iproc==0) then
           write(*,*) '------------------'
           write(*,*) 'Serial Calculation'
           write(*,"(1x,a,3(1pe20.12))") "eht, exc, vxc:",eh,exc,vxc
           write(*,'(a,3(i0,1x))') '  Max diff at: ',i1_max,i2_max,i3_max
           write(*,"(1x,a,1pe20.12)") '     Max diff:',diff_ser,&
                '       result:',rhopot(i1_max,i2_max,i3_max),&
                '     original:',potential(i1_max,i2_max,i3_max)
        end if

        !Maximum difference, parallel-serial
        call compare(n01,n02,ncomp,rhopot(1,1,i3sd+i3xcsh),density(1,1,i3sd+i3xcsh),&
             i1_max,i2_max,i3_max,max_diff)

!!     print *,'max_diff,i1_max,i2_max,i3_max,i3s,i3xcsh,n3p',max_diff,i1_max,i2_max,i3_max,&
!!          i3s,i3xcsh,n3p

        if (nproc > 1) then
           !extract the max
           call MPI_ALLREDUCE(max_diff,diff_parser,1,MPI_double_precision,  &
                MPI_MAX,MPI_COMM_WORLD,ierr)
        else
           diff_parser=max_diff
        end if

        if (iproc==0) then
           write(*,*) '------------------'
           write(*,'(a,3(i0,1x))')&
                'difference parallel-serial, at',i1_max,i2_max,i3_max
           write(*,"(1x,a,1pe12.4)")&
                '    Max diff:',diff_parser,&
                '    parallel:',density(i1_max,i2_max,i3_max),&
                '      serial:',rhopot(i1_max,i2_max,i3_max)
           write(*,"(1x,a,3(1pe12.4))")&
                "energy_diffs:",ehartree-eh,eexcu-exc,vexcu-vxc
        end if
     end if
  end if
  if (iproc==0 .and. .not. onlykernel) then

     call regroup_data(geocode,hx,hy,hz)

     i2=i2_max
     do i3=1,n03
        do i1=1,n01
           j1=n01/2+1-abs(n01/2+1-i1)
           j2=n02/2+1-abs(n02/2+1-i2)
           j3=n03/2+1-abs(n03/2+1-i3)
           write(11,*)i1,i3,rhopot(i1,i2,i3),potential(i1,i2,i3),&
                density(i1,i2,i3)
        end do
     end do
     i3=i3_max
     do i2=1,n02
        do i1=1,n01
           j1=n01/2+1-abs(n01/2+1-i1)
           j2=n02/2+1-abs(n02/2+1-i2)
           j3=n03/2+1-abs(n03/2+1-i3)
           write(12,*)i1,i2,rhopot(i1,i2,i3),potential(i1,i2,i3),&
                density(i1,i2,i3)
        end do
     end do

  end if

  if (.not. onlykernel) then
     call f_free(density)
     call f_free(rhopot)
     call f_free(potential)
     call f_free(pot_ion)
  end if

  call f_lib_finalize()

  call MPI_FINALIZE(ierr)


contains


subroutine regroup_data(geocode,hx,hy,hz)
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  real(kind=8), intent(in) :: hx,hy,hz
  !local variables
  real(kind=8) :: hgrid
  hgrid=max(hx,hy,hz)
  if (geocode == 'S') hgrid=hy

  !open(unit=60,file='time.par',status='unknown')
  !read(60,*)
  !read(60,*)string,tcp1,tcp2,pcp
  !read(60,*)string,tcm1,tcm2,pcm
  !read(60,*)string,tk1,tk2,pk
  !read(60,*)string,txc1,txc2,pxc
  !close(60)
  !tcp=tcp2
  !tcm=tcm2
  !tk=tk2
  !txc=txc2
  
!  write(99,'(a2,3(i4),1pe9.2,1pe10.3,4(1pe9.2),4(0pf5.1),1pe9.2)')&
 !      geocode,n01,n02,n03,hgrid,max_diff,tcp,tcm,tk,txc,pcp,pcm,pk,pxc,diff_parser
  
end subroutine regroup_data


subroutine compare(n01,n02,n03,potential,density,i1_max,i2_max,i3_max,max_diff)
  implicit none
  integer, intent(in) :: n01,n02,n03
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,density
  integer, intent(out) :: i1_max,i2_max,i3_max
  real(kind=8), intent(out) :: max_diff

  !local variables
  integer :: i1,i2,i3
  real(kind=8) :: factor
  max_diff = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  do i3=1,n03
     do i2=1,n02 
        do i1=1,n01
           factor=abs(potential(i1,i2,i3)-density(i1,i2,i3))
           if (max_diff < factor) then
              max_diff = factor
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do
end subroutine compare


!> This subroutine builds some analytic functions that can be used for 
!! testing the poisson solver.
!! The default choice is already well-tuned for comparison.
!! WARNING: not all the test functions can be used for all the boundary conditions of
!! the poisson solver, in order to have a reliable analytic comparison.
!! The parameters of the functions must be adjusted in order to have a sufficiently localized
!! function in the isolated direction and an explicitly periodic function in the periodic ones.
!! Beware of the high-frequency components that may falsify the results when hgrid is too high.
subroutine test_functions(geocode,ixc,n01,n02,n03,acell,a_gauss,hx,hy,hz,&
     density,potential,rhopot,pot_ion,mu0,alpha,beta,gamma)
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  integer, intent(in) :: n01,n02,n03,ixc
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz,mu0
  !triclinic lattice
  real(kind=8), intent(in) :: alpha, beta, gamma
  real(kind=8), dimension(n01,n02,n03), intent(out) :: density,potential,rhopot,pot_ion

  !local variables
  integer :: i1,i2,i3,ifx,ify,ifz
  real(kind=8) :: x,x1,x2,x3,y,z,length,denval,pi,a2,derf,factor,r,r2,r0
  real(kind=8) :: fx,fx1,fx2,fy,fy1,fy2,fz,fz1,fz2,a,ax,ay,az,bx,by,bz,tt,potion_fac
  real(kind=8) :: monopole
  real(kind=8), dimension(3) :: dipole
 
  !non-orthorhombic lattice
  real(kind=8), dimension(3,3) :: gu,gd
  real(kind=8) :: detg


  !triclinic cell
  !covariant metric
  gd(1,1) = 1.0_dp
  gd(1,2) = dcos(alpha)
  gd(1,3) = dcos(beta)
  gd(2,2) = 1.0_dp
  gd(2,3) = dcos(gamma)
  gd(3,3) = 1.0_dp

  gd(2,1) = gd(1,2)
  gd(3,1) = gd(1,3)
  gd(3,2) = gd(2,3)
  !
  detg = 1.0_dp - dcos(alpha)**2 - dcos(beta)**2 - dcos(gamma)**2 + 2.0_dp*dcos(alpha)*dcos(beta)*dcos(gamma)

  write(*,*) 'detg =', detg
  !
  !contravariant metric
  gu(1,1) = (dsin(gamma)**2)/detg
  gu(1,2) = (dcos(beta)*dcos(gamma)-dcos(alpha))/detg
  gu(1,3) = (dcos(alpha)*dcos(gamma)-dcos(beta))/detg
  gu(2,2) = (dsin(beta)**2)/detg
  gu(2,3) = (dcos(alpha)*dcos(beta)-dcos(gamma))/detg
  gu(3,3) = (dsin(alpha)**2)/detg
  !
  gu(2,1) = gu(1,2)
  gu(3,1) = gu(1,3)
  gu(3,2) = gu(2,3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ixc==0) denval=0.d0

  if (trim(geocode) == 'P') then
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     !test functions in the three directions
     ifx=6
     ify=6
     ifz=6
     !parameters of the test functions
     ax=length
     ay=length
     az=length
     
     !the following b's are not used, actually
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0
     
   
     ! !original version

     !plot of the functions used
     do i1=1,n03
        x = hx*real(i1-n01/2-1,kind=8)!valid if hy=hz
        y = hz*real(i1-n03/2-1,kind=8) 
        call functions(x,ax,bx,fx,fx1,fx2,ifx)
        call functions(y,az,bz,fz,fz1,fz2,ifz)
        write(20,*)i1,fx,fx2,fz,fz2
     end do

     !Initialization of density and potential
     denval=0.d0 !value for keeping the density positive
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz1,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy1,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx1,fx2,ifx)
              potential(i1,i2,i3) =  -16.d0*datan(1.d0)*fx*fy*fz
              !density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2-mu0**2*fx*fy*fz
              !triclinic lattice
              density(i1,i2,i3) = -mu0**2*fx*fy*fz
              density(i1,i2,i3) = density(i1,i2,i3) + gu(1,1)*fx2*fy*fz+gu(2,2)*fx*fy2*fz+gu(3,3)*fx*fy*fz2
              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(gu(1,2)*fx1*fy1*fz+gu(1,3)*fx1*fy*fz1+gu(2,3)*fx*fy1*fz1)
              denval=max(denval,-density(i1,i2,i3))
           end do
        end do
     end do


     ! !tweaked version: for debugging the solver for non-orthorhombic cells

     ! pi = 4.d0*atan(1.d0)
     ! a2 = a_gauss**2/2
     ! !mu0 = 1.e0_dp

     ! !Normalization
     ! !factor = a_gauss*sqrt(pi)/2.0_dp
     ! factor = 2.0_dp
     ! !gaussian function
     ! do i3=1,n03
     !    !x3 = hz*real(i3-n03/2,kind=8)
     !    do i2=1,n02
     !       !x2 = hy*real(i2-n02/2,kind=8)
     !       do i1=1,n01
     !          x1 = hx*real(i1-n01/2,kind=8)+hy*real(i2-n02/2,kind=8)*dcos(alpha)+hz*real(i3-n03/2,kind=8)*dcos(beta)
     !          x2 = hy*real(i2-n02/2,kind=8)*dsin(alpha) + & 
     !               & hz*real(i3-n03/2,kind=8)*(-dcos(alpha)*dcos(beta)+dcos(gamma))/dsin(alpha)
     !          x3 = hz*real(i3-n03/2,kind=8)*sqrt(detg)/dsin(alpha)
     !          !r2 = x1*x1+x2*x2+x3*x3
     !          !triclinic lattice:
     !          !r2 = gd(1,1)*x1*x1+gd(2,2)*x2*x2+gd(3,3)*x3*x3+2.0_dp*(gd(1,2)*x1*x2+gd(1,3)*x1*x3+gd(2,3)*x2*x3)
     !          r2 = x1*x1+x2*x2+x3*x3
     !          !density(i1,i2,i3) = factor*exp(-r2/a2)
     !          r = sqrt(r2)
     !          !Potential from a gaussian
     !          potential(i1,i2,i3) = dexp(-r2/a2)
     !          density(i1,i2,i3) = 4.0_dp*r2/a2**2-6.0_dp/a2
     !          density(i1,i2,i3) = -potential(i1,i2,i3)*density(i1,i2,i3)/16.0_dp/datan(1.0_dp)
     !       end do
     !    end do
     ! end do


     
     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           !j1=n01/2+1-abs(n01/2+1-i1)
           !j2=n02/2+1-abs(n02/2+1-i2)
           !j3=n03/2+1-abs(n03/2+1-i3)
           write(200,*) i1,i3,density(i1,i2,i3),potential(i1,i2,i3)               
        end do
     end do

!plane capacitor oriented along the y direction
!!     do i2=1,n02
!!        if (i2==n02/4) then
!!           do i3=1,n03
!!              do i1=1,n01
!!                 density(i1,i2,i3)=1.d0!real(i2,kind=8)
!!              end do
!!           end do
!!        else if (i2==3*n02/4) then
!!           do i3=1,n03
!!              do i1=1,n01
!!                 density(i1,i2,i3)=-1.d0!real(i2,kind=8)
!!              end do
!!           end do
!!        else
!!           do i3=1,n03
!!              do i1=1,n01
!!                 density(i1,i2,i3)=0.d0
!!              end do
!!           end do
!!        end if
!!     end do
!!     denval=0.d0

     if (ixc==0) denval=0.d0



  else if (trim(geocode) == 'S') then
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     !test functions in the three directions
     ifx=5
     ifz=1
     !non-periodic dimension
     ify=6
     !parameters of the test functions
     ax=length
     ay=length
     az=length
     !b's are not used, actually
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0!real(nu,kind=8)

     
     
     !non-periodic dimension
     !ay=length!4.d0*a

     density(:,:,:) = 0.d0!1d-20 !added

     !plot of the functions used
     do i1=1,n02
        x = hx*real(i1-n02/2-1,kind=8)!valid if hy=hz
        y = hy*real(i1-n02/2-1,kind=8) 
        call functions(x,ax,bx,fx,fx1,fx2,ifx)
        call functions(y,ay,by,fy,fy1,fy2,ify)
        write(20,*)i1,fx,fx2,fy,fy2
     end do

     !Initialisation of density and potential
     !Normalisation
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz1,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fz1,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx1,fx2,ifx)
              potential(i1,i2,i3) =  -16.d0*datan(1.d0)*fx*fy*fz
              density(i1,i2,i3) = -mu0**2*fx*fy*fz
              density(i1,i2,i3) = density(i1,i2,i3) + gu(1,1)*fx2*fy*fz+gu(2,2)*fx*fy2*fz+gu(3,3)*fx*fy*fz2
              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(gu(1,2)*fx1*fy1*fz+gu(1,3)*fx1*fy*fz1+gu(2,3)*fx*fy1*fz1)
              !old:
              !density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2 - mu0**2*fx*fy*fz
              denval=max(denval,-density(i1,i2,i3))
           end do
        end do
     end do

     
     ! !plane capacitor oriented along the y direction
     ! do i2=1,n02
     !    if (i2==n02/4) then
     !       do i3=1,n03
     !          do i1=1,n01
     !             density(i1,i2,i3)=1.d0!real(i2,kind=8)
     !          end do
     !       end do
     !    else if (i2==3*n02/4) then
     !       do i3=1,n03
     !          do i1=1,n01
     !             density(i1,i2,i3)=-1.d0!real(i2,kind=8)
     !          end do
     !       end do
     !    else
     !       do i3=1,n03
     !          do i1=1,n01
     !             density(i1,i2,i3)=0.d0
     !          end do
     !       end do
     !    end if
     ! end do



     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           !j1=n01/2+1-abs(n01/2+1-i1)
           !j2=n02/2+1-abs(n02/2+1-i2)
           !j3=n03/2+1-abs(n03/2+1-i3)
           write(200,*) i1,i3,density(i1,i2,i3),potential(i1,i2,i3)               
        end do
     end do


     if (ixc==0) denval=0.d0

  
!   else if (trim(geocode) == 'F') then

!      !grid for the free BC case
!      !hgrid=max(hx,hy,hz)

!      pi = 4.d0*atan(1.d0)
!      a2 = a_gauss**2

!      !Normalization
!      factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
!      !gaussian function
!      do i3=1,n03
!         x3 = hz*real(i3-n03/2,kind=8)
!         do i2=1,n02
!            x2 = hy*real(i2-n02/2,kind=8)
!            do i1=1,n01
!               x1 = hx*real(i1-n01/2,kind=8)
!               r2 = x1*x1+x2*x2+x3*x3
!               density(i1,i2,i3) = factor*exp(-r2/a2)
!               r = sqrt(r2)
!               !Potential from a gaussian
!               if (r == 0.d0) then
!                  potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
!               else
!                  potential(i1,i2,i3) = derf(r/a_gauss)/r
!               end if
!            end do
!         end do
!      end do

!      i2=n02/2
!      do i3=1,n03
!         do i1=1,n01
!            !j1=n01/2+1-abs(n01/2+1-i1)
!            !j2=n02/2+1-abs(n02/2+1-i2)
!            !j3=n03/2+1-abs(n03/2+1-i3)
!            write(200,*) i1,i3,density(i1,i2,i3),potential(i1,i2,i3)               
!         end do
!      end do

   
! !plane capacitor oriented along the y direction
! !!     do i2=1,n02
! !!        if (i2==n02/4) then
! !!           do i3=1,n03
! !!              do i1=1,n01
! !!                 density(i1,i2,i3)=1.d0!real(i2,kind=8)
! !!              end do
! !!           end do
! !!        else if (i2==3*n02/4) then
! !!           do i3=1,n03
! !!              do i1=1,n01
! !!                 density(i1,i2,i3)=-1.d0!real(i2,kind=8)
! !!              end do
! !!           end do
! !!        else
! !!           do i3=1,n03
! !!              do i1=1,n01
! !!                 density(i1,i2,i3)=0.d0
! !!              end do
! !!           end do
! !!        end if
! !!     end do
     
!      denval=0.d0


  else if (trim(geocode) == 'H' .or. trim(geocode) == 'F') then

     !hgrid=max(hx,hy,hz)

     pi = 4.d0*atan(1.d0)
     a2 = a_gauss**2
     !mu0 = 1.e0_dp

     !Normalization
     !factor = a_gauss*sqrt(pi)/2.0_dp
     factor = 2.0_dp-a_gauss*dexp(a2*mu0**2/4.0_dp)*sqrt(pi)*mu0*derfc(mu0*a_gauss/2.0_dp)
     !gaussian function
     do i3=1,n03
        x3 = hz*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hy*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hx*real(i1-n01/2,kind=8)
              !r2 = x1*x1+x2*x2+x3*x3
              !triclinic lattice:
              r2 = x1*x1+x2*x2+x3*x3
              !density(i1,i2,i3) = factor*exp(-r2/a2)
              r = sqrt(r2)
              !Potential from a gaussian
              if (r == 0.d0) then
                 potential(i1,i2,i3) = 1.0_dp
              else
                 potential(i1,i2,i3) = -2+derfc(r/a_gauss-a_gauss*mu0/2.0_dp)
                 potential(i1,i2,i3) = potential(i1,i2,i3)+dexp(2.0_dp*r*mu0)*derfc(r/a_gauss+a_gauss*mu0/2.0_dp)
                 potential(i1,i2,i3) = potential(i1,i2,i3)*a_gauss*dexp(-mu0*r)*sqrt(pi)/(-2*r*factor*dexp(-a2*mu0**2/4.0_dp))
              end if
              !density(i1,i2,i3) = exp(-r2/a2)/4.0_dp/factor**2 + 0.1_dp**2/(4*pi)*potential(i1,i2,i3)
              density(i1,i2,i3) = dexp(-r2/a2)/factor/(a2*pi)
           end do
        end do
     end do

     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           !j1=n01/2+1-abs(n01/2+1-i1)
           !j2=n02/2+1-abs(n02/2+1-i2)
           !j3=n03/2+1-abs(n03/2+1-i3)
           write(200,*) i1,i3,density(i1,i2,i3),potential(i1,i2,i3)               
        end do
     end do

        
     denval=0.d0


  else if (trim(geocode) == 'W') then
     !parameters for the test functions
     length=acell
     !a=0.5d0/a_gauss**2
     !test functions in the three directions
     !isolated directions
     ifx=6
     ify=6
     !periodic direction
     ifz=5
     !parameters of the test functions
     
     ax = length
     ay = length
     az = length
     
     bx = 2.d0
     by = 2.d0
     bz = 2.d0
  

     density(:,:,:) = 0.d0!1d-20 !added
     factor = 2.0d0


     !plot of the functions used
     do i1=1,min(n01,n03)
        x = hx*real(i1-n01/2-1,kind=8)!isolated
        z = hz*real(i1-n03/2-1,kind=8)!periodic
        call functions(x,ax,bx,fx,fx1,fx2,ifx)
        call functions(z,az,bz,fz,fz1,fz2,ifz)
        write(20,*) i1,fx,fx2,fz,fz2
     end do


     !Initialization of density and potential
   

     select case(mode)
     case ("monopolar")
        ! Gaussian Density Distribution in (x,y)
        ! this is the configuration yielding non-zero monopole
        do i3=1,n03
           x3 = hz*real(i3-n03/2-1,kind=8)
           call functions(x3,az,bz,fz,fz1,fz2,1)
           do i2=1,n02
              x2 = hy*real(i2-n02/2-1,kind=8)
              !call functions(x2,ay,by,fy,fy1,fy2,ify)
              do i1=1,n01
                 x1 = hx*real(i1-n01/2-1,kind=8)
                 r2 = x1*x1+x2*x2
                 r = sqrt(r2)
                 !call functions(x1,ax,bx,fx,fx1,fx2,ifx)
                 if  (r == 0.d0) then
                    !EulerGamma = 0.5772156649015328d
                    density(i1,i2,i3) = dexp(-factor*r2)
                    potential(i1,i2,i3) = (-EulerGamma - dlog(factor))/(4.0d0*factor)
                 else
                    call e1xb(factor*r2,e1)
                    density(i1,i2,i3) = dexp(-factor*r2)
                    potential(i1,i2,i3) = (e1+dlog(r2))/(4.0d0*factor)
                 end if
                 !note that in this case we cannot account for the screening in the following way,
                 !density(i1,i2,i3) = density(i1,i2,i3) + mu0**2*potential(i1,i2,i3)
                 !because the extra-term, proportional to the potential, is not localized
                 !in the non-periodic directions
                 potential(i1,i2,i3) = -16.0d0*datan(1.0d0)*potential(i1,i2,i3)
              end do
           end do
        end do

     case("zigzag_model_wire")

        density = 0.d0
        potential = 0.d0

!!$        density(n01/4,n02/2,n03/4) = -1.0d0
!!$        density(3*n01/4,n02/2,3*n03/4) = 1.0d0
        density(n01/2,n02/2,n03/2) = 1.0d0

        !factor=(16.d0/acell)**2
        !r0 = acell/4.d0
        !the following is evaluated analytically by imposing that
        !\int_0^\infty r*(-exp(-factor(r-r0)^2)+denval*exp(-factor*r^2)) = 0
        !denval=sqrt(4.d0*datan(1.d0)*factor)*r0*(1.d0+derf(sqrt(factor)*r0))
        !do i3=1,n03
        ! x3 = hz*real(i3-n03/2-1,kind=8)
        ! do i2=1,n02
        ! x2 = hy*real(i2-n02/2-1,kind=8)
        ! do i1=1,n01
        ! x1 = hx*real(i1-n01/2-1,kind=8)
        ! !r2 = x1*x1+x2*x2+x3*x3
        ! !r = sqrt(r2)
        ! !in this configuration denval is used so as to achieve zero monopole
        ! density(i1,i2,i3) = -1.d0*dexp(-factor*(x1-r0)**2)*dexp(-factor*x2**2)*dexp(-factor*(x3-r0)**2) &
        ! + 1.0d0*dexp(-factor*(x1+r0)**2)*dexp(-factor*x2**2)*dexp(-factor*(x3+r0)**2)
        ! density(i1,i2,i3) = density(i1,i2,i3)*(factor/4.d0/datan(1.d0))**(3.d0/2.d0)
        ! end do
        ! end do
        !end do

     case ("charged_thin_wire")
        do i3=1,n03
           do i2=1,n02
              do i1=1,n01
                 if (i1 == n01/2+1 .and. i2 == n02/2+1) density(i1,i2,i3) = 1.0d0
              end do
           end do
        end do
     case("cylindrical_capacitor")
        !mimicked by two Gaussian charge distributions,
        !one localized around r = 0,
        !the other around r0 =acell/4
        factor=3.d0*acell
        r0 = acell/4.d0
        !the following is evaluated analytically by imposing that
        !\int_0^\infty r*(-exp(-factor(r-r0)^2)+denval*exp(-factor*r^2)) = 0
        denval=sqrt(4.d0*datan(1.d0)*factor)*r0*(1.d0+derf(sqrt(factor)*r0))
        do i3=1,n03
           x3 = hz*real(i3-n03/2-1,kind=8)
           do i2=1,n02
              x2 = hy*real(i2-n02/2-1,kind=8)
              do i1=1,n01
                 x1 = hx*real(i1-n01/2-1,kind=8)
                 r2 = x1*x1+x2*x2
                 r = sqrt(r2)
                 !in this configuration denval is used so as to achieve zero monopole
                 density(i1,i2,i3) = density(i1,i2,i3) + denval*dexp(-factor*r2) - dexp(-factor*(r-r0)**2)
              end do
           end do
        end do
     case default
        denval=0.d0 !value for keeping the density positive
        do i3=1,n03
           x3 = hz*real(i3-n03/2-1,kind=8)
           call functions(x3,az,bz,fz,fz1,fz2,ifz)
           do i2=1,n02
              x2 = hy*real(i2-n02/2-1,kind=8)
              call functions(x2,ay,by,fy,fy1,fy2,ify)
              do i1=1,n01
                 x1 = hx*real(i1-n01/2-1,kind=8)
                 call functions(x1,ax,bx,fx,fx1,fx2,ifx)
                 potential(i1,i2,i3) = -fx*fy*fz              
                 density(i1,i2,i3) = (fx2*fy*fz+fx*fy2*fz+fx*fy*fz2+mu0**2*potential(i1,i2,i3))/(16.d0*datan(1.d0))
                 denval=max(denval,-density(i1,i2,i3))
              end do
           end do
        end do
     end select

  
     ! !acerioni: r = sqrt(x**2+y**2); V(x,y,z) = ArcTan(a*r)*f(z)/(a*r)
     ! do i3=1,n03
     !    x3 = hz*real(i3-n03/2-1,kind=8)
     !    call functions(x3,az,bz,fz,fz,1fz2,ifz)
     !    do i2=1,n02
     !       x2 = hy*real(i2-n02/2-1,kind=8)
     !       !call functions(x2,ay,by,fy,fy1,fy2,ify)
     !       do i1=1,n01
     !          x1 = hx*real(i1-n01/2-1,kind=8)
     !          r2 = x1*x1+x2*x2
     !          r = sqrt(r2)
     !          fxy = datan(factor*r)/(factor*r)
     !          !call functions(x1,ax,bx,fx,fx1,fx2,ifx)
     !          if (r == 0.d0) then
     !             potential(i1,i2,i3) = potential(i1,i2,i3) + 1.d0*fz
     !             density(i1,i2,i3) = density(i1,i2,i3) - fz*4.d0/3.d0*factor**2
     !             density(i1,i2,i3) = density(i1,i2,i3) + 1.d0*fz2
     !          else
     !             density(i1,i2,i3) = density(i1,i2,i3) + & 
     !                  fz*(-3.d0*factor**2/(1+factor**2*r2)**2 - 1.d0/r2/(1+factor**2*r2)**2 + fxy/r2)
     !             density(i1,i2,i3) = density(i1,i2,i3) + fxy*fz2
     !             !denval=max(denval,-density(i1,i2,i3))
     !             potential(i1,i2,i3) = potential(i1,i2,i3) + fxy*fz
     !          end if
     !          density(i1,i2,i3) = density(i1,i2,i3) / (-16.d0*datan(1.d0))
     !          !density(i1,i2,i3) = -density(i1,i2,i3)
     !       end do
     !    end do
     ! end do
     ! !acerioni

     ! !acerioni: density = delta(x,y)*"constant = 1 along z"
     ! do i3=1,n03
     !    x3 = hz*real(i3-n03/2-1,kind=8)
     !    call functions(x3,az,bz,fz,fz1,fz2,ifz)
     !    do i2=1,n02
     !       x2 = hy*real(i2-n02/2-1,kind=8)
     !       !call functions(x2,ay,by,fy,fy1,fy2,ify)
     !       do i1=1,n01
     !          x1 = hx*real(i1-n01/2-1,kind=8)
     !          r2 = x1*x1+x2*x2
     !          r = sqrt(r2)
     !          fxy = datan(factor*r)/(factor*r)
     !          !call functions(x1,ax,bx,fx,fx1,fx2,ifx)
     !          if (r == 0.d0) then
     !             potential(i1,i2,i3) = 0.d0*fz
     !             density(i1,i2,i3) = 1.d0/(hx*hy)
     !          else
     !             potential(i1,i2,i3) = - 2.d0*log(r)
     !          end if
     !          !density(i1,i2,i3) = density(i1,i2,i3) / (-16.d0*datan(1.d0))
     !       end do
     !    end do
     ! end do
     ! !acerioni

     

     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           !j1=n01/2+1-abs(n01/2+1-i1)
           !j2=n02/2+1-abs(n02/2+1-i2)
           !j3=n03/2+1-abs(n03/2+1-i3)
           write(200,*) i1,i3,density(i1,i2,i3),potential(i1,i2,i3)               
        end do
     end do



     if (ixc==0) denval=0.d0

  else

     print *,'geometry code not admitted',geocode
     stop

  end if



  !!! evaluation of the monopolar contribution !!!
  monopole = 0.d0
  dipole=0.0d0
  do i3 = 1, n03
     do i2 = 1, n02
        do i1 = 1, n01
           monopole = monopole + density(i1,i2,i3)
           dipole(1) = dipole(1) + density(i1,i2,i3)*real(i1-n01/2*hx,kind=8)
           dipole(2) = dipole(2) + density(i1,i2,i3)*real(i2-n02/2*hy,kind=8)
           dipole(3) = dipole(3) + density(i1,i2,i3)*real(i3-n03/2*hz,kind=8)
        end do
     end do
  end do
  write(*,*) 'monopole = ', monopole
  write(*,*) 'dipole = ', dipole
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! For ixc/=0 the XC potential is added to the solution, and an analytic comparison is no more
  ! possible. In that case the only possible comparison is between the serial and the parallel case
  ! To ease the comparison between the serial and the parallel case we add a random pot_ion
  ! to the potential.


  if (ixc==0) then
     potion_fac=0.d0
  else
     potion_fac=1.d0
  end if

  rhopot(:,:,:) = density(:,:,:) + denval
     do i3=1,n03
        do i2=1,n02
           do i1=1,n01
              call random_number(tt)
              !tt=0.d0!1.d0
              pot_ion(i1,i2,i3)=tt
              potential(i1,i2,i3)=potential(i1,i2,i3)+potion_fac*tt
!!              !for the ixc/=0 case
!!              call random_number(tt)
!!              rhopot(i1,i2,i3)=abs(tt)
           end do
        end do
     end do
     if (denval /= 0.d0) density=rhopot

end subroutine test_functions


subroutine functions(x,a,b,f,f1,f2,whichone)
  implicit none
  integer, intent(in) :: whichone
  real(kind=8), intent(in) :: x,a,b
  real(kind=8), intent(out) :: f,f1,f2
  !local variables
  real(kind=8) :: r,r2,y,yp,ys,factor,pi,g,h,g1,g2,h1,h2
  real(kind=8) :: length,frequency,nu,sigma,agauss

  !f1 = 0.0_dp

  pi = 4.d0*datan(1.d0)
  select case(whichone)
  case(1)
     !constant
     f=1.d0
     f1=0.d0
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
     f1=-dsin(r)*frequency*pi/length
     f2=-(frequency*pi/length)**2*dcos(r)
  case(5)
     !exp of a cosine, a=length
     nu=2.d0
     r=pi*nu/a*x
     y=dcos(r)
     yp=-dsin(r)
     f=dexp(y)/dexp(1.0_dp)
     factor=(pi*nu/a)**2*(-y+yp**2)
     f1 = f*pi*nu/a*yp
     f2 = factor*f
  case(6)
     !gaussian times "shrinked" gaussian, sigma=length/10
     length=1.d0*a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     g=dexp(-y**2)
     g1=-2.d0*y*yp*g
     g2=factor*dexp(-y**2)
     sigma=length/10.0d0
     agauss=0.5d0/sigma**2
     r2=agauss*x**2
     h=dexp(-r2)
     h1=-2.d0*agauss*x*h
     h2=(-2.d0*agauss+4.d0*agauss*r2)*dexp(-r2)
     f=g*h
     f1=g1*h+g*h1
     f2=g2*h+g*h2+2.d0*g1*h1
  case(7)
     !sine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dsin(r)
     f2=-(frequency*pi/length)**2*dsin(r)
  case(8)
     !atan with a=length, b=frequency
     length=a
     nu = length
     f=(datan(nu*x/length))**2
     f2=2.0d0*nu**2*length*(length-2.0d0*nu*x*f)/(length**2+nu**2*x**2)**2
  end select

end subroutine functions


!> Purpose: Compute exponential integral E1(x)
!! Output:  
subroutine e1xb(x,e1)
  implicit none
  !Arguments
  real(kind=8), intent(in) :: x   !< x  Argument of E1(x)
  real(kind=8), intent(out) :: e1 !< E1 --- E1(x)  ( x > 0 )
  !Local variables
  real(kind=8), parameter :: ga=0.5772156649015328d0 !< EulerGamma
  real(kind=8) :: r,t0,t
  integer :: k,m

  if (x.eq.0.0) then
     e1=1.0d+300
  else if (x.le.1.0) then
     e1=1.0d0
     r=1.0d0
     do k=1,25
        r=-r*k*x/(k+1.0d0)**2
        e1=e1+r
        if (abs(r) <= abs(e1)*1.0d-15) then 
           exit
        end if
     end do
      e1=-ga-dlog(x)+x*e1
   else
        m=20+int(80.0/x)
        t0=0.0d0
        do k=m,1,-1
           t0=k/(1.0d0+k/(x+t0))
        end do
           t=1.0d0/(x+t0)
           e1=dexp(-x)*t
        endif
        
end subroutine e1xb
      
end program PSolver_Program
