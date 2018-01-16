!> @file
!!  Program test for the stress tensor in Poisson solver
!!  May work either in parallel or in serial case
!!  And for different geometries
!! @author
!!    Copyright (C) 2006-2012 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test program for the stress tensor in Poisson solver
program PS_StressCheck
  use Poisson_Solver
  use wrapper_mpi
  use time_profiling
  use dynamic_memory
  use yaml_output
  use dictionaries
  use yaml_parse
  use yaml_strings
  use numerics
  use box
  use f_utils
  implicit none
  !include 'mpif.h'
  !Order of interpolating scaling function
  !integer, parameter :: itype_scf=8
  character(len=*), parameter :: subname='Poisson_Solver'
  real(kind=8), parameter :: a_gauss = 1.0d0, a2 = a_gauss**2
  !Length of the box
  real(kind=8), parameter :: acell = 10.0d0
  real(kind=8), parameter :: EulerGamma = 0.5772156649015328d0
  type(cell) :: mesh
  integer, parameter :: nord = 16

  !Type of function
  integer, parameter :: FUNC_CONSTANT = 1
  integer, parameter :: FUNC_GAUSSIAN = 2
  integer, parameter :: FUNC_GAUSSIAN_SHRINKED = 3
  integer, parameter :: FUNC_COSINE = 4
  integer, parameter :: FUNC_EXP_COSINE = 5
  integer, parameter :: FUNC_SHRINK_GAUSSIAN = 6
  integer, parameter :: FUNC_SINE = 7
  integer, parameter :: FUNC_ATAN = 8

  character(len=*), parameter :: inputs=&
       "- {name: ndim, shortname: n, default: 30,"//&
       "  help_string: Size of the simulation domain,"//&
       "  help_dict: {Allowed values: list of integers}}"//f_cr//&
       
       "- {name: geocode,shortname: g, default: P,"//&
       "  help_string: Boundary conditions,"//&
       "  help_dict: {Usage: set the boundary conditions of the run,"//&
       "              Allowed values: [F, S , W, P]}}"//f_cr//&
       
       "- {name: angdeg, shortname: d, default: 90.0,"//&
       "  help_string: Degrees of the angles between the directions,"//&
       "  help_dict: {Allowed values: arrays of floats}}"//f_cr//&
       
       "- {name: nstress, shortname: s, default: 1,"//&
       "  help_string: Divide the intervall (acell,2*acell) in nstress points in each direction,"//&
       "  help_dict: {Allowed values: one integer, nstress-1 SHOULD BE a divisor of n01,n02,n03}}"//f_cr//&
       
       "- {name: volstress, shortname: v, default: none,"//&
       "  help_string: Enable stress calculation varying x,y,z concurrently,"//&
       "  help_dict: {Allowed values: logical}}"//f_cr//&
       
       "- {name: input, shortname: i, default: None,"//&
       "  help_string: Inpufile of Poisson Solver,"//&
       "  help_dict: {Allowed values: dictionary in yaml format (mapping)}}"


  character(len=1) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  character(len=1) :: datacode
  character(len=30) :: mode
  real(kind=8), dimension(:,:,:), allocatable :: density,rhopot,potential,pot_ion
  real(kind=8), dimension(:), allocatable :: ene_acell,dene
  real(kind=8), dimension(:,:), allocatable :: stress_ana
  real(kind=8), dimension(:,:,:), allocatable :: stress_ps
  logical :: volstress 
  type(coulomb_operator) :: karray
  integer, dimension(3) :: ndims
  real(f_double), dimension(3) :: angdeg,angrad,hgrids
  type(dictionary), pointer :: dict
  real(kind=8) :: hx,hy,hz,max_diff,eh,exc,vxc,hgrid,diff_parser,offset,x,y
  real(kind=8) :: ehartree,eexcu,vexcu,diff_par,diff_ser,e1,acell_var,da,tr
  integer :: n01,n02,n03,unit3,unit14,unit15
  integer :: i1,i2,i3,j1,j2,j3,i1_max,i2_max,i3_max,iproc,nproc,ierr,i3sd,ncomp
  integer :: n_cell,ixc,i3xcsh,i3s,is,i,ii,ii_fin
  ! nstress-1 SHOULD BE a divisor of n01, n02 and n03 concurrently;
  ! It divided the intervall [acell,2*acell] in nstress points in each
  ! direction.
  integer :: nstress
  !> Stress tensor: Add the stress tensor part from the Hartree potential
  real(dp), dimension(6) :: stresst
  logical :: alsoserial,onlykernel
  integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
  logical :: wrtfiles=.false.
  !triclinic lattice
  real(kind=8) :: alpha,beta,gamma,detg
  real(kind=8), dimension(:,:,:,:), pointer :: rhocore_fake
  type(dictionary), pointer :: options,input
  external :: gather_timings  
  nullify(rhocore_fake)


  !mode = "charged_thin_wire"
  !mode="cylindrical_capacitor"
  !mode="monopolar"
  mode="zigzag_model_wire"
    

  call f_lib_initialize()

  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()
  call f_malloc_set_status(iproc=iproc)

  nullify(dict,input)
  call yaml_argparse(options,inputs)
  if (iproc==0) then
     call yaml_new_document()
     call yaml_map('Commandline options provided',options)
  end if
  ndims=options//'ndim'
  n01=ndims(1)
  n02=ndims(2)
  n03=ndims(3)
  ixc=0 !not needed anymore
  geocode=options//'geocode'
  angdeg=options//'angdeg'
  nstress=options//'nstress'
  volstress=options//'volstress'
  input=options .get. 'input'
  call dict_copy(dict,input) !null if absent

  call dict_free(options)

  alpha = angdeg(1)/180.0_f_double*pi!2.0_dp*datan(1.0_dp) !to be modified
  beta  = angdeg(2)/180.0_f_double*pi!2.0_dp*datan(1.0_dp)
  gamma = angdeg(3)/180.0_f_double*pi!2.0_dp*datan(1.0_dp)
  angrad(1) = angdeg(1)/180.0_f_double*pi
  angrad(2) = angdeg(2)/180.0_f_double*pi
  angrad(3) = angdeg(3)/180.0_f_double*pi

  detg = 1.0_dp - dcos(alpha)**2 - dcos(beta)**2 - dcos(gamma)**2 + 2.0_dp*dcos(alpha)*dcos(beta)*dcos(gamma)

  !write(*,*) 'mu0 =', mu0

  !perform also the comparison with the serial case
  alsoserial=.false.
  onlykernel=.false.
  !code for the Poisson Solver in the parallel case
  !datacode='G'

  hx=acell/real(n01,kind=8)
  hy=acell/real(n02,kind=8)
  hz=acell/real(n03,kind=8)
  hgrids=(/hx,hy,hz/)
  !grid for the free BC case
  hgrid=max(hx,hy,hz)

  mesh=cell_new(geocode,ndims,hgrids,alpha_bc=alpha,beta_ac=beta,gamma_ab=gamma)
   if (iproc==0) then
    call yaml_map('Angles',[alpha,beta,gamma]*180.0_dp*oneopi)
    call yaml_map('Contravariant Metric',mesh%gu)
    call yaml_map('Covariant Metric',mesh%gd)
    call yaml_map('Product of the two',matmul(mesh%gu,mesh%gd))
    call yaml_map('Covariant determinant',mesh%detgd)
   end if

   !Allocation of the energy vs acell vector
   ene_acell = f_malloc((/ nstress /),id='ene_acell')
   dene = f_malloc((/ nstress /),id='dene')
   stress_ana = f_malloc((/ nstress, 3 /),id='stress_ana')
   stress_ps = f_malloc((/ nstress, 6, 3 /),id='stress_ps')
   unit3=205
   unit14=214
   unit15=215
   if (wrtfiles) call f_open_file(unit=unit3,file='func_ene_acell.dat')

   da=1.d0
   if (nstress.gt.1)  da=acell/real(nstress-1,kind=8)

! Start of the stress code
  if (volstress) then
   ii_fin=1
  else
   ii_fin=3
  end if

  do ii=1,ii_fin ! loop on the three x, y, z components.

  do is=1,nstress

   if (iproc==0) then
     call yaml_comment('Stress itetation',hfill='-')
     call yaml_mapping_open('PS stress input')
     call yaml_map('PS stress iteration', is)
   end if

   acell_var=acell+real((is-1),kind=8)*da
   if (volstress) then
    hx=acell_var/real(n01,kind=8)
    hy=acell_var/real(n02,kind=8)
    hz=acell_var/real(n03,kind=8)
   else
    if (ii.eq.1) then
     hx=acell_var/real(n01,kind=8)
     hy=acell/real(n02,kind=8)
     hz=acell/real(n03,kind=8)
    else if (ii.eq.2) then
     hx=acell/real(n01,kind=8)
     hy=acell_var/real(n02,kind=8)
     hz=acell/real(n03,kind=8)
    else if (ii.eq.3) then
     hx=acell/real(n01,kind=8)
     hy=acell/real(n02,kind=8)
     hz=acell_var/real(n03,kind=8)
    end if
   end if
   
   hgrids=(/hx,hy,hz/)
   !grid for the free BC case
   hgrid=max(hx,hy,hz)
   !mesh=cell_new(geocode,ndims,hgrids,angrad) 


  select case(geocode)
  
  case('P')

     !if (iproc==0) print *,"PSolver, periodic BC: ",n01,n02,n03,'processes',nproc
     if (iproc==0) then
        call yaml_map('PSolver, periodic BC', (/ n01,n02,n03 /))
        call yaml_map('processes',nproc)
     end if
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
  
  case('S')

     !if (iproc==0) print *,"PSolver for surfaces: ",n01,n02,n03,'processes',nproc
     if (iproc==0) then
        call yaml_map('PSolver, surface BC', (/ n01,n02,n03/))
        call yaml_map('processes',nproc)
     end if
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  
  case('F')

     !if (iproc==0) print *,"PSolver, free BC: ",n01,n02,n03,'processes',nproc
     if (iproc==0) then
        call yaml_map('PSolver, free BC', (/ n01,n02,n03 /))
        call yaml_map('processes',nproc)
     end if
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  
  case('W')
    
     !if (iproc==0) print *,"PSolver, wires BC: ",n01,n02,n03,'processes',nproc
     if (iproc==0) then
        call yaml_map('PSolver, wire BC', (/ n01,n02,n03 /))
        call yaml_map('processes',nproc)
     end if
     call W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
     
  case('H')
   
     !if (iproc==0) print *,"PSolver, Helmholtz Equation Solver: ",n01,n02,n03,'processes',nproc
     if (iproc==0) then
        call yaml_map('PSolver, Helmholtz Equation Solver', (/ n01,n02,n03 /))
        call yaml_map('processes',nproc)
     end if
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
  
  end select

  if (iproc==0) call yaml_mapping_close()

  !write(*,*) n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3

  !initialize memory counting
  !call memocc(0,iproc,'count','start')

  !Step size
  n_cell = max(n01,n02,n03)

  !we must choose properly a test case with a positive density
  !itype_scf=16

  !write(*,'(a12,i4)') ' itype_scf = ', itype_scf 
  !call yaml_map('itype_scf',itype_scf)

  call f_timing_reset(filename='time.yaml',master=iproc==0)
  !call timing(nproc,'time.prc','IN')

  !dict => yaml_load('{kernel: {screening:'//mu0//', isf_order:'//itype_scf//'}}')
  

  karray=pkernel_init(iproc,nproc,dict,&
       geocode,(/n01,n02,n03/),(/hx,hy,hz/),&
       alpha_bc=beta,beta_ac=alpha,gamma_ab=gamma)

  mesh=karray%mesh

   if (iproc==0) then
    call yaml_map('box size',acell_var)
    call yaml_map('hgrids',hgrids)
    call yaml_map('volume element',mesh%volume_element)
   end if


  call pkernel_set(karray,verbose=.true.)

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

     if (iproc==0) call yaml_map('ixc',ixc)

     call test_functions(geocode,ixc,n01,n02,n03,acell,acell_var,a_gauss,hx,hy,hz,&
          density,potential,rhopot,pot_ion,0.0_dp,alpha,gamma,ii,volstress) !onehalf*pi,onehalf*pi,onehalf*pi)!

     !calculate expected hartree enegy
     if (wrtfiles) then
      i2=n02/2
      do i3=1,n03
         do i1=1,n01
            j1=n01/2+1-abs(n01/2+1-i1)
            j2=n02/2+1-abs(n02/2+1-i2)
            j3=n03/2+1-abs(n03/2+1-i3)
            write(110,'(2(1x,I8),2(1x,e22.15))')i1,i3,rhopot(i1,i2,i3),potential(i1,i2,i3)               
         end do
      end do
     end if

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
        offset=offset*hx*hy*hz*sqrt(mesh%detgd) ! /// to be fixed ///
        !write(*,*) 'offset = ',offset
        if (iproc==0) call yaml_map('offset',offset)
     end if

!!$     !dimension needed for allocations
!!$     call PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,(ixc>10),.false.,0,n3d,n3p,n3pi,i3xcsh,i3s)
!!$
!!$     !dimension for comparison in the global or distributed poisson solver
!!$     if (datacode == 'G') then
!!$        i3sd=1
!!$        ncomp=n03
!!$     else if (datacode == 'D') then
!!$        i3sd=i3s
!!$        ncomp=n3p
!!$     end if

     if (karray%opt%datacode=='G') then
        i3sd=1
        ncomp=n03
     else
        i3sd=karray%grid%istart+1
        ncomp=karray%grid%n3p
     end if
     i3s=i3sd
     i3xcsh=0

!!  print *,'iproc,i3xcsh,i3s',iproc,i3xcsh,i3s

!!$     print *,'density',density(25,25,25)
!!$     print *,'potential',potential(25,25,25)
!!$     print *,'rhopot',rhopot(25,25,25)
!!$     print *,'pot_ion',pot_ion(25,25,25)

     !apply the Poisson Solver (case with distributed potential)
     eexcu=0.0_gp
     vexcu=0.0_gp
     !offset=0.0_gp
     call H_potential(datacode,karray,density(1,1,i3sd),pot_ion(1,1,i3s+i3xcsh),ehartree,offset,.false.,stress_tensor=stresst)
!!$     call PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
!!$          density(1,1,i3sd),karray%kernel,pot_ion(1,1,i3s+i3xcsh),ehartree,eexcu,vexcu,offset,.true.,1,alpha,beta,gamma)

     !print *,'potential integral',sum(density)
     !this has to be corrected with the volume element of mesh
     do i=1,6
      stress_ps(is,i,ii)=stresst(i)
     end do
     eexcu=sum(density)*hx*hy*hx*sqrt(mesh%detgd)
     if (nproc >1) call mpiallred(eexcu,1,op=MPI_SUM)
     if (iproc==0) call yaml_map('potential integral',eexcu)

     if (wrtfiles) then
      i3=n03/2
      do i2=1,n02
         do i1=1,n01
            !j1=n01/2+1-abs(n01/2+1-i1)
            !j2=n02/2+1-abs(n02/2+1-i2)
            !j3=n03/2+1-abs(n03/2+1-i3)
            write(111,'(2(1x,i6),3(1x,1pe25.16e3))') i1,i2,rhopot(i1,i2,i3),potential(i1,i2,i3),density(i1,i2,i3)
            !write(111,*) i1*hx+hy*i2*dcos(alpha)+i3*hz*dcos(beta), &
            !     i2*hy*dsin(alpha)+i3*hz*(-dcos(alpha)*dcos(beta)+dcos(gamma))/dsin(alpha), &
            !     rhopot(i1,i2,i3),potential(i1,i2,i3), &
            !     density(i1,i2,i3)
         end do
      end do

      if (wrtfiles) call f_open_file(unit=unit14,file='output_x_line.dat')
      if (wrtfiles) call f_open_file(unit=unit15,file='output_y_line.dat')

      i3=n03/2+1
      i2=n02/2+1
      do i1=1,n01
       x = hx*real(i1-n01/2-1,kind=8)
       write(unit14,'(1(1x,i6),4(1x,1pe25.16e3))') i1,x,rhopot(i1,i2,i3),potential(i1,i2,i3),density(i1,i2,i3)
      end do
      i3=n03/2+1
      i1=n01/2+1
      do i2=1,n02
       y = hy*real(i2-n02/2-1,kind=8)
       write(unit15,'(1(1x,i6),4(1x,1pe25.16e3))') i2,y,rhopot(i1,i2,i3),potential(i1,i2,i3),density(i1,i2,i3)
      end do

      if (wrtfiles) call f_close(unit14)
      if (wrtfiles) call f_close(unit15)
 
      i2=n02/2
      do i3=1,n03
         !do i2=1,n02
            do i1=1,n01
               !j1=n01/2+1-abs(n01/2+1-i1)
               !j2=n02/2+1-abs(n02/2+1-i2)
               !j3=n03/2+1-abs(n03/2+1-i3)
               write(112,'(2(1x,i6),3(1x,1pe25.16e3))')i1,i3,rhopot(i1,i2,i3),potential(i1,i2,i3),density(i1,i2,i3)
            end do
         !end do
      end do
!!$     print *,'density2',density(25,25,25)
!!$     print *,'potential2',potential(25,25,25)
!!$     print *,'rhopot2',rhopot(25,25,25)
!!$     print *,'pot_ion2',pot_ion(25,25,25)
     end if
  end if

  
  call f_timing_stop(mpi_comm=karray%mpi_env%mpi_comm,nproc=karray%mpi_env%nproc,gather_routine=gather_timings)
  call pkernel_free(karray)
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
        call yaml_mapping_open('Parallel calculation')
        call yaml_map('Ehartree',ehartree)
        call yaml_map('Stress tensor',stresst)
        call yaml_map('Max diff at',[i1_max,i2_max,i3_max])
        call yaml_map('Max diff',diff_par)
        call yaml_map('Result',density(i1_max,i2_max,i3_max))
        call yaml_map('Original',potential(i1_max,i2_max,i3_max))
        call yaml_mapping_close()
     end if

     ene_acell(is)=ehartree

  end if

!  call mpibarrier()

  !Serial case
  if (alsoserial) then
     call f_timing_reset(filename='time_serial.yaml',master=iproc==0)
     !call timing(0,'             ','IN')

     karray=pkernel_init(0,1,dict,&
          geocode,(/n01,n02,n03/),(/hx,hy,hz/),&
          alpha_bc=alpha,beta_ac=beta,gamma_ab=gamma)
     call dict_free(dict)

     call pkernel_set(karray,verbose=.true.)

     call pkernel_free(karray)

     call f_timing_stop(mpi_comm=mpiworld())

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

    if (wrtfiles) then
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
  end if

  if (.not. onlykernel) then
     call f_free(density)
     call f_free(rhopot)
     call f_free(potential)
     call f_free(pot_ion)
  end if

  end do
! End of the stress loop

!post-processing of stress calculation
  call fssnord1DmatNabla('F',nstress,da,ene_acell,dene,nord)

   do is=1,nstress
    acell_var=acell+real((is-1),kind=8)*da
    if (volstress) then
     stress_ana(is,ii)=-dene(is)/(acell_var*acell_var)
    else
     stress_ana(is,ii)=-dene(is)/(acell*acell)
    end if
    if (wrtfiles) write(unit3,'(1(1x,i8),5(1x,1pe26.14e3))')is,acell_var,ene_acell(is),dene(is),&
                  stress_ana(is,ii),stress_ps(is,ii,ii)
   end do

  end do ! loop external to the stress one, for the three x,y,z directions.

  if (iproc==0) then
   call yaml_map('PS stress iteration', is)
  end if
  if (volstress) then
   tr=stress_ps(nstress/2,1,1)+stress_ps(nstress/2,2,1)+stress_ps(nstress/2,3,1)
   if (iproc == 0) then
    call yaml_comment('Stress post-processing',hfill='-')
    call yaml_mapping_open('Comparison between analytical vs psolver varing x,y,z concurrently')
    call yaml_map('Comparison at nstress/2',nstress/2)
    call yaml_map('stress analytical ',stress_ana(nstress/2,1))
    call yaml_map('stress psolver trace',tr)
    call yaml_mapping_close()
   end if
  else  
   if (iproc == 0) then
    call yaml_comment('Stress post-processing',hfill='-')
    call yaml_mapping_open('Comparison between analytical vs psolver varing x,y,z individully')
    call yaml_map('Comparison at nstress/2',nstress/2)
    call yaml_map('stress analytical x',stress_ana(nstress/2,1))
    call yaml_map('stress psolver x',stress_ps(nstress/2,1,1))
    call yaml_map('stress analytical y',stress_ana(nstress/2,2))
    call yaml_map('stress psolver y',stress_ps(nstress/2,2,2))
    call yaml_map('stress analytical z',stress_ana(nstress/2,3))
    call yaml_map('stress psolver z',stress_ps(nstress/2,3,3))
    call yaml_mapping_close()
   end if
  end if

  call f_free(ene_acell)
  call f_free(dene)
  call f_free(stress_ps)
  call f_free(stress_ana)
  call dict_free(dict)
  if (wrtfiles) call f_close(unit3)

  call mpifinalize()
  call f_lib_finalize()

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
subroutine test_functions(geocode,ixc,n01,n02,n03,acell,acell_var,a_gauss,hx,hy,hz,&
     density,potential,rhopot,pot_ion,mu0,alpha,gamma,ii,volstress)
  use yaml_output
  use f_utils
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  integer, intent(in) :: n01,n02,n03,ixc
  real(kind=8), intent(in) :: acell,acell_var,a_gauss,hx,hy,hz,mu0
  !triclinic lattice
  real(kind=8), intent(in) :: alpha,gamma
  integer, intent(in) :: ii
  real(kind=8), dimension(n01,n02,n03), intent(out) :: density,potential,rhopot,pot_ion
  logical, intent(in) :: volstress

  !local variables
  integer :: i1,i2,i3,ifx,ify,ifz,unit,unit2,unit3,unit4,unit5
  real(kind=8) :: x,x1,x2,x3,y,z,length,denval,a2,derf,factor,r,r2,r0,erfc_yy,erf_yy
  real(kind=8) :: fx,fx1,fx2,fy,fy1,fy2,fz,fz1,fz2,a,ax,ay,az,bx,by,bz,tt,potion_fac
  real(kind=8) :: monopole,kx,ky,kz,k
  real(kind=8), dimension(3) :: dipole
  logical :: gauss_dens=.false.
 
!!  !non-orthorhombic lattice
!!  real(kind=8), dimension(3,3) :: gu,gd
!!  real(kind=8) :: detg
!!
!!
!!  !triclinic cell
!!  !covariant metric
!!  gd(1,1) = 1.0_dp
!!  gd(1,2) = dcos(alpha)
!!  gd(1,3) = dcos(beta)
!!  gd(2,2) = 1.0_dp
!!  gd(2,3) = dcos(gamma)
!!  gd(3,3) = 1.0_dp
!!
!!  gd(2,1) = gd(1,2)
!!  gd(3,1) = gd(1,3)
!!  gd(3,2) = gd(2,3)
!!  !
!!  detg = 1.0_dp - dcos(alpha)**2 - dcos(beta)**2 - dcos(gamma)**2 + 2.0_dp*dcos(alpha)*dcos(beta)*dcos(gamma)
!!
!!  !write(*,*) 'detg =', detg
!!  if (iproc==0) call yaml_map('detg',detg)
!!  !
!!  !contravariant metric
!!  gu(1,1) = (dsin(gamma)**2)/detg
!!  gu(1,2) = (dcos(beta)*dcos(gamma)-dcos(alpha))/detg
!!  gu(1,3) = (dcos(alpha)*dcos(gamma)-dcos(beta))/detg
!!  gu(2,2) = (dsin(beta)**2)/detg
!!  gu(2,3) = (dcos(alpha)*dcos(beta)-dcos(gamma))/detg
!!  gu(3,3) = (dsin(alpha)**2)/detg
!!  !
!!  gu(2,1) = gu(1,2)
!!  gu(3,1) = gu(1,3)
!!  gu(3,2) = gu(2,3)
!!
!!  !gu=gd !test
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  unit=200
  unit2=201
  unit3=202
  unit4=203
  unit5=204
  if (wrtfiles) then
   call f_open_file(unit=unit,file='references_xz_vac.dat')
   call f_open_file(unit=unit2,file='references_xy_vac.dat')
   call f_open_file(unit=unit3,file='references_line_y_vac.dat')
   call f_open_file(unit=unit4,file='references_line_x_vac.dat')
   call f_open_file(unit=unit5,file='ref_funcs.dat')
  end if
  kx=1.d0
  ky=1.d0
  kz=1.d0

  if (ixc==0) denval=0.d0

  if (trim(geocode) == 'P') then
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     if (.not. gauss_dens) then
      !test functions in the three directions
      ifx=FUNC_COSINE !FUNC_SHRINK_GAUSSIAN
      ify=FUNC_COSINE !FUNC_SHRINK_GAUSSIAN
      ifz=FUNC_COSINE !FUNC_SHRINK_GAUSSIAN
      !parameters of the test functions
      ax=length
      ay=length
      az=length
!      if (ii.eq.1) then
!       ax=acell_var
!      else if (ii.eq.2) then
!       ay=acell_var
!      else if (ii.eq.3) then
!       az=acell_var
!      end if
     else
      !test functions in the three directions
      ifx=FUNC_GAUSSIAN!_SHRINKED
      ify=FUNC_GAUSSIAN!_SHRINKED
      ifz=FUNC_GAUSSIAN!_SHRINKED
      !parameters of the test functions
      ax=acell*0.05d0
      ay=acell*0.05d0
      az=acell*0.05d0
!      if (ii.eq.1) then
!       ax=acell_var*0.05d0
!      else if (ii.eq.2) then
!       ay=acell_var*0.05d0
!      else if (ii.eq.3) then
!       az=acell_var*0.05d0
!      end if
     end if

!     ifx=FUNC_EXP_COSINE
!     ify=FUNC_EXP_COSINE
!     ifz=FUNC_EXP_COSINE
!     ax=length
!     ay=length
!     az=length

     !test functions in the three directions
     ifx=FUNC_GAUSSIAN!_SHRINKED
     ify=FUNC_GAUSSIAN!_SHRINKED
     ifz=FUNC_GAUSSIAN!_SHRINKED
     !parameters of the test functions
     ax=acell*0.05d0
     ay=acell*0.05d0
     az=acell*0.05d0
     if (volstress) then
      kx=acell/acell_var
      ky=acell/acell_var
      kz=acell/acell_var
     else
      if (ii.eq.1) then
       kx=acell/acell_var
      else if (ii.eq.2) then
       ky=acell/acell_var
      else if (ii.eq.3) then
       kz=acell/acell_var
      end if
     end if
 
     !the following b's are not used, actually
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0
   
     !original version
     if (wrtfiles) then
      do i1=1,n01
         x = hx*real(i1-n01/2-1,kind=8)!valid if hy=hz
         y = hy*real(i1-n02/2-1,kind=8) 
         call functions(x,ax,kx,fx,fx1,fx2,ifx)
         call functions(y,ay,ky,fy,fy1,fy2,ify)
         write(unit5,'(1x,I8,6(1x,1pe26.14e3))') i1,x,fx,fx2,y,fy,fy2
      end do
     end if
     !Initialization of density and potential
     denval=0.d0 !value for keeping the density positive
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,kz,fz,fz1,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,ky,fy,fy1,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,kx,fx,fx1,fx2,ifx)
              potential(i1,i2,i3) =  -16.d0*datan(1.d0)*fx*fy*fz
              !density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2-mu0**2*fx*fy*fz
              !triclinic lattice
!!              density(i1,i2,i3) = -mu0**2*fx*fy*fz
!!              density(i1,i2,i3) = density(i1,i2,i3) + gu(1,1)*fx2*fy*fz+gu(2,2)*fx*fy2*fz+gu(3,3)*fx*fy*fz2
!!              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(gu(1,2)*fx1*fy1*fz+gu(1,3)*fx1*fy*fz1+gu(2,3)*fx*fy1*fz1)
              density(i1,i2,i3) = -mu0**2*fx*fy*fz
              density(i1,i2,i3) = density(i1,i2,i3) + mesh%gu(1,1)*fx2*fy*fz+mesh%gu(2,2)*fx*fy2*fz+&
                                  mesh%gu(3,3)*fx*fy*fz2
              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(mesh%gu(1,2)*fx1*fy1*fz+&
                                  mesh%gu(1,3)*fx1*fy*fz1+mesh%gu(2,3)*fx*fy1*fz1)
              if (gauss_dens) density(i1,i2,i3) = fx*fy*fz
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
     ifx=FUNC_EXP_COSINE
     ifz=FUNC_EXP_COSINE !FUNC_CONSTANT
     !non-periodic dimension
     ify=FUNC_SHRINK_GAUSSIAN
     !parameters of the test functions
     ax=length
     ay=length
     az=length
!!!     !test functions in the three directions
!!!     ifx=FUNC_GAUSSIAN !FUNC_SHRINK_GAUSSIAN
!!!     ify=FUNC_GAUSSIAN !FUNC_SHRINK_GAUSSIAN
!!!     ifz=FUNC_GAUSSIAN !FUNC_SHRINK_GAUSSIAN
!!!     !parameters of the test functions
!!!     ax=length*0.05d0
!!!     ay=length*0.05d0
!!!     az=length*0.05d0

     !b's are not used, actually
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0!real(nu,kind=8)

     call f_assert(alpha-onehalf*pi,id='Alpha angle invalid')
     call f_assert(gamma-onehalf*pi,id='Gamma angle invalid for S BC')
     
     !non-periodic dimension
     !ay=length!4.d0*a

     density(:,:,:) = 0.d0!1d-20 !added

     if (wrtfiles) then
      do i1=1,n02
         x = hx*real(i1-n02/2-1,kind=8)!valid if hy=hz
         y = hy*real(i1-n02/2-1,kind=8) 
         z = hz*real(i1-n02/2-1,kind=8)
         call functions(x,ax,bx,fx,fx1,fx2,ifx)
         call functions(y,ay,by,fy,fy1,fy2,ify)
         call functions(z,az,bz,fz,fz1,fz2,ifz)
         write(unit5,'(1x,I8,6(1x,e22.15))') i1,fx,fx2,fy,fy2,fz,fz2
      end do
     end if
     !Initialisation of density and potential
     !Normalisation
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz1,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy1,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx1,fx2,ifx)
              potential(i1,i2,i3) =  -fourpi*fx*fy*fz
!!              density(i1,i2,i3) = -mu0**2*fx*fy*fz
!!              density(i1,i2,i3) = density(i1,i2,i3) + gu(1,1)*fx2*fy*fz+gu(2,2)*fx*fy2*fz+gu(3,3)*fx*fy*fz2
!!              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(gu(1,2)*fx1*fy1*fz+gu(1,3)*fx1*fy*fz1+gu(2,3)*fx*fy1*fz1)
              density(i1,i2,i3) = -mu0**2*fx*fy*fz
              density(i1,i2,i3) = density(i1,i2,i3) + mesh%gu(1,1)*fx2*fy*fz+mesh%gu(2,2)*fx*fy2*fz+&
                                  mesh%gu(3,3)*fx*fy*fz2
              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(mesh%gu(1,2)*fx1*fy1*fz+&
                                  mesh%gu(1,3)*fx1*fy*fz1+mesh%gu(2,3)*fx*fy1*fz1)
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

     if (ixc==0) denval=0.d0

  
   else if (trim(geocode) == 'F') then

      !grid for the free BC case
      !hgrid=max(hx,hy,hz)

!      pi = 4.d0*atan(1.d0)
      a2 = a_gauss**2
    if (.false.) then

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
               density(i1,i2,i3) = factor*exp(-r2/a2)
               r = sqrt(r2)
               !Potential from a gaussian
               if (r == 0.d0) then
                  potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
               else
                  potential(i1,i2,i3) = derf(r/a_gauss)/r
               end if
            end do
         end do
      end do
    else
     !parameters for the test functions
     length=acell
     !test functions in the three directions
     ifx=FUNC_GAUSSIAN
     ifz=FUNC_GAUSSIAN
     !non-periodic dimension
     ify=FUNC_GAUSSIAN
     !parameters of the test functions
     ax=length*0.05d0
     ay=length*0.05d0
     az=length*0.05d0
     !b's are not used, actually
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0!real(nu,kind=8)

     call f_assert(alpha-onehalf*pi,id='Alpha angle invalid')
     call f_assert(gamma-onehalf*pi,id='Gamma angle invalid for S BC')
     
     !non-periodic dimension
     !ay=length!4.d0*a

     density(:,:,:) = 0.d0!1d-20 !added

     if (wrtfiles) then
      do i1=1,n02
         x = hx*real(i1-n02/2-1,kind=8)!valid if hy=hz
         y = hy*real(i1-n02/2-1,kind=8) 
         z = hz*real(i1-n02/2-1,kind=8)
         call functions(x,ax,bx,fx,fx1,fx2,ifx)
         call functions(y,ay,by,fy,fy1,fy2,ify)
         call functions(z,az,bz,fz,fz1,fz2,ifz)
         write(unit5,'(1x,I8,6(1x,1pe26.14e3))') i1,fx,fx2,fy,fy2,fz,fz2
      end do
     end if
     !Initialisation of density and potential
     !Normalisation
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz1,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy1,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx1,fx2,ifx)
              potential(i1,i2,i3) =  -fourpi*fx*fy*fz
              density(i1,i2,i3) = -mu0**2*fx*fy*fz
              density(i1,i2,i3) = density(i1,i2,i3) + mesh%gu(1,1)*fx2*fy*fz+mesh%gu(2,2)*fx*fy2*fz+&
                                  mesh%gu(3,3)*fx*fy*fz2
              density(i1,i2,i3) = density(i1,i2,i3) + 2.0_dp*(mesh%gu(1,2)*fx1*fy1*fz+&
                                  mesh%gu(1,3)*fx1*fy*fz1+mesh%gu(2,3)*fx*fy1*fz1)
              !old:
              !density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2 - mu0**2*fx*fy*fz
              denval=max(denval,-density(i1,i2,i3))
           end do
        end do
     end do

    end if
   
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
     
      denval=0.d0


  else if (trim(geocode) == 'H' .or. trim(geocode) == 'F') then

     !hgrid=max(hx,hy,hz)

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
                 call derf_local(erf_yy,r/a_gauss-a_gauss*mu0/2.0_dp)
                 erfc_yy=1.0_dp-erf_yy
                 potential(i1,i2,i3) = -2.0_dp+erfc_yy 
                 !potential(i1,i2,i3) = -2+derfc(r/a_gauss-a_gauss*mu0/2.0_dp)
                 call derf_local(erf_yy,r/a_gauss+a_gauss*mu0/2.0_dp)
                 erfc_yy=1.0_dp-erf_yy
                 potential(i1,i2,i3) = potential(i1,i2,i3)+dexp(2.0_dp*r*mu0)*erfc_yy
                 !potential(i1,i2,i3) = potential(i1,i2,i3)+dexp(2.0_dp*r*mu0)*derfc(r/a_gauss+a_gauss*mu0/2.0_dp)
                 potential(i1,i2,i3) = potential(i1,i2,i3)*a_gauss*dexp(-mu0*r)*sqrt(pi)/(-2*r*factor*dexp(-a2*mu0**2/4.0_dp))
              end if
              !density(i1,i2,i3) = exp(-r2/a2)/4.0_dp/factor**2 + 0.1_dp**2/(4*pi)*potential(i1,i2,i3)
              density(i1,i2,i3) = safe_exp(-r2/a2)/factor/(a2*pi)
           end do
        end do
     end do

     denval=0.d0

  else if (trim(geocode) == 'W') then
     !parameters for the test functions
     length=acell
     !a=0.5d0/a_gauss**2
     !test functions in the three directions
     !isolated directions
     ifx=FUNC_SHRINK_GAUSSIAN
     ify=FUNC_SHRINK_GAUSSIAN
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

     if (wrtfiles) then
      do i1=1,min(n01,n03)
        x = hx*real(i1-n01/2-1,kind=8)!isolated
        z = hz*real(i1-n03/2-1,kind=8)!periodic
        call functions(x,ax,bx,fx,fx1,fx2,ifx)
        call functions(z,az,bz,fz,fz1,fz2,ifz)
        write(unit5,*) i1,fx,fx2,fz,fz2
      end do
     end if
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

     

     if (ixc==0) denval=0.d0

  else

     !print *,'geometry code not admitted',geocode
     !stop
     call f_err_throw('geometry code not admitted "'//geocode//'"')

  end if

  k=1.d0
  if (gauss_dens) then
  k=0.d0
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           k=k+density(i1,i2,i3)
        end do
     end do
  end do
  k=k*mesh%volume_element
  end if

  !!! evaluation of the monopolar contribution !!!
  monopole = 0.d0
  dipole=0.0d0
  do i3 = 1, n03
     do i2 = 1, n02
        do i1 = 1, n01
           if (gauss_dens) density(i1,i2,i3)=density(i1,i2,i3)/k
           monopole = monopole + density(i1,i2,i3)
           dipole(1) = dipole(1) + density(i1,i2,i3)*real(i1-n01/2*hx,kind=8)
           dipole(2) = dipole(2) + density(i1,i2,i3)*real(i2-n02/2*hy,kind=8)
           dipole(3) = dipole(3) + density(i1,i2,i3)*real(i3-n03/2*hz,kind=8)
        end do
     end do
  end do
  !write(*,*) 'monopole = ', monopole
  !write(*,*) 'dipole = ', dipole
  monopole=monopole*mesh%volume_element
  dipole=dipole*mesh%volume_element

 
  if (iproc==0) then
     call yaml_map('monopole',monopole)
     call yaml_map('dipole',dipole)
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (wrtfiles) then
     i2=n02/2
     do i3=1,n03
        do i1=1,n01
           write(unit,'(2(1x,I8),2(1x,1pe26.14e3))') i1,i3,potential(i1,i2,i3),density(i1,i2,i3)
        end do
     end do
     i3=n03/2
     do i2=1,n02
        do i1=1,n01
           write(unit2,'(2(1x,I8),2(1x,1pe26.14e3))') i1,i2,potential(i1,i2,i3),density(i1,i2,i3)
        end do
     end do
     i1=n01/2
     i3=n03/2
     do i2=1,n02
      write(unit3,'(1x,I8,2(1x,1pe26.14e3))') i2,potential(i1,i2,i3),density(i1,i2,i3)
     end do
     i2=n02/2
     i3=n03/2
     do i1=1,n01
      write(unit4,'(1x,I8,2(1x,1pe26.14e3))') i1,potential(i1,i2,i3),density(i1,i2,i3)
     end do
    call f_close(unit)
    call f_close(unit2)
    call f_close(unit3)
    call f_close(unit4)
    call f_close(unit5)
  end if
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


!> Purpose: Compute exponential integral E1(x)
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
      
end program PS_StressCheck


!> Define the test functions
subroutine functions(x,a,b,f,f1,f2,whichone)
  use futile, dp => f_double
  use numerics
  implicit none
  integer, intent(in) :: whichone   !< Choose the function
  real(kind=8), intent(in) :: x     !< Argument of the function
  real(kind=8), intent(in) :: a,b   !< Parameter of the functions
  real(kind=8), intent(out) :: f    !< The value of the function
  real(kind=8), intent(out) :: f1   !< The value of the first derivative
  real(kind=8), intent(out) :: f2   !< The value of the second derivative
  !local variables
  !Type of function
  integer, parameter :: FUNC_CONSTANT = 1
  integer, parameter :: FUNC_GAUSSIAN = 2
  integer, parameter :: FUNC_GAUSSIAN_SHRINKED = 3
  integer, parameter :: FUNC_COSINE = 4
  integer, parameter :: FUNC_EXP_COSINE = 5
  integer, parameter :: FUNC_SHRINK_GAUSSIAN = 6
  integer, parameter :: FUNC_SINE = 7
  integer, parameter :: FUNC_ATAN = 8
  integer, parameter :: FUNC_ERF = 9

  real(kind=8) :: r,r2,y,yp,ys,factor,g,h,g1,g2,h1,h2,c
  real(kind=8) :: length,frequency,nu,sigma,agauss,derf

  select case(whichone)
  case(FUNC_CONSTANT)
     !constant
     f=1.d0
     f1=0.d0
     f2=0.d0
  case(FUNC_GAUSSIAN)
!!     !gaussian of sigma s.t. a=1/(2*sigma^2)
!!     r2=a*x**2
!!     f=dexp(-r2) !<checed
!!     f1=-2.d0*a*x*f !<checked
!!     f2=(-2.d0*a+4.d0*a*r2)*f !<checked
     !gaussian of sigma s.t. a=1/(2*sigma^2)
     factor=1.0d0!/(a*dsqrt(2.0d0*pi))
     c=(b/a)**2
     r2=x**2
     f=factor*safe_exp(-0.5d0*c*r2) !<checed
     f1=-f*x*c !<checked
     f2=-f*c-x*c*f1 !<checked
  case(FUNC_GAUSSIAN_SHRINKED)
     !gaussian "shrinked" with a=length of the system
     length=a
     r=b*pi*x/length
     y=tan(r)
!!$     yp=pi/length*1.d0/(dcos(r))**2
!!$     ys=2.d0*pi/length*y*yp
!!$     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
!!$     !!!!here we still need the first derivative
!!$     f2=factor*dexp(-y**2)
     f=dexp(-y**2) !<checked
     f1=-2.d0*pi*f*y/(length*cos(r)**2) !<checked
     f2=2.d0*pi**2*(2.d0*y**6 + y**4 - 2.d0*y**2 - 1.d0)/length**2*f !<checked
  case(FUNC_COSINE)
     !cosine with a=length, b=frequency
     length=a
     frequency=b*2.d0
     r=frequency*pi*x/length
     f=dcos(r) !<checked
     f1=-dsin(r)*frequency*pi/length !<checked
     f2=-(frequency*pi/length)**2*dcos(r) !<checked
  case(FUNC_EXP_COSINE)
     !exp of a cosine, a=length
     nu=2.d0
     r=b*pi*nu/a*x
     y=cos(r)
     yp=-sin(r)
     f=safe_exp(y) !<checked /dexp(1.0_dp) !<to be checked
     factor=(pi*nu/a)**2*(-y+yp**2)
     f1 = f*pi*nu/a*yp !<checked
     f2 = factor*f !<checked
  case(FUNC_SHRINK_GAUSSIAN)
     !gaussian times "shrinked" gaussian, sigma=length/10
     length=1.d0*a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     g=dexp(-y**2) !<checked
     g1=-2.d0*y*yp*g !<checked
     !g2=factor*dexp(-y**2)
     g2=2.d0*pi**2*(2.d0*y**6 + y**4 - 2.d0*y**2 - 1.d0)/length**2*g !<check
     sigma=length/10.0d0
     agauss=0.5d0/sigma**2
     r2=agauss*x**2
     h=dexp(-r2) !<checked
     h1=-2.d0*agauss*x*h !<checked
     h2=(-2.d0*agauss+4.d0*agauss*r2)*h !<checked
     f=g*h !<checked
     f1=g1*h+g*h1 !<checked
     f2=g2*h+g*h2+2.d0*g1*h1 !<checked
  case(FUNC_SINE)
     !sine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dsin(r) !<checked
     f1=frequency*pi*cos(r)/length !<checked
     f2=-(frequency*pi/length)**2*sin(r) !<checked
  case(FUNC_ATAN)
     r=a*x
     factor=r**2+1.d0
     f=atan(r) !<checked
     f1=a/factor !<checked
     f2=-2.d0*r*a**2/factor**2 !<checked
!!$     !atan with a=length, b=frequency
!!$     length=a
!!$     nu = length
!!$     f=(datan(nu*x/length))**2 !<checked
!!$     !!here first derivative is lacking
!!$     f2=2.0d0*nu**2*length*(length-2.0d0*nu*x*f)/(length**2+nu**2*x**2)**2 
  case(FUNC_ERF)
     !error function with a=sigma
     factor=sqrt(2.d0/pi)/a
     r=x
     y=x/(sqrt(2.d0)*a)
     if (abs(x)<=1.d-15) then
        f=factor
        f1=0.d0 !<checked
        f2=-sqrt(2.d0/pi)/(3.d0*a**3) !<checked
     else
        f=derf(y)/r
        y=x*x
        y=y/(2.d0*a**2)
        g=dexp(-y)
        h=1.d0/a**2+2.d0/x**2
        f1=-f/x+factor*g/x !<checked
        f2=-factor*g*h+2.d0*f/x**2  !<checked
     end if
  case default
     !print *,"Unknow function:",whichone
     !stop
     call f_err_throw('Unknown function '//trim(yaml_toa(whichone)))
  end select

end subroutine functions

!> Error function in double precision
subroutine derf_local(derf_yy,yy)
  use PSbase, only: dp
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

subroutine fssnord1DmatNabla(geocode,n01,hx,u,du,nord)
      implicit none
!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      character(len=1), intent(in) :: geocode
      integer, intent(in) :: n01,nord
      real(kind=8), intent(in) :: hx
      real(kind=8), dimension(n01) :: u
      real(kind=8), dimension(n01) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,i1,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c1DF
      logical :: perx

      n = nord+1
      m = nord/2
      n_cell = n01

      !buffers associated to the geocode
      !conditions for periodicity
      perx=(geocode /= 'F')

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
        c1DF(i,j)=0.d0
       end do
      end do

      include 'FiniteDiffCorff.inc'

      do i1=1,n01
   
       du(i1) = 0.0d0
   
       if (i1.le.m) then
        if (perx) then
         do j=-m,m
          ii=modulo(i1 + j + n01 - 1, n01 ) + 1
          du(i1) = du(i1) + c1D(j,0)*u(ii)
         end do
        else
         do j=-m,m
          du(i1) = du(i1) + c1D(j,i1-m-1)*u(j+m+1)
         end do
        end if
       else if (i1.gt.n01-m) then
        if (perx) then
         do j=-m,m
          ii=modulo(i1 + j - 1, n01 ) + 1
          du(i1) = du(i1) + c1D(j,0)*u(ii)
         end do
        else
         do j=-m,m
          du(i1) = du(i1) + c1D(j,i1-n01+m)*u(n01 + j - m)
         end do
        end if
       else
        do j=-m,m
         du(i1) = du(i1) + c1D(j,0)*u(i1 + j)
        end do
       end if
       du(i1)=du(i1)/hx
   
      end do

end subroutine fssnord1DmatNabla
