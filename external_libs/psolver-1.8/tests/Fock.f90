!> @file
!!  Program test for Fock Operator
!!  May work either in parallel or in serial case
!!  And for different geometries
!! @author
!!    Copyright (C) 2016-2016 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!> Test program for the Poisson Solver
program Fock_Operator_Program
  use Poisson_Solver
  use wrapper_mpi
  use futile, dict_set=>set
  use overlap_point_to_point
  use box
  use PSbox
  use f_blas, only: f_dot
  implicit none
  !Order of interpolating scaling function
  real(kind=8), parameter :: a_gauss = 1.0, a2 = a_gauss**2
  real(kind=8), parameter :: acell = 10.d0
  logical :: symmetric,usegpu
  character(len=1) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  real(f_double) :: offset,eexctX,tel,sfac,trm,ehartree_exp
  integer :: iproc,nproc,norb,norbu,norbd,norbpj,nspin,prc
  integer :: i_stat,igpu,iorb,jproc,ndim,ngroup,isorb
  integer, dimension(3) :: nxyz
  real(f_double), dimension(3) :: hgrids
  type(OP2P_data) :: OP2P
  type(OP2P_iterator) :: iter
  type(cell) :: mesh
  type(coulomb_operator) :: pkernel
  type(dictionary), pointer :: dict_input,options
  integer, dimension(:,:), allocatable :: nobj_par
  real(f_double), dimension(:), allocatable :: occup,spinsgn,rp_ij
  real(f_double), dimension(:,:,:), allocatable :: density,rhopot
  real(f_double), dimension(:,:,:), allocatable :: potential,pot_ion
  real(f_double), dimension(:,:,:,:), allocatable :: psir,dpsir
  external :: gather_timings  

  call f_lib_initialize() 

  !read command line
  call parse_command_line_options(options)

  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()

  !initialize categories for the Poisson Solver
  call PS_initialize_timing_categories()
  !tell to f_malloc who is the master
  call f_malloc_set_status(memory_limit=0.e0,iproc=iproc)
  call f_routine(id='Fock_Operator_Program')
  if (iproc ==0) then
     call yaml_set_stream(record_length=92,tabbing=30)!unit=70,filename='log.yaml')
     call yaml_new_document()
     call PSolver_logo()
  end if
  !initialize memory counting and timings
  call f_timing_reset(filename='time.yaml',master=iproc==0)
  nxyz=options//'ndim'
  geocode=options//'geocode'
  usegpu = options // 'accel'

  call dict_init(dict_input)
  if (usegpu) then 
    call dict_set(dict_input//'setup'//'accel','CUDA')
    call dict_set(dict_input//'setup'//'use_gpu_direct',.true.)
  end if
  !lower the verbosity for the kernel initialization
  call dict_set(dict_input //'setup'//'verbose', .false.)
  if ('input' .in. options) &
       call dict_copy(dest=dict_input,src=options//'input')

  !need norbu,norbd to derive occup and spinsgn
  nspin=dict_len(options//'norb')
  if (nspin==0) then
     nspin=1
     norbu=options//'norb'
  else if (nspin==1) then
     norbu=options//'norb'//0
     norbd=0
  else if (nspin==2) then
     norbu=options//'norb'//0
     norbd=options//'norb'//1
  end if
  norb=norbu+norbd
  call dict_free(options)
 
  !initialize the Poisson Solver main structure
  hgrids=acell/nxyz
  pkernel=pkernel_init(0,1,dict_input,geocode,nxyz,hgrids)
  call dict_free(dict_input)
  call pkernel_set(pkernel,verbose=iproc==0)

  !fill the psir functions, first find the test functions
  density = f_malloc(nxyz,id='density')
  !Density then potential
  rhopot = f_malloc(nxyz,id='rhopot')
  potential = f_malloc(nxyz,id='potential')
  !ionic potential
  pot_ion = f_malloc(nxyz,id='pot_ion')

  mesh=cell_new(geocode,nxyz,hgrids)
  call test_functions_new(mesh,1,a_gauss,& !_box
       density,potential,rhopot,pot_ion,offset)
  !calculate the expected hartree energy
  ehartree_exp=0.5_dp*f_dot(potential,density)*mesh%volume_element


  !then fill for all the functions the proposed density
  !the psi should be transformed in real space, do it within the orbital basis iterators
  psir = f_malloc0([nxyz(1),nxyz(2),nxyz(3),norbp(norb,nproc,iproc)],id='psir')
  do iorb=1,norbp(norb,nproc,iproc)
     call f_memcpy(n=product(nxyz),src=density(1,1,1),&
          dest=psir(1,1,1,iorb))
  end do

  call f_free(density)
  call f_free(potential)
  call f_free(rhopot)
  call f_free(pot_ion)
  dpsir = f_malloc0([nxyz(1),nxyz(2),nxyz(3),norbp(norb,nproc,iproc)],id='dpsir')

  occup=f_malloc(norb,id='occup')
  spinsgn=f_malloc(norb,id='spinsgn')

  occup=1.0_f_double
  spinsgn(1:norbu)=1.0_f_double
  spinsgn(norbu+1:norb)=-1.0_f_double
  !starting new approach for the exact exchange
  if (nspin==2) then
     sfac=1.0d0
     ngroup=2
  else
     sfac=0.5d0
     ngroup=1
  end if
  !construct the OP2P scheme and test it
  !use temporarily tyhe nvrct_par array
  nobj_par = f_malloc0((/ 0.to.nproc-1, 1.to.ngroup /),id='nobj_par')
  isorb=0
  do jproc=0,nproc-1
     norbpj=norbp(norb,nproc,jproc)
     !transition region
     if (isorb+norbpj > norbu .and. isorb < norbu) then
        nobj_par(jproc,1)=norbu-isorb
        if (ngroup==2) nobj_par(jproc,2)=isorb+norbpj-norbu
     else if (isorb >= norbu .and. ngroup==2) then
        nobj_par(jproc,2)=norbpj
     else
        nobj_par(jproc,1)=norbpj
     end if
     isorb=isorb+norbpj
  end do
  ndim=product(nxyz)

  symmetric=.true.
  call f_zero(dpsir)
  if (iproc==0) call yaml_map('Orbital repartition',nobj_par)
  call OP2P_unitary_test(mpiworld(),iproc,nproc,ngroup,ndim,nobj_par,symmetric)
  !this part should go inside the kernel initialization
  if(pkernel%igpu==1 .and. pkernel%initCufftPlan==0) then
     igpu=0
     if (iproc==0) call yaml_warning("not enough memory to allocate cuFFT plans on the GPU - no GPUDirect either")
  else
     igpu=.if. pkernel%use_gpu_direct .then. pkernel%igpu .else. 0
  end if
  
  call initialize_OP2P_data(OP2P,mpiworld(),iproc,nproc,ngroup,ndim,nobj_par,igpu,symmetric)

  !this part is also inaesthetic
  if(igpu==1 .and. OP2P%gpudirect==1) pkernel%stay_on_gpu=1
  !allocate work array for the internal exctx calculation
  rp_ij = f_malloc(ndim,id='rp_ij')
  eexctX=0.0d0
  !initialize the OP2P descriptor for the communication
  call set_OP2P_iterator(iproc,OP2P,iter,norbp(norb,nproc,iproc),psir,dpsir)
  !main loop
  if (iproc == 0) call yaml_newline()
  OP2P_exctx_loop: do
     call OP2P_communication_step(iproc,OP2P,iter)
     if(igpu==1) call synchronize() !this can be moved inside the communication step
     if (iter%event == OP2P_EXIT) exit OP2P_exctx_loop
     call internal_calculation_exctx(iter%istep,sfac,pkernel,norb,&
          occup,spinsgn,&
          iter%remote_result,iter%nloc_i,iter%nloc_j,iter%isloc_i,iter%isloc_j,&
          iter%phi_i,iter%phi_j,eexctX,rp_ij)
     if (iproc == 0) then
        !this part can be replaced by a progress bar
        call OP2P_info(iter,OP2P,prc,tel,trm)
        call yaml_comment('Exact exchange calculation: '+prc**'(i3)'+&
             '%; Time (s): Elapsed='+tel**'(1pg12.2)'&
             +', Remaining='+trm**'(1pg12.2)')
     end if
  end do OP2P_exctx_loop

  !we need to get the result back from the card (and synchronize with it, finally)
  if(pkernel%igpu==1 .and. pkernel%stay_on_gpu==1) then
     !here wrappers should be used
     call get_gpu_data(1, eexctX, pkernel%w%eexctX_GPU )
     !f_zero should go here
     call cudamemset(pkernel%w%eexctX_GPU,0,1,i_stat)
     if (i_stat /= 0) call f_err_throw('error cudamalloc eexctX_GPU (GPU out of memory ?) ')
     pkernel%stay_on_gpu=0
  end if
  call free_OP2P_data(OP2P)
  if (nproc > 1) call mpiallred(eexctX,1,op=MPI_SUM)
  if (iproc == 0) then
     call yaml_map('Exact Exchange Energy',eexctX,fmt='(1pe18.11)')
     call yaml_map('Expected Exchange Energy',ehartree_exp,fmt='(1pe18.11)')
  end if
  call f_free(nobj_par)
  call f_free(rp_ij)
  call f_free(dpsir)
  call f_free(psir)
  call f_free(occup,spinsgn)

  call f_timing_stop(mpi_comm=mpiworld(),nproc=nproc,gather_routine=gather_timings)

  call pkernel_free(pkernel)
  call f_release_routine()
  if (iproc==0) then
     call yaml_release_document()
     call yaml_close_all_streams()
  end if

  call mpifinalize()
  call f_lib_finalize()

contains

  pure function norbp(norb,nproc,jproc)
    implicit none
    integer, intent(in) :: jproc,norb,nproc
    integer :: norbp

    norbp=norb/nproc !integer division
    norbp=min(max(norb-jproc*norbp,0),norbp)

  end function norbp

  !>identify the options from command line
  !! and write the result in options dict
  subroutine parse_command_line_options(options)
    use yaml_parse
    use dictionaries
    implicit none
    !> dictionary of the options of the run
    !! on entry, it contains the options for initializing
    !! on exit, it contains in the key "BigDFT", a list of the 
    !! dictionaries of each of the run that the local instance of BigDFT
    !! code has to execute
    type(dictionary), pointer :: options
    !local variables
    type(yaml_cl_parse) :: parser !< command line parser

    !define command-line options
    parser=yaml_cl_parse_null()
    !between these lines, for another executable using BigDFT as a blackbox,
    !other command line options can be specified
    !then the bigdft options can be specified
    call PS_check_options(parser)
    !extra options of the Fock test
    call Fock_test_options(parser)
    !parse command line, and retrieve arguments
    call yaml_cl_parse_cmd_line(parser,args=options)
    !free command line parser information
    call yaml_cl_parse_free(parser)

  end subroutine parse_command_line_options

  subroutine Fock_test_options(parser)
  use yaml_parse
  use dictionaries, only: dict_new,operator(.is.)
  implicit none
  type(yaml_cl_parse), intent(inout) :: parser

  call yaml_cl_parse_option(parser,&
       '{name: norb,'//&
       'shortname: o,'//&
       'default: 10,'//&
       'help_string: Number of orbitals for each spin,'//&
       'help_dict: {Allowed values: list of integers}}')

end subroutine Fock_test_options

end program Fock_Operator_Program
