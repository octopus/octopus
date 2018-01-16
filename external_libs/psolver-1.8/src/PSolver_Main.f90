!> @file
!!    Main routine to perform Poisson solver calculation
!! @author
!!    Creation date: February 2007<br/>
!!    Luigi Genovese
!!    Copyright (C) 2002-2017 BigDFT group<br/>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculate the Hartree potential by solving the Poisson equation 
!! @f$\nabla^2 V(x,y,z)=-4 \pi \rho(x,y,z)@f$
!! from a given @f$\rho@f$, 
!! for different boundary conditions an for different data distributions.
!! Following the boundary conditions, it applies the Poisson Kernel previously calculated.
!!    
!! @warning
!!    The dimensions of the arrays must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
!!    Moreover, for the cases with the exchange and correlation the density must be initialised
!!    to 10^-20 and not to zero.
!!
!! @todo
!!    Wire boundary condition is missing
subroutine Electrostatic_Solver(kernel,rhov,energies,pot_ion,rho_ion,ehartree)
  use FDder
  implicit none
  !> kernel of the coulomb operator, it also contains metadata about the parallelisation scheme
  !! and the data distributions in the grid.
  !! can be also used to gather the distributed arrays for data processing or poltting purposes
  type(coulomb_operator), intent(inout) :: kernel
  !> on input, density of the (G)Pe. On output, electrostatic potential, possibly corrected with extra term in
  !! the case of rho-dependent cavity when the suitable variable of the options datatype is set. 
  !!The latter correction term is useful to define a KS DFT potential for the definition of the Hamiltonian out of 
  !!the Electrostatic environment defined from rho
  !! the last dimension is either kernel%ndims(3) or kernel%grid%n3p, depending if kernel%opt%datacode is 'G' or 'D' respectively
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),*), intent(inout) :: rhov
  !> Datatype containing the energies and th stress tensor.
  !! the components are filled accordin to the coulomb operator set ans the options given to the solver.
  type(PSolver_energies), intent(out), optional :: energies
  !> Additional external potential that is added to the output, if present.
  !! Usually represents the potential of the ions that is needed to define the full electrostatic potential of a Vacuum Poisson Equation
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%grid%n3p), intent(inout), optional, target :: pot_ion
  !> Additional external density that is added to the output input, if present.
  !! The treatment of the Poisson Equation is done with the sum of the two densities whereas the rho-dependent cavity and some components
  !! of the energies are calculated only with the input rho.
  real(dp), dimension(kernel%ndims(1), kernel%ndims(2), kernel%grid%n3p), intent(inout), optional, target :: rho_ion
  !> Electrostatic Energy of the system given as @f$\int \rho \cdot V @f$, where @f$\rho@f$ and @f$V@f$ correspond to the 
  !! values of rhov array in input and output, respectively. This value is already reduced such that each of the 
  !! MPI tasks have the same value
  real(dp), intent(out), optional :: ehartree
  !local variables
  real(dp), parameter :: max_ratioex_PB = 1.0d2
  logical :: sum_pi,sum_ri,build_c,is_vextra,plot_cavity,wrtmsg
  logical :: cudasolver,poisson_boltzmann,calc_nabla2pot,needmem
  integer :: i3start,n1,n23,i3s,i23s,i23sd2,i3sd2,i_PB,i3s_pot_pb
  real(dp) :: IntSur,IntVol,e_static,norm_nonvac,ehartreeLOC,res_PB
  type(PSolver_energies) :: energs
  real(dp), dimension(:,:), allocatable :: depsdrho,dsurfdrho
  real(dp), dimension(:,:,:), allocatable :: rhopot_full,nabla2_rhopot,delta_rho,cc_rho
  real(dp), dimension(:,:,:,:), allocatable :: nabla_rho
  real(dp), dimension(:,:,:), pointer :: pot_ion_eff,vextra_eff
  !character(len=3) :: quiet

  call f_routine(id='Electrostatic_Solver')

   energs=PSolver_energies_null()

   select case(kernel%opt%datacode)
   case('G')
      !starting address of rhopot in the case of global i/o
      i3start=kernel%grid%istart*kernel%ndims(1)*kernel%ndims(2)+1
      i23s=kernel%grid%istart*kernel%ndims(2)+1
      i3s=kernel%grid%istart+1
      if (kernel%grid%n3p == 0) then
         i3start=1
         i23s=1
         i3s=1
      end if
   case('D')
      !distributed i/o
      i3start=1
      i23s=1
      i3s=1
   case default
      call f_err_throw('datacode ("'//kernel%opt%datacode//&
           '") not admitted in PSolver')
   end select
   IntSur=0.0_gp
   IntVol=0.0_gp

   !aliasing
   n23=kernel%grid%m3*kernel%grid%n3p
   n1=kernel%grid%m1

   i3sd2=kernel%grid%istart+1
   i23sd2=kernel%grid%istart*kernel%ndims(2)+1
   if (kernel%grid%n3p==0) then
      i3sd2=1
      i23sd2=1
   end if

   poisson_boltzmann=.not. (kernel%method .hasattr. PS_PB_NONE_ENUM)
  

   wrtmsg=.true.
  select case(kernel%opt%verbosity_level)
  case(0)
     !call f_strcpy(quiet,'YES')
     wrtmsg=.false.
  case(1)
     !call f_strcpy(quiet,'NO')
     wrtmsg=.true.
  end select
  !in any case only master proc writes messages
  wrtmsg=wrtmsg .and. kernel%mpi_env%iproc==0 .and. kernel%mpi_env%igroup==0

  !decide what to do 
  sum_pi=present(pot_ion) .and. n23 > 0
  sum_ri=present(rho_ion) .and. n23 > 0
  build_c=(kernel%method .hasattr. PS_SCCS_ENUM) .and. kernel%opt%update_cavity
  calc_nabla2pot=build_c .or. kernel%opt%final_call
  is_vextra=sum_pi .or. build_c
  plot_cavity=kernel%opt%cavity_info .and. (kernel%method /= PS_VAC_ENUM)
  cudasolver= (kernel%igpu==1 .and. .not. kernel%opt%calculate_strten)

  !check of the input variables, if needed
!!$  if (kernel%method /= PS_VAC_ENUM .and. .not. present(rho_ion)) then
!!$     call f_err_throw('Error in Electrostatic_Solver, for a cavity treatment the array rho_ion is needed')
!!$  end if
  if (kernel%method == PS_PI_ENUM) then
     if (.not. associated(kernel%w%oneoeps) .or. .not. associated(kernel%w%dlogeps))&
          call f_err_throw('The PI method needs the arrays of dlogeps and oneoeps,'//&
          ' use pkernel_set_epsilon routine to set them')
  else if (kernel%method == PS_PCG_ENUM) then
     if (.not. associated(kernel%w%oneoeps) .or. .not. associated(kernel%w%corr))&
          call f_err_throw('The PCG method needs the arrays corr and oneosqrteps,'//&
          ' use pkernel_set_epsilon routine to set them')
  end if

  if (wrtmsg) then
     call yaml_mapping_open('Poisson Solver')
     select case(kernel%geocode)
     case('F')
        call yaml_map('BC','Free')
     case('P')
        call yaml_map('BC','Periodic')
     case('S')
        call yaml_map('BC','Surface')
     case('W')
        call yaml_map('BC','Wires')
     end select
     call yaml_map('Box',kernel%ndims,fmt='(i5)')
     call yaml_map('MPI tasks',kernel%mpi_env%nproc,fmt='(i5)')
     if (cudasolver) call yaml_map('GPU acceleration',.true.)
  end if
  
  !in the case of SC cavity, gather the full density and determine the depsdrho
  !here the if statement for the SC cavity should be put
  !print *,'method',trim(char(kernel%method)),associated(kernel%method%family),trim(char(kernel%method%family))
  if (calc_nabla2pot) then
     rhopot_full=f_malloc(kernel%ndims,id='rhopot_full')
     nabla2_rhopot=f_malloc(kernel%ndims,id='nabla2_rhopot')
  end if

  if (build_c) then
     delta_rho=f_malloc(kernel%ndims,id='delta_rho')
     cc_rho=f_malloc(kernel%ndims,id='cc_rho')
     nabla_rho=f_malloc([kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),3],id='nabla_rho')
     depsdrho=f_malloc([n1,n23],id='depsdrho')
     dsurfdrho=f_malloc([n1,n23],id='dsurfdrho')
     !useless for datacode= G 
     if (kernel%opt%datacode == 'D') then
        call PS_gather(rhov(1,1,i3s),kernel,dest=rhopot_full)
     else
        call f_memcpy(n=product(kernel%ndims),src=rhov(1,1,1),dest=rhopot_full(1,1,1))
     end if
     !call pkernel_build_epsilon(kernel,work_full,eps0,depsdrho,dsurfdrho)
     call rebuild_cavity_from_rho(rhopot_full,nabla_rho,nabla2_rhopot,delta_rho,cc_rho,depsdrho,dsurfdrho,&
          kernel,IntSur,IntVol)    
     call f_free(nabla_rho)
     call f_free(delta_rho)
     call f_free(cc_rho)
  end if

  if (kernel%method == PS_PI_ENUM .and. calc_nabla2pot) call f_free(rhopot_full)

  !add the ionic density to the potential, calculate also the integral
  !between the rho and pot_ion and the extra potential if present
  e_static=0.0_dp
  if (sum_ri) then
     if (sum_pi) then
        call finalize_hartree_results(.true.,cudasolver,kernel,rho_ion,&
             kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
             kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
             pot_ion,rhov(1,1,i3s),rhov(1,1,i3s),e_static)
     else
        call axpy(n1*n23,1.0_dp,rho_ion(1,1,1),1,rhov(1,1,i3s),1)
     end if
     if (wrtmsg) call yaml_map('Summing up ionic density',.true.)
  end if
   
  !inform about the quantity of charge "outside" the cavity
  if (kernel%method /= PS_VAC_ENUM) then
     call f_zero(norm_nonvac)
     if (n23 >0) call nonvacuum_projection(n1,n23,rhov(1,1,i3s),kernel%w%oneoeps,norm_nonvac)
     norm_nonvac=norm_nonvac*product(kernel%hgrids)
     call PS_reduce(norm_nonvac,kernel)
     if (wrtmsg) call yaml_map('Integral of the density in the nonvacuum region',norm_nonvac)
  end if
  
  !the allocation of the rho array is maybe not needed
  call PS_allocate_lowlevel_workarrays(poisson_boltzmann,cudasolver,&
       rhov(1,1,i3s),kernel)

  i3s_pot_pb=0  !in the PCG case
  if (kernel%method == PS_PI_ENUM) i3s_pot_pb=i23sd2-1


  !switch between neutral and ionic solution (GPe or PBe)
  if (poisson_boltzmann) then
     if (kernel%opt%use_pb_input_guess) then
        call PB_iteration(n1,n23,i3s_pot_pb,1.0_dp,kernel%cavity,rhov(1,1,i3s),kernel%w%pot,kernel%w%eps,&
             kernel%w%rho_ions,kernel%w%rho_pb,res_PB)
        if (.not. kernel%opt%use_input_guess .and. kernel%method == PS_PCG_ENUM ) call f_zero(kernel%w%pot)
     else
        call f_zero(kernel%w%rho_ions)
     end if
     if (wrtmsg) call yaml_sequence_open('Poisson Boltzmann solver')
     loop_pb: do i_PB=1,kernel%max_iter_PB
        if (wrtmsg) then
           call yaml_sequence(advance='no')
           call yaml_mapping_open()
        end if
        call Parallel_GPS(kernel,cudasolver,kernel%opt%potential_integral,energs%strten,&
             wrtmsg,kernel%w%rho_pb,kernel%opt%use_input_guess)
        
        call PB_iteration(n1,n23,i3s_pot_pb,kernel%PB_eta,kernel%cavity,rhov(1,1,i3s),kernel%w%pot,kernel%w%eps,&
             kernel%w%rho_ions,kernel%w%rho_pb,res_PB)
        if (kernel%method == PS_PCG_ENUM) call f_memcpy(src=kernel%w%rho_pb,dest=kernel%w%res)
        call PS_reduce(res_PB,kernel)
        res_PB=sqrt(res_PB/product(kernel%ndims))
        if (wrtmsg) then
           call yaml_newline()
           call EPS_iter_output(i_PB,0.0_dp,res_PB,0.0_dp,0.0_dp,0.0_dp)
           call yaml_mapping_close()
        end if
        if (res_PB < kernel%minres_PB .or. res_PB > max_ratioex_PB) exit loop_pb
     end do loop_pb
     if (wrtmsg) call yaml_sequence_close()

  else
     !call the Generalized Poisson Solver
     call Parallel_GPS(kernel,cudasolver,kernel%opt%potential_integral,energs%strten,&
          wrtmsg,rhov(1,1,i3s),kernel%opt%use_input_guess)
  end if

  !this part is not important now, to be fixed later
  if (plot_cavity) then
     if (kernel%method == PS_PCG_ENUM) then
        needmem=.not. allocated(rhopot_full)
        if (needmem) rhopot_full=f_malloc(kernel%ndims,id='rhopot_full')
        call PS_gather(src=kernel%w%pot,dest=rhopot_full,kernel=kernel)
        call polarization_charge(kernel,rhopot_full,rhov(1,1,i3s))
        if (needmem) call f_free(rhopot_full)
     end if
  end if
     
!!$  !here we should extract the information on the cavity
!!$  !--------------------------------------
!!$  ! Polarization charge, to be calculated on-demand
!!$  call pol_charge(kernel,pot_full,rho,kernel%w%pot)
!!$  !--------------------------------------

   call PS_release_lowlevel_workarrays(kernel,keep_rhopol=plot_cavity)

  !the external ionic potential is referenced if present
  if (sum_pi) then
     pot_ion_eff=>pot_ion
  else
     !point to a temporary array
     pot_ion_eff=>kernel%w%zf !<however unused
  end if
  if (is_vextra) then
     if (build_c) then
        if(.not. associated(kernel%w%zf)) then 
           kernel%w%zf = f_malloc_ptr([kernel%grid%md1, kernel%grid%md3, 2*kernel%grid%md2/kernel%mpi_env%nproc],id='zf')
        end if
        vextra_eff=>kernel%w%zf
     else
        vextra_eff=>pot_ion_eff
     end if
  else
     vextra_eff=>kernel%w%zf !hovever unused
  end if

  ehartreeLOC=0.0_gp
  select case(trim(str(kernel%method)))
  case('VAC')
     call finalize_hartree_results(present(pot_ion),cudasolver,kernel,&
          pot_ion_eff,&
          kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
          kernel%grid%md1,kernel%grid%md3,2*(kernel%grid%md2/kernel%mpi_env%nproc),&
          rhov(1,1,i3s),kernel%w%zf,rhov(1,1,i3s),ehartreeLOC)
  case('PI')
     !if statement for SC cavity
     if (calc_nabla2pot) then
        !in the PI method the potential is allocated as a full array
        call nabla_u_square(kernel%geocode,kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
             kernel%w%pot,nabla2_rhopot,kernel%nord,kernel%hgrids)
!!$        call add_Vextra(n1,n23,nabla2_rhopot(1,1,i3sd2),&
!!$             depsdrho,dsurfdrho,kernel%cavity,kernel%opt%only_electrostatic,&
!!$             sum_pi,pot_ion_eff,vextra_eff)
       end if

!!$     !here the harteee energy can be calculated and the ionic potential
!!$       !added
!!$       call finalize_hartree_results(is_vextra,cudasolver,kernel,&
!!$          vextra_eff,& 
!!$          kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
!!$          kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
!!$          rhov(1,1,i3s),kernel%w%pot(1,i23sd2),rhov(1,1,i3s),ehartreeLOC)

  case('PCG')
     !only useful for gpu, bring back the x array
     call update_pot_from_device(cudasolver, kernel, kernel%w%pot)

     if (calc_nabla2pot) then
        call PS_gather(src=kernel%w%pot,dest=rhopot_full,kernel=kernel)
        call nabla_u_square(kernel%geocode,kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
             rhopot_full,nabla2_rhopot,kernel%nord,kernel%hgrids)
        call f_free(rhopot_full)
!!$        call add_Vextra(n1,n23,nabla2_rhopot(1,1,i3sd2),&
!!$             depsdrho,dsurfdrho,kernel%cavity,kernel%opt%only_electrostatic,&
!!$             sum_pi,pot_ion_eff,vextra_eff)
     end if

!!$     !here the harteee energy can be calculated and the ionic potential
!!$     !added
!!$     call finalize_hartree_results(is_vextra,cudasolver,kernel,&
!!$          vextra_eff,&
!!$          kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
!!$          kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
!!$          rhov(1,1,i3s),kernel%w%pot,rhov(1,1,i3s),ehartreeLOC)

  end select

  if (kernel%method /= PS_VAC_ENUM) then
     !here the harteee energy can be calculated and the ionic potential
     !added

     if (build_c) then
        call add_Vextra(n1,n23,nabla2_rhopot(1,1,i3sd2),&
             depsdrho,dsurfdrho,kernel%cavity,kernel%opt%only_electrostatic,&
             sum_pi,pot_ion_eff,vextra_eff)
     end if

     call finalize_hartree_results(is_vextra,cudasolver,kernel,&
          vextra_eff,& 
          kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
          kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
          rhov(1,1,i3s),kernel%w%pot(1,i3s_pot_pb+1),rhov(1,1,i3s),ehartreeLOC)
  end if

  if (calc_nabla2pot) then
     !destroy oneoverepsilon array in the case of the forces for the rigid case
     if (kernel%opt%final_call) call nabla2pot_epsm1(n1,n23,kernel%w%eps,nabla2_rhopot(1,1,i3sd2),kernel%w%oneoeps)
     call f_free(nabla2_rhopot)
  end if


  if (build_c) then
     call f_free(depsdrho)
     call f_free(dsurfdrho)
  end if
  nullify(vextra_eff,pot_ion_eff)

  call release_PS_potential(kernel%keepzf,kernel%w,kernel%opt%use_input_guess &
       .or. (kernel%opt%use_pb_input_guess .and. poisson_boltzmann))
  
  !gather the full result in the case of datacode = G
  if (kernel%opt%datacode == 'G') call PS_gather(rhov,kernel)


  if (build_c) then
   kernel%IntSur=IntSur
   kernel%IntVol=IntVol
  end if
  if (kernel%opt%only_electrostatic .or. kernel%method == PS_VAC_ENUM) then
   kernel%IntSur=0.0_gp
   kernel%IntVol=0.0_gp
  end if
  !evaluating the total ehartree + e_static if needed
  !also cavitation energy can be given
  energs%hartree=ehartreeLOC*0.5_dp*kernel%mesh%volume_element!product(kernel%hgrids)
  energs%eVextra=e_static*kernel%mesh%volume_element!product(kernel%hgrids)
  energs%cavitation=(kernel%cavity%gammaS+kernel%cavity%alphaS)*kernel%IntSur+&
       kernel%cavity%betaV*kernel%IntVol

  call PS_reduce(energs,kernel)

  if (present(energies)) then
     energies=energs
  end if

  if (present(ehartree)) then
     ehartree=energs%hartree
  end if

  if (wrtmsg) then
     if (kernel%cavity%gammaS*kernel%IntSur /= 0.0_gp) &
          call yaml_map('Cavitation energy',kernel%cavity%gammaS*kernel%IntSur)
     if (kernel%cavity%alphaS*kernel%IntSur /= 0.0_gp) &
          call yaml_map('Repulsion energy',kernel%cavity%alphaS*kernel%IntSur)
     if (kernel%cavity%betaV*kernel%IntVol /= 0.0_gp) &
          call yaml_map('Dispersion energy',kernel%cavity%betaV*kernel%IntVol)
     if (energs%cavitation /= 0.0_gp) &
          call yaml_map('Non-eletrostatic energy',energs%cavitation)
  end if

  if (wrtmsg) call yaml_mapping_close()

  call f_release_routine()
end subroutine Electrostatic_Solver


!>Generalized Poisson Solver working with parallel data distribution.
!!should not allocate memory as the memory handling is supposed to be done 
!!outside. the only exception for the moment is represented by
!! apply_kernel as the inner poisson solver still allocates heap memory.
subroutine Parallel_GPS(kernel,cudasolver,offset,strten,wrtmsg,rho_dist,use_input_guess)
  use FDder!, only: update_rhopol
  implicit none
  logical, intent(in) :: cudasolver,wrtmsg
  real(dp), intent(in) :: offset
  type(coulomb_operator), intent(inout) :: kernel
  !>density in input, distributed version. Not modified even if written as inout
  real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(inout) :: rho_dist
  real(dp), dimension(6), intent(out) :: strten !<stress tensor
  logical, intent(in) :: use_input_guess
  !local variables
  integer, parameter :: NEVER_DONE=0,DONE=1
  real(dp), parameter :: max_ratioex = 1.0e10_dp !< just to avoid crazy results
  real(dp), parameter :: no_ig_minres=1.0d-2 !< just to neglect wrong inputguess
  integer :: n1,n23,i1,i23,ip,i23s,iinit
  real(dp) :: rpoints,rhores2,beta,ratio,normr,normb,alpha,q
  !aliasings
  call f_timing(TCAT_PSOLV_COMPUT,'ON')
  rpoints=product(real(kernel%ndims,dp))

  n23=kernel%grid%m3*kernel%grid%n3p
  n1=kernel%grid%m1

  !now switch the treatment according to the method used
  select case(trim(str(kernel%method)))
  case('VAC')
     !initalise to zero the zf array 
     !call f_zero(kernel%w%zf)
     !core psolver routine
     call apply_kernel(cudasolver,kernel,rho_dist,offset,strten,kernel%w%zf,.false.)
  case('PI')

     if (wrtmsg) &
          call yaml_sequence_open('Embedded PSolver Polarization Iteration Method')

     if (use_input_guess) then
        !gathering the data to obtain the distribution array
        !call PS_gather(kernel%w%pot,kernel) !not needed as in PI the W%pot array is global
        call update_rhopol(kernel%mesh,kernel%w%pot,kernel%nord,1.0_dp,kernel%w%eps,&
             kernel%w%dlogeps,kernel%w%rho,rhores2)
     end if

     ip=0
     iinit=NEVER_DONE
     pi_loop: do while (ip <= kernel%max_iter)
      ip=ip+1
!     pi_loop: do ip=1,kernel%max_iter
        !update the needed part of rhopot array
        !irho=1
        !i3s=kernel%grid%istart
        i23s=kernel%grid%istart*kernel%grid%m3
        do i23=1,n23
           do i1=1,n1
              kernel%w%pot(i1,i23+i23s)=& !this is full
              !rhopot(irho)=&
                   kernel%w%oneoeps(i1,i23)*rho_dist(i1,i23)+&
                   kernel%w%rho(i1,i23+i23s) !this is full
              kernel%w%rho_pol(i1,i23)=kernel%w%pot(i1,i23+i23s)-rho_dist(i1,i23)
              !irho=irho+1
           end do
        end do

        !initalise to zero the zf array 
        !call f_zero(kernel%w%zf)
        call apply_kernel(cudasolver,kernel,kernel%w%pot(1,i23s+1),&
             offset,strten,kernel%w%zf,.true.)
        !gathering the data to obtain the distribution array
        call PS_gather(kernel%w%pot,kernel)

        !update rhopol and calculate residue
        !reduction of the residue not necessary
        call update_rhopol(kernel%mesh,kernel%w%pot,kernel%nord,kernel%PI_eta,kernel%w%eps,&
             kernel%w%dlogeps,kernel%w%rho,rhores2)

        rhores2=sqrt(rhores2/rpoints)

        if (wrtmsg) then 
           call yaml_newline()
           call yaml_sequence(advance='no')
           call EPS_iter_output(ip,0.0_dp,rhores2,0.0_dp,0.0_dp,0.0_dp)
        end if

        if (rhores2 < kernel%minres) exit pi_loop

        if ((rhores2 > no_ig_minres).and.use_input_guess .and. iinit==NEVER_DONE) then
           call f_zero(kernel%w%rho)
!!$           i23s=kernel%grid%istart*kernel%grid%m3
!!$           do i23=1,n23
!!$              do i1=1,n1
!!$                 kernel%w%rho(i1,i23+i23s)=0.0_dp !this is full
!!$              end do
!!$           end do
           if (kernel%mpi_env%iproc==0) &
                call yaml_warning('Input guess not used due to residual norm >'+no_ig_minres)
           iinit=DONE
        end if
     end do pi_loop
     if (wrtmsg) call yaml_sequence_close()
  case('PCG')
     if (wrtmsg) then
        call yaml_newline()
        call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
     end if

     normb=dot(n1*n23,kernel%w%res(1,1),1,kernel%w%res(1,1),1)
     call PS_reduce(normb,kernel)
     normb=sqrt(normb/rpoints)
     normr=normb
     
     iinit=1
     if (use_input_guess) then
      iinit=2
 
      !$omp parallel do default(shared) private(i1,i23,q)
      do i23=1,n23
         do i1=1,n1
            q=kernel%w%pot(i1,i23)*kernel%w%corr(i1,i23)
            kernel%w%q(i1,i23)=kernel%w%pot(i1,i23)
            kernel%w%z(i1,i23)=(kernel%w%res(i1,i23)-q)*&
                                 kernel%w%oneoeps(i1,i23)
         end do
      end do
      !$omp end parallel do
 
      call apply_kernel(cudasolver,kernel,kernel%w%z,offset,strten,kernel%w%zf,.true.)

      normr=0.0_dp
      !$omp parallel do default(shared) private(i1,i23) &
      !$omp reduction(+:normr)
      do i23=1,n23
         do i1=1,n1
            kernel%w%pot(i1,i23)=kernel%w%z(i1,i23)*kernel%w%oneoeps(i1,i23)
            kernel%w%z(i1,i23)=kernel%w%res(i1,i23)
            kernel%w%res(i1,i23)=(kernel%w%q(i1,i23) - kernel%w%pot(i1,i23))*kernel%w%corr(i1,i23)
!!$            if (kernel%w%res(i1,i23) /= 0.0_dp) then
!!$               print *,'facs',kernel%w%q(i1,i23)
!!$               print *,'facs2',kernel%w%pot(i1,i23)
!!$               print *,'facs3',kernel%w%corr(i1,i23)
!!$               print *,'here',kernel%w%res(i1,i23)
!!$               print *,'there',kernel%w%res(i1,i23)**2
!!$            end if
            kernel%w%q(i1,i23)=0.d0
            normr=normr+kernel%w%res(i1,i23)**2
         end do
      end do
      !$omp end parallel do
      !print *,'here',sum(kernel%w%res)
      !normr=dot(n1*n23,kernel%w%res(1,1),1,kernel%w%res(1,1),1)
      call PS_reduce(normr,kernel)
      normr=sqrt(normr/rpoints)
      ratio=normr/normb
 
      if (wrtmsg) then
       call yaml_newline()
       call yaml_sequence(advance='no')
       call EPS_iter_output(1,0.0_dp,normr,0.0_dp,0.0_dp,0.0_dp)
      end if
      if (normr < kernel%minres) iinit=kernel%max_iter+10
 
      if (normr > no_ig_minres ) then
 
       !wipe out the IG information as it turned out to be useless
       !$omp parallel do default(shared) private(i1,i23)
       do i23=1,n23
          do i1=1,n1
             kernel%w%res(i1,i23)=kernel%w%z(i1,i23)
             kernel%w%pot(i1,i23)=0.d0
          end do
       end do
       !$omp end parallel do
 
       iinit=1
       !call yaml_warning('Input guess not used due to residual norm > 1')
       if (kernel%mpi_env%iproc==0) &
            call yaml_warning('Input guess not used due to residual norm >'+no_ig_minres)
      end if

     end if

     beta=1.d0
     ratio=1.d0

     !$omp parallel do default(shared) private(i1,i23)
     do i23=1,n23
        do i1=1,n1
           kernel%w%z(i1,i23)=kernel%w%res(i1,i23)*kernel%w%oneoeps(i1,i23)
        end do
     end do
     !$omp end parallel do

     PCG_loop: do ip=iinit,kernel%max_iter

        !initalise to zero the zf array 
        !call f_zero(kernel%w%zf)
        !  Apply the Preconditioner
        call apply_kernel(cudasolver,kernel,kernel%w%z,offset,strten,kernel%w%zf,.true.)

        call apply_reductions(ip,cudasolver,kernel,&
             kernel%w%res,kernel%w%pot,kernel%w%p,kernel%w%q,kernel%w%z,&
             alpha,beta,normr)

        normr=sqrt(normr/rpoints)

        ratio=normr/normb
        if (wrtmsg) then
           call yaml_newline()
           call yaml_sequence(advance='no')
           !call EPS_iter_output(ip,normb,normr,ratio,alpha,beta)
           call EPS_iter_output(ip,0.0_dp,normr,0.0_dp,0.0_dp,0.0_dp)
        end if
        if (normr < kernel%minres .or. normr > max_ratioex) exit PCG_loop
     end do PCG_loop
     if (wrtmsg) call yaml_sequence_close()
  end select
  call f_timing(TCAT_PSOLV_COMPUT,'OF')  
end subroutine Parallel_GPS

subroutine H_potential(datacode,kernel,rhopot,pot_ion,eh,offset,sumpion,&
     quiet,rho_ion,stress_tensor) !optional argument
  use FDder
   implicit none

   !> kernel of the Poisson equation. It is provided in distributed case, with
   !! dimensions that are related to the output of the PS_dim4allocation routine
   !! it MUST be created by following the same geocode as the Poisson Solver.
   !! it is declared as inout as in the cavity treatment it might be modified
   type(coulomb_operator), intent(inout) :: kernel

   !> @copydoc poisson_solver::doc::datacode
   character(len=1), intent(in) :: datacode

   !> Logical value which states whether to sum pot_ion to the final result or not
   !!   .true.  rhopot will be the Hartree potential + pot_ion
   !!           pot_ion will be untouched
   !!   .false. rhopot will be only the Hartree potential
   !!           pot_ion will be ignored
   logical, intent(in) :: sumpion

   !> Total integral on the supercell of the final potential on output
   real(dp), intent(in) :: offset

   !> Hartree Energy (Hartree)
   real(gp), intent(out) :: eh

   !> On input, it represents the density values on the grid points
   !! On output, it is the Hartree potential (and maybe also pot_ion)
   real(dp), dimension(*), intent(inout) :: rhopot

   !> Additional external potential that is added to the output, 
   !! when the XC parameter ixc/=0 and sumpion=.true.
   !! When sumpion=.true., it is always provided in the distributed form,
   !! clearly without the overlapping terms which are needed only for the XC part
   real(dp), dimension(*), intent(inout) :: pot_ion

   !> Optional argument to avoid output writings
   character(len=3), intent(in), optional :: quiet
   
   !> Additional external density added to the input, regardless of the value of sumpion,
   !! without the overlapping terms which are needed only for the XC part
   real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(inout), optional :: rho_ion


   !> Stress tensor: Add the stress tensor part from the Hartree potential
   real(dp), dimension(6), intent(out), optional :: stress_tensor

   !local variables
   character(len=*), parameter :: subname='H_potential'
   real(dp), parameter :: max_ratioex = 1.0e10_dp,eps0=78.36d0 !to be inserted in pkernel
   logical :: wrtmsg,cudasolver,global,verb
   integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
   integer :: ierr,ind,ind2,ind3,indp,ind2p,ind3p,i
   integer :: i1,i2,i3,j2,istart,iend,i3start,jend,jproc,i3xcsh
   integer :: nxc,istden,istglo,ip,irho,i3s,i23,i23s,n23
   integer(f_integer) :: ierr_4
   real(dp) :: ehartreeLOC,pot,rhores2,beta
   real(dp) :: alpha,ratio,normb,normr,norm_nonvac,e_static,rpoints
   !real(dp) :: scal
   real(dp), dimension(6) :: strten
   real(dp), dimension(:,:), allocatable :: rho,rhopol,q,p,z,depsdrho,dsurfdrho
   real(dp), dimension(:,:,:), allocatable :: work_full,pot_full
   !integer, dimension(:,:), allocatable :: gather_arr
   type(PSolver_energies) :: energies

   !kernel%opt=PSolver_options_null()
   global=.false.
   global=datacode=='G'
   verb=.true.
   if (present(quiet)) then
      if ((quiet .eqv. 'yes') .and. kernel%method == PS_VAC_ENUM) &
           verb=.false.
   end if

   call PS_set_options(kernel,global_data=global,&
        calculate_strten=present(stress_tensor),verbose=verb,&
        update_cavity=kernel%method .hasattr. PS_SCCS_ENUM,&
        potential_integral=offset,&
        final_call=(kernel%method .hasattr. 'rigid') .and.&
        present(stress_tensor))

   if (sumpion .and. present(rho_ion) .and. kernel%method /= PS_VAC_ENUM) then
      call Electrostatic_Solver(kernel,rhopot,energies,pot_ion,rho_ion)
   else if (sumpion) then
      call Electrostatic_Solver(kernel,rhopot,energies,pot_ion)
   else if (present(rho_ion)) then
      call Electrostatic_Solver(kernel,rhopot,energies,rho_ion=rho_ion)
   else
      call Electrostatic_Solver(kernel,rhopot,energies)
   end if

   !retrieve the energy and stress
   !eh=energies%hartree+energies%eVextra
   eh=energies%hartree+energies%eVextra+energies%cavitation
   if (present(stress_tensor)) stress_tensor=energies%strten

END SUBROUTINE H_potential

subroutine extra_sccs_potential(kernel,work_full,depsdrho,dsurfdrho,pot,eps0)
  implicit none
  type(coulomb_operator), intent(in) :: kernel
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(out) :: work_full
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(inout) :: depsdrho
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(in) :: dsurfdrho
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p) :: pot !intent in
  real(dp), intent(in) :: eps0

  !first gather the potential to calculate the derivative
  if (kernel%mpi_env%nproc > 1) then
     call mpiallgather(pot,recvbuf=work_full,recvcounts=kernel%counts,&
          displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
  else
     call f_memcpy(src=pot,dest=work_full)
  end if

  !then calculate the extra potential and add it to pot
  call sccs_extra_potential(kernel,work_full,depsdrho,dsurfdrho,eps0)
  
end subroutine extra_sccs_potential

subroutine pol_charge(kernel,pot_full,rho,pot)
  implicit none
  type(coulomb_operator), intent(inout) :: kernel
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(out) :: pot_full
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(inout) :: rho
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p) :: pot !intent in

  !first gather the potential to calculate the derivative
  if (kernel%mpi_env%nproc > 1) then
     call mpiallgather(pot,recvbuf=pot_full,recvcounts=kernel%counts,&
          displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
  else
     call f_memcpy(src=pot,dest=pot_full)
  end if

  !calculate the extra potential and add it to pot
  call polarization_charge(kernel,pot_full,rho)
  
end subroutine pol_charge

subroutine apply_reductions(ip, gpu, kernel, r, x, p, q, z, alpha, beta, normr)
  use f_utils, only: f_zero
  use time_profiling, only: f_timing
  use yaml_output
  implicit none
  logical, intent(in) :: gpu !< logical variable controlling the gpu acceleration
  type(coulomb_operator), intent(in) :: kernel 
  real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(inout) :: r,x,p,q,z
  real(dp), intent(inout) :: beta
  real(dp), intent(out) :: alpha,normr
  !local variables
  integer :: n1,n23,i_stat,ierr,i23,i1,size1, ip
  real(dp) :: zeta, rval, beta0, epsc, kappa, pval, qval

  n23=kernel%grid%m3*kernel%grid%n3p
  n1=kernel%grid%m1
  beta0 = beta
  beta=0.d0
  if (kernel%gpuPCGRed == 0) then !CPU case
         !$omp parallel do default(shared) private(i1,i23,rval,zeta) &
         !$omp reduction(+:beta)
         do i23=1,n23
            do i1=1,n1
               zeta=z(i1,i23)
               zeta=zeta*kernel%w%oneoeps(i1,i23)
               rval=r(i1,i23)
               rval=rval*zeta
               beta=beta+rval
               z(i1,i23)=zeta
            end do
         end do
         !$omp end parallel do
         call PS_reduce(beta,kernel)
         kappa=0.d0
         !$omp parallel do default(shared) private(i1,i23,epsc,zeta)&
         !$omp private(pval,qval,rval) reduction(+:kappa)
         do i23=1,n23
            do i1=1,n1
               zeta=z(i1,i23)
               epsc=kernel%w%corr(i1,i23)
               pval=p(i1,i23)
               qval=q(i1,i23)
               rval=r(i1,i23)
               pval = zeta+(beta/beta0)*pval
               qval = zeta*epsc+rval+(beta/beta0)*qval
               p(i1,i23) = pval
               q(i1,i23) = qval
               rval=pval*qval
               kappa = kappa+rval 
            end do
         end do
         !$omp end parallel do
         call PS_reduce(kappa,kernel)

         alpha = beta/kappa

         normr=0.d0
         !$omp parallel do default(shared) private(i1,i23,rval) &
         !$omp reduction(+:normr)
         do i23=1,n23
            do i1=1,n1
               x(i1,i23) = x(i1,i23) + alpha*p(i1,i23)
               r(i1,i23) = r(i1,i23) - alpha*q(i1,i23)
               z(i1,i23) = r(i1,i23)*kernel%w%oneoeps(i1,i23)
               rval=r(i1,i23)*r(i1,i23)
               normr=normr+rval
            end do
         end do
         !$omp end parallel do
         call PS_reduce(normr,kernel)

  else
!naive method, with allocations/free at each time .. 
!may need to store more pointers inside kernel
  size1=n1*n23
  if (kernel%keepGPUmemory == 0) then
    call cudamalloc(size1,kernel%w%z_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%w%r_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%w%oneoeps_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%w%p_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%w%q_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%w%x_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%w%corr_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat


    call cudamalloc(sizeof(alpha),kernel%w%alpha_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(sizeof(beta),kernel%w%beta_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(sizeof(beta0),kernel%w%beta0_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(sizeof(kappa),kernel%w%kappa_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
  end if

  if (ip == 1) then 
    call reset_gpu_data(size1,r,kernel%w%r_GPU)
    call reset_gpu_data(size1,kernel%w%oneoeps,kernel%w%oneoeps_GPU)
    call cudamemset(kernel%w%p_GPU, 0, size1,i_stat)
    if (i_stat /= 0) print *,'error cudamemset p',i_stat
    call cudamemset(kernel%w%q_GPU, 0, size1,i_stat)
    if (i_stat /= 0) print *,'error cudamemset q',i_stat
    call cudamemset(kernel%w%x_GPU, 0, size1,i_stat)
    if (i_stat /= 0) print *,'error cudamemset x',i_stat
  end if

  call reset_gpu_data(size1,z,kernel%w%z_GPU)
  call cudamemset(kernel%w%beta_GPU, 0, sizeof(beta),i_stat)
  call first_reduction_kernel(n1,n23,kernel%w%p_GPU,kernel%w%q_GPU,kernel%w%r_GPU,&
kernel%w%x_GPU,kernel%w%z_GPU,kernel%w%corr_GPU, kernel%w%oneoeps_GPU, kernel%w%alpha_GPU,&
 kernel%w%beta_GPU, kernel%w%beta0_GPU, kernel%w%kappa_GPU, kernel%w%reduc_GPU, beta)

  call cudamemset(kernel%w%kappa_GPU, 0, sizeof(kappa),i_stat)

!TODO : gpudirect.
  call PS_reduce(beta,kernel)
         kappa=0.d0

  call reset_gpu_data(sizeof(beta),beta, kernel%w%beta_GPU)
  call reset_gpu_data(sizeof(beta0),beta0, kernel%w%beta0_GPU)

  if (ip == 1) then 
    call reset_gpu_data(size1,kernel%w%corr,kernel%w%corr_GPU)
  end if

  call second_reduction_kernel(n1,n23,kernel%w%p_GPU,kernel%w%q_GPU,kernel%w%r_GPU,&
kernel%w%x_GPU,kernel%w%z_GPU,kernel%w%corr_GPU, kernel%w%oneoeps_GPU, kernel%w%alpha_GPU,&
 kernel%w%beta_GPU, kernel%w%beta0_GPU, kernel%w%kappa_GPU, kernel%w%reduc_GPU, kappa)

!  call get_gpu_data(size1,p,kernel%p_GPU)
!  call get_gpu_data(size1,q,kernel%q_GPU)

  call PS_reduce(kappa,kernel)

  alpha = beta/kappa
  call reset_gpu_data(sizeof(alpha),alpha,kernel%w%alpha_GPU)
  normr=0.d0

  call third_reduction_kernel(n1,n23,kernel%w%p_GPU,kernel%w%q_GPU,kernel%w%r_GPU,&
kernel%w%x_GPU,kernel%w%z_GPU,kernel%w%corr_GPU, kernel%w%oneoeps_GPU, kernel%w%alpha_GPU,&
 kernel%w%beta_GPU, kernel%w%beta0_GPU, kernel%w%kappa_GPU, kernel%w%reduc_GPU, normr)

  call PS_reduce(normr,kernel)
  
  call get_gpu_data(size1,z,kernel%w%z_GPU)


  if (kernel%keepGPUmemory == 0) then
    call cudafree(kernel%w%z_GPU)
    call cudafree(kernel%w%r_GPU)
    call cudafree(kernel%w%oneoeps_GPU)
    call cudafree(kernel%w%p_GPU)
    call cudafree(kernel%w%q_GPU)
    call cudafree(kernel%w%x_GPU)
    call cudafree(kernel%w%corr_GPU)
    call cudafree(kernel%w%alpha_GPU)
    call cudafree(kernel%w%beta_GPU)
    call cudafree(kernel%w%beta0_GPU)
    call cudafree(kernel%w%kappa_GPU)
  end if

  end if

end subroutine apply_reductions

!at the end of the loop, we have to synchronize potential from gpu
subroutine update_pot_from_device(gpu, kernel,x)
  type(coulomb_operator), intent(in) :: kernel 
  real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(inout) :: x
  logical, intent(in) :: gpu !< logical variable controlling the gpu acceleration
  integer size1
!  if (.false.) then !CPU case
  if (kernel%gpuPCGRed==1) then !CPU case
    size1=kernel%grid%m3*kernel%grid%n3p*kernel%grid%m1
    call get_gpu_data(size1,x,kernel%w%x_GPU)
  end if 
end subroutine update_pot_from_device


subroutine EPS_iter_output(iter,normb,normr,ratio,alpha,beta)
  implicit none
  integer, intent(in) :: iter
  real(dp), intent(in) :: normb,normr,ratio,beta,alpha
  !local variables
  character(len=*), parameter :: vrb='(1pe25.17)'!'(1pe16.4)'

  !call yaml_newline()
  call yaml_mapping_open('Iteration quality',flow=.true.)
  if (beta /= 0.0_dp) call yaml_comment('GPS Iteration '+iter,hfill='_')
  !write the PCG iteration
  call yaml_map('iter',iter,fmt='(i4)')
  !call yaml_map('rho_norm',normb)
  if (normr/=0.0_dp) call yaml_map('res',normr,fmt=vrb)
  if (ratio /= 0.0_dp) call yaml_map('ratio',ratio,fmt=vrb)
  if (alpha /= 0.0_dp) call yaml_map('alpha',alpha,fmt=vrb)
  if (beta /= 0.0_dp) call yaml_map('beta',beta,fmt=vrb)

  call yaml_mapping_close()
end subroutine EPS_iter_output

!> Plot in different files the available information related to the coulomb operator.
!! The results is adapted to the actual setup of the calculation.
subroutine PS_dump_coulomb_operator(kernel,prefix)
  use IObox
  use yaml_output
  implicit none
  type(coulomb_operator), intent(in) :: kernel
  character(len=*), intent(in) :: prefix
  !local variables
  logical :: master
  real(dp), dimension(:,:,:), allocatable :: global_arr

  master=kernel%mpi_env%iproc==0 .and. kernel%mpi_env%igroup==0

  !if available plot the dielectric function
  if (kernel%method /= 'VAC') then
     if (master) call yaml_map('Writing dielectric cavity in file','dielectric_cavity')
     global_arr = f_malloc(kernel%ndims,id='global_arr')
     call PS_gather(src=kernel%w%eps,dest=global_arr,kernel=kernel)
     if (kernel%method .hasattr. PS_RIGID_ENUM) then
        !we might add the atoms in the case of a rigid cavity, they are in
        call dump_field(trim(prefix)//'dielectric_cavity',&
             kernel%geocode,kernel%ndims,kernel%hgrids,1,&
             global_arr,rxyz=kernel%w%rxyz)
     else
        !charge dependent case
        call dump_field(trim(prefix)//'dielectric_cavity',&
             kernel%geocode,kernel%ndims,kernel%hgrids,1,global_arr)
     end if
     !now check if the polarization charge is available
     if (associated(kernel%w%rho_pol)) then
        call PS_gather(src=kernel%w%rho_pol,dest=global_arr,kernel=kernel)
        call dump_field(trim(prefix)//'polarization_charge',&
             kernel%geocode,kernel%ndims,kernel%hgrids,1,global_arr)
     end if
     !now check if the ionic charge of the Poisson boltzmann charge is available
     if (associated(kernel%w%rho_ions)) then
        call PS_gather(src=kernel%w%rho_ions,dest=global_arr,kernel=kernel)
        call dump_field(trim(prefix)//'solvent_density',&
             kernel%geocode,kernel%ndims,kernel%hgrids,1,global_arr)
     end if
     call f_free(global_arr)
  end if
  
end subroutine PS_dump_coulomb_operator

!> Calculate the dimensions needed for the allocation of the arrays 
!! related to the Poisson Solver
!!
!! @warning
!!    The XC enlarging due to GGA part is not present for surfaces and 
!!    periodic boundary condition. This is related to the fact that the calculation of the
!!    gradient and the White-Bird correction are not yet implemented for non-isolated systems
!! @author Luigi Genovese
!! @date February 2007
subroutine PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,use_gradient,use_wb_corr,&
      igpu,n3d,n3p,n3pi,i3xcsh,i3s)
   implicit none
   character(len=1), intent(in) :: geocode  !< @copydoc poisson_solver::doc::geocode
   character(len=1), intent(in) :: datacode !< @copydoc poisson_solver::doc::datacode
   integer, intent(in) :: iproc        !< Process Id
   integer, intent(in) :: nproc        !< Number of processes
   integer, intent(in) :: n01,n02,n03  !< Dimensions of the real space grid to be hit with the Poisson Solver
   logical, intent(in) :: use_gradient !< .true. if functional is using the gradient.
   logical, intent(in) :: use_wb_corr  !< .true. if functional is using WB corrections.
   integer, intent(in) :: igpu         !< Is GPU enabled ?
   !> Third dimension of the density. For distributed data, it takes into account 
   !! the enlarging needed for calculating the XC functionals.
   !! For global data it is simply equal to n03. 
   !! When there are too many processes and there is no room for the density n3d=0
   integer, intent(out) :: n3d
   !> Third dimension for the potential. The same as n3d, but without 
   !! taking into account the enlargment for the XC part. For non-GGA XC, n3p=n3d.
   integer, intent(out) :: n3p
   !> Dimension of the pot_ion array, always with distributed data. 
   !! For distributed data n3pi=n3p
   integer, intent(out) :: n3pi
   !> Shift of the density that must be performed to enter in the 
   !! non-overlapping region. Useful for recovering the values of the potential
   !! when using GGA XC functionals. If the density starts from rhopot(1,1,1),
   !! the potential starts from rhopot(1,1,i3xcsh+1). 
   !! For non-GGA XCs and for global distribution data i3xcsh=0
   integer, intent(out) :: i3xcsh
   !> Starting point of the density effectively treated by each processor 
   !! in the third direction.
   !! It takes into account also the XC enlarging. The array rhopot will correspond
   !! To the planes of third coordinate from i3s to i3s+n3d-1. 
   !! The potential to the planes from i3s+i3xcsh to i3s+i3xcsh+n3p-1
   !! The array pot_ion to the planes from i3s+i3xcsh to i3s+i3xcsh+n3pi-1
   !! For global disposition i3s is equal to distributed case with i3xcsh=0.
   integer, intent(out) :: i3s

   !local variables
   !n(c) integer, parameter :: nordgr=4
   integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
   integer :: istart,iend,nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr
   integer :: n3pr1,n3pr2
   
   
   !calculate the dimensions wrt the geocode
   select case(geocode)
   case('P')
      call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
       if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   case('S')
      call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,igpu,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   case('F','H')
      call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,igpu,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2/2) then
           md2=(md2/nproc+1)*nproc 
        endif
      endif
   case('W')
      call W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,igpu,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   case default
      !write(*,*) geocode
      !stop 'PS_dim4allocation: geometry code not admitted'
      call f_err_throw('Geometry code "'//geocode//'" not admitted in PS_dim4allocation')
   end select
   
   !formal start and end of the slice
   istart=iproc*(md2/nproc)
   iend=min((iproc+1)*md2/nproc,m2)
   
   select case(datacode)
   case('D')
      call xc_dimensions(geocode,use_gradient,use_wb_corr,istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh)
   
      nwb=nxcl+nxc+nxcr-2
      nxt=nwbr+nwb+nwbl
   
      n3p=nxc
      n3d=nxt
      n3pi=n3p
   case('G')
      n3d=n03
      n3p=n03
      i3xcsh=0
      i3s=min(istart,m2-1)+1
      n3pi=max(iend-istart,0)
   case default
      !print *,datacode
      !stop 'PS_dim4allocation: data code not admitted'
      call f_err_throw('data code "'//datacode//'" not admitted in PS_dim4allocation')
   end select

!!  print *,'P4,iproc',iproc,'nxc,ncxl,ncxr,nwbl,nwbr',nxc,nxcl,nxcr,nwbl,nwbr,&
!!       'ixc,n3d,n3p,i3xcsh,i3s',ixc,n3d,n3p,i3xcsh,i3s

END SUBROUTINE PS_dim4allocation


!> Calculate the dimensions to be used for the XC part, taking into account also
!! the White-bird correction which should be made for some GGA functionals
!! @warning It is imperative that iend <=m2
subroutine xc_dimensions(geocode,use_gradient,use_wb_corr,&
     & istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh)
  implicit none

  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: use_gradient     !< .true. if functional is using the gradient.
  logical, intent(in) :: use_wb_corr      !< .true. if functional is using WB corrections.
  integer, intent(in) :: istart,iend      
  integer, intent(in) :: m2               !< dimension to be parallelised
  integer, intent(out) :: nxc             !< size of the parallelised XC potential
  integer, intent(out) :: nxcl,nxcr       !< left and right buffers for calculating the WB correction after call drivexc
  integer, intent(out) :: nwbl,nwbr       !< left and right buffers for calculating the gradient to pass to drivexc    
  integer, intent(out) :: i3s             !< starting addres of the distributed dimension
  integer, intent(out) :: i3xcsh          !< shift to be applied to i3s for having the striting address of the potential
  !local variables
  integer, parameter :: nordgr=4 !the order of the finite-difference gradient (fixed)

  if (istart <= m2-1) then
     nxc=iend-istart
     if (use_gradient .and. geocode == 'F') then
        if (.not. use_wb_corr) then
           !now the dimension of the part required for the gradient
           nwbl=min(istart,nordgr)
           nwbr=min(m2-iend,nordgr) !always m2 < iend
           nxcl=1
           nxcr=1
        else
           !now the dimension of the part required for the gradient
           if(istart<=nordgr) then
              nxcl=istart+1
              nwbl=0
           else
              nxcl=nordgr+1
              nwbl=min(nordgr,istart-nordgr)
           end if
           if(iend>=m2-nordgr+1) then
              nxcr=m2-iend+1
              nwbr=0
           else
              nxcr=nordgr+1
              nwbr=min(nordgr,m2-nordgr-iend)
           end if
        end if
     else if (geocode /= 'F' .and. use_gradient .and. nxc /= m2) then
        if (.not. use_wb_corr) then
           !now the dimension of the part required for the gradient
           nwbl=nordgr
           nwbr=nordgr
           nxcl=1
           nxcr=1
        else
           nxcl=nordgr+1
           nwbl=nordgr
           nxcr=nordgr+1
           nwbr=nordgr
        end if
     !this case is also considered below
     !else if (geocode /= 'F' .and. use_gradient .and. nxc == m2) then
     else 
        nwbl=0
        nwbr=0
        nxcl=1
        nxcr=1
     end if
     i3xcsh=nxcl+nwbl-1
     i3s=istart+1-i3xcsh
  else
     nwbl=0
     nwbr=0
     nxcl=1
     nxcr=1
     nxc=0
     i3xcsh=0
     i3s=m2
  end if
END SUBROUTINE xc_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the periodic system
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the other geometries of the Poisson Solver.
!!    The dimensions 2 and 3 are exchanged.
!! @author Luigi Genovese
!! @date October 2006
subroutine P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,enlarge_md2)
 implicit none
 !Arguments
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03   !< original real dimensions (input)
 integer, intent(in) :: nproc
 integer, intent(out) :: m1,m2,m3     !< original real dimension, with m2 and m3 exchanged
 integer, intent(out) :: n1,n2,n3     !< the first FFT dimensions (even for the moment - the medium point being n/2+1)
 integer, intent(out) :: md1,md2,md3  !< the n1,n2,n3 dimensions. They contain the real unpadded space.
                                      !! md2 is further enlarged to be a multiple of nproc
 integer, intent(out) :: nd1,nd2,nd3  !< Fourier dimensions for which the kernel is injective,
                                      !! formally 1/8 of the fourier grid. Here the dimension nd3 is
                                      !! enlarged to be a multiple of nproc
 !Local variables
 integer :: l1,l2,l3
 
 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 

 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3

 call fourier_dim(l1,n1)
 if (n1 /= m1) then
    call f_err_throw('The FFT in the x direction is not allowed, n01 dimension '//n01)
 end if

 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    call f_err_throw('The FFT in the z direction is not allowed, n03 dimension '//n03)
 end if
 
 call fourier_dim(l3,n3)
 if (n3 /= m3) then
    call f_err_throw('The FFT in the y direction is not allowed, n02 dimension '//n02)
 end if

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3

 !enlarge the md2 dimension to be compatible with MPI_ALLTOALL communication
 do while(nproc*(md2/nproc) < n2)
    md2=md2+1
 end do
 
 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1 
 nd2=n2/2+1
 nd3=n3/2+1
 
 !enlarge the md2 dimension to be compatible with MPI_ALLTOALL communication
 do while(modulo(nd3,nproc) /= 0)
    nd3=nd3+1
 end do

END SUBROUTINE P_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the surface system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                which has been properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!! @author Luigi Genovese
!! @date October 2006
subroutine S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2,non_ortho)
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc,gpu
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3
 logical, intent(in), optional :: non_ortho
 !local variables
 logical :: north

 north=.false.
 if (present(non_ortho)) north=non_ortho
 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
 else
  l3=2*m3
 endif

 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3

 call fourier_dim(l1,n1)
 if (n1 /= m1) then
    print *,'the FFT in the x direction is not allowed'
    print *,'n01 dimension',n01
    stop
 end if
 
 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    print *,'the FFT in the z direction is not allowed'
    print *,'n03 dimension',n03
    stop
 end if

 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2)
    !151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    !goto 151
    !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 !note: for non-orthorhombic cells the dimensions of the kernel 
 !should become double either in nd1 or in nd2

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1
 if (north) nd2=n2
 nd3=n3/2+1
 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !goto 250
    !endif
 end do

END SUBROUTINE S_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the Wires BC system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                which has been properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!! @author Luigi Genovese
!! @date October 2006
subroutine W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2)
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc,gpu
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
 else
  l3=2*m3
 endif

 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
       exit
    end if
    l1=l1+1
 end do

 
 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    print *,'the FFT in the z direction is not allowed'
    print *,'n03 dimension',n03
    stop
 end if

 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do

 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2)
    !151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    !goto 151
    !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1
 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !goto 250
    !endif
 end do

END SUBROUTINE W_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! zero-padded convolution
!!
!! @warning
!!    The dimension m2 and m3 correspond to n03 and n02 respectively
!!    this is needed since the convolution routine manage arrays of dimension
!!    (md1,md3,md2/nproc)
!! @author Luigi Genovese
!! @date February 2006
subroutine F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2)
 implicit none
 !Arguments
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03  !< Original real dimensions
 integer, intent(in) :: nproc,gpu
 integer, intent(out) :: n1,n2       !< The first FFT even dimensions greater that 2*m1, 2*m2
 integer, intent(out) :: n3          !< The double of the first FFT even dimension greater than m3
                                     !! (improved for the HalFFT procedure)
 integer, intent(out) :: md1,md2,md3 !< Half of n1,n2,n3 dimension. They contain the real unpadded space,
                                     !! which has been properly enlarged to be compatible with the FFT dimensions n_i.
                                     !! md2 is further enlarged to be a multiple of nproc
 integer, intent(out) :: nd1,nd2,nd3 !< Fourier dimensions for which the kernel FFT is injective,
                                     !! formally 1/8 of the fourier grid. Here the dimension nd3 is
                                     !! enlarged to be a multiple of nproc
 integer, intent(out) :: m1,m2,m3
 !Local variables
 integer :: l1,l2,l3, mul3

 !dimensions of the density in the real space, inverted for convenience
 m1=n01
 m2=n03
 m3=n02
 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=2*m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
  mul3=2
 else
  l3=2*m3
  mul3=4 ! in GPU we still need this dimension's size to be multiple of 4
 endif
 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3
 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
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
    if (modulo(n3,mul3) == 0) then
       exit
    end if
    l3=l3+1
 end do
 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2/2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2/2)
   !151 if (nproc*(md2/nproc).lt.n2/2) then
    md2=md2+1
   !goto 151
   !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1

 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !    goto 250
    ! endif
 end do

END SUBROUTINE F_FFT_dimensions
