!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
#include "global.h"

module floquet_oct_m
  use blas_oct_m
  use calc_mode_par_oct_m
  use comm_oct_m
  use density_oct_m
  use derivatives_oct_m
  use distributed_oct_m
  use eigensolver_oct_m
  use energy_calc_oct_m
  use excited_states_oct_m
  use forces_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_base_oct_m
  use hamiltonian_oct_m
  use ion_dynamics_oct_m
  use io_oct_m
  use iso_c_binding
  use kick_oct_m
  use kpoints_oct_m
  use lalg_adv_oct_m
  use lasers_oct_m
  use loct_math_oct_m
  use loct_oct_m
  use magnetic_oct_m
  use math_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mix_oct_m
  use mpi_lib_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use output_oct_m
  use parser_oct_m
  use partial_charges_oct_m
  use pert_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scf_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use sort_oct_m
  use species_oct_m
  use states_calc_oct_m
  use states_dim_oct_m
  use states_io_oct_m
  use states_oct_m
  use states_restart_oct_m
  use subspace_oct_m
  use system_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use write_iter_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                       &
       floquet_init,              &
       floquet_hamiltonians_init, &
       floquet_hamiltonian_update,&
       floquet_hamiltonian_solve, &
       floquet_init_dressed_wfs,  &
       floquet_hamiltonian_run_solver, &
       floquet_restart_dressed_st, &
       floquet_load_td_hamiltonians, &
       floquet_td_hamiltonians_sample, &
       floquet_calc_occupations,  &
       floquet_td_state,          &
       floquet_FBZ_filter,        &
!        zfloquet_FBZ_subspace_diag, &
!        dfloquet_FBZ_subspace_diag, &
       floquet_calc_FBZ_coefficients, & 
       floquet_init_td

  integer, public, parameter ::    &
       FLOQUET_NONE            = 0, &      
       FLOQUET_NON_INTERACTING = 1, &
       FLOQUET_FROZEN_PHONON   = 2, &
       FLOQUET_INTERACTING     = 3

contains
  
  !------------------------------------------------------------------------  
  subroutine floquet_init(sys,this,dim)
    type(system_t), intent(in)       :: sys
    type(floquet_t),    intent(out)  :: this
    integer,              intent(in) :: dim ! the standard dimension of the groundstate

    type(block_t)     :: blk
    integer :: ia, idir, idim, it, count

    FLOAT :: time_step

    PUSH_SUB(floquet_init)

    ! Default values
    this%floquet_apply = .false.
    this%sample = .false.
    this%FBZ_solver = .false.
    this%sample_one_only  = .true.
    this%lambda = M_ZERO


    !%Variable FloquetBoson
    !%Type flag
    !%Default classical
    !%Section Floquet
    !%Description
    !% Dont use it now!!
    !%Option classical_floquet 0
    !%
    !%Option qed_photon 1
    !% 
    !%Option qed_phonon 2
    !% 
    !%End
    call parse_variable('FloquetBoson', OPTION__FLOQUETBOSON__CLASSICAL_FLOQUET, this%boson)


    !%Variable FloquetMode
    !%Type flag
    !%Default non_interacting
    !%Section Floquet
    !%Description
    !% Types of Floquet analysis.
    !%Option none 0
    !%
    !%Option non_interacting 1
    !% 
    !%Option frozen_phonon 2
    !% 
    !%Option interacting 3
    !%
    !%End

    call parse_variable('FloquetMode', FLOQUET_NONE, this%mode)

    !%Variable FloquetTimePropagation
    !%Type logical
    !%Default no
    !%Section Floquet
    !%Description
    !% Perform time propagation of Floquet states. Makes sense only for FloquetBoson=qed_photon.
    !%End
    call parse_variable('FloquetTimePropagation', .false., this%propagate)



    if (this%mode == FLOQUET_NONE) then
      POP_SUB(floquet_init)
      return
    end if 
    
    ! Init mpi communicators
    this%is_parallel = multicomm_strategy_is_parallel(sys%mc, P_STRATEGY_OTHER)
    if(this%is_parallel) then
      call mpi_grp_init(this%mpi_grp, sys%mc%group_comm(P_STRATEGY_OTHER))
    else
      call mpi_grp_init(this%mpi_grp, -1)
    end if
    
    
    if (this%boson == OPTION__FLOQUETBOSON__CLASSICAL_FLOQUET) then
      call messages_print_stress(stdout, 'Floquet')
    else
      call messages_print_stress(stdout, 'Floquet QED')
    end if  
      

    call messages_print_var_option(stdout, 'FloquetBoson',  this%boson)

    call messages_print_var_option(stdout, 'FloquetMode',  this%mode)
    
    call messages_print_var_value(stdout,  'Allows time propagation',  this%propagate)
    

    if(this%mode==FLOQUET_FROZEN_PHONON) then
       if(.not.parse_is_defined('IonsTimeDependentDisplacements')) then
          write(message(1),'(a)') 'Please specify IonsTimeDependentDisplacements in order'
          write(message(2),'(a)') 'to use FloquetMode=frozen_phonon'
          call messages_fatal(2) 
       end if
    end if
    
    if (this%mode == FLOQUET_INTERACTING)  this%sample = .true.
    
    !%Variable FloquetInitialization
    !%Type flag
    !%Default gs
    !%Section Floquet
    !%Description
    !% Floquet diagonalization intial guess.
    !%Option f_gs 0
    !% Use ground state wavefunctions.
    !%Option f_rand  1
    !% Use random wavefunctions.
    !%End
    call parse_variable('FloquetInitialization', OPTION__FLOQUETINITIALIZATION__F_GS, this%init)
    call messages_print_var_option(stdout, 'FloquetInitialization',  this%init)
    
    !%Variable FloquetDimension
    !%Type integer
    !%Default -1
    !%Section Floquet
    !%Description
    !% Order of Floquet Hamiltonian. If negative number is given, downfolding
    !%is performed.
    !% If entered as a block allows to specify asymmetric dimensions 
    !%
    !% <tt>%FloquetDimension
    !% <br>&nbsp;&nbsp; mindim| maxdim 
    !% <br>%</tt>
    !%
    !%End
    this%order(:)=-1
    
    if(parse_block('FloquetDimension', blk) == 0) then
      if(parse_block_cols(blk,0) < 2) call messages_input_error('FloquetDimension')
      do idim = 1, 2
        call parse_block_integer(blk, 0, idim - 1, this%order(idim))
      end do
    else
      call parse_variable('FloquetDimension',-1,this%order(2))
      this%order(1)= -this%order(2)
    end if
    
    if(this%order(2) >= 0) then
      write(message(1),'(a,i5,a,i5,a)') 'Info: Floque-Hamiltonian tensor dimension [min,max]: [',&
                                         this%order(1),',', this%order(2),']'
      call messages_info(1)
      !Dimension of multiphoton Floquet-Hamiltonian
      this%floquet_dim = this%order(2)-this%order(1)+1

    else
      message(1) = 'Floquet-Hamiltonian downfolding not implemented for interacting propagation.'
      call messages_fatal(1)
      !this%downfolding = .true.
      !this%Forder = 1
      !this%Fdim = 3
    endif


    !%Variable FloquetModeCalcOccupations
    !%Type logical
    !%Default yes
    !%Section Floquet
    !%Description
    !% Calculate occupations of Floquet states.
    !%End
    call parse_variable('FloquetModeCalcOccupations', .true., this%calc_occupations)
    call messages_print_var_value(stdout,'Calculate occupations',  this%calc_occupations)

    !%Variable FloquetOccupationsCutDim
    !%Type integer
    !%Default 0
    !%Section Floquet
    !%Description
    !% Number of outermost subdimension with occupations set to zero
    !%End
    call parse_variable('FloquetOccupationsCutDim', 0, this%occ_cut)
    call messages_print_var_value(stdout,'Cut occupation of outer subdimensions',  this%occ_cut)

    !%Variable FloquetFrequency
    !%Type float
    !%Default 0
    !%Section Floquet
    !%Description
    !% Frequency for the Floquet analysis, this should be the carrier
    !%frequency or integer multiples of it.
    !% Other options will work, but likely be nonsense.
    !%
    !%End
    call parse_variable('FloquetFrequency', M_ZERO, this%omega, units_inp%energy)
    call messages_print_var_value(stdout,'Frequency used for Floquet analysis', this%omega)
    if(this%omega==M_ZERO .and. this%boson/=OPTION__FLOQUETBOSON__QED_PHOTON) then
      message(1) = "Please give a non-zero value for FloquetFrequency"
      call messages_fatal(1)
    endif

    ! get time of one cycle
    this%Tcycle=M_TWO*M_PI/this%omega

    !%Variable FloquetMaximumSolverIterations
    !%Type integer
    !%Default 35
    !%Section Floquet
    !%Description
    !% Maximumn Number of calls to eigensolver for solving the Floquet Hamiltonian
    !%
    !%End
    call parse_variable('FloquetMaximumSolverIterations ', 35 ,this%max_solve_iter)
    call messages_print_var_value(stdout,'Maximum eigensolver iterations', this%max_solve_iter)


    if (this%boson == OPTION__FLOQUETBOSON__CLASSICAL_FLOQUET) then

      !%Variable FloquetBZSolver
      !%Type logical
      !%Default no
      !%Section Floquet
      !%Description
      !% Use dedicated Floquet solver for the Floquet Brillouin zone
      !%End
      call parse_variable('FloquetBZSolver', .false., this%FBZ_solver)
      if(this%FBZ_solver) then
        message(1) = "Info: Using Floquet-BZ solver"
        call messages_info(1)
      endif

      !%Variable FloquetFBZStartDownfoldingScfSteps
      !%Type integer
      !%Default 200
      !%Section Floquet
      !%Description
      !% Expert option. 
      !%
      !%End    
      call parse_variable('FloquetFBZStartDownfoldingScfSteps', 200, this%cf_nsteps)
      call messages_print_var_value(stdout,'FloquetFBZStartDownfoldingScfSteps', this%cf_nsteps)


      !%Variable FloquetModeSampleOneCycle
      !%Type logical
      !%Default yes
      !%Section Floquet
      !%Description
      !% Stop sampling Floquet Hamiltoninans after the first cycle.
      !%End
      call parse_variable('FloquetModeSampleOneCycle', .true., this%sample_one_only)
      call messages_print_var_value(stdout,'FloquetModeSampleOneCycle',  this%sample_one_only)

      !%Variable FloquetSample
      !%Type integer
      !%Default FloquetDimension*3
      !%Section Floquet
      !%Description
      !% Number of points on which one Floquet cycle is sampled in the
      !%time-integral of the Floquet analysis.
      !%
      !%End
      this%nt = maxval(abs(this%order(:)))*3
      if (maxval(abs(this%order(:))) == 0 ) this%nt = 3 
    
      call parse_variable('FloquetSample', this%nt, this%nt)
      call messages_print_var_value(stdout,'Number of Floquet time-sampling points', this%nT)
      this%dt = this%Tcycle/real(this%nT)


      this%count = 1

      ! re-read time step from input
      call parse_variable('TDTimeStep', M_ZERO, time_step, unit = units_inp%time)
      if(time_step == M_ZERO) then
         message(1) = 'Did not find time-step in Floquet init, please give a value for TDTimeStep'
        call messages_fatal(1)
      end if
      this%interval = int(this%dt/time_step)
      this%ncycle = this%interval*this%nT


      call messages_print_var_value(stdout,'Steps in Floquet time-sampling interval',  this%interval)
      call messages_print_var_value(stdout,'Steps in Floquet time-sampling cycle',  this%ncycle)

    end if
    
    if (this%boson == OPTION__FLOQUETBOSON__QED_PHOTON) then
      
      this%mode = FLOQUET_NON_INTERACTING
      this%nt = 1
      
      !%Variable FloquetCavityLambda
      !%Type float
      !%Default 1
      !%Section Floquet
      !%Description
      !% Dont use it
      !%
      !%End
      call parse_variable('FloquetCavityLambda', M_ONE, this%lambda)
      call messages_print_var_value(stdout,'Cavity lambda', this%lambda)
      
      !%Variable FloquetCavityModePolarization
      !%Type float
      !%Default (1,0,0)
      !%Section Floquet
      !%Description
      !%
      !% <tt>%FloquetCavityModePolarization
      !% <br>&nbsp;&nbsp; x | y | z
      !% <br>%</tt>
      !%
      !%End
      this%pol(1:3)=(/M_z1,M_z0,M_z0/)
    
      if(parse_block('FloquetCavityModePolarization', blk) == 0) then
        if(parse_block_cols(blk,0) < 3) call messages_input_error('FloquetCavityModePolarization')
        do idim = 1, 3
          call parse_block_cmplx(blk, 0, idim - 1, this%pol(idim))
        end do
      end if
      
      this%pol(:)=this%pol(:)/sqrt(sum(abs(this%pol(:))**2))

      write(message(1),'(3x,a,3(a1,f7.4,a1,f7.4,a1))') 'Info: Cavity polarization direction: ', &
            '(', real(this%pol(1)), ',', aimag(this%pol(1)), '), ', &
            '(', real(this%pol(2)), ',', aimag(this%pol(2)), '), ', &
            '(', real(this%pol(3)), ',', aimag(this%pol(3)), ')'


!       write(message(1),'(a,e8.2,a,e8.2,a,e8.2,a,e8.2,a,e8.2,a,e8.2,a)') 'Info: Cavity polarization direction : (',&
!                                          this%pol(1),',', this%pol(2),',',this%pol(3),')'
      call messages_info(1)

      !Read the scf threshold for the relative density
      call parse_variable('ConvRelDens', CNST(1e-5), this%conv_rel_dens)
      
      
    end if

    this%spindim = dim


    ! distribute the stacked loop of the Floquet-Hamiltonian
    call distributed_init(this%flat_idx,this%nT*this%floquet_dim,sys%mc%group_comm(P_STRATEGY_OTHER),"Floquet-flat-index")
    SAFE_ALLOCATE(this%idx_map(this%nT*this%floquet_dim,2))
    count = 0
    do it=1,this%nT
       do idim=this%order(1),this%order(2)
          count = count + 1
          this%idx_map(count,1)=it
          this%idx_map(count,2)=idim
       end do
    end do
    
        
    call messages_print_stress(stdout)

   POP_SUB(floquet_init)

  end subroutine floquet_init


  subroutine floquet_init_td(sys, hm, dressed_st)
    type(system_t),      intent(in)    :: sys
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(inout) :: dressed_st

    type(grid_t),   pointer :: gr

    PUSH_SUB(floquet_init_td)

    gr => sys%gr

    if(hm%F%propagate) then
      
      ! Override undressed states initialization
      ! First: end states
      call states_deallocate_wfns(dressed_st)
      call states_end(dressed_st)
            
      ! Second: reinitialize the state object with the Floquet dimension
      dressed_st%floquet_dim = hm%F%floquet_dim
      dressed_st%floquet_FBZ = hm%F%FBZ_solver
      call states_init(dressed_st, gr, hm%geo)
      call kpoints_distribute(dressed_st%d,sys%mc)
      call states_distribute_nodes(dressed_st,sys%mc)
      call states_exec_init(dressed_st, sys%mc)
      call states_allocate_wfns(dressed_st,gr%der%mesh, wfs_type = TYPE_CMPLX)
      call states_densities_init(dressed_st, sys%gr, sys%geo)

            
      ! set dimension of Floquet Hamiltonian                                                    
      hm%d%dim = dressed_st%d%dim

      ! Switch on the Floquet Hamiltonian 
      hm%F%floquet_apply = .true.
        
    end if
    
    POP_SUB(floquet_init_td)
  end subroutine floquet_init_td



  subroutine floquet_hamiltonians_init(this, gr, st, sys)
    type(hamiltonian_t), intent(inout) :: this ! this is not great, as everyhting should be within the floquet_t
    type(grid_t),      intent(inout)   :: gr
    type(states_t),    intent(inout)   :: st !< at iter=0 this is the ggroundstate
    type(system_t),    intent(inout)   :: sys

    CMPLX, allocatable ::  temp_state1(:,:), temp_state2(:,:)
    FLOAT, allocatable :: eigenval(:), bands(:,:)
    character(len=80) :: filename
    integer :: it, ik, nst,ip, idim, ispin, ist, jst, ierr
    type(mesh_t) :: mesh
    type(states_t) :: gs_st
    FLOAT :: time_step, time
    type(scf_t) :: scf ! used for frozen_phonon
    integer :: ia, space_dim
    type(ion_dynamics_t) :: ions
    FLOAT, allocatable :: frozen_bands(:,:)

    CMPLX, allocatable :: psi_i(:,:), psi_j(:,:), hpsi1(:,:), hpsi2(:,:), gmat(:,:) 
    type(restart_t) :: restart
!     FLOAT  :: tmp(4,2)

    PUSH_SUB(floquet_hamiltonian_init)



    mesh = gr%der%mesh
    nst = st%nst

    ! the Hamiltonain gets assigned an array of td-Hamiltonians
    ! this is a bit recursive, so maybe there should be a Flqoeut moduel or something
    nullify(this%td_hm)
    SAFE_ALLOCATE(this%td_hm(1:this%F%nT))

    if(this%F%mode == FLOQUET_FROZEN_PHONON) then
       SAFE_ALLOCATE(frozen_bands(1:st%nst,1:gr%sb%kpoints%reduced%npoints))
       frozen_bands(1:st%nst,1:gr%sb%kpoints%reduced%npoints) = M_ZERO
       call ion_dynamics_init(ions,this%geo)
    end if

    ! initialize the instances of the Hamiltonians
    do it=1,this%F%nT
       this%td_hm(it)%F%Tcycle = this%F%Tcycle
       this%td_hm(it)%F%omega = this%F%omega
       this%td_hm(it)%F%dt = this%F%dt

       SAFE_ALLOCATE(this%td_hm(it)%geo)

       select case(this%F%mode)

       case(FLOQUET_FROZEN_PHONON) 
         time = it*this%F%dt
         space_dim = this%geo%space%dim
         
         call ion_dynamics_propagate(ions, gr%sb, this%geo, time, M_ZERO)
         call geometry_copy(this%td_hm(it)%geo, this%geo)
         this%td_hm(it)%geo%skip_species_pot_init = .true.

         call hamiltonian_init(this%td_hm(it), gr, this%td_hm(it)%geo, st, sys%ks%theory_level, &
                               sys%ks%xc_family,  family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
         call hamiltonian_epot_generate(this%td_hm(it), gr, this%td_hm(it)%geo, st, time=M_ZERO)

         call scf_init(scf,gr,this%td_hm(it)%geo,st,sys%mc,this%td_hm(it),sys%ks, conv_force = CNST(1e-8))
         call scf_run(scf,sys%mc,gr,this%td_hm(it)%geo,st,sys%ks,this%td_hm(it),sys%outp, gs_run=.false.)
         call scf_end(scf)

         frozen_bands(1:st%nst,1:gr%sb%kpoints%reduced%npoints) = &
              frozen_bands(1:st%nst,1:gr%sb%kpoints%reduced%npoints) + &
                     M_ONE/this%F%nT*st%eigenval(1:st%nst,1:gr%sb%kpoints%reduced%npoints)

! tmp(1:4,it) = st%eigenval(1:4,1)

         if(simul_box_is_periodic(gr%sb)) then
           write(filename,'(I5)') it
           filename = 'BO_bands_'//trim(adjustl(filename))
           call states_write_bandstructure(FLOQUET_DIR, st%nst, st, gr%sb, filename)
         end if

         call floquet_save_td_hamiltonian(this%td_hm(it), sys, it)

       case(FLOQUET_NON_INTERACTING)
         call geometry_copy(this%td_hm(it)%geo, this%geo)
          
         ! set flag to prevent species types to be touched, because
         ! hm%geo is only a pointer to the global geo instance
         this%td_hm(it)%geo%skip_species_pot_init = .true.

         call hamiltonian_init(this%td_hm(it), gr, this%td_hm(it)%geo, st, sys%ks%theory_level, &
                               sys%ks%xc_family,  family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
         time =this%F%Tcycle+ it*this%F%dt ! offset in time to catch switchon cycle
         do ispin=1,this%d%nspin
            forall (ip = 1:gr%mesh%np) this%td_hm(it)%vhxc(ip,ispin) = this%vhxc(ip, ispin)
         end do
         forall (ip = 1:gr%mesh%np) this%td_hm(it)%ep%vpsl(ip)= this%ep%vpsl(ip)
         call hamiltonian_epot_generate(this%td_hm(it), gr, this%td_hm(it)%geo, st, time=time)
         call hamiltonian_update(this%td_hm(it), gr%der%mesh, gr%der%boundaries,time=time)

        case(FLOQUET_INTERACTING)
           ! init is on the fly           
        end select

     enddo


     if(this%F%mode == FLOQUET_FROZEN_PHONON) then
        st%eigenval(1:st%nst,1:gr%sb%kpoints%reduced%npoints) = frozen_bands(1:st%nst,1:gr%sb%kpoints%reduced%npoints)
        if(simul_box_is_periodic(gr%sb)) then
          filename = 'frozen_bands'
          call states_write_bandstructure(FLOQUET_DIR, st%nst, st, gr%sb, filename)
        end if


        call states_init(gs_st, gr, sys%geo)
        call kpoints_distribute(gs_st%d,sys%mc)
        call states_distribute_nodes(gs_st,sys%mc)
        call states_exec_init(gs_st, sys%mc)
        call states_densities_init(gs_st, sys%gr, sys%geo)
        
        
        call states_allocate_wfns(gs_st, gr%mesh, wfs_type = TYPE_CMPLX)
      
        call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh, exact=.true.)

        if(ierr == 0) call states_load(restart, gs_st, gr, ierr, label = ": gs")
        if (ierr /= 0) then
          message(1) = 'Unable to read ground-state wavefunctions.'
          call messages_fatal(1)
        end if
        call restart_end(restart)
        
        ! Calculate the deformation potential matrix elements 
        SAFE_ALLOCATE(gmat(1:gs_st%nst, 1:gs_st%d%nik))
        SAFE_ALLOCATE(psi_i(1:gr%der%mesh%np_part,gs_st%d%dim)) 
        SAFE_ALLOCATE(psi_j(1:gr%der%mesh%np_part,gs_st%d%dim)) 
        SAFE_ALLOCATE(hpsi1(1:gr%der%mesh%np_part,gs_st%d%dim)) 
        SAFE_ALLOCATE(hpsi2(1:gr%der%mesh%np_part,gs_st%d%dim)) 

        do ist=gs_st%st_start,gs_st%st_end
          gmat = M_z0

          do ik=gs_st%d%kpt%start,gs_st%d%kpt%end

            call states_get_state(gs_st,gr%der%mesh,ist,ik,psi_i)

            call zhamiltonian_apply(this%td_hm(1), gr%der, psi_i, hpsi1, ist, ik)            
            call zhamiltonian_apply(this%td_hm(2), gr%der, psi_i, hpsi2, ist, ik)                        


            do jst=gs_st%st_start,gs_st%st_end
              call states_get_state(gs_st,gr%der%mesh,jst,ik,psi_j)
              
              gmat(jst,ik) = zmf_dotp(gr%mesh, gs_st%d%dim, psi_j, hpsi1)
              gmat(jst,ik) = gmat(jst,ik) - zmf_dotp(gr%mesh, st%d%dim, psi_j, hpsi2)
            
!               print *, ist, jst, gmat(jst,ik), tmp(jst,2)-tmp(jst,1), gmat(jst,ik)/(tmp(jst,2)-tmp(jst,1))
            end do
            
          enddo
          
          if(gs_st%parallel_in_states .or. gs_st%d%kpt%parallel) then
            call comm_allreduce(st%st_kpt_mpi_grp%comm, gmat(1:gs_st%nst, 1:gs_st%d%nik))
          end if
          
          if(simul_box_is_periodic(gr%sb)) then
            write(filename,'(a,i6)') 'Re_deformation_pot_i-', ist
            call states_write_bandstructure(FLOQUET_DIR, gs_st%nst, gs_st, gr%sb, filename, vec = real(gmat))
            write(filename,'(a,i6)') 'Im_deformation_pot_i-', ist
            call states_write_bandstructure(FLOQUET_DIR, gs_st%nst, gs_st, gr%sb, filename, vec = aimag(gmat))
          end if
        enddo          

        SAFE_DEALLOCATE_A(gmat)
        SAFE_DEALLOCATE_A(psi_i)
        SAFE_DEALLOCATE_A(psi_j)
        SAFE_DEALLOCATE_A(hpsi1)
        SAFE_DEALLOCATE_A(hpsi2)
        
        call states_deallocate_wfns(gs_st)
        
        call exit()
     end if


     POP_SUB(floquet_hamiltonian_init)
        
  end subroutine floquet_hamiltonians_init

  !--------------------------------------------
  subroutine floquet_td_hamiltonians_sample(hm,sys,iter)
  type(hamiltonian_t), intent(inout) :: hm
  type(system_t),      intent(in)    :: sys
  integer,               intent (in) :: iter

  integer :: it


  PUSH_SUB(floquet_td_hamiltonians_sample)   

  it = mod(iter/hm%F%interval,hm%F%nT)
  if(it==0) it=hm%F%nT

  call floquet_save_td_hamiltonian(hm, sys, it)

  POP_SUB(floquet_td_hamiltonians_sample)   

  end subroutine floquet_td_hamiltonians_sample

   !--------------------------------------------
   ! this is only called if F%mode=interacting
   subroutine floquet_hamiltonian_update(hm,st,gr,sys,iter)
     type(hamiltonian_t), intent(inout) :: hm
     type(states_t)      , intent(inout):: st
     type(grid_t),      intent(inout)   :: gr
     type(system_t),      intent(in)    :: sys
     integer,               intent (in) :: iter

     integer :: it

     integer :: ip, ispin
     FLOAT :: time

     PUSH_SUB(floquet_hamiltonian_update)


     it = mod(iter/hm%F%interval,hm%F%nT)
     if(it==0) it=hm%F%nT
     
     time = iter/hm%F%interval*hm%F%dt

     ! set the geometry of the td-hamiltonian
     call geometry_copy(hm%td_hm(it)%geo, hm%geo)

     ! set flag to prevent species types to be touched, because                                                   
     ! hm%geo is only a pointer to the global geo instance                                                        
     hm%td_hm(it)%geo%skip_species_pot_init = .true.

     call hamiltonian_init(hm%td_hm(it), gr, hm%td_hm(it)%geo, st, sys%ks%theory_level, &
                           sys%ks%xc_family, family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
     
     call hamiltonian_update(hm%td_hm(it), gr%der%mesh, gr%der%boundaries,time=time)

     ! fill time-dependent hamiltonian structure with scf-fields at this time
     do ispin=1,hm%d%nspin
       forall (ip = 1:gr%mesh%np) hm%td_hm(it)%vhxc(ip,ispin) = hm%vhxc(ip, ispin)
     end do
     forall (ip = 1:gr%mesh%np) hm%td_hm(it)%ep%vpsl(ip) = hm%ep%vpsl(ip)

     call floquet_save_td_hamiltonian(hm%td_hm(it), sys, it)


     call hamiltonian_epot_generate(hm%td_hm(it), gr, hm%td_hm(it)%geo, st, time=time)

     POP_SUB(floquet_hamiltonian_update)

   end subroutine floquet_hamiltonian_update

    !--------------------------------------------
    subroutine floquet_save_td_hamiltonian(hm, sys, it)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t),      intent(in)    :: sys
      integer,             intent(in)    :: it

      type(restart_t) :: restart   
      type(grid_t),   pointer :: gr
      integer :: ierr
      character(len=128) :: filename

      PUSH_SUB(floquet_save_td_hamiltonian)

      gr => sys%gr      
      
      call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_DUMP, sys%mc, ierr, gr%der%mesh)

      if (hm%theory_level /= INDEPENDENT_PARTICLES) then
        call hamiltonian_dump_vhxc(restart, hm, gr%mesh, ierr, idx = it)
      end if
        
      
      if (mpi_grp_is_root(mpi_world)) then      
        write(filename, fmt='(a,a,i6.6)') trim(restart_dir(restart)),'/geometry', it
        call geometry_write_xyz(hm%geo,trim(filename))
      end if


      call restart_end(restart)
      
      POP_SUB(floquet_save_td_hamiltonian)
    end subroutine floquet_save_td_hamiltonian

    !--------------------------------------------
    subroutine floquet_load_td_hamiltonians(hm, sys, ierr)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t),      intent(in)    :: sys
      integer,             intent(out)   :: ierr

      type(restart_t) :: restart   
      type(grid_t),   pointer :: gr
      type(states_t), pointer :: st        
      integer :: it
      character(len=128) :: filename
      FLOAT :: time

      PUSH_SUB(floquet_load_td_hamiltonians)

      gr => sys%gr      
      st => sys%st      
      
      call messages_print_stress(stdout, 'Load Floquet Hamitonians')
      call messages_print_stress(stdout)
      
      
      call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_LOAD, sys%mc, ierr, gr%der%mesh)
      
      if (ierr /=0 ) then
        return
        POP_SUB(floquet_load_td_hamiltonians)
      end if

      do it=1, hm%F%nT
        time = it * hm%F%dt

        write(message(1),'(a)') ' '
        write(message(2),'(a,i6,a,f9.3,a)') 'Info: Floquet td-hamiltonian: ', it,' (time =',time,')'
        write(message(3),'(a)') ' '
        call messages_info(3)        
        

        ! Read the geometry
        write(filename, fmt='(a,a,i6.6, a)') trim(restart_dir(restart)),'/geometry', it,'.xyz'
        call geometry_init(hm%td_hm(it)%geo, sys%space,  xyzfname = filename)
        
        
        call hamiltonian_init(hm%td_hm(it), gr, hm%td_hm(it)%geo, st, sys%ks%theory_level, &
                              sys%ks%xc_family,family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
 

        if (hm%theory_level /= INDEPENDENT_PARTICLES) then
          call hamiltonian_load_vhxc(restart, hm%td_hm(it), gr%mesh, ierr, idx = it)  
        end if

        if (ierr /=0 ) then
          return
          POP_SUB(floquet_load_td_hamiltonians)
        end if
        
        call hamiltonian_epot_generate(hm%td_hm(it), gr, hm%td_hm(it)%geo, st, time=time)
        
      end do

      call restart_end(restart)

      call messages_print_stress(stdout, 'Done')
      call messages_print_stress(stdout)
      
      POP_SUB(floquet_load_td_hamiltonians)
    end subroutine floquet_load_td_hamiltonians

    !--------------------------------------------
    subroutine floquet_hamiltonian_solve(hm,gr,sys,st, fromScratch)
      type(hamiltonian_t), intent(inout) :: hm
      type(grid_t),        intent(inout) :: gr
      type(system_t),      intent(inout) :: sys
      type(states_t),         intent(in) :: st
      logical, optional,      intent(in) :: fromScratch  

      type(states_t) :: dressed_st
      PUSH_SUB(floquet_hamiltonian_solve)

      call floquet_init_dressed_wfs(hm, sys, dressed_st, optional_default(fromScratch, .false.))  

      call floquet_hamiltonian_run_solver(hm, sys, dressed_st)
         
      call states_end(dressed_st)
      
      POP_SUB(floquet_hamiltonian_solve)
    end subroutine floquet_hamiltonian_solve

    !--------------------------------------------
    subroutine floquet_restart_dressed_st(hm, sys, dressed_st, ierr, fromScratch)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t),         intent(in) :: sys
      type(states_t),      intent(inout) :: dressed_st
      integer,               intent(out) :: ierr
      logical, optional,      intent(in) :: fromScratch  
      
      integer :: nik, dim, nst, iter
      type(grid_t),   pointer :: gr
      type(restart_t) :: restart        
      
      
      PUSH_SUB(floquet_restart_dressed_st)
      
      gr => sys%gr
      
      ! solver iteration
      iter = 0
      
      call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_LOAD, &
                       sys%mc, ierr, gr%der%mesh)

      if(ierr == 0 .and. .not. optional_default(fromScratch, .false.)) then 
        call states_look(restart, nik, dim, nst, ierr)
        if(dim==dressed_st%d%dim .and. nik==gr%sb%kpoints%reduced%npoints &
            .and. nst==dressed_st%nst) then

          call states_load(restart, dressed_st, gr, ierr, iter)
        else
          write(message(1),'(a)') 'Floquet restart structure not commensurate.'
          call messages_warning(1)
        end if
      end if
      
      call restart_end(restart)
      
      hm%F%iter = iter
      
      
      POP_SUB(floquet_restart_dressed_st)
    end subroutine floquet_restart_dressed_st


    !--------------------------------------------
    subroutine floquet_init_dressed_wfs(hm, sys, dressed_st, fromScratch)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t),      intent(inout) :: sys
      type(states_t),        intent(out) :: dressed_st
      logical,                intent(in) :: fromScratch
      
      type(restart_t) :: restart
      type(grid_t),   pointer :: gr
      type(states_t), pointer :: st
      CMPLX, allocatable :: temp_state1(:,:), temp_state2(:,:)
      integer            :: iter , ik, in, im, ist, idim, ierr, nik, dim, nst
      
      
      PUSH_SUB(floquet_init_dressed_wfs)
      
      gr => sys%gr
      st => sys%st
      
      ! initialize a state object with the Floquet dimension
      dressed_st%floquet_dim = hm%F%floquet_dim
      dressed_st%floquet_FBZ = hm%F%FBZ_solver
      call states_init(dressed_st, gr, hm%geo)
      call kpoints_distribute(dressed_st%d,sys%mc)
      call states_distribute_nodes(dressed_st,sys%mc)
      call states_exec_init(dressed_st, sys%mc)
      call states_allocate_wfns(dressed_st,gr%der%mesh, wfs_type = TYPE_CMPLX)
      call states_densities_init(dressed_st, sys%gr, sys%geo)

      call floquet_restart_dressed_st(hm, sys, dressed_st, ierr, fromScratch)
      
      if(ierr == 0 .and. .not. fromScratch) then
        continue
      else
        
        ierr = 0 
        if (hm%F%init == OPTION__FLOQUETINITIALIZATION__F_GS ) then
          call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)

          if(ierr == 0) call states_load(restart, sys%st, sys%gr, ierr, label = ": gs")
          if (ierr /= 0) then
            message(1) = 'Unable to read ground-state wavefunctions.'
            message(2) = 'Floquet states will be randomly initialized.'
            call messages_warning(2)
          end if
          call restart_end(restart)
        end if

        ! only use gs states for init, if they are not distributed 
        if(.not. st%parallel_in_states .and. ierr == 0 &
           .and.  hm%F%init == OPTION__FLOQUETINITIALIZATION__F_GS ) then
          message(1) ='Info: Initialize Floquet states with gs wavefunctions'
          call messages_info(1)
          ! initialize floquet states from scratch
          SAFE_ALLOCATE(temp_state1(1:gr%der%mesh%np,st%d%dim))
          SAFE_ALLOCATE(temp_state2(1:gr%der%mesh%np,hm%F%floquet_dim*st%d%dim))
          
          if(hm%F%FBZ_solver) then
            do ik=st%d%kpt%start,st%d%kpt%end
              do ist=st%st_start,st%st_end
                call states_get_state(st,gr%der%mesh,ist,ik,temp_state1)
                temp_state2(:,:) = M_ZERO
                do idim=1,st%d%dim
                  temp_state2(1:gr%der%mesh%np,(-hm%F%order(1)+1)*st%d%dim+idim) = temp_state1(1:gr%der%mesh%np,idim)
                end do
                call states_set_state(dressed_st,gr%der%mesh, ist, ik,temp_state2)
              enddo
            enddo
          else
             do ik=st%d%kpt%start,st%d%kpt%end
              do in=1,hm%F%floquet_dim
                do ist=st%st_start,st%st_end
                  call states_get_state(st,gr%der%mesh,ist,ik,temp_state1)
                  temp_state2(:,:) = M_ZERO
                  do idim=1,st%d%dim
                    temp_state2(1:gr%der%mesh%np,(in-1)*st%d%dim+idim) = temp_state1(1:gr%der%mesh%np,idim)
                  end do
                  call states_set_state(dressed_st,gr%der%mesh, (in-1)*st%nst+ist, ik,temp_state2)
                enddo
              enddo
            enddo
          end if
          SAFE_DEALLOCATE_A(temp_state1)
          SAFE_DEALLOCATE_A(temp_state2)
        else
          message(1) ='Info: Initialize Floquet states with random wavefunctions'
          call messages_info(1)
          ! randomize
          SAFE_ALLOCATE(temp_state1(1:gr%der%mesh%np,st%d%dim))
          SAFE_ALLOCATE(temp_state2(1:gr%der%mesh%np,hm%F%floquet_dim*st%d%dim))
          do ik=st%d%kpt%start,st%d%kpt%end
            do ist=dressed_st%st_start,dressed_st%st_end
                
              do in=1,hm%F%floquet_dim*st%d%dim  
                if(dressed_st%randomization == PAR_INDEPENDENT) then
                  call zmf_random(gr%der%mesh, temp_state2(:,in), gr%der%mesh%vp%xlocal-1, &
                             seed=mpi_world%rank, normalized= .true.)
                else
                  call zmf_random(gr%der%mesh, temp_state2(:, in),seed=mpi_world%rank,normalized = .true.)
                end if
              enddo
                
              call states_set_state(dressed_st,gr%der%mesh, ist, ik,temp_state2)
            end do
          enddo
          SAFE_DEALLOCATE_A(temp_state1)
          SAFE_DEALLOCATE_A(temp_state2)
        end if

    end if

      
      POP_SUB(floquet_init_dressed_wfs)
      
    end subroutine floquet_init_dressed_wfs


    !--------------------------------------------
    subroutine floquet_hamiltonian_run_solver(hm, sys, dressed_st)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t), intent(inout)      :: sys
      type(states_t), intent(inout)      :: dressed_st
        

      logical :: converged, converged_elec, have_gs
      integer :: iter , maxiter, ik, in, im, ist, idim, ierr, nik, dim, nst, is
      type(eigensolver_t) :: eigens, eigens_elec
      type(restart_t) :: restart
      character(len=80) :: filename
      character(len=80) :: iterstr
      character(len=40) :: msg
      integer :: file
      integer :: iunit, ii, iatom, idir, iter_elec
      real(8) :: itime, etime, rel_dens

      type(grid_t),   pointer :: gr
      type(states_t), pointer :: st
      type(states_t)          :: FBZ_st
      type(mix_t)             :: smix
        
      FLOAT, allocatable :: spect(:,:), me(:,:), rhoin(:,:,:), rhoout(:,:,:), mixrho(:,:,:), dressed_rhoold(:,:,:), tmp(:)
 
      PUSH_SUB(floquet_hamiltonian_run_solver)
 
      gr => sys%gr
      st => sys%st

      call messages_print_stress(stdout, 'Floquet diagonalization')

      iter = hm%F%iter

      if(.not.hm%F%FBZ_solver) hm%F%floquet_apply = .true.

      ! set dimension of Floquet Hamiltonian                                                    
      hm%d%dim = dressed_st%d%dim
      
      call eigensolver_init(eigens, gr, dressed_st, sys%ks%xc)
      eigens%converged(:) = 0

      ! the subspace diagonalization is hardcoded (and only applied once)
!      if(hm%F%FBZ_solver) eigens%sdiag%method = OPTION__SUBSPACEDIAGONALIZATION__FLOQUET_SS
      if(hm%F%FBZ_solver) eigens%sdiag%method = OPTION__SUBSPACEDIAGONALIZATION__NONE

      have_gs= .false.
      if (hm%F%calc_occupations) then
        message(1) = 'Info: Read ground-state wavefunctions to calculate occupations.'
        call messages_info(1)
        
        call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)

        if(ierr == 0) call states_load(restart, sys%st, sys%gr, ierr, label = ": gs")
        if (ierr /= 0) then
          message(1) = 'Unable to read ground-state wavefunctions.'
          message(2) = 'Occupations will not be calculated.'
          call messages_warning(2)
        else  
          have_gs = .true.
        end if
        call restart_end(restart)
      end if

      call mix_init(smix, gr%fine%der, gr%fine%mesh%np, 1, st%d%nspin)
      SAFE_ALLOCATE(rhoin(gr%mesh%np,1:1,st%d%nspin))
      SAFE_ALLOCATE(rhoout(gr%mesh%np,1:1,st%d%nspin))
      SAFE_ALLOCATE(mixrho(gr%mesh%np,1:1,st%d%nspin))
      SAFE_ALLOCATE(dressed_rhoold(gr%mesh%np,1:1,st%d%nspin))
      SAFE_ALLOCATE(tmp(1:gr%fine%mesh%np))

      converged=.false.
      maxiter = hm%F%max_solve_iter
      itime = loct_clock()
      dressed_rhoold = M_ZERO
      rel_dens = M_HUGE

      do while(.not.converged .and. iter <= maxiter) 
         if (hm%F%boson==OPTION__FLOQUETBOSON__QED_PHOTON  .and. rel_dens<hm%F%conv_rel_dens) then
            exit
         end if

         write(msg,'(a,i5)') 'Iter #', iter         
         call messages_print_stress(stdout, trim(msg))
        
         if(hm%F%FBZ_solver) then
            call floquet_FBZ_eigensolver_run(eigens, sys, dressed_st, hm, iter,converged)
         else
           call eigensolver_run(eigens, gr, dressed_st, hm, 1,converged)
         end if

         
         
         if (hm%F%calc_occupations) then
           call floquet_calc_occupations(hm, sys, dressed_st, st)
           call density_calc(dressed_st, sys%gr,dressed_st%rho)

           rhoin(1:gr%fine%mesh%np, 1, 1:st%d%nspin) = st%rho(1:gr%fine%mesh%np, 1:st%d%nspin)
           rhoout(1:gr%fine%mesh%np, 1, 1:st%d%nspin) = dressed_st%rho(1:gr%fine%mesh%np, 1:st%d%nspin)

           if (hm%F%boson==OPTION__FLOQUETBOSON__QED_PHOTON) then
              !Convergence check on the relative density
              rel_dens = M_ZERO
              do is = 1, st%d%nspin
                tmp = abs(rhoout(1:gr%fine%mesh%np, 1, is) - dressed_rhoold(1:gr%fine%mesh%np, 1, is))
                rel_dens = rel_dens + dmf_integrate(gr%fine%mesh, tmp)
              end do
              rel_dens = rel_dens / st%qtot
              write(message(1),'(a,es15.8,4a)') 'Relative Density: ', rel_dens, ' [au]'
              call messages_info(1)
              dressed_rhoold(1:gr%fine%mesh%np, 1, 1:st%d%nspin) = rhoout(1:gr%fine%mesh%np, 1, 1:st%d%nspin)
   
              call dmixing(smix, rhoin, rhoout, mixrho)
              st%rho(1:gr%fine%mesh%np, 1:st%d%nspin) = mixrho(1:gr%fine%mesh%np, 1, 1:st%d%nspin)
              dressed_st%rho(1:gr%fine%mesh%np, 1:st%d%nspin) = mixrho(1:gr%fine%mesh%np, 1, 1:st%d%nspin)
           end if

           call v_ks_calc(sys%ks, hm, dressed_st, sys%geo)
           call energy_calc_total(hm, gr, dressed_st, iunit = 0)
           write(message(1),'(a,es15.8,4a)') 'Energy tot: ', units_from_atomic(units_out%energy, hm%energy%total) &
                                                           , ' [',  trim(units_abbrev(units_out%energy)), ']'
           call messages_info(1)
         end if 


         call smear_find_fermi_energy(dressed_st%smear, dressed_st%eigenval, dressed_st%occ, dressed_st%qtot, &
           dressed_st%d%nik, dressed_st%nst, dressed_st%d%kweights)
               
         ! compute Floquet-spin expectation values
         if(dressed_st%d%ispin == SPINORS) then
            call floquet_calc_spin(hm%F,gr%mesh,dressed_st)
         end if

         write(message(1),'(a,i6)') 'Converged eigenvectors: ', sum(eigens%converged(1:st%d%nik))
         call messages_info(1)
         call states_write_eigenvalues(stdout, dressed_st%nst, dressed_st, gr%sb, eigens%diff, &
                                       compact = .true.)         
         
         etime = loct_clock() - itime
         itime = etime + itime
         write(message(1),'(a,i5,a,f14.2)') 'Elapsed time for iter # ', iter,':', etime
         call messages_info(1)

         
         write(iterstr,'(I5)') iter !hm%F_count
         if (simul_box_is_periodic(gr%sb) .and. kpoints_have_zero_weight_path(gr%sb%kpoints)) then
           filename = 'bands'!//trim(adjustl(iterstr))

           call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, gr%sb, filename)
           
           if (hm%F%calc_occupations .and. have_gs) then                                 
             filename = 'bands_occ'!//trim(adjustl(iterstr))
             call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, gr%sb, filename, vec = dressed_st%occ)
           end if

           if(dressed_st%d%ispin == SPINORS) then
             filename = 'bands_Sx'
             call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, gr%sb, &
                                             filename, vec = dressed_st%spin(1,:,:))
             filename = 'bands_Sy'
             call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, gr%sb, &
                                             filename, vec = dressed_st%spin(2,:,:))
             filename = 'bands_Sz'
             call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, gr%sb, &
                                             filename, vec = dressed_st%spin(3,:,:))
           end if
            
         else
                      
           filename = FLOQUET_DIR//'/eigenvalues'!//trim(adjustl(iterstr))
           
           call states_write_eigenvalues(filename, dressed_st%nst, dressed_st, gr%sb, eigens%diff)
         end if
         
         call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_DUMP, &
                           sys%mc, ierr, gr%der%mesh)
         call states_dump(restart, dressed_st, gr, ierr, iter)
         call states_dump_rho(restart,dressed_st, gr, ierr, iter=iter)
         call restart_end(restart)

!         ! Ground State Update for QED
!         if (hm%F%boson==OPTION__FLOQUETBOSON__QED_PHOTON) then
!            !switch off floquet hamiltonian                                                           
!            hm%F%floquet_apply = .false.
!            hm%d%dim = st%d%dim
!
!            call eigensolver_init(eigens_elec, gr, st)
!            eigens_elec%converged(:) = 0
!
!            converged_elec = .false.
!            iter_elec = 0 
!            do while(.not.converged_elec.and.iter_elec <= maxiter)
!
!               call eigensolver_run(eigens_elec, gr, st, hm, 1, converged_elec)
!
!               iter_elec = iter_elec +1 
!
!            enddo
!            
!            call states_fermi(st, gr%mesh)
!
!            hm%d%dim = dressed_st%d%dim
!            hm%F%floquet_apply = .true.
!            call eigensolver_end(eigens_elec)
!         endif


         iter = iter +1
      end do

      call eigensolver_end(eigens)
      call mix_end(smix)

      ! filter the FBZ
      if(.not.hm%F%FBZ_solver .and. hm%F%boson == OPTION__FLOQUETBOSON__CLASSICAL_FLOQUET) then
         ! this triggers FBZ allocation 
         hm%F%FBZ_solver = .true.
         call floquet_init_dressed_wfs(hm, sys, FBZ_st, fromScratch=.true.)
         ! reset flag
         hm%F%FBZ_solver = .false.
         call floquet_FBZ_filter(sys, hm, dressed_st, FBZ_st)
         if (simul_box_is_periodic(sys%gr%sb) .and. kpoints_have_zero_weight_path(sys%gr%sb%kpoints)) then
           filename = 'FBZ_filtered_bands'
           call states_write_bandstructure(FLOQUET_DIR, FBZ_st%nst, FBZ_st, sys%gr%sb, filename, vec = FBZ_st%occ)
         else
           filename = trim(FLOQUET_DIR)//'FBZ_filtered_eigenvalues'
           call states_write_eigenvalues(filename, FBZ_st%nst, FBZ_st, sys%gr%sb)
         end if
         call states_end(FBZ_st)
      end if

      if(mpi_grp_is_root(mpi_world)) then
        iunit = io_open(trim(FLOQUET_DIR) // "/info" , action='write')
      else
        iunit = 0
      end if
      
      call energy_calc_total(hm, gr, dressed_st, iunit, full = .true.)
      call forces_calculate(gr, sys%geo, hm, dressed_st, sys%ks)

      if(mpi_grp_is_root(mpi_world)) then
        write(iunit, '(3a)') 'Energy [', trim(units_abbrev(units_out%energy)), ']:'
        write(iunit,'(1x)')
        

        write(iunit,'(3a)') 'Floquet-forces on the ions [', trim(units_abbrev(units_out%force)), "]"
        write(iunit,'(a,10x,99(14x,a))') ' Ion', (index2axis(idir), idir = 1, gr%sb%dim)
        do iatom = 1, sys%geo%natoms
          write(iunit,'(i4,a10,10f15.6)') iatom, trim(species_label(sys%geo%atom(iatom)%species)), &
               (units_from_atomic(units_out%force, sys%geo%atom(iatom)%f(idir)), idir=1, gr%sb%dim)
        end do
        write(iunit,'(1x,100a1)') ("-", ii = 1, 13 + gr%sb%dim * 15)
        write(iunit,'(a14, 10f15.6)') " Max abs force", &
            (units_from_atomic(units_out%force, maxval(abs(sys%geo%atom(1:sys%geo%natoms)%f(idir)))), idir=1, gr%sb%dim)
        write(iunit,'(a14, 10f15.6)') " Total force", &
            (units_from_atomic(units_out%force, sum(sys%geo%atom(1:sys%geo%natoms)%f(idir))), idir=1, gr%sb%dim)

        call io_close(iunit)        
      end if
      
      call output_all(sys%outp, sys%gr, sys%geo, dressed_st, hm, sys%ks, trim(FLOQUET_DIR))
      
      !switch off floquet hamiltonian                                                           
      hm%F%floquet_apply = .false.                                                               
  
      ! reset dimension
      hm%d%dim = hm%F%spindim
      
      hm%F%count=hm%F%count + 1
      
      SAFE_DEALLOCATE_A(spect)
      SAFE_DEALLOCATE_A(me)
      SAFE_DEALLOCATE_A(rhoin)
      SAFE_DEALLOCATE_A(rhoout)
      SAFE_DEALLOCATE_A(mixrho)
      SAFE_DEALLOCATE_A(dressed_rhoold)
      SAFE_DEALLOCATE_A(tmp)
      
      POP_SUB(floquet_hamiltonian_run_solver)

    end subroutine floquet_hamiltonian_run_solver

    !--------------------------------
    subroutine floquet_td_state(F,mesh,F_psi,F_eval,time,psi_t)
  
    type(floquet_t) :: F
    type(mesh_t)    :: mesh
    CMPLX   :: F_psi(:,:)
    FLOAT   :: F_eval
    FLOAT   :: time
    CMPLX   :: psi_t(:,:)
  
    integer :: idim, im, imm
  
    psi_t = M_z0
    do im= F%order(1),F%order(2)
      imm = im - F%order(1) + 1
      do idim=1,F%spindim
        psi_t(1:mesh%np,idim)  =  psi_t(1:mesh%np,idim) + &
                                  exp(M_zI*im*F%omega*time)*F_psi(1:mesh%np,(imm-1)*F%spindim+idim)
      end do
    end do
    psi_t(1:mesh%np,:) = exp(-M_zI*F_eval*time)* psi_t(1:mesh%np,:)

    ! normalize td-state
!     psi_t(1:mesh%np,1:F%spindim)  = psi_t(1:mesh%np,1:F%spindim)/zmf_nrm2(mesh,F%spindim,psi_t)

    end subroutine floquet_td_state


    !--------------------------------
    subroutine floquet_calc_spin(F,mesh,st)

    type(floquet_t) :: F
    type(mesh_t)    :: mesh
    type(states_t)  :: st
    CMPLX, allocatable   :: psi(:,:), temp(:,:)

    integer :: ik, ist,idim, im, imm

    PUSH_SUB(floquet_calc_spin)

    SAFE_ALLOCATE(psi(1:mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(temp(1:mesh%np,1:F%spindim))
    
    st%spin = M_ZERO

    do ik=st%d%kpt%start,st%d%kpt%end
      do ist=st%st_start,st%st_end
        call states_get_state(st,mesh,ist,ik,psi)

        do im= F%order(1),F%order(2)
          imm = im - F%order(1) + 1
          temp = M_ZERO
          do idim=1,F%spindim
            temp(1:mesh%np, idim) = psi(1:mesh%np,(imm-1)*F%spindim+idim)
          end do
          st%spin(1:3,ist,ik) = st%spin(1:3,ist,ik) + state_spin(mesh,temp)
        end do
      end do
    end do

    if(st%d%kpt%parallel) then
       call comm_allreduce(st%d%kpt%mpi_grp%comm,  st%spin(:,:,:))
    end if

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(temp)

    POP_SUB(floquet_calc_spin)
   end subroutine floquet_calc_spin


    !--------------------------------------------
    ! Occupation strategy wrapper routine 
    !--------------------------------------------
    subroutine floquet_calc_occupations(hm, sys, dressed_st, st)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t), intent(inout)      :: sys
      type(states_t), intent(inout)      :: dressed_st
      type(states_t), intent(in)         :: st

      PUSH_SUB(floquet_calc_occupations)


      select case (hm%F%boson)
      case (OPTION__FLOQUETBOSON__CLASSICAL_FLOQUET)

        if(hm%F%FBZ_solver) then
          call floquet_calc_FBZ_coefficients(hm, sys, dressed_st, st, M_ZERO)
        else
          call floquet_calc_occupations_sudden(hm, sys, dressed_st)
        end if

      case (OPTION__FLOQUETBOSON__QED_PHOTON)
        call floquet_calc_occupations_sudden(hm, sys, dressed_st)

      case (OPTION__FLOQUETBOSON__QED_PHONON)

      end select
      
      POP_SUB(floquet_calc_occupations)
    end subroutine floquet_calc_occupations

    !--------------------------------------------
    subroutine floquet_calc_occupations_norms(hm, sys, st)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t), intent(inout)      :: sys
      type(states_t), intent(inout)      :: st

      integer :: Fdim, spindim, ik, ist, im, idim
      integer, allocatable :: idx(:,:)
      FLOAT, allocatable :: norms(:), evals(:,:), occs(:,:)
      CMPLX, allocatable :: state(:,:), state_reshape(:,:,:)

      PUSH_SUB(floquet_calc_occupations_norms)

      ASSERT(st%parallel_in_states .eqv. .false.)

      Fdim = hm%F%floquet_dim
      spindim = hm%F%spindim
      SAFE_ALLOCATE(norms(1:st%nst))
      SAFE_ALLOCATE(idx(1:st%nst,1:st%d%nik))
      SAFE_ALLOCATE(evals(1:st%nst,1:st%d%nik))
      SAFE_ALLOCATE(occs(1:st%nst,1:st%d%nik))
      SAFE_ALLOCATE(state(1:sys%gr%mesh%np,1:st%d%dim))
      SAFE_ALLOCATE(state_reshape(1:sys%gr%mesh%np,1:spindim,1:Fdim))
      
      idx(:,:) = 0
      
      do ik=st%d%kpt%start,st%d%kpt%end
         norms = M_ZERO
         do ist=1,st%nst
           call states_get_state(st,sys%gr%mesh, ist, ik, state)
           state_reshape = reshape(state,(/sys%gr%mesh%np,spindim,Fdim/))
           do idim=1,spindim
             norms(ist) = norms(ist) + zmf_nrm2(sys%gr%mesh, state_reshape(:,idim,1))
           end do
         end do
         call sort(norms, idx(:,ik))
         
         do ist=1,st%nst
           evals(ist,ik) = st%eigenval(idx(ist,ik),ik)
         end do
       end do
      
       if(st%d%kpt%parallel) then
         call comm_allreduce(st%d%kpt%mpi_grp%comm,  evals(:,:))
       end if
       
      
       call smear_find_fermi_energy(st%smear, evals, occs, st%qtot, st%d%nik, st%nst, st%d%kweights)

       call smear_fill_occupations(st%smear, evals, occs, st%d%nik, st%nst)

       do ik=st%d%kpt%start,st%d%kpt%end
         do ist=1,st%nst
           st%occ(idx(ist,ik),ik) = occs(ist,ik)
         end do
       end do

       SAFE_DEALLOCATE_A(norms)
       SAFE_DEALLOCATE_A(idx)
       SAFE_DEALLOCATE_A(evals)
       SAFE_DEALLOCATE_A(occs)
       SAFE_DEALLOCATE_A(state)
       SAFE_DEALLOCATE_A(state_reshape)
      
      POP_SUB(floquet_calc_occupations_norms)
    end subroutine floquet_calc_occupations_norms


    !--------------------------------------------
    subroutine floquet_calc_occupations_sudden(hm, sys, dressed_st)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t), intent(inout)      :: sys
      type(states_t), intent(inout)      :: dressed_st

      CMPLX, allocatable :: u_ma(:,:),  psi(:,:), tmp(:)
      FLOAT :: omega, dt, sum_dr, sum_gs
      FLOAT, allocatable :: occs_on_subdim(:)
      integer :: idx, im, it, ist, idim, nT, ik, ia, Fdim(2), imm
      
      type(states_t), pointer :: gs_st
      type(mesh_t),   pointer :: mesh


      PUSH_SUB(floquet_calc_occupations_sudden)

      mesh  => sys%gr%der%mesh
      gs_st => sys%st
      
      dt=hm%F%dt
      nT=hm%F%nT
      omega=hm%F%omega
      Fdim(:)=hm%F%order(:)


      SAFE_ALLOCATE(psi(1:mesh%np,gs_st%d%dim))
      SAFE_ALLOCATE(u_ma(1:mesh%np,hm%F%floquet_dim*hm%F%spindim))
      SAFE_ALLOCATE(tmp(gs_st%d%dim))
      SAFE_ALLOCATE(occs_on_subdim(hm%F%floquet_dim))

      dressed_st%occ(:,:) = M_ZERO

      do ik=dressed_st%d%kpt%start,dressed_st%d%kpt%end

        do ia=dressed_st%st_start,dressed_st%st_end
          call states_get_state(dressed_st, mesh, ia, ik, u_ma)
          dressed_st%occ(ia,ik) = M_ZERO 

          do ist=gs_st%st_start,gs_st%st_end
            call states_get_state(gs_st, mesh, ist, ik, psi)
            
            tmp(:) = M_ZERO
            if (hm%F%boson==OPTION__FLOQUETBOSON__QED_PHOTON) then 

               do idim=1,gs_st%d%dim
                 tmp(idim)  =  tmp(idim) + zmf_dotp(mesh, psi(1:mesh%np,idim),u_ma(1:mesh%np,idim))
               end do

            else

               do idx=hm%F%flat_idx%start, hm%F%flat_idx%end
                  it = hm%F%idx_map(idx,1)
                  im = hm%F%idx_map(idx,2)
                  imm = im - Fdim(1) + 1
                  do idim=1,gs_st%d%dim
                    tmp(idim)  =  tmp(idim) + exp(-M_zI*im*omega*it*dt)/nT * zmf_dotp(mesh, psi(1:mesh%np,idim), &
                                                                             u_ma(1:mesh%np,(imm-1)*gs_st%d%dim+idim))
                  end do
                  
               end do

            end if

            if(hm%F%is_parallel) call comm_allreduce(hm%F%mpi_grp%comm, tmp(:))
            
            dressed_st%occ(ia,ik) = dressed_st%occ(ia,ik) + abs(sum(tmp(:)))**2 * gs_st%occ(ist,ik)
!             print *, ist, "gs_st%occ(ist,:) =", gs_st%occ(ist,:)
!             print *, tmp(:)

          enddo
          ! cut out occupations depending on sub-space dimensions
          if(hm%F%occ_cut > 0) then
            occs_on_subdim = M_ZERO
            do imm=1,hm%F%floquet_dim
              do idim=1,gs_st%d%dim
                occs_on_subdim(imm) = occs_on_subdim(imm) + zmf_dotp(mesh,u_ma(1:mesh%np,(imm-1)*gs_st%d%dim+idim), &
                                                                          u_ma(1:mesh%np,(imm-1)*gs_st%d%dim+idim))
              end do
            end do
            print *, maxloc(occs_on_subdim,dim=1)
            if(maxloc(occs_on_subdim,dim=1) <= hm%F%occ_cut .or. &
              maxloc(occs_on_subdim,dim=1) >= hm%F%floquet_dim-hm%F%occ_cut)  dressed_st%occ(ia,ik) = M_ZERO
         end if

        enddo
        
      enddo
      
      
      if(dressed_st%parallel_in_states .or. dressed_st%d%kpt%parallel) then
        call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  dressed_st%occ(:,:))
      end if

! variant using reshape, probably wonrg, but woudl be faster
!      SAFE_ALLOCATE(occs_on_subdim(1:gs_st%nst,hm%F%floquet_dim,1:dressed_st%d%kpt%nglobal))
!      occs_on_subdim = reshape(dressed_st%occ(:,:),(/gs_st%nst,hm%F%floquet_dim,dressed_st%d%kpt%nglobal/))
!      occs_on_subdim(1:gs_st%nst,1:hm%F%occ_cut,1:dressed_st%d%kpt%nglobal) = M_ZERO
!      occs_on_subdim(1:gs_st%nst,hm%F%floquet_dim-hm%F%occ_cut:hm%F%floquet_dim,1:dressed_st%d%kpt%nglobal) = M_ZERO
!      dressed_st%occ(:,:) = reshape(occs_on_subdim,(/dressed_st%nst,dressed_st%d%kpt%nglobal/))
      SAFE_DEALLOCATE_A(occs_on_subdim)

      ! occupations checksum 
      do ik=1, dressed_st%d%kpt%nglobal
        sum_dr = sum(dressed_st%occ(:,ik))
        sum_gs = sum(gs_st%occ(:,ik))
        if (mpi_grp_is_root(mpi_world)) then
          if( abs(sum_dr/sum_gs - M_ONE) > 1E-5) then
            call messages_write('Occupations checksum failed for kpoint = ')
            call messages_write(ik, fmt = '(i6)')
            call messages_write(':   gs_occ =  ')
            call messages_write(sum_gs, fmt ='(f12.6)')
            call messages_write('   floquet_occ =  ')
            call messages_write(sum_dr, fmt ='(f12.6)')
            call messages_warning()
          end if
        end if
        dressed_st%occ(:,ik) = dressed_st%occ(:,ik)*sum_gs/sum_dr
      end do
      
      
      SAFE_DEALLOCATE_A(tmp)
      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(u_ma)
    
      POP_SUB(floquet_calc_occupations_sudden)
    
    end subroutine floquet_calc_occupations_sudden


    !--------------------------------------------
    subroutine floquet_calc_FBZ_coefficients(hm, sys, FBZ_st, ref_st, time)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t), intent(inout)      :: sys
      type(states_t), intent(inout)      :: FBZ_st
      type(states_t), intent(in)         :: ref_st
      FLOAT,          intent(in)         :: time            

      CMPLX, allocatable ::  psi(:,:), u_m(:,:), matba(:,:)
      CMPLX, allocatable ::  Fpsi_a(:,:), Fpsi_b(:,:), rhs(:,:), sol(:,:)    
      integer :: ib, ia , ik, ist, idim            
      type(mesh_t),   pointer :: mesh

      PUSH_SUB(floquet_calc_FBZ_coefficients)

      ASSERT(FBZ_st%nst==ref_st%nst)

      mesh  => sys%gr%der%mesh
      
      SAFE_ALLOCATE(matba(1:FBZ_st%nst,1:FBZ_st%nst))
      SAFE_ALLOCATE(rhs(1,1:FBZ_st%nst))
      SAFE_ALLOCATE(sol(1,1:FBZ_st%nst))

      SAFE_ALLOCATE(psi(1:mesh%np,ref_st%d%dim))
      SAFE_ALLOCATE(Fpsi_a(1:mesh%np,FBZ_st%d%dim))
      SAFE_ALLOCATE(Fpsi_b(1:mesh%np,FBZ_st%d%dim))
      SAFE_ALLOCATE(u_m(1:mesh%np,hm%F%floquet_dim*hm%F%spindim))
      
      if (.not. associated(FBZ_st%coeff)) then
        SAFE_ALLOCATE(FBZ_st%coeff(1:FBZ_st%nst,FBZ_st%d%nik))
      end if

      FBZ_st%coeff(:,:) = M_z0

!       do ik=FBZ_st%d%kpt%start,FBZ_st%d%kpt%end
!         rhs(:,:) = M_z0
!         matba(:,:) = M_z0
!
!         do ib=FBZ_st%st_start,FBZ_st%st_end
!           !print *,'test info', ik, ib, FBZ_st%d%kpt%start, FBZ_st%d%kpt%end,FBZ_st%st_start,FBZ_st%st_end
!           call states_get_state(FBZ_st, mesh, ib, ik, u_m)
!           call floquet_td_state(hm%F,mesh,u_m,FBZ_st%eigenval(ib,ik),time,Fpsi_b)
!
!           do ia=1,FBZ_st%nst
!             call states_get_state(FBZ_st, mesh, ia, ik, u_m)
!             call floquet_td_state(hm%F,mesh,u_m,FBZ_st%eigenval(ia,ik),time,Fpsi_a)
!             do idim=1,hm%F%spindim
!               matba(ib,ia) = matba(ib,ia) + zmf_dotp(mesh, Fpsi_b(1:mesh%np,idim), Fpsi_a(1:mesh%np,idim))
!             end do
!
!           end do !ia
!
!           do ist=1,ref_st%nst
!             call states_get_state(ref_st, mesh, ist, ik, psi)
!             do idim=1,hm%F%spindim
!               rhs(1,ib) = rhs(1,ib) + ref_st%occ(ist,ik)*zmf_dotp(mesh, Fpsi_b(1:mesh%np,idim), psi(1:mesh%np,idim))
!             end do
!
!           end do! ist
!         end do !ib
!         if(FBZ_st%parallel_in_states)  then
!           call comm_allreduce(FBZ_st%mpi_grp%comm, rhs)
!           call comm_allreduce(FBZ_st%mpi_grp%comm, matba)
!         end if
!
!
!         ! solve coefficient system
!         call lalg_linsyssolve(FBZ_st%nst, 1, matba, rhs, sol)
!         FBZ_st%coeff(:,ik) = sol(1,:)
!
!       end do! ik

      !TD Floquet states are orthogonal 
      ! \phi_i(t)= \sum_a f_a^i\psi_a(t)
      ! fa=sum_i f_a^i
      do ik=FBZ_st%d%kpt%start,FBZ_st%d%kpt%end
  
        do ia=FBZ_st%st_start,FBZ_st%st_end
          if (FBZ_st%eigenval(ia,ik) == M_HUGE) cycle ! skip garbage
          !print *,'test info', ik, ia, FBZ_st%d%kpt%start, FBZ_st%d%kpt%end,FBZ_st%st_start,FBZ_st%st_end, FBZ_st%eigenval(ia,ik)
          call states_get_state(FBZ_st, mesh, ia, ik, u_m)
          call floquet_td_state(hm%F,mesh,u_m,FBZ_st%eigenval(ia,ik),time,Fpsi_a)
    
          do ist=1,ref_st%nst
            call states_get_state(ref_st, mesh, ist, ik, psi)
            do idim=1,hm%F%spindim
              FBZ_st%coeff(ia,ik) = FBZ_st%coeff(ia,ik) + ref_st%occ(ist,ik)*zmf_dotp(mesh, Fpsi_a(1:mesh%np,idim), psi(1:mesh%np,idim))
            end do           

          end do! ist
        end do !ia
  
      end do! ik



      
      
      if(FBZ_st%parallel_in_states .or. FBZ_st%d%kpt%parallel) then
        call comm_allreduce(FBZ_st%st_kpt_mpi_grp%comm,  FBZ_st%coeff)
      end if

      FBZ_st%occ(:,:) =  abs(FBZ_st%coeff(:,:))
      
      
      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(u_m)
      SAFE_DEALLOCATE_A(Fpsi_a)
      SAFE_DEALLOCATE_A(Fpsi_b)
      
      SAFE_DEALLOCATE_A(matba)
      SAFE_DEALLOCATE_A(rhs)
      SAFE_DEALLOCATE_A(sol)
    
      POP_SUB(floquet_calc_FBZ_coefficients)    
    end subroutine floquet_calc_FBZ_coefficients


    
  subroutine floquet_FBZ_eigensolver_run(eigens, sys, st, hm, iter, conv)
   type(eigensolver_t),  intent(inout) :: eigens
   type(system_t),         intent(in)  :: sys
   type(states_t),       intent(inout) :: st
   type(hamiltonian_t),  intent(inout) :: hm
   integer,              intent(in)    :: iter
   logical,    optional, intent(out)   :: conv

   type(states_t)        :: sub_st
   type(mesh_t), pointer :: mesh
   integer               :: ik, ist, idim, im
   FLOAT                 :: diff
   CMPLX, allocatable    :: full_state(:,:), sub_state(:,:), full_state2(:,:)
   character(len=80) :: filename

   PUSH_SUB(floquet_FBZ_eigensolver_run)

   ASSERT(st%nst == sys%st%nst )

   mesh => sys%gr%der%mesh

   eigens%folded_spectrum = .true.

   if(iter ==0) then
     SAFE_ALLOCATE(eigens%spectrum_shift(1:st%nst,st%d%nik))

     st%eigenval = M_ZERO
     do ik=st%d%kpt%start,st%d%kpt%end
       call zfloquet_FBZ_subspace_diag(sys%gr%der, st, hm, ik, start=.true.)
     end do
     call comm_allreduce(st%st_kpt_mpi_grp%comm,  st%eigenval(:,:))

     if (simul_box_is_periodic(sys%gr%sb) .and. kpoints_have_zero_weight_path(sys%gr%sb%kpoints)) then
       filename = 'FBZ_start_bands'
       call states_write_bandstructure(FLOQUET_DIR, st%nst, st, sys%gr%sb, filename)
     else                    
       filename = FLOQUET_DIR//'/FBZ_start_eigenvalues'
       call states_write_eigenvalues(filename, st%nst, st, sys%gr%sb)
     end if

     eigens%spectrum_shift(1:st%nst,1:st%d%nik) = st%eigenval(1:st%nst,1:st%d%nik)

!SAFE_ALLOCATE(full_state(1:mesh%np,1:st%d%dim))             
!SAFE_ALLOCATE(full_state2(1:mesh%np,1:st%d%dim))
!do ik=st%d%kpt%start,st%d%kpt%end 
!  do ist=1,st%nst
!    call states_get_state(st,sys%gr%der%mesh, ist, ik, full_state)
!    call zhamiltonian_apply(hm, sys%gr%der, full_state, full_state2, ist, ik)
!  !  print *, ik, ist, zstates_residue(sys%gr%mesh, st%d%dim, full_state2, st%eigenval(ist, ik),full_state)
!  print *, st%eigenval(ist, ik),  zmf_dotp(sys%gr%mesh, st%d%dim, full_state, full_state2, reduce = .false.)
!  end do
!end do

   end if

   hm%F%floquet_apply = .true.
   call eigensolver_run(eigens, sys%gr, st, hm, iter,conv)

   ! put here floquet-dimension analysis for convergence etc. ..

   hm%F%floquet_apply = .false.
   POP_SUB(floquet_FBZ_eigensolver_run)

  end subroutine floquet_FBZ_eigensolver_run


  subroutine floquet_FBZ_filter(sys, hm, st_full, st_filtered)
    type(system_t),    intent(in)    :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    type(states_t),     intent(in) :: st_full
    type(states_t),     intent(inout):: st_filtered
  
    integer :: ik,ist,jst, Fdim, spindim, count, idim, im
    FLOAT, allocatable :: norms(:,:)
    CMPLX, allocatable :: state(:,:), state_reshape(:,:,:), state0(:,:)

    PUSH_SUB(floquet_FBZ_filter)

    ASSERT(st_full%d%kpt%start == st_filtered%d%kpt%start)
    ASSERT(st_full%d%kpt%end == st_filtered%d%kpt%end)
    ASSERT(st_full%d%dim == st_filtered%d%dim)

    Fdim = hm%F%floquet_dim
    spindim = hm%F%spindim
    SAFE_ALLOCATE(norms(1:st_full%nst,1:Fdim))
    SAFE_ALLOCATE(state(1:sys%gr%mesh%np,1:st_full%d%dim))
    SAFE_ALLOCATE(state_reshape(1:sys%gr%mesh%np,1:spindim,1:Fdim))
    st_filtered%eigenval(:,:) = M_HUGE !set to some nonsense value to rise flags

    do ik=st_full%d%kpt%start,st_full%d%kpt%end
       norms = M_ZERO
       do ist=1,st_full%nst
          call states_get_state(st_full,sys%gr%mesh, ist, ik, state)
          state_reshape = reshape(state,(/sys%gr%mesh%np,spindim,Fdim/))
          do im=1,Fdim
            do idim=1,spindim
              norms(ist,im) = norms(ist,im)+ zmf_nrm2(sys%gr%mesh, state_reshape(:,idim,im))
            end do
          end do
       end do
     
       count = 0
       do ist=1,st_full%nst
          if(maxloc(norms(ist,:),dim=1) == Fdim/2+1) then
            count = count + 1
            if(count .le. st_filtered%nst ) then
              call states_get_state(st_full    ,sys%gr%mesh,   ist, ik, state)
              call states_set_state(st_filtered,sys%gr%mesh, count, ik, state)
              st_filtered%eigenval(count,ik) = st_full%eigenval(ist,ik)
            else
              call messages_write('FBZ filtering is ambigous for kpoint = ')
              call messages_write(ik, fmt = '(i6)')
              call messages_warning()
            end if
          end if
       end do

    end do

    ! verify the FBZ contains only unique states by checking sectorwise orthogonality 
    SAFE_ALLOCATE(state0(1:sys%gr%mesh%np,1:spindim))
    do ik=st_filtered%d%kpt%start,st_filtered%d%kpt%end
      do ist=1,st_filtered%nst
        call states_get_state(st_filtered,sys%gr%mesh, ist, ik, state)
        state_reshape = reshape(state,(/sys%gr%mesh%np,spindim,Fdim/))
        state0(1:sys%gr%mesh%np,1:spindim) =  state_reshape(1:sys%gr%mesh%np,1:spindim,Fdim/2+1)
        do jst=1,ist-1
           call states_get_state(st_filtered,sys%gr%mesh, jst, ik, state)
           state_reshape = reshape(state,(/sys%gr%mesh%np,spindim,Fdim/))
           do im=1,Fdim
              do idim=1,spindim
                 ! this should be zero, but under no circumstance can it be close to the norm of state0
                 if(abs(zmf_dotp(sys%gr%mesh, state0(:,idim),state_reshape(:,idim,im))) .gt. M_HALF*zmf_nrm2(sys%gr%mesh, state0(:,idim))) then
                    message(1) = 'FBZ filtering failed'
                    call messages_warning(1)
                 end if
              end do
           end do
         end do
       end do
     end do

    SAFE_DEALLOCATE_A(norms)
    SAFE_DEALLOCATE_A(state)
    SAFE_DEALLOCATE_A(state0)
    SAFE_DEALLOCATE_A(state_reshape)

    POP_SUB(floquet_FBZ_filter)
  end subroutine floquet_FBZ_filter
        

end module floquet_oct_m
