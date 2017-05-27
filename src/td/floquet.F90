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
  use iso_c_binding
  use calc_mode_par_oct_m
  use comm_oct_m
  use distributed_oct_m
  use eigensolver_oct_m
  use excited_states_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use hamiltonian_base_oct_m
  use output_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use kick_oct_m
  use kpoints_oct_m
  use lasers_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
  use loct_math_oct_m
  use magnetic_oct_m
  use math_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use multicomm_oct_m
  use parser_oct_m
  use partial_charges_oct_m
  use pert_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scf_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_dim_oct_m
  use states_io_oct_m
  use states_restart_oct_m
  use system_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
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
       floquet_td_hamiltonians_sample

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
    


    !%Variable TDFloquetMode
    !%Type flag
    !%Default non_interacting
    !%Section Floquet
    !%Description
    !% Types of Floquet analysis performed when TDOutput=td_floquet
    !%Option none 0
    !%
    !%Option non_interacting 1
    !% 
    !%Option frozen_phonon 2
    !% 
    !%Option interacting 3
    !%
    !%End

    call parse_variable('TDFloquetMode', FLOQUET_NONE, this%mode)


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
    
    

    call messages_print_stress(stdout, 'Floquet')
    call messages_print_var_option(stdout, 'TDFloquetMode',  this%mode)

    if(this%mode==FLOQUET_FROZEN_PHONON) then
       if(.not.parse_is_defined('IonsTimeDependentDisplacements')) then
          write(message(1),'(a)') 'Please specify IonsTimeDependentDisplacements in order'
          write(message(2),'(a)') 'to use TDFloquetMode=frozen_phonon'
          call messages_fatal(2) 
       end if
    end if
    
    if (this%mode == FLOQUET_INTERACTING)  this%sample = .true.
    
    !%Variable TDFloquetInitialization
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
    call parse_variable('TDFloquetInitialization', OPTION__TDFLOQUETINITIALIZATION__F_GS, this%init)
    call messages_print_var_option(stdout, 'TDFloquetInitialization',  this%init)
    
    !%Variable TDFloquetModeSampleOneCycle
    !%Type logical
    !%Default yes
    !%Section Floquet
    !%Description
    !% Stop sampling Floquet Hamiltoninans after the first cycle.
    !%End
    call parse_variable('TDFloquetModeSampleOneCycle', .true., this%sample_one_only)
    call messages_print_var_value(stdout,'TDFloquetModeSampleOneCycle',  this%sample_one_only)


    !%Variable TDFloquetModeCalcOccupations
    !%Type logical
    !%Default yes
    !%Section Floquet
    !%Description
    !% Calculate occupations of Floquet states.
    !%End
    call parse_variable('TDFloquetModeCalcOccupations', .true., this%calc_occupations)
    call messages_print_var_value(stdout,'Calculate occupations',  this%calc_occupations)


    !%Variable TDFloquetFrequency
    !%Type float
    !%Default 0
    !%Section Floquet
    !%Description
    !% Frequency for the Floquet analysis, this should be the carrier
    !%frequency or integer multiples of it.
    !% Other options will work, but likely be nonsense.
    !%
    !%End
    call parse_variable('TDFloquetFrequency', M_ZERO, this%omega, units_inp%energy)
    call messages_print_var_value(stdout,'Frequency used for Floquet analysis', this%omega)
    if(this%omega==M_ZERO) then
      message(1) = "Please give a non-zero value for TDFloquetFrequency"
      call messages_fatal(1)
    endif

    ! get time of one cycle
    this%Tcycle=M_TWO*M_PI/this%omega


    !%Variable TDFloquetMaximumSolverIterations
    !%Type integer
    !%Default 35
    !%Section Floquet
    !%Description
    !% Maximumn Number of calls to eigensolver for solving the Floquet Hamiltonian
    !%
    !%End
    call parse_variable('TDFloquetMaximumSolverIterations ', 35 ,this%max_solve_iter)
    call messages_print_var_value(stdout,'Maximum eigensolver iterations', this%max_solve_iter)


    !%Variable TDFloquetDimension
    !%Type integer
    !%Default -1
    !%Section Floquet
    !%Description
    !% Order of Floquet Hamiltonian. If negative number is given, downfolding
    !%is performed.
    !% If entered as a block allows to specify asymmetric dimensions 
    !%
    !% <tt>%TDFloquetDimension
    !% <br>&nbsp;&nbsp; mindim| maxdim 
    !% <br>%</tt>
    !%
    !%End
    this%order(:)=-1
    
    if(parse_block('TDFloquetDimension', blk) == 0) then
      if(parse_block_cols(blk,0) < 2) call messages_input_error('TDFloquetDimension')
      do idim = 1, 2
        call parse_block_integer(blk, 0, idim - 1, this%order(idim))
      end do
    else
      call parse_variable('TDFloquetDimension',-1,this%order(2))
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

    !%Variable TDFloquetSample
    !%Type integer
    !%Default TDFloquetDimension*3
    !%Section Floquet
    !%Description
    !% Number of points on which one Floquet cycle is sampled in the
    !%time-integral of the Floquet analysis.
    !%
    !%End
    call parse_variable('TDFloquetSample',maxval(abs(this%order(:)))*3 ,this%nt)
    call messages_print_var_value(stdout,'Number of Floquet time-sampling points', this%nT)
    this%dt = this%Tcycle/real(this%nT)


    this%count = 1
    this%spindim = dim

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

  subroutine floquet_hamiltonians_init(this, gr, st, sys)
    type(hamiltonian_t), intent(inout) :: this ! this is not great, as everyhting should be within the floquet_t
    type(grid_t),      intent(inout)   :: gr
    type(states_t),    intent(inout)   :: st !< at iter=0 this is the ggroundstate
    type(system_t),    intent(inout)   :: sys

    CMPLX, allocatable ::  temp_state1(:,:), temp_state2(:,:)
    FLOAT, allocatable :: eigenval(:), bands(:,:)
    character(len=80) :: filename
    integer :: it, ik, nst,ip, idim, ispin
    type(mesh_t) :: mesh
    type(states_t) :: hm_st
    FLOAT :: time_step, time
    type(scf_t) :: scf ! used for frozen_phonon
    integer :: ia, space_dim
    type(ion_dynamics_t) :: ions
    FLOAT, allocatable :: frozen_bands(:,:)


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
                               sys%ks%xc_family,sys%ks%xc_flags,  family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
         call hamiltonian_epot_generate(this%td_hm(it), gr, this%td_hm(it)%geo, st, time=M_ZERO)

         call scf_init(scf,gr,this%td_hm(it)%geo,st,this%td_hm(it))
         call scf_run(scf,sys%mc,gr,this%td_hm(it)%geo,st,sys%ks,this%td_hm(it),sys%outp, gs_run=.false.)
         call scf_end(scf)

         frozen_bands(1:st%nst,1:gr%sb%kpoints%reduced%npoints) = &
              frozen_bands(1:st%nst,1:gr%sb%kpoints%reduced%npoints) + &
                     M_ONE/this%F%nT*st%eigenval(1:st%nst,1:gr%sb%kpoints%reduced%npoints)

         write(filename,'(I5)') it
         filename = 'BO_bands_'//trim(adjustl(filename))
         call states_write_bandstructure(FLOQUET_DIR, st%nst, st, gr%sb, filename)

         call floquet_save_td_hamiltonian(this%td_hm(it), sys, it)

       case(FLOQUET_NON_INTERACTING)
         call geometry_copy(this%td_hm(it)%geo, this%geo)
          
         ! set flag to prevent species types to be touched, because
         ! hm%geo is only a pointer to the global geo instance
         this%td_hm(it)%geo%skip_species_pot_init = .true.

         call hamiltonian_init(this%td_hm(it), gr, this%td_hm(it)%geo, st, sys%ks%theory_level, &
                               sys%ks%xc_family,sys%ks%xc_flags,  family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
         time =this%F%Tcycle+ it*this%F%dt ! offset in time to catch switchon cycle
         do ispin=1,this%d%nspin
            forall (ip = 1:gr%mesh%np) this%td_hm(it)%vhxc(ip,ispin) = this%vhxc(ip, ispin)
         end do
         forall (ip = 1:gr%mesh%np) this%td_hm(it)%ep%vpsl(ip)= this%ep%vpsl(ip)
         call hamiltonian_epot_generate(this%td_hm(it), gr, this%td_hm(it)%geo, st, time=time)
         call hamiltonian_update(this%td_hm(it), gr%der%mesh,time=time)

        case(FLOQUET_INTERACTING)
           ! init is on the fly           
        end select

     enddo


     if(this%F%mode == FLOQUET_FROZEN_PHONON) then
        st%eigenval(1:st%nst,1:gr%sb%kpoints%reduced%npoints) = frozen_bands(1:st%nst,1:gr%sb%kpoints%reduced%npoints)
        filename = 'frozen_bands'
        call states_write_bandstructure(FLOQUET_DIR, st%nst, st, gr%sb, filename)
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
                           sys%ks%xc_family,sys%ks%xc_flags, family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
     
     call hamiltonian_update(hm%td_hm(it), gr%der%mesh,time=time)

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
                              sys%ks%xc_family,sys%ks%xc_flags, family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
 

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
      call states_init(dressed_st, gr, hm%geo,floquet_dim=hm%F%floquet_dim)
      call kpoints_distribute(dressed_st%d,sys%mc)
      call states_distribute_nodes(dressed_st,sys%mc)
      call states_exec_init(dressed_st, sys%mc)
      call states_allocate_wfns(dressed_st,gr%der%mesh, wfs_type = TYPE_CMPLX)

  
      call floquet_restart_dressed_st(hm, sys, dressed_st, ierr, fromScratch)
      
      if(ierr == 0 .and. .not. fromScratch) then
         call states_berry_connection(FLOQUET_DIR,'floquet_berry_connection',dressed_st, gr,gr%sb)
      else
        
        ierr = 0 
        if (hm%F%init == OPTION__TDFLOQUETINITIALIZATION__F_GS ) then
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
           .and.  hm%F%init == OPTION__TDFLOQUETINITIALIZATION__F_GS) then
           message(1) ='Info: Initialize Floquet states with gs wavefunctions'
           call messages_info(1)
            ! initialize floquet states from scratch
            SAFE_ALLOCATE(temp_state1(1:gr%der%mesh%np,st%d%dim))
            SAFE_ALLOCATE(temp_state2(1:gr%der%mesh%np,hm%F%floquet_dim*st%d%dim))
       
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
        

      logical :: converged, have_gs
      integer :: iter , maxiter, ik, in, im, ist, idim, ierr, nik, dim, nst
      type(eigensolver_t) :: eigens
      type(restart_t) :: restart
      character(len=80) :: filename
      character(len=80) :: iterstr
      character(len=40) :: msg
      integer :: file
      real(8) :: itime, etime

      type(grid_t),   pointer :: gr
      type(states_t), pointer :: st
        
      FLOAT, allocatable :: spect(:,:), me(:,:)  
 
      PUSH_SUB(floquet_hamiltonian_run_solver)
 
      gr => sys%gr
      st => sys%st


      call messages_print_stress(stdout, 'Floquet diagonalization')

      iter = hm%F%iter

      hm%F%floquet_apply = .true.
      ! set dimension of Floquet Hamiltonian                                                    
      hm%d%dim = dressed_st%d%dim
      
      call eigensolver_init(eigens, gr, dressed_st)
      eigens%converged(:) = 0

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


      converged=.false.
      maxiter = hm%F%max_solve_iter
      itime = loct_clock()
      do while(.not.converged.and.iter <= maxiter)
         write(msg,'(a,i5)') 'Iter #', iter         
         call messages_print_stress(stdout, trim(msg))
        
         call eigensolver_run(eigens, gr, dressed_st, hm, 1,converged)
         
         if (hm%F%calc_occupations .and. have_gs) then
           call floquet_calc_occupations(hm, sys, dressed_st)
         end if 
         
         call smear_find_fermi_energy(dressed_st%smear, dressed_st%eigenval, dressed_st%occ, dressed_st%qtot, &
           dressed_st%d%nik, dressed_st%nst, dressed_st%d%kweights)
                 

         write(message(1),'(a,i6)') 'Converged eigenvectors: ', sum(eigens%converged(1:st%d%nik))
         call messages_info(1)
         call states_write_eigenvalues(stdout, dressed_st%nst, dressed_st, gr%sb, eigens%diff, &
                                       compact = .true.)         
         
         etime = loct_clock() - itime
         itime = etime + itime
         write(message(1),'(a,i5,a,f14.2)') 'Elapsed time for iter # ', iter,':', etime
         call messages_info(1)

         
         write(iterstr,'(I5)') iter !hm%F_count
         if (simul_box_is_periodic(gr%sb)) then
           filename = 'floquet_multibands_'//trim(adjustl(iterstr))

           call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, gr%sb, filename)
           
           if (hm%F%calc_occupations .and. have_gs) then                                 
             filename = 'floquet_multibands_occ_'//trim(adjustl(iterstr))
             call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, gr%sb, filename, vec = dressed_st%occ)
           end if
                     
         else
                      
           filename = FLOQUET_DIR//'/floquet_eigenvalues_'//trim(adjustl(iterstr))
           call states_write_eigenvalues(filename, dressed_st%nst, dressed_st, gr%sb, eigens%diff)
         end if
         
         call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_DUMP, &
                           sys%mc, ierr, gr%der%mesh)
         call states_dump(restart, dressed_st, gr, ierr, iter)
         call restart_end(restart)


         iter = iter +1
      end do

      call eigensolver_end(eigens)
      !
      
      !switch off floquet hamiltonian                                                           
      hm%F%floquet_apply = .false.                                                               
      ! reset dimension
      hm%d%dim = hm%F%spindim
      
      hm%F%count=hm%F%count + 1
      
      SAFE_DEALLOCATE_A(spect)
      SAFE_DEALLOCATE_A(me)
      
      POP_SUB(floquet_hamiltonian_run_solver)

    end subroutine floquet_hamiltonian_run_solver


    !--------------------------------------------
    subroutine floquet_calc_occupations(hm, sys, dressed_st)
      type(hamiltonian_t), intent(inout) :: hm
      type(system_t), intent(inout)      :: sys
      type(states_t), intent(inout)      :: dressed_st

      CMPLX, allocatable :: u_ma(:,:),  psi(:,:), tmp(:)
      FLOAT :: omega, dt, sum_dr, sum_gs   
      integer :: idx, im, it, ist, idim, nT, ik, ia, Fdim(2), imm
      
      type(states_t), pointer :: gs_st
      type(mesh_t),   pointer :: mesh


      PUSH_SUB(floquet_calc_occupations)

      mesh  => sys%gr%der%mesh
      gs_st => sys%st
      
      dt=hm%F%dt
      nT=hm%F%nT
      omega=hm%F%omega
      Fdim(:)=hm%F%order(:)


      SAFE_ALLOCATE(psi(1:mesh%np,gs_st%d%dim))
      SAFE_ALLOCATE(u_ma(1:mesh%np,hm%F%floquet_dim*hm%F%spindim))
      SAFE_ALLOCATE(tmp(gs_st%d%dim))


      do ik=dressed_st%d%kpt%start,dressed_st%d%kpt%end

        do ia=dressed_st%st_start,dressed_st%st_end
          call states_get_state(dressed_st, mesh, ia, ik, u_ma)
          dressed_st%occ(ia,ik) = M_ZERO 
          
          do ist=gs_st%st_start,gs_st%st_end
            call states_get_state(gs_st, mesh, ist, ik, psi)
            
            tmp(:) = M_ZERO
            do idx=hm%F%flat_idx%start, hm%F%flat_idx%end
               it = hm%F%idx_map(idx,1)
               im = hm%F%idx_map(idx,2)
               imm = im - Fdim(1) + 1
               do idim=1,gs_st%d%dim
                 tmp(idim)  =  tmp(idim) + exp(-M_zI*im*omega*it*dt)/nT * zmf_dotp(mesh, psi(1:mesh%np,idim), &
                                                                          u_ma(1:mesh%np,(imm-1)*gs_st%d%dim+idim))
               end do
               
            end do 
            if(hm%F%is_parallel) call comm_allreduce(hm%F%mpi_grp%comm, tmp(:))   
            
            dressed_st%occ(ia,ik) = dressed_st%occ(ia,ik) + abs(sum(tmp(:)))**2 * gs_st%occ(ist,ik)
!             print *, ist, "gs_st%occ(ist,:) =", gs_st%occ(ist,:)
!             print *, tmp(:)
            
          enddo
        enddo
        
        ! occupations checksum 
        sum_dr = sum(dressed_st%occ(:,ik))
        sum_gs = sum(gs_st%occ(:,ik))
        if( abs(sum_dr -sum_gs) > 1E-6) then
          call messages_write('Occupations checksum failed for kpoint = ')
          call messages_write(ik, fmt = '(i6)')
          call messages_write(':   gs_occ =  ')
          call messages_write(sum_gs, fmt ='(f12.6)')
          call messages_write('   floquet_occ =  ')
          call messages_write(sum_dr, fmt ='(f12.6)')
          call messages_warning()
        end if
      enddo
      
!       print *, "occsum = ", sum(dressed_st%occ(:,1))
!       print *, "gs_occsum = ", sum(gs_st%occ(:,1))


      
      SAFE_DEALLOCATE_A(tmp)
      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(u_ma)
    
      POP_SUB(floquet_calc_occupations)
    
    end subroutine floquet_calc_occupations
    


    
    
    

end module floquet_oct_m
