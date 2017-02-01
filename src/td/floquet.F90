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
  use comm_oct_m
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
  use parser_oct_m
  use partial_charges_oct_m
  use pert_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scf_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_dim_oct_m
  use states_restart_oct_m
  use system_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::                       &
       floquet_init,              &
       floquet_hamiltonians_init, &
       floquet_hamiltonian_update,&
       floquet_hamiltonian_solve  

  integer, public, parameter ::    &
       FLOQUET_NON_INTERACTING = 1, &
       FLOQUET_FROZEN_PHONON   = 2, &
       FLOQUET_INTERACTING     = 3

contains

  subroutine floquet_init(this,geo,dim)
    type(floquet_t),    intent(out)  :: this
    type(geometry_t), intent(in)     :: geo
    integer :: dim ! the standard dimension of the groundstate
    type(block_t)     :: blk
    integer :: ia, idir

    FLOAT :: time_step

    PUSH_SUB(floquet_init)

    !%Variable TDFloquetMode
    !%Type flag
    !%Default non_interacting
    !%Section Time-Dependent::TD Output
    !%Description
    !% Types of Floquet analysis performed when TDOutput=td_floquet
    !%Option non_interacting 1
    !% 
    !%Option frozen_phonon 2
    !% 
    !%Option interacting 3
    !%
    !%End

    call parse_variable('TDFloquetMode', 1,this%mode)

    if(this%mode==FLOQUET_FROZEN_PHONON) then
      
      !%Variable TDFloquetFrozenDistortion
      !%Type block
      !%Section Time-Dependent::TD Output
      !%Description
      !%
      !%End
      if(.not.parse_block('TDFloquetFrozenDistortion', blk) == 0) then
        write(message(1),'(a)') 'Internal error while reading TDFloquetFrozenDistortion.'
        call messages_fatal(1)
      end if
      if(parse_block_n(blk) /= geo%natoms) then
         write(message(1),'(a)') 'Please provide in then TDFloquetFrozenDistortion block an'
         write(message(2),'(a,i3,a)') 'initial distortion for all ',geo%natoms, ' atoms.'
         call messages_fatal(2)
      end if

      SAFE_ALLOCATE(this%frozen_distortion(1:geo%natoms,1:geo%space%dim))
      do ia=1,geo%natoms
        do idir = 1, geo%space%dim
          call parse_block_float(blk, ia-1, idir-1,this%frozen_distortion(ia,idir))
        end do
      end do

      write(message(1),'(a)') 'Initial distortions given in input file:'
      call messages_info(1)
      ! TODO: print here the values also
    end if

    !%Variable TDFloquetFrequency
    !%Type float
    !%Default 0
    !%Section Time-Dependent::TD Output
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

    !%Variable TDFloquetSample
    !%Type integer
    !%Default 20
    !%Section Time-Dependent::TD Output
    !%Description
    !% Number of points on which one Floquet cycle is sampled in the
    !%time-integral of the Floquet analysis.
    !%
    !%End
    call parse_variable('TDFloquetSample',20 ,this%nt)
    call messages_print_var_value(stdout,'Number of Floquet time-sampling points', this%nT)
    this%dt = this%Tcycle/real(this%nT)

    !%Variable TDFloquetMaximumSolverIterations
    !%Type integer
    !%Default 35
    !%Section Time-Dependent::TD Output
    !%Description
    !% Maximumn Number of calls to eigensolver for solving the Floquet Hamiltonian
    !%
    !%End
    call parse_variable('TDFloquetMaximumSolverIterations ', 35 ,this%max_solve_iter)
    call messages_print_var_value(stdout,'Maximum eigensolver iterations', this%max_solve_iter)


    !%Variable TDFloquetDimension
    !%Type integer
    !%Default -1
    !%Section Time-Dependent::TD Output
    !%Description
    !% Order of Floquet Hamiltonian. If negative number is given, downfolding
    !%is performed.
    !%End
    call parse_variable('TDFloquetDimension',-1,this%order)
    if(this%order.ge.0) then
      call messages_print_var_value(stdout,'Order of multiphoton Floquet-Hamiltonian', this%order)
      !Dimension of multiphoton Floquet-Hamiltonian
      this%floquet_dim = 2*this%order+1
    else
      message(1) = 'Floquet-Hamiltonian downfolding not implemented for interacting propagation.'
      call messages_fatal(1)
      !this%downfolding = .true.
      !this%Forder = 1
      !this%Fdim = 3
    endif

    this%count = 1
    this%spindim = dim

    ! re-read time step from input
    call parse_variable('TDTimeStep', M_ZERO, time_step, unit = units_inp%time)
    if(time_step == M_ZERO) then
       message(1) = 'Did not find time-step in Floquet init, plase give a value for TDTimeStep'
      call messages_fatal(1)
    end if
    this%interval = int(this%dt/time_step)
    this%ncycle = this%interval*this%nT

    call messages_print_var_value(stdout,'Steps in Floquet time-sampling interval',  this%interval)
    call messages_print_var_value(stdout,'Steps in Floquet time-sampling cycle',  this%ncycle)

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
integer :: nik, ist
FLOAT, allocatable :: frozen_bands(:,:)


    PUSH_SUB(floquet_hamiltonian_init)

    mesh = gr%der%mesh
    nst = st%nst

    !for now no domain distributionallowed
    ASSERT(mesh%np == mesh%np_global)

    ! the Hamiltonain gets assigned an array of td-Hamiltonians
    ! this is a bit recursive, so maybe there should be a Flqoeut moduel or something
    nullify(this%td_hm)
    SAFE_ALLOCATE(this%td_hm(1:this%F%nT))

    if(this%F%mode == FLOQUET_FROZEN_PHONON) then
       SAFE_ALLOCATE(frozen_bands(st%nst,gr%sb%kpoints%reduced%npoints))
    end if

    ! initialize the instances of the Hamiltonians
    do it=1,this%F%nT
       this%td_hm(it)%F%Tcycle = this%F%Tcycle
       this%td_hm(it)%F%omega = this%F%omega
       this%td_hm(it)%F%dt = this%F%dt

       SAFE_ALLOCATE(this%td_hm(it)%geo)
       if(.not.this%F%mode==FLOQUET_INTERACTING) then
          call geometry_copy(this%td_hm(it)%geo, this%geo)

          ! set flag to prevent species types to be touched, because
          ! hm%geo is only a pointer to the global geo instance
          this%td_hm(it)%geo%skip_species_pot_init = .true.
       end if

       select case(this%F%mode)

       case(FLOQUET_FROZEN_PHONON) 
         time = it*this%F%dt
         space_dim = this%geo%space%dim
         do ia=1,this%geo%natoms
           this%td_hm(it)%geo%atom(ia)%x(1:space_dim) = this%td_hm(it)%geo%atom(ia)%x(1:space_dim) + &
                                                this%F%frozen_distortion(ia,1:space_dim)*cos(time*this%F%omega)
         end do

         call hamiltonian_init(this%td_hm(it), gr, this%td_hm(it)%geo, st, &
                                     sys%ks%theory_level, sys%ks%xc_family,sys%ks%xc_flags)
         call hamiltonian_epot_generate(this%td_hm(it), gr, this%td_hm(it)%geo, st, time=M_ZERO)

         call scf_init(scf,gr,this%td_hm(it)%geo,st,this%td_hm(it))
         call scf_run(scf,sys%mc,gr,this%td_hm(it)%geo,st,sys%ks,this%td_hm(it),sys%outp, gs_run=.false.)
         call scf_end(scf)

         frozen_bands(st%nst,gr%sb%kpoints%reduced%npoints) = &
              frozen_bands(st%nst,gr%sb%kpoints%reduced%npoints) + M_ONE/this%F%nT*st%eigenval(st%nst,gr%sb%kpoints%reduced%npoints)

         write(filename,'(I5)') it
         filename = 'BO_bands_'//trim(adjustl(filename))
         open(unit=98765,file=filename)
         nik=gr%sb%kpoints%nik_skip
         do ik=gr%sb%kpoints%reduced%npoints-nik+1,gr%sb%kpoints%reduced%npoints
            do ist=1,st%nst
               write(98765,'(e12.6, 1x)',advance='no') st%eigenval(ist, ik)
            end do
            write(98765,'(1x)')
         end do
         close(98765)


       case(FLOQUET_NON_INTERACTING)
         call hamiltonian_init(this%td_hm(it), gr, this%td_hm(it)%geo, st, &
                                    sys%ks%theory_level, sys%ks%xc_family,sys%ks%xc_flags)
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
        open(unit=98765,file='frozen_bands')
        nik=gr%sb%kpoints%nik_skip
        do ik=gr%sb%kpoints%reduced%npoints-nik+1,gr%sb%kpoints%reduced%npoints
           do ist=1,st%nst
              write(98765,'(e12.6, 1x)',advance='no') frozen_bands(ist, ik)
           end do
           write(98765,'(1x)')
        end do
        close(98765)
        
     end if


     POP_SUB(floquet_hamiltonian_init)
        
   end subroutine floquet_hamiltonians_init

   !--------------------------------------------
   ! this is only called if F%mode=interacting
   subroutine floquet_hamiltonian_update(hm,st,gr,sys,iter)
     type(hamiltonian_t), intent(inout) :: hm
     type(states_t)      , intent(inout):: st
     type(grid_t),      intent(inout)   :: gr
     type(system_t),      intent(in)    :: sys
     integer :: iter
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

     call hamiltonian_init(hm%td_hm(it), gr, hm%td_hm(it)%geo, st, &
           sys%ks%theory_level, sys%ks%xc_family,sys%ks%xc_flags)
     
     call hamiltonian_update(hm%td_hm(it), gr%der%mesh,time=time)

     ! fill time-dependent hamiltonian structure with scf-fields at this time
     do ispin=1,hm%d%nspin
       forall (ip = 1:gr%mesh%np) hm%td_hm(it)%vhxc(ip,ispin) = hm%vhxc(ip, ispin)
     end do
     forall (ip = 1:gr%mesh%np) hm%td_hm(it)%ep%vpsl(ip) = hm%ep%vpsl(ip)

     call hamiltonian_epot_generate(hm%td_hm(it), gr, hm%td_hm(it)%geo, st, time=time)

     PUSH_SUB(floquet_hamiltonian_update)

    end subroutine floquet_hamiltonian_update

    !--------------------------------------------
    subroutine floquet_hamiltonian_solve(out_floquet,hm,gr,sys,st)
      type(c_ptr),       intent(inout)   :: out_floquet
      type(hamiltonian_t), intent(inout) :: hm
      type(grid_t),      intent(inout)   :: gr
      type(system_t), intent(inout)      :: sys
      type(states_t), intent(in)         :: st

      logical :: converged
      integer :: iter , maxiter, ik, in, im, ist, idim, ierr, nik, dim, nst
      CMPLX, allocatable :: temp_state1(:,:), temp_state2(:,:)
      type(eigensolver_t) :: eigens
      type(states_t) :: dressed_st
      type(restart_t) :: restart

! temporary variables for outputting 
character(len=80) :: filename
integer :: file

      ! initialize a state object with the Floquet dimension
      call states_init(dressed_st, gr, hm%geo,floquet_dim=hm%F%floquet_dim)
      call kpoints_distribute(dressed_st%d,sys%mc)
      call states_distribute_nodes(dressed_st,sys%mc)
      call states_allocate_wfns(dressed_st,gr%der%mesh)

      ! solver iteration
      iter = 0

      call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_LOAD, &
                        dressed_st%dom_st_kpt_mpi_grp, ierr, gr%der%mesh)
      if(ierr == 0) then
         call states_look(restart, nik, dim, nst, ierr)
         if(dim==dressed_st%d%dim .and. nik==gr%sb%kpoints%reduced%npoints &
              .and. nst==dressed_st%nst) then
            call states_load(restart, dressed_st, gr, ierr, iter)
         else
            ! this doesnt need to be fatal
            write(message(1),'(a)') 'Floquet restart structure not commensurate.'
            call messages_fatal(1)
         end if
      else
         ! initialize floquet states from scratch
         SAFE_ALLOCATE(temp_state1(1:gr%der%mesh%np,st%d%dim))
         SAFE_ALLOCATE(temp_state2(1:gr%der%mesh%np,hm%F%floquet_dim))
         
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
      end if

      call restart_end(restart)

      hm%F%floquet_apply = .true.
      ! set dimension of Floquet Hamiltonian                                                    
      hm%d%dim = dressed_st%d%dim
      
      call eigensolver_init(eigens, gr, dressed_st)
      ! no subspace diag implemented yet
      eigens%sdiag%method = OPTION__SUBSPACEDIAGONALIZATION__NONE

      converged=.false.
      maxiter = hm%F%max_solve_iter
      do while(.not.converged.and.iter <= maxiter)
         call eigensolver_run(eigens, gr, dressed_st, hm, 1,converged)

! this will be replaced with a a call to plot_bandstructure
if(mpi_world%rank==0) then
   file = 987654
   write(filename,'(I5)') iter !hm%F_count
   filename = 'floquet_multibands_'//trim(adjustl(filename))
   open(unit=file,file=filename)
   ! we are only interested in k-points with zero weight
   nik=gr%sb%kpoints%nik_skip
   do ik=gr%sb%kpoints%reduced%npoints-nik+1,gr%sb%kpoints%reduced%npoints
      do ist=1,dressed_st%nst
         write(file,'(e12.6, 1x)',advance='no') dressed_st%eigenval(ist, ik)
      end do
      write(file,'(1x)')
   end do
   close(file)
endif

         iter = iter +1
      end do

      call eigensolver_end(eigens)
      !
      
      !switch off floquet hamiltonian                                                           
      hm%F%floquet_apply = .false.                                                               
      ! reset dimension
      hm%d%dim = hm%F%spindim
      
      hm%F%count=hm%F%count + 1
      
      ! here we might want to do other things with the states...

      ! write states
      call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_DUMP, &
                        dressed_st%dom_st_kpt_mpi_grp, ierr, gr%der%mesh)
      iter = iter -1
      call states_dump(restart, dressed_st, gr, ierr, iter)

      call states_end(dressed_st)



    end subroutine floquet_hamiltonian_solve

    !---------------------------------------
    !subroutine floquet_end()

end module floquet_oct_m
