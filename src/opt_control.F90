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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module opt_control_m
  use excited_states_m
  use datasets_m
  use varinfo_m
  use global_m
  use filter_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lasers_m
  use lib_oct_parser_m
  use lib_oct_m
  use messages_m
  use mesh_m
  use mesh_function_m
  use output_m
  use geometry_m
  use states_m
  use states_output_m
  use string_m
  use system_m
  use td_rti_m
  use td_write_m
  use timedep_m
  use units_m
  use v_ks_m
  use external_pot_m
  use restart_m
  use tdf_m

  implicit none

  private
  public :: opt_control_run

  integer, parameter, private  ::  &
    oct_is_groundstate   = 1,      &
    oct_is_excited       = 2,      &
    oct_is_superposition = 3,      &
    oct_is_userdefined   = 4         
    
  integer, parameter, private  ::  &
    oct_tg_groundstate   = 1,      &
    oct_tg_excited       = 2,      &
    oct_tg_superposition = 3,      &
    oct_tg_userdefined   = 4,      &
    oct_tg_local         = 5,      &        
    oct_tg_td_local      = 6   

  integer, parameter, private  ::  &
    oct_algorithm_zbr98 = 1,       &
    oct_algorithm_zr98  = 2,       &
    oct_algorithm_wg05  = 3       
 
  integer, parameter, private  ::  &
    oct_targetmode_static = 0,     &
    oct_targetmode_td     = 1
  
  integer, parameter, private  ::  &
    oct_tgtype_state = 1,          &
    oct_tgtype_local = 2
  
  integer, parameter, private  ::  &
    oct_ftype_gauss = 1,           &
    oct_ftype_step  = 2

  type oct_t
    FLOAT   :: targetfluence
    logical :: mode_fixed_fluence
    integer :: targetmode
    integer :: filtermode
    integer :: totype
    integer :: algorithm_type
    logical :: td_tg_state 
    logical :: oct_double_check
    integer :: istype
    logical :: dump_intermediate
  end type oct_t

  type td_target_t
    integer :: type       ! local or state 
    FLOAT   :: par 
    FLOAT   :: width      ! width parameter
    FLOAT   :: weight     ! relative weight of filter
    character(len=1024) :: expression(MAX_DIM) ! expression of user def func
    integer             :: ftype               ! which function: gaussian etc
    CMPLX, pointer :: tdshape(:,:) ! evaluate user defined functions
  end type td_target_t

  type td_target_set_t
    integer :: no_tdtargets
    type(td_target_t), pointer :: tdtg(:)
    FLOAT, pointer :: td_fitness(:)
    CMPLX, pointer :: tdtarget(:)
  end type td_target_set_t

  type oct_iterator_t
    FLOAT              :: eps
    integer            :: ctr_iter_max
    integer            :: ctr_iter
    FLOAT              :: functional
    FLOAT              :: old_functional
    FLOAT              :: overlap
    FLOAT, pointer     :: convergence(:,:)
    FLOAT              :: bestJ, bestJ1, bestJ_fluence, bestJ1_fluence
    FLOAT              :: bestJ_J1, bestJ1_J
    integer            :: bestJ_ctr_iter, bestJ1_ctr_iter
  end type oct_iterator_t

  type oct_control_parameters_t
    integer :: no_parameters
    FLOAT   :: dt
    integer :: ntiter
    type(tdf_t), pointer :: f(:)
    CMPLX, pointer  :: laser_pol(:, :)
  end type oct_control_parameters_t

  type oct_penalty_t 
    logical          :: mode_tdpenalty
    FLOAT, pointer   :: a_penalty(:)
    type(tdf_t), pointer :: td_penalty(:)
  end type oct_penalty_t

contains

  ! ---------------------------------------------------------
  subroutine opt_control_run(sys, h)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h

    type(oct_t)                :: oct
    type(oct_penalty_t)        :: penalty
    type(oct_iterator_t)       :: iterator
    type(td_t)                 :: td
    type(grid_t),      pointer :: gr   ! some shortcuts
    type(states_t),    pointer :: st
    type(states_t)             :: chi, psi, target_st, initial_st, psi2
    type(filter_t),    pointer :: f(:)
    type(td_target_set_t)      :: tdt
    type(oct_control_parameters_t) :: par, par_tmp
    integer :: i, ierr
    character(len=80)  :: filename

    call push_sub('opt_control.opt_control_run')

    ! Checks that the run is actually possible with the current settings
    call check_runmode_constrains(sys, h)

    call io_mkdir('opt-control')

    ! some shortcuts
    gr  => sys%gr
    st  => sys%st

    call td_init(sys, h, td)
    call states_allocate_wfns(st, gr%m, M_CMPLX)

    ! Initialize a bunch of states: initial, target, auxiliary.
    psi        = st
    psi2       = st
    chi        = st
    initial_st = st
    target_st  = st

    call oct_read_inp(oct)

    ! Initial guess for the laser: read from the input file.
    call laser_init(h%ep%no_lasers, h%ep%lasers, gr%m)

    do i = 1, h%ep%no_lasers
      call laser_to_numerical(h%ep%lasers(i), td%dt, td%max_iter)
    end do

    call parameters_init(par, h%ep%no_lasers, td)
    call parameters_init(par_tmp, h%ep%no_lasers, td)
    call parameters_set(par, h%ep)
    call parameters_set(par_tmp, h%ep)
    call parameters_write('opt-control/initial_laser', par)

    write(message(1),'(a,f14.8)') 'Input: Fluence of Initial laser ', laser_fluence(par)
    call write_info(1)

    call oct_iterator_init(iterator, oct)

    call penalty_init(penalty, oct, par, td%max_iter, iterator%ctr_iter_max, td%dt)

    if(oct%filtermode.gt.0) then
      nullify(f)
      call def_filter(td%max_iter, td%dt, oct%filtermode, f)
      call filter_write(f, NDIM, td%dt, td%max_iter)
    end if

    call def_istate(oct, gr, sys%geo, initial_st)
    call def_toperator(oct, gr, sys%geo, target_st)
    call tdtargetset_init(oct, gr, td, tdt)

    call check_faulty_runmodes(oct, penalty)

    call states_output(initial_st, gr, 'opt-control/initial1', sys%outp)
    call states_output(target_st, gr, 'opt-control/target1', sys%outp)

    ! mode switcher
    select case(oct%algorithm_type)
      case(oct_algorithm_zbr98) ;  call scheme_zbr98
      case(oct_algorithm_zr98)  ;  call scheme_zr98
      case(oct_algorithm_wg05)  ;  call scheme_wg05
    case default
      write(message(1),'(a)') "The OCT algorithm has not been implemented yet."
      write(message(2),'(a)') "Choose: ZR98, ZBR98, WG05."
      call write_fatal(2)
    end select

    ! Some informative output.
    call output(oct, iterator)
    call states_output(target_st, gr, 'opt-control/target', sys%outp)
    call states_output(psi, gr, 'opt-control/final', sys%outp)

    ! do final test run: propagate initial state with optimal field
    call parameters_to_h(par, h%ep)
    call oct_finalcheck(oct, initial_st, target_st, sys, h, td, tdt)

    ! clean up
    call parameters_end(par)
    call parameters_end(par_tmp)
    call oct_iterator_end(iterator)
    call tdtargetset_end(tdt)
    call filter_end(f)
    call td_end(td)
    call pop_sub()

  contains

    
    ! ---------------------------------------------------------
    subroutine scheme_zr98
      integer :: ierr
      call push_sub('opt_control.scheme_ZR98')

      message(1) = "Info: Starting Optimal Control Run  using Scheme: ZR98"
      call write_info(1)
      
      ! first propagate chi to ti
      message(1) = "Info: Initial forward propagation"
      call write_info(1)

      psi = initial_st
      call parameters_to_h(par, h%ep)
      call propagate_forward(oct, sys, h, td, tdt, psi) 
      
      if(oct%dump_intermediate) call states_output(psi, gr, 'opt-control/prop1', sys%outp)
        
      call states_output(initial_st, gr, 'opt-control/zr98_istate1', sys%outp)

      ctr_loop: do
        
         ! defines chi
        call target_calc(oct, gr, target_st, psi, chi)
        
        if(iteration_manager(oct, penalty, gr, tdt%td_fitness, par, td, psi, target_st, iterator)) exit ctr_loop
        
        call bwd_step(oct_algorithm_zr98)
        
        if(oct%dump_intermediate) then
          write(filename,'(a,i3.3)') 'opt-control/b_laser.', iterator%ctr_iter
          call parameters_write(filename, par_tmp)
        end if
        
        ! forward propagation
        psi = initial_st
        psi2 = initial_st
        call states_output(psi, gr, 'opt-control/last_bwd', sys%outp)
        call fwd_step(oct_algorithm_zr98)
        write(filename,'(a,i3.3)') 'opt-control/PsiT.', iterator%ctr_iter
        call states_output(psi, gr, filename, sys%outp)
        
        if(oct%dump_intermediate) then
          write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
          call parameters_write(filename, par)
        end if
        
      end do ctr_loop
      
      call pop_sub()
    end subroutine scheme_zr98
    

    !---------------------------------------
    subroutine scheme_wg05
      integer :: ierr
      FLOAT :: fluence
      call push_sub('opt_control.scheme_WG05')
      
      message(1) = "Info: Starting OCT iteration using scheme: WG05"
      call write_info(1)
      
      ctr_loop: do
         
        ! first propagate chi to ti
        message(1) = "Info: Initial forward propagation"
        call write_info(1)
        
        psi = initial_st
        call parameters_to_h(par, h%ep)
        call propagate_forward(oct, sys, h, td, tdt, psi) 
        
        if(oct%dump_intermediate) then
          write(filename,'(a,i3.3)') 'opt-control/PsiT.', iterator%ctr_iter
          call states_output(psi, gr, 'opt-control/last_bwd', sys%outp)
        end if


        ! define target state
        call target_calc(oct, gr, target_st, psi, chi) ! defines chi
        
        if(iteration_manager(oct, penalty, gr, tdt%td_fitness, par, td, psi, target_st, iterator)) exit ctr_loop
        
        call bwd_step(oct_algorithm_wg05)
        !!!WARNING: this probably shoudl not be here?
        fluence = laser_fluence(par_tmp)

        ! WARNING: This assumes that there is only one control parameter!
        ! WARNING: Untested.
        if (oct%filtermode.gt.0) call apply_filter(td%max_iter, f, par_tmp%f(1))

        ! recalc field
        if (oct%mode_fixed_fluence) then
          !!!WARNING: this probably shoudl not be here?
          fluence = laser_fluence(par_tmp)

          penalty%a_penalty(iterator%ctr_iter + 1) = &
            sqrt(fluence * penalty%a_penalty(iterator%ctr_iter)**2 / oct%targetfluence )

          if(oct%dump_intermediate) then
            write (6,*) 'actual penalty', penalty%a_penalty(iterator%ctr_iter)
            write (6,*) 'next penalty', penalty%a_penalty(iterator%ctr_iter + 1)
          end if

          !!!!!WARNING: This is a very strange statement?
          call tdf_set_numerical(penalty%td_penalty(1), &
                                 spread(penalty%a_penalty(iterator%ctr_iter + 1), 1, td%max_iter+1) )

          call tdf_scalar_multiply( ( penalty%a_penalty(iterator%ctr_iter) / penalty%a_penalty(iterator%ctr_iter + 1) ), &
                                    par_tmp%f(1) )

          fluence = laser_fluence(par_tmp)
        end if

        par%f(1) = par_tmp%f(1)
        
        ! dump here: since the fwd_step is missing
        write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
        call parameters_write(filename, par_tmp)
        
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_wg05

 
    ! ---------------------------------------------------------
    subroutine scheme_zbr98
      integer :: ierr
      call push_sub('opt_control.scheme_zbr98')
      
      message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
      call write_info(1)
      
      ! first propagate chi to T
      message(1) = "Info: Initial backward propagation"
      call write_info(1)
      
      call target_calc(oct, gr, target_st, psi, chi)

      call parameters_to_h(par, h%ep)
      call propagate_backward(sys, h, td, chi)

      if(oct%dump_intermediate) call states_output(chi, gr, 'opt-control/initial_propZBR98', sys%outp)

      ctr_loop: do

        ! check for stop file and delete file
        message(1) = "Info: Setup forward"
        call write_info(1)
        
        ! forward propagation
        psi = initial_st
        
        call fwd_step(oct_algorithm_zbr98)

        write(filename,'(a,i3.3)') 'opt-control/PsiT.', iterator%ctr_iter
        call states_output(psi, gr, filename, sys%outp)
        
        if(iteration_manager(oct, penalty, gr, tdt%td_fitness, par, td, psi, target_st, iterator)) exit ctr_loop
        
         ! and now backward
        call target_calc(oct, gr, target_st, psi, chi)
        call bwd_step(oct_algorithm_zbr98)
        call states_output(psi, gr, 'opt-control/last_bwd', sys%outp)
        
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_zbr98


    ! ----------------------------------------------------------
    subroutine fwd_step(method)
      integer, intent(in) :: method
      integer :: i
      FLOAT, allocatable :: dens_tmp(:,:)

      call push_sub('opt_control.fwd_step')

      ALLOCATE(dens_tmp(gr%m%np_part, st%d%nspin), gr%m%np_part*st%d%nspin) 
      
      ! setup forward propagation
      call states_densities_init(psi, gr, sys%geo)
      call states_calc_dens(psi, NP_PART, dens_tmp)
      psi%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
      call td_rti_run_zero_iter(h, td%tr)

      message(1) = "Info: Propagating forward"
      call write_info(1)

      do i = 1, td%max_iter
        call prop_iter_fwd(i, method)
        if(oct%targetmode==oct_targetmode_td) &
          call calc_tdfitness(tdt, gr, psi, tdt%td_fitness(i))     
      end do
      
      ! dump new laser field
      !if(oct%dump_intermediate) then
      !  write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
      !  call parameters_write(filename, par)
      !endif

      deallocate(dens_tmp)      
      call pop_sub()
    end subroutine fwd_step
    
    
    ! --------------------------------------------------------
    subroutine bwd_step(method) 
      integer, intent(in) :: method
      integer :: i
      FLOAT, allocatable :: dens_tmp(:,:)

      call push_sub('opt_control.bwd_step')

      ALLOCATE(dens_tmp(gr%m%np_part, st%d%nspin), gr%m%np_part*st%d%nspin) 
      ! setup backward propagation
      call states_calc_dens(chi, NP_PART, dens_tmp)
      chi%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
      call td_rti_run_zero_iter(h, td%tr)

      message(1) = "Info: Propagating backward"
      call write_info(1)
      !call loct_progress_bar(td%max_iter-1, td%max_iter-1)
      td%dt = -td%dt
      do i = td%max_iter-1, 0, -1
        call prop_iter_bwd(i,method)
        !call loct_progress_bar(td%max_iter-1-i, td%max_iter-1)
      end do
      td%dt = -td%dt

      ! This is not the place to write anything.
      !if(oct%dump_intermediate) then
      !  write(filename,'(a,i3.3)') 'opt-control/b_laser.', iterator%ctr_iter
      !  call parameters_write(filename, par_tmp)
      !end if

      deallocate(dens_tmp)
      call pop_sub()
    end subroutine bwd_step


    ! ---------------------------------------------------------
    subroutine prop_iter_fwd(iter,method)
      integer, intent(in) :: method
      integer,           intent(in) :: iter
      FLOAT, allocatable :: dens_tmp(:,:)
      
      call push_sub('opt_control.prop_iter_fwd')
 
      ALLOCATE(dens_tmp(gr%m%np_part, st%d%nspin), gr%m%np_part*st%d%nspin) 
      
      ! chi(0) --> chi(T) [laser_tmp]
      !         |------------laser     
      ! psi(0) --> psi(T) 
      
      ! new electric field
      if(oct%targetmode==oct_targetmode_td) then
        ! psi(0) --> psi(T) with old field (for chi)
        call calc_inh(psi2, gr, tdt, iter, td%max_iter, td%dt, chi)
        ! psi2
        call parameters_to_h(par, h%ep)
        call states_calc_dens(psi2, NP_PART, dens_tmp)
        psi2%rho = dens_tmp
        call v_ks_calc(gr, sys%ks, h, psi2, calc_eigenval=.true.)
        call td_rti_dt(sys%ks, h, gr, psi2, td%tr, abs(iter*td%dt), abs(td%dt))
      end if

      call update_field(oct, penalty, iter-1, par, gr, td, psi, chi)
      
      ! chi
      call parameters_to_h(par_tmp, h%ep)
      call states_calc_dens(chi, NP_PART, dens_tmp)
      chi%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.) 
      call td_rti_dt(sys%ks, h, gr, chi, td%tr, abs(iter*td%dt), abs(td%dt))
      
      ! psi
      call parameters_to_h(par, h%ep)
      call states_calc_dens(psi, NP_PART, dens_tmp)
      psi%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, abs(iter*td%dt), abs(td%dt))

      deallocate(dens_tmp)     
      call pop_sub()
    end subroutine prop_iter_fwd
   
   
    ! ---------------------------------------------------------
    ! do backward propagation step with update of field
    ! works for time-independent and td targets
    subroutine prop_iter_bwd(iter, method)
      integer, intent(in) :: iter
      integer, intent(in) :: method
      FLOAT, allocatable :: dens_tmp(:,:)
      
      call push_sub('opt_control.prop_iter_bwd')
 
      ALLOCATE(dens_tmp(gr%m%np_part, st%d%nspin), gr%m%np_part*st%d%nspin) 
 
      ! chi(T) --> chi(0)
      !         |------------laser_tmp
      ! psi(T) --> psi(0) [laser]
      
      ! calc inh first, then update field
      if(oct%targetmode==oct_targetmode_td) then
        call calc_inh(psi, gr, tdt, iter, td%max_iter, td%dt, chi)
      end if
 
      ! new electric field
      call update_field(oct, penalty, iter+1, par_tmp, gr, td, psi, chi)
           
      ! chi
      call parameters_to_h(par_tmp, h%ep)
      call states_calc_dens(chi, NP_PART, dens_tmp)
      chi%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, chi, td%tr, abs(iter*td%dt),     td%dt )
      
      ! psi
      call  parameters_to_h(par, h%ep)
      call states_calc_dens(psi, NP_PART, dens_tmp)
      psi%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, abs(iter*td%dt),     td%dt )
 
      deallocate(dens_tmp)     
      call pop_sub()
    end subroutine prop_iter_bwd
 
  end subroutine opt_control_run


  ! ---------------------------------------------------------
  ! This subroutine just stops the run if some of the settings
  ! are not compatible with the run mode, either because it is
  ! meaningless, or because the run mode is still not fully
  ! implemented.
  !
  ! TODO: Right now just a couple of checks are made, but there are many other constrains.
  subroutine check_runmode_constrains(sys, h)
    type(system_t), target, intent(in) :: sys
    type(hamiltonian_t),    intent(in) :: h

    integer :: no_electrons

    ! Only dipole approximation in length gauge.
    if(h%gauge.ne.LENGTH) then
      write(message(1),'(a)') "So far only length gauge is supported in optimal control runs."
      call write_fatal(1)
    end if

    ! Only single-electron calculations.
    no_electrons = nint(sum(sys%st%occ(1,:)))
    if(no_electrons .ne. 1) then
      write(message(1),'(a,i4)') 'Number of electrons (currently ',no_electrons,') should be just one.'
      write(message(2),'(a)')    'Optimal control theory for many electron systems not yet developed!'
      call write_fatal(2)
    end if

  end subroutine check_runmode_constrains

#include "opt_control_read.F90"
#include "opt_control_aux.F90"
#include "opt_control_defstates.F90"
#include "opt_control_propagation.F90"
#include "opt_control_finalcheck.F90"
#include "opt_control_iter.F90"
#include "opt_control_output.F90"
#include "opt_control_parameters.F90"
#include "opt_control_penalty.F90"
#include "opt_control_tdtarget.F90"

end module opt_control_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
