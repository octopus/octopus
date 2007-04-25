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
  use states_m
  use states_output_m
  use string_m
  use system_m
  use td_rti_m
  use td_write_m
  use timedep_m
  use units_m
  use v_ks_m

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
    FLOAT :: eps
    integer :: ctr_iter_max
    FLOAT :: penalty
    FLOAT, pointer :: a_penalty(:)
    FLOAT :: targetfluence
    logical :: mode_fixed_fluence
    FLOAT, pointer :: tdpenalty(:, :)
    integer :: targetmode
    integer :: filtermode

    CMPLX              :: laser_pol(MAX_DIM)
    integer :: dof

    integer :: algorithm_type

    logical :: oct_double_check
  end type oct_t


  type td_target_t
    integer :: type       ! local or state 
    FLOAT   :: par 
    FLOAT   :: width      ! width parameter
    FLOAT   :: weight     ! relative weight of filter
    
    character(len=1024) :: expression(MAX_DIM) ! expression of user def func
    integer             :: ftype               ! which function: gaussian etc
    
    FLOAT, pointer :: spatial(:)   ! x,y,z dependence
    CMPLX, pointer :: tdshape(:,:) ! evaluate user defined functions

  end type td_target_t

contains

  ! ---------------------------------------------------------
  subroutine opt_control_run(sys, h)
    implicit none

    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h

    type(oct_t)                :: oct
    type(td_t)                 :: td
    type(grid_t),      pointer :: gr   ! some shortcuts
    type(states_t),    pointer :: st
    type(states_t) :: psi_i
    type(states_t)             :: chi, psi, target_st, initial_st, psi2
    type(filter_t),    pointer :: f(:)
    type(td_target_t), pointer :: td_tg(:)
    type(td_write_t)           :: write_handler

    integer            :: No_electrons
    FLOAT, pointer     :: v_old_i(:,:,:), v_old_f(:,:,:)
    FLOAT, pointer     :: laser_tmp(:,:), laser(:,:)
    FLOAT, allocatable :: dens_tmp(:,:)
    integer            :: i, ctr_iter
    integer            :: iunit
    FLOAT              :: overlap, functional, old_functional
    FLOAT              :: fluence, u, t
    FLOAT, allocatable :: convergence(:,:), field(:), td_fitness(:)
    CMPLX, allocatable :: tdtarget(:)
    
    FLOAT              :: bestJ, bestJ1, bestJ_fluence, bestJ1_fluence
    FLOAT              :: bestJ_J1, bestJ1_J
    integer            :: bestJ_ctr_iter, bestJ1_ctr_iter

    character(len=80)  :: filename

    integer            :: istype, totype
    
    logical            :: dump_intermediate, mode_fixed_fluence
    logical            :: mode_tdpenalty, td_tg_state 


    call push_sub('opt_control.opt_control_run')

    ! Initially nullify the pointers
    nullify(f)
    nullify(td_tg)
    nullify(v_old_i)
    nullify(v_old_f)
    nullify(laser_tmp)
    nullify(laser)


    ! TODO: internal debug config, give some power to the user later
    dump_intermediate  = .TRUE.  ! dump laser fields during iteration
                                 ! this might create a lot of data

    mode_tdpenalty     = .FALSE. ! if TD penalty is used

    td_tg_state        = .FALSE. ! check if state target exist

    ! CHECK:: only dipole approximation in length gauge
    !         only single particle (yet)
    if(h%gauge.ne.1) then
      write(message(1),'(a)') "So far only length gauge is supported..."
      call write_fatal(1)
    end if

    ! initialize oct run, defines initial laser... initial state
    call init_()

    ! mode switcher
    select case(oct%algorithm_type)
      
    case(oct_algorithm_zbr98)    
      ! (0)  rabitz zbr98 type algorithm (states only)    
      call scheme_ZBR98(oct%algorithm_type) 

    case(oct_algorithm_zr98)     
      ! (1)  rabitz zr98  type algorithm 
      call scheme_ZR98(oct%algorithm_type)  
      ! TODO: o generalize elements
      !       o time-dependent targets
      
    case(oct_algorithm_wg05)  
      ! (2) Werschnik, Gross Type: allows fixed fluence and filtering
      call scheme_WG05(oct%algorithm_type)  
      ! (11) td targets rabitz type
      ! (12) td targets wg     type (with filter)
      
    ! TODO: Krotov scheme
    ! MT03 scheme: with additional parameters
   
    case DEFAULT
      write(message(1),'(a)') "Unknown choice of OCT algorithm."
      write(message(2),'(a)') "Choose: ZR98, ZBR98, WG05 ."
      call write_fatal(2)
      
    end select
    
    ! output: States, Lasers, Convergence
    call output()

    ! do final test run: propagate initial state with optimal field
    if(oct%oct_double_check) then
      message(1) = "Info: Optimization finished...checking the field"
      call write_info(1)
      
      call states_copy(psi, initial_st)
      
      call propagate_forward(psi)

      if(oct%targetmode==oct_targetmode_td) then
        overlap = SUM(td_fitness) * abs(td%dt)
        write(message(1), '(a,f14.8)') " => overlap:", overlap
        call write_info(1)
      else
        overlap = abs(zstates_mpdotp(gr%m, psi, target_st))
        write(message(1), '(6x,a,f14.8)') " => overlap:", overlap
        call write_info(1)
      end if
    end if

    ! clean up
    td%tr%v_old => v_old_i
    nullify(h%ep%lasers(1)%numerical)
    deallocate(laser)
    deallocate(laser_tmp)

    deallocate(convergence)
!! PROBLEMS WITH DEALLOCATION
    !write(6,*) 'DEBUG5'
    !nullify(v_old_f)
    !deallocate(v_old_f)
!! ----
    call states_end(psi_i)
    call end_()
    
    call pop_sub()

  contains

    
    ! ---------------------------------------------------------
    subroutine scheme_zr98(method)
      integer, intent(in) :: method
      
      integer :: ierr
      logical :: stoploop
     
      call push_sub('opt_control.scheme_ZR98')

      message(1) = "Info: Starting Optimal Control Run  using Scheme: ZR98"
      call write_info(1)
      
      stoploop = .FALSE.
      ! first propagate chi to ti
      message(1) = "Info: Initial forward propagation"
      call write_info(1)

      call states_copy(psi, initial_st)
      call propagate_forward(psi) 
      
      if(in_debug_mode) call zoutput_function(sys%outp%how, &
        'opt-control', 'prop1', gr%m, gr%sb, psi%zpsi(:,1,1,1), M_ONE, ierr)
        
      call zoutput_function(sys%outp%how,'opt-control','zr98_istate1',gr%m,gr%sb,initial_st%zpsi(:,1,1,1),M_ONE,ierr)
      old_functional = -CNST(1e10)
      ctr_iter = 0
      
      ctr_loop: do
        
        ! check for stop file and delete file
        
        ! define target state
        call target_calc(method, target_st, psi, chi) ! defines chi
        
        call iteration_manager(stoploop)
        if (stoploop)  exit ctr_loop
        
        call bwd_step(method)
        
        if(dump_intermediate) then
          write(filename,'(a,i3.3)') 'opt-control/b_laser.', ctr_iter
          call write_field(filename, laser_tmp, td%max_iter, td%dt)
        end if
        
        ! forward propagation
        call states_copy(psi, initial_st)
        call states_copy(psi2, initial_st) ! for td target
        call zoutput_function(sys%outp%how,'opt-control','last_bwd',gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr)
        call fwd_step(method)
        write(filename,'(a,i3.3)') 'PsiT.', ctr_iter
        call zoutput_function(sys%outp%how,'opt-control',filename,gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr)
        
        if(dump_intermediate) then
          write(filename,'(a,i3.3)') 'opt-control/laser.', ctr_iter
          call write_field(filename, laser, td%max_iter, td%dt)
        end if
        
      end do ctr_loop
      
      call pop_sub()
    end subroutine scheme_zr98
    

    !---------------------------------------
    subroutine scheme_wg05(method)
      integer, intent(in) :: method
      
      integer :: ierr
      logical :: stoploop

      call push_sub('opt_control.scheme_WG05')
      
      stoploop = .FALSE.
    
      message(1) = "Info: Starting OCT iteration using scheme: WG05"
      call write_info(1)
      
      old_functional = -CNST(1e10)
      ctr_iter = 0

      ctr_loop: do
         
        ! first propagate chi to ti
        message(1) = "Info: Initial forward propagation"
        call write_info(1)
        
        call states_copy(psi, initial_st)
        call propagate_forward(psi) 
        
        if(in_debug_mode) then
          write(filename,'(a,i3.3)') 'PsiT.', ctr_iter
          call zoutput_function(sys%outp%how,'opt-control',filename,gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr)
        end if

        ! define target state
        call target_calc(method,target_st,psi,chi) ! defines chi
        
        call iteration_manager(stoploop)
        if (stoploop)  exit ctr_loop
        
        call bwd_step(method)
        call calc_fluence(laser_tmp, fluence)      
        
        ! filter stuff / alpha stuff
        if (oct%filtermode.gt.0) & 
          call apply_filter(gr, td%max_iter, oct%filtermode, f, laser_tmp)

        ! recalc field
        if (oct%mode_fixed_fluence) then
          call calc_fluence(laser_tmp, fluence)      

          if(in_debug_mode) write (6,*) 'actual fluence', fluence

          oct%a_penalty(ctr_iter + 1) = &
            sqrt(fluence * oct%a_penalty(ctr_iter)**2 / oct%targetfluence )

          if(in_debug_mode) then
            write (6,*) 'actual penalty', oct%a_penalty(ctr_iter)
            write (6,*) 'next penalty', oct%a_penalty(ctr_iter + 1)
          end if

          oct%tdpenalty = oct%a_penalty(ctr_iter + 1)
          laser_tmp = laser_tmp * oct%a_penalty(ctr_iter) / oct%a_penalty(ctr_iter + 1)
          call calc_fluence(laser_tmp, fluence)

          if(in_debug_mode) write (6,*) 'renormalized', fluence, oct%targetfluence

        end if

        laser = laser_tmp
        
        ! dump here: since the fwd_step is missing
        write(filename,'(a,i3.3)') 'opt-control/laser.', ctr_iter
        call write_field(filename, laser, td%max_iter, td%dt)   
        
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_wg05

 
    ! ---------------------------------------------------------
    subroutine scheme_zbr98(method)
      integer, intent(in)   :: method

      integer :: ierr
      logical :: stoploop
      
      call push_sub('opt_control.scheme_zbr98')
      
      message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
      call write_info(1)
      
      stoploop = .FALSE.

      ! first propagate chi to T
      message(1) = "Info: Initial backward propagation"
      call write_info(1)
      
      call target_calc(method, target_st, psi, chi)
      call propagate_backward(chi) 

      if(in_debug_mode) call zoutput_function(sys%outp%how, &
        'opt-control', 'initial_propZBR98', gr%m, gr%sb, chi%zpsi(:,1,1,1), M_ONE, ierr)

      old_functional = -CNST(1e10)
      ctr_iter = 0

      ctr_loop: do

        ! check for stop file and delete file
        message(1) = "Info: Setup forward"
        call write_info(1)
        
        ! forward propagation
        call states_copy(psi, initial_st)  
        
        call fwd_step(method)

        if(in_debug_mode) then
          write(6,*) 'norm', zstates_dotp(gr%m, psi%d%dim,  psi%zpsi(:,:, 1, 1), psi%zpsi(:,:, 1, 1))
          write(6,*) 'norm', zstates_dotp(gr%m, psi%d%dim,  chi%zpsi(:,:, 1, 1), chi%zpsi(:,:, 1, 1))
        end if

        write(filename,'(a,i3.3)') 'PsiT.', ctr_iter
        call zoutput_function(sys%outp%how,'opt-control',filename,gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE, ierr) 
        
        if(in_debug_mode) write(6,*) 'penalty ',oct%tdpenalty(:,2)

        ! iteration managament ! make own subroutine
        call iteration_manager(stoploop)
        if (stoploop)  exit ctr_loop
        
         ! and now backward
        call target_calc(method,target_st,psi,chi)
        call bwd_step(method)
        call zoutput_function(sys%outp%how,'opt-control','last_bwd',gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE, ierr) 
        
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_zbr98


    ! ----------------------------------------------------------
    subroutine fwd_step(method)
      integer, intent(in) :: method

      call push_sub('opt_control.fwd_step')
      
      ! setup forward propagation
      call states_densities_init(psi, gr, sys%geo)
      call states_calc_dens(psi, NP_PART, dens_tmp)

      ! psi%rho = M_ZERO

      psi%rho = dens_tmp

      call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
      do i = 1, st%d%nspin
        v_old_i(:, i, 2) = h%vhxc(:, i)
      end do
      v_old_i(:, :, 3) = v_old_i(:, :, 2)

      ! new:
      v_old_i(:, :, 0) = v_old_i(:, :, 2)
      v_old_i(:, :, 1) = v_old_i(:, :, 2)
      do i = 1, 3
        v_old_f(:,:,i) = v_old_f(:,:,0) ! this one comes from the previous propagation
      end do

      message(1) = "Info: Propagating forward"
      call write_info(1)
!      call loct_progress_bar(-1, td%max_iter-1)

      h%ep%lasers(1)%dt = td%dt

      do i = 1, td%max_iter
  
        call prop_iter_fwd(i,method)
        if(oct%targetmode==oct_targetmode_td) &
          call calc_tdfitness(td_tg, tdtarget, td_fitness(i))     
        ! if td_target
        ! call loct_progress_bar(i-1, td%max_iter-1)
      end do
      
      ! dump new laser field
      if(dump_intermediate) then
        write(filename,'(a,i3.3)') 'opt-control/laser.', ctr_iter
        call write_field(filename, laser, td%max_iter, td%dt)
      endif
      
      call pop_sub()
    end subroutine fwd_step
    
    
    ! --------------------------------------------------------
    subroutine bwd_step(method) 
      integer, intent(in) :: method

      call push_sub('opt_control.bwd_step')

      ! setup backward propagation
      call states_calc_dens(chi, NP_PART, dens_tmp)
      chi%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
      do i = 1, st%d%nspin
        v_old_f(:, i, 2) = h%vhxc(:, i)
      end do
      v_old_f(:, :, 3) = v_old_f(:, :, 2)
      
      do i = 2, 3
        v_old_i(:,:,i) = v_old_i(:,:,1) ! this one comes the previous propagation
      end do

      message(1) = "Info: Propagating backward"
      call write_info(1)
      !call loct_progress_bar(td%max_iter-1, td%max_iter-1)
      td%dt = -td%dt
      h%ep%lasers(1)%dt = td%dt
      do i = td%max_iter-1, 0, -1
        call prop_iter_bwd(i,method)
        !call loct_progress_bar(td%max_iter-1-i, td%max_iter-1)
      end do
      td%dt = -td%dt

      ! dump new laser field
      if(dump_intermediate) then
        write(filename,'(a,i3.3)') 'opt-control/b_laser.', ctr_iter
        call write_field(filename, laser_tmp, td%max_iter, td%dt)
      end if

      call pop_sub()
    end subroutine bwd_step


    ! ---------------------------------------------------------
    subroutine iteration_manager(stoploop)
      logical, intent(out) :: stoploop

      call push_sub('opt_control.iteration_manager')
      
      stoploop = .false.
      ! calculate overlap 
      message(1) = "Info: Calculate overlap"
      call write_info(1)
      call calc_J() ! psi mit target_st
      
      message(1) = "Info: Loop control"
      call write_info(1)
      
      ! TODO:: check for STOP FILE AND delete it

      if((ctr_iter .eq. oct%ctr_iter_max).or.(oct%eps>M_ZERO.and.abs(functional-old_functional) < oct%eps)) then
        if((ctr_iter .eq. oct%ctr_iter_max)) then
          message(1) = "Info: Maximum number of iterations reached"
          call write_info(1)
        endif
        if(oct%eps > M_ZERO .and. abs(functional-old_functional) < oct%eps ) then
          message(1) = "Info: Convergence threshold reached"
          call write_info(1)
        endif
        
        stoploop = .TRUE.
      end if

      write(message(1), '(a,i3)') 'Info: Optimal control iteration #', ctr_iter
      call write_info(1)

      if(oct%mode_fixed_fluence) then
        write(message(1), '(6x,a,f10.5,a,f10.5,a,f10.5,a,f10.5)') &
          " => J1:", overlap, "   J: " , functional,  "  I: " , fluence, &
          " penalty: ", oct%a_penalty(ctr_iter)
      else
        write(message(1), '(6x,a,f14.8,a,f20.8,a,f14.8)') &
          " => J1:", overlap, "   J: " , functional,  "  I: " , fluence
      end if
      call write_info(1)

      ! store field with best J
      if(functional.gt.bestJ) then
        bestJ          = functional
        bestJ_J1       = overlap
        bestJ_fluence  = fluence
        bestJ_ctr_iter = ctr_iter
        ! dump to disc
        write(filename,'(a)') 'opt-control/laser.bestJ'
        call write_field(filename, laser, td%max_iter, td%dt)
      end if
      
      ! store field with best J1
      if(overlap.gt.bestJ1) then
        bestJ1          = overlap
        bestJ1_J        = functional
        bestJ1_fluence  = fluence       
        bestJ1_ctr_iter = ctr_iter
        ! dump to disc
        write(filename,'(a)') 'opt-control/laser.bestJ1'
        call write_field(filename, laser, td%max_iter, td%dt)
      end if
      
      ctr_iter = ctr_iter + 1
      old_functional = functional
      
      call pop_sub()
    end subroutine iteration_manager

    
    ! ---------------------------------------------------------
    subroutine calc_tdfitness(tg, tdtarget, merit)
      type(td_target_t), intent(in) :: tg(:)
      CMPLX,             intent(in) :: tdtarget(:)
      FLOAT,             intent(out):: merit
      integer             :: jj, ik, p, dim

      call push_sub('opt_control.calc_tdfitness')

      merit = M_ZERO
      do jj=1, size(tg)
        if(tg(jj)%type.eq.oct_tgtype_local) then
          do ik = 1, psi_i%d%nik
            do p  = psi_i%st_start, psi_i%st_end
              do dim = 1, psi_i%d%dim
                merit = merit + &
                  zmf_integrate(gr%m, tdtarget(:)* &
                  abs(psi%zpsi(:,dim,ik,p))**2)
              end do
            end do
          end do
        else
          do ik = 1, psi_i%d%nik
            do p  = psi_i%st_start, psi_i%st_end
              do dim = 1, psi_i%d%dim
                merit = merit + &
                  abs(zmf_integrate(gr%m, tdtarget(:)* &
                  conjg(psi%zpsi(:,dim,ik,p))))**2
              end do
            end do
          end do
        end if
      end do

      call pop_sub()
    end subroutine calc_tdfitness

    ! ---------------------------------------------------------
    ! calculate chi = \hat{O} psi
    ! do loop <target_st|Psi> for all States
    subroutine target_calc(method, targetst, psi_in, chi_out)
      integer,           intent(in)  :: method
      type(states_t),    intent(in)  :: targetst, psi_in
      type(states_t),    intent(inout) :: chi_out
      
      CMPLX   :: olap
      integer :: ik, p, dim
    
      call push_sub('opt_control.target_calc')

      if(oct%targetmode==oct_targetmode_static) then
        if(totype.eq.oct_tg_local) then ! only zr98 and wg05
          do ik = 1, psi_i%d%nik
            do p  = psi_i%st_start, psi_i%st_end
              do dim = 1, psi_i%d%dim
                ! multiply orbtials with local operator
                ! FIXME: for multiple particles 1,1,1 -> dim,p,ik
                chi_out%zpsi(:,dim,p,ik) = targetst%zpsi(:, 1, 1 , 1)*psi_in%zpsi(:, dim, p, ik)
              end do
            end do
          end do
        else ! totype nonlocal (all other totypes)
          do ik = 1, psi_i%d%nik
            do p  = psi_i%st_start, psi_i%st_end
              olap = M_z0
              do dim = 1, psi_i%d%dim
                olap = olap + zmf_integrate(gr%m,conjg(targetst%zpsi(:, dim, p, ik))*psi_in%zpsi(:, dim, p, ik))
              end do
              if(method == oct_algorithm_zr98) &
                chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:, :, p, ik)
              if(method == oct_algorithm_zbr98) then
                chi_out%zpsi(:,:,p,ik) = targetst%zpsi(:, :, p, ik)
              end if
              if(method == oct_algorithm_wg05) &
                chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:, :, p, ik)
            end do
          end do
        end if
      else
        ! time-dependent target
        message(1) = 'Info: Time-dependent Target selected'
        call write_info(1)
        chi_out%zpsi = M_z0
      end if

!      do ik = 1, psi%d%nik
!        do p  = psi%st_start, psi%st_end
!          olap = M_z0     
!          olap = zstates_dotp(gr%m, psi_in%d%dim,  targetst%zpsi(:,:, p, ik), psi_in%zpsi(:,:, p, ik))
!          if(method == 'ZR98') &
!            chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:,:,p,ik)
!          if(method == 'ZBR98') &
!            chi_out%zpsi(:,:,p,ik) = targetst%zpsi(:,:,p,ik)
!          if(method == 'WG05') &
!            chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:,:,p,ik)
!        end do
!      end do
      
      call pop_sub()
   end subroutine target_calc

   
   ! ---------------------------------------------------------
   subroutine prop_iter_fwd(iter,method)
     integer, intent(in) :: method
     integer,           intent(in) :: iter
     
     call push_sub('opt_control.prop_iter_fwd')
     
     ! chi(0) --> chi(T) [laser_tmp]
     !         |------------laser     
     ! psi(0) --> psi(T) 
     
     ! new electric field
     if(oct%targetmode==oct_targetmode_td) then
       ! psi(0) --> psi(T) with old field (for chi)
       call calc_inh(psi2,td_tg,iter,&
         td%dt,chi)
       ! psi2
       h%ep%lasers(1)%numerical => laser_tmp
       td%tr%v_old => v_old_i
       h%ep%lasers(1)%numerical => laser
       call states_calc_dens(psi2, NP_PART, dens_tmp)
       psi2%rho = dens_tmp
       call v_ks_calc(gr, sys%ks, h, psi2, calc_eigenval=.true.)
       call td_rti_dt(sys%ks, h, gr, psi2, td%tr, abs(iter*td%dt), abs(td%dt))
     end if
     
     call update_field(iter-1, laser, method)
     
     ! chi
     td%tr%v_old => v_old_f
     h%ep%lasers(1)%numerical => laser_tmp
     call states_calc_dens(chi, NP_PART, dens_tmp)
     chi%rho = dens_tmp
     call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.) 
     call td_rti_dt(sys%ks, h, gr, chi, td%tr, abs(iter*td%dt), abs(td%dt))
     
     ! psi
     td%tr%v_old => v_old_i
     h%ep%lasers(1)%numerical => laser
     call states_calc_dens(psi, NP_PART, dens_tmp)
     psi%rho = dens_tmp
     call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
     call td_rti_dt(sys%ks, h, gr, psi, td%tr, abs(iter*td%dt), abs(td%dt))
     ! write(6,*) iter, psi%zpsi(45:50,1,1,1)
     ! write(6,*) laser(2*iter,1)
     
     call pop_sub()
   end subroutine prop_iter_fwd
   
   
   ! ---------------------------------------------------------
   ! do backward propagation step with update of field
   ! works for time-independent and td targets
   subroutine prop_iter_bwd(iter, method)
     integer, intent(in) :: iter
     integer, intent(in) :: method
     
     call push_sub('opt_control.prop_iter_bwd')
     ! chi(T) --> chi(0)
     !         |------------laser_tmp
     ! psi(T) --> psi(0) [laser]
     
     ! calc inh first, then update field
     if(oct%targetmode==oct_targetmode_td) &
       call calc_inh(psi,td_tg,iter, &
       td%dt,chi)
     ! new electric field
     call update_field(iter+1, laser_tmp, method)!, tdpenalty(:,iter*2))
          
     ! chi
     td%tr%v_old => v_old_f
     h%ep%lasers(1)%numerical => laser_tmp   
     call states_calc_dens(chi, NP_PART, dens_tmp)
     chi%rho = dens_tmp
     call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
     call td_rti_dt(sys%ks, h, gr, chi, td%tr, abs(iter*td%dt),     td%dt )
     
     ! psi
     td%tr%v_old => v_old_i
     h%ep%lasers(1)%numerical => laser
     call states_calc_dens(psi, NP_PART, dens_tmp)
     psi%rho = dens_tmp
     call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
     call td_rti_dt(sys%ks, h, gr, psi, td%tr, abs(iter*td%dt),     td%dt )
     
     call pop_sub()
   end subroutine prop_iter_bwd
   
   
   ! ---------------------------------------------------------------
   subroutine calc_inh(psi_n, tg, iter, dt, chi_n)
     FLOAT,           intent(in)     :: dt
     integer,         intent(in)     :: iter
     type(states_t),  intent(in)     :: psi_n
     type(states_t),  intent(inout)  :: chi_n
     type(td_target_t), intent(in)   :: tg(:)
     
     CMPLX               :: tgt(NP_PART)
     integer             :: jj, dim
     CMPLX               :: olap

     call push_sub('opt_control.calc_inh')
     
     ! tdtarget is defined globally
     !      CMPLX               :: tdop(gr%m%np)
     
     ! distinguish between local and wavefunction td target
     ! TODO: change to right dimensions
     !       build target operator
     
     ! FIXME: more flexibility when defining multiple target operators
     ! usually one has only one specie of operators, i.e., only local or only non-local, 
     ! when mixing these, this routine has to be improved.
     if(tg(1)%type.eq.oct_tgtype_local) then
       ! local td operator
       ! build operator
       !do ip = 1, gr%m%np
       !   x = gr%m%x(ip, :)
       !   r = sqrt(sum(x(:)**2))
       
       !   call loct_parse_expression(psi_re, psi_im, &
       !        x(1), x(2), x(3), r, tt, td_tg(jj)%expression)
       
       !   ! need tdtarget for fitness calculation
       ! tdtarget(:) = (psi_re+m_zI*psi_im)
       do jj=1, size(tg)       
         call build_tdtarget(tg(jj),tgt, iter) ! tdtarget is build
         tdtarget(:) = tdtarget(:) + tgt(:) 
       end do
       do dim=1, psi_i%d%dim
         chi_n%zpsi(:,dim,1,1) = chi_n%zpsi(:,dim,1,1) &
           - sign(M_ONE,dt)/real(td%max_iter)*tdtarget(:)* &
           psi_n%zpsi(:,dim,1,1)
       enddo
       ! call zoutput_function(sys%outp%how,'opt-control','inh',gr%m,gr%sb,tdtarget,M_ONE, ierr) 
     else 
       ! build operator 
       ! non-local td operator
       !do ip = 1, gr%m%np
       !   x = gr%m%x(ip, :)
       !   r = sqrt(sum(x(:)**2))
       !   call loct_parse_expression(psi_re, psi_im, &
       !        x(1), x(2), x(3), r, tt, td_tg(jj)%expression)
       !   tdtarget(ip) = (psi_re+m_zI*psi_im)  
       ! need tdtarget for fitness calculation
       do jj=1, size(tg) 
         call build_tdtarget(tg(jj), tgt, iter)
         tdtarget = tgt
         olap= m_z0
         do dim=1, psi_i%d%dim
           olap = zmf_integrate(gr%m,psi_n%zpsi(:,dim,1,1)*&
             conjg(tdtarget(:)))
           chi_n%zpsi(:,dim,1,1) = chi_n%zpsi(:,dim,1,1) &
             - sign(M_ONE,dt)/real(td%max_iter)*olap*tdtarget(:)
         end do
       end do
     end if
     
     call pop_sub()
   end subroutine calc_inh
   

   ! ---------------------------------------------------------
   subroutine build_tdtarget(tdtg, tgt, iter)
     type(td_target_t), intent(in)   :: tdtg
     integer,           intent(in)   :: iter
     CMPLX  ,           intent(out)  :: tgt(:)

     integer :: pol
     call push_sub('opt_control.build_tdtarget')
     
     select case(tdtg%ftype)
     case(1)
       ! gauss
       ! Product x = gr%m%x(ip, :) 
       tgt = M_ONE
       do pol=1, NDIM
         tgt = tgt * exp( - (gr%m%x(:, pol) - &
           real(tdtg%tdshape(pol,iter),REAL_PRECISION ))**2/(M_TWO*tdtg%width**2))
       enddo
       
       !case(2)
       ! tgt = M-zero
       !   ! step
       !   Nwidth = tg(jj)%width
       !   Nwmin  = tg(jj)%numerical
       !   Nwmax  = 
       !   do kk=1,  gr%m%np  tg(jj)%width  tg(jj)%tdshape
       !   
       !   end do
     case default
       message(1)="Unknown TD Target Operator: Choose Gaussian(1) or Step(2)"
       call write_fatal(1)
     end select
     
      call pop_sub()
    end subroutine build_tdtarget


    ! ---------------------------------------------------------
    subroutine update_field(iter, l, method)
      integer, intent(in)      :: iter
      FLOAT,   intent(inout)   :: l(1:NDIM,0:2*td%max_iter)
      integer, intent(in)      :: method
      
      CMPLX :: d1
      CMPLX :: d2(NDIM)
      integer :: ik, p, dim, i, pol
      
      call push_sub('opt_control.update_field')
      
      ! case SWITCH
      ! check dipole moment function - > velocity gauge
      
      d1 = M_z0 
      d2 = M_z0 
      do ik = 1, psi_i%d%nik
        do p  = psi_i%st_start, psi_i%st_end
          do dim = 1, psi_i%d%dim
            do pol=1, NDIM
              d2(pol) = d2(pol) + zmf_integrate(gr%m,&
                conjg(chi%zpsi(:, dim, p, ik))*oct%laser_pol(pol)&
                *gr%m%x(:,pol)*psi%zpsi(:, dim, p, ik))
            enddo
            if(method.eq.oct_algorithm_zbr98) then
              d1 = zmf_integrate(gr%m,conjg(psi%zpsi(:, dim, p, ik))&
                *chi%zpsi(:, dim, p, ik))
              
            else
              d1 = M_z1
            end if
          end do
        end do
      end do
      ! Q: How to distinguish between the cases ?
      !if((NDIM-dof).eq.1) then
      !l(1:NDIM, 2*iter) = M_ONE/tdpenalty(1:NDIM,2*iter)*real(m_z1/(m_z2*m_zI) * &
      !(laser_pol(1:NDIM)*(d1*d2(1:NDIM)+conjg(d1)*conjg(d2(1:NDIM)))))
      !endif
      !if((NDIM-dof).eq.0) then
      l(1:NDIM, 2*iter) = aimag(d1*d2(1:NDIM))/oct%tdpenalty(1:NDIM,2*iter)
      !endif
      ! extrapolate to t+-dt/2
      i = int(sign(M_ONE, td%dt))
      if(iter==0.or.iter==td%max_iter) then
        l(1:NDIM, 2*iter+  i) = l(1:NDIM, 2*iter)
        l(1:NDIM, 2*iter+2*i) = l(1:NDIM, 2*iter)
      else
        l(1:NDIM, 2*iter+  i) = M_HALF*(M_THREE*l(1:NDIM, 2*iter) -       l(1:NDIM, 2*iter-2*i))
        l(1:NDIM, 2*iter+2*i) = M_HALF*( M_FOUR*l(1:NDIM, 2*iter) - M_TWO*l(1:NDIM, 2*iter-2*i))
      end if

      call pop_sub()
    end subroutine update_field


    ! ---------------------------------------------------------
    subroutine propagate_backward(psi_n) ! give initial chi and laser
      type(states_t), intent(inout)  :: psi_n
      call push_sub('opt_control.propagate_backward')
      
      message(1) = "Info: Backward propagating Chi"
      call write_info(1)
      ! setup the hamiltonian
      call states_calc_dens(psi_n, NP_PART, dens_tmp)
      psi_n%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)

      ! setup start of the propagation
      do i = 1, st%d%nspin
        v_old_f(:, i, 2) = h%vhxc(:, i)
      end do
      v_old_f(:, :, 3) = v_old_f(:, :, 2)

      h%ep%lasers(1)%numerical => laser
      td%tr%v_old => v_old_f

      td%dt = -td%dt
      h%ep%lasers(1)%dt = td%dt

      do i = td%max_iter-1, 0, -1
        ! time iterate wavefunctions
        call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), td%dt)
        ! update
        call states_calc_dens(psi_n, NP_PART, dens_tmp)
        psi_n%rho = dens_tmp
        call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
        
      end do
      td%dt = -td%dt
      
      call pop_sub()
    end subroutine propagate_backward


    ! ---------------------------------------------------------
    subroutine propagate_forward(psi_n) ! give initial psi and laser
      type(states_t), intent(inout)  :: psi_n

      integer :: ierr, ii
      FLOAT   :: etime

      call push_sub('opt_control.propagate_forward')

      message(1) = "Info: Forward propagating Psi"
      call write_info(1)

      ! setup the hamiltonian
      call states_calc_dens(psi_n, NP_PART, dens_tmp)
      psi_n%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
      ! setup start of the propagation
      do i = 1, st%d%nspin
        v_old_i(:, i, 2) = h%vhxc(:, i)
      end do
      v_old_f(:, :, 3) = v_old_i(:, :, 2)

      h%ep%lasers(1)%numerical => laser
      td%tr%v_old => v_old_i

      h%ep%lasers(1)%dt = td%dt      
      !call loct_progress_bar(-1, td%max_iter-1)

      ii = 1

      etime = loct_clock()
      do i = 1, td%max_iter 
        ! time iterate wavefunctions
        call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), abs(td%dt))

        ! if td_target
        if(oct%targetmode==oct_targetmode_td) &
          call calc_tdfitness(td_tg, tdtarget, td_fitness(i))
        
        ! update
        call states_calc_dens(psi_n, NP_PART, dens_tmp)
        psi_n%rho = dens_tmp
        call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
        call hamiltonian_energy(h, sys%gr, sys%geo, psi_n, -1)

        ! only write in final run
        if(oct%oct_double_check) then

          write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
            i*td%dt       / units_out%time%factor, &
            (h%etot + sys%geo%kinetic_energy) / units_out%energy%factor, &
            loct_clock() - etime
          call write_info(1)
          etime = loct_clock()

          ii = ii + 1 
          if(ii==sys%outp%iter+1 .or. i == td%max_iter) then ! output 
            if(i == td%max_iter) sys%outp%iter = ii - 1 
            ii = 1 
            call td_write_iter(write_handler, gr, psi_n, h, sys%geo, td%kick, td%dt, i)
            call td_write_data(write_handler, gr, psi_n, h, sys%outp, sys%geo, td%dt, i) 
          end if

        end if

        ! call loct_progress_bar(i-1, td%max_iter-1)
      end do
      
      call zoutput_function(sys%outp%how,'opt-control','last_fwd',gr%m,gr%sb,psi_n%zpsi(:,1,1,1),M_ONE, ierr) 
      
      message(1) = ""
      call write_info(1)
      call pop_sub()
    end subroutine propagate_forward


    ! ---------------------------------------------------------
    subroutine calc_fluence(laserin, fluenceout)
      FLOAT, intent(in) :: laserin(:,:)
      FLOAT, intent(out):: fluenceout
      call push_sub('opt_control.calc_fluence')

      fluenceout = SUM(laserin**2) * abs(td%dt)/M_TWO      

      call pop_sub()
    end subroutine calc_fluence


    ! ---------------------------------------------------------
    subroutine calc_J()
      FLOAT :: J2
      call push_sub('opt_control.calc_j')

      call calc_fluence(laser,fluence)
      J2 = SUM(oct%tdpenalty * laser**2) * abs(td%dt)/M_TWO
      if(oct%targetmode==oct_targetmode_td) then
         ! 1/T * int(<Psi| O | Psi>)
         overlap = SUM(td_fitness) / real(td%max_iter, REAL_PRECISION) 
      else
        overlap = abs(zstates_mpdotp(gr%m, psi, target_st))
        write(message(1), '(6x,a,f14.8)') " => overlap:", overlap
        call write_info(1)
      end if
      functional = overlap - J2

      convergence(1,ctr_iter) = functional
      convergence(2,ctr_iter) = overlap
      convergence(3,ctr_iter) = fluence
      convergence(4,ctr_iter) = oct%tdpenalty(1,1)

      call pop_sub()
    end subroutine calc_J


    ! ---------------------------------------------------------
    subroutine read_state(st, m, filename)
      type(states_t),   intent(inout) :: st
      type(mesh_t),     intent(in)  :: m
      character(len=*), intent(in)  :: filename

      integer :: ierr

      call push_sub('opt_control.read_state')

      call zinput_function('tmp/restart_gs/'//trim(filename), m, st%zpsi(:, 1, 1, 1), ierr, is_tmp=.TRUE.)
      ! if we do not succeed try obf
      if(ierr>0) call zinput_function ('tmp/restart_gs/'//trim(filename)//'.obf', &
        m, st%zpsi(:, 1, 1, 1), ierr, is_tmp=.TRUE.)
      ! if we do not succeed try NetCDF
      if(ierr>0) call zinput_function('tmp/restart_gs/'//trim(filename)//'.ncdf', &
        m, st%zpsi(:, 1, 1, 1), ierr, is_tmp=.TRUE.)

      if(ierr > 0) then
         message(1) = "Unsuccesfull read of states in 'tmp/restart_gs/" // trim(filename) // "'"
        call write_fatal(1)
      end if

      call pop_sub()
    end subroutine read_state

    
    ! ---------------------------------------------------------
    subroutine write_field(filename, las, steps, dt)
      integer,          intent(in) :: steps
      character(len=*), intent(in) :: filename
      FLOAT,            intent(in) :: las(1:NDIM,0:2*steps), dt ! ndim, step
      integer :: i, iunit
      FLOAT   :: tgrid(0:2*steps)

      call push_sub('opt_control.write_field')
      
      call t_lookup(2*steps+1,dt/real(2,REAL_PRECISION),tgrid)
      iunit = io_open(filename, action='write')
      do i = 0, 2*steps, 2
         write(iunit, '(4es30.16e4)') tgrid(i), las(:, i)
      end do
      call io_close(iunit)

      call pop_sub()
    end subroutine write_field


    ! ---------------------------------------------------------
    subroutine write_fieldw(filename, las)
      ! in w=(2pi f) space
      character(len=*), intent(in) :: filename
      FLOAT,            intent(in) :: las(1:NDIM,0:2*td%max_iter)
    
      integer :: i, iunit
      FLOAT   :: wgrid(0:2*td%max_iter)

      call push_sub('opt_control.write_fieldw')

      call w_lookup(2*td%max_iter+1,td%dt,wgrid)

      iunit = io_open(filename, action='write')
      do i = 0, 2*td%max_iter
         write(iunit, '(4es30.16e4)') wgrid(i), las(:, i)
      end do
      call io_close(iunit)
      
      call pop_sub()
    end subroutine write_fieldw


    ! ---------------------------------------------------------
    subroutine output()
      integer :: iunit, loop, ierr

      call push_sub('opt_control.output')
      
      iunit = io_open('opt-control/info', action='write')
      write(iunit, '(a,i4)')    'Total Iterations = ', ctr_iter
      write(iunit, '(a,f14.8)') 'Last Overlap    = ', overlap
      write(iunit, '(a,f14.8)') 'Last Functional = ', functional
      write(iunit, '(a)') 
      write(iunit, '(a)')       'Best value of functional'
      write(iunit, '(a,i4)')    'Iteration  = ', bestJ_ctr_iter
      write(iunit, '(a,f14.8)') 'Overlap    = ', bestJ_J1
      write(iunit, '(a,f14.8)') 'Functional = ', bestJ
      write(iunit, '(a,f14.8)') 'Fluence = ',    bestJ_fluence
      write(iunit, '(a)') 
      write(iunit, '(a)')       'Best value of target functional'
      write(iunit, '(a,i4)')    'Iteration  = ', bestJ1_ctr_iter
      write(iunit, '(a,f14.8)') 'Overlap    = ', bestJ1
      write(iunit, '(a,es18.8)') 'Functional = ', bestJ1_J
      write(iunit, '(a,f14.8)') 'Fluence = ',    bestJ1_fluence
      call io_close(iunit)
      message(1) = "Info: Output States"
      call write_info(1)
      ! should output wavefunctions ;)
      
      ! dump convergence: J,P,fluence,penalty
      iunit = io_open('opt-control/convergence', action='write')
      ! header
      write(iunit, '(4(a))') '# iteration ','functional ','overlap ','penalty '
      ! data
      do loop=1,oct%ctr_iter_max
         write(iunit, '(i6,3f18.8,es20.10)') loop, convergence(1,loop), convergence(2,loop), &
           convergence(3,loop),convergence(4,loop)
      end do
      call io_close(iunit)
      !do loop=1,ctr_iter_max
      !   write(6,'(f14.8,f14.8)') convergence(1,loop),convergence(2,loop)
      !end do
      ! assign final wave function to psi_i
      psi_i%zpsi = psi%zpsi


      ! check: sth goes wrong here: if psi_i replaced by psi 
      ! probably something is ill-defined
      call states_output(psi_i, gr, 'opt-control', sys%outp)
      call zoutput_function(sys%outp%how,'opt-control','target',gr%m,gr%sb,target_st%zpsi(:,1,1,1),M_ONE,ierr) 
      call zoutput_function(sys%outp%how,'opt-control','final',gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr) 

      call pop_sub()
    end subroutine output


    !------------------------------------------------------------------------
    subroutine def_istate()
      integer           :: kk, no_c, state, no_blk
      C_POINTER         :: blk
      integer           :: p, ik, ib, idim, inst, inik
      integer           :: id ,is, ip, ierr, no_states, isize
      character(len=10) :: fname
      type(states_t)    :: tmp_st 
      FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
      CMPLX             :: c_weight
      
      call push_sub('opt_control.def_istate')

      ! allocate tmp_state
      call states_init(tmp_st, gr, sys%geo) 
      if(h%ep%nvnl > 0) then
         ALLOCATE(tmp_st%rho_core(NP_PART), NP_PART)
         tmp_st%rho_core(NP_PART) = psi_i%rho_core(NP_PART)
      end if
      tmp_st%st_start = psi_i%st_start
      tmp_st%st_end = psi_i%st_end
      tmp_st%d%wfs_type = M_CMPLX
      isize = NP_PART*tmp_st%d%dim*(tmp_st%st_end-tmp_st%st_start+1)*tmp_st%d%nik
      ALLOCATE(tmp_st%zpsi(NP_PART, tmp_st%d%dim, tmp_st%st_start:tmp_st%st_end, tmp_st%d%nik), isize)

      !%Variable OCTInitialState
      !%Type integer
      !%Section Optimal Control
      !%Description
      !% The string OCTInitialState describes the initial state of the quantum system
      !% Possible arguments are:
      !%Option oct_is_groundstate 1
      !% start in the ground state 
      !%Option oct_is_excited 2
      !% start in the excited state given by OCTISnumber
      !% (ordered by energy)
      !%Option oct_is_superposition 3
      !% start in a superposition of states defined by the block OCTISsuperposition)
      !%Option oct_is_userdefined 4
      !% start in a userdefined state 
      !%End

      ! parse input
      call loct_parse_int(check_inp('OCTInitialState'), oct_is_groundstate, istype)
     
      ! make it work for single particle
      select case(istype)
      case(oct_is_groundstate) 
        message(1) =  'Info: Using Ground State for InitialState'
        call write_info(1)
        ! TODO: more particles
        call read_state(initial_st, gr%m, "0000000001")      
      case(oct_is_excited)  
        message(1) =  'Info: Using Excited State for InitialState'
        call write_info(1)

        !%Variable OCTTargetStateNumber
        !%Type integer
        !%Section Optimal Control
        !%Description
        !% Specify the target state, ordered by energy
        !% 1 corresponds to the ground state
        !% default is 2
        !%End
        call loct_parse_int(check_inp('OCTInitialStateNumber'), 2, state)
        
        write(fname,'(i10.10)') state
        call read_state(initial_st, gr%m, fname)
      case(oct_is_superposition)   
        message(1) =  'Info: Using Superposition of States for InitialState'
        call write_info(1)

        !%Variable OCTInitialSuperposition
        !%Type block
        !%Section Optimal Control
        !%Description
        !% The syntax of each line is, then:
        !%
        !% <tt>%OCTInitialSuperposition
        !% <br>&nbsp;&nbsp;state1 | weight1 | state2 | weight2 | ... 
        !% <br>%</tt>
        !%
        !% Note that weight can be complex, to produce current carrying superpositions.
        !% Example:
        !%
        !% <tt>%OCTInitialSuperposition
        !% <br>&nbsp;&nbsp;1 | 1 | 2 | -1 | ... 
        !% <br>%</tt>
        !%
        !%End 
        if(loct_parse_block(check_inp('OCTInitialSuperposition'),blk)==0) then
          no_blk = loct_parse_block_n(blk)
          ! TODO: for each orbital            
          do i=1, no_blk
            ! The structure of the block is:
            ! domain | function_type | center | width | weight 
            no_c = loct_parse_block_cols(blk, i-1)
            initial_st%zpsi(:,:,i,1) = M_z0
            do kk=1, no_c, 2
              call loct_parse_block_int(blk, i-1, kk-1, state)
              call loct_parse_block_cmplx(blk, i-1, kk, c_weight)
              write(fname,'(i10.10)') state
              !write(6,*) 'pars ', state, c_weight
              call read_state(tmp_st, gr%m, fname)
              initial_st%zpsi(:,:,1,1) = initial_st%zpsi(:,:,1,1) + c_weight * tmp_st%zpsi(:,:,1,1)
            enddo
          end do
          ! normalize state
          do ik = 1, psi%d%nik
            do p  = psi%st_start, psi%st_end
              call zstates_normalize_orbital(gr%m, psi_i%d%dim, &
                initial_st%zpsi(:,:, p, ik))
            enddo
          enddo
          !            end do
          
        end if
        
      case(oct_is_userdefined) 
        message(1) =  'Info: Building userdefined InitialState'
        call write_info(1)
        
        !%Variable OCTInitialUserdefined
        !%Type block
        !%Section Optimal Control
        !%Description
        !% 
        !% Example:
        !%
        !% <tt>%UserDefinedStates
        !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
        !% <br>%</tt>
        !%  
        !%End
        if(loct_parse_block(check_inp('OCTInitialUserdefined'),blk)==0) then
          
          no_states = loct_parse_block_n(blk)
          do ib = 1, no_states
            write(6,*) 'HELLO USERDEF' 
            call loct_parse_block_int(blk, ib-1, 0, idim)
            call loct_parse_block_int(blk, ib-1, 1, inst)
            call loct_parse_block_int(blk, ib-1, 2, inik)
            write(6,*) ' DEBUG: ', idim,inst,inik
            ! read formula strings and convert to C strings
            do id = 1, psi_i%d%dim
              do is = 1, psi_i%nst
                do ik = 1, psi_i%d%nik   
                  
                  ! does the block entry match and is this node responsible?
                  if(.not.(id.eq.idim .and. is.eq.inst .and. ik.eq.inik    &
                    .and. psi_i%st_start.le.is .and. psi_i%st_end.ge.is) ) cycle
                  
                  ! parse formula string
                  call loct_parse_block_string(                            &
                    blk, ib-1, 3, psi_i%user_def_states(id, is, ik))
                  ! convert to C string
                  call conv_to_C_string(psi_i%user_def_states(id, is, ik))
                  
                  do ip = 1, gr%m%np
                    x = gr%m%x(ip, :)
                    r = sqrt(sum(x(:)**2))
                    
                    ! parse user defined expressions
                    call loct_parse_expression(psi_re, psi_im,             &
                      x(1), x(2), x(3), r, M_ZERO, initial_st%user_def_states(id, is, ik))
                    ! fill state
                    !write(6,*) psi_re, psi_im
                    initial_st%zpsi(ip, id, is, ik) = psi_re + M_zI*psi_im
                  end do
                  ! normalize orbital
                  call zstates_normalize_orbital(gr%m, psi_i%d%dim, initial_st%zpsi(:,:, is, ik))
                end do
              end do
            enddo
          end do
          call loct_parse_block_end(blk)
          call zoutput_function(sys%outp%how,'opt-control','is_userdef',gr%m,gr%sb,initial_st%zpsi(:,1,1,1),u,ierr)
        else
          message(1) = '"UserDefinedStates" has to be specified as block.'
          call write_fatal(1)
        end if
        
      case default
        write(message(1),'(a)') "No valid initial state defined."
        write(message(2),'(a)') "Choosing the ground state."
        call write_info(2)
      end select
      
      if (in_debug_mode) call zoutput_function(sys%outp%how, &
        'opt-control', 'istate', gr%m, gr%sb, initial_st%zpsi(:,1,1,1), u, ierr)
      
      call pop_sub()
    end subroutine def_istate
    

    ! ----------------------------------------------------------------------
    subroutine def_toperator()
      integer           :: no_tds, no_c, state
      C_POINTER         :: blk
      integer           :: ierr, isize
      character(len=10) :: fname
      integer           :: ik, ib
      integer           :: no_states, kk, p, ip
      integer           :: no_blk, pol
      FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
      CMPLX             :: c_weight      
      type(states_t)    :: tmp_st
      character(len=1024) :: expression

      call push_sub('opt_control.def_toperator')

      ! allocate tmp_state
      call states_init(tmp_st, gr, sys%geo) 
      if(h%ep%nvnl > 0) then
         ALLOCATE(tmp_st%rho_core(NP_PART), NP_PART)
         tmp_st%rho_core(NP_PART) = psi_i%rho_core(NP_PART)
      end if

      tmp_st%st_start = psi_i%st_start
      tmp_st%st_end   = psi_i%st_end
      tmp_st%d%wfs_type = M_CMPLX
      isize = NP_PART*tmp_st%d%dim*(tmp_st%st_end-tmp_st%st_start+1)*tmp_st%d%nik
      ALLOCATE(tmp_st%zpsi(NP_PART, tmp_st%d%dim, tmp_st%st_start:tmp_st%st_end, tmp_st%d%nik), isize)

      ! prepare targetoperator 
      !%Variable OCTTargetOperator
      !%Type integer
      !%Section Optimal Control
      !%Description
      !% The string OCTTargetOperator describes the initial state of the quantum system
      !% Possible arguments are:
      !%Option *oct_tg_groundstate* 1 
      !% Targetoperator is a projection operator on the ground state
      !%Option *oct_tg_excited* 2
      !% Targetoperator is a projection operator on the excited state given by OCTTOnumber
      !%Option *oct_tg_superposition* 3
      !% Targetoperator is a projection operator on a superposition of states defined by the block OCTTOsuperposition)
      !%Option *oct_tg_userdefined* 4
      !% Targetoperator is a projection operator on a user defined state
      !%Option *oct_tg_local* 5
      !% Targetoperator is a local operator. Specify the shape and position within the block OCTLocalTarget.
      !%Option *oct_tg_td_local* 6
      !% Target operator is time-dependent, please specify block OCTTdTarget
      !%End
      call loct_parse_int(check_inp('OCTTargetOperator'),oct_tg_excited, totype)
      select case(totype)
      case(oct_tg_groundstate)
        message(1) =  'Info: Using Ground State for TargetOperator'
        call write_info(1)
        ! Target State is ground state
        ! rarely used: photo association
        ! maybe other cooling stuff
        ! TODO: expand to many particles
        call read_state(target_st, gr%m, "0000000001")
        
      case(oct_tg_excited) 
        message(1) =  'Info: Using Excited State for TargetOperator'
        call write_info(1)

        ! read in excited state
        !%Variable OCTTargetStateNumber
        !%Type integer
        !%Section Optimal Control
        !%Description
        !% Specify the target state, ordered by energy
        !% 1 corresponds to the ground state
        !% Default is 2
        !%End
        call loct_parse_int(check_inp('OCTTargetStateNumber'), 2, state)
        write(fname,'(i10.10)') state
        call read_state(target_st, gr%m, fname)
      case(oct_tg_superposition)  
        message(1) =  'Info: Using Superposition of States for TargetOperator'
        call write_info(1)

        !%Variable OCTTargetSuperposition
        !%Type block
        !%Section Optimal Control
        !%Description
        !% The syntax of each line is, then:
        !%
        !% <tt>%OCTTargetSuperposition
        !% <br>&nbsp;&nbsp;state1 | weight1 | state2 | weight2 | ... 
        !% <br>%</tt>
        !% 
        !% Example:
        !%
        !% <tt>%OCTTargetSuperposition
        !% <br>&nbsp;&nbsp;1 | 1 | 2 | -1 | ... 
        !% <br>%</tt>
        !%
        !%End
        if(loct_parse_block(check_inp('OCTTargetSuperposition'),blk)==0) then
          no_blk = loct_parse_block_n(blk)
          ! TODO: for each orbital            
          do i=1, no_blk
            ! The structure of the block is:
            ! domain | function_type | center | width | weight 
            no_c = loct_parse_block_cols(blk, i-1)
            target_st%zpsi(:,:,i,1) = M_z0
            do kk=1, no_c, 2
              call loct_parse_block_int(blk, i-1, kk-1, state)
              call loct_parse_block_cmplx(blk, i-1, kk, c_weight)
              write(fname,'(i10.10)') state
              !write(6,*) 'pars ', state, c_weight
              call read_state(tmp_st, gr%m, fname)
              target_st%zpsi(:,:,1,1) = target_st%zpsi(:,:,1,1) + c_weight * tmp_st%zpsi(:,:,1,1)
            enddo
          end do
          ! normalize state
          do ik = 1, psi%d%nik
            do p  = psi%st_start, psi%st_end
              call zstates_normalize_orbital(gr%m, psi_i%d%dim, &
                target_st%zpsi(:,:, p, ik))
            enddo
          enddo
          !            end do
          
        end if
        

      case(oct_tg_userdefined) 
        message(1) =  'Info: Using userdefined state for Targetoperator'
        call write_info(1)
        !%Variable OCTInitialUserdefined
        !%Type block
        !%Section Optimal Control
        !%Description
        !% 
        !% Example:
        !%
        !% <tt>%UserDefinedStates
        !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
        !% <br>%</tt>
        !%  
        !%End
      case(oct_tg_local) 
        message(1) =  'Info: Using Local Target'
        call write_info(1)
        !%Variable OCTLocalTarget
        !%Type block
        !%Section Optimal Control
        !%Description
        !% This block defines the shape and position of the local target operator. 
        !% It is possible to define several local operators which are then summed up.
        !% The syntax of each line is, then:
        !%
        !% <tt>%OCTLocalTarget
        !% <br>&nbsp;&nbsp; "exp(-r^2)*exp(-i*0.2*x)"
        !% <br>%</tt>
        !%  
        !% Example
        !%
        !% <tt>%OCTLocalTarget
        !% <br>&nbsp;&nbsp; "exp(-r^2)*exp(-i*0.2*x)"
        !% <br>%</tt>
        !%
        !%End
        ! parse input
        
        if(loct_parse_block(check_inp('OCTLocalTarget'),blk)==0) then
          no_states = loct_parse_block_n(blk)
          target_st%zpsi(:, 1, 1, 1) = m_z0
          do ib = 1, no_states
            ! parse formula string
            call loct_parse_block_string(blk, ib-1, 0, expression)
            ! convert to C string
            call conv_to_C_string(expression)
            
            do ip = 1, gr%m%np
              x = gr%m%x(ip, :)
              r = sqrt(sum(x(:)**2))
              
              ! parse user defined expressions
              call loct_parse_expression(psi_re, psi_im, &
                x(1), x(2), x(3), r, M_ZERO, expression)
              
              ! fill state
              target_st%zpsi(ip, 1, 1, 1) = target_st%zpsi(ip, 1, 1, 1) &
                + psi_re + M_zI*psi_im
            end do
            
            ! normalize orbital
            call zstates_normalize_orbital(gr%m, psi_i%d%dim, &
              target_st%zpsi(:,:, 1, 1))
          end do
          call loct_parse_block_end(blk)
          call zoutput_function(sys%outp%how,'opt-control','tg_local',&
            gr%m,gr%sb,target_st%zpsi(:,1,1,1),u,ierr)
        else
          message(1) = '"OCTLocalTarget" has to be specified as block.'
          call write_fatal(1)
        end if
      case(oct_tg_td_local)
        target_st%zpsi(:,1,1,1) = m_z0
        
      case default
        write(message(1),'(a)') "Target Operator not properly defined."
        call write_fatal(1)
      end select
      

      !%Variable OCTTdTarget
      !%Type block
      !%Section Optimal Control
      !%Description
      !% octopus features also time-dependent targets, i.e., one wants to optimize a laser 
      !% field that achieves a predefined time-dependent target. An example, could be the 
      !% evolution of occupation numbers in time. A time-dependent target consists of two 
      !% parts, i.e., the operator itself (a projection or a local operator) and its 
      !% time-dependence. Both are specified in one row of the block. You may enter as many 
      !% rows as you like. For time-dependent occupation targets OCTOPUS takes care of the 
      !% normalization.
      !% 
      !% The structure of the block is as follows:
      !%
      !% <tt>%OCTTdTarget
      !% <br>&nbsp;&nbsp;type | ftype | width | weight | g_x(t) | g_y(t) | g_z(t) 
      !% <br>%</tt>
      !%  
      !% Type:
      !% Choose betweem local and state target. 
      !% *oct_tgtype_local*
      !% *oct_tgtype_state*
      !%
      !% Ftype: 
      !% Spatial function of target operator.
      !% *oct_ftype_gauss*
      !%
      !% width:
      !% width of the Gaussian 
      !%
      !% weight:
      !% If multiple operators ae defined weight the importance between them
      !% g_x(t), g_y(t), g_z(t):
      !% describe the time-dependence
      !% 
      !% Example: 
      !% R0=1.0
      !% <tt>%OCTTdTarget
      !% <br>&nbsp;&nbsp;oct_tgtype_local | oct_ftype_gauss | 20 | 1 | "R0*sin(1.0*t)" | "R0*cos(1.0*t)" | "0"
      !% <br>%</tt>
      !% 
      !% In this example we define a narrow Gaussian which circles with radius R0 around the center. 
      !% 
      !%End
      if((oct%targetmode==oct_targetmode_td) &
        .AND.(loct_parse_block(check_inp('OCTTdTarget'),blk)==0)) then
        
        no_tds = loct_parse_block_n(blk)
        ALLOCATE(td_tg(no_tds), no_tds)
        do i=1, no_tds
          do pol=1, MAX_DIM
            td_tg(i)%expression(pol) = " "
          end do
          ! The structure of the block is:
          ! domain | function_type | center | width | weight 
          no_c = loct_parse_block_cols(blk, i-1)
          !td_tg(i)%type = oct_tgtype_local
          call loct_parse_block_int(blk, i-1, 0, td_tg(i)%type)
          call loct_parse_block_int(blk, i-1, 1, td_tg(i)%ftype)
          call loct_parse_block_float(blk, i-1, 2, td_tg(i)%width)
          call loct_parse_block_float(blk, i-1, 3, td_tg(i)%weight)
          do pol=1, NDIM
            call loct_parse_block_string(blk, i-1, 3+pol, &
              td_tg(i)%expression(pol))
          end do
          !
          ALLOCATE(td_tg(i)%tdshape(NDIM,0:2*td%max_iter), NDIM*(2*td%max_iter+1))
          call build_tdshape(td_tg(i), td%max_iter, td%dt)
          
        end do
        
      end if
      ! calc norm, give warning when to small
      
      call pop_sub()
    end subroutine def_toperator


    !----------------------------------------------------------
    subroutine build_tdshape(td_tg, steps, dt)
      type(td_target_t), intent(inout) :: td_tg
      integer,           intent(in)    :: steps
      FLOAT,             intent(in)    :: dt

      FLOAT   :: tgrid(0:2*steps)
      FLOAT   :: f_re, f_im
      integer :: kk, pol

      call push_sub('opt_control.build_tdshape')

      call t_lookup(2*steps+1,dt/real(2,REAL_PRECISION),tgrid)
      do kk=0, 2*steps
        do pol=1,NDIM
          f_re = M_ZERO
          f_im = M_ZERO
          !call loct_parse_expression(f_re, f_im, "t", tgrid(kk), &
          !     td_tg%expression(pol))
          call loct_parse_expression(f_re, f_im, & 	
            M_ZERO, M_ZERO, M_ZERO, M_ZERO, tgrid(kk), &
            td_tg%expression(pol)) 
          td_tg%tdshape(pol,kk) = f_re + M_zI*f_im
        end do
      end do
      call pop_sub()
    end subroutine build_tdshape
    
    
    !---------------------------------------------------------------
    subroutine init_tdtarget()
      integer             :: jj
      FLOAT               :: mxloc(MAX_DIM)
      integer             :: tt, pos(1)
      FLOAT               :: t
      CMPLX               :: tgt(NP_PART)
      call push_sub('opt_control.init_tdtarget')

      if(.not.associated(td_tg)) then
         call pop_sub(); return
      end if

      ! td target fitness
      ALLOCATE(td_fitness(0:td%max_iter), td%max_iter+1)
      td_fitness = M_ZERO
      ! to avoid double parsing and build up tdtarget array
      ALLOCATE(tdtarget(1:NP_PART), NP_PART)
      tdtarget   = M_z0

      ! generate some output to analyse the defined target
      do jj=1, size(td_tg)
         write(message(1),'(a,i2.2)') 'Info: Dumping td target ', jj
         call write_info(1)
         write(filename,'(a,i2.2)') 'opt-control/td_target_',jj
         iunit = io_open(filename, action='write')
         ! build target_function and calc inhomogeneity
         do tt=0, td%max_iter, 20 
            t = real(tt,REAL_PRECISION)*td%dt
            !write(6,*) tt, trim(td_tg(jj)%expression(1))
            !write(6,*) tt, trim(td_tg(jj)%expression(2))
            !write(6,*) tt, trim(td_tg(jj)%expression(3))
            call build_tdtarget(td_tg(jj),tgt,tt)
            tdtarget = tgt
            pos = maxloc(abs(tdtarget))
            !write(6,*) pos
            mxloc(:) = gr%m%x(pos(1),:) 
            !write(6,*) mxloc
            write(iunit, '(4es30.16e4)') t, mxloc, td_tg(jj)%tdshape(1,tt)
         end do
         close(iunit)
      end do
    
      call pop_sub()
    end subroutine init_tdtarget


    ! ---------------------------------------------------------
    subroutine init_()
      integer :: ierr, kk, jj

      call push_sub('opt_control.init')  

      call io_mkdir('opt-control')
      ! some shortcuts
      gr  => sys%gr
      st  => sys%st

      ! prepare unit factor u for output
      u  = M_ONE/units_out%length%factor**NDIM

      call td_init(gr, td, sys%st, sys%outp)
      call td_write_init(write_handler, gr, sys%st, sys%geo, (td%move_ions>0), h%ep%with_gauge_field, td%iter, td%dt)
      call states_allocate_wfns(st, gr%m, M_CMPLX)

      ALLOCATE(dens_tmp(gr%m%np_part, st%d%nspin), gr%m%np_part*st%d%nspin) 

      ! psi_i is initialized in system_init
      psi_i = st

      No_electrons = SUM(psi_i%occ(1,:))
      write(message(1),'(a,i4)') 'Info: Number of electrons: ', &
        No_electrons
      call write_info(1)

      ! FIXME: is this check obsolete ?
      ! psi_i never really used
      !psi_i should have complex wavefunctions
      if (psi_i%d%wfs_type /= M_CMPLX) then
        message(1) = "error in init_.opt_control"
        call write_fatal(1)
      end if

      ! call write_info(2)
      v_old_i => td%tr%v_old

      psi = st
      psi2 = st
      chi = st
      initial_st = st
      target_st = st

      ! allocate memory
      ALLOCATE(v_old_f(NP, chi%d%nspin, 0:3), NP_PART*chi%d%nspin*(3+1))
      v_old_f = M_ZERO
      ALLOCATE(laser(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter+1))
      laser = M_ZERO
      ALLOCATE(laser_tmp(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter+1))
      laser_tmp = M_ZERO
      ALLOCATE(field(1:NDIM), NDIM)
      field = M_ZERO

      ! use same laser as time-dependent run mode
      call laser_init(h%ep%no_lasers,h%ep%lasers,gr%m)

      !
      bestJ           = M_ZERO
      bestJ1          = M_ZERO
      bestJ_ctr_iter  = M_ZERO
      bestJ1_ctr_iter = M_ZERO

      do jj = 0, 2*td%max_iter
        t = td%dt*(jj-1)/M_TWO
        !i = int(abs(M_TWO*t/l(1)%dt) + M_HALF)

        do kk=1, h%ep%no_lasers
          call laser_field(gr%sb,h%ep%no_lasers,h%ep%lasers,t,field)
          laser(:,jj) = laser(:,jj) + field
        end do
      enddo
      write(filename,'(a)') 'opt-control/initial_laser.1'
      call write_field(filename, laser, td%max_iter, td%dt)
      ! - deallocate multiple lasers
      call laser_end(h%ep%no_lasers,h%ep%lasers)
      ! - NOW: allocate just one laser for optimal control
      h%ep%no_lasers = 1
      ALLOCATE(h%ep%lasers(1), 1)
      ! ALLOCATE(h%ep%lasers(1)%numerical(1:NDIM,0:2*td%max_iter),NDIM*(2*td%max_iter +1) )
      h%ep%lasers(1)%envelope = 99 ! internal type
      h%ep%lasers(1)%dt = td%dt
      h%ep%lasers(1)%numerical => laser   
      h%ep%lasers(1)%numerical => laser_tmp
      write(filename,'(a)') 'opt-control/initial_laser.2'
      call write_field(filename, laser, td%max_iter, td%dt)
      write(message(1),'(a,f14.8)') 'Input: Fluence of Initial laser ', sum(laser**2)*abs(td%dt)
      call write_info(1)
      ! initial laser definition end

      ! call after laser is setup ! do not change order
      call oct_read_inp(oct)
      
      ! This should be temporal

      h%ep%lasers(1)%pol(:) = oct%laser_pol(:)

      ALLOCATE(oct%tdpenalty(NDIM,0:2*td%max_iter),NDIM*(2*td%max_iter+1))
      oct%tdpenalty = M_ONE
      call def_tdpenalty(gr, oct%tdpenalty, td%max_iter, td%dt, mode_tdpenalty)
      oct%tdpenalty = oct%tdpenalty * oct%penalty

      if(in_debug_mode) then
        write(filename,'(a)') 'opt-control/td_penalty'
        call write_field(filename, oct%tdpenalty, td%max_iter, td%dt)
      end if

      if(oct%filtermode.gt.0) then
        call def_filter(gr, td%max_iter, td%dt, oct%filtermode, f)
        !! DEBUG
        !write(6,*) 'filter defined', size(f)
        do kk=1, size(f)
          write(filename,'(a,i2.2)') 'opt-control/filter', kk
          if(f(kk)%domain.eq.1) then
            call write_fieldw(filename, real(f(kk)%numerical(:,:), REAL_PRECISION))
          else
            !write(6,*) f(kk)%numerical(:,0)
            iunit = io_open(filename, action='write')
            do i = 0, 2*td%max_iter
              write(iunit, '(4ES30.16E4)') i*td%dt*M_HALF, f(kk)%numerical(:,i)
            end do
            call io_close(iunit)
            !call write_field(filename, real(f(kk)%numerical(:,:), REAL_PRECISION), td%max_iter, td%dt)
          end if
          !!
        end do
      end if

      ! global fitness
      ALLOCATE(convergence(4,0:oct%ctr_iter_max),(oct%ctr_iter_max+1)*4)

      ! initial state
      ! target operator
      call def_istate()
      call def_toperator()
      call init_tdtarget()

      call check_faulty_runmodes()

      call zoutput_function(sys%outp%how,'opt-control','initial1',gr%m,gr%sb,initial_st%zpsi(:,1,1,1),u,ierr)
      call zoutput_function(sys%outp%how,'opt-control','target1',gr%m,gr%sb,target_st%zpsi(:,1,1,1),u,ierr) 

      call pop_sub()
    end subroutine init_
    
    
    ! ---------------------------------------------------------
    subroutine check_faulty_runmodes()

      integer :: jj
      call push_sub('opt_control.check_faulty_runmodes')

      ! try to avoid ill defined combinations of run modes
      ! be careful with the order !!
      
      ! FixedFluence and Filter only work with WG05
      if((oct%mode_fixed_fluence).and.(oct%algorithm_type.ne.oct_algorithm_wg05)) then
        write(message(1),'(a)') "Cannot optimize to a given fluence with the chosen algorithm."
        write(message(2),'(a)') "Switching to scheme WG05."         
        call write_info(2)
        oct%algorithm_type = oct_algorithm_wg05
      end if
      
      ! Filters work only with WG05
      if((oct%filtermode.gt.0).and.(oct%algorithm_type.ne.oct_algorithm_wg05)) then
        write(message(1),'(a)') "Warning: Cannot use filters with the chosen algorithm."
        write(message(2),'(a)') "Warning: Switching to scheme WG05."
        call write_info(2)
        oct%algorithm_type = oct_algorithm_wg05
      end if
      
      ! tdpenalty and fixed fluence do not work !
      if((mode_tdpenalty).AND.(oct%mode_fixed_fluence)) then
        write(message(1),'(a)') "Warning: Cannot use fixed fluence and" &
          //" td penalty."
        write(message(2),'(a)') "Warning: Disabling td penalty."
        call write_info(2)
        oct%tdpenalty = oct%penalty
      end if
      
      ! local targets only in ZR98 and WG05
      if((totype.eq.oct_tg_local) & 
        .AND.(oct%algorithm_type.eq.oct_algorithm_zbr98)) then
        write(message(1),'(a)') "Warning: Local targets work" &
          // " only with ZR98 and WG05."
        write(message(2),'(a)') "Warning: Switching to ZR98."
        call write_info(2)
        oct%algorithm_type = oct_algorithm_zr98
      end if
      
      ! tdtargets only in ZR98 and WG05
      if((oct%targetmode.eq.oct_targetmode_td) & 
        .AND.(oct%algorithm_type.eq.oct_algorithm_zbr98)) then
        write(message(1),'(a)') "Warning: Time-dependent targets work" &
          // " only with ZR98 and WG05."
        write(message(2),'(a)') "Warning: Please change algorithm type."
        call write_fatal(2)
      end if
      
      ! td target with states works only for 1 electron
      ! TODO: HOW TO RETRIEVE NUMBER OF ELECTRONS
      !       UNCOMMENT LAST LINES
      if(associated(td_tg)) then
        do jj=1, size(td_tg)
          if(td_tg(jj)%type.eq.oct_tgtype_state) &
            td_tg_state = .TRUE.
        end do
      end if
      ! occ in states_t
      !write(6,*) 'DEBUG: occ: ', shape(psi_i%occ), size(psi_i%occ)
      if( (td_tg_state) .AND. (SUM(psi_i%occ(1,:)).gt.1) ) then
      !if(td_tg_state) then
      !     if(psi_i%occ.gt.1) then
        write(message(1),'(a)') "Warning: Time-dependent state targets" &
          // " work only with one electron."
        write(message(2),'(a)') "Warning: Please change input file."
        call write_fatal(2)
      end if

      call pop_sub()      
    end subroutine check_faulty_runmodes


    ! ---------------------------------------------------------
    subroutine end_()
      call push_sub('opt_control.end_')
      ! TODO:: deallocate filters (if there were any)
      call filter_end(f)

      call td_end(td)
      call pop_sub()
    end subroutine end_

  end subroutine opt_control_run

#include "opt_control_read.F90"

end module opt_control_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
