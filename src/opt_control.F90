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
    FLOAT            :: eps
    integer          :: ctr_iter_max
    FLOAT            :: penalty
    FLOAT, pointer   :: a_penalty(:)
    FLOAT, pointer   :: tdpenalty(:, :)
    logical          :: mode_tdpenalty
    FLOAT            :: targetfluence
    logical :: mode_fixed_fluence
    integer :: targetmode
    integer :: filtermode
    integer :: totype
    CMPLX   :: laser_pol(MAX_DIM)
    integer :: dof
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

  type oct_iterator_t
    integer            :: ctr_iter
    FLOAT              :: functional
    FLOAT              :: old_functional
    FLOAT              :: overlap
    FLOAT, pointer     :: convergence(:,:)
    FLOAT              :: bestJ, bestJ1, bestJ_fluence, bestJ1_fluence
    FLOAT              :: bestJ_J1, bestJ1_J
    integer            :: bestJ_ctr_iter, bestJ1_ctr_iter
  end type oct_iterator_t

contains

  ! ---------------------------------------------------------
  subroutine opt_control_run(sys, h)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h

    type(oct_t)                :: oct
    type(oct_iterator_t)       :: iterator
    type(td_t)                 :: td
    type(grid_t),      pointer :: gr   ! some shortcuts
    type(states_t),    pointer :: st
    type(states_t)             :: chi, psi, target_st, initial_st, psi2
    type(filter_t),    pointer :: f(:)
    type(td_target_t), pointer :: td_tg(:)
    FLOAT, pointer     :: laser_tmp(:,:), laser(:,:)
    FLOAT, allocatable :: dens_tmp(:,:)
    FLOAT, allocatable :: td_fitness(:)
    CMPLX, allocatable :: tdtarget(:)
    integer :: ierr
    character(len=80)  :: filename

    call push_sub('opt_control.opt_control_run')

    ! Checks that the run is actually possible with the current settings
    call check_runmode_constrains(sys, h)

    ! initialize oct run, defines initial laser... initial state
    call init_()

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
    
    ! Output some info to stdout
    call output(oct, iterator)

    call states_output(psi, gr, 'opt-control', sys%outp)
    call zoutput_function(sys%outp%how, 'opt-control', 'target', gr%m, gr%sb, target_st%zpsi(:,1,1,1), M_ONE, ierr) 
    call zoutput_function(sys%outp%how, 'opt-control', 'final', gr%m, gr%sb, psi%zpsi(:,1,1,1), M_ONE, ierr) 

    ! do final test run: propagate initial state with optimal field
    call oct_finalcheck(oct, initial_st, target_st, sys, h, td, laser, td_tg, tdtarget)

    ! clean up
    nullify(h%ep%lasers(1)%numerical)
    deallocate(laser)
    deallocate(laser_tmp)
    call oct_iterator_end(iterator)
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
      call propagate_forward(oct, sys, h, td, laser, td_tg, tdtarget, td_fitness, psi) 
      
      if(in_debug_mode) call zoutput_function(sys%outp%how, &
        'opt-control', 'prop1', gr%m, gr%sb, psi%zpsi(:,1,1,1), M_ONE, ierr)
        
      call zoutput_function(sys%outp%how,'opt-control','zr98_istate1',gr%m,gr%sb,initial_st%zpsi(:,1,1,1),M_ONE,ierr)

      call oct_iterator_init(iterator, oct)
      
      ctr_loop: do
        
        ! check for stop file and delete file
        
        ! define target state
        call target_calc(oct, gr, target_st, psi, chi) ! defines chi
        
        if(iteration_manager(oct, gr, td_fitness, laser, td, psi, target_st, iterator)) exit ctr_loop
        
        call bwd_step(oct_algorithm_zr98)
        
        if(oct%dump_intermediate) then
          write(filename,'(a,i3.3)') 'opt-control/b_laser.', iterator%ctr_iter
          call write_field(filename, td%max_iter, NDIM, laser_tmp, td%dt)
        end if
        
        ! forward propagation
        psi = initial_st
        psi2 = initial_st
        call zoutput_function(sys%outp%how,'opt-control','last_bwd',gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr)
        call fwd_step(oct_algorithm_zr98)
        write(filename,'(a,i3.3)') 'PsiT.', iterator%ctr_iter
        call zoutput_function(sys%outp%how,'opt-control',filename,gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr)
        
        if(oct%dump_intermediate) then
          write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
          call write_field(filename, td%max_iter, NDIM, laser, td%dt)
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
      
      call oct_iterator_init(iterator, oct)

      ctr_loop: do
         
        ! first propagate chi to ti
        message(1) = "Info: Initial forward propagation"
        call write_info(1)
        
        psi = initial_st
        call propagate_forward(oct, sys, h, td, laser, td_tg, tdtarget, td_fitness, psi) 
        
        if(in_debug_mode) then
          write(filename,'(a,i3.3)') 'PsiT.', iterator%ctr_iter
          call zoutput_function(sys%outp%how,'opt-control',filename,gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr)
        end if

        ! define target state
        call target_calc(oct, gr, target_st, psi, chi) ! defines chi
        
        if(iteration_manager(oct, gr, td_fitness, laser, td, psi, target_st, iterator)) exit ctr_loop
        
        call bwd_step(oct_algorithm_wg05)
        !!!WARNING: this probably shoudl not be here?
        fluence = laser_fluence(laser_tmp, td%dt)
        
        ! filter stuff / alpha stuff
        if (oct%filtermode.gt.0) & 
          call apply_filter(gr, td%max_iter, oct%filtermode, f, laser_tmp)

        ! recalc field
        if (oct%mode_fixed_fluence) then
          !!!WARNING: this probably shoudl not be here?
          fluence = laser_fluence(laser_tmp, td%dt)
          if(in_debug_mode) write (6,*) 'actual fluence', fluence

          oct%a_penalty(iterator%ctr_iter + 1) = &
            sqrt(fluence * oct%a_penalty(iterator%ctr_iter)**2 / oct%targetfluence )

          if(in_debug_mode) then
            write (6,*) 'actual penalty', oct%a_penalty(iterator%ctr_iter)
            write (6,*) 'next penalty', oct%a_penalty(iterator%ctr_iter + 1)
          end if

          oct%tdpenalty = oct%a_penalty(iterator%ctr_iter + 1)
          laser_tmp = laser_tmp * oct%a_penalty(iterator%ctr_iter) / oct%a_penalty(iterator%ctr_iter + 1)
          fluence = laser_fluence(laser_tmp, td%dt)

          if(in_debug_mode) write (6,*) 'renormalized', fluence, oct%targetfluence

        end if

        laser = laser_tmp
        
        ! dump here: since the fwd_step is missing
        write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
        call write_field(filename, td%max_iter, NDIM, laser, td%dt)   
        
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_wg05

 
    ! ---------------------------------------------------------
    subroutine scheme_zbr98!(method)
      integer :: ierr
      call push_sub('opt_control.scheme_zbr98')
      
      message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
      call write_info(1)
      
      ! first propagate chi to T
      message(1) = "Info: Initial backward propagation"
      call write_info(1)
      
      call target_calc(oct, gr, target_st, psi, chi)
      call propagate_backward(sys, h, td, laser, chi)

      if(in_debug_mode) call zoutput_function(sys%outp%how, &
        'opt-control', 'initial_propZBR98', gr%m, gr%sb, chi%zpsi(:,1,1,1), M_ONE, ierr)

      call oct_iterator_init(iterator, oct)

      ctr_loop: do

        ! check for stop file and delete file
        message(1) = "Info: Setup forward"
        call write_info(1)
        
        ! forward propagation
        psi = initial_st
        
        call fwd_step(oct_algorithm_zbr98)

        if(in_debug_mode) then
          write(6,*) 'norm', zstates_dotp(gr%m, psi%d%dim,  psi%zpsi(:,:, 1, 1), psi%zpsi(:,:, 1, 1))
          write(6,*) 'norm', zstates_dotp(gr%m, psi%d%dim,  chi%zpsi(:,:, 1, 1), chi%zpsi(:,:, 1, 1))
        end if

        write(filename,'(a,i3.3)') 'PsiT.', iterator%ctr_iter
        call zoutput_function(sys%outp%how,'opt-control',filename,gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE, ierr) 
        
        if(in_debug_mode) write(6,*) 'penalty ',oct%tdpenalty(:,2)

        if(iteration_manager(oct, gr, td_fitness, laser, td, psi, target_st, iterator)) exit ctr_loop
        
         ! and now backward
        call target_calc(oct, gr, target_st, psi, chi)
        call bwd_step(oct_algorithm_zbr98)
        call zoutput_function(sys%outp%how,'opt-control','last_bwd',gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE, ierr) 
        
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_zbr98


    ! ----------------------------------------------------------
    subroutine fwd_step(method)
      integer, intent(in) :: method
      integer :: i

      call push_sub('opt_control.fwd_step')
      
      ! setup forward propagation
      call states_densities_init(psi, gr, sys%geo)
      call states_calc_dens(psi, NP_PART, dens_tmp)

      ! psi%rho = M_ZERO

      psi%rho = dens_tmp

      call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
      call td_rti_run_zero_iter(h, td%tr)

      message(1) = "Info: Propagating forward"
      call write_info(1)
!      call loct_progress_bar(-1, td%max_iter-1)

      h%ep%lasers(1)%dt = td%dt

      do i = 1, td%max_iter
  
        call prop_iter_fwd(i,method)
        if(oct%targetmode==oct_targetmode_td) &
          call calc_tdfitness(td_tg, gr, psi, tdtarget, td_fitness(i))     
        ! if td_target
        ! call loct_progress_bar(i-1, td%max_iter-1)
      end do
      
      ! dump new laser field
      if(oct%dump_intermediate) then
        write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
        call write_field(filename, td%max_iter, NDIM, laser, td%dt)
      endif
      
      call pop_sub()
    end subroutine fwd_step
    
    
    ! --------------------------------------------------------
    subroutine bwd_step(method) 
      integer, intent(in) :: method
      integer :: i

      call push_sub('opt_control.bwd_step')

      ! setup backward propagation
      call states_calc_dens(chi, NP_PART, dens_tmp)
      chi%rho = dens_tmp
      call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
      call td_rti_run_zero_iter(h, td%tr)

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
      if(oct%dump_intermediate) then
        write(filename,'(a,i3.3)') 'opt-control/b_laser.', iterator%ctr_iter
        call write_field(filename, td%max_iter, NDIM, laser_tmp, td%dt)
      end if

      call pop_sub()
    end subroutine bwd_step


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
       h%ep%lasers(1)%numerical => laser
       call states_calc_dens(psi2, NP_PART, dens_tmp)
       psi2%rho = dens_tmp
       call v_ks_calc(gr, sys%ks, h, psi2, calc_eigenval=.true.)
       call td_rti_dt(sys%ks, h, gr, psi2, td%tr, abs(iter*td%dt), abs(td%dt))
     end if
     
     call update_field(iter-1, laser, method)
     
     ! chi
     h%ep%lasers(1)%numerical => laser_tmp
     call states_calc_dens(chi, NP_PART, dens_tmp)
     chi%rho = dens_tmp
     call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.) 
     call td_rti_dt(sys%ks, h, gr, chi, td%tr, abs(iter*td%dt), abs(td%dt))
     
     ! psi
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
     h%ep%lasers(1)%numerical => laser_tmp   
     call states_calc_dens(chi, NP_PART, dens_tmp)
     chi%rho = dens_tmp
     call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
     call td_rti_dt(sys%ks, h, gr, chi, td%tr, abs(iter*td%dt),     td%dt )
     
     ! psi
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
       do dim=1, chi_n%d%dim
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
         do dim=1, psi_n%d%dim
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
      do ik = 1, psi%d%nik
        do p  = psi%st_start, psi%st_end
          do dim = 1, psi%d%dim
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
    subroutine init_tdtarget(td_tg)
      type(td_target_t), pointer :: td_tg(:)

      integer             :: i, jj, iunit
      FLOAT               :: mxloc(MAX_DIM)
      integer             :: tt, pos(1)
      FLOAT               :: t
      CMPLX               :: tgt(NP_PART)

      integer :: no_c, no_tds, pol
      C_POINTER         :: blk

      call push_sub('opt_control.init_tdtarget')

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
      no_tds = 0
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
  
      if(no_tds.eq.0) then
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
      integer :: ierr, kk, jj, iunit, i
      FLOAT, allocatable :: field(:)
      FLOAT :: t

      call push_sub('opt_control.init')  

      call io_mkdir('opt-control')

      ! Initially nullify the pointers
      nullify(f)
      nullify(td_tg)
      nullify(laser_tmp)
      nullify(laser)

      ! some shortcuts
      gr  => sys%gr
      st  => sys%st

      call td_init(gr, td, sys%st, sys%outp)
      call states_allocate_wfns(st, gr%m, M_CMPLX)

      ALLOCATE(dens_tmp(gr%m%np_part, st%d%nspin), gr%m%np_part*st%d%nspin) 

      ! Initialize a bunch of states: initial, target, auxiliary.
      psi        = st
      psi2       = st
      chi        = st
      initial_st = st
      target_st  = st

      ! allocate memory
      ALLOCATE(laser(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter+1))
      laser = M_ZERO
      ALLOCATE(laser_tmp(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter+1))
      laser_tmp = M_ZERO
      ALLOCATE(field(1:NDIM), NDIM)
      field = M_ZERO

      ! use same laser as time-dependent run mode
      call laser_init(h%ep%no_lasers, h%ep%lasers,gr%m)

      !
      do jj = 0, 2*td%max_iter
        t = td%dt*(jj-1)/M_TWO
        !i = int(abs(M_TWO*t/l(1)%dt) + M_HALF)

        do kk=1, h%ep%no_lasers
          call laser_field(gr%sb,h%ep%no_lasers,h%ep%lasers,t,field)
          laser(:,jj) = laser(:,jj) + field
        end do
      enddo
      write(filename,'(a)') 'opt-control/initial_laser.1'
      call write_field(filename, td%max_iter, NDIM, laser, td%dt)
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
      call write_field(filename, td%max_iter, NDIM, laser, td%dt)
      write(message(1),'(a,f14.8)') 'Input: Fluence of Initial laser ', sum(laser**2)*abs(td%dt)
      call write_info(1)
      deallocate(field)
      ! initial laser definition end

      ! call after laser is setup ! do not change order
      call oct_read_inp(oct)
      
      ! This should be temporal

      h%ep%lasers(1)%pol(:) = oct%laser_pol(:)

      ALLOCATE(oct%tdpenalty(NDIM,0:2*td%max_iter),NDIM*(2*td%max_iter+1))
      oct%tdpenalty = M_ONE
      call def_tdpenalty(gr, oct%tdpenalty, td%max_iter, td%dt, oct%mode_tdpenalty)
      oct%tdpenalty = oct%tdpenalty * oct%penalty

      if(in_debug_mode) then
        write(filename,'(a)') 'opt-control/td_penalty'
        call write_field(filename, td%max_iter, NDIM, oct%tdpenalty, td%dt)
      end if

      if(oct%filtermode.gt.0) then
        call def_filter(gr, td%max_iter, td%dt, oct%filtermode, f)
        !! DEBUG
        !write(6,*) 'filter defined', size(f)
        do kk=1, size(f)
          write(filename,'(a,i2.2)') 'opt-control/filter', kk
          if(f(kk)%domain.eq.1) then
            call write_fieldw(filename, NDIM, td%max_iter, real(f(kk)%numerical(:,:), REAL_PRECISION), td%dt)
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

      ! initial state
      ! target operator
      call def_istate(oct, gr, sys%outp, initial_st)
      call def_toperator(oct, gr, sys%outp, target_st)
      call init_tdtarget(td_tg)

      call check_faulty_runmodes(oct)

      call zoutput_function(sys%outp%how,'opt-control','initial1', gr%m, gr%sb, initial_st%zpsi(:,1,1,1), &
                            M_ONE/units_out%length%factor**NDIM, ierr)
      call zoutput_function(sys%outp%how,'opt-control','target1', gr%m, gr%sb, target_st%zpsi(:,1,1,1), &
                            M_ONE/units_out%length%factor**NDIM, ierr) 

      call pop_sub()
    end subroutine init_
    
    
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

end module opt_control_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
