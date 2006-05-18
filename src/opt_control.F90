!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
  use output_m
  use units_m
!  use filter_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use mesh_m
  use grid_m
  use states_m
  use system_m
  use restart_m
  use v_ks_m
  use hamiltonian_m
  use timedep_m
  use td_rti_m
  use mesh_function_m
  use lasers_m
  implicit none

  private
  public :: opt_control_run

  integer, parameter, private  ::  &
    oct_istate_gs = 1,             &
    oct_istate_ex = 2,             &
    oct_istate_sp = 3,             &
    oct_istate_ud = 4             
  
  integer, parameter, private  ::  &
    oct_algorithm_zbr98 = 1,       &
    oct_algorithm_zr98  = 2,       &
    oct_algorithm_wg05  = 3       


contains

  ! ---------------------------------------------------------
  subroutine opt_control_run(sys, h)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h

    type(td_t)                :: td
    type(grid_t),     pointer :: gr   ! some shortcuts
    type(states_t),   pointer :: st, psi_i
    type(states_t)            :: chi, psi, target_st, initial_st

    FLOAT, pointer :: v_old_i(:,:,:), v_old_f(:,:,:)
    FLOAT, pointer :: laser_tmp(:,:), laser(:,:)
    integer :: i, ctr_iter, ctr_iter_max
    integer :: targetmode, filtermode
    FLOAT :: eps, penalty, overlap, functional, old_functional
    FLOAT :: fluence, u, t
    FLOAT, allocatable :: convergence(:,:) , field(:)
    
    character(len=80) :: filename
    character(len=5)  :: algtype, totype
    integer           :: istype, algorithm_type

    call push_sub('opt_control.opt_control_run')

    ! CHECK:: only dipole approximation in length gauge
    !         only single particle (yet)
    if(h%gauge.ne.1) then
       write(message(1),'(a)') "So far only length gauge is supported..."
       call write_fatal(1)
    end if

    ! call write_debug_marker(1)
    
    ! CHECK if run mode is defined
    ! define feature_vector - > check if all options match

    ! initialize oct run, defines initial laser... initial state
    call init_()

    ! mode switcher
    select case(algorithm_type)
      
    case(oct_algorithm_zbr98)    ! alternative name: BFB
      call scheme_ZBR98(algtype) ! (0)  rabitz zbr98 type algorithm (states only)    

    case(oct_algorithm_zr98)     ! fbf
      call scheme_ZR98(algtype)  ! (1)  rabitz zr98  type algorithm 
      ! TODO: o generalize elements
      !       o time-dependent targets

    case(oct_algorithm_wg05)  ! FB
      call scheme_WG05(algtype)  ! (2) Werschnik, Gross Type: allows fixed fluence and filtering
      ! (11) td targets rabitz type
      ! (12) td targets wg     type (with filter)

    case DEFAULT
      write(message(1),'(a)') "Unknown choice of OCT algorithm."
      write(message(2),'(a)') "Choose: ZR98, ZBR98, WG05 ."
      call write_fatal(2)

    end select

    ! output: States, Lasers, Convergence
    call output()

    ! done ... clean up
    write(message(1), '(a,i3)') 'Info: Cleaning up..#', ctr_iter
    call write_info(1)
    ! clean up
    td%tr%v_old => v_old_i
    nullify(h%ep%lasers(1)%numerical)
    deallocate(laser)
    ! write(6,*) 'DEBUG3'    
    deallocate(laser_tmp)
    ! write(6,*) 'DEBUG4'

!deallocate(convergence)
!! PROBLEMS WITH DEALLOCATION

    !nullify(v_old_f)
    !deallocate(v_old_f)
!! ----
    ! write(6,*) 'DEBUG5'
    call states_end(psi_i)
    ! write(6,*) 'DEBUG6'
    call end_()
    ! write(6,*) 'DEBUG7'
    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine scheme_zr98(method)
      character (len=*), intent(in) :: method

      integer :: ierr
      logical :: stoploop
     
      call push_sub('opt_control.scheme_ZR98')

      message(1) = "Info: Starting OCT iteration using scheme: ZR98"
      call write_info(1)
      
      stoploop = .FALSE.
      ! first propagate chi to ti
      message(1) = "Info: Initial forward propagation"
      call write_info(1)
      
      psi = initial_st
      call prop_fwd(psi) 
!! DEBUG
      call zoutput_function(sys%outp%how,'opt-control','prop1',gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr)
!!stop

      old_functional = -CNST(1e10)
      ctr_iter = 0

      ctr_loop: do

         ! check for stop file and delete file

         ! define target state
         call target_calc('ZR98', target_st, psi, chi) ! defines chi

         call iteration_manager(stoploop)
         if (stoploop)  exit ctr_loop

         message(1) = "Info: Setup backward"
         call write_info(1)
!! CHECK: why do we need st%rho ?
         ! setup backward propagation
         call states_calc_dens(chi, NP_PART, st%rho)
         call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
         do i = 1, st%d%nspin
            v_old_f(:, i, 2) = h%vhxc(:, i)
         end do
         v_old_f(:, :, 3) = v_old_f(:, :, 2)
         
         do i = 2, 3
            v_old_i(:,:,i) = v_old_i(:,:,1) ! this one comes the previous propagation
         end do
         ! and now backward
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
         ! write(stdout, '(1x)')
         write(filename,'(a,i3.3)') 'opt-control/b_laser.', ctr_iter
         call write_field(filename, laser_tmp)

         ! setup forward propagation
         psi = initial_st  
         call states_calc_dens(psi, NP_PART, st%rho)
         call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
         do i = 1, st%d%nspin
            v_old_i(:, i, 2) = h%vhxc(:, i)
         end do
         v_old_i(:, :, 3) = v_old_i(:, :, 2)
         
         do i = 2, 3
            v_old_f(:,:,i) = v_old_f(:,:,1) ! this one comes from the previous propagation
         end do
         message(1) = "Info: Propagating forward"
         call write_info(1)
         !call loct_progress_bar(-1, td%max_iter-1)
         h%ep%lasers(1)%dt = td%dt
         do i = 1, td%max_iter
            call prop_iter_fwd(i,method)
            !call loct_progress_bar(i-1, td%max_iter-1)
         end do
         ! write(stdout, '(1x)')
         
         ! write new field to file
         write(filename,'(a,i3.3)') 'opt-control/laser.', ctr_iter
         call write_field(filename, laser)
      
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_zr98


    !---------------------------------------
    subroutine scheme_wg05(method)
      character (len=*), intent(in) :: method
      
      integer :: ierr
      logical :: stoploop

      call push_sub('opt_control.scheme_WG05')
      
      stoploop = .FALSE.
    
      message(1) = "Info: Starting OCT iteration using scheme: WG05"
      call write_info(1)
      
      old_functional = -CNST(1e10)
      ctr_iter = 0

      ! SELECT case(controlmode)
      ! CONTROL MODE
      !   case("penalty")
      ! FILTER (FREQ / TIME)
      ! FIXED FLUENCE
      !      case DEFAULT
      !         ! penalty
      !      END SELECT

      ctr_loop: do
         
         ! first propagate chi to ti
         message(1) = "Info: Initial forward propagation"
         call write_info(1)
         
         psi = initial_st
         call prop_fwd(psi) 
         call zoutput_function(sys%outp%how,'opt-control','prop1',gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr)
         
         ! define target state
         call target_calc('ZR98',target_st,psi,chi) ! defines chi

         call iteration_manager(stoploop)
         if (stoploop)  exit ctr_loop

         message(1) = "Info: Setup backward"
         call write_info(1)

         ! setup backward propagation
         call states_calc_dens(chi, NP_PART, st%rho)
         call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
         do i = 1, st%d%nspin
            v_old_f(:, i, 2) = h%vhxc(:, i)
         end do
         v_old_f(:, :, 3) = v_old_f(:, :, 2)
         
         do i = 2, 3
            v_old_i(:,:,i) = v_old_i(:,:,1) ! this one comes the previous propagation
         end do
         ! and now backward
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
         ! write(stdout, '(1x)')
         write(filename,'(a,i3.3)') 'opt-control/b_laser.', ctr_iter
         call write_field(filename, laser_tmp)

         !! FILTER STUFF / ALPHA STUFF

         ! RECALC FIELD

         laser = laser_tmp
         !!
         
         ! write new field to file
         write(filename,'(a,i3.3)') 'opt-control/laser.', ctr_iter
         call write_field(filename, laser)
         
      end do ctr_loop

      call pop_sub()
    end subroutine scheme_wg05


    ! ---------------------------------------------------------
    subroutine scheme_zbr98(method)
      character (len=*), intent(in)   :: method

      integer :: ierr
      logical :: stoploop
      
      call push_sub('opt_control.scheme_zbr98')
      
      message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
      call write_info(1)
      
      stoploop = .FALSE.

      ! first propagate chi to T
      message(1) = "Info: Initial backward propagation"
      call write_info(1)
      
      call target_calc('ZBR98', target_st, psi, chi)
      call prop_bwd(chi) 
!! DEBUG
      call zoutput_function(sys%outp%how, 'opt-control', 'initial_propZBR98', gr%m, gr%sb, psi%zpsi(:,1,1,1), M_ONE, ierr)
!! stop

      old_functional = -CNST(1e10)
      ctr_iter = 0

      ctr_loop: do

         ! check for stop file and delete file

         message(1) = "Info: Setup forward"
         call write_info(1)

!! CHECK: WHY DO WE NEED st%rho ?
         call states_calc_dens(chi, NP_PART, st%rho)
         call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
         do i = 1, st%d%nspin
            v_old_f(:, i, 2) = h%vhxc(:, i)
         end do
         v_old_f(:, :, 3) = v_old_f(:, :, 2)
         
         do i = 2, 3
            v_old_i(:,:,i) = v_old_i(:,:,1) ! this one comes the previous propagation
         end do

         ! setup forward propagation
         psi = initial_st  
         call states_calc_dens(psi, NP_PART, st%rho)
         call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
         do i = 1, st%d%nspin
            v_old_i(:, i, 2) = h%vhxc(:, i)
         end do
         v_old_i(:, :, 3) = v_old_i(:, :, 2)
         
         do i = 2, 3
            v_old_f(:,:,i) = v_old_f(:,:,1) ! this one comes from the previous propagation
         end do
         message(1) = "Info: Propagating forward"
         call write_info(1)

         h%ep%lasers(1)%dt = td%dt
         do i = 1, td%max_iter
            call prop_iter_fwd(i,method)
         end do
         ! write(stdout, '(1x)')
         
         ! write new field to file
         write(filename,'(a,i3.3)') 'opt-control/laser.', ctr_iter
         call write_field(filename, laser)


         ! iteration managament ! make own subroutine
         call iteration_manager(stoploop)
         if (stoploop)  exit ctr_loop

         ! and now backward
         call target_calc('ZBR98',target_st,psi,chi)
         message(1) = "Info: Propagating backward"
         call write_info(1)
         td%dt = -td%dt
         h%ep%lasers(1)%dt = td%dt
         do i = td%max_iter-1, 0, -1
            call prop_iter_bwd(i,method)
         end do
         td%dt = -td%dt
         ! write(stdout, '(1x)')
         write(filename,'(a,i3.3)') 'opt-control/b_laser.', ctr_iter
         call write_field(filename, laser_tmp)

      end do ctr_loop

      call pop_sub()
    end subroutine scheme_zbr98


    ! ---------------------------------------------------------
    subroutine iteration_manager(stoploop)
      logical, intent(out) :: stoploop
      
      ! calculate overlap 
      message(1) = "Info: Calculate overlap"
      call write_info(1)
      call calc_J() ! psi mit target_st
      
      message(1) = "Info: Loop control"
      call write_info(1)
      
      ! TODO:: check for STOP FILE AND delete it

      if((ctr_iter .eq. ctr_iter_max).or.(eps>M_ZERO.and.abs(functional-old_functional) < eps)) then
         if((ctr_iter .eq. ctr_iter_max)) then
            message(1) = "Info: Maximum number of iterations reached"
            call write_info(1)
         endif
         if(eps > M_ZERO .and. abs(functional-old_functional) < eps ) then
            message(1) = "Info: Convergence threshold reached"
            call write_info(1)
         endif
 
         stoploop = .TRUE.

      end if
      write(message(1), '(a,i3)') 'Info: Optimal control iteration #', ctr_iter
      call write_info(1)
      ctr_iter = ctr_iter + 1
      old_functional = functional
      
    end subroutine iteration_manager

    
    ! ---------------------------------------------------------
    subroutine target_calc(method, targetst, psi_in, chi_out)
      ! calculate chi = \hat{O} psi
      ! do loop <target_st|Psi> for all States
      character (len=*), intent(in)  :: method
      type(states_t),    intent(in)  :: targetst, psi_in
      type(states_t),    intent(out) :: chi_out

      CMPLX   :: olap
      integer :: ik, p, dim
      
      do ik = 1, psi_i%d%nik
        do p  = psi_i%st_start, psi_i%st_end
           olap = M_z0
           do dim = 1, psi_i%d%dim
             olap = olap + zmf_integrate(gr%m,conjg(targetst%zpsi(:, dim, p, ik))*psi_in%zpsi(:, dim, p, ik))
          end do
          if(method == 'ZR98') &
               chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:, :, p, ik)
          if(method == 'ZBR98') &
               chi_out%zpsi(:,:,p,ik) = targetst%zpsi(:, :, p, ik)
          if(method == 'WG05') &
               chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:, :, p, ik)
       end do
    end do
    
!      do ik = 1, psi%d%nik
!        do p  = psi%st_start, psi%st_end
!           olap = M_z0     
!           olap = zstates_dotp(gr%m, psi_in%d%dim,  targetst%zpsi(:,:, p, ik), psi_in%zpsi(:,:, p, ik))
    !       if(method == 'ZR98') &
    !            chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:,:,p,ik)
    !       if(method == 'ZBR98') &
    !            chi_out%zpsi(:,:,p,ik) = targetst%zpsi(:,:,p,ik)
    !       if(method == 'WG05') &
    !            chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:,:,p,ik)
!! DEBUG
!           write(6,*) 'overlap: ', olap
!!
     !   end do
     !end do

   end subroutine target_calc


   ! ---------------------------------------------------------
   subroutine prop_iter_fwd(iter,method)
     character (len=*), intent(in) :: method
     integer,           intent(in) :: iter

     call push_sub('opt_control.prop_iter_fwd')
     
     ! chi(0) --> chi(T) [laser_tmp]
     !         |------------laser     
     ! psi(0) --> psi(T) 
     
      ! new electric field
      call update_field(iter-1, laser, method)!
    
      ! chi
      td%tr%v_old => v_old_f
      h%ep%lasers(1)%numerical => laser_tmp
      call states_calc_dens(chi, NP_PART, st%rho)
      call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, chi, td%tr, abs(iter*td%dt), abs(td%dt))

      ! psi
      td%tr%v_old => v_old_i
      h%ep%lasers(1)%numerical => laser
      call states_calc_dens(psi, NP_PART, st%rho)
      call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, abs(iter*td%dt), abs(td%dt))
      !write(6,*) iter, psi%zpsi(45:50,1,1,1)
      !write(6,*) laser(2*iter,1)

      call pop_sub()
    end subroutine prop_iter_fwd


    ! ---------------------------------------------------------
    subroutine prop_iter_bwd(iter, method)
      integer, intent(in)           :: iter
      character (len=*), intent(in) :: method

      call push_sub('opt_control.prop_iter_bwd')
      ! chi(T) --> chi(0)
      !         |------------laser_tmp
      ! psi(T) --> psi(0) [laser]
 
      ! new electric field
      call update_field(iter+1, laser_tmp, method)

      ! psi
      td%tr%v_old => v_old_i
      h%ep%lasers(1)%numerical => laser
      call states_calc_dens(psi_i, NP_PART, st%rho)
      call v_ks_calc(gr, sys%ks, h, psi_i, calc_eigenval=.true.)
      ! abs
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, abs(iter*td%dt), td%dt)

      ! chi
      td%tr%v_old => v_old_f
      h%ep%lasers(1)%numerical => laser_tmp
      call states_calc_dens(chi, NP_PART, st%rho)
      call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
      ! abs
      call td_rti_dt(sys%ks, h, gr, chi, td%tr, abs(iter*td%dt), td%dt)
!      write(6,*) laser_tmp(2*iter,1)

      call pop_sub()
    end subroutine prop_iter_bwd


    ! ---------------------------------------------------------
    subroutine update_field(iter, l, method)
      integer, intent(in)           :: iter
      FLOAT,   intent(inout)        :: l(0:2*td%max_iter, NDIM)
      character (len=*), intent(in) :: method
      
      CMPLX :: d1
      CMPLX :: d2(MAX_DIM)
      integer :: ik, p, dim, i, pol

      ! case SWITCH
      ! CHECK DIPOLE MOMENT FUNCTION - > VELOCITY GAUGE

!! REPLACE gr%m%x(:,pol) by Dipol: mu_x = x +y , mu_x = x+iy

      d1 = M_z0; d2 = M_z0 
      do ik = 1, psi_i%d%nik
        do p  = psi_i%st_start, psi_i%st_end
          do dim = 1, psi_i%d%dim
             do pol=1, NDIM
                !write(6,*) ik,p,dim,pol,d2
                d2(pol) = d2(pol) + zmf_integrate(gr%m,conjg(chi%zpsi(:, dim, p, ik))*gr%m%x(:,pol)*psi%zpsi(:, dim, p, ik))
             enddo
             d1 = d1 + zmf_integrate(gr%m,conjg(psi%zpsi(:, dim, p, ik))*chi%zpsi(:, dim, p, ik))
          end do
        end do
      end do

      if(method=="ZBR98") l(2*iter, 1:NDIM) = -aimag(d1*d2(1:NDIM))/penalty
      if(method=="ZR98")  l(2*iter, 1:NDIM) = -aimag(d2(1:NDIM))/penalty
      if(method=="WG05")  l(2*iter, 1:NDIM) = -aimag(d2(1:NDIM))/penalty

      !write(6,*) d2
      ! penalty als block einlesen !
      !l(2*iter, 1:NDIM) = -aimag(d1*d2(1:NDIM))/penalty
      l(2*iter, 1:NDIM) = -aimag(d2(1:NDIM)/penalty)
      !fluence =  sum(l(2*iter, 1:NDIM)**2)*abs(td%dt)

      ! extrapolate to t+-dt/2
      i = int(sign(M_ONE, td%dt))
      if(iter==0.or.iter==td%max_iter) then
        l(2*iter+  i, 1:NDIM) = l(2*iter, 1:NDIM)
        l(2*iter+2*i, 1:NDIM) = l(2*iter, 1:NDIM)
      else
        l(2*iter+  i, 1:NDIM) = M_HALF*(M_THREE*l(2*iter, 1:NDIM) -       l(2*iter-2*i, 1:NDIM))
        l(2*iter+2*i, 1:NDIM) = M_HALF*( M_FOUR*l(2*iter, 1:NDIM) - M_TWO*l(2*iter-2*i, 1:NDIM))
      end if
      !write(6,*) l(2*iter,1)
    end subroutine update_field


    ! ---------------------------------------------------------
    subroutine prop_bwd(psi_n) ! give initial chi and laser
      type(states_t), intent(inout)  :: psi_n

      message(1) = "Info: Backward propagating Chi"
      call write_info(1)

      ! setup the hamiltonian
      call states_calc_dens(psi_n, NP_PART, st%rho)
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
        call states_calc_dens(psi_n, NP_PART, st%rho)
        call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)

      end do
      td%dt = -td%dt
      message(1) = ""
      call write_info(1)
    end subroutine prop_bwd


    ! ---------------------------------------------------------
    subroutine prop_fwd(psi_n) ! give initial psi and laser
      type(states_t), intent(inout)  :: psi_n

      integer :: ierr

      message(1) = "Info: Forward propagating Psi"
      call write_info(1)

      ! setup the hamiltonian
      call states_calc_dens(psi_n, NP_PART, st%rho)
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
      
      do i = 1, td%max_iter 
        ! time iterate wavefunctions
        call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), abs(td%dt))

        ! update
        call states_calc_dens(psi_n, NP_PART, st%rho)
        call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
      
        !call loct_progress_bar(i-1, td%max_iter-1)
      end do

      call zoutput_function(sys%outp%how,'opt-control','target',gr%m,gr%sb,psi_n%zpsi(:,1,1,1),M_ONE, ierr) 

      message(1) = ""
      call write_info(1)
    end subroutine prop_fwd


    ! ---------------------------------------------------------
    subroutine calc_J()

      fluence = SUM(laser**2)*abs(td%dt)
    
      call calc_overlap()

      convergence(1,ctr_iter) = functional
      convergence(2,ctr_iter) = overlap
      convergence(3,ctr_iter) = fluence
      convergence(4,ctr_iter) = penalty

    end subroutine calc_J


    ! ---------------------------------------------------------
    subroutine calc_overlap()
      integer :: ik, p, dim

      overlap = M_z0;
      do ik = 1, psi%d%nik
        do p  = psi%st_start, psi%st_end
           do dim = 1, psi_i%d%dim
          ! WARNING gives garbage when calculated through zstates_dotp
          !overlap = zstates_dotp(gr%m, psi%d%dim, psi%zpsi(:,:, p, ik), target_st%zpsi(:,:, p, ik))
              !write(6,*) p, psi%zpsi(1:3, dim, p, ik)
              !write(6,*) ik,target_st%zpsi(1:3, dim, p, ik)
              overlap= overlap + abs(zmf_integrate(gr%m,conjg(psi%zpsi(:, dim, p, ik))*target_st%zpsi(:, dim, p, ik)))
              functional = overlap - penalty*fluence
              write(message(1), '(6x,i3,1x,i3,a,f14.8,a,f14.8,a,f14.8)') &
                ik, p, " => J1:", overlap, "   J: " , functional,  "  I: " , fluence
              call write_info(1)
           end do
        end do
     end do
   
   end subroutine calc_overlap


    ! ---------------------------------------------------------
    subroutine read_state(st, m, filename)
      type(states_t),   intent(out) :: st
      type(mesh_t),     intent(in)  :: m
      character(len=*), intent(in)  :: filename

      FLOAT   :: phi(1:m%np_part)
      integer :: ierr

      !write(6,*) 'mesh_change: ', mesh_change
      !call zrestart_read_function('tmp/restart_gs',filename, gr%m,phi, ierr)
      call dinput_function('tmp/restart_gs/'//trim(filename), m, phi(:), ierr, is_tmp=.TRUE.)
      !call zrestart_read('tmp/restart_gs', st,gr, ierr)

      st%zpsi(:,1,1,1) = phi
      if(ierr.ne.0) then
         message(1) = "Unsuccesfull read of states in 'tmp/restart_gs/" // trim(filename) // "'"
        call write_fatal(1)
      end if
    end subroutine read_state

    
    ! ---------------------------------------------------------
    subroutine write_field(filename, las)
      character(len=*), intent(in) :: filename
      FLOAT,            intent(in) :: las(1:NDIM,0:2*td%max_iter)

      integer :: i, iunit

      iunit = io_open(filename, action='write')
      do i = 0, 2*td%max_iter
        write(iunit, '(4es20.12)') i*td%dt/M_TWO, las(:, i)
      end do
      call io_close(iunit)

    end subroutine write_field
    

    ! ---------------------------------------------------------
    subroutine output()
      integer :: iunit, loop, ierr

      call push_sub('opt_control.output')
      
      iunit = io_open('opt-control/info', action='write')
      write(iunit, '(a,i4)')    'Iterations = ', ctr_iter
      write(iunit, '(a,f14.8)') 'Overlap    = ', overlap
      write(iunit, '(a,f14.8)') 'Functional = ', functional
      call io_close(iunit)
      message(1) = "Info: Output States"
      call write_info(1)
      ! should output wavefunctions ;)
      
      ! dump convergence: J,P,fluence,penalty
      iunit = io_open('opt-control/convergence', action='write')
      ! HEADER
      write(iunit, '(4(a))') '# iteration ','functional ','overlap ','penalty '
      ! DATA
      do loop=1,ctr_iter_max
         write(iunit, '(i6,4(f14.8))') loop, convergence(1,loop), convergence(2,loop), convergence(3,loop),convergence(4,loop)
      end do
      call io_close(iunit)
      do loop=1,ctr_iter_max
         write(6,'(f14.8,f14.8)') convergence(1,loop),convergence(2,loop)
      end do
      ! assign final wave function to psi_i
      psi_i%zpsi = psi%zpsi


      ! CHECK: STH GOES WRONG HERE: IF psi_i REPLACED BY PSI 
      ! probably something is ill-defined
      call states_output(psi_i, gr, 'opt-control', sys%outp)
      call zoutput_function(sys%outp%how,'opt-control','target',gr%m,gr%sb,target_st%zpsi(:,1,1,1),M_ONE,ierr) 
      call zoutput_function(sys%outp%how,'opt-control','final',gr%m,gr%sb,psi%zpsi(:,1,1,1),M_ONE,ierr) 
      call pop_sub()
    end subroutine output

    !subroutine dipole_moment_function()
    ! define dipole moment function 
    ! free in x- and y
    ! mu_x = x 
    ! mu_y = y
    ! or force elliptic / linear : use only first component
    ! mu_x = x + iy
    !

    !end subroutine dipole_moment_function


    subroutine def_istate()

      ! prepare initialstate
      !%Variable OCTInitialState
      !%Type String
      !%Section Optimal Control
      !%Description
      !% The string OCTInitialState describes the initial state of the quantum system
      !% Possible arguments are:
      !%Option istate_gs 1
      !% start in the ground state 
      !%Option istate_ex 2
      !% start in the excited state given by OCTISnumber
      !% (ordered by energy)
      !%Option istate_sp 3
      !% start in a superposition of states defined by the block OCTISsuperposition)
      !%Option istate_ud 4
      !% start in a userdefined state 
      !%End

      ! parse input
      call loct_parse_int(check_inp('OCTInitialState'), oct_istate_gs, istype)
     

      ! make it work for single particle
      call read_state(initial_st, gr%m, "0000000001")
      ! check options in parameter file

      ! initial superposition
      ! call read_superposition

      ! userdefined state
      ! call user defined state
      
    end subroutine def_istate
    

    ! ----------------------------------------------------------------------
    subroutine def_toperator()

      ! prepare targetoperator      
      !%Variable OCTTargetState
      !%Type String
      !%Section Optimal Control
      !%Description
      !% The string OCTTargetOperator describes the initial state of the quantum system
      !% Possible arguments are:
      !% state_gs - Targetoperator is a projection operator on the ground state 
      !% state_ex - Targetoperator is a projection operator on the excited state given by OCTTOnumber
      !% (ordered by energy)
      !% state_sp - Targetoperator is a projection operator on a superposition of states defined by the block OCTTOsuperposition)
      !% state_ud - Targetoperator is a projection operator on a user defined state
      !% local_ud - Targetoperator is a local operator
      !%End

      ! parse input
      call loct_parse_string(check_inp('OCTTargetOperator'),'gs', totype)
      call read_state(target_st, gr%m, "0000000002")

    end subroutine def_toperator


    ! ---------------------------------------------------------
    subroutine init_()
      integer            :: ierr, kk, jj

      call push_sub('opt_control.init_')  

      call io_mkdir('opt-control')
      ! some shortcuts
      gr  => sys%gr
      st  => sys%st

      ! prepare unit factor u for output
      u  = M_ONE/units_out%length%factor**NDIM

      call td_init(gr, td, sys%st, sys%outp)

      call states_allocate_wfns(st, gr%m, M_CMPLX)

      ! psi_i is initialized in system_init
      psi_i => st

      ! call write_info(2)
      v_old_i => td%tr%v_old

      ! now we initialize chi. This will repeat some stuff
      call states_init(psi, gr)
      call states_init(chi, gr)
      call states_init(initial_st, gr)
      call states_init(target_st, gr)
      if(h%ep%nvnl > 0) then
        ALLOCATE(psi%rho_core(NP_PART), NP_PART)
        psi%rho_core(NP_PART) = psi_i%rho_core(NP_PART)  
        
        ALLOCATE(initial_st%rho_core(NP_PART), NP_PART)
        initial_st%rho_core(NP_PART) = psi_i%rho_core(NP_PART)  
        
        ALLOCATE(chi%rho_core(NP_PART), NP_PART)
        chi%rho_core(NP_PART) = psi_i%rho_core(NP_PART) 
        
        ALLOCATE(target_st%rho_core(NP_PART), NP_PART)
        target_st%rho_core(NP_PART) = psi_i%rho_core(NP_PART)
      end if 

      psi%st_start = psi_i%st_start
      psi%st_end = psi_i%st_end
      initial_st%st_start = psi_i%st_start
      initial_st%st_end = psi_i%st_end
      target_st%st_start = psi_i%st_start
      target_st%st_end = psi_i%st_end
      chi%st_start = psi_i%st_start
      chi%st_end = psi_i%st_end
      call states_allocate_wfns(psi,        gr%m, M_CMPLX)
      call states_allocate_wfns(chi,        gr%m, M_CMPLX)
      call states_allocate_wfns(initial_st, gr%m, M_CMPLX)
      call states_allocate_wfns(target_st,  gr%m, M_CMPLX)

      ALLOCATE(v_old_f(NP_PART, chi%d%nspin, 3), NP_PART*chi%d%nspin*3)


      ! INITIAL LASER FIELD 

      ! allocate memory
      ALLOCATE(laser(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter+1))
      ALLOCATE(laser_tmp(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter+1))
      ALLOCATE(field(1:NDIM), NDIM)

      ! use same laser as time-dependent run mode
      call laser_init(h%ep%no_lasers,h%ep%lasers,gr%m)
      !! - built one numerical laser
      laser_tmp = M_ZERO  
      laser     = M_ZERO

      do jj = 1, td%max_iter
         t = td%dt*(jj-1)/M_TWO
         !i = int(abs(M_TWO*t/l(1)%dt) + M_HALF)
         ! 
         do kk=1, h%ep%no_lasers
            call laser_field(gr%sb,h%ep%no_lasers,h%ep%lasers,t,field)
            laser(:,jj) = laser(:,jj) + field
         end do
      enddo

      write(filename,'(a)') 'opt-control/initial_laser.1'

      call write_field(filename, laser)
      !! - deallocate multiple lasers
      call laser_end(h%ep%no_lasers,h%ep%lasers)
      !! - NOW: allocate just one laser for optimal control
      h%ep%no_lasers = 1
      ALLOCATE(h%ep%lasers(1), 1)
      h%ep%lasers(1)%envelope = 99 ! internal type
      h%ep%lasers(1)%dt = td%dt
      h%ep%lasers(1)%numerical => laser   
      h%ep%lasers(1)%numerical => laser_tmp

      write(filename,'(a)') 'opt-control/initial_laser.2'
      call write_field(filename, laser)
      ! write(6,*) 'LASER STRENGTH', SUM(laser**2)*abs(td%dt)

      ! initial laser definition end
      
   
      !! read in oct parameters
      ! description here
      call loct_parse_float(check_inp('OCTPenalty'), M_ONE, penalty)
      call loct_parse_float(check_inp('OCTEps'), CNST(1e-3), eps)
      call loct_parse_int(check_inp('OCTMaxIter'), 10, ctr_iter_max)
      call loct_parse_int(check_inp('OCTTargetMode'),0,targetmode)
      call loct_parse_int(check_inp('OCControlFilterMode'),0,filtermode)
      call loct_parse_int(check_inp('OCTScheme'), oct_algorithm_zr98, algorithm_type)


      if(ctr_iter_max < 0.and.eps<M_ZERO) then
        message(1) = "OptControlMaxIter and OptControlEps can not be both <0"
        call write_fatal(1)
      end if

      if(ctr_iter_max < 0) ctr_iter_max = huge(ctr_iter_max)

      ALLOCATE(convergence(4,ctr_iter_max),ctr_iter_max*4)
      ! initial state
      ! replace by extra routine
      call def_istate()
      call def_toperator()

      ! write(6,*) size(psi_i%zpsi)
      ! write(6,*) shape(psi_i%zpsi)

      call zoutput_function(sys%outp%how,'opt-control','initial1',gr%m,gr%sb,initial_st%zpsi(:,1,1,1),u,ierr)
      call zoutput_function(sys%outp%how,'opt-control','target1',gr%m,gr%sb,target_st%zpsi(:,1,2,1),u,ierr) 

     call pop_sub()
    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      call td_end(td)
    end subroutine end_

  end subroutine opt_control_run

end module opt_control_m
