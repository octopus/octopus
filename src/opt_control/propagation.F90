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
!! $Id: opt_control.F90 2873 2007-04-29 22:05:29Z acastro $

#include "global.h"

module opt_control_propagation_m
  use datasets_m
  use varinfo_m
  use global_m
  use lib_oct_m
  use messages_m
  use units_m
  use grid_m
  use mesh_function_m
  use states_m
  use excited_states_m
  use hamiltonian_m
  use system_m
  use timedep_m
  use td_rti_m
  use td_write_m
  use opt_control_constants_m
  use opt_control_tdtarget_m
  use opt_control_parameters_m
  use v_ks_m
  use tdf_m

  private
  public :: propagate_forward,  &
            propagate_backward, &
            fwd_step,           &
            bwd_step,           &
            calc_tdfitness,     &
            update_field


  contains

  ! ---------------------------------------------------------
  subroutine propagate_forward(targetmode, sys, h, td, tdt, psi_n, write_iter)
    integer, intent(in) :: targetmode
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(td_t),          intent(inout) :: td
    type(td_target_set_t), intent(inout)  :: tdt
    type(states_t),      intent(inout) :: psi_n
    logical, optional, intent(in)      :: write_iter

    integer :: ierr, ii, i
    logical :: write_iter_ = .false.
    FLOAT   :: etime
    FLOAT, allocatable :: dens(:,:)
    type(grid_t),  pointer :: gr
    type(td_write_t)           :: write_handler

    call push_sub('opt_control.propagate_forward')

    message(1) = "Info: Forward propagating Psi"
    call write_info(1)

    write_iter_ = .false.
    if(present(write_iter)) write_iter_ = write_iter

    gr => sys%gr

    if(write_iter_) then
      call td_write_init(write_handler, gr, sys%st, sys%geo, (td%move_ions>0), h%ep%with_gauge_field, td%iter, td%dt)
    end if

    ALLOCATE(dens(NP_PART, psi_n%d%nspin), NP_PART*psi_n%d%nspin) 

    ! setup the hamiltonian
    call states_calc_dens(psi_n, NP_PART, dens)
    psi_n%rho = dens
    call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
    call td_rti_run_zero_iter(h, td%tr)

    ii = 1

    etime = loct_clock()
    do i = 1, td%max_iter 
      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), abs(td%dt), td%max_iter)

      ! if td_target
      if(targetmode==oct_targetmode_td) &
        call calc_tdfitness(tdt, gr, psi_n, tdt%td_fitness(i))
      
      ! update
      call states_calc_dens(psi_n, NP_PART, dens)
      psi_n%rho = dens
      call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
      call hamiltonian_energy(h, sys%gr, sys%geo, psi_n, -1)

      ! only write in final run
      if(write_iter_) then

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

    end do
    
    if(write_iter_) call td_write_end(write_handler)
    deallocate(dens)
    call pop_sub()
  end subroutine propagate_forward


  ! ---------------------------------------------------------
  subroutine propagate_backward(sys, h, td, psi_n)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(td_t),          intent(inout) :: td
    type(states_t), intent(inout)      :: psi_n

    integer :: i
    FLOAT, allocatable :: dens(:,:)
    type(grid_t),  pointer :: gr

    call push_sub('opt_control.propagate_backward')
    
    message(1) = "Info: Backward propagating Chi"
    call write_info(1)

    gr => sys%gr

    ALLOCATE(dens(NP_PART, psi_n%d%nspin), NP_PART*psi_n%d%nspin) 

    ! setup the hamiltonian
    call states_calc_dens(psi_n, NP_PART, dens)
    psi_n%rho = dens
    call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
    call td_rti_run_zero_iter(h, td%tr)

    td%dt = -td%dt
    do i = td%max_iter-1, 0, -1
      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), td%dt, td%max_iter)
      ! update
      call states_calc_dens(psi_n, NP_PART, dens)
      psi_n%rho = dens
      call v_ks_calc(gr, sys%ks, h, psi_n, calc_eigenval=.true.)
      
    end do
    td%dt = -td%dt
    
    deallocate(dens)
    call pop_sub()
  end subroutine propagate_backward



  ! ---------------------------------------------------------
  ! chi(0) --> chi(T) [laser_tmp]
  !         |------------laser     
  ! psi(0) --> psi(T) 
  ! ---------------------------------------------------------
  subroutine prop_iter_fwd(iter, method, targetmode, gr, ks, h, tdt, td, par, par_tmp, psi, chi, psi2)
    integer,           intent(in) :: iter
    integer, intent(in) :: method
    integer, intent(in) :: targetmode
    type(grid_t), intent(inout) :: gr
    type(v_ks_t), intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: h
    type(td_target_set_t), intent(inout)  :: tdt
    type(td_t), intent(inout) :: td
    type(oct_control_parameters_t), intent(inout) :: par, par_tmp
    type(states_t), intent(inout) :: psi, chi, psi2


    FLOAT, allocatable :: dens_tmp(:,:)
    integer :: nspin
    
    call push_sub('opt_control.prop_iter_fwd')

    nspin = psi%d%nspin

    ALLOCATE(dens_tmp(NP_PART, nspin), NP_PART*nspin) 

    ! new electric field
    if(targetmode==oct_targetmode_td) then
      ! psi(0) --> psi(T) with old field (for chi)
      call calc_inh(psi2, gr, tdt, iter, td%max_iter, td%dt, chi)
      ! psi2
      call parameters_to_h(par, h%ep)
      call states_calc_dens(psi2, NP_PART, dens_tmp)
      psi2%rho = dens_tmp
      call v_ks_calc(gr, ks, h, psi2, calc_eigenval=.true.)
      call td_rti_dt(ks, h, gr, psi2, td%tr, abs(iter*td%dt), abs(td%dt), td%max_iter)
    end if

    call update_field(method, iter-1, par, gr, td, psi, chi)
    
    ! chi
    call parameters_to_h(par_tmp, h%ep)
    call states_calc_dens(chi, NP_PART, dens_tmp)
    chi%rho = dens_tmp
    call v_ks_calc(gr, ks, h, chi, calc_eigenval=.true.)
    call td_rti_dt(ks, h, gr, chi, td%tr, abs(iter*td%dt), abs(td%dt), td%max_iter)
    
    ! psi
    call parameters_to_h(par, h%ep)
    call states_calc_dens(psi, NP_PART, dens_tmp)
    psi%rho = dens_tmp
    call v_ks_calc(gr, ks, h, psi, calc_eigenval=.true.)
    call td_rti_dt(ks, h, gr, psi, td%tr, abs(iter*td%dt), abs(td%dt), td%max_iter)

    deallocate(dens_tmp)     
    call pop_sub()
  end subroutine prop_iter_fwd


  ! ----------------------------------------------------------
  subroutine fwd_step(method, targetmode, sys, td, h, tdt, par, par_tmp, chi, psi)
    integer, intent(in)                           :: method
    integer, intent(in)                           :: targetmode
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(td_target_set_t), intent(inout)          :: tdt
    type(oct_control_parameters_t), intent(inout) :: par, par_tmp
    type(states_t), intent(inout)                 :: chi
    type(states_t), intent(inout)                 :: psi

    integer :: i, nspin
    FLOAT, allocatable :: dens_tmp(:,:)
    type(states_t) :: psi2

    type(grid_t), pointer :: gr

    call push_sub('opt_control.fwd_step')

    gr => sys%gr
    nspin = psi%d%nspin

    ALLOCATE(dens_tmp(NP_PART, nspin), NP_PART*nspin) 

    psi2 = psi
    
    ! setup forward propagation
    call states_densities_init(psi, gr, sys%geo)
    call states_calc_dens(psi, NP_PART, dens_tmp)
    psi%rho = dens_tmp
    call v_ks_calc(gr, sys%ks, h, psi, calc_eigenval=.true.)
    call td_rti_run_zero_iter(h, td%tr)

    message(1) = "Info: Propagating forward"
    call write_info(1)

    do i = 1, td%max_iter
      call prop_iter_fwd(i, method, targetmode, gr, sys%ks, h, tdt, td, par, par_tmp, psi, chi, psi2)
      if(targetmode==oct_targetmode_td) call calc_tdfitness(tdt, gr, psi, tdt%td_fitness(i))     
    end do

    call states_end(psi2)
    nullify(gr)
    deallocate(dens_tmp)      
    call pop_sub()
  end subroutine fwd_step


  ! --------------------------------------------------------
  subroutine bwd_step(method, targetmode, sys, td, h, tdt, par, par_tmp, chi, psi) 
    integer, intent(in) :: method
    integer, intent(in) :: targetmode
    type(system_t), intent(inout) :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(td_target_set_t), intent(inout)          :: tdt
    type(oct_control_parameters_t), intent(inout) :: par, par_tmp
    type(states_t), intent(inout)                 :: chi
    type(states_t), intent(inout)                 :: psi

    integer :: i, nspin
    FLOAT, allocatable :: dens_tmp(:,:)
    type(grid_t), pointer :: gr

    call push_sub('opt_control.bwd_step')

    gr => sys%gr
    nspin = psi%d%nspin

    ALLOCATE(dens_tmp(NP_PART, nspin), NP_PART*nspin) 

    ! setup backward propagation
    call states_calc_dens(chi, NP_PART, dens_tmp)
    chi%rho = dens_tmp
    call v_ks_calc(gr, sys%ks, h, chi, calc_eigenval=.true.)
    call td_rti_run_zero_iter(h, td%tr)

    message(1) = "Info: Propagating backward"
    call write_info(1)

    td%dt = -td%dt
    do i = td%max_iter-1, 0, -1
      call prop_iter_bwd(i, method, targetmode, gr, sys%ks, h, tdt, td, par, par_tmp, psi, chi)
    end do
    td%dt = -td%dt

    deallocate(dens_tmp)
    call pop_sub()
  end subroutine bwd_step


  ! ---------------------------------------------------------
  ! do backward propagation step with update of field
  ! works for time-independent and td targets
  ! chi(T) --> chi(0)
  !         |------------laser_tmp
  ! psi(T) --> psi(0) [laser]
  ! ---------------------------------------------------------
  subroutine prop_iter_bwd(iter, method, targetmode, gr, ks, h, tdt, td, par, par_tmp, psi, chi)
    integer, intent(in) :: iter
    integer, intent(in) :: method
    integer, intent(in) :: targetmode
    type(grid_t), intent(inout) :: gr
    type(v_ks_t), intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: h
    type(td_target_set_t), intent(inout)  :: tdt
    type(td_t), intent(inout) :: td
    type(oct_control_parameters_t), intent(inout) :: par, par_tmp
    type(states_t), intent(inout) :: psi, chi
    FLOAT, allocatable :: dens_tmp(:,:)

    integer :: nspin
    
    call push_sub('opt_control.prop_iter_bwd')

    nspin = psi%d%nspin

    ALLOCATE(dens_tmp(NP_PART, nspin), NP_PART*nspin) 
     
    ! calc inh first, then update field
    if(targetmode==oct_targetmode_td) then
      call calc_inh(psi, gr, tdt, iter, td%max_iter, td%dt, chi)
    end if

    ! new electric field
    call update_field(method, iter+1, par_tmp, gr, td, psi, chi)
         
    ! chi
    call parameters_to_h(par_tmp, h%ep)
    call states_calc_dens(chi, NP_PART, dens_tmp)
    chi%rho = dens_tmp
    call v_ks_calc(gr, ks, h, chi, calc_eigenval=.true.)
    call td_rti_dt(ks, h, gr, chi, td%tr, abs(iter*td%dt), td%dt, td%max_iter)
    
    ! psi
    call  parameters_to_h(par, h%ep)
    call states_calc_dens(psi, NP_PART, dens_tmp)
    psi%rho = dens_tmp
    call v_ks_calc(gr, ks, h, psi, calc_eigenval=.true.)
    call td_rti_dt(ks, h, gr, psi, td%tr, abs(iter*td%dt), td%dt, td%max_iter)

    deallocate(dens_tmp)     
    call pop_sub()
  end subroutine prop_iter_bwd


  ! ---------------------------------------------------------
  subroutine calc_tdfitness(tdt, gr, psi, merit)
    type(td_target_set_t), intent(inout) :: tdt
    type(grid_t),      intent(in) :: gr
    type(states_t),    intent(in) :: psi
    FLOAT,             intent(out):: merit
    integer             :: jj, ik, p, dim

    call push_sub('opt_control.calc_tdfitness')

    merit = M_ZERO
    do jj = 1, tdt%no_tdtargets
      if(tdt%tdtg(jj)%type.eq.oct_tgtype_local) then
        do ik = 1, psi%d%nik
          do p  = psi%st_start, psi%st_end
            do dim = 1, psi%d%dim
              merit = merit + &
                zmf_integrate(gr%m, tdt%tdtarget(:)* &
                abs(psi%zpsi(:,dim,ik,p))**2)
            end do
          end do
        end do
      else
        do ik = 1, psi%d%nik
          do p  = psi%st_start, psi%st_end
            do dim = 1, psi%d%dim
              merit = merit + &
                abs(zmf_integrate(gr%m, tdt%tdtarget(:)* &
                conjg(psi%zpsi(:,dim,ik,p))))**2
            end do
          end do
        end do
      end if
    end do

    call pop_sub()
  end subroutine calc_tdfitness


  ! ---------------------------------------------------------
  subroutine update_field(algorithm_type, iter, cp, gr, td, psi, chi)
    integer, intent(in) :: algorithm_type
    integer, intent(in)        :: iter
    type(oct_control_parameters_t), intent(inout) :: cp
    type(grid_t), intent(in)   :: gr
    type(td_t), intent(in)     :: td
    type(states_t), intent(in) :: psi
    type(states_t), intent(in) :: chi
    
    CMPLX :: d1
    CMPLX :: d2(NDIM)
    integer :: ik, p, dim, i, pol
    FLOAT :: value

    CMPLX, allocatable :: rpsi(:, :)
    
    call push_sub('opt_control.update_field')
    
    ALLOCATE(rpsi(gr%m%np_part, psi%d%dim), gr%m%np_part*psi%d%dim)

    ! TODO This should be a product between Slater determinants.
    d2 = M_z0 
    do ik = 1, psi%d%nik
      do p  = psi%st_start, psi%st_end
        do pol = 1, NDIM
          do dim = 1, psi%d%dim
            rpsi(:, dim) = psi%zpsi(:, dim, p, ik)*cp%laser_pol(pol, 1)*gr%m%x(:, pol)
          end do
          d2(pol) = zstates_dotp(gr%m, psi%d%dim, chi%zpsi(:, :, p, ik), rpsi)
        end do
      end do
    end do
    deallocate(rpsi)

    d1 = M_z1
    if(algorithm_type .eq. oct_algorithm_zbr98) d1 = zstates_mpdotp(gr%m, psi, chi)

    value = aimag(d1*d2(1))/tdf(cp%td_penalty(1), iter+1)
    call tdf_set_numerical(cp%f(1), iter+1, value)
    i = int(sign(M_ONE, td%dt))
    if(iter==0.or.iter==td%max_iter) then
      call tdf_set_numerical(cp%f(1), iter+1, value)
    else
      value = M_HALF*( M_FOUR*tdf(cp%f(1), iter+1) - M_TWO*tdf(cp%f(1), iter-i+1))
      call tdf_set_numerical(cp%f(1), iter+i+1, value)
    end if

    call pop_sub()
  end subroutine update_field

end module opt_control_propagation_m
