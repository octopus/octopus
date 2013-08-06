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
!! $Id: target.F90 2870 2007-04-28 06:26:47Z acastro $

#include "global.h"

module target_m
  use datasets_m
  use density_m
  use derivatives_m
  use epot_m
  use excited_states_m
  use fft_m
  use forces_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use ion_dynamics_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use opt_control_global_m
  use opt_control_state_m
  use output_m
  use parser_m
  use profiling_m
  use restart_m
  use species_m
  use species_pot_m
  use spectrum_m
  use states_m
  use states_calc_m
  use states_dim_m
  use string_m
  use td_calc_m
  use td_m
  use types_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public :: target_t,                &
            target_get_state,        &
            target_init,             &
            target_end,              &
            target_output,           &
            target_tdcalc,           &
            target_inh,              &
            target_mode,             &
            target_type,             &
            target_j1,               &
            target_chi,              &
            target_move_ions,        &
            target_curr_functional,  &
            target_init_propagation


  integer, public, parameter ::       &
    oct_tg_groundstate      = 1,      &
    oct_tg_excited          = 2,      &
    oct_tg_gstransformation = 3,      &
    oct_tg_userdefined      = 4,      &
    oct_tg_jdensity         = 5,      &        
    oct_tg_local            = 6,      &
    oct_tg_td_local         = 7,      &
    oct_tg_exclude_state    = 8,      &
    oct_tg_hhg              = 9,      &
    oct_tg_velocity         = 10,     &
    oct_tg_hhgnew           = 12,     &
    oct_tg_classical        = 13

  integer, public, parameter ::       &
    oct_targetmode_static = 0,        &
    oct_targetmode_td     = 1

  integer, public, parameter ::       &
    oct_no_curr              = 0,     &
    oct_curr_square          = 1,     &
    oct_max_curr_ring        = 2,     &
    oct_curr_square_td       = 3

  type target_t
    private
    integer :: type
    type(states_t) :: st
    type(excited_states_t) :: est
    FLOAT, pointer :: rho(:) => null()
    FLOAT, pointer :: td_fitness(:) => null()
    character(len=200) :: td_local_target
    character(len=80) :: excluded_states_list
    character(len=4096) :: vel_input_string
    character(len=4096) :: classical_input_string
    character(len=1024), pointer :: vel_der_array(:,:) => null()
    character(len=1024), pointer :: mom_der_array(:,:) => null()
    character(len=1024), pointer :: pos_der_array(:,:) => null()
    FLOAT, pointer :: grad_local_pot(:,:,:) => null()
    logical :: move_ions
    integer :: hhg_nks
    integer, pointer :: hhg_k(:) => null()
    FLOAT,   pointer :: hhg_alpha(:) => null()
    FLOAT,   pointer :: hhg_a(:) => null()
    FLOAT   :: hhg_w0
    FLOAT   :: dt
    integer :: curr_functional
    FLOAT   :: density_weight
    FLOAT   :: curr_weight
    integer :: strt_iter_curr_tg
    FLOAT, pointer :: spatial_curr_wgt(:) => null()
    character(len=1000) :: plateau_string
    CMPLX, pointer :: acc(:, :)
    CMPLX, pointer :: vel(:, :)
    CMPLX, pointer :: gvec(:, :)
    FLOAT, pointer :: alpha(:)
    type(fft_t) :: fft_handler
  end type target_t


contains

  ! ---------------------------------------------------------
  !> This routine performs all the things that must be initialized
  !! prior to a forward evolution, regarding the target. Right now
  !! some of those initializations are not done here, and should
  !! be moved.
  subroutine target_init_propagation(tg)
    type(target_t), intent(inout)    :: tg
    PUSH_SUB(target_init_propagation)

    select case(tg%type)
    case(oct_tg_hhgnew)
      tg%vel = M_z0
      tg%gvec = M_z0
      tg%acc = M_z0
    end select

    POP_SUB(target_init_propagation)
  end subroutine target_init_propagation


  ! ----------------------------------------------------------------------
  !> This just copies the states_t variable present in target, into st.
  subroutine target_get_state(tg, st)
    type(target_t), intent(in)    :: tg
    type(states_t), intent(inout) :: st

    PUSH_SUB(target_get_state)
    call states_copy(st, tg%st)

    POP_SUB(target_get_state)
  end subroutine target_get_state
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> The target is initialized, mainly by reading from the inp file.
  subroutine target_init(gr, geo, qcs, td, w0, tg, oct, ep)
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(in)    :: geo
    type(opt_control_state_t),   intent(inout) :: qcs
    type(td_t),       intent(in)    :: td
    FLOAT,            intent(in)    :: w0
    type(target_t),   intent(inout) :: tg
    type(oct_t),      intent(in)    :: oct
    type(epot_t),     intent(inout) :: ep

    type(states_t), pointer :: stin

    PUSH_SUB(target_init)

    stin => opt_control_point_qs(qcs)

    !%Variable OCTTargetOperator
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default oct_tg_gstransformation
    !%Description
    !% The variable <tt>OCTTargetOperator</tt> prescribes which kind of target functional is
    !% to be used.
    !%Option oct_tg_groundstate 1 
    !% The target operator is a projection operator on the ground state, <i>i.e.</i> the
    !% objective is to populate the ground state as much as possible.
    !%Option oct_tg_excited 2
    !% The target operator is an "excited state". This means that the target operator
    !% is a linear combination of Slater determinants, each one formed by replacing
    !% in the ground-state Slater determinant one occupied state with one excited
    !% state (<i>i.e.</i> "single excitations"). The description of which excitations are
    !% used, and with which weights, should be given in a file called
    !% <tt>oct-excited-state-target</tt>. This is still in very preliminary, experimental
    !% phase. See the documentation of subroutine <tt>excited_states_init</tt> in the source
    !% code in order to use this feature.
    !%Option oct_tg_gstransformation 3
    !% The target operator is a projection operator on a transformation of the ground-state 
    !% orbitals defined by the block <tt>OCTTargetTransformStates</tt>.
    !%Option oct_tg_userdefined 4
    !% Allows to define target state by using <tt>OCTTargetUserdefined</tt>.
    !%Option oct_tg_jdensity 5
    !% EXPERIMENTAL: 
    !%Option oct_tg_local 6
    !% The target operator is a local operator.
    !%Option oct_tg_td_local 7
    !% The target operator is a time-dependent local operator.
    !%Option oct_tg_exclude_state 8
    !% Target operator is the projection onto the complement of a given state, given by the
    !% block <tt>OCTTargetTransformStates</tt>. This means that the target operator is the unity
    !% operator minus the projector onto that state.
    !%Option oct_tg_hhg 9
    !% The target is the optimization of the HHG yield. You must supply the OCTOptimizeHarmonicSpectrum
    !% block, and it attempts to optimize te maximum of the spectrum around each harmonic peak. You may
    !% use only one of the gradient-less optimization schemes.
    !%Option oct_tg_velocity 10
    !% The target is a function of the velocities of the nuclei at the end of the influence of
    !% the external field, defined by <tt>OCTVelocityTarget</tt>
    !%Option oct_tg_hhgnew 12
    !% EXPERIMENTAL: The  target is the optimization of the HHG yield. You must supply the
    !% OCTHarmonicWeigth string. It attempts to optimized the integral of the harmonic spectrum multiplied
    !% by some user defined weight function.
    !%Option oct_tg_classical 13
    !% EXPERIMENTAL
    !%End
    call parse_integer(datasets_check('OCTTargetOperator'), oct_tg_gstransformation, tg%type)
    if(.not.varinfo_valid_option('OCTTargetOperator', tg%type)) &
      call input_error('OCTTargetOperator')

    call states_copy(tg%st, stin)
    call states_deallocate_wfns(tg%st)
    call states_allocate_wfns(tg%st, gr%mesh, TYPE_CMPLX)
    nullify(tg%td_fitness)

    select case(tg%type)
    case(oct_tg_groundstate)
      call target_init_groundstate(gr, tg)
    case(oct_tg_excited)
      call messages_experimental('OCTTargetOperator = oct_tg_excited')
      call target_init_excited(gr, tg)
    case(oct_tg_exclude_state)
      call target_init_exclude(gr, tg)
    case(oct_tg_gstransformation)
      call target_init_gstransformation(gr, tg)
    case(oct_tg_userdefined) 
      call target_init_userdefined(gr, tg)
    case(oct_tg_jdensity)
      call target_init_density(gr, tg, stin, td)
    case(oct_tg_local)
      call target_init_local(gr, tg)
    case(oct_tg_td_local)
      call target_init_tdlocal(gr, tg, td)
    case(oct_tg_hhg)
      call target_init_hhg(tg, td, w0)
    case(oct_tg_hhgnew)
      call messages_experimental('OCTTargetOperator = oct_tg_hhgnew')
      call target_init_hhgnew(gr, tg, td, geo, ep)
    case(oct_tg_velocity)
      call target_init_velocity(gr, geo, tg, oct, td, ep)
    case(oct_tg_classical)
      call messages_experimental('OCTTargetOperator = oct_tg_classical')
      call target_init_classical(geo, tg, td)
    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call messages_fatal(1)
    end select

    nullify(stin)
    POP_SUB(target_init)
  end subroutine target_init
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_end(tg, oct)
    type(target_t), intent(inout) :: tg
    type(oct_t), intent(in)       :: oct

    PUSH_SUB(target_end)

    call states_end(tg%st)

    select case(tg%type)
    case(oct_tg_groundstate)
      call target_end_groundstate()
    case(oct_tg_excited)
      call target_end_excited()
    case(oct_tg_exclude_state)
      call target_end_exclude()
    case(oct_tg_gstransformation)
      call target_end_gstransformation()
    case(oct_tg_userdefined) 
      call target_end_userdefined()
    case(oct_tg_jdensity)
      call target_end_density(tg)
    case(oct_tg_local)
      call target_end_local(tg)
    case(oct_tg_td_local)
      call target_end_tdlocal(tg)
    case(oct_tg_hhg)
      call target_end_hhg(tg)
    case(oct_tg_hhgnew)
      call target_end_hhgnew(tg, oct)
    case(oct_tg_velocity)
      call target_end_velocity(tg, oct)
    case(oct_tg_classical)
      call target_end_classical(tg)
    end select

    POP_SUB(target_end)
  end subroutine target_end
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_output(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    PUSH_SUB(target_output)

    select case(tg%type)
    case(oct_tg_groundstate)
      call target_output_groundstate(tg, gr, dir, geo, outp)
    case(oct_tg_excited)
      call target_output_excited(tg, gr, dir, geo, outp)
    case(oct_tg_exclude_state)
      call target_output_exclude(tg, gr, dir, geo, outp)
    case(oct_tg_gstransformation)
      call target_output_gstransformation(tg, gr, dir, geo, outp)
    case(oct_tg_userdefined) 
      call target_output_userdefined(tg, gr, dir, geo, outp)
    case(oct_tg_jdensity)
      call target_output_density(tg, gr, dir, geo, outp)
    case(oct_tg_local)
      call target_output_local(tg, gr, dir, geo, outp)
    case(oct_tg_td_local)
      call target_output_tdlocal(tg, gr, dir, geo, outp)
    case(oct_tg_hhg)
      call target_output_hhg(tg, gr, dir, geo, outp)
    case(oct_tg_hhgnew)
      call target_output_hhg(tg, gr, dir, geo, outp)
    case(oct_tg_velocity)
      call target_output_velocity(tg, gr, dir, geo, outp)
    case(oct_tg_classical)
      call target_output_classical
    end select
    
    POP_SUB(target_output)
  end subroutine target_output
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculates, at a given point in time marked by the integer
  !! index, the integrand of the target functional:
  !! <Psi(t)|\hat{O}(t)|Psi(t)>.
  subroutine target_tdcalc(tg, hm, gr, geo, psi, time, max_time)
    type(target_t),      intent(inout) :: tg
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(states_t),      intent(inout) :: psi
    integer,             intent(in)    :: time
    integer,             intent(in)    :: max_time

    if(target_mode(tg)  /= oct_targetmode_td) return

    PUSH_SUB(target_tdcalc)

    tg%td_fitness(time) = M_ZERO

    select case(tg%type)
    case(oct_tg_hhgnew)
      call target_tdcalc_hhgnew(tg, gr, psi, time, max_time)
    case(oct_tg_velocity)
      call target_tdcalc_velocity(tg, hm, gr, geo, psi, time, max_time)
    case(oct_tg_td_local)
      call target_tdcalc_tdlocal(tg, gr, psi, time)
    case(oct_tg_hhg)
      call target_tdcalc_hhg(tg, hm, gr, geo, psi, time)
    case(oct_tg_jdensity)
      call target_tdcalc_density(tg, gr, psi, time)
    case default
      message(1) = 'Error in target.target_tdcalc: default.'
      call messages_fatal(1)
    end select

    POP_SUB(target_tdcalc)
  end subroutine target_tdcalc
  ! ----------------------------------------------------------------------



  ! ---------------------------------------------------------------
  !> Calculates the inhomogeneous term that appears in the equation
  !! for chi, and places it into inh.
  subroutine target_inh(psi, gr, tg, time, inh, iter)
    type(states_t),    intent(inout)     :: psi
    type(grid_t),      intent(in)        :: gr
    type(target_t),    intent(inout)     :: tg
    FLOAT,             intent(in)        :: time
    type(states_t),    intent(inout)     :: inh
    integer,           intent(in)        :: iter
 
    integer :: ik, ist, ip, idim
    CMPLX :: gvec(MAX_DIM)

    PUSH_SUB(target_inh)

    select case(tg%type)
    case(oct_tg_td_local)
      call target_build_tdlocal(tg, gr, time)
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
        inh%zpsi(ip, idim, ist, ik) = - psi%occ(ist, ik) * tg%rho(ip) * psi%zpsi(ip, idim, ist, ik)
      end forall

    case(oct_tg_hhgnew)
      gvec(1:gr%sb%dim) = real(tg%gvec(iter+1, 1:gr%sb%dim), REAL_PRECISION)
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
        inh%zpsi(ip, idim, ist, ik) = &
           - psi%occ(ist, ik) * M_TWO * sum(tg%grad_local_pot(1, ip, 1:gr%sb%dim) * gvec(1:gr%sb%dim)) * &
           psi%zpsi(ip, idim, ist, ik)
      end forall

    case(oct_tg_velocity)
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
         inh%zpsi(ip, idim, ist, ik) = - psi%occ(ist, ik) * tg%rho(ip) * psi%zpsi(ip, idim, ist, ik)
      end forall
   
    case(oct_tg_jdensity)
      if (abs(nint(time/tg%dt)) >= tg%strt_iter_curr_tg) then
        inh%zpsi =  -chi_current(tg, gr, psi)
      else
        inh%zpsi = M_ZERO
      end if     

    case default
      write(message(1),'(a)') 'Internal error in target_inh'
      call messages_fatal(1)
  
    end select

    POP_SUB(target_inh)
  end subroutine target_inh
  !----------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculates the J1 functional, i.e.:
  !! <Psi(T)|\hat{O}|Psi(T) in the time-independent
  !! case, or else \int_0^T dt <Psi(t)|\hat{O}(t)|Psi(t) in 
  !! the time-dependent case.
  FLOAT function target_j1(tg, gr, qcpsi, geo) result(j1)
    type(target_t), intent(inout)   :: tg
    type(grid_t),   intent(inout)   :: gr
    type(opt_control_state_t), intent(inout)   :: qcpsi
    type(geometry_t), intent(in), optional :: geo

    type(states_t), pointer :: psi

    psi => opt_control_point_qs(qcpsi)

    PUSH_SUB(target_j1)

    select case(tg%type)
    case(oct_tg_groundstate)
      j1 = target_j1_groundstate(tg, gr, psi)
    case(oct_tg_excited)
      j1 = target_j1_excited(tg, gr, psi)
    case(oct_tg_gstransformation)
      j1 = target_j1_gstransformation(tg, gr, psi)
    case(oct_tg_userdefined)
      j1 = target_j1_userdefined(tg, gr, psi)
    case(oct_tg_jdensity)
      j1 = target_j1_density(gr, tg, psi)
    case(oct_tg_local)
      j1 = target_j1_local(gr, tg, psi)
    case(oct_tg_td_local)
      j1 = target_j1_tdlocal(tg)
    case(oct_tg_exclude_state)
      j1 = target_j1_exclude(gr, tg, psi)
    case(oct_tg_hhg)
      j1 = target_j1_hhg(tg)
    case(oct_tg_hhgnew)
      j1 = target_j1_hhgnew(gr, tg)
    case(oct_tg_velocity)
      j1 = target_j1_velocity(tg, geo)
    case(oct_tg_classical)
      j1 = target_j1_classical(tg, qcpsi)
    end select

    nullify(psi)
    POP_SUB(target_j1)
  end function target_j1
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculate |chi(T)> = \hat{O}(T) |psi(T)>
  subroutine target_chi(tg, gr, qcpsi_in, qcchi_out, geo)
    type(target_t),    intent(inout) :: tg
    type(grid_t),      intent(inout) :: gr
    type(opt_control_state_t), target, intent(inout) :: qcpsi_in
    type(opt_control_state_t), target, intent(inout) :: qcchi_out
    type(geometry_t),  intent(in)    :: geo

    type(states_t), pointer :: psi_in, chi_out
    PUSH_SUB(target_chi)

    psi_in => opt_control_point_qs(qcpsi_in)
    chi_out => opt_control_point_qs(qcchi_out)

    select case(tg%type)
    case(oct_tg_groundstate)
      call target_chi_groundstate(tg, gr, psi_in, chi_out)
    case(oct_tg_excited) 
      call target_chi_excited(tg, gr, psi_in, chi_out)
    case(oct_tg_gstransformation)
      call target_chi_gstransformation(tg, gr, psi_in, chi_out)
    case(oct_tg_userdefined)
      call target_chi_userdefined(tg, gr, psi_in, chi_out)
    case(oct_tg_jdensity)
      call target_chi_density(tg, gr, psi_in, chi_out)
    case(oct_tg_local)
      call target_chi_local(tg, gr, psi_in, chi_out)
    case(oct_tg_td_local)
      call target_chi_tdlocal(gr, chi_out)
    case(oct_tg_exclude_state)
      call target_chi_exclude(tg, gr, psi_in, chi_out)
    case(oct_tg_hhg)
      call target_chi_hhg(gr, chi_out)
    case(oct_tg_hhgnew)
      call target_chi_hhg(gr, chi_out)
    case(oct_tg_velocity)
      call target_chi_velocity(gr, tg, chi_out, geo)
    case(oct_tg_classical)
      call target_chi_classical(tg, qcpsi_in, qcchi_out, geo)
    end select

    nullify(psi_in)
    nullify(chi_out)
    POP_SUB(target_chi)
  end subroutine target_chi


  ! ----------------------------------------------------------------------
  integer pure function target_mode(tg)
    type(target_t), intent(in) :: tg

    select case(tg%type)
    case(oct_tg_td_local, oct_tg_hhg, oct_tg_velocity, oct_tg_hhgnew)
      target_mode = oct_targetmode_td
    case default
      target_mode = oct_targetmode_static
    end select

    ! allow specific current functionals to be td
    ! Attention: yet combined with static density target,
    ! the total target is considered td.
    select case(tg%curr_functional)
    case(oct_curr_square_td) 
      target_mode = oct_targetmode_td
    end select

  end function target_mode


  ! ----------------------------------------------------------------------
  integer pure function target_type(tg)
    type(target_t), intent(in) :: tg
    target_type = tg%type
  end function target_type
  ! ----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  integer pure function target_curr_functional(tg)
    type(target_t), intent(in) :: tg
    target_curr_functional = tg%curr_functional
  end function target_curr_functional
  !-----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  logical pure function target_move_ions(tg)
    type(target_t), intent(in) :: tg
    target_move_ions = tg%move_ions
  end function target_move_ions
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  logical pure function is_spatial_curr_wgt(tg)
    type(target_t), intent(in) :: tg
    is_spatial_curr_wgt = associated(tg%spatial_curr_wgt)
  end function is_spatial_curr_wgt
  ! ----------------------------------------------------------------------

#include "target_density_inc.F90"
#include "target_velocity_inc.F90"
#include "target_hhg_inc.F90"
#include "target_groundstate_inc.F90"
#include "target_excited_inc.F90"
#include "target_gstransformation_inc.F90"
#include "target_exclude_inc.F90"
#include "target_userdefined_inc.F90"
#include "target_local_inc.F90"
#include "target_tdlocal_inc.F90"
#include "target_classical_inc.F90"

end module target_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
