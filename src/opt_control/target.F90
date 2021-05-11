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

module target_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use density_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use excited_states_oct_m
  use fft_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use opt_control_global_oct_m
  use opt_control_state_oct_m
  use output_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use species_oct_m
  use spectrum_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use string_oct_m
  use td_calc_oct_m
  use td_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

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
    oct_tg_classical        = 13,     &
    oct_tg_spin             = 14

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
    type(states_elec_t) :: st
    type(excited_states_t) :: est
    FLOAT, allocatable :: rho(:)
    FLOAT, allocatable :: td_fitness(:)
    character(len=200) :: td_local_target
    character(len=80) :: excluded_states_list
    character(len=4096) :: vel_input_string
    character(len=4096) :: classical_input_string
    character(len=1024), allocatable :: vel_der_array(:,:)
    character(len=1024), allocatable :: mom_der_array(:,:)
    character(len=1024), allocatable :: pos_der_array(:,:)
    FLOAT, allocatable :: grad_local_pot(:,:,:)
    logical :: move_ions
    integer :: hhg_nks
    integer, allocatable :: hhg_k(:)
    FLOAT,   allocatable :: hhg_alpha(:)
    FLOAT,   allocatable :: hhg_a(:)
    FLOAT   :: hhg_w0
    FLOAT   :: dt
    integer :: curr_functional
    FLOAT   :: density_weight
    FLOAT   :: curr_weight
    integer :: strt_iter_curr_tg
    FLOAT, allocatable :: spatial_curr_wgt(:)
    character(len=1000) :: plateau_string
    CMPLX, allocatable :: acc(:, :)
    CMPLX, allocatable :: vel(:, :)
    CMPLX, allocatable :: gvec(:, :)
    FLOAT, allocatable :: alpha(:)
    CMPLX :: spin_matrix(2, 2)
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
  !> This just copies the states_elec_t variable present in target, into st.
  subroutine target_get_state(tg, st)
    type(target_t),      intent(in)    :: tg
    type(states_elec_t), intent(inout) :: st

    PUSH_SUB(target_get_state)
    call states_elec_copy(st, tg%st)

    POP_SUB(target_get_state)
  end subroutine target_get_state
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> The target is initialized, mainly by reading from the inp file.
  subroutine target_init(gr, kpoints, namespace, space, ions, qcs, td, w0, tg, oct, ep, mc)
    type(grid_t),                intent(in)    :: gr
    type(kpoints_t),             intent(in)    :: kpoints
    type(namespace_t),           intent(in)    :: namespace
    type(space_t),               intent(in)    :: space
    type(ions_t),                intent(in)    :: ions
    type(opt_control_state_t),   intent(inout) :: qcs
    type(td_t),                  intent(in)    :: td
    FLOAT,                       intent(in)    :: w0
    type(target_t),              intent(inout) :: tg
    type(oct_t),                 intent(in)    :: oct
    type(epot_t),                intent(inout) :: ep
    type(multicomm_t),           intent(in)    :: mc

    integer :: ierr
    type(states_elec_t), pointer :: stin
    type(restart_t) :: restart

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
    !% (Experimental) The target operator is an "excited state". This means that the target operator
    !% is a linear combination of Slater determinants, each one formed by replacing
    !% in the ground-state Slater determinant one occupied state with one excited
    !% state (<i>i.e.</i> "single excitations"). The description of which excitations are
    !% used, and with which weights, should be given in a file called
    !% <tt>oct-excited-state-target</tt>.
    !% See the documentation of subroutine <tt>excited_states_elec_init</tt> in the source
    !% code in order to use this feature.
    !%Option oct_tg_gstransformation 3
    !% The target operator is a projection operator on a transformation of the ground-state 
    !% orbitals defined by the block <tt>OCTTargetTransformStates</tt>.
    !%Option oct_tg_userdefined 4
    !% (Experimental) Allows to define target state by using <tt>OCTTargetUserdefined</tt>.
    !%Option oct_tg_jdensity 5
    !% (Experimental)
    !%Option oct_tg_local 6
    !% (Experimental) The target operator is a local operator.
    !%Option oct_tg_td_local 7
    !% (Experimental) The target operator is a time-dependent local operator.
    !%Option oct_tg_exclude_state 8
    !% (Experimental) Target operator is the projection onto the complement of a given state, given by the
    !% block <tt>OCTTargetTransformStates</tt>. This means that the target operator is the unity
    !% operator minus the projector onto that state.
    !%Option oct_tg_hhg 9
    !% (Experimental) The target is the optimization of the HHG yield. You must supply the <tt>OCTOptimizeHarmonicSpectrum</tt>
    !% block, and it attempts to optimize the maximum of the spectrum around each harmonic peak. You may
    !% use only one of the gradient-less optimization schemes.
    !%Option oct_tg_velocity 10
    !% (Experimental) The target is a function of the velocities of the nuclei at the end of the influence of
    !% the external field, defined by <tt>OCTVelocityTarget</tt>
    !%Option oct_tg_hhgnew 12
    !% (Experimental) The target is the optimization of the HHG yield. You must supply the
    !% <tt>OCTHarmonicWeight</tt> string. It attempts to optimize the integral of the harmonic spectrum multiplied
    !% by some user-defined weight function.
    !%Option oct_tg_classical 13
    !% (Experimental)
    !%Option oct_tg_spin 14
    !% (Experimental)
    !%End
    call parse_variable(namespace, 'OCTTargetOperator', oct_tg_gstransformation, tg%type)
      if(tg%type == oct_tg_excited) call messages_experimental('OCTTargetOperator = oct_tg_excited')
      if(tg%type == oct_tg_userdefined) call messages_experimental('OCTTargetOperator = oct_tg_userdefined')
      if(tg%type == oct_tg_jdensity) call messages_experimental('OCTTargetOperator = oct_tg_jdensity')
      if(tg%type == oct_tg_local) call messages_experimental('OCTTargetOperator = oct_tg_local')
      if(tg%type == oct_tg_td_local) call messages_experimental('OCTTargetOperator = oct_tg_td_local')
      if(tg%type == oct_tg_exclude_state) call messages_experimental('OCTTargetOperator = oct_tg_exclude_state')
      if(tg%type == oct_tg_hhg) call messages_experimental('OCTTargetOperator = oct_tg_hhg')
      if(tg%type == oct_tg_velocity) call messages_experimental('OCTTargetOperator = oct_tg_velocity')
      if(tg%type == oct_tg_hhgnew) call messages_experimental('OCTTargetOperator = oct_tg_hhgnew')
      if(tg%type == oct_tg_classical) call messages_experimental('OCTTargetOperator = oct_tg_classical')
      if(tg%type == oct_tg_spin) call messages_experimental('OCTTargetOperator = oct_tg_spin')


    if(.not.varinfo_valid_option('OCTTargetOperator', tg%type)) &
      call messages_input_error(namespace, 'OCTTargetOperator')

    call states_elec_copy(tg%st, stin)
    call states_elec_deallocate_wfns(tg%st)
    call states_elec_allocate_wfns(tg%st, gr%mesh, TYPE_CMPLX)
    call restart_init(restart, namespace, RESTART_GS, RESTART_TYPE_LOAD, mc, ierr, mesh=gr%mesh, exact=.true.)
    if(ierr /= 0) then
      message(1) = "Could not read gs for OCTTargetOperator."
      call messages_fatal(1)
    end if

    select case(tg%type)
    case(oct_tg_groundstate)
      call target_init_groundstate(gr%mesh, namespace, space, tg, td, restart, kpoints)
    case(oct_tg_excited)
      call messages_experimental('OCTTargetOperator = oct_tg_excited')
      call target_init_excited(gr%mesh, namespace, space, tg, td, restart, kpoints)
    case(oct_tg_exclude_state)
      call target_init_exclude(gr%mesh, namespace, space, tg, td, restart, kpoints)
    case(oct_tg_gstransformation)
      call target_init_gstransformation(gr, namespace, space, tg, td, restart, kpoints)
    case(oct_tg_userdefined) 
      call target_init_userdefined(gr, namespace, tg, td)
    case(oct_tg_jdensity)
      call target_init_density(gr, kpoints, namespace, space, tg, stin, td, restart)
    case(oct_tg_local)
      call target_init_local(gr, namespace, tg, td)
    case(oct_tg_td_local)
      call target_init_tdlocal(gr, namespace, tg, td)
    case(oct_tg_hhg)
      call target_init_hhg(tg, namespace, td, w0)
    case(oct_tg_hhgnew)
      call messages_experimental('OCTTargetOperator = oct_tg_hhgnew')
      call target_init_hhgnew(gr, namespace, tg, td, ions, ep)
    case(oct_tg_velocity)
      call target_init_velocity(gr, namespace, ions, tg, oct, td, ep)
    case(oct_tg_classical)
      call messages_experimental('OCTTargetOperator = oct_tg_classical')
      call target_init_classical(ions, namespace, tg, td, oct)
    case(oct_tg_spin)
      call messages_experimental('OCTTargetOperator = oct_tg_spin')
      call target_init_spin(tg, namespace)
    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call messages_fatal(1)
    end select

    call restart_end(restart)

    nullify(stin)
    POP_SUB(target_init)
  end subroutine target_init
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_end(tg, oct)
    type(target_t), intent(inout) :: tg
    type(oct_t), intent(in)       :: oct

    PUSH_SUB(target_end)

    call states_elec_end(tg%st)

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
  subroutine target_output(tg, namespace, space, gr, dir, ions, hm, outp)
    type(target_t),           intent(inout) :: tg
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr
    character(len=*),         intent(in)    :: dir
    type(ions_t),             intent(in)    :: ions
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(output_t),           intent(in)    :: outp

    PUSH_SUB(target_output)

    select case(tg%type)
    case(oct_tg_groundstate)
      call target_output_groundstate(tg, namespace, space, gr, dir, ions, hm, outp)
    case(oct_tg_excited)
      call target_output_excited(tg, namespace, space, gr, dir, ions, hm, outp)
    case(oct_tg_exclude_state)
      call target_output_exclude(tg, namespace, space, gr, dir, ions, hm, outp)
    case(oct_tg_gstransformation)
      call target_output_gstransformation(tg, namespace, space, gr, dir, ions, hm, outp)
    case(oct_tg_userdefined) 
      call target_output_userdefined(tg, namespace, space, gr, dir, ions, hm, outp)
    case(oct_tg_jdensity)
      call target_output_density(tg, namespace, space, gr%mesh, dir, ions, outp)
    case(oct_tg_local)
      call target_output_local(tg, namespace, space, gr%mesh, dir, ions, outp)
    case(oct_tg_td_local)
      call target_output_tdlocal(tg, namespace, space, gr, dir, ions, outp)
    case(oct_tg_hhg)
      call target_output_hhg(tg, namespace, space, gr, dir, ions, hm, outp)
    case(oct_tg_hhgnew)
      call target_output_hhg(tg, namespace, space, gr, dir, ions, hm, outp)
    case(oct_tg_velocity)
      call target_output_velocity(tg, namespace, space, gr, dir, ions, hm, outp)
    case(oct_tg_classical)
      call target_output_classical()
    end select
    
    POP_SUB(target_output)
  end subroutine target_output
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculates, at a given point in time marked by the integer
  !! index, the integrand of the target functional:
  !! <Psi(t)|\hat{O}(t)|Psi(t)>.
  subroutine target_tdcalc(tg, namespace, space, hm, gr, ions, psi, time, max_time)
    type(target_t),           intent(inout) :: tg
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(inout) :: ions
    type(states_elec_t),      intent(inout) :: psi
    integer,                  intent(in)    :: time
    integer,                  intent(in)    :: max_time

    if(target_mode(tg)  /= oct_targetmode_td) return

    PUSH_SUB(target_tdcalc)

    tg%td_fitness(time) = M_ZERO

    select case(tg%type)
    case(oct_tg_hhgnew)
      call target_tdcalc_hhgnew(tg, gr, psi, time, max_time)
    case(oct_tg_velocity)
      call target_tdcalc_velocity(tg, hm, gr, ions, psi, time, max_time)
    case(oct_tg_td_local)
      call target_tdcalc_tdlocal(tg, gr, psi, time)
    case(oct_tg_hhg)
      call target_tdcalc_hhg(tg, namespace, space, hm, gr, ions, psi, time)
    case(oct_tg_jdensity)
      call target_tdcalc_density(tg, gr, hm%kpoints, psi, time)
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
  subroutine target_inh(psi, gr, kpoints, tg, time, inh, iter)
    type(states_elec_t), intent(inout)     :: psi
    type(grid_t),        intent(in)        :: gr
    type(kpoints_t),     intent(in)        :: kpoints
    type(target_t),      intent(inout)     :: tg
    FLOAT,               intent(in)        :: time
    type(states_elec_t), intent(inout)     :: inh
    integer,             intent(in)        :: iter
 
    integer :: ik, ist, ip, idim, ib
    CMPLX, allocatable :: zpsi(:)
    CMPLX :: gvec(MAX_DIM)

    PUSH_SUB(target_inh)

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np))
    
    select case(tg%type)
    case(oct_tg_td_local)

      call target_build_tdlocal(tg, gr, time)

      do ik = inh%d%kpt%start, inh%d%kpt%end
        do ist = inh%st_start, inh%st_end
          do idim = 1, inh%d%dim
            call states_elec_get_state(psi, gr%mesh, idim, ist, ik, zpsi)
            zpsi(1:gr%mesh%np) = -psi%occ(ist, ik)*tg%rho(1:gr%mesh%np)*zpsi(1:gr%mesh%np)
            call states_elec_set_state(inh, gr%mesh, idim, ist, ik, zpsi)
          end do
        end do
      end do
      
    case(oct_tg_hhgnew)
      gvec(1:gr%sb%dim) = TOFLOAT(tg%gvec(iter + 1, 1:gr%sb%dim))

      do ik = inh%d%kpt%start, inh%d%kpt%end
        do ist = inh%st_start, inh%st_end
          do idim = 1, inh%d%dim
            call states_elec_get_state(psi, gr%mesh, idim, ist, ik, zpsi)
            do ip = 1, gr%mesh%np
              zpsi(ip) = -psi%occ(ist, ik)*M_TWO*sum(tg%grad_local_pot(1, ip, 1:gr%sb%dim)*gvec(1:gr%sb%dim))*zpsi(ip)
            end do
            call states_elec_set_state(inh, gr%mesh, idim, ist, ik, zpsi)
          end do
        end do
      end do

    case(oct_tg_velocity)

      do ik = inh%d%kpt%start, inh%d%kpt%end
        do ist = inh%st_start, inh%st_end
          do idim = 1, inh%d%dim
            call states_elec_get_state(psi, gr%mesh, idim, ist, ik, zpsi)
            do ip = 1, gr%mesh%np
              zpsi(ip) = -psi%occ(ist, ik)*tg%rho(ip)*zpsi(ip)
            end do
            call states_elec_set_state(inh, gr%mesh, idim, ist, ik, zpsi)
          end do
        end do
      end do
   
    case(oct_tg_jdensity)

      do ik = inh%d%kpt%start, inh%d%kpt%end
        do ib = inh%group%block_start, inh%group%block_end
          call batch_set_zero(inh%group%psib(ib, ik))
        end do
      end do
        
      if (abs(nint(time/tg%dt)) >= tg%strt_iter_curr_tg) then
        call chi_current(tg, gr, kpoints, CNST(-1.0), psi, inh)
      end if     

    case default
      write(message(1),'(a)') 'Internal error in target_inh'
      call messages_fatal(1)
  
    end select

    SAFE_DEALLOCATE_A(zpsi)
    
    POP_SUB(target_inh)
  end subroutine target_inh
  !----------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculates the J1 functional, i.e.:
  !! <Psi(T)|\hat{O}|Psi(T) in the time-independent
  !! case, or else \int_0^T dt <Psi(t)|\hat{O}(t)|Psi(t) in 
  !! the time-dependent case.
  FLOAT function target_j1(tg, namespace, gr, kpoints, qcpsi, ions) result(j1)
    type(target_t),             intent(inout)   :: tg
    type(namespace_t),          intent(in)      :: namespace
    type(grid_t),               intent(in)      :: gr
    type(kpoints_t),            intent(in)      :: kpoints
    type(opt_control_state_t),  intent(inout)   :: qcpsi
    type(ions_t),     optional, intent(in)      :: ions

    type(states_elec_t), pointer :: psi

    psi => opt_control_point_qs(qcpsi)

    PUSH_SUB(target_j1)

    select case(tg%type)
    case(oct_tg_groundstate)
      j1 = target_j1_groundstate(tg, gr, psi)
    case(oct_tg_excited)
      j1 = target_j1_excited(tg, namespace, gr, psi)
    case(oct_tg_gstransformation)
      j1 = target_j1_gstransformation(tg, gr, psi)
    case(oct_tg_userdefined)
      j1 = target_j1_userdefined(tg, gr, psi)
    case(oct_tg_jdensity)
      j1 = target_j1_density(gr, kpoints, tg, psi)
    case(oct_tg_local)
      j1 = target_j1_local(gr%mesh, tg, psi)
    case(oct_tg_td_local)
      j1 = target_j1_tdlocal(tg)
    case(oct_tg_exclude_state)
      j1 = target_j1_exclude(gr, tg, psi)
    case(oct_tg_hhg)
      j1 = target_j1_hhg(tg, namespace)
    case(oct_tg_hhgnew)
      j1 = target_j1_hhgnew(gr, tg)
    case(oct_tg_velocity)
      j1 = target_j1_velocity(tg, ions)
    case(oct_tg_classical)
      j1 = target_j1_classical(tg, qcpsi)
    case(oct_tg_spin)
      j1 = target_j1_spin(tg, gr, psi)
    end select

    nullify(psi)
    POP_SUB(target_j1)
  end function target_j1
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculate |chi(T)> = \hat{O}(T) |psi(T)>
  subroutine target_chi(tg, namespace, gr, kpoints, qcpsi_in, qcchi_out, ions)
    type(target_t),                    intent(inout) :: tg
    type(namespace_t),                 intent(in)    :: namespace
    type(grid_t),                      intent(in)    :: gr
    type(kpoints_t),                   intent(in)    :: kpoints
    type(opt_control_state_t), target, intent(inout) :: qcpsi_in
    type(opt_control_state_t), target, intent(inout) :: qcchi_out
    type(ions_t),                      intent(in)    :: ions

    FLOAT, pointer :: q(:, :), p(:, :)
    type(states_elec_t), pointer :: psi_in, chi_out
    PUSH_SUB(target_chi)

    psi_in => opt_control_point_qs(qcpsi_in)
    chi_out => opt_control_point_qs(qcchi_out)

    select case(tg%type)
    case(oct_tg_groundstate)

      call target_chi_groundstate(tg, gr, psi_in, chi_out)
    case(oct_tg_excited) 
      call target_chi_excited(tg, namespace, gr, psi_in, chi_out)
    case(oct_tg_gstransformation)
      call target_chi_gstransformation(tg, gr, psi_in, chi_out)
    case(oct_tg_userdefined)
      call target_chi_userdefined(tg, gr, psi_in, chi_out)
    case(oct_tg_jdensity)
      call target_chi_density(tg, gr, kpoints, psi_in, chi_out)
    case(oct_tg_local)
      call target_chi_local(tg, gr%mesh, psi_in, chi_out)
    case(oct_tg_td_local)
      call target_chi_tdlocal(chi_out)
    case(oct_tg_exclude_state)
      call target_chi_exclude(tg, gr, psi_in, chi_out)
    case(oct_tg_hhg)
      call target_chi_hhg(chi_out)
    case(oct_tg_hhgnew)
      call target_chi_hhg(chi_out)
    case(oct_tg_velocity)
      call target_chi_velocity(gr, tg, chi_out, ions)
    case(oct_tg_classical)
      call target_chi_classical(tg, qcpsi_in, qcchi_out, ions)
    case(oct_tg_spin)
      call target_chi_spin(tg, gr, psi_in, chi_out)
    end select

    ! Unless the target is "classical", the co-state classical variables are zero at time t=T.
    if(tg%type .ne. oct_tg_classical ) then
      q => opt_control_point_q(qcchi_out)
      p => opt_control_point_p(qcchi_out)
      q = M_ZERO
      p = M_ZERO
      nullify(q)
      nullify(p)
    end if

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

    is_spatial_curr_wgt = allocated(tg%spatial_curr_wgt)

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
#include "target_spin_inc.F90"

end module target_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
