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
  use excited_states_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
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
  use states_oct_m
  use states_calc_oct_m
  use states_dim_oct_m
  use states_restart_oct_m
  use td_calc_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::                      &
       floquet_init,             &
       floquet_hamiltonians_init

contains

  subroutine floquet_init(this,dim)
    type(floquet_t),    intent(out)    :: this
    integer :: dim ! the standard dimension of the groundstate

    FLOAT :: time_step

    PUSH_SUB(floquet_init)

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

    ! re-read time stepfrom input
    call parse_variable('TDTimeStep', M_ONE, time_step, unit = units_inp%time)
    this%interval = int(this%dt/time_step)
    this%ncycle = this%interval*this%nT

    call messages_print_var_value(stdout,'Steps in Floquet time-sampling interval',  this%interval)
    call messages_print_var_value(stdout,'Steps in Floquet time-sampling cycle',  this%ncycle)

   POP_SUB(floquet_init)

  end subroutine floquet_init

  subroutine floquet_hamiltonians_init(this,gr, st, ks, iter)
    type(hamiltonian_t), intent(inout) :: this ! this is not great, as everyhting should be within the floquet_t
    type(grid_t),      intent(inout)   :: gr
    type(states_t),    intent(inout)   :: st !< at iter=0 this is the ggroundstate
    type(v_ks_t),      intent(in)      :: ks
    integer,           intent(in)      :: iter

    CMPLX, allocatable :: hmss(:,:), psi(:,:,:), hpsi(:,:,:), temp_state1(:,:), temp_state2(:,:)
    CMPLX, allocatable :: HFloquet(:,:,:), HFloq_eff(:,:), temp(:,:)
    FLOAT, allocatable :: eigenval(:), bands(:,:)
    character(len=80) :: filename
    integer :: it, nT, ik, ist, jst, in, im, inm, file, idim, nik, ik_count
    integer ::  m0, n0, n1, nst, ii, jj, lim_nst
    type(mesh_t) :: mesh
    type(states_t) :: hm_st
    FLOAT :: time_step

    PUSH_SUB(floquet_hamiltonian_init)

    mesh = gr%der%mesh
    nst = st%nst

    !for now no domain distributionallowed
    ASSERT(mesh%np == mesh%np_global)

   ! this is used to initialize the hpsi (more effiecient ways?)
    call states_copy(hm_st, st)

    ! the Hamiltonain gets assigned an array of td-Hamiltonians
    ! this is a bit recursive, so maybe there should be a Flqoeut moduel or something
    nullify(this%td_hm)
    SAFE_ALLOCATE(this%td_hm(1:this%F%nT))
    
    POP_SUB(floquet_hamiltonian_init)
    
  end subroutine floquet_hamiltonians_init

end module floquet_oct_m
