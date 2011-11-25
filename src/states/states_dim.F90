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

module states_dim_m
  use blas_m
  use calc_mode_m
  use datasets_m
  use distributed_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use kpoints_m
  use lalg_basic_m
  use loct_m
  use math_m
  use messages_m
  use mesh_m
  use mesh_function_m
  use mpi_m
  use mpi_lib_m
  use multicomm_m
  use parser_m
  use profiling_m
  use simul_box_m
  use species_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private

  public ::                           &
    states_dim_t,                     &
    states_dim_null,                  &
    states_dim_copy,                  &
    states_dim_end,                   &
    is_spin_down,                     &
    is_spin_up,                       &
    states_choose_kpoints,            &
    states_dim_get_spin_index,        &
    states_dim_get_kpoint_index,      &
    kpoints_distribute

  ! Parameters...
  integer, public, parameter :: &
    UNPOLARIZED    = 1,         &
    SPIN_POLARIZED = 2,         &
    SPINORS        = 3

  ! Spin-polarized k-indices for non-periodic systems.
  integer, public, parameter :: &
    SPIN_DOWN = 1,              &
    SPIN_UP   = 2

  ! Orthogonalization methods
  integer, public, parameter ::   &
    ORTH_GS     = 1,               &
    ORTH_PAR_GS = 2,               &
    ORTH_MGS    = 3,               &
    ORTH_QR     = 4

  type states_dim_t
    integer :: dim                  !< Dimension of the state (one or two for spinors)
    integer :: nik                  !< Number of irreducible subspaces
    integer :: nik_axis(MAX_DIM)    !< Number of k-points per axis
    integer :: ispin                !< spin mode (unpolarized, spin-polarized, spinors)
    integer :: nspin                !< dimension of rho (1, 2 or 4)
    integer :: spin_channels        !< 1 or 2, whether spin is or not considered.
    logical :: cdft                 !< Are we using current-DFT or not?
    FLOAT, pointer :: kweights(:)   !< weights for the k-point integrations
    type(distributed_t) :: kpt
    integer :: block_size
    integer :: orth_method
    logical :: pack_states
    FLOAT   :: cl_states_mem
  end type states_dim_t

contains

  ! ---------------------------------------------------------
  subroutine states_dim_null(this)
    type(states_dim_t), intent(inout) :: this

    PUSH_SUB(states_dim_null)

    nullify(this%kweights)
    call distributed_nullify(this%kpt)

    POP_SUB(states_dim_null)
  end subroutine states_dim_null


  ! ---------------------------------------------------------
  subroutine states_dim_copy(dout, din)
    type(states_dim_t), intent(out) :: dout
    type(states_dim_t), intent(in)  :: din

    PUSH_SUB(states_dim_copy)

    dout%dim            = din%dim
    dout%nik            = din%nik
    dout%nik_axis       = din%nik_axis
    dout%ispin          = din%ispin
    dout%nspin          = din%nspin
    dout%spin_channels  = din%spin_channels
    dout%cdft           = din%cdft
    dout%block_size     = din%block_size
    dout%pack_states    = din%pack_states
    dout%cl_states_mem  = din%cl_states_mem

    call loct_pointer_copy(dout%kweights, din%kweights)

    call distributed_copy(din%kpt, dout%kpt)

    POP_SUB(states_dim_copy)
  end subroutine states_dim_copy


  ! ---------------------------------------------------------
  subroutine states_dim_end(dim)
    type(states_dim_t), intent(inout) :: dim

    PUSH_SUB(states_dim_end)

    call distributed_end(dim%kpt)

    SAFE_DEALLOCATE_P(dim%kweights)

    POP_SUB(states_dim_end)
  end subroutine states_dim_end


  ! ---------------------------------------------------------
  ! Returns true if k-point ik denotes spin-up, in spin-polarized case.
  logical pure function is_spin_up(ik)
    integer, intent(in) :: ik

    is_spin_up = even(ik)

  end function is_spin_up


  ! ---------------------------------------------------------
  ! Returns true if k-point ik denotes spin-down, in spin-polarized case.
  logical pure function is_spin_down(ik)
    integer, intent(in) :: ik

    is_spin_down = odd(ik)

  end function is_spin_down

  ! ---------------------------------------------------------
  integer pure function states_dim_get_spin_index(this, iq) result(index)
    type(states_dim_t), intent(in) :: this
    integer,            intent(in) :: iq

    if(this%ispin == SPIN_POLARIZED) then
      index = 1 + mod(iq - 1, 2)
    else
      index = 1
    end if

  end function states_dim_get_spin_index


  ! ---------------------------------------------------------
  integer pure function states_dim_get_kpoint_index(this, iq) result(index)
    type(states_dim_t), intent(in) :: this
    integer,            intent(in) :: iq
    
    if(this%ispin == SPIN_POLARIZED) then
      index = 1 + (iq - 1)/2
    else
      index = iq
    end if
    
  end function states_dim_get_kpoint_index


  ! ---------------------------------------------------------
  subroutine kpoints_distribute(this, mc)
    type(states_dim_t), intent(inout) :: this
    type(multicomm_t),  intent(in)    :: mc

    PUSH_SUB(kpoints_distribute)
    call distributed_init(this%kpt, this%nik, mc%group_comm(P_STRATEGY_KPOINTS), "k-points")

    POP_SUB(kpoints_distribute)
  end subroutine kpoints_distribute
  
  ! ---------------------------------------------------------
  subroutine states_choose_kpoints(dd, sb, geo)
    type(states_dim_t), intent(inout) :: dd
    type(simul_box_t),  intent(in)    :: sb
    type(geometry_t),   intent(in)    :: geo

    integer :: ik, iq

    PUSH_SUB(states_choose_kpoints)

    dd%nik = kpoints_number(sb%kpoints)

    if (dd%ispin == SPIN_POLARIZED) dd%nik = 2*dd%nik

    SAFE_ALLOCATE(dd%kweights(1:dd%nik))

    do iq = 1, dd%nik
      ik = states_dim_get_kpoint_index(dd, iq)
      dd%kweights(iq) = kpoints_get_weight(sb%kpoints, ik)
    end do

    if(in_debug_mode) call print_kpoints_debug
    POP_SUB(states_choose_kpoints)

  contains
    subroutine print_kpoints_debug
      integer :: iunit

      PUSH_SUB(states_choose_kpoints.print_kpoints_debug)
      
      call io_mkdir('debug/')
      iunit = io_open('debug/kpoints', action = 'write')
      call kpoints_write_info(sb%kpoints, iunit)      
      call io_close(iunit)

      POP_SUB(states_choose_kpoints.print_kpoints_debug)
    end subroutine print_kpoints_debug

  end subroutine states_choose_kpoints

end module states_dim_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
