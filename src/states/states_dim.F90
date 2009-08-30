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
  use crystal_m
  use datasets_m
  use distributed_m
  use geometry_m
  use global_m
  use grid_m
  use io_function_m
  use io_m
  use lalg_basic_m
  use loct_m
  use loct_parser_m
  use math_m
  use messages_m
  use mesh_function_m
  use mesh_m
  use mpi_m
  use mpi_lib_m
  use multicomm_m
  use profiling_m
  use simul_box_m
  use units_m
  use varinfo_m
  use species_m

  implicit none

  private

  public ::                           &
    states_dim_t,                     &
    states_dim_copy,                  &
    states_dim_end,                   &
    is_spin_down,                     &
    is_spin_up,                       &
    states_choose_kpoints,            &
    states_dim_get_spin_index,        &
    kpoints_write_info,               &
    kpoints_distribute,               &
    kpoint_is_gamma,                  &
    kpoint_index

  ! Parameters...
  integer, public, parameter :: &
    UNPOLARIZED    = 1,         &
    SPIN_POLARIZED = 2,         &
    SPINORS        = 3

  ! Spin-polarized k-indices for non-periodic systems.
  integer, public, parameter :: &
    SPIN_DOWN = 1,              &
    SPIN_UP   = 2
  
  type states_dim_t
    integer :: dim                  ! Dimension of the state (one or two for spinors)
    integer :: nik                  ! Number of irreducible subspaces
    integer :: nik_axis(MAX_DIM)    ! Number of kpoints per axis
    integer :: ispin                ! spin mode (unpolarized, spin polarized, spinors)
    integer :: nspin                ! dimension of rho (1, 2 or 4)
    integer :: spin_channels        ! 1 or 2, whether spin is or not considered.
    logical :: cdft                 ! Are we using Current-DFT or not?
    FLOAT, pointer :: kpoints(:,:)  ! obviously the kpoints
    FLOAT, pointer :: kweights(:)   ! weights for the kpoint integrations
    type(distributed_t) :: kpt
    integer :: block_size
    integer :: window_size
  end type states_dim_t

contains

  subroutine states_dim_copy(dout, din)
    type(states_dim_t), intent(out) :: dout
    type(states_dim_t), intent(in)  :: din

    call push_sub('states_dim.states_dim_copy')

    dout%dim            = din%dim
    dout%nik            = din%nik
    dout%nik_axis       = din%nik_axis
    dout%ispin          = din%ispin
    dout%nspin          = din%nspin
    dout%spin_channels  = din%spin_channels
    dout%cdft           = din%cdft
    dout%block_size     = din%block_size
    dout%window_size    = din%window_size

    call loct_pointer_copy(dout%kpoints, din%kpoints)
    call loct_pointer_copy(dout%kweights, din%kweights)

    call distributed_copy(din%kpt, dout%kpt)

    call pop_sub()
  end subroutine states_dim_copy

  ! ---------------------------------------------------------
  subroutine states_dim_end(d)
    type(states_dim_t), intent(inout) :: d

    call distributed_end(d%kpt)

    SAFE_DEALLOCATE_P(d%kpoints)
    SAFE_DEALLOCATE_P(d%kweights)

  end subroutine states_dim_end


  ! ---------------------------------------------------------
  ! Returns if k-point ik denotes spin-up, in spin-polarized case.
  logical function is_spin_up(ik)
    integer, intent(in) :: ik

    call push_sub('states_dim.is_spin_up')

    is_spin_up = even(ik)

    call pop_sub()
  end function is_spin_up


  ! ---------------------------------------------------------
  ! Returns if k-point ik denotes spin-down, in spin-polarized case.
  logical function is_spin_down(ik)
    integer, intent(in) :: ik

    call push_sub('states_dim.is_spin_down')

    is_spin_down = odd(ik)

    call pop_sub()
  end function is_spin_down

  integer function states_dim_get_spin_index(this, iq) result(index)
    type(states_dim_t), intent(in) :: this
    integer,            intent(in) :: iq
    
    if(this%ispin == SPIN_POLARIZED) then
      index = 1 + mod(iq - 1, 2)
    else
      index = 1
    end if
    
  end function states_dim_get_spin_index
  
  subroutine kpoints_distribute(this, mc)
    type(states_dim_t), intent(inout) :: this
    type(multicomm_t),  intent(in)    :: mc

    call distributed_init(this%kpt, this%nik, mc, P_STRATEGY_KPOINTS, "K points")

  end subroutine kpoints_distribute
  
#include "states_kpoints_inc.F90"

end module states_dim_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
