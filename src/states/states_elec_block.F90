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

!> This module contains routines to multiply blocks of states by
!! blocks of states and general dense matrices.
!! THESE ROUTINES ARE DEPRECATED, THEY SHOULD NOT BE USED BY NEW CODE.
!! They are slow (too many data copies) and they are not maintained.

module states_elec_block_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use multicomm_oct_m
  use profiling_oct_m
  use states_elec_oct_m

  implicit none

  private

  public ::                    &
    states_elec_blockt_mul,         &
    states_elec_block_matr_mul_add


!   ! This type encapsulates an array of wavefunctions with numbers
!   ! st_start to st_end. It is used to pass arrays of states around
!   ! between subroutines because array in types do not get their
!   ! index boundaries reset. In the bottom line, this avoids the need to
!   ! pass array boundaries as additional arguments.
!   ! The fields st_start, st_end, and nst are for convenience, the following
!   ! equalities should hold:
!   !   lbound(stb, 3)    = st_start
!   !   ubound(stb, 3)    = st_end
!   !   st_end-st_start+1 = nst
!   ! with stb of type(states_elec_block_t).
!   type states_elec_block_t
!     integer        :: st_start
!     integer        :: st_end
!     integer        :: nst
!     FLOAT, pointer :: dpsi(:, :, :) ! dpsi(gr%mesh%np_part, DIM, st_start:st_end)
!     CMPLX, pointer :: zpsi(:, :, :) ! zpsi(gr%mesh%np_part, DIM, st_start:st_end)
!   end type states_elec_block_t

  interface states_elec_blockt_mul
    module procedure dstates_elec_blockt_mul, zstates_elec_blockt_mul
  end interface states_elec_blockt_mul

  interface states_elec_block_matr_mul_add
    module procedure dstates_elec_block_matr_mul_add, zstates_elec_block_matr_mul_add
  end interface states_elec_block_matr_mul_add

  type(profile_t), save, public :: &
    C_PROFILING_BLOCKT,            &
    C_PROFILING_BLOCKT_AR,         &
    C_PROFILING_BLOCKT_MM,         &
    C_PROFILING_BLOCKT_CP,         &
    C_PROFILING_BLOCK_MATR,        &
    C_PROFILING_BLOCK_MATR_CP,     &
    C_PROFILING_BLOCK_MATR_MM
  
contains

  ! ---------------------------------------------------------
  !> From the global state number index set global_idx(1:length)
  !! create the local sets idx(1:cnt(rank), rank), for rank=0, ..., comm_size-1
  !! using the states distribution stored in st.
  subroutine states_elec_block_local_idx(st, global_idx, length, cnt, idx)
    type(states_elec_t), intent(in) :: st
    integer,             intent(in) :: global_idx(:)
    integer,             intent(in) :: length
    integer,             pointer    :: cnt(:)
    integer,             pointer    :: idx(:, :)

    integer :: size, node, ist, i

    PUSH_SUB(states_elec_block_local_idx)

    size = st%mpi_grp%size
    SAFE_ALLOCATE(cnt(0:size-1))

    ! Count the how many vectors each node has.
    cnt = 0
    do i = 1, length
      cnt(st%node(global_idx(i))) = cnt(st%node(global_idx(i))) + 1
    end do
    ! Allocate space, it is a bit more than really required but makes the code simpler.
    SAFE_ALLOCATE(idx(1:maxval(cnt), 0:size-1))

    ! Now set up the index sets.
    cnt = 0
    idx = 0
    do ist = 1, st%nst
      node = st%node(ist)
      ! A state ist is only included if its in the global index set.
      if(member(ist, global_idx)) then
        cnt(node)            = cnt(node) + 1
        idx(cnt(node), node) = ist
      end if
    end do

    POP_SUB(states_elec_block_local_idx)
  end subroutine states_elec_block_local_idx


#include "undef.F90"
#include "real.F90"
#include "states_elec_block_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_elec_block_inc.F90"
#include "undef.F90"
end module states_elec_block_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
