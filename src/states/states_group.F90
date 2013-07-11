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
!! $Id: states.F90 10697 2013-05-23 10:56:09Z acastro $

#include "global.h"

module states_group_m
  use batch_m
  use batch_ops_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::                           &
    states_group_t,                   &
    states_group_null

  type states_group_t
    type(batch_t), pointer   :: psib(:, :)            !< A set of wave-functions blocks
    integer                  :: nblocks               !< The number of blocks
    integer                  :: block_start           !< The lowest index of local blocks
    integer                  :: block_end             !< The highest index of local blocks
    integer, pointer         :: iblock(:, :)          !< A map, that for each state index, returns the index of block containing it
    integer, pointer         :: block_range(:, :)     !< Each block contains states from block_range(:, 1) to block_range(:, 2)
    integer, pointer         :: block_size(:)         !< The number of states in each block.
    logical, pointer         :: block_is_local(:, :)  !< It is true if the block is in this node.
    logical                  :: block_initialized     !< For keeping track of the blocks to avoid memory leaks
  end type states_group_t

contains

  ! ---------------------------------------------------------
  subroutine states_group_null(this)
    type(states_group_t), intent(inout) :: this

    PUSH_SUB(states_group_null)

    nullify(this%psib)
    nullify(this%iblock)
    nullify(this%block_range)
    nullify(this%block_size)
    nullify(this%block_is_local)

    this%block_initialized = .false.

    POP_SUB(states_group_null)
  end subroutine states_group_null

end module states_group_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
