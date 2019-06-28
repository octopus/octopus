!! Copyright (C) 2002-20069 M. Marques, X. Andrade
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

module transfer_table_oct_m

  implicit none

  private
  public ::                         &
    transfer_table_t,               &
    transfer_table_nullify


  type transfer_table_t
    ! Components are public by default
    integer          ::  n_coarse
    integer, pointer :: to_coarse(:)

    integer          ::  n_fine
    integer          ::  n_fine1,       n_fine2,       n_fine4,       n_fine8
    integer, pointer :: to_fine1(:,:), to_fine2(:,:), to_fine4(:,:), to_fine8(:,:)
    integer, pointer :: fine_i(:)
  end type transfer_table_t

contains

  ! ---------------------------------------------------------
  elemental subroutine transfer_table_nullify(this)
    type(transfer_table_t), intent(out) :: this

    this%n_coarse = 0
    nullify(this%to_coarse)
    this%n_fine = 0
    this%n_fine1 = 0
    this%n_fine2 = 0
    this%n_fine4 = 0
    this%n_fine8 = 0
    nullify(this%to_fine1, this%to_fine2, this%to_fine4, this%to_fine8)
    nullify(this%fine_i)

  end subroutine transfer_table_nullify

end module transfer_table_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
