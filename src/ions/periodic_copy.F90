!! Copyright (C) 2007 X. Andrade, M. Marques
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

module periodic_copy_oct_m
  use global_oct_m
  use io_oct_m
  use lattice_vectors_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private

  public ::                   &
    periodic_copy_t,          &
    periodic_copy_init,       &
    periodic_copy_end,        &
    periodic_copy_position,   &
    periodic_copy_num

  type periodic_copy_t
    private
    integer :: num
    FLOAT :: pos(1:MAX_DIM)
    FLOAT :: pos_chi(1:MAX_DIM)
    FLOAT :: range
    integer :: nbmax(1:MAX_DIM), nbmin(1:MAX_DIM)
    integer, allocatable :: icell(:, :) !< (dim, num)
  end type periodic_copy_t

contains

  subroutine periodic_copy_init(this, space, latt, lsize, pos, range)
    type(periodic_copy_t),    intent(out) :: this
    type(space_t),            intent(in)  :: space
    type(lattice_vectors_t),  intent(in)  :: latt
    FLOAT,                 intent(in)  :: lsize(:)
    FLOAT,                 intent(in)  :: pos(:) !< (dim)
    FLOAT,                 intent(in)  :: range

    integer :: pd, jj, kk, idir

    PUSH_SUB(periodic_copy_init)

    ASSERT(range >= M_ZERO)

    this%range = range
    this%pos(1:space%dim) = pos(1:space%dim)

    if(.not. space%is_periodic()) then
      this%num = 1
      this%nbmin = 0
      this%nbmax = 0

      POP_SUB(periodic_copy_init)
      return
    end if

    pd = space%periodic_dim

    !convert the position to the orthogonal space
    this%pos_chi(1:pd) = matmul(pos(1:pd), latt%klattice_primitive(1:pd, 1:pd))

    this%nbmin(1:pd) = -nint(-(this%pos_chi(1:pd) - range)/(M_TWO*lsize(1:pd)) + M_HALF)
    this%nbmax(1:pd) = nint((this%pos_chi(1:pd) + range)/(M_TWO*lsize(1:pd)) + M_HALF)
    ! no copies in non-periodic directions

    this%num = product(this%nbmax(1:space%periodic_dim) - this%nbmin(1:space%periodic_dim) + 1)
    SAFE_ALLOCATE(this%icell(1:space%periodic_dim, 1:this%num))

    do jj = 1, this%num
      kk = jj - 1
      do idir = space%periodic_dim, 1, -1
        this%icell(idir, jj) = mod(kk, this%nbmax(idir) - this%nbmin(idir) + 1) + this%nbmin(idir)
        if(idir > 1) &
          kk = kk / (this%nbmax(idir) - this%nbmin(idir) + 1)
      end do
    end do

    POP_SUB(periodic_copy_init)
  end subroutine periodic_copy_init

  ! ----------------------------------------------------------------

  subroutine periodic_copy_end(this)
    type(periodic_copy_t), intent(inout) :: this

    PUSH_SUB(periodic_copy_end)

    SAFE_DEALLOCATE_A(this%icell)

    this%nbmin = 0
    this%nbmax = 0

    POP_SUB(periodic_copy_end)
  end subroutine periodic_copy_end

  ! ----------------------------------------------------------------

  integer pure function periodic_copy_num(this) result(num)
    type(periodic_copy_t), intent(in)    :: this

    ! no push_sub allowed in pure function
    num = this%num

  end function periodic_copy_num
  
  ! ----------------------------------------------------------------

  pure function periodic_copy_position(this, space, latt, lsize, ii) result(pcopy)
    type(periodic_copy_t),   intent(in)  :: this
    type(space_t),           intent(in)  :: space
    type(lattice_vectors_t), intent(in)  :: latt
    FLOAT,                   intent(in)  :: lsize(:)
    integer,                 intent(in)  :: ii
    FLOAT                                :: pcopy(space%dim)
    
    integer :: pd

    pd = space%periodic_dim

    if(.not. space%is_periodic()) then
      pcopy(1:space%dim) = this%pos(1:space%dim)
      return
    end if

    pcopy(1:pd) = this%pos_chi(1:pd) - M_TWO*lsize(1:pd)*this%icell(1:pd, ii)
    pcopy(1:pd) = matmul(latt%rlattice_primitive(1:pd, 1:pd), pcopy(1:pd))
    pcopy(pd + 1:space%dim) = this%pos(pd+1:space%dim)

  end function periodic_copy_position

end module periodic_copy_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
