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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: simul_box.F90 3479 2007-11-09 08:36:10Z xavier $

#include "global.h"

module solids_m
  use global_m
  use simul_box_m

  implicit none

  private

  public ::                   &
    periodic_copy_t, &
    periodic_copy_init,       &
    periodic_copy_end,        &
    periodic_copy_position,   &
    periodic_copy_num
  
  type :: periodic_copy_t
    private
    FLOAT :: pos(1:MAX_DIM)
    FLOAT :: range
    integer :: nbmax(MAX_DIM), nbmin(MAX_DIM)
  end type periodic_copy_t

contains

  subroutine periodic_copy_init(this, sb, pos, range)
    type(periodic_copy_t), intent(out)   :: this
    type(simul_box_t),       intent(in)    :: sb
    FLOAT,                   intent(in)    :: pos(1:MAX_DIM)
    FLOAT,                   intent(in)    :: range
    
    this%range = range

    !convert the position to the orthogonal space
    this%pos = matmul(pos, sb%klattice_unitary)

    this%nbmin = int((this%pos - range)/sb%lsize - 0.5)
    this%nbmax = int((this%pos + range)/sb%lsize + 0.5)

    !there are no copies in non-periodic directions
    this%nbmin(sb%periodic_dim + 1:sb%dim) = 0
    this%nbmax(sb%periodic_dim + 1:sb%dim) = 0

  end subroutine periodic_copy_init

  subroutine periodic_copy_end(this)
    type(periodic_copy_t), intent(out) :: this

    this%nbmin = 0
    this%nbmax = 0

  end subroutine periodic_copy_end

  integer function periodic_copy_num(this) result(num)
    type(periodic_copy_t), intent(out) :: this
    
    num = product(this%nbmax - this%nbmin + 1)
    
  end function periodic_copy_num

  function periodic_copy_position(this, sb, ii) result(pcopy)
    type(periodic_copy_t), intent(in)  :: this
    type(simul_box_t),       intent(in)  :: sb
    integer, intent(in)                  :: ii
    FLOAT                                :: pcopy(MAX_DIM) 
    
    integer :: icell1, icell2, icell3, jj
    
    jj = 1
    do icell1 = this%nbmin(1), this%nbmax(1)
      do icell2 = this%nbmin(2), this%nbmax(2)
        do icell3 = this%nbmin(3), this%nbmax(3)
          
          if (jj == ii) then 
            pcopy(:) = this%pos(:) - M_TWO*sb%lsize(:)*(/icell1, icell2, icell3/)
            pcopy(:) = matmul(sb%rlattice, pcopy(:))
            return
          end if
          jj = jj + 1

        end do
      end do
    end do

  end function periodic_copy_position

end module solids_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
