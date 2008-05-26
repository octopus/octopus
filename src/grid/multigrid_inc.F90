!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

  ! ---------------------------------------------------------
  subroutine X(multigrid_coarse2fine)(level, f_coarse, f_fine)
    type(multigrid_level_t), intent(in)  :: level
    R_TYPE,                  intent(in)  :: f_coarse(:)
    R_TYPE,                  intent(out) :: f_fine(:)

    FLOAT, pointer :: vol_pp(:)

    integer :: i, i1, i2, i4, i8
    integer :: jj, j(8)
    FLOAT   :: vol_total

    call push_sub('multigrid.multigrid_coarse2fine')

    vol_pp => level%m%vol_pp

    i1 = 0;  i2 = 0;  i4 = 0;  i8 = 0;
    do i = 1, level%n_fine
      select case(level%fine_i(i))
      case(1)
        i1 = i1 + 1
        j(1:1) = level%to_fine1(1:1, i1)
      case(2)
        i2 = i2 + 1
        j(1:2) = level%to_fine2(1:2, i2)
      case(4)
        i4 = i4 + 1
        j(1:4) = level%to_fine4(1:4, i4)
      case(8)
        i8 = i8 + 1
        j(1:8) = level%to_fine8(1:8, i8)
      end select

      f_fine(i) = M_ZERO
      vol_total = M_ZERO
      do jj = 1, level%fine_i(i)
        f_fine(i) = f_fine(i) + vol_pp(j(jj))*f_coarse(j(jj))
        vol_total = vol_total + vol_pp(j(jj))
      end do
      f_fine(i) = f_fine(i)/vol_total

    end do

    call pop_sub()
  end subroutine X(multigrid_coarse2fine)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
