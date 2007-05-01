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
!! $Id: opt_control.F90 2875 2007-04-30 16:54:15Z acastro $

#include "global.h"


  ! ---------------------------------------------------------
  subroutine oct_iterator_init(iterator, oct)
    type(oct_iterator_t), intent(inout) :: iterator
    type(oct_t), intent(in)             :: oct

    call push_sub('opt_control_iter.oct_iter_init')

    iterator%old_functional = -CNST(1e10)
    iterator%ctr_iter       = 0

    ALLOCATE(iterator%convergence(4,0:oct%ctr_iter_max),(oct%ctr_iter_max+1)*4)

    iterator%bestJ           = M_ZERO
    iterator%bestJ1          = M_ZERO
    iterator%bestJ_ctr_iter  = M_ZERO
    iterator%bestJ1_ctr_iter = M_ZERO

    call pop_sub()
  end subroutine oct_iterator_init


  ! ---------------------------------------------------------
  subroutine oct_iterator_end(iterator)
    type(oct_iterator_t), intent(inout) :: iterator

    call push_sub('opt_control_iter.oct_iter_end')

    deallocate(iterator%convergence)
    nullify(iterator%convergence)

    call pop_sub()
  end subroutine oct_iterator_end


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
