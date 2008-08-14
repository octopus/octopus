!! Copyright (C) 2008 X. Andrade
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
!! $Id: batch_inc.F90 4298 2008-06-18 15:03:00Z dstrubbe $

subroutine X(batch_init_contiguous)(this, st_start, st_end, psi)
  type(batch_t),  intent(out)   :: this
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  R_TYPE, target, intent(in)    :: psi(:, :, st_start:)

  integer :: ist

  call batch_init_empty(this, st_end - st_start + 1)

  do ist = st_start, st_end
    call X(batch_add_state)(this, ist, psi(:, :, ist))
  end do

end subroutine X(batch_init_contiguous)

subroutine X(batch_add_state)(this, ist, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: ist
  R_TYPE, target, intent(in)    :: psi(:, :)

  ASSERT(this%current <= this%nst)

  this%states(this%current)%ist = ist
  this%states(this%current)%X(psi) => psi

  this%current = this%current  + 1

end subroutine X(batch_add_state)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
