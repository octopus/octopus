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

subroutine X(batch_init_contiguous)(this, st_start, st_end, ik, psi)
  type(batch_t),  intent(out)   :: this
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  integer,        intent(in)    :: ik
  R_TYPE, target, intent(in)    :: psi(:, :, st_start:)

  integer :: place, ist

  call batch_init_empty(this, st_end - st_start + 1)

  place = 1
  do ist = st_start, st_end
    call X(batch_add_state)(this, place, ist, ik, psi(:, :, ist))
    place = place + 1
  end do

end subroutine X(batch_init_contiguous)

subroutine X(batch_add_state)(this, place, ist, ik, psi)
  type(batch_t),  intent(out)   :: this
  integer,        intent(in)    :: place
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: ik
  R_TYPE, target, intent(in)    :: psi(:, :)

  ASSERT(place <= this%nst)

  this%states(place)%ist = ist
  this%states(place)%ik  = ik
  this%states(place)%X(psi) => psi

end subroutine X(batch_add_state)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
