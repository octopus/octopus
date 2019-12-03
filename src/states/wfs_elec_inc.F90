!! Copyright (C) 2019 M. Oliveira
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

!--------------------------------------------------------------
subroutine X(wfs_elec_init_contiguous)(this, dim, st_start, st_end, psi, ik)
  type(wfs_elec_t), intent(out)   :: this
  integer,          intent(in)    :: dim
  integer,          intent(in)    :: st_start
  integer,          intent(in)    :: st_end
  integer,          intent(in)    :: ik
  R_TYPE,   target, intent(in)    :: psi(:, :, st_start:)

  integer :: ist

  PUSH_SUB(X(wfs_elec_init_contiguous))

  this%ik = ik
  this%has_phase = .false.
  call batch_init(this%batch_t, dim,  st_start, st_end, psi)

  POP_SUB(X(wfs_elec_init_contiguous))
end subroutine X(wfs_elec_init_contiguous)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
