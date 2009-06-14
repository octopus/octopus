! Copyright (C) 2009 X. Andrade
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
!! $Id: scissor_inc.F90 4866 2009-01-14 23:21:40Z xavier $

subroutine X(scissor_init)(this,st, gap, psi)
  type(scissor_t),        intent(inout) :: this
  type(states_t), target, intent(in)    :: st
  FLOAT,                  intent(in)    :: gap
  R_TYPE, target,         intent(in)    :: psi(:, :, :, :)

  ASSERT(.not. this%apply)

  this%apply = .true.
  this%gap = gap
  this%X(psi) => psi
  this%st => st
end subroutine X(scissor_init)

subroutine X(scissor_apply)(this, mesh, ik, psi, spsi)
  type(scissor_t), intent(in)    :: this
  type(mesh_t),    intent(in)    :: mesh
  integer,         intent(in)    :: ik
  R_TYPE,          intent(in)    :: psi(:, :)
  R_TYPE,          intent(inout) :: spsi(:, :)

  integer :: ist, idim
  R_TYPE  :: dot

  ASSERT(this%apply)
  ASSERT(associated(this%X(psi)))
  ASSERT(associated(this%st))

  do idim = 1, this%st%d%dim
    call lalg_axpy(mesh%np, this%gap, psi(:, idim), spsi(:, idim))
  end do
  do ist = 1, this%st%nst
    if(this%st%occ(ist, ik) < CNST(0.001)) cycle
    dot = X(mf_dotp)(mesh, this%st%d%dim, this%X(psi)(:, :, ist, ik), psi)
    do idim = 1, this%st%d%dim
      call lalg_axpy(mesh%np, -this%gap*dot, this%X(psi)(:, idim, ist, ik), spsi(:, idim))
    end do
  end do
  
end subroutine X(scissor_apply)
