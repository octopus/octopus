!! Copyright (C) 2009 X. Andrade
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
!! $Id: em_field_inc.F90 3988 2008-03-31 15:06:50Z fnog $

subroutine X(em_field_apply_batch)(this, mesh, psib, vpsib)
  type(em_field_t),    intent(in)    :: this
  type(mesh_t),        intent(in)    :: mesh
  type(batch_t),       intent(in)    :: psib
  type(batch_t),       intent(inout) :: vpsib

  integer :: ist, idim, ipmax, ip, ip2, bs

  if(associated(this%potential)) then

    bs = hardware%dblock_size

    do ip = 1, mesh%np, bs
      ipmax = min(mesh%np, ip + bs - 1)

      forall (ist = 1:psib%nst, idim = 1:psib%dim, ip2 = ip:ipmax)
        vpsib%states(ist)%X(psi)(ip2, idim) = this%potential(ip2)*psib%states(ist)%X(psi)(ip2, idim)
      end forall
      
    end do

    call profiling_count_operations((R_ADD + R_MUL*psib%nst)*mesh%np)
    call profiling_count_transfers(mesh%np, M_ONE)
    call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))
    
  end if

  ASSERT(.not. associated(this%vector_potential))
  ASSERT(.not. associated(this%uniform_magnetic_field))

end subroutine X(em_field_apply_batch)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
