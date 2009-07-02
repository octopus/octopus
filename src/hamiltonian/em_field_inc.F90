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

  integer :: ist, idim, ipmax, ip

  if(associated(this%potential)) then
    !$omp parallel do private(ipmax, idim, ip)
    do ist = 1, psib%nst
      forall (ip=1:mesh%np, idim = 1:psib%dim)
        vpsib%states(ist)%X(psi)(ip, idim) = this%potential(ip)*psib%states(ist)%X(psi)(ip, idim)
      end forall
    end do
    !$omp end parallel do

    call profiling_count_operations((R_MUL*psib%nst)*mesh%np)
    call profiling_count_transfers(mesh%np, M_ONE)
    call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))
    
  end if

  ASSERT(.not. associated(this%vector_potential))
  ASSERT(.not. associated(this%uniform_magnetic_field))

end subroutine X(em_field_apply_batch)

subroutine X(em_field_apply)(this, std, mesh, time, psi, vpsi)
  type(em_field_t),    intent(in)    :: this
  type(states_dim_t),  intent(in)    :: std
  type(mesh_t),        intent(in)    :: mesh
  FLOAT,               intent(in)    :: time
  R_TYPE, target,      intent(inout) :: psi(:,:)
  R_TYPE,              intent(out)   :: vpsi(:,:)

  type(batch_t) :: psib, vpsib

  call batch_init(psib, std%dim, 1)
  call batch_add_state(psib, 1, psi)
  call batch_init(vpsib, std%dim, 1)
  call batch_add_state(vpsib, 1, vpsi)

  call X(em_field_apply_batch)(this, mesh, psib, vpsib)

  call batch_end(psib)
  call batch_end(vpsib)
end subroutine X(em_field_apply)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
