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
!! $Id: epot.F90 2648 2007-01-09 19:08:10Z lorenzen $


!------------------------------------------------------------------------------
! X(kb_project) calculates the action of the projector kb_p on the psi 
! wavefunction. The action of the projector kb_p is defined as:
! \hat{kb_p} |psi> = \sum_{i}^kb_p%n_c p%e(i) |kb_p%p(:, i)><kb_p%p(:, i)|psi>
! The result is summed up to ppsi.
!------------------------------------------------------------------------------
subroutine X(kb_project)(mesh, sm, kb_p, dim, psi, ppsi)
  type(mesh_t),         intent(in)  :: mesh
  type(submesh_t),      intent(in)  :: sm
  type(kb_projector_t), intent(in)  :: kb_p
  integer,              intent(in)  :: dim
  R_TYPE,               intent(in)  :: psi(:, :)  ! psi(kb%n_s, dim)
  R_TYPE,               intent(out) :: ppsi(:, :) ! ppsi(kb%n_s, dim)

  integer :: i, idim, n_s
  R_TYPE :: uvpsi

  call push_sub('kb_projector_inc.kb_project')

  n_s = kb_p%n_s
  ppsi = R_TOTYPE(M_ZERO)

  do idim = 1, dim
    do i = 1, kb_p%n_c
      if (kb_p%e(i) == M_ZERO) cycle

      uvpsi = X(dsm_integrate_prod)(mesh, sm, psi(1:n_s, idim), kb_p%p(1:n_s, i))

#ifdef R_TCOMPLEX
      ppsi(1:n_s, idim) = ppsi(1:n_s, idim) + kb_p%e(i)*uvpsi*kb_p%p(1:n_s, i)
#else
      call lalg_axpy(n_s, kb_p%e(i)*uvpsi, kb_p%p(:, i), ppsi(:, idim))
#endif
      
    end do
  end do

  call pop_sub()
end subroutine X(kb_project)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
