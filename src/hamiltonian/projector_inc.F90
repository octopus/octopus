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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

!------------------------------------------------------------------------------
!> X(projector_matrix_element) calculates <psia|projector|psib>
R_TYPE function X(projector_matrix_element)(pj, dim, ik, psia, psib) result(apb)
  type(projector_t), target, intent(in)    :: pj
  integer,                   intent(in)    :: dim
  integer,                   intent(in)    :: ik
  R_TYPE,                    intent(in)    :: psia(:, :)  !< psia(1:mesh%np, dim)
  R_TYPE,                    intent(in)    :: psib(:, :)  !< psib(1:mesh%np, dim)

  integer ::  ns, idim, is
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)
  type(mesh_t), pointer :: mesh

  PUSH_SUB(X(projector_matrix_element))

  ns = pj%sphere%np

  ASSERT(associated(pj%sphere%mesh))
  mesh => pj%sphere%mesh

  SAFE_ALLOCATE(lpsi(1:ns, 1:dim))
  SAFE_ALLOCATE(plpsi(1:ns, 1:dim))

  do idim = 1, dim
    if(associated(pj%phase)) then
      do is = 1, ns
        lpsi(is, idim) = psib(pj%sphere%map(is), idim)*pj%phase(is, ik)
      end do
    else
      do is = 1, ns
        lpsi(is, idim) = psib(pj%sphere%map(is), idim)
      end do
    end if
  end do

  call X(project_sphere)(mesh, pj, dim, lpsi, plpsi)

  apb = M_ZERO
  do idim = 1, dim
    if(associated(pj%phase)) then
      do is = 1, ns
        plpsi(is, idim) = R_CONJ(psia(pj%sphere%map(is), idim))*plpsi(is, idim)*conjg(pj%phase(is, ik))
      end do
    else
      do is = 1, ns
        plpsi(is, idim) = R_CONJ(psia(pj%sphere%map(is), idim))*plpsi(is, idim)
      end do
    end if

    if(ns > 0) then
      apb = apb + X(sm_integrate)(mesh, pj%sphere, plpsi(1:ns, idim))
    else
      apb = apb + X(sm_integrate)(mesh, pj%sphere)
    end if
  end do

  SAFE_DEALLOCATE_A(lpsi)
  SAFE_DEALLOCATE_A(plpsi)
  
  POP_SUB(X(projector_matrix_element))
end function X(projector_matrix_element)

!------------------------------------------------------------------------------

subroutine X(project_sphere)(mesh, pj, dim, psi, ppsi)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:ns, dim)
  R_TYPE,            intent(out)   :: ppsi(:, :)  ! ppsi(1:ns, dim)

  integer :: ll, mm

  PUSH_SUB(X(project_sphere))

  ppsi = M_ZERO

  do ll = 0, pj%lmax
    if (ll == pj%lloc) cycle
    do mm = -ll, ll
      call X(kb_project)(mesh, pj%sphere, pj%kb_p(ll, mm), dim, psi, ppsi)
    end do
  end do

  POP_SUB(X(project_sphere))
end subroutine X(project_sphere)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
