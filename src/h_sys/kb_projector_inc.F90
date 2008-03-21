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

  R_TYPE :: uvpsi(1:2)
#if defined(HAVE_MPI)
  R_TYPE :: uvpsi_tmp(1:2)
#endif

  call push_sub('kb_projector_inc.kb_project')

  call X(kb_project_bra)(mesh, sm, kb_p, dim, psi, uvpsi)

#if defined(HAVE_MPI)
  if(mesh%parallel_in_domains) then
    call MPI_Allreduce(uvpsi, uvpsi_tmp, 2, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    uvpsi = uvpsi_tmp
  end if
#endif

  uvpsi(1:2) = kb_p%e(1:2)*uvpsi(1:2)
  call X(kb_project_ket)(mesh, sm, kb_p, dim, uvpsi, ppsi)

  call pop_sub()
end subroutine X(kb_project)

subroutine X(kb_project_bra)(mesh, sm, kb_p, dim, psi, uvpsi)
  type(mesh_t),         intent(in)  :: mesh
  type(submesh_t),      intent(in)  :: sm
  type(kb_projector_t), intent(in)  :: kb_p
  integer,              intent(in)  :: dim
  R_TYPE,               intent(in)  :: psi(:, :)  ! psi(1:ns, 1:dim)
  R_TYPE,               intent(out) :: uvpsi(1:2)

  integer :: ic, idim, ns, is
  FLOAT, allocatable :: kp(:, :)

  call push_sub('kb_projector_inc.kb_project')

  ns = kb_p%n_s

  ALLOCATE(kp(1:ns, 1:2), ns*2)

  kp = M_ZERO

  if(mesh%use_curvlinear) then
    do ic = 1, kb_p%n_c
      kp(1:ns, ic) = kb_p%p(1:ns, ic)*mesh%vol_pp(sm%jxyz(1:ns))
    end do
  else
    do ic = 1, kb_p%n_c
      kp(1:ns, ic) = kb_p%p(1:ns, ic)*mesh%vol_pp(1)
    end do
  end if

  call profiling_count_operations(ns*kb_p%n_c)

  uvpsi = M_ZERO

  do idim = 1, dim
    do is = 1, ns
      uvpsi(1) = uvpsi(1) + psi(is, idim)*kp(is, 1)
      uvpsi(2) = uvpsi(2) + psi(is, idim)*kp(is, 2)
    end do
  end do

  call profiling_count_operations(ns*dim*kb_p%n_c*(R_ADD + R_MUL))

  deallocate(kp)

  call pop_sub()
end subroutine X(kb_project_bra)

subroutine X(kb_project_ket)(mesh, sm, kb_p, dim, uvpsi, psi)
  type(mesh_t),         intent(in)    :: mesh
  type(submesh_t),      intent(in)    :: sm
  type(kb_projector_t), intent(in)    :: kb_p
  integer,              intent(in)    :: dim
  R_TYPE,               intent(in)    :: uvpsi(1:2)
  R_TYPE,               intent(out)   :: psi(:, :) ! psi(1:ns, 1:dim)

  integer :: ic, idim, ns, is

  call push_sub('kb_projector_inc.kb_project')

  psi = M_ZERO

  ns = kb_p%n_s

  do idim = 1, dim
    do is = 1, ns
      do ic = 1, kb_p%n_c
        psi(is, idim) = psi(is, idim) + uvpsi(ic)*kb_p%p(is, ic)
      end do
    end do
  end do

  call profiling_count_operations(ns*dim*kb_p%n_c*2*R_ADD)

  call pop_sub()
end subroutine X(kb_project_ket)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
