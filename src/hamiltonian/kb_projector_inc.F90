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
!> X(kb_project) calculates the action of the projector kb_p on the psi 
!! wavefunction. The action of the projector kb_p is defined as:
!! \f[
!! \hat{kb_p} |psi> = \sum_{i}^kb_p%n_c p%e(i) |kb_p%p(:, i)><kb_p%p(:, i)|psi>
!! \f]
!! The result is summed up to ppsi.
!------------------------------------------------------------------------------
subroutine X(kb_project)(mesh, sm, kb_p, dim, psi, ppsi)
  type(mesh_t),         intent(in)    :: mesh
  type(submesh_t),      intent(in)    :: sm
  type(kb_projector_t), intent(in)    :: kb_p
  integer,              intent(in)    :: dim
  R_TYPE,               intent(in)    :: psi(:, :)  !< psi(kb%n_s, dim)
  R_TYPE,               intent(inout) :: ppsi(:, :) !< ppsi(kb%n_s, dim)

  R_TYPE :: uvpsi(1:2)
#if defined(HAVE_MPI)
  R_TYPE :: uvpsi_tmp(1:2)
#endif

  PUSH_SUB(X(kb_project))

  call X(kb_project_bra)(mesh, sm, kb_p, dim, psi, uvpsi)

#if defined(HAVE_MPI)
  if(mesh%parallel_in_domains) then
    call MPI_Allreduce(uvpsi, uvpsi_tmp, 2, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    uvpsi = uvpsi_tmp
  end if
#endif

  call X(kb_project_ket)(mesh, sm, kb_p, dim, uvpsi, ppsi)

  POP_SUB(X(kb_project))

end subroutine X(kb_project)

!--------------------------------------------------------------
!> THREADSAFE
subroutine X(kb_project_bra)(mesh, sm, kb_p, dim, psi, uvpsi)
  type(mesh_t),         intent(in)  :: mesh
  type(submesh_t),      intent(in)  :: sm
  type(kb_projector_t), intent(in)  :: kb_p
  integer,              intent(in)  :: dim
  R_TYPE,               intent(in)  :: psi(:, :)  ! psi(1:ns, 1:dim)
  R_TYPE,               intent(out) :: uvpsi(1:kb_p%n_c, 1:dim)

  integer :: ic, idim, ns, is

#ifndef HAVE_OPENMP
  PUSH_SUB(X(kb_project_bra))
#endif

  ns = kb_p%n_s

  uvpsi(1:kb_p%n_c, 1:dim) = M_ZERO

  if(mesh%use_curvilinear) then

    do idim = 1, dim
      do ic = 1, kb_p%n_c
        do is = 1, ns
          uvpsi(ic, idim) = uvpsi(ic, idim) + (kb_p%p(is, ic)*psi(is, idim))*mesh%vol_pp(sm%map(is))
        end do
      end do
    end do

  else

    do idim = 1, dim
      do ic = 1, kb_p%n_c
        do is = 1, ns
          uvpsi(ic, idim) = uvpsi(ic, idim) + psi(is, idim)*kb_p%p(is, ic)
        end do
      end do
    end do

    uvpsi(1:kb_p%n_c, 1:dim) = uvpsi(1:kb_p%n_c, 1:dim)*mesh%vol_pp(1)
  end if

#ifndef HAVE_OPENMP
  POP_SUB(X(kb_project_bra))
#endif
end subroutine X(kb_project_bra)

!--------------------------------------------------------------
!> THREADSAFE
subroutine X(kb_project_ket)(mesh, sm, kb_p, dim, uvpsi, psi)
  type(mesh_t),         intent(in)    :: mesh
  type(submesh_t),      intent(in)    :: sm
  type(kb_projector_t), intent(in)    :: kb_p
  integer,              intent(in)    :: dim
  R_TYPE,               intent(inout) :: uvpsi(1:kb_p%n_c, 1:dim)
  R_TYPE,               intent(inout) :: psi(:, :) ! psi(1:ns, 1:dim)

  integer :: ic, idim, ns, is
#ifndef HAVE_OPENMP
  PUSH_SUB(X(kb_project_ket))
#endif

  ns = kb_p%n_s

  call X(kb_mul_energies)(kb_p, dim, uvpsi)

  do idim = 1, dim
    do ic = 1, kb_p%n_c
      do is = 1, ns
        psi(is, idim) = psi(is, idim) + uvpsi(ic, idim)*kb_p%p(is, ic)
      end do
    end do
  end do

#ifndef HAVE_OPENMP
  POP_SUB(X(kb_project_ket))
#endif
end subroutine X(kb_project_ket)

!--------------------------------------------------------------
!> THREADSAFE
subroutine X(kb_mul_energies)(kb_p, dim, uvpsi)
  type(kb_projector_t), intent(in)    :: kb_p
  integer,              intent(in)    :: dim
  R_TYPE,               intent(inout) :: uvpsi(1:kb_p%n_c, 1:dim)
  
  integer :: idim

  do idim = 1, dim
    uvpsi(1:kb_p%n_c, idim) = uvpsi(1:kb_p%n_c, idim)*kb_p%e(1:kb_p%n_c)
  end do
end subroutine X(kb_mul_energies)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
