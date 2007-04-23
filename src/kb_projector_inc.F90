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
subroutine X(kb_project)(mesh, kb_p, dim, psi, ppsi, phases)
  type(mesh_t),         intent(in)  :: mesh
  type(kb_projector_t), intent(in)  :: kb_p
  integer,              intent(in)  :: dim
  R_TYPE,               intent(in)  :: psi(:, :)  ! psi(kb%n_s, dim)
  R_TYPE,               intent(out) :: ppsi(:, :) ! ppsi(kb%n_s, dim)
  CMPLX, optional,      intent(in)  :: phases(:)

  integer :: i, idim, n_s
  R_TYPE :: uvpsi
#if defined(HAVE_MPI)
  R_TYPE :: tmp
#endif

  call push_sub('kb_projector_inc.kb_project')

  n_s = kb_p%n_s
  ppsi = R_TOTYPE(M_ZERO)

  do idim = 1, dim
    do i = 1, kb_p%n_c
      if (kb_p%e(i) == M_ZERO) cycle
      uvpsi = sum(psi(1:n_s, idim)*kb_p%p(1:n_s, i))

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(uvpsi, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
        uvpsi = tmp
      end if
#endif

      if (present(phases)) then
        ppsi(1:n_s, idim) = ppsi(1:n_s, idim) + &
             kb_p%e(i) * uvpsi * kb_p%p(1:n_s, i) * R_CONJ(phases(1:n_s))
      else
        ppsi(1:n_s, idim) = ppsi(1:n_s, idim) + &
             kb_p%e(i) * uvpsi * kb_p%p(1:n_s, i)
      end if
    end do
  end do

  call pop_sub()
end subroutine X(kb_project)


!------------------------------------------------------------------------------
! X(kb_dproject) calculates:
!  \sum_{i}^kb_p%n_c\sum{k}^3 p%e(i) <psi|kb_p%dp(:, k, i)><kb_p%p(:, i)|psi>
!------------------------------------------------------------------------------
function X(kb_dproject)(mesh, kb_p, dim, psi, phases) result(res)
  type(mesh_t),         intent(in)    :: mesh
  type(kb_projector_t), intent(in) :: kb_p
  integer,              intent(in) :: dim
  R_TYPE,               intent(in) :: psi(:, :)  ! psi(kb%n_s, dim)
  CMPLX, optional,      intent(in) :: phases(:)
  R_TYPE :: res(3)

  integer :: n_s, i, k, idim
  R_TYPE :: uvpsi
  R_TYPE, allocatable :: ppsi(:)
#if defined(HAVE_MPI)
  R_TYPE :: tmp
  R_TYPE :: tmp2(3)
#endif

  call push_sub('kb_projector_inc.kb_dproject')

  res = R_TOTYPE(M_ZERO)
  n_s = kb_p%n_s
  ALLOCATE(ppsi(n_s), n_s)

  do idim = 1, dim
    do i = 1, kb_p%n_c
      if (kb_p%e(i) == M_ZERO) cycle
      uvpsi = sum(psi(1:n_s, idim)*kb_p%p(1:n_s, i))
#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(uvpsi, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
        uvpsi = tmp
      end if
#endif

      do k = 1, 3
        if (present(phases)) then
          ppsi(1:n_s) = kb_p%e(i) * uvpsi * kb_p%dp(1:n_s, k, i) * R_CONJ(phases(1:n_s))
        else
          ppsi(1:n_s) = kb_p%e(i) * uvpsi * kb_p%dp(1:n_s, k, i)
        end if
        res(k) = res(k) + sum(R_CONJ(psi(1:n_s, idim)) * ppsi(1:n_s))
      end do
    end do
  end do

  deallocate(ppsi)

#if defined(HAVE_MPI)
  if(mesh%parallel_in_domains) then
    call MPI_Allreduce(res(1), tmp2(1), 3, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    res = tmp2
  end if
#endif

  call pop_sub()
end function X(kb_dproject)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
