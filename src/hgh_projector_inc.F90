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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: epot.F90 2648 2007-01-09 19:08:10Z lorenzen $


!------------------------------------------------------------------------------
! X(hgh_project) calculates the action of the projector hgh_p on the psi 
! wavefunction. The action of the projector hgh_p is defined as:
! \hat{hgh_p} |psi> = \sum_{ij}^3 p%h(i,j) |hgh_p%p(:, i)><hgh_p%p(:, j)|psi>
! The result is summed up to ppsi.
!
! If including the spin-orbit coupling there is another term to be added to ppsi:
! \sum_{ij}^3\sum{k}^3 p%k(i,j) |hgh_p%p(:, i)><hgh_p%lp(:, k, j)|\hat{S(k)}|psi>
!------------------------------------------------------------------------------
subroutine X(hgh_project)(mesh, hgh_p, dim, psi, ppsi, reltype, phases)
  type(mesh_t),          intent(in)  :: mesh
  type(hgh_projector_t), intent(in)  :: hgh_p
  integer,               intent(in)  :: dim
  R_TYPE,                intent(in)  :: psi(:, :)  ! psi(hgh%n_s, dim)
  R_TYPE,                intent(out) :: ppsi(:, :) ! ppsi(hgh%n_s, dim)
  integer,               intent(in)  :: reltype
  CMPLX, optional,       intent(in)  :: phases(:)

  integer :: n_s, i, j, k, idim
  R_TYPE :: uvpsi
#ifdef R_TCOMPLEX
  CMPLX, allocatable :: lp_psi(:, :, :)
#endif
#if defined(HAVE_MPI)
  R_TYPE :: tmp
#endif

  n_s = hgh_p%n_s
  ppsi = R_TOTYPE(M_ZERO)

  do idim = 1, dim
    do j = 1, 3
      if (all(hgh_p%h(:, j) == M_ZERO)) cycle
      uvpsi = sum(psi(1:n_s, idim)*hgh_p%p(1:n_s, j))

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(uvpsi, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
        uvpsi = tmp
      end if
#endif
      do i = 1, 3
        if (hgh_p%h(i, j) == M_ZERO) cycle
        if (present(phases)) then
          ppsi(1:n_s, idim) = ppsi(1:n_s, idim) + &
               hgh_p%h(i, j) * uvpsi * hgh_p%p(1:n_s, i) * R_CONJ(phases(1:n_s))
        else
          ppsi(1:n_s, idim) = ppsi(1:n_s, idim) + &
               hgh_p%h(i, j) * uvpsi * hgh_p%p(1:n_s, i)
        end if
      end do
    end do
  end do

  if (reltype == 1) then
#ifdef R_TCOMPLEX
    ! Add spin-orbit term
    ASSERT(dim == 2)

    ALLOCATE(lp_psi(n_s, 3, dim), n_s*3*dim)
    lp_psi = M_Z0

    do idim = 1, dim
      do k = 1, 3
        do j = 1, 3
          if (all(hgh_p%k(:, j) == M_ZERO)) cycle
          uvpsi = sum(psi(1:n_s, idim)*hgh_p%lp(1:n_s, k, j))
#if defined(HAVE_MPI)
          if(mesh%parallel_in_domains) then
            call MPI_Allreduce(uvpsi, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
            uvpsi = tmp
          end if
#endif
          
          do i = 1, 3
            if (hgh_p%k(i, j) == M_ZERO) cycle
            if (present(phases)) then
              lp_psi(1:n_s, k, idim) = lp_psi(1:n_s, k, idim) + &
                   hgh_p%k(i, j) * uvpsi * hgh_p%p(1:n_s, i) * R_CONJ(phases(1:n_s))
            else
              lp_psi(1:n_s, k, idim) = lp_psi(1:n_s, k, idim) + &
                   hgh_p%k(i, j) * uvpsi * hgh_p%p(1:n_s, i)
            end if
          end do
        end do
      end do
    end do

    ppsi(1:n_s, 1) = ppsi(1:n_s, 1) + M_zI*M_HALF*( lp_psi(1:n_s, 3, 1) + lp_psi(1:n_s, 1, 2) - M_zI*lp_psi(1:n_s, 2, 2))
    ppsi(1:n_s, 2) = ppsi(1:n_s, 2) + M_zI*M_HALF*(-lp_psi(1:n_s, 3, 2) + lp_psi(1:n_s, 1, 1) + M_zI*lp_psi(1:n_s, 2, 1))

    deallocate(lp_psi)

#endif
  end if

end subroutine X(hgh_project)


!------------------------------------------------------------------------------
! X(hgh_dproject) calculates:
!  \sum_{ij}^3\sum{k}^3 p%h(i,j) <psi|hgh_p%dp(:, k, i)><hgh_p%p(:, j)|psi>
! Here the spin-orbit coupling will always be ignored.
!------------------------------------------------------------------------------
function X(hgh_dproject)(mesh, hgh_p, dim, psi, phases) result(res)
  type(mesh_t),          intent(in)    :: mesh
  type(hgh_projector_t), intent(in) :: hgh_p
  integer,               intent(in) :: dim
  R_TYPE,                intent(in) :: psi(:, :)  ! psi(hgh%n_s, dim)
  CMPLX, optional,       intent(in) :: phases(:)
  R_TYPE :: res(3)

  integer :: n_s, i, j, k, idim
  R_TYPE :: uvpsi
  R_TYPE, allocatable :: ppsi(:)
#if defined(HAVE_MPI)
  R_TYPE :: tmp
#endif

  res = R_TOTYPE(M_ZERO)
  n_s = hgh_p%n_s
  ALLOCATE(ppsi(n_s), n_s)

  do idim = 1, dim
    do j = 1, 3
      if (all(hgh_p%h(:, j) == M_ZERO)) cycle
      uvpsi = sum(psi(1:n_s, idim)*hgh_p%p(1:n_s, j))
#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(uvpsi, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
        uvpsi = tmp
      end if
#endif
      do i = 1, 3
        if (hgh_p%h(i, j) == M_ZERO) cycle
        do k = 1, 3
          if (present(phases)) then
            ppsi(1:n_s) = hgh_p%h(i, j) * uvpsi * hgh_p%dp(1:n_s, k, i) * R_CONJ(phases(1:n_s))
          else
            ppsi(1:n_s) = hgh_p%h(i, j) * uvpsi * hgh_p%dp(1:n_s, k, i)
          end if
          res(k) = res(k) + sum(R_CONJ(psi(1:n_s, idim)) * ppsi(1:n_s))
        end do
      end do
    end do
  end do

  deallocate(ppsi)

#if defined(HAVE_MPI)
  if(mesh%parallel_in_domains) then
    call MPI_Allreduce(res, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    res = tmp
  end if
#endif

end function X(hgh_dproject)
