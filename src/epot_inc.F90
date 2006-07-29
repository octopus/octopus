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
!! $Id$

!------------------------------------------------------------------------------
! X(project) calculates the action of the sum of the np projectors p(1:np) on
! the psi wavefunction. The result is summed up to ppsi:
! |ppsi> = |ppsi> + \sum_{ip=1}^{np} \hat{p}(ip) |psi>
! The action of the projector p is defined as:
! \hat{p} |psi> = \sum_{ij} p%uvu(i,j) |p%ket(:, i)><p%bra(:, j)|psi>
!------------------------------------------------------------------------------
subroutine X(project)(mesh, p, n_projectors, psi, ppsi, periodic, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: p(:)
  integer,           intent(in)    :: n_projectors
  R_TYPE,            intent(in)    :: psi(:)   ! psi(1:mesh%np)
  R_TYPE,            intent(inout) :: ppsi(:)  ! ppsi(1:mesh%np)
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik

  integer :: i, j, n_s, ip, k
  R_TYPE, allocatable :: lpsi(:), plpsi(:)
  R_TYPE :: uvpsi
#if defined(HAVE_MPI)
  R_TYPE :: tmp
#endif

  call push_sub('epot_inc.project')

  ! index labels the atom
  k = p(1)%iatom - 1 ! This way I make sure that k is not equal to p(1)%iatom

  do ip = 1, n_projectors

    if(p(ip)%iatom .ne. k) then
      if(ip.ne.1) ppsi(p(ip-1)%jxyz(1:n_s)) = ppsi(p(ip-1)%jxyz(1:n_s)) + plpsi(1:n_s)
      n_s = p(ip)%n_points_in_sphere

      deallocate(lpsi, plpsi, stat = j)
      ALLOCATE( lpsi(n_s), n_s)
      ALLOCATE(plpsi(n_s), n_s)

      lpsi(1:n_s)  = psi(p(ip)%jxyz(1:n_s))*mesh%vol_pp(p(ip)%jxyz(1:n_s))
      if(periodic) lpsi(1:n_s)  = lpsi(1:n_s) * p(ip)%phases(1:n_s, ik)
      plpsi(1:n_s) = R_TOTYPE(M_ZERO)

      k = p(ip)%iatom
    end if

    do j = 1, p(ip)%n_channels
      uvpsi = sum(lpsi(1:n_s)*p(ip)%bra(1:n_s, j))
#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(uvpsi, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
        uvpsi = tmp
      end if
#endif
      do i = 1, p(ip)%n_channels
        if(periodic) then
          plpsi(1:n_s) = plpsi(1:n_s) + &
            p(ip)%uvu(i, j) * uvpsi * p(ip)%ket(1:n_s, i) * R_CONJ(p(ip)%phases(1:n_s, ik))
        else
          plpsi(1:n_s) = plpsi(1:n_s) + p(ip)%uvu(i, j) * uvpsi * p(ip)%ket(1:n_s, i)
        end if
      end do
    end do

  end do

  ppsi(p(n_projectors)%jxyz(1:n_s)) = ppsi(p(n_projectors)%jxyz(1:n_s)) + plpsi(1:n_s)

  deallocate(plpsi, lpsi)
  call pop_sub()
end subroutine X(project)


!------------------------------------------------------------------------------
! X(psiprojectpsi) calculates the expectation of the psi wavefuction taken over
! the sum of the projectors p(1:np).
! The action of the projector p is defined as above:
! \hat{p} |psi> = \sum_{ij} p%uvu(i,j) |p%ket(:, i)><p%bra(:, j)|psi>.
!------------------------------------------------------------------------------
R_TYPE function X(psiprojectpsi)(mesh, p, n_projectors, psi, periodic, ik) result(res)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: p(:)
  integer,           intent(in)    :: n_projectors
  R_TYPE,            intent(in)    :: psi(:)   ! psi(1:mesh%np)
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik

  integer :: i, j, n_s, ip, k
  R_TYPE, allocatable :: lpsi(:), plpsi(:)
  R_TYPE :: uvpsi
#if defined(HAVE_MPI)
  R_TYPE :: tmp
#endif

  call push_sub('epot_inc.project')

  res = R_TOTYPE(M_ZERO)

  ! index labels the atom
  k = p(1)%iatom - 1 ! This way I make sure that k is not equal to p(1)%iatom

  do ip = 1, n_projectors

    if(p(ip)%iatom .ne. k) then
      n_s = p(ip)%n_points_in_sphere
      deallocate(lpsi, plpsi, stat = j)
      ALLOCATE( lpsi(n_s), n_s)
      ALLOCATE(plpsi(n_s), n_s)

      lpsi(1:n_s)  = psi(p(ip)%jxyz(1:n_s))*mesh%vol_pp(p(ip)%jxyz(1:n_s))
      if(periodic) lpsi(1:n_s)  = lpsi(1:n_s) * p(ip)%phases(1:n_s, ik)

      k = p(ip)%iatom
    end if

    do j = 1, p(ip)%n_channels
      uvpsi = sum(lpsi(1:n_s)*p(ip)%bra(1:n_s, j))
#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(uvpsi, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
        uvpsi = tmp
      end if
#endif
      do i = 1, p(ip)%n_channels
        if(periodic) then
          plpsi(1:n_s) = p(ip)%uvu(i, j) * uvpsi * p(ip)%ket(1:n_s, i) * R_CONJ(p(ip)%phases(1:n_s, ik))
        else
          plpsi(1:n_s) = p(ip)%uvu(i, j) * uvpsi * p(ip)%ket(1:n_s, i)
        end if
        res = res + sum(R_CONJ(lpsi(1:n_s)) * plpsi(1:n_s))
      end do
    end do

  end do

#if defined(HAVE_MPI)
  if(mesh%parallel_in_domains) then
    call MPI_Allreduce(res, tmp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    res = tmp
  end if
#endif

  deallocate(plpsi, lpsi)
  call pop_sub()
end function X(psiprojectpsi)
