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
! X(project) calculates the action of the sum of the projectors p(1:n_projectors)
! on the psi wavefunction. The result is summed up to ppsi
subroutine X(project)(mesh, p, n_projectors, dim, psi, ppsi, reltype, periodic, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: p(:)
  integer,           intent(in)    :: n_projectors
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: ppsi(:, :)  ! ppsi(1:mesh%np, dim)
  integer,           intent(in)    :: reltype
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik

  integer :: n_s, k, ip, idim
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)

  call push_sub('epot_inc.project')

  k = p(1)%iatom - 1 ! This way I make sure that k is not equal to p(1)%iatom
  do ip = 1, n_projectors

    if(p(ip)%iatom .ne. k) then
      n_s = p(ip)%n_s
      if(allocated(lpsi))   deallocate(lpsi)
      if(allocated(plpsi))  deallocate(plpsi)
      ALLOCATE(lpsi(n_s, dim),  n_s*dim)
      ALLOCATE(plpsi(n_s, dim), n_s*dim)

      do idim = 1, dim
        lpsi(1:n_s, idim)  = psi(p(ip)%jxyz(1:n_s), idim)*mesh%vol_pp(p(ip)%jxyz(1:n_s))
        if(periodic) lpsi(1:n_s, idim)  = lpsi(1:n_s, idim) * p(ip)%phases(1:n_s, ik)
      end do

      k = p(ip)%iatom
    end if

    select case (p(ip)%type)
    case (M_HGH)
      if (periodic) then
        call X(hgh_project)(mesh, p(ip)%hgh_p, dim, lpsi, plpsi, reltype, p(ip)%phases(:, ik))
      else
        call X(hgh_project)(mesh, p(ip)%hgh_p, dim, lpsi, plpsi, reltype)
      end if
    case (M_KB)
      if (periodic) then
        call X(kb_project)(mesh, p(ip)%kb_p, dim, lpsi, plpsi, p(ip)%phases(:, ik))
      else
        call X(kb_project)(mesh, p(ip)%kb_p, dim, lpsi, plpsi)
      end if
    case (M_RKB)
#ifdef R_TCOMPLEX
      !This can only be aplied to complex spinor wave-functions
      if (periodic) then
        call rkb_project(mesh, p(ip)%rkb_p, lpsi, plpsi, p(ip)%phases(:, ik))
      else
        call rkb_project(mesh, p(ip)%rkb_p, lpsi, plpsi)
      end if
#endif
    end select

    do idim = 1, dim
      ppsi(p(ip)%jxyz(1:n_s), idim) = ppsi(p(ip)%jxyz(1:n_s), idim) + plpsi(1:n_s, idim)
    end do

  end do

  deallocate(plpsi, lpsi)

  call pop_sub()
end subroutine X(project)


!------------------------------------------------------------------------------
! X(psidprojectpsi) is used to calculate the contribution of the non-local part
! to the force acting on the ions.
!------------------------------------------------------------------------------
function X(psidprojectpsi)(mesh, p, n_projectors, dim, psi, periodic, ik) result(res)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: p(:)
  integer,           intent(in)    :: n_projectors
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:mesh%np, 1:dim)
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik
  R_TYPE :: res(3)


  integer ::  n_s, k, ip, idim
  R_TYPE, allocatable :: lpsi(:, :)
#if defined(HAVE_MPI)
  R_TYPE :: tmp
#endif

  call push_sub('epot_inc.psidprojectpsi')

  res = R_TOTYPE(M_ZERO)

  ! index labels the atom
  k = p(1)%iatom - 1 ! This way I make sure that k is not equal to p(1)%iatom

  do ip = 1, n_projectors

    if(p(ip)%iatom .ne. k) then
      n_s = p(ip)%n_s
      if(allocated(lpsi))  deallocate(lpsi)
      ALLOCATE( lpsi(n_s, dim), n_s*dim)

      do idim = 1, dim
        lpsi(1:n_s, idim)  = psi(p(ip)%jxyz(1:n_s), idim)*mesh%vol_pp(p(ip)%jxyz(1:n_s))
        if(periodic) lpsi(1:n_s, idim)  = lpsi(1:n_s, idim) * p(ip)%phases(1:n_s, ik)
      end do

      k = p(ip)%iatom
    end if

    select case (p(ip)%type)
    case (M_HGH)
      if (periodic) then
        res = res + X(hgh_dproject)(mesh, p(ip)%hgh_p, dim, lpsi, p(ip)%phases(:, ik))
      else
        res = res + X(hgh_dproject)(mesh, p(ip)%hgh_p, dim, lpsi)
      end if
    case (M_KB)
      if (periodic) then
        res = res + X(kb_dproject)(mesh, p(ip)%kb_p, dim, lpsi, p(ip)%phases(:, ik))
      else
        res = res + X(kb_dproject)(mesh, p(ip)%kb_p, dim, lpsi)
      end if
    case (M_RKB)
#ifdef R_TCOMPLEX
      !This can only be aplied to complex spinor wave-functions
      if (periodic) then
        res = res + rkb_dproject(mesh, p(ip)%rkb_p, lpsi, p(ip)%phases(:, ik))
      else
        res = res + rkb_dproject(mesh, p(ip)%rkb_p, lpsi)
      end if
#endif
    end select

  end do

  if(allocated(lpsi)) deallocate(lpsi)

  call pop_sub()
end function X(psidprojectpsi)
