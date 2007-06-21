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
!! $Id$

!------------------------------------------------------------------------------
! X(project_psi) calculates the action of a projector on the psi wavefunction.
! The result is summed up to ppsi
subroutine X(project_psi)(mesh, pj, dim, psi, ppsi, reltype, periodic, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: ppsi(:, :)  ! ppsi(1:mesh%np, dim)
  integer,           intent(in)    :: reltype
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik

  integer :: n_s, idim
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)

  call push_sub('projector_inc.project')

  n_s = pj%sphere%ns
  
  ALLOCATE(lpsi(n_s, dim),  n_s*dim)
  ALLOCATE(plpsi(n_s, dim), n_s*dim)
  
  do idim = 1, dim
    lpsi(1:n_s, idim)  = psi(pj%sphere%jxyz(1:n_s), idim)
  end do

  call X(project_sphere)(mesh, pj, dim, lpsi, plpsi, reltype, periodic, ik)

  do idim = 1, dim
    ppsi(pj%sphere%jxyz(1:n_s), idim) = ppsi(pj%sphere%jxyz(1:n_s), idim) + plpsi(1:n_s, idim)
  end do

  if(allocated(lpsi))   deallocate(lpsi)
  if(allocated(plpsi))  deallocate(plpsi)

  call pop_sub()
end subroutine X(project_psi)

!------------------------------------------------------------------------------
! X(psia_project_psib) calculates <psia|projector|psib>
R_TYPE function X(psia_project_psib)(mesh, pj, dim, psia, psib, reltype, periodic, ik) result(apb)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psia(:, :)  ! psia(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: psib(:, :)  ! psib(1:mesh%np, dim)
  integer,           intent(in)    :: reltype
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik

  integer ::  n_s, idim
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)

  call push_sub('projector_inc.psia_project_psib')

  n_s = pj%sphere%ns

  ALLOCATE(lpsi(n_s, dim),  n_s*dim)
  ALLOCATE(plpsi(n_s, dim), n_s*dim)
  
  do idim = 1, dim
    lpsi(1:n_s, idim)  = psib(pj%sphere%jxyz(1:n_s), idim)
  end do

  call X(project_sphere)(mesh, pj, dim, lpsi, plpsi, reltype, periodic, ik)

  apb = M_ZERO
  do idim = 1, dim
    plpsi(1:n_s, 1) = R_CONJ(psia(pj%sphere%jxyz(1:n_s), idim)) * plpsi(1:n_s, idim)
  end do
  apb = apb + X(sm_integrate)(mesh, pj%sphere, plpsi(1:n_s, 1))

  deallocate(lpsi, plpsi)

  call pop_sub()
end function X(psia_project_psib)

subroutine X(project_sphere)(mesh, pj, dim, psi, ppsi, reltype, periodic, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:n_s, dim)
  R_TYPE,            intent(inout) :: ppsi(:, :)  ! ppsi(1:n_s, dim)
  integer,           intent(in)    :: reltype
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik

  integer :: n_s, idim
  R_TYPE, allocatable :: ph_psi(:, :)

  call push_sub('projector_inc.project')

  n_s = pj%sphere%ns

  if(pj%type == M_LOCAL) then

    do idim = 1, dim
      ppsi(1:n_s, idim) = pj%local_p%v(1:n_s) * psi(1:n_s, idim)
    end do

  else

    if (periodic) then
      ALLOCATE(ph_psi(n_s, dim),  n_s*dim)

      do idim = 1, dim
        ph_psi(1:n_s, idim)  = psi(1:n_s, idim) * pj%phases(1:n_s, ik)
      end do
    end if

    select case (pj%type)
    case (M_HGH)
      if (periodic) then
        call X(hgh_project)(mesh, pj%sphere, pj%hgh_p, dim, ph_psi, ppsi, reltype, pj%phases(:, ik))
      else
        call X(hgh_project)(mesh, pj%sphere, pj%hgh_p, dim, psi, ppsi, reltype)
      end if
    case (M_KB)
      if (periodic) then
        call X(kb_project)(mesh, pj%sphere, pj%kb_p, dim, ph_psi, ppsi, pj%phases(:, ik))
      else
        call X(kb_project)(mesh, pj%sphere, pj%kb_p, dim, psi, ppsi)
      end if
    case (M_RKB)
#ifdef R_TCOMPLEX
      !This can only be aplied to complex spinor wave-functions
      if (periodic) then
        call rkb_project(mesh, pj%sphere, pj%rkb_p, ph_psi, ppsi, pj%phases(:, ik))
      else
        call rkb_project(mesh, pj%sphere, pj%rkb_p, psi, ppsi)
      end if
#endif
    end select

  end if

  if(periodic) deallocate(ph_psi)

  call pop_sub()
end subroutine X(project_sphere)

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


  integer ::  n_s, ip, idim
  R_TYPE, allocatable :: lpsi(:, :)

  call push_sub('projector_inc.psidprojectpsi')

  res = R_TOTYPE(M_ZERO)

  ! index labels the atom
  do ip = 1, n_projectors
    if(p(ip)%type == M_LOCAL) cycle

    n_s = p(ip)%sphere%ns
    if(allocated(lpsi))  deallocate(lpsi)
    ALLOCATE( lpsi(n_s, dim), n_s*dim)
    
    do idim = 1, dim
      lpsi(1:n_s, idim)  = psi(p(ip)%sphere%jxyz(1:n_s), idim)
      if(periodic) lpsi(1:n_s, idim)  = lpsi(1:n_s, idim) * p(ip)%phases(1:n_s, ik)
    end do
    
    select case (p(ip)%type)
    case (M_HGH)
      if (periodic) then
        res = res + X(hgh_dproject)(mesh, p(ip)%sphere, p(ip)%hgh_p, dim, lpsi, p(ip)%phases(:, ik))
      else
        res = res + X(hgh_dproject)(mesh, p(ip)%sphere, p(ip)%hgh_p, dim, lpsi)
      end if
    case (M_KB)
      if (periodic) then
        res = res + X(kb_dproject)(mesh, p(ip)%sphere, p(ip)%kb_p, dim, lpsi, p(ip)%phases(:, ik))
      else
        res = res + X(kb_dproject)(mesh, p(ip)%sphere, p(ip)%kb_p, dim, lpsi)
      end if
    case (M_RKB)
#ifdef R_TCOMPLEX
      !This can only be aplied to complex spinor wave-functions
      if (periodic) then
        res = res + rkb_dproject(mesh, p(ip)%sphere, p(ip)%rkb_p, lpsi, p(ip)%phases(:, ik))
      else
        res = res + rkb_dproject(mesh, p(ip)%sphere, p(ip)%rkb_p, lpsi)
      end if
#endif
    end select

  end do

  if(allocated(lpsi)) deallocate(lpsi)

  call pop_sub()
end function X(psidprojectpsi)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
