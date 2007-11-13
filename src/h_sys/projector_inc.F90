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
subroutine X(project_psi)(mesh, pj, dim, psi, ppsi, reltype)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: ppsi(:, :)  ! ppsi(1:mesh%np, dim)
  integer,           intent(in)    :: reltype

  integer :: n_s, idim
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)

  call push_sub('projector_inc.project_psi')

  n_s = pj%sphere%ns
  
  ALLOCATE(lpsi(n_s, dim),  n_s*dim)
  ALLOCATE(plpsi(n_s, dim), n_s*dim)
  
  do idim = 1, dim
    lpsi(1:n_s, idim)  = psi(pj%sphere%jxyz(1:n_s), idim)
  end do

  call X(project_sphere)(mesh, pj, dim, lpsi, plpsi, reltype)

  do idim = 1, dim
    ppsi(pj%sphere%jxyz(1:n_s), idim) = ppsi(pj%sphere%jxyz(1:n_s), idim) + plpsi(1:n_s, idim)
  end do

  if(allocated(lpsi))   deallocate(lpsi)
  if(allocated(plpsi))  deallocate(plpsi)

  call pop_sub()
end subroutine X(project_psi)

!------------------------------------------------------------------------------
! X(psia_project_psib) calculates <psia|projector|psib>
R_TYPE function X(psia_project_psib)(mesh, pj, dim, psia, psib, reltype) result(apb)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psia(:, :)  ! psia(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: psib(:, :)  ! psib(1:mesh%np, dim)
  integer,           intent(in)    :: reltype

  integer ::  n_s, idim
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)

  call push_sub('projector_inc.psia_project_psib')

  n_s = pj%sphere%ns

  ALLOCATE(lpsi(n_s, dim),  n_s*dim)
  ALLOCATE(plpsi(n_s, dim), n_s*dim)
  
  do idim = 1, dim
    lpsi(1:n_s, idim)  = psib(pj%sphere%jxyz(1:n_s), idim)
  end do

  call X(project_sphere)(mesh, pj, dim, lpsi, plpsi, reltype)

  apb = M_ZERO
  do idim = 1, dim
    plpsi(1:n_s, idim) = R_CONJ(psia(pj%sphere%jxyz(1:n_s), idim)) * plpsi(1:n_s, idim)
    apb = apb + X(sm_integrate)(mesh, pj%sphere, plpsi(1:n_s, idim))
  end do

  deallocate(lpsi, plpsi)

  call pop_sub()
end function X(psia_project_psib)

subroutine X(project_sphere)(mesh, pj, dim, psi, ppsi, reltype)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:n_s, dim)
  R_TYPE,            intent(inout) :: ppsi(:, :)  ! ppsi(1:n_s, dim)
  integer,           intent(in)    :: reltype

  integer :: n_s, idim

  call push_sub('projector_inc.project_sphere')

  n_s = pj%sphere%ns

  if(pj%type == M_LOCAL) then

    do idim = 1, dim
      ppsi(1:n_s, idim) = pj%local_p%v(1:n_s) * psi(1:n_s, idim)
    end do

  else

    select case (pj%type)
    case (M_HGH)
      call X(hgh_project)(mesh, pj%sphere, pj%hgh_p, dim, psi, ppsi, reltype)
    case (M_KB)
      call X(kb_project)(mesh, pj%sphere, pj%kb_p, dim, psi, ppsi)
    case (M_RKB)
#ifdef R_TCOMPLEX
      !This can only be aplied to complex spinor wave-functions
      call rkb_project(mesh, pj%sphere, pj%rkb_p, psi, ppsi)
#endif
    end select
    
  end if

  call pop_sub()
end subroutine X(project_sphere)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
