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
! X(hgh_project) calculates the action of the projector hgh_p on the psi 
! wavefunction. The action of the projector hgh_p is defined as:
! \hat{hgh_p} |psi> = \sum_{ij}^3 p%h(i,j) |hgh_p%p(:, i)><hgh_p%p(:, j)|psi>
! The result is summed up to ppsi.
!
! If including the spin-orbit coupling there is another term to be added to ppsi:
! \sum_{ij}^3\sum{k}^3 p%k(i,j) |hgh_p%p(:, i)><hgh_p%lp(:, k, j)|\hat{S(k)}|psi>
!------------------------------------------------------------------------------
subroutine X(hgh_project)(mesh, sm, hgh_p, dim, psi, ppsi, reltype)
  type(mesh_t),          intent(in)  :: mesh
  type(submesh_t),       intent(in)  :: sm
  type(hgh_projector_t), intent(in)  :: hgh_p
  integer,               intent(in)  :: dim
  R_TYPE,                intent(in)  :: psi(:, :)  ! psi(hgh%n_s, dim)
  R_TYPE,                intent(out) :: ppsi(:, :) ! ppsi(hgh%n_s, dim)
  integer,               intent(in)  :: reltype

  integer :: n_s, i, j, idim
#ifdef R_TCOMPLEX
  integer :: k
#endif
  R_TYPE :: uvpsi
#ifdef R_TCOMPLEX
  CMPLX, allocatable :: lp_psi(:, :, :)
#endif

  n_s = hgh_p%n_s
  ppsi = R_TOTYPE(M_ZERO)

  do idim = 1, dim
    do j = 1, 3
      if (all(hgh_p%h(:, j) == M_ZERO)) cycle
      uvpsi = X(dsm_integrate_prod)(mesh, sm, psi(1:n_s, idim), hgh_p%p(1:n_s, j))

      do i = 1, 3
        if (hgh_p%h(i, j) == M_ZERO) cycle
#ifdef R_TCOMPLEX
        ppsi(1:n_s, idim) = ppsi(1:n_s, idim) + hgh_p%h(i, j)*uvpsi*hgh_p%p(1:n_s, i)
#else
        call lalg_axpy(n_s, hgh_p%h(i, j)*uvpsi, hgh_p%p(:, i), ppsi(:, idim))
#endif
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
          uvpsi = zdsm_integrate_prod(mesh, sm, psi(1:n_s, idim), hgh_p%lp(1:n_s, k, j))
          
          do i = 1, 3
            if (hgh_p%k(i, j) == M_ZERO) cycle
            lp_psi(1:n_s, k, idim) = lp_psi(1:n_s, k, idim) + &
                 hgh_p%k(i, j) * uvpsi * hgh_p%p(1:n_s, i)
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

subroutine X(hgh_project_bra)(mesh, sm, hgh_p, dim, reltype, psi, uvpsi)
  type(mesh_t),          intent(in)  :: mesh
  type(submesh_t),       intent(in)  :: sm
  type(hgh_projector_t), intent(in)  :: hgh_p
  integer,               intent(in)  :: dim
  integer,               intent(in)  :: reltype
  R_TYPE,                intent(in)  :: psi(:, :)
  R_TYPE,                intent(out) :: uvpsi(1:2, 1:4, 1:3)

  integer :: n_s, jj, idim, kk
  R_TYPE, allocatable :: bra(:, :, :)

#ifndef R_TCOMPLEX
  ASSERT(reltype == 0)
#endif

  n_s = hgh_p%n_s
  uvpsi = M_ZERO

  ALLOCATE(bra(1:n_s, 1:4, 1:3), n_s*4*3)

  bra = M_ZERO

  if(mesh%use_curvlinear) then
    do jj = 1, 3
      if(reltype == 1) then
        do kk = 1, 3
          bra(1:n_s, kk, jj) = hgh_p%lp(1:n_s, kk, jj)*mesh%vol_pp(sm%jxyz(1:n_s))
        end do
      end if
      bra(1:n_s, 4, jj) = hgh_p%p(1:n_s, jj)*mesh%vol_pp(sm%jxyz(1:n_s))
    end do
  else
    if(reltype == 1) bra(1:n_s, 1:3, 1:3) = hgh_p%lp(1:n_s, 1:3, 1:3)*mesh%vol_pp(1)
    bra(1:n_s, 4, 1:3) = hgh_p%p(1:n_s, 1:3)*mesh%vol_pp(1)
  end if

  do idim = 1, dim
    do kk = 1, 4
      do jj = 1, 3
        uvpsi(idim, kk, jj) = sum(psi(1:n_s, idim)*bra(1:n_s, kk, jj))
      end do
    end do
  end do

  deallocate(bra)

end subroutine X(hgh_project_bra)

subroutine X(hgh_project_ket)(mesh, sm, hgh_p, dim, reltype, uvpsi, ppsi)
  type(mesh_t),          intent(in)    :: mesh
  type(submesh_t),       intent(in)    :: sm
  type(hgh_projector_t), intent(in)    :: hgh_p
  integer,               intent(in)    :: dim
  integer,               intent(in)    :: reltype
  R_TYPE,                intent(in)    :: uvpsi(1:2, 1:4, 1:3)
  R_TYPE,                intent(inout) :: ppsi(:, :)

  integer :: n_s, ii, jj, idim
  integer :: kk
  CMPLX, allocatable :: lp_psi(:, :, :)

  n_s = hgh_p%n_s

  do idim = 1, dim
    do jj = 1, 3
      do ii = 1, 3
        ppsi(1:n_s, idim) = ppsi(1:n_s, idim) + hgh_p%h(ii, jj)*uvpsi(idim, 4, jj)*hgh_p%p(1:n_s, ii)
      end do
    end do
  end do
  
  if (reltype == 1) then

    ALLOCATE(lp_psi(n_s, 3, dim), n_s*3*dim)
    lp_psi = M_Z0

    do idim = 1, dim
      do kk = 1, 3
        do jj = 1, 3
          do ii = 1, 3
            lp_psi(1:n_s, kk, idim) = lp_psi(1:n_s, kk, idim) + hgh_p%k(ii, jj)*uvpsi(idim, kk, jj)*hgh_p%p(1:n_s, ii)
          end do
        end do
      end do
    end do

    ppsi(1:n_s, 1) = ppsi(1:n_s, 1) + M_zI*M_HALF*( lp_psi(1:n_s, 3, 1) + lp_psi(1:n_s, 1, 2) - M_zI*lp_psi(1:n_s, 2, 2))
    ppsi(1:n_s, 2) = ppsi(1:n_s, 2) + M_zI*M_HALF*(-lp_psi(1:n_s, 3, 2) + lp_psi(1:n_s, 1, 1) + M_zI*lp_psi(1:n_s, 2, 1))
    
    deallocate(lp_psi)
   
  end if
  
end subroutine X(hgh_project_ket)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
