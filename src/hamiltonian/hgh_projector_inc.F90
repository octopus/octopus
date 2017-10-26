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
!> X(hgh_project) calculates the action of the projector hgh_p on the psi 
!! wavefunction. The action of the projector hgh_p is defined as:
!! \f[
!! \hat{hgh_p} |psi> = \sum_{ij}^3 p\%h(i,j) |hgh_p\%p(:, i)><hgh_p\%p(:, j)|psi>
!! \f]
!! The result is summed up to ppsi.
!!
!! If including the spin-orbit coupling there is another term to be added to ppsi:
!! \f[
!! \sum_{ij}^3\sum{k}^3 p\%k(i,j) |hgh_p\%p(:, i)><hgh_p\%lp(:, k, j)|\hat{S(k)}|psi>
!! \f]
!------------------------------------------------------------------------------
subroutine X(hgh_project)(mesh, sm, hgh_p, dim, psi, ppsi, reltype)
  type(mesh_t),          intent(in)    :: mesh
  type(submesh_t),       intent(in)    :: sm
  type(hgh_projector_t), intent(in)    :: hgh_p
  integer,               intent(in)    :: dim
  R_TYPE,                intent(in)    :: psi(:, :)  !< (hgh\%n_s, dim)
  R_TYPE,                intent(inout) :: ppsi(:, :) !< (hgh\%n_s, dim)
  integer,               intent(in)    :: reltype

  R_TYPE :: uvpsi(1:2, 1:12)

  call X(hgh_project_bra)(mesh, sm, hgh_p, dim, reltype, psi, uvpsi)

  if(mesh%parallel_in_domains) call comm_allreduce(mesh%vp%comm, uvpsi)

  call X(hgh_project_ket)(hgh_p, dim, reltype, uvpsi, ppsi)
  
end subroutine X(hgh_project)

!-------------------------------------------------------------------------
!> THREADSAFE
subroutine X(hgh_project_bra)(mesh, sm, hgh_p, dim, reltype, psi, uvpsi)
  type(mesh_t),          intent(in)  :: mesh
  type(submesh_t),       intent(in)  :: sm
  type(hgh_projector_t), intent(in)  :: hgh_p
  integer,               intent(in)  :: dim
  integer,               intent(in)  :: reltype
  R_TYPE,                intent(in)  :: psi(:, :)
  R_TYPE,                intent(out) :: uvpsi(:,:) !< (dim, 12)

  integer :: n_s, jj, idim, kk
  R_TYPE, allocatable :: bra(:, :, :)
  type(profile_t), save :: prof
  integer :: block_size, sp, ep

  call profiling_in(prof, "HGH_PROJECT_BRA")

  ! This routine uses blocking to optimize cache usage. One block of
  ! |phi> is loaded in cache L1 and then then we calculate the dot
  ! product of it with the corresponding blocks of |psi_k>, next we
  ! load another block and do the same. This way we only have to load
  ! |psi> from the L2 or memory.
  block_size = hardware%X(block_size)


#ifndef R_TCOMPLEX
  ASSERT(reltype == 0)
#endif

  n_s = hgh_p%n_s
  uvpsi = M_ZERO

  if(mesh%use_curvilinear) then

    SAFE_ALLOCATE(bra(1:n_s, 1:4, 1:3))
    bra = M_ZERO

    do jj = 1, 3
      if(reltype == 1) then
        do kk = 1, 3
          bra(1:n_s, kk, jj) = hgh_p%lp(1:n_s, kk, jj)*mesh%vol_pp(sm%map(1:n_s))
        end do
      end if
      bra(1:n_s, 4, jj) = hgh_p%p(1:n_s, jj)*mesh%vol_pp(sm%map(1:n_s))
    end do

    do idim = 1, dim
      do kk = 1, 4 
        do jj = 1, 3
          uvpsi(idim, hgh_index(kk, jj)) = sum(psi(1:n_s, idim)*bra(1:n_s, kk, jj))
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(bra)

  else

    do idim = 1, dim
      do sp = 1, n_s, block_size
        ep = sp - 1 + min(block_size, n_s - sp + 1)
        do jj = 1, 3
          uvpsi(idim, hgh_index(4, jj)) = uvpsi(idim, hgh_index(4, jj)) + &
                            sum(psi(sp:ep, idim)*hgh_p%p(sp:ep, jj))*mesh%volume_element
          if(reltype == 1) then
            do kk = 1, 3
              uvpsi(idim, hgh_index(kk, jj)) = uvpsi(idim, hgh_index(kk, jj)) + &
                  sum(psi(sp:ep, idim)*hgh_p%lp(sp:ep, kk, jj))*mesh%volume_element
            end do
          end if
        end do
      end do
    end do
  end if

  call profiling_out(prof)

end subroutine X(hgh_project_bra)

!-------------------------------------------------------------------------
!> THREADSAFE
subroutine X(hgh_project_ket)(hgh_p, dim, reltype, uvpsi, ppsi)
  type(hgh_projector_t), intent(in)    :: hgh_p
  integer,               intent(in)    :: dim
  integer,               intent(in)    :: reltype
  R_TYPE,                intent(in)    :: uvpsi(:,:) !< (dim, 12)
  R_TYPE,                intent(inout) :: ppsi(:, :)

  integer :: n_s, ii, jj, idim
  integer :: kk
  CMPLX, allocatable :: lp_psi(:, :, :)
  R_TYPE :: weight
  CMPLX  :: zweight

  type(profile_t), save :: prof

  call profiling_in(prof, "HGH_PROJECT_KET")

  n_s = hgh_p%n_s

  do idim = 1, dim
    do ii = 1, 3
      weight = R_TOTYPE(M_ZERO)
      do jj = 1, 3
        weight = weight + hgh_p%h(ii, jj)*uvpsi(idim, hgh_index(4, jj))
      end do
      ppsi(1:n_s, idim) = ppsi(1:n_s, idim) + weight*hgh_p%p(1:n_s, ii)
    end do
  end do
  
  if (reltype == 1) then

    SAFE_ALLOCATE(lp_psi(1:n_s, 1:3, 1:dim))
    lp_psi = M_Z0

    do idim = 1, dim
      do kk = 1, 3
        do ii = 1, 3
          zweight = M_z0
          do jj = 1, 3
            zweight = zweight + hgh_p%k(ii, jj)*uvpsi(idim, hgh_index(kk, jj))
          end do
          lp_psi(1:n_s, kk, idim) = lp_psi(1:n_s, kk, idim) + zweight*hgh_p%p(1:n_s, ii)
        end do
      end do
    end do

    ppsi(1:n_s, 1) = ppsi(1:n_s, 1) + M_zI*M_HALF*( lp_psi(1:n_s, 3, 1) + lp_psi(1:n_s, 1, 2) - M_zI*lp_psi(1:n_s, 2, 2))
    ppsi(1:n_s, 2) = ppsi(1:n_s, 2) + M_zI*M_HALF*(-lp_psi(1:n_s, 3, 2) + lp_psi(1:n_s, 1, 1) + M_zI*lp_psi(1:n_s, 2, 1))
    
    SAFE_DEALLOCATE_A(lp_psi)
   
  end if

  call profiling_out(prof)  

end subroutine X(hgh_project_ket)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
