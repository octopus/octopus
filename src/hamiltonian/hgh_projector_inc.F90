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
subroutine X(hgh_project)(mesh, sm, hgh_p, ll, lmax, dim, psi, ppsi, reltype)
  type(mesh_t),          intent(in)    :: mesh
  type(submesh_t),       intent(in)    :: sm
  type(hgh_projector_t), intent(in)    :: hgh_p(-lmax:)
  integer,               intent(in)    :: ll
  integer,               intent(in)    :: lmax
  integer,               intent(in)    :: dim
  R_TYPE,                intent(in)    :: psi(:, :)  !< (hgh\%n_s, dim)
  R_TYPE,                intent(inout) :: ppsi(:, :) !< (hgh\%n_s, dim)
  integer,               intent(in)    :: reltype

  integer :: mm
  R_TYPE :: uvpsi(1:dim, 1:3, -ll:ll)

  do mm = -ll,ll
    call X(hgh_project_bra)(mesh, sm, hgh_p(mm), dim, reltype, psi, uvpsi(1:dim, 1:3, mm))
  end do

  if(mesh%parallel_in_domains) call comm_allreduce(mesh%vp%comm, uvpsi)

  call X(hgh_project_ket)(hgh_p, ll, lmax, dim, reltype, uvpsi, ppsi)
  
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
  R_TYPE,               intent(out)  :: uvpsi(:,:) !< (dim, 3)

  integer :: n_s, jj, idim
  R_TYPE, allocatable :: bra(:, :)
  type(profile_t), save :: prof
  integer :: block_size, sp, ep, size

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

    SAFE_ALLOCATE(bra(1:n_s, 1:3))
    bra = M_ZERO

    do jj = 1, 3
      if(reltype == 1) then
        bra(1:n_s, jj) = R_CONJ(hgh_p%X(p)(1:n_s, jj))*mesh%vol_pp(sm%map(1:n_s))
      else
        bra(1:n_s, jj) = hgh_p%dp(1:n_s, jj)*mesh%vol_pp(sm%map(1:n_s))
      endif
    end do

    do idim = 1, dim
      do jj = 1, 3
        uvpsi(idim, jj) = sum(psi(1:n_s, idim)*bra(1:n_s, jj))
      end do
    end do

    SAFE_DEALLOCATE_A(bra)

  else

    if(reltype == 1) then
      do sp = 1, n_s, block_size
        ep = sp - 1 + min(block_size, n_s - sp + 1)
        size = min(block_size, n_s - sp + 1)
        do jj = 1, 3
          do idim = 1, dim
#ifdef R_TCOMPLEX
            uvpsi(idim, jj) = uvpsi(idim, jj) + blas_dot(size, hgh_p%zp(sp, jj), 1, psi(sp, idim), 1)*mesh%volume_element
#endif
          end do
        end do
      end do
    else
      do sp = 1, n_s, block_size
        ep = sp - 1 + min(block_size, n_s - sp + 1)
        size = min(block_size, n_s - sp + 1)
        do jj = 1, 3
          do idim = 1, dim
            uvpsi(idim, jj) = uvpsi(idim, jj) + sum(psi(sp:ep, idim)*hgh_p%dp(sp:ep, jj))*mesh%volume_element
          end do
        end do
      end do
    end if

  end if

  call profiling_out(prof)

end subroutine X(hgh_project_bra)

!-------------------------------------------------------------------------
!> THREADSAFE
subroutine X(hgh_project_ket)(hgh_p, ll, lmax, dim, reltype, uvpsi, ppsi)
  type(hgh_projector_t), intent(in)    :: hgh_p(-lmax:)
  integer,               intent(in)    :: ll
  integer,               intent(in)    :: lmax
  integer,               intent(in)    :: dim
  integer,               intent(in)    :: reltype
  R_TYPE,                intent(in)    :: uvpsi(:,:,-ll:) !< (dim, 3, 2*ll+1)
  R_TYPE,                intent(inout) :: ppsi(:, :)

  integer :: n_s, ii, jj, idim, mm
  R_TYPE :: weight(1:3,1:dim, -ll:ll)
  CMPLX  :: zweight(1:3,1:dim, -ll:ll)

  integer :: block_size, sp, ep, size
  type(profile_t), save :: prof

  call profiling_in(prof, "HGH_PROJECT_KET")

  ! This routine uses blocking to optimize cache usage. One block of
  ! |phi> is loaded in cache L1 and then then we calculate the dot
  ! product of it with the corresponding blocks of |psi_k>, next we
  ! load another block and do the same. This way we only have to load
  ! |psi> from the L2 or memory.
  block_size = hardware%X(block_size)

  n_s = hgh_p(0)%n_s

  do mm = -ll,ll
    weight(1:3, 1:dim, mm) = R_TOTYPE(M_ZERO)

    !We first compute for each value of ii and idim the weight of the projector hgh_p%p(1:n_s, ii)
    !Doing that we need to only apply once the each projector
    do idim = 1, dim
      do ii = 1, 3
        do jj = 1, 3
          weight(ii,idim,mm) = weight(ii,idim,mm) + hgh_p(mm)%h(ii, jj)*uvpsi(idim, jj, mm)
        end do
      end do
    end do
  
    if (reltype == 1) then

      !We compute for each value of ii, kk and idim the weight of the projector hgh_p%p(1:n_s, ii)
      !Doing that we need to only apply once the each projector 
      zweight(1:3,1:dim,mm) = M_z0
      
      do ii = 1, 3
        do jj = 1, 3

          zweight(ii,1,mm) =  zweight(ii,1,mm) + mm*hgh_p(mm)%k(ii, jj)*uvpsi(1, jj, mm)
          zweight(ii,2,mm) =  zweight(ii,2,mm) - mm*hgh_p(mm)%k(ii, jj)*uvpsi(2, jj, mm)

          if(mm < ll) then 
            zweight(ii,1,mm) =  zweight(ii,1,mm) + hgh_p(mm)%k(ii, jj)*sqrt(real(ll*(ll+1)-mm*(mm+1)))&
               * uvpsi(2, jj, mm+1) 
          end if

          if(-mm < ll) then
            zweight(ii,2,mm) =  zweight(ii,2,mm) + hgh_p(mm)%k(ii, jj)*sqrt(real(ll*(ll+1)-mm*(mm-1)))&
               * uvpsi(1, jj, mm-1)
          end if

        end do
      end do
    end if
  end do


  if (reltype == 1) then
    do sp = 1, n_s, block_size
      size = min(block_size, n_s - sp + 1) 
      !We now apply the projectors
      do mm=-ll, ll
        do idim = 1, dim
          do ii = 1, 3
            !If we have SOC, we can only have complex wfns 
#ifdef R_TCOMPLEX
            call blas_axpy(size, (M_HALF*zweight(ii,idim,mm)+weight(ii,idim,mm)), hgh_p(mm)%zp(sp, ii), 1, ppsi(sp, idim), 1)
#endif 
          end do
        end do
      end do
     end do
  else

    do sp = 1, n_s, block_size
      size = min(block_size, n_s - sp + 1)
      ep = sp -1 + size
      !We now apply the projectors
      do mm=-ll, ll 
        do idim = 1, dim
          do ii = 1, 3
           !In case of complex wfns, we cannot use blas
#ifdef R_TCOMPLEX
            ppsi(sp:ep, idim) = ppsi(sp:ep, idim) + weight(ii,idim,mm)*hgh_p(mm)%dp(sp:ep, ii)
#else
            call blas_axpy(size, weight(ii,idim,mm), hgh_p(mm)%dp(sp, ii), 1, ppsi(sp, idim), 1)
#endif
          end do
        end do
      end do
    end do
  end if


  call profiling_out(prof)  

end subroutine X(hgh_project_ket)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
