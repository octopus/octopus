!! Copyright (C) 2009 M. Marques, A. Castro, M. Verstraete
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

!------------------------------------------------------------
! This routine calculates the one-body density matrix gamma 
! for particle ikeeppart, used in higher dimensional model
! hamiltonian calculations (MJV, NH) 
!------------------------------------------------------------
subroutine X(mf_calculate_gamma)(ikeeppart, mb_1part, nparticles_densmat, &
     mesh, psi, gamma)
  integer, intent(in)      :: ikeeppart
  integer, intent(in)      :: nparticles_densmat
  type(modelmb_1part_t), intent(in) :: mb_1part
  type(mesh_t), intent(in) :: mesh
  R_TYPE, intent(in)       :: psi(:)
  R_TYPE, intent(out)       :: gamma(:, :)

  integer :: icoord, icoordp, icoord_diff
  integer :: jdim, ip_global, ip, ipp_global
  integer, allocatable :: ix(:), ix_1part(:), ixp(:)
  integer, allocatable :: forward_map_gamma(:)
  integer, allocatable :: icoord_map(:)
  FLOAT :: volume_element
  R_TYPE, allocatable :: psi_p(:,:,:)
  type(batch_t) :: wfbatch

  PUSH_SUB(X(mf_calculate_gamma))

  SAFE_ALLOCATE(ix(1:MAX_DIM))
  SAFE_ALLOCATE(ixp(1:MAX_DIM))
  SAFE_ALLOCATE(psi_p(1:mesh%np,1,1))
  SAFE_ALLOCATE(forward_map_gamma(1:mesh%np_global))
  SAFE_ALLOCATE(icoord_map(1:mesh%np))

  volume_element = 1.0d0
  do jdim = 1, MAX_DIM
    if (mesh%spacing(jdim) > 1.e-10) volume_element=volume_element*mesh%spacing(jdim)
  end do
  do jdim = (ikeeppart - 1)*mb_1part%ndim1part + 1, ikeeppart*mb_1part%ndim1part
    if (mesh%spacing(jdim) > 1.e-10) volume_element = volume_element/mesh%spacing(jdim)
  end do

  ASSERT (ubound(gamma, dim=1) == mb_1part%npt)
  ASSERT (ubound(gamma, dim=2) == mb_1part%npt)
  ASSERT (ubound(psi,   dim=1) >= mesh%np)

  gamma = R_TOTYPE(M_ZERO)

  SAFE_ALLOCATE(ix_1part(1:mb_1part%ndim1part))

  ! loop over the points of psi we have locally
  do ip = 1, mesh%np
    ! find global index
    ip_global = ip
    if (mesh%parallel_in_domains) ip_global = mesh%vp%local(ip + mesh%vp%xlocal(mesh%vp%partno) - 1)

    ! find coordinates of present point in full MAX_DIM space
    call index_to_coords(mesh%idx, mesh%sb%dim, ip_global, ix)
    
    ! find index of present coordinates for particle ikeeppart
    ix_1part = ix((ikeeppart - 1)*mb_1part%ndim1part + 1:ikeeppart*mb_1part%ndim1part)
    call hypercube_x_to_i(mb_1part%hypercube_1part, mb_1part%ndim1part, &
      mb_1part%nr_1part, mb_1part%enlarge_1part(1), ix_1part, icoord)

    icoord_map(ip) = icoord
  end do

  ! loop over the difference between
  !  * the x` position of the kept particle, which we will impose below
  !  * and the x position of the local particle
  do icoord_diff = 1, mb_1part%npt

    ! make global map of all points to their image with x`=icoord + icoord_diff (modulus npt of course)
    ! this map will be the same on all processors
    do ip_global = 1, mesh%np_global

      ! find coordinates of present point in full MAX_DIM space
      call index_to_coords(mesh%idx, mesh%sb%dim, ip_global, ix)
      
      ! prime position will be identical to ix, apart from the ikeeppart particle
      ixp = ix
        
      ! find index of present coordinates for particle ikeeppart
      ix_1part = ix((ikeeppart - 1)*mb_1part%ndim1part + 1:ikeeppart*mb_1part%ndim1part)
      call hypercube_x_to_i(mb_1part%hypercube_1part, mb_1part%ndim1part, &
        mb_1part%nr_1part, mb_1part%enlarge_1part(1), ix_1part, icoord)

      icoordp = icoord + icoord_diff
      if (icoordp > mb_1part%npt) icoordp = icoordp - mb_1part%npt

      ! find equivalent position of particle prime
      call hypercube_i_to_x(mb_1part%hypercube_1part, mb_1part%ndim1part, &
        mb_1part%nr_1part, mb_1part%enlarge_1part(1), icoordp, ix_1part)
        
      ! change coordinates of particle ikeeppart only 
      ixp((ikeeppart - 1)*mb_1part%ndim1part + 1:ikeeppart*mb_1part%ndim1part) = ix_1part
        
      ! find new index for general point prime
      ipp_global = index_from_coords(mesh%idx, mesh%sb%dim, ixp)
      forward_map_gamma(ipp_global) = ip_global
    end do

    ! use map to recover the corresponding np points for psi
    if (mesh%parallel_in_domains) then
      psi_p(:,1,1) = psi(1:mesh%np)
      call batch_init (wfbatch, 1, 1, 1, psi_p)
      call X(mesh_batch_exchange_points) (mesh, wfbatch, forward_map=forward_map_gamma)
      call batch_end(wfbatch)
    else
      psi_p(forward_map_gamma(1:mesh%np),1,1) = psi(1:mesh%np)
    end if

    ! accumulate in gamma each pair of positions local processor now has
    do ip = 1, mesh%np
      icoordp = icoord_map(ip) + icoord_diff
      if (icoordp > mb_1part%npt) icoordp = icoordp - mb_1part%npt

      gamma(icoord_map(ip), icoordp) = gamma(icoord_map(ip), icoordp) + &
         nparticles_densmat*volume_element*psi(ip)*R_CONJ(psi_p (ip,1,1))
    end do
  end do

  if(mesh%parallel_in_domains) call comm_allreduce(mesh%mpi_grp%comm, gamma, &
&      dim = (/mb_1part%npt, mb_1part%npt/))

  SAFE_DEALLOCATE_A(forward_map_gamma)
  SAFE_DEALLOCATE_A(icoord_map)
  SAFE_DEALLOCATE_A(psi_p)
  SAFE_DEALLOCATE_A(ix)
  SAFE_DEALLOCATE_A(ixp)
  SAFE_DEALLOCATE_A(ix_1part)

  POP_SUB(X(mf_calculate_gamma))
end subroutine X(mf_calculate_gamma)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
