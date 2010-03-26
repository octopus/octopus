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
  CMPLX, intent(out)       :: gamma(:, :)

  integer :: imesh, imeshp, icoord, icoordp, jdim
  integer, allocatable :: ix(:), ix_1part(:), ixp(:)
  FLOAT :: volume_element
  R_TYPE, allocatable :: psi_global(:)

  call push_sub('modelmb_density_matrix_inc.Xmf_calculate_gamma')

  ! In case of running parallel in domains, we need to operate psi_global, which 
  ! contains the full wavefunction after "gathering" all the domains.
  SAFE_ALLOCATE(psi_global(1:mesh%np_part_global))
  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    call X(vec_allgather)(mesh%vp, psi_global, psi)
#endif
  else
    psi_global(1:mesh%np_part_global) = psi(1:mesh%np_part_global)
  end if

  SAFE_ALLOCATE(ix(1:MAX_DIM))
  SAFE_ALLOCATE(ixp(1:MAX_DIM))

  volume_element = 1.0d0
  do jdim = 1, MAX_DIM
    if (mesh%spacing(jdim) > 1.e-10) volume_element=volume_element*mesh%spacing(jdim)
  end do
  do jdim = (ikeeppart - 1)*mb_1part%ndim1part + 1, ikeeppart*mb_1part%ndim1part
    if (mesh%spacing(jdim) > 1.e-10) volume_element = volume_element/mesh%spacing(jdim)
  end do

  gamma = R_TOTYPE(M_ZERO)

  SAFE_ALLOCATE(ix_1part(1:mb_1part%ndim1part))
  xloop: do imesh = 1, mesh%np_global
    ! find coordinates of present point in full MAX_DIM space
     call index_to_coords(mesh%idx, mesh%sb%dim, imesh, ix)
     
     ! find index of present coordinates for particle ikeeppart
     ix_1part = ix((ikeeppart - 1)*mb_1part%ndim1part + 1:ikeeppart*mb_1part%ndim1part)
     call hypercube_x_to_i(mb_1part%hypercube_1part, mb_1part%ndim1part, &
       mb_1part%nr_1part, mb_1part%enlarge_1part(1), ix_1part, icoord)
     
     ! run over all possible position indices of particle prime
     ixp = ix
     xprimeloop: do icoordp = 1, mb_1part%npt_part
       
       ! find equivalent position of particle prime
       call hypercube_i_to_x(mb_1part%hypercube_1part, mb_1part%ndim1part, &
         mb_1part%nr_1part, mb_1part%enlarge_1part(1), icoordp, ix_1part)
       
       ! change coordinates of particle ikeeppart only 
       ixp((ikeeppart - 1)*mb_1part%ndim1part + 1:ikeeppart*mb_1part%ndim1part) = ix_1part
       
       ! find new index for general point prime
       imeshp = index_from_coords(mesh%idx, mesh%sb%dim, ixp)
       
       ! accumulate into density matrix
       gamma(icoord, icoordp) = gamma(icoord, icoordp) + &
         nparticles_densmat*volume_element*psi_global(imesh)*R_CONJ(psi_global (imeshp))
     end do xprimeloop
  end do xloop

  SAFE_DEALLOCATE_A(psi_global)
  SAFE_DEALLOCATE_A(ix)
  SAFE_DEALLOCATE_A(ixp)
  SAFE_DEALLOCATE_A(ix_1part)

  call pop_sub('modelmb_density_matrix_inc.Xmf_calculate_gamma')
end subroutine X(mf_calculate_gamma)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
