!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Verstraete
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
! This routine calculates the density rho for particle 
! ikeeppart, used in higher-dimensional model
! hamiltonian calculations (MJV, NH) 
!------------------------------------------------------------
subroutine X(modelmb_density_calculate)(ikeeppart, mb_1part, nparticles_dens, &
     mesh, psi, rho)
  integer, intent(in)      :: ikeeppart
  integer, intent(in)      :: nparticles_dens
  type(mesh_t), intent(in) :: mesh
  type(modelmb_1part_t), intent(in) :: mb_1part
  R_TYPE, intent(in)       :: psi(:)
  FLOAT, intent(out)       :: rho(:)

  integer :: imesh, icoord, jdim, ip_global
  integer, allocatable :: ix(:), ix_1part(:), ixp(:)
  FLOAT :: volume_element

  PUSH_SUB(X(modelmb_density_calculate))

  SAFE_ALLOCATE(ix(1:MAX_DIM))
  SAFE_ALLOCATE(ixp(1:MAX_DIM))

  volume_element = 1.0d0
  do jdim = 1, MAX_DIM
    if (mesh%spacing(jdim) > 1.e-10) volume_element=volume_element*mesh%spacing(jdim)
  end do
  do jdim = (ikeeppart - 1)*mb_1part%ndim1part + 1, ikeeppart*mb_1part%ndim1part
    if (mesh%spacing(jdim) > 1.e-10) volume_element = volume_element/mesh%spacing(jdim)
  end do

  rho = R_TOTYPE(M_ZERO)

  SAFE_ALLOCATE(ix_1part(1:mb_1part%ndim1part))
  xloop: do imesh = 1, mesh%np
! go from local to global point index

    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
      ip_global = mesh%vp%local(imesh + mesh%vp%xlocal(mesh%vp%partno) - 1)
#endif
    else
      ip_global = imesh
    end if

! find coordinates of present point in full dimensional space
    call index_to_coords(mesh%idx, mesh%sb%dim, ip_global, ix)

! find index of present coordinates for particle ikeeppart
    ix_1part = ix((ikeeppart - 1)*mb_1part%ndim1part + 1:ikeeppart*mb_1part%ndim1part)
    call hypercube_x_to_i(mb_1part%hypercube_1part, mb_1part%ndim1part, mb_1part%nr_1part, mb_1part%enlarge_1part(1), &
      ix_1part, icoord)

! accumulate into scalar density
    rho(icoord) = rho(icoord) + &
      nparticles_dens*volume_element*TOFLOAT(psi(imesh)*R_CONJ(psi(imesh)))
  end do xloop

! need to collect all the rho() here in parallel case
  if(mesh%parallel_in_domains) call comm_allreduce(mesh%mpi_grp%comm, rho, dim = mb_1part%npt_part)

  SAFE_DEALLOCATE_A(ix)
  SAFE_DEALLOCATE_A(ixp)
  SAFE_DEALLOCATE_A(ix_1part)

  POP_SUB(X(modelmb_density_calculate))

end subroutine X(modelmb_density_calculate)

subroutine X(modelmb_2partdensity_calculate)(ikeeppart1, ikeeppart2, mb_1part, & 
        nparticles_dens, mesh, psi, rho2)
  integer, intent(in)      :: ikeeppart1, ikeeppart2
  integer, intent(in)      :: nparticles_dens
  type(mesh_t), intent(in) :: mesh
  type(modelmb_1part_t), intent(in) :: mb_1part
  R_TYPE, intent(in)       :: psi(:)
  FLOAT, intent(out)       :: rho2(:,:)

  integer :: imesh, icoord1, jdim, icoord2
  integer, allocatable :: ix(:), ix_1part(:), ixp(:), ixp_1part(:)
  FLOAT :: volume_element
  R_TYPE, allocatable :: psi_global(:)

  PUSH_SUB(X(modelmb_2partdensity_calculate))

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
  do jdim = (ikeeppart1 - 1)*mb_1part%ndim1part + 1, ikeeppart1*mb_1part%ndim1part
    if (mesh%spacing(jdim) > 1.e-10) volume_element = volume_element/mesh%spacing(jdim)
  end do
  do jdim = (ikeeppart2 - 1)*mb_1part%ndim1part + 1, ikeeppart2*mb_1part%ndim1part
    if (mesh%spacing(jdim) > 1.e-10) volume_element = volume_element/mesh%spacing(jdim)
  end do
  
  rho2 = R_TOTYPE(M_ZERO)

  SAFE_ALLOCATE(ix_1part(1:mb_1part%ndim1part))
  SAFE_ALLOCATE(ixp_1part(1:mb_1part%ndim1part))
  
  do imesh = 1, mesh%np_global
    ! find coordinates of present point in full MAX_DIM space
    call index_to_coords(mesh%idx, mesh%sb%dim, imesh, ix)
    ! find index of present coordinates for particle ikeeppart1
    ix_1part = ix((ikeeppart1 - 1)*mb_1part%ndim1part + 1:ikeeppart1*mb_1part%ndim1part)
    call hypercube_x_to_i(mb_1part%hypercube_1part, mb_1part%ndim1part, mb_1part%nr_1part, &
      mb_1part%enlarge_1part(1), ix_1part, icoord1)
    ! find index of present coordinates for particle ikeeppart2
    ixp_1part = ixp((ikeeppart2 - 1)*mb_1part%ndim1part + 1:ikeeppart2*mb_1part%ndim1part)
    call hypercube_x_to_i(mb_1part%hypercube_1part, mb_1part%ndim1part, mb_1part%nr_1part, & 
      mb_1part%enlarge_1part(1), ixp_1part, icoord2)
    ! accumulate into two-particle density
    rho2(icoord1, icoord2) = rho2(icoord1, icoord2) + nparticles_dens*(nparticles_dens-1)* &
      volume_element*psi_global(imesh)*R_CONJ(psi_global(imesh))
  end do 

  SAFE_DEALLOCATE_A(psi_global)
  SAFE_DEALLOCATE_A(ix)
  SAFE_DEALLOCATE_A(ixp)
  SAFE_DEALLOCATE_A(ix_1part)
  SAFE_DEALLOCATE_A(ixp_1part)
  POP_SUB(X(modelmb_2partdensity_calculate))

end subroutine X(modelmb_2partdensity_calculate)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
