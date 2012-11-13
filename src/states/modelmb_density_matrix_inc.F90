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


! ---------------------------------------------------------
subroutine X(modelmb_density_matrix_write)(gr, st, wf, mm, denmat)
  type(grid_t),           intent(in) :: gr
  type(states_t),         intent(in) :: st
  R_TYPE,                 intent(in) :: wf(1:gr%mesh%np_part)
  integer,                intent(in) :: mm
  type(modelmb_denmat_t), intent(in) :: denmat

  integer :: jj, ll, j, err_code, iunit, ndims, ndim1part
  integer :: ikeeppart, idir
  integer :: idensmat, nparticles
  integer, allocatable :: npoints(:)
  integer, allocatable :: ix_1part(:), ix_1part_p(:)
  logical :: bof
  character(len=200) :: filename
  R_TYPE, allocatable :: densmatr(:, :), evectors(:, :)
  R_TYPE, allocatable :: densmatr_tmp(:, :)
  FLOAT, allocatable :: evalues(:), density(:)

  type(modelmb_1part_t) :: mb_1part
  FLOAT, allocatable :: dipole_moment(:)

  PUSH_SUB(X(modelmb_density_matrix_write))


  ! The algorithm should consider how many dimensions the wavefunction has (ndims),
  ! and how many (and which) dimensions should be integrated away.
  ndims = gr%sb%dim

  ndim1part=st%modelmbparticles%ndim

  call modelmb_1part_nullify(mb_1part)
  SAFE_ALLOCATE(  ix_1part(1:ndim1part))
  SAFE_ALLOCATE(ix_1part_p(1:ndim1part))
  SAFE_ALLOCATE(dipole_moment(1:ndim1part))

  ! Allocation of the arrays that store the limiting indices for each direction
  SAFE_ALLOCATE(npoints(1:ndims))
  do j = 1, ndims
    npoints(j) = gr%mesh%idx%ll(j)
  end do


  ! loop over desired density matrices
  densmat_loop: do idensmat = 1, denmat%ndensmat_to_calculate
    ikeeppart = denmat%particle_kept(idensmat)
    nparticles = st%modelmbparticles%nparticles_per_type(st%modelmbparticles%particletype(ikeeppart))

    call modelmb_1part_init(mb_1part, gr%mesh, ikeeppart, ndim1part, gr%sb%box_offset)

    SAFE_ALLOCATE(densmatr(1:mb_1part%npt, 1:mb_1part%npt))
    SAFE_ALLOCATE(evectors(1:mb_1part%npt, 1:mb_1part%npt))
    SAFE_ALLOCATE(evalues(1:mb_1part%npt))
    SAFE_ALLOCATE(density(1:mb_1part%npt))

    
    densmatr  = R_TOTYPE(M_ZERO)

    !   calculate the 1-particle density matrix for this many-body state, and for the chosen
    !   particle being the free coordinate
    call X(mf_calculate_gamma)(ikeeppart, mb_1part, nparticles, &
          gr%mesh, wf, densmatr)

    ! Only node zero writes.
    ! mjv 14/3/2009: is this still at the right place in the file? None of
    ! this works in parallel yet...
    if(.not. mpi_grp_is_root(mpi_world)) cycle

    !Diagonalize the density matrix
    bof=.true.
    SAFE_ALLOCATE(densmatr_tmp(1:mb_1part%npt, 1:mb_1part%npt))
    densmatr_tmp=densmatr
    ! CHECK: should we only be diagonalizing the main grid points, as
    ! opposed to the full mb_1part%npt?
    evectors = densmatr_tmp
    call lalg_eigensolve(mb_1part%npt, evectors, evalues, bof, err_code)
    SAFE_DEALLOCATE_A(densmatr_tmp)
  
    !NOTE: The highest eigenvalues are the last ones not the first!!!
    !      Writing is therefore in reverse order
    evectors = evectors/sqrt(mb_1part%vol_elem_1part)
    evalues  = evalues*mb_1part%vol_elem_1part

    !Write everything into files
    write(filename,'(a,i3.3,a,i2.2)') trim(denmat%dirname)//'/occnumb_ip',ikeeppart,'_imb',mm
    iunit = io_open(trim(filename), action='write')

    do jj = mb_1part%npt, 1, -1
      write(iunit,'(i4.4,es11.3)') mb_1part%npt-jj+1, evalues(jj)
    end do
        
    call io_close(iunit)

    do jj = mb_1part%npt-denmat%nnatorb_prt(idensmat)+1, mb_1part%npt
      write(filename,'(a,i3.3,a,i2.2,a,i4.4)') trim(denmat%dirname)//'/natorb_ip', &
            ikeeppart,'_imb', mm, '_', mb_1part%npt-jj+1
      iunit = io_open(filename, action='write')
      do ll = 1, mb_1part%npt
        call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, &
             mb_1part%enlarge_1part(1), ll, ix_1part)
        do idir=1,ndim1part
          write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
        end do
        write(iunit,'(es11.3,es11.3)') evectors(ll,jj) 
        !), aimag(evectors(ll,jj)) ! format is too long for real wf case, but should be ok for most compilers
      end do
    call io_close(iunit)
    end do

    write(filename,'(a,i3.3,a,i2.2)') trim(denmat%dirname)//'/densmatr_ip', ikeeppart,'_imb', mm
    iunit = io_open(filename,action='write')
    do jj = 1, mb_1part%npt
      call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, &
           mb_1part%enlarge_1part(1), jj, ix_1part)
      do ll = 1, mb_1part%npt
        call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, &
             mb_1part%enlarge_1part(1), ll, ix_1part_p)
        do idir=1,ndim1part
          write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
        end do
        do idir=1,ndim1part
          write(iunit,'(es11.3)', ADVANCE='no') ix_1part_p(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
        end do
        write(iunit,'(es11.3,es11.3)') densmatr(jj,ll)
        !), aimag(densmatr(jj,ll)) ! format is too long for real wf case, but should be ok for most compilers
      end do
      write(iunit,*)
    end do
    call io_close(iunit)

    write(filename,'(a,i3.3,a,i2.2)') trim(denmat%dirname)//'/density_ip', ikeeppart,'_imb', mm
    iunit = io_open(filename,action='write')
    do jj = 1, mb_1part%npt
      call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, &
           mb_1part%enlarge_1part(1), jj, ix_1part)
      do idir=1,ndim1part
        write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
      end do
      write(iunit,'(es18.10)') real(densmatr(jj,jj))
    end do
    call io_close(iunit)


    ! calculate dipole moment from density for this particle
    dipole_moment(:) = M_ZERO
    do jj = 1,mb_1part%npt
      call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, &
           mb_1part%enlarge_1part(1), jj, ix_1part)
      dipole_moment = dipole_moment+(ix_1part(:)*mb_1part%h_1part(:)+mb_1part%origin(:))&
                    *TOFLOAT(densmatr(jj,jj))&
                    *st%modelmbparticles%charge_particle(ikeeppart)
    end do
    ! note: for eventual multiple particles in 4D (eg 8D total) this would fail to give the last values of dipole_moment
    write (message(1),'(a,I6,a,I6,a,I6)') 'For particle ', ikeeppart, ' of mb state ', mm
    write (message(2),'(a,3E20.10)') 'The dipole moment is (in a.u. = e bohr):     ', dipole_moment(1:min(3,ndim1part))
    write (message(3),'(a,E15.3)') '     with intrinsic numerical error usually <= ', 1.e-6*mb_1part%npt
    call messages_info(3)

    SAFE_DEALLOCATE_A(evectors)
    SAFE_DEALLOCATE_A(evalues)
    SAFE_DEALLOCATE_A(density)
    SAFE_DEALLOCATE_A(densmatr)
  
    call modelmb_1part_end(mb_1part)
  
  end do densmat_loop ! loop over densmats to output


  SAFE_DEALLOCATE_A(ix_1part)
  SAFE_DEALLOCATE_A(ix_1part_p)
  SAFE_DEALLOCATE_A(npoints)
  SAFE_DEALLOCATE_A(dipole_moment)

  POP_SUB(X(modelmb_density_matrix_write))
end subroutine X(modelmb_density_matrix_write)
! ---------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
