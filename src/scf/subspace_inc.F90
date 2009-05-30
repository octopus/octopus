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

! ---------------------------------------------------------
! This routine diagonalises the Hamiltonian in the subspace defined by the states.
subroutine X(subspace_diag)(gr, st, hm, ik, diff)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  integer,             intent(in)    :: ik
  FLOAT, optional,     intent(out)   :: diff(:)

  R_TYPE, allocatable :: h_subspace(:, :), vec(:, :), f(:, :, :), psi(:, :, :)
  integer             :: ist, ist2, jst, size, idim
  FLOAT               :: nrm2
  type(profile_t),     save    :: diagon_prof
  type(batch_t) :: psib, hpsib, whole_psib

  call push_sub('subspace_inc.Xsubspace_diag')
  call profiling_in(diagon_prof, "SUBSPACE_DIAG")

#ifdef HAVE_MPI
  if(st%parallel_in_states) then
    if(present(diff)) then
      call X(subspace_diag_par_states)(gr, st, hm, ik, diff)
    else
      call X(subspace_diag_par_states)(gr, st, hm, ik)
    end if
  else
#endif

    SAFE_ALLOCATE(h_subspace(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(       vec(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(f(1:gr%mesh%np, 1:st%d%dim, 1:st%d%block_size))

    if(st%np_size) then
      SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%block_size))
    end if

    ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
    do ist = st%st_start, st%st_end, st%d%block_size
      size = min(st%d%block_size, st%st_end - ist + 1)

      if(st%np_size) then
        call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi)
        call batch_set(psib, gr%mesh%np, st%X(psi)(:, :, ist:, ik))
      else
        call batch_init(psib, hm%d%dim, ist, ist + size - 1, st%X(psi)(:, :, ist:, ik))
      end if

      call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, f)

      call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik)

      call batch_end(psib)

      call batch_init(whole_psib, hm%d%dim, ist, st%nst, st%X(psi)(:, :, ist:st%nst, ik))
      call X(mesh_batch_dotp_matrix)(gr%mesh, hpsib, whole_psib, h_subspace)
      call batch_end(whole_psib)
      call batch_end(hpsib)
      
    end do

    do ist = st%st_start, st%st_end
      do jst = st%st_start, st%st_end
        h_subspace(jst, ist) = R_CONJ(h_subspace(ist, jst))
      end do
    end do

    ! Diagonalize the hamiltonian in the subspace.
    call lalg_eigensolve(st%nst, h_subspace, vec, st%eigenval(:, ik))

    ! Calculate the new eigenfunctions as a linear combination of the
    ! old ones.
    call batch_init(whole_psib, hm%d%dim, 1, st%nst, st%X(psi)(:, :, :, ik))
    call X(mesh_batch_rotate)(gr%mesh, whole_psib, vec)
    call batch_end(whole_psib)

    ! Renormalize.
    do ist = st%st_start, st%st_end
      nrm2 = X(mf_nrm2)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik))
      do idim = 1, st%d%dim
        call lalg_scal(gr%mesh%np, M_ONE/nrm2, st%X(psi)(:, idim, ist, ik))
      end do
    end do

    ! Recalculate the residues if requested by the diff argument.
    if(present(diff)) then 
      
      do ist = st%st_start, st%st_end, st%d%block_size
        size = min(st%d%block_size, st%st_end - ist + 1)

        if(st%np_size) then
          call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi)
          call batch_set(psib, gr%mesh%np, st%X(psi)(:, :, ist:, ik))
        else
          call batch_init(psib, hm%d%dim, ist, ist + size - 1, st%X(psi)(:, :, ist:, ik))
        end if
        call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, f)
        
        call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik)

        call batch_end(psib)
        call batch_end(hpsib)
        
        do ist2 = ist, ist + size - 1
          diff(ist2) = X(states_residue)(gr%mesh, st%d%dim, f(:, :, ist2 - ist + 1), st%eigenval(ist2, ik), st%X(psi)(:, :, ist2, ik))
        end do
      end do

    end if

    if(st%np_size) then
      SAFE_DEALLOCATE_A(psi)
    end if

    SAFE_DEALLOCATE_A(f)
    SAFE_DEALLOCATE_A(h_subspace)
    SAFE_DEALLOCATE_A(vec)

#ifdef HAVE_MPI
  end if
#endif

  call profiling_out(diagon_prof)
  call pop_sub()

end subroutine X(subspace_diag)

#ifdef HAVE_MPI
! --------------------------------------------------------- 
! This routine diagonalises the Hamiltonian in the subspace defined by
! the states, this version is aware of parallelization in states but
! consumes more memory.
!
subroutine X(subspace_diag_par_states)(gr, st, hm, ik, diff)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  integer,             intent(in)    :: ik
  FLOAT, optional,     intent(out)   :: diff(1:st%nst)

  R_TYPE, allocatable :: h_subspace(:,:), vec(:,:), f(:,:,:)
  integer             :: i
  FLOAT               :: nrm2
#if defined(HAVE_MPI)
  integer             :: tmp
  FLOAT               :: ldiff(st%lnst)
#endif

  call push_sub('subspace_inc.Xsubspace_diag_par_states')

  SAFE_ALLOCATE(h_subspace(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(       vec(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(f(1:gr%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))

  ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
  do i = st%st_start, st%st_end
    call X(hamiltonian_apply)(hm, gr, st%X(psi)(:, :, i, ik), f(:, :, i), i, ik)
  end do
  call states_blockt_mul(gr%mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
       st%X(psi)(:, :, :, ik), f, h_subspace, symm=.true.)

  ! Diagonalize the hamiltonian in the subspace.
  call lalg_eigensolve(st%nst, h_subspace, vec, st%eigenval(:, ik))

  ! The new states are the given by the eigenvectors of the matrix.
  f(1:gr%mesh%np, 1:st%d%dim, st%st_start:st%st_end) = st%X(psi)(1:gr%mesh%np, 1:st%d%dim, st%st_start:st%st_end, ik)

  call states_block_matr_mul(gr%mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
       f, vec, st%X(psi)(:, :, :, ik))

  ! Renormalize.
  do i = st%st_start, st%st_end
    nrm2 = X(mf_nrm2)(gr%mesh, st%d%dim, st%X(psi)(:, :, i, ik))
    st%X(psi)(1:gr%mesh%np, 1:st%d%dim, i, ik) = st%X(psi)(1:gr%mesh%np, 1:st%d%dim, i, ik)/nrm2
  end do

  ! Recalculate the residues if requested by the diff argument.
  if(present(diff)) then 
    do i = st%st_start, st%st_end
      call X(hamiltonian_apply)(hm, gr, st%X(psi)(:, :, i, ik) , f(:, :, st%st_start), i, ik)
      diff(i) = X(states_residue)(gr%mesh, st%d%dim, f(:, :, st%st_start), st%eigenval(i, ik), &
           st%X(psi)(:, :, i, ik))
    end do

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      ldiff = diff(st%st_start:st%st_end)
      call lmpi_gen_allgatherv(st%lnst, ldiff, tmp, diff(:), st%mpi_grp)
    end if
#endif
  end if
  
  SAFE_DEALLOCATE_A(f)
  SAFE_DEALLOCATE_A(h_subspace)
  SAFE_DEALLOCATE_A(vec)
  
  call pop_sub()
  
end subroutine X(subspace_diag_par_states) 
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
