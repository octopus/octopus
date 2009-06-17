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
subroutine X(subspace_diag)(gr, st, hm, ik, eigenval, psi, diff)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st
  type(hamiltonian_t), intent(inout) :: hm
  integer,             intent(in)    :: ik
  FLOAT,               intent(out)   :: eigenval(:)
  R_TYPE,              intent(inout) :: psi(:, :, :)
  FLOAT, optional,     intent(out)   :: diff(:)

  R_TYPE, allocatable :: h_subspace(:, :), f(:, :, :), psi2(:, :, :)
  integer             :: ist, ist2, jst, size, idim
  FLOAT               :: nrm2
  type(profile_t),     save    :: diagon_prof
  type(batch_t) :: psib, hpsib, whole_psib

  call push_sub('subspace_inc.Xsubspace_diag')
  call profiling_in(diagon_prof, "SUBSPACE_DIAG")

#ifdef HAVE_MPI
  if(st%parallel_in_states) then
    if(present(diff)) then
      call X(subspace_diag_par_states)(gr, st, hm, ik, psi, diff)
    else
      call X(subspace_diag_par_states)(gr, st, hm, ik, psi)
    end if
  else
#endif

    SAFE_ALLOCATE(h_subspace(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(f(1:gr%mesh%np, 1:st%d%dim, 1:st%d%block_size))

    if(st%np_size) then
      SAFE_ALLOCATE(psi2(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%block_size))
    end if

    ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
    do ist = st%st_start, st%st_end, st%d%block_size
      size = min(st%d%block_size, st%st_end - ist + 1)

      if(st%np_size) then
        call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi2)
        call batch_set(psib, gr%mesh%np, psi(:, :, ist:))
      else
        call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi(:, :, ist:))
      end if

      call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, f)

      call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik)

      call batch_end(psib)

      call batch_init(whole_psib, hm%d%dim, ist, st%nst, psi(:, :, ist:st%nst))
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
    call lalg_eigensolve(st%nst, h_subspace, eigenval(:))

    ! Calculate the new eigenfunctions as a linear combination of the
    ! old ones.
    call batch_init(whole_psib, hm%d%dim, 1, st%nst, psi(:, :, :))
    call X(mesh_batch_rotate)(gr%mesh, whole_psib, h_subspace)
    call batch_end(whole_psib)

    ! Renormalize.
    do ist = st%st_start, st%st_end
      nrm2 = X(mf_nrm2)(gr%mesh, st%d%dim, psi(:, :, ist))
      do idim = 1, st%d%dim
        call lalg_scal(gr%mesh%np, M_ONE/nrm2, psi(:, idim, ist))
      end do
    end do

    ! Recalculate the residues if requested by the diff argument.
    if(present(diff)) then 
      
      do ist = st%st_start, st%st_end, st%d%block_size
        size = min(st%d%block_size, st%st_end - ist + 1)

        if(st%np_size) then
          call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi2)
          call batch_set(psib, gr%mesh%np, psi(:, :, ist:))
        else
          call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi(:, :, ist:))
        end if
        call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, f)
        
        call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik)

        call batch_end(psib)
        call batch_end(hpsib)
        
        do ist2 = ist, ist + size - 1
          diff(ist2) = X(states_residue)(gr%mesh, st%d%dim, f(:, :, ist2 - ist + 1), eigenval(ist2), psi(:, :, ist2))
        end do
      end do

    end if

    if(st%np_size) then
      SAFE_DEALLOCATE_A(psi2)
    end if

    SAFE_DEALLOCATE_A(f)
    SAFE_DEALLOCATE_A(h_subspace)

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
subroutine X(subspace_diag_par_states)(gr, st, hm, ik, psi, diff)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: psi(:, :, :)
  FLOAT, optional,     intent(out)   :: diff(1:st%nst)

  R_TYPE, allocatable :: h_subspace(:,:), f(:,:,:)
  integer             :: i
  FLOAT               :: nrm2
#if defined(HAVE_MPI)
  integer             :: tmp
  FLOAT               :: ldiff(st%lnst)
#endif

  call push_sub('subspace_inc.Xsubspace_diag_par_states')

  SAFE_ALLOCATE(h_subspace(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(f(1:gr%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))

  ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
  do i = st%st_start, st%st_end
    call X(hamiltonian_apply)(hm, gr, psi(:, :, i), f(:, :, i), i, ik)
  end do
  call states_blockt_mul(gr%mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
       psi(:, :, :), f, h_subspace, symm=.true.)

  ! Diagonalize the hamiltonian in the subspace.
  call lalg_eigensolve(st%nst, h_subspace, eigenval(:))

  ! The new states are the given by the eigenvectors of the matrix.
  f(1:gr%mesh%np, 1:st%d%dim, st%st_start:st%st_end) = psi(1:gr%mesh%np, 1:st%d%dim, st%st_start:st%st_end)

  call states_block_matr_mul(gr%mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
       f, h_subspace, psi(:, :, :))

  ! Renormalize.
  do i = st%st_start, st%st_end
    nrm2 = X(mf_nrm2)(gr%mesh, st%d%dim, psi(:, :, i))
    psi(1:gr%mesh%np, 1:st%d%dim, i) = psi(1:gr%mesh%np, 1:st%d%dim, i)/nrm2
  end do

  ! Recalculate the residues if requested by the diff argument.
  if(present(diff)) then 
    do i = st%st_start, st%st_end
      call X(hamiltonian_apply)(hm, gr, psi(:, :, i) , f(:, :, st%st_start), i, ik)
      diff(i) = X(states_residue)(gr%mesh, st%d%dim, f(:, :, st%st_start), eigenval(i), &
           psi(:, :, i))
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
  
  call pop_sub()
  
end subroutine X(subspace_diag_par_states) 
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
